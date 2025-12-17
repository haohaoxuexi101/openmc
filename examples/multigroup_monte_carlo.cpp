#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __CUDACC__
#include <cuda_runtime.h>
#endif

namespace {

constexpr int kGroups = 7;
constexpr int kMaterials = 2; // 0 fuel, 1 moderator

#ifdef __CUDACC__
#define HD __host__ __device__
#define H __host__
#define D __device__
#else
#define HD
#define H
#define D
#endif

struct Vec3 {
  double x {0.0};
  double y {0.0};
  double z {0.0};
};

struct Particle {
  Vec3 r;          // position [cm]
  Vec3 u;          // direction (unit vector)
  int g {0};       // energy group (0-based)
  double w {1.0};  // weight
};

struct Mgxs {
  double total[kGroups];
  double absorption[kGroups];
  double fission[kGroups];
  double nu_fission[kGroups];
  double chi[kGroups];
  double scatter[kGroups * kGroups]; // row-major: g_out major, g_in minor (GxG)
};

struct MgxsView {
  const double* total;
  const double* absorption;
  const double* fission;
  const double* nu_fission;
  const double* chi;
  const double* scatter;
};

struct Tallies {
  std::vector<double> flux_tracklength; // size: materials * groups
  std::vector<double> nu_fission;
  int groups {0};

  explicit Tallies(int materials, int groups_)
    : flux_tracklength(materials * groups_, 0.0), nu_fission(materials * groups_, 0.0),
      groups(groups_) {}

  void accumulate_flux(int material, int g, double track) { flux_tracklength[material * groups + g] += track; }
  void accumulate_nu_fission(int material, int g, double contrib) { nu_fission[material * groups + g] += contrib; }
};

struct TalliesView {
  double* flux_tracklength;
  double* nu_fission;
};

struct Rng {
  // Counter-based SplitMix64 for deterministic parallel streams.
  uint64_t state;
  HD explicit Rng(uint64_t seed) : state(seed + 0x9E3779B97F4A7C15ull) {}

  HD uint64_t next_uint64()
  {
    uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }

  HD double uniform() { return (next_uint64() >> 11) * (1.0 / 9007199254740992.0); }
};

HD Vec3 sample_isotropic(Rng& rng)
{
  const double mu = 2.0 * rng.uniform() - 1.0;
  const double phi = 2.0 * M_PI * rng.uniform();
  const double sin_t = std::sqrt(fmax(0.0, 1.0 - mu * mu));
  return {sin_t * std::cos(phi), sin_t * std::sin(phi), mu};
}

HD int sample_energy_group(const double chi[kGroups], Rng& rng)
{
  const double xi = rng.uniform();
  double cdf = 0.0;
  for (int g = 0; g < kGroups; ++g) {
    cdf += chi[g];
    if (xi <= cdf || g + 1 == kGroups) return g;
  }
  return kGroups - 1;
}

HD int sample_scatter_group(const MgxsView& xs, int g_in, Rng& rng)
{
  const double xi = rng.uniform();
  double cdf = 0.0;
  for (int g_out = 0; g_out < kGroups; ++g_out) {
    cdf += xs.scatter[g_out * kGroups + g_in] / xs.total[g_in];
    if (xi <= cdf || g_out + 1 == kGroups) return g_out;
  }
  return kGroups - 1;
}

Mgxs make_uo2()
{
  Mgxs xs {};
  const double total[kGroups] {0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678, 0.5644058};
  const double absorption[kGroups] {8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02, 3.0020E-02, 1.1126E-01, 2.8278E-01};
  const double fission[kGroups] {7.21206E-03, 8.19301E-04, 6.45320E-03, 1.85648E-02, 1.78084E-02, 8.30348E-02, 2.16004E-01};
  const double nu_fission[kGroups] {2.005998E-02, 2.027303E-03, 1.570599E-02, 4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01};
  const double chi[kGroups] {5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0, 0.0, 0.0};
  const double scatter[kGroups * kGroups] {
    0.1275370, 0.0423780, 0.0000094, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.3244560, 0.0016314, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4509400, 0.0026792, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.4525650, 0.0055664, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0001253, 0.2714010, 0.0102550, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0012968, 0.2658020, 0.0168090,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0085458, 0.2730800};
  for (int g = 0; g < kGroups; ++g) {
    xs.total[g] = total[g];
    xs.absorption[g] = absorption[g];
    xs.fission[g] = fission[g];
    xs.nu_fission[g] = nu_fission[g];
    xs.chi[g] = chi[g];
    for (int gin = 0; gin < kGroups; ++gin)
      xs.scatter[g * kGroups + gin] = scatter[g * kGroups + gin];
  }
  return xs;
}

Mgxs make_water()
{
  Mgxs xs {};
  const double total[kGroups] {0.15920605, 0.412969593, 0.59030986, 0.58435, 0.718, 1.2544497, 2.650379};
  const double absorption[kGroups] {6.0105E-04, 1.5793E-05, 3.3716E-04, 1.9406E-03, 5.7416E-03, 1.5001E-02, 3.7239E-02};
  const double chi[kGroups] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double fission[kGroups] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double nu_fission[kGroups] {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  const double scatter[kGroups * kGroups] {
    0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0, 0.0,
    0.0, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010,
    0.0, 0.0, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034,
    0.0, 0.0, 0.0, 0.0910284, 0.4155100, 0.0637320, 0.0121390,
    0.0, 0.0, 0.0, 0.0000714, 0.1391380, 0.5118200, 0.0612290,
    0.0, 0.0, 0.0, 0.0, 0.0022157, 0.6999130, 0.5373200,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.1324400, 2.4807000};
  for (int g = 0; g < kGroups; ++g) {
    xs.total[g] = total[g];
    xs.absorption[g] = absorption[g];
    xs.fission[g] = fission[g];
    xs.nu_fission[g] = nu_fission[g];
    xs.chi[g] = chi[g];
    for (int gin = 0; gin < kGroups; ++gin)
      xs.scatter[g * kGroups + gin] = scatter[g * kGroups + gin];
  }
  return xs;
}

struct SurfaceHit { double distance; int axis; };
struct FuelHit { double distance; bool leaving_fuel; };

HD Vec3 lattice_index(const Vec3& r, const Vec3& lower, const Vec3& pitch)
{
  return {(r.x - lower.x) / pitch.x, (r.y - lower.y) / pitch.y, (r.z - lower.z) / pitch.z};
}

HD int universe_id_from_lattice(const Vec3& r)
{
  // 3x3x3 lattice; layout follows earlier XML example.
  const Vec3 lower {-1.89, -1.89, -150.0};
  const Vec3 pitch {1.26, 1.26, 100.0};
  Vec3 idx = lattice_index(r, lower, pitch);
  int i = static_cast<int>(std::floor(idx.x));
  int j = static_cast<int>(std::floor(idx.y));
  int k = static_cast<int>(std::floor(idx.z));
  if (i < 0 || i >= 3 || j < 0 || j >= 3 || k < 0 || k >= 3) return -1;

  // z-planes: [k=0 -> moderator], [k=1 -> center fuel], [k=2 -> moderator]
  if (k == 1 && i == 1 && j == 1) return 1; // central fuel pin
  if (k == 1 && (i == 1 || j == 1)) return 2; // water surrounding central column/row
  return 2; // blanket moderator elsewhere
}

HD bool inside_fuel_pin(const Vec3& r)
{
  // Central cylindrical fuel pellet of radius 0.54 cm centered in its pin cell.
  const Vec3 lower {-1.89, -1.89, -150.0};
  const Vec3 pitch {1.26, 1.26, 100.0};
  Vec3 idx = lattice_index(r, lower, pitch);
  const int i = static_cast<int>(std::floor(idx.x));
  const int j = static_cast<int>(std::floor(idx.y));
  const int k = static_cast<int>(std::floor(idx.z));
  if (i != 1 || j != 1 || k != 1) return false;
  const Vec3 center {lower.x + (i + 0.5) * pitch.x, lower.y + (j + 0.5) * pitch.y, r.z};
  const double dx = r.x - center.x;
  const double dy = r.y - center.y;
  return dx * dx + dy * dy <= 0.54 * 0.54;
}

HD int material_from_position(const Vec3& r)
{
  const int uni = universe_id_from_lattice(r);
  if (uni == 1) {
    return inside_fuel_pin(r) ? 0 : 1; // 0: fuel, 1: moderator
  }
  return uni == 2 ? 1 : -1;
}

HD SurfaceHit distance_to_surfaces(const Particle& p)
{
  // Reflective box: [-1.89,1.89] in x,y; [-150,150] in z.
  const double bounds[3][2] {{-1.89, 1.89}, {-1.89, 1.89}, {-150.0, 150.0}};
  double min_d = std::numeric_limits<double>::infinity();
  int axis = -1;
  const double dir[3] {p.u.x, p.u.y, p.u.z};
  const double pos[3] {p.r.x, p.r.y, p.r.z};
  for (int a = 0; a < 3; ++a) {
    if (std::abs(dir[a]) < 1e-12) continue;
    for (int s = 0; s < 2; ++s) {
      const double d = (bounds[a][s] - pos[a]) / dir[a];
      if (d > 1e-10 && d < min_d) {
        min_d = d;
        axis = a;
      }
    }
  }
  return {min_d, axis};
}

HD void reflect(Particle& p, int axis)
{
  if (axis == 0) p.u.x *= -1.0;
  if (axis == 1) p.u.y *= -1.0;
  if (axis == 2) p.u.z *= -1.0;
}

HD FuelHit distance_to_fuel_boundary(const Particle& p)
{
  // Only meaningful inside the central pin cell; returns distance and whether exiting fuel (true) or entering fuel (false).
  const Vec3 lower {-1.89, -1.89, -150.0};
  const Vec3 pitch {1.26, 1.26, 100.0};
  Vec3 idx = lattice_index(p.r, lower, pitch);
  const int i = static_cast<int>(std::floor(idx.x));
  const int j = static_cast<int>(std::floor(idx.y));
  const int k = static_cast<int>(std::floor(idx.z));
  if (i != 1 || j != 1 || k != 1) return {std::numeric_limits<double>::infinity(), false};

  const Vec3 center {lower.x + 1.5 * pitch.x, lower.y + 1.5 * pitch.y, p.r.z};
  const double dx = p.r.x - center.x;
  const double dy = p.r.y - center.y;
  const double b = dx * p.u.x + dy * p.u.y;
  const double c = dx * dx + dy * dy - 0.54 * 0.54;
  const double disc = b * b - c;
  if (disc <= 0.0) return {std::numeric_limits<double>::infinity(), false};
  const double sqrt_disc = std::sqrt(disc);
  double d1 = -b - sqrt_disc;
  double d2 = -b + sqrt_disc;
  double d_exit = std::numeric_limits<double>::infinity();
  if (d1 > 1e-10) d_exit = d1;
  else if (d2 > 1e-10) d_exit = d2;
  const bool leaving_fuel = inside_fuel_pin(p.r);
  return {d_exit, leaving_fuel};
}

HD MgxsView view_xs(const Mgxs& xs)
{
  return {xs.total, xs.absorption, xs.fission, xs.nu_fission, xs.chi, xs.scatter};
}

HD Particle sample_source(Rng& rng, const Mgxs& fuel)
{
  // Sample uniformly in fuel region only (central pin) to keep source fissionable.
  Particle p;
  p.g = sample_energy_group(fuel.chi, rng);
  p.u = sample_isotropic(rng);

  const Vec3 lower {-1.89, -1.89, -150.0};
  const Vec3 pitch {1.26, 1.26, 100.0};
  Vec3 center {lower.x + 1.5 * pitch.x, lower.y + 1.5 * pitch.y, 0.0};
  // Rejection sample inside cylinder within central cell.
  while (true) {
    const double rx = (rng.uniform() - 0.5) * pitch.x;
    const double ry = (rng.uniform() - 0.5) * pitch.y;
    if (rx * rx + ry * ry <= 0.54 * 0.54) {
      p.r.x = center.x + rx;
      p.r.y = center.y + ry;
      break;
    }
  }
  p.r.z = lower.z + pitch.z * (1.0 + rng.uniform()); // always middle lattice plane
  return p;
}

void transport_history(Particle& p, const Mgxs* materials, Tallies& tallies, Rng& rng,
  double& produced_neutrons, std::vector<Particle>& fission_bank)
{
  const MgxsView xs_set[kMaterials] {view_xs(materials[0]), view_xs(materials[1])};
  while (true) {
    const int m = material_from_position(p.r);
    if (m < 0) return; // leaked (should not happen with reflective boundaries)
    const MgxsView xs = xs_set[m];

    const double sig_t = xs.total[p.g];
    const double free_path = -std::log(rng.uniform()) / sig_t;
    SurfaceHit surf = distance_to_surfaces(p);
    FuelHit fuel = distance_to_fuel_boundary(p);
    const double min_surface = std::min(surf.distance, fuel.distance);

    if (free_path < min_surface) {
      // Collision occurs before any surface.
      const double track = free_path;
      p.r.x += track * p.u.x;
      p.r.y += track * p.u.y;
      p.r.z += track * p.u.z;
      tallies.accumulate_flux(m, p.g, track * p.w);

      const double xi = rng.uniform();
      const double scatter_prob = (xs.total[p.g] - xs.absorption[p.g]) / xs.total[p.g];
      const double fission_prob = xs.fission[p.g] / xs.total[p.g];

      if (xi < scatter_prob) {
        // Scatter
        p.g = sample_scatter_group(xs, p.g, rng);
        p.u = sample_isotropic(rng);
        continue;
      }

      if (xi < scatter_prob + fission_prob && xs.nu_fission[p.g] > 0.0) {
        const double nu = xs.nu_fission[p.g];
        const int num_new = static_cast<int>(std::floor(nu + rng.uniform()));
        tallies.accumulate_nu_fission(m, p.g, nu * p.w);
        produced_neutrons += nu * p.w;
        for (int n = 0; n < num_new; ++n) {
          Particle q;
          q.r = p.r;
          q.u = sample_isotropic(rng);
          q.g = sample_energy_group(xs.chi, rng);
          q.w = p.w;
          fission_bank.push_back(q);
        }
        return; // particle absorbed in fission
      }

      // Absorption without fission
      return;
    }

    // Surface interaction happens before collision.
    const double track = min_surface;
    p.r.x += track * p.u.x;
    p.r.y += track * p.u.y;
    p.r.z += track * p.u.z;
    tallies.accumulate_flux(m, p.g, track * p.w);

    if (surf.distance < fuel.distance) {
      reflect(p, surf.axis);
    }
  }
}

void run_eigenvalue_cpu(int batches, int inactive, int histories, uint64_t seed)
{
  Mgxs materials[kMaterials] {make_uo2(), make_water()};

  Tallies tallies(kMaterials, kGroups);

  std::vector<Particle> source_bank(histories);
  Rng seeder(seed);
  for (auto& p : source_bank) p = sample_source(seeder, materials[0]);

  for (int b = 0; b < batches; ++b) {
    Tallies batch_tallies(kMaterials, kGroups);
    double produced = 0.0;
    std::vector<Particle> next_bank;
    next_bank.reserve(histories);

#pragma omp parallel
    {
      std::vector<Particle> local_bank;
      double local_produced = 0.0;
      Tallies local_tallies(kMaterials, kGroups);

#pragma omp for schedule(static)
      for (int i = 0; i < histories; ++i) {
        Rng rng(seed + static_cast<uint64_t>(b) * histories + i);
        Particle p = source_bank[i];
        transport_history(p, materials, local_tallies, rng, local_produced, local_bank);
      }

#pragma omp critical
      {
        produced += local_produced;
        for (auto& q : local_bank) next_bank.push_back(q);
        for (size_t i = 0; i < local_tallies.flux_tracklength.size(); ++i)
          batch_tallies.flux_tracklength[i] += local_tallies.flux_tracklength[i];
        for (size_t i = 0; i < local_tallies.nu_fission.size(); ++i)
          batch_tallies.nu_fission[i] += local_tallies.nu_fission[i];
      }
    }

    const double keff = produced / static_cast<double>(histories);
    if (next_bank.empty()) throw std::runtime_error("All particles died; adjust cross sections or source.");
    // Resample to fixed population.
    std::vector<Particle> resampled;
    resampled.reserve(histories);
    Rng bank_rng(seed + static_cast<uint64_t>(b) * 7919);
    for (int i = 0; i < histories; ++i) {
      const size_t pick = static_cast<size_t>(bank_rng.uniform() * next_bank.size());
      resampled.push_back(next_bank[pick]);
    }
    source_bank.swap(resampled);

    if (b >= inactive) {
      for (size_t i = 0; i < tallies.flux_tracklength.size(); ++i)
        tallies.flux_tracklength[i] += batch_tallies.flux_tracklength[i];
      for (size_t i = 0; i < tallies.nu_fission.size(); ++i)
        tallies.nu_fission[i] += batch_tallies.nu_fission[i];
    }

    std::cout << "Batch " << std::setw(3) << b + 1 << " keff ~ " << std::setprecision(6) << keff <<
      " bank size " << next_bank.size() << std::endl;
  }

  std::cout << "\nTallied track-length flux by material/group:\n";
  for (int m = 0; m < kMaterials; ++m) {
    std::cout << (m == 0 ? "Fuel   : " : "Water  : ");
    for (int g = 0; g < kGroups; ++g) {
      const double val = tallies.flux_tracklength[m * kGroups + g];
      std::cout << std::setw(10) << std::setprecision(4) << val;
    }
    std::cout << "\n";
  }
}

#ifdef __CUDACC__
inline void check_cuda(cudaError_t err, const char* msg)
{
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
  }
}

struct DeviceBuffers {
  Particle* source {nullptr};
  Particle* next_bank {nullptr};
  int* bank_count {nullptr};
  double* produced_sum {nullptr};
  double* flux {nullptr};
  double* nu_fission {nullptr};
};

struct DeviceXsBuffers {
  double* total {nullptr};
  double* absorption {nullptr};
  double* fission {nullptr};
  double* nu_fission {nullptr};
  double* chi {nullptr};
  double* scatter {nullptr};
};

void allocate_device_xs(const Mgxs& host, DeviceXsBuffers& dev, MgxsView& view)
{
  check_cuda(cudaMalloc(&dev.total, sizeof(double) * kGroups), "alloc xs total");
  check_cuda(cudaMalloc(&dev.absorption, sizeof(double) * kGroups), "alloc xs absorption");
  check_cuda(cudaMalloc(&dev.fission, sizeof(double) * kGroups), "alloc xs fission");
  check_cuda(cudaMalloc(&dev.nu_fission, sizeof(double) * kGroups), "alloc xs nu_fission");
  check_cuda(cudaMalloc(&dev.chi, sizeof(double) * kGroups), "alloc xs chi");
  check_cuda(cudaMalloc(&dev.scatter, sizeof(double) * kGroups * kGroups), "alloc xs scatter");
  check_cuda(cudaMemcpy(dev.total, host.total, sizeof(double) * kGroups, cudaMemcpyHostToDevice), "copy xs total");
  check_cuda(cudaMemcpy(dev.absorption, host.absorption, sizeof(double) * kGroups, cudaMemcpyHostToDevice), "copy xs absorption");
  check_cuda(cudaMemcpy(dev.fission, host.fission, sizeof(double) * kGroups, cudaMemcpyHostToDevice), "copy xs fission");
  check_cuda(cudaMemcpy(dev.nu_fission, host.nu_fission, sizeof(double) * kGroups, cudaMemcpyHostToDevice), "copy xs nu_fission");
  check_cuda(cudaMemcpy(dev.chi, host.chi, sizeof(double) * kGroups, cudaMemcpyHostToDevice), "copy xs chi");
  check_cuda(cudaMemcpy(dev.scatter, host.scatter, sizeof(double) * kGroups * kGroups, cudaMemcpyHostToDevice), "copy xs scatter");
  view = {dev.total, dev.absorption, dev.fission, dev.nu_fission, dev.chi, dev.scatter};
}

D void atomic_accumulate(double* arr, int idx, double val)
{
  atomicAdd(&arr[idx], val);
}

D void atomic_increment(int* counter, int val)
{
  atomicAdd(counter, val);
}

D void append_fission(Particle* bank, int capacity, int* counter, const Particle& q, int copies)
{
  const int slot = atomicAdd(counter, copies);
  for (int n = 0; n < copies; ++n) {
    if (slot + n < capacity) bank[slot + n] = q;
  }
}

D void transport_history_device(Particle p, MgxsView fuel, MgxsView water, TalliesView tallies, Rng rng,
  int bank_capacity, Particle* bank, int* bank_counter, double* produced)
{
  MgxsView xs_set[kMaterials] {fuel, water};
  while (true) {
    const int m = material_from_position(p.r);
    if (m < 0) return;
    const MgxsView xs = xs_set[m];

    const double sig_t = xs.total[p.g];
    const double free_path = -log(rng.uniform()) / sig_t;
    SurfaceHit surf = distance_to_surfaces(p);
    FuelHit fuel_hit = distance_to_fuel_boundary(p);
    const double min_surface = fmin(surf.distance, fuel_hit.distance);

    if (free_path < min_surface) {
      const double track = free_path;
      p.r.x += track * p.u.x;
      p.r.y += track * p.u.y;
      p.r.z += track * p.u.z;
      atomic_accumulate(tallies.flux_tracklength, m * kGroups + p.g, track * p.w);

      const double xi = rng.uniform();
      const double scatter_prob = (xs.total[p.g] - xs.absorption[p.g]) / xs.total[p.g];
      const double fission_prob = xs.fission[p.g] / xs.total[p.g];

      if (xi < scatter_prob) {
        p.g = sample_scatter_group(xs, p.g, rng);
        p.u = sample_isotropic(rng);
        continue;
      }

      if (xi < scatter_prob + fission_prob && xs.nu_fission[p.g] > 0.0) {
        const double nu = xs.nu_fission[p.g];
        const int num_new = static_cast<int>(floor(nu + rng.uniform()));
        atomic_accumulate(tallies.nu_fission, m * kGroups + p.g, nu * p.w);
        atomicAdd(produced, nu * p.w);
        for (int n = 0; n < num_new; ++n) {
          Particle q;
          q.r = p.r;
          q.u = sample_isotropic(rng);
          q.g = sample_energy_group(xs.chi, rng);
          q.w = p.w;
          append_fission(bank, bank_capacity, bank_counter, q, 1);
        }
        return;
      }

      return;
    }

    const double track = min_surface;
    p.r.x += track * p.u.x;
    p.r.y += track * p.u.y;
    p.r.z += track * p.u.z;
    atomic_accumulate(tallies.flux_tracklength, m * kGroups + p.g, track * p.w);

    if (surf.distance < fuel_hit.distance) {
      reflect(p, surf.axis);
    }
  }
}

D void transport_kernel(const Particle* source, int histories, MgxsView fuel, MgxsView water,
  TalliesView tallies, Particle* bank, int bank_capacity, int* bank_counter, double* produced,
  uint64_t seed, int batch)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= histories) return;
  Rng rng(seed + static_cast<uint64_t>(batch) * histories + tid);
  Particle p = source[tid];
  transport_history_device(p, fuel, water, tallies, rng, bank_capacity, bank, bank_counter, produced);
}

void run_eigenvalue_cuda(int batches, int inactive, int histories, uint64_t seed)
{
  Mgxs fuel = make_uo2();
  Mgxs water = make_water();
  MgxsView fuel_v {};
  MgxsView water_v {};
  DeviceXsBuffers fuel_dev {};
  DeviceXsBuffers water_dev {};
  allocate_device_xs(fuel, fuel_dev, fuel_v);
  allocate_device_xs(water, water_dev, water_v);

  const int tally_size = kMaterials * kGroups;
  DeviceBuffers buffers;
  check_cuda(cudaMalloc(&buffers.source, sizeof(Particle) * histories), "alloc source");
  const int max_bank = histories * 8; // generous cap for demo purposes
  check_cuda(cudaMalloc(&buffers.next_bank, sizeof(Particle) * max_bank), "alloc bank");
  check_cuda(cudaMalloc(&buffers.bank_count, sizeof(int)), "alloc bank count");
  check_cuda(cudaMalloc(&buffers.produced_sum, sizeof(double)), "alloc produced");
  check_cuda(cudaMalloc(&buffers.flux, sizeof(double) * tally_size), "alloc flux");
  check_cuda(cudaMalloc(&buffers.nu_fission, sizeof(double) * tally_size), "alloc nu fission");

  std::vector<Particle> source_bank(histories);
  Rng seeder(seed);
  for (auto& p : source_bank) p = sample_source(seeder, fuel);

  std::vector<double> flux_host(tally_size, 0.0);
  std::vector<double> nu_fission_host(tally_size, 0.0);
  std::vector<double> flux_total(tally_size, 0.0);
  std::vector<double> nu_fission_total(tally_size, 0.0);

  for (int b = 0; b < batches; ++b) {
    check_cuda(cudaMemcpy(buffers.source, source_bank.data(), sizeof(Particle) * histories, cudaMemcpyHostToDevice), "copy source");
    check_cuda(cudaMemset(buffers.flux, 0, sizeof(double) * tally_size), "zero flux");
    check_cuda(cudaMemset(buffers.nu_fission, 0, sizeof(double) * tally_size), "zero nu fission");
    check_cuda(cudaMemset(buffers.bank_count, 0, sizeof(int)), "zero bank count");
    check_cuda(cudaMemset(buffers.produced_sum, 0, sizeof(double)), "zero produced");

    TalliesView tallies_v {buffers.flux, buffers.nu_fission};

    const int threads = 256;
    const int blocks = (histories + threads - 1) / threads;
    transport_kernel<<<blocks, threads>>>(buffers.source, histories, fuel_v, water_v, tallies_v,
      buffers.next_bank, max_bank, buffers.bank_count, buffers.produced_sum, seed, b);
    check_cuda(cudaDeviceSynchronize(), "run kernel");

    int bank_size = 0;
    check_cuda(cudaMemcpy(&bank_size, buffers.bank_count, sizeof(int), cudaMemcpyDeviceToHost), "copy bank size");
    std::vector<Particle> next_bank(bank_size);
    if (bank_size > 0) {
      check_cuda(cudaMemcpy(next_bank.data(), buffers.next_bank, sizeof(Particle) * bank_size, cudaMemcpyDeviceToHost), "copy bank");
    }
    double produced = 0.0;
    check_cuda(cudaMemcpy(&produced, buffers.produced_sum, sizeof(double), cudaMemcpyDeviceToHost), "copy produced");
    check_cuda(cudaMemcpy(flux_host.data(), buffers.flux, sizeof(double) * tally_size, cudaMemcpyDeviceToHost), "copy flux");
    check_cuda(cudaMemcpy(nu_fission_host.data(), buffers.nu_fission, sizeof(double) * tally_size, cudaMemcpyDeviceToHost), "copy nu fission");

    const double keff = produced / static_cast<double>(histories);
    if (next_bank.empty()) throw std::runtime_error("All particles died; adjust cross sections or source.");

    std::vector<Particle> resampled;
    resampled.reserve(histories);
    Rng bank_rng(seed + static_cast<uint64_t>(b) * 7919);
    for (int i = 0; i < histories; ++i) {
      const size_t pick = static_cast<size_t>(bank_rng.uniform() * next_bank.size());
      resampled.push_back(next_bank[pick]);
    }
    source_bank.swap(resampled);

    if (b >= inactive) {
      for (int i = 0; i < tally_size; ++i) {
        flux_total[i] += flux_host[i];
        nu_fission_total[i] += nu_fission_host[i];
      }
    }

    std::cout << "[CUDA] Batch " << std::setw(3) << b + 1 << " keff ~ " << std::setprecision(6) << keff <<
      " bank size " << next_bank.size() << std::endl;
  }

  std::cout << "\nTallied track-length flux by material/group (CUDA):\n";
  for (int m = 0; m < kMaterials; ++m) {
    std::cout << (m == 0 ? "Fuel   : " : "Water  : ");
    for (int g = 0; g < kGroups; ++g) {
      const double val = flux_total[m * kGroups + g];
      std::cout << std::setw(10) << std::setprecision(4) << val;
    }
    std::cout << "\n";
  }

  cudaFree(buffers.source);
  cudaFree(buffers.next_bank);
  cudaFree(buffers.bank_count);
  cudaFree(buffers.produced_sum);
  cudaFree(buffers.flux);
  cudaFree(buffers.nu_fission);
  cudaFree(fuel_dev.total);
  cudaFree(fuel_dev.absorption);
  cudaFree(fuel_dev.fission);
  cudaFree(fuel_dev.nu_fission);
  cudaFree(fuel_dev.chi);
  cudaFree(fuel_dev.scatter);
  cudaFree(water_dev.total);
  cudaFree(water_dev.absorption);
  cudaFree(water_dev.fission);
  cudaFree(water_dev.nu_fission);
  cudaFree(water_dev.chi);
  cudaFree(water_dev.scatter);
}
#endif

} // namespace

int main()
{
  try {
    const int batches = 20;
    const int inactive = 5;
    const int histories = 5000;
    const uint64_t seed = 12345;
#ifdef __CUDACC__
    run_eigenvalue_cuda(batches, inactive, histories, seed);
#else
    run_eigenvalue_cpu(batches, inactive, histories, seed);
#endif
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
    return 1;
  }
  return 0;
}
