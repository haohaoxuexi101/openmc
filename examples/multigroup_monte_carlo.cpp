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

// A self-contained multigroup Monte Carlo eigenvalue driver using the
// canonical C5G7 benchmark materials and a 17x17 PWR assembly geometry. The
// code supports both CPU (OpenMP) and CUDA execution paths while sharing
// identical transport, sampling, and cross-section data.

namespace {

// --- Global constants -------------------------------------------------------

constexpr int kGroups = 7;
constexpr int kMaterials = 7; // UO2, MOX4.3, MOX7, MOX8.7, guide tube, fission chamber, moderator

#ifdef __CUDACC__
#define HD __host__ __device__
#define H __host__
#define D __device__
#else
#define HD
#define H
#define D
#endif

// --- Basic types ------------------------------------------------------------

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
  double scatter[kGroups * kGroups]; // row-major g_out major, g_in minor
};

struct MgxsView {
  const double* total {nullptr};
  const double* absorption {nullptr};
  const double* fission {nullptr};
  const double* nu_fission {nullptr};
  const double* chi {nullptr};
  const double* scatter {nullptr};
};

struct Tallies {
  std::vector<double> flux_tracklength; // size materials * groups
  std::vector<double> nu_fission;       // size materials * groups
  int groups {0};

  explicit Tallies(int materials, int groups_)
    : flux_tracklength(materials * groups_, 0.0), nu_fission(materials * groups_, 0.0),
      groups(groups_) {}

  void accumulate_flux(int material, int g, double track) { flux_tracklength[material * groups + g] += track; }
  void accumulate_nu_fission(int material, int g, double contrib) { nu_fission[material * groups + g] += contrib; }
};

struct TalliesView {
  double* flux_tracklength {nullptr};
  double* nu_fission {nullptr};
};

// --- RNG utilities ----------------------------------------------------------

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

// --- Cross sections (C5G7) --------------------------------------------------

template<int N>
void copy_matrix(const double (&src)[N], double* dst)
{
  for (int i = 0; i < N; ++i) dst[i] = src[i];
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
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  copy_matrix(fission, xs.fission);
  copy_matrix(nu_fission, xs.nu_fission);
  copy_matrix(chi, xs.chi);
  copy_matrix(scatter, xs.scatter);
  return xs;
}

Mgxs make_mox43()
{
  Mgxs xs {};
  const double total[kGroups] {0.2067141, 0.3720612, 0.5461590, 0.7072748, 0.8303936, 1.0163840, 1.3887890};
  const double absorption[kGroups] {8.4339E-03, 8.2496E-03, 1.1531E-02, 3.3469E-02, 1.1760E-02, 1.0684E-01, 2.0088E-01};
  const double fission[kGroups] {7.62741E-03, 8.62873E-04, 5.82028E-03, 1.59407E-02, 1.31554E-02, 6.22917E-02, 1.56182E-01};
  const double nu_fission[kGroups] {2.25352E-02, 2.54835E-03, 1.73080E-02, 4.74910E-02, 3.93030E-02, 1.83715E-01, 4.25602E-01};
  const double chi[kGroups] {5.6443E-01, 4.3138E-01, 3.9706E-04, 1.2780E-07, 0.0, 0.0, 0.0};
  const double scatter[kGroups * kGroups] {
    0.1391880, 0.0636751, 0.0000071, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.2842863, 0.0035455, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4411510, 0.0037315, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.4556427, 0.0070729, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0001902, 0.2722965, 0.0122263, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0017017, 0.2711640, 0.0203723,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0100599, 0.2733531};
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  copy_matrix(fission, xs.fission);
  copy_matrix(nu_fission, xs.nu_fission);
  copy_matrix(chi, xs.chi);
  copy_matrix(scatter, xs.scatter);
  return xs;
}

Mgxs make_mox70()
{
  Mgxs xs {};
  const double total[kGroups] {0.2329149, 0.4090928, 0.5995840, 0.7846960, 0.9396500, 1.1071730, 1.4978630};
  const double absorption[kGroups] {1.0276E-02, 8.9707E-03, 1.2927E-02, 3.8145E-02, 1.2623E-02, 1.2568E-01, 2.4186E-01};
  const double fission[kGroups] {9.31340E-03, 9.41510E-04, 6.93910E-03, 1.79229E-02, 1.46028E-02, 7.45250E-02, 1.84205E-01};
  const double nu_fission[kGroups] {2.75466E-02, 2.78006E-03, 1.99786E-02, 5.37516E-02, 4.39470E-02, 2.24475E-01, 5.51603E-01};
  const double chi[kGroups] {5.8435E-01, 4.1454E-01, 3.8911E-04, 1.2526E-07, 0.0, 0.0, 0.0};
  const double scatter[kGroups * kGroups] {
    0.1482810, 0.0740283, 0.0000060, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.3078030, 0.0049670, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4485700, 0.0043862, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.4631300, 0.0084785, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0002533, 0.2800380, 0.0151235, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0019209, 0.2748570, 0.0245007,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0103100, 0.2742880};
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  copy_matrix(fission, xs.fission);
  copy_matrix(nu_fission, xs.nu_fission);
  copy_matrix(chi, xs.chi);
  copy_matrix(scatter, xs.scatter);
  return xs;
}

Mgxs make_mox87()
{
  Mgxs xs {};
  const double total[kGroups] {0.2656125, 0.4562789, 0.6582474, 0.8521908, 1.0586080, 1.1294030, 1.6162980};
  const double absorption[kGroups] {1.3244E-02, 9.2792E-03, 1.3550E-02, 4.0259E-02, 1.2711E-02, 1.4070E-01, 2.8130E-01};
  const double fission[kGroups] {1.17681E-02, 9.89892E-04, 7.83352E-03, 1.95275E-02, 1.51236E-02, 8.05580E-02, 2.01081E-01};
  const double nu_fission[kGroups] {3.42542E-02, 2.92601E-03, 2.16224E-02, 5.64131E-02, 4.38663E-02, 2.38571E-01, 5.94373E-01};
  const double chi[kGroups] {5.9232E-01, 4.0698E-01, 5.6992E-04, 1.9398E-07, 0.0, 0.0, 0.0};
  const double scatter[kGroups * kGroups] {
    0.1607090, 0.0826723, 0.0000084, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.3385210, 0.0065672, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4629600, 0.0050405, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.4783570, 0.0095188, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0003205, 0.2834870, 0.0170060, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0020084, 0.2783430, 0.0279491,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0108236, 0.2753690};
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  copy_matrix(fission, xs.fission);
  copy_matrix(nu_fission, xs.nu_fission);
  copy_matrix(chi, xs.chi);
  copy_matrix(scatter, xs.scatter);
  return xs;
}

Mgxs make_guide_tube()
{
  Mgxs xs {};
  const double total[kGroups] {0.1280039, 0.2204697, 0.4110794, 0.6739530, 1.3565440, 2.5091700, 3.0501960};
  const double absorption[kGroups] {2.9257E-04, 3.3750E-04, 6.1563E-04, 1.3425E-03, 2.5154E-03, 5.5440E-03, 1.29040E-02};
  const double scatter[kGroups * kGroups] {
    0.0872789, 0.0212057, 0.0000010, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.1991640, 0.0099369, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.3815910, 0.0275374, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.6121010, 0.0602343, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0003300, 0.8741970, 0.1393720, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0026673, 1.8327600, 0.6796830,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0654613, 2.9664220};
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  for (double& v : xs.fission) v = 0.0;
  for (double& v : xs.nu_fission) v = 0.0;
  for (int g = 0; g < kGroups; ++g) xs.chi[g] = 0.0;
  copy_matrix(scatter, xs.scatter);
  return xs;
}

Mgxs make_fission_chamber()
{
  Mgxs xs {};
  const double total[kGroups] {0.0395845, 0.1618260, 0.3832426, 0.6366331, 1.3060330, 2.3291780, 2.8384730};
  const double absorption[kGroups] {1.6434E-03, 1.2911E-03, 1.5146E-03, 2.4366E-03, 9.1171E-04, 2.5078E-02, 9.8135E-02};
  const double scatter[kGroups * kGroups] {
    0.0268930, 0.0100952, 0.0000016, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.1427090, 0.0174997, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.3473940, 0.0346911, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.5685810, 0.0639203, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0004130, 0.9197070, 0.1365960, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0020517, 1.8284900, 0.5975180,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0499530, 2.6750040};
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  for (double& v : xs.fission) v = 0.0;
  for (double& v : xs.nu_fission) v = 0.0;
  for (int g = 0; g < kGroups; ++g) xs.chi[g] = 0.0;
  copy_matrix(scatter, xs.scatter);
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
  copy_matrix(total, xs.total);
  copy_matrix(absorption, xs.absorption);
  copy_matrix(fission, xs.fission);
  copy_matrix(nu_fission, xs.nu_fission);
  copy_matrix(chi, xs.chi);
  copy_matrix(scatter, xs.scatter);
  return xs;
}

HD MgxsView view_xs(const Mgxs& xs)
{
  return {xs.total, xs.absorption, xs.fission, xs.nu_fission, xs.chi, xs.scatter};
}

// --- Geometry helpers (C5G7 assembly) --------------------------------------

constexpr double kPitch = 1.26;          // cm pin pitch
constexpr double kFuelRadius = 0.54;     // cm fuel pellet radius
constexpr double kAssemblyHalf = 17 * kPitch * 0.5; // half-width of 17x17 lattice

// 17x17 lattice map (row-major y,x) using material IDs defined above.
// Pattern follows the C5G7 benchmark MOX assembly with a central fission
// chamber (FC) and 24 guide tubes (GT).
const int kAssemblyMap[17 * 17] = {
  6,5,1,1,1,1,1,5,6,5,1,1,1,1,1,5,6,
  5,4,3,3,3,3,3,4,5,4,3,3,3,3,3,4,5,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  5,4,3,3,3,3,3,4,5,4,3,3,3,3,3,4,5,
  6,5,1,1,1,1,1,5,5,5,1,1,1,1,1,5,6,
  5,4,3,3,3,3,3,4,5,4,3,3,3,3,3,4,5,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  1,3,2,2,2,2,2,3,1,3,2,2,2,2,2,3,1,
  5,4,3,3,3,3,3,4,5,4,3,3,3,3,3,4,5,
  6,5,1,1,1,1,1,5,6,5,1,1,1,1,1,5,6
};

// Mapping for which pins contain fissionable fuel (UO2 or MOX) to bias the
// initial source distribution.
inline bool pin_is_fuel(int mat_id) { return mat_id == 1 || mat_id == 2 || mat_id == 3; }

HD inline bool in_bounds(const Vec3& r)
{
  return std::abs(r.x) < kAssemblyHalf && std::abs(r.y) < kAssemblyHalf && std::abs(r.z) < 200.0;
}

HD std::tuple<int, int, int> lattice_index(const Vec3& r)
{
  const double ix = (r.x + kAssemblyHalf) / kPitch;
  const double iy = (r.y + kAssemblyHalf) / kPitch;
  return {static_cast<int>(std::floor(ix)), static_cast<int>(std::floor(iy)), 0};
}

HD int assembly_material(const Vec3& r)
{
  auto [i, j, _] = lattice_index(r);
  if (i < 0 || i >= 17 || j < 0 || j >= 17) return -1;
  return kAssemblyMap[j * 17 + i];
}

HD bool inside_fuel_pin(const Vec3& r)
{
  auto [i, j, _] = lattice_index(r);
  const int mat = assembly_material(r);
  if (mat < 1 || mat > 3) return false;
  const double cx = -kAssemblyHalf + (i + 0.5) * kPitch;
  const double cy = -kAssemblyHalf + (j + 0.5) * kPitch;
  const double dx = r.x - cx;
  const double dy = r.y - cy;
  return dx * dx + dy * dy <= kFuelRadius * kFuelRadius;
}

HD int material_from_position(const Vec3& r)
{
  const int id = assembly_material(r);
  if (id < 0) return -1;
  if (id == 5 || id == 6) return id; // guide tube or moderator water region
  // Fuel pins split between fuel pellet and water moderator
  if (inside_fuel_pin(r)) return id;
  return 6; // outside pellet is moderator
}

struct SurfaceHit { double distance; int axis; };

HD SurfaceHit distance_to_surfaces(const Particle& p)
{
  const double bounds[3][2] {{-kAssemblyHalf, kAssemblyHalf}, {-kAssemblyHalf, kAssemblyHalf}, {-200.0, 200.0}};
  double min_d = std::numeric_limits<double>::infinity();
  int axis = -1;
  const double dir[3] {p.u.x, p.u.y, p.u.z};
  const double pos[3] {p.r.x, p.r.y, p.r.z};
  for (int a = 0; a < 3; ++a) {
    if (std::abs(dir[a]) < 1e-12) continue;
    for (int s = 0; s < 2; ++s) {
      const double d = (bounds[a][s] - pos[a]) / dir[a];
      if (d > 1e-10 && d < min_d) { min_d = d; axis = a; }
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

struct FuelHit { double distance; bool leaving_fuel; };

HD FuelHit distance_to_fuel_boundary(const Particle& p)
{
  auto [i, j, _] = lattice_index(p.r);
  const int mat = assembly_material(p.r);
  if (mat < 1 || mat > 3) return {std::numeric_limits<double>::infinity(), false};
  const double cx = -kAssemblyHalf + (i + 0.5) * kPitch;
  const double cy = -kAssemblyHalf + (j + 0.5) * kPitch;
  const double dx = p.r.x - cx;
  const double dy = p.r.y - cy;
  const double b = dx * p.u.x + dy * p.u.y;
  const double c = dx * dx + dy * dy - kFuelRadius * kFuelRadius;
  const double disc = b * b - c;
  if (disc <= 0.0) return {std::numeric_limits<double>::infinity(), false};
  const double sqrt_disc = std::sqrt(disc);
  const double d1 = -b - sqrt_disc;
  const double d2 = -b + sqrt_disc;
  double d_exit = std::numeric_limits<double>::infinity();
  if (d1 > 1e-10) d_exit = d1;
  else if (d2 > 1e-10) d_exit = d2;
  const bool leaving_fuel = inside_fuel_pin(p.r);
  return {d_exit, leaving_fuel};
}

// --- Source sampling --------------------------------------------------------

Particle sample_source(Rng& rng, const std::vector<int>& fuel_map)
{
  Particle p;
  // Choose a random fuel pin location uniformly.
  const size_t idx = static_cast<size_t>(rng.uniform() * fuel_map.size());
  const int pin = fuel_map[idx];
  const int px = pin % 17;
  const int py = pin / 17;
  const double cx = -kAssemblyHalf + (px + 0.5) * kPitch;
  const double cy = -kAssemblyHalf + (py + 0.5) * kPitch;

  // Rejection sample within pellet.
  while (true) {
    const double rx = (rng.uniform() - 0.5) * kPitch;
    const double ry = (rng.uniform() - 0.5) * kPitch;
    if (rx * rx + ry * ry <= kFuelRadius * kFuelRadius) {
      p.r.x = cx + rx;
      p.r.y = cy + ry;
      break;
    }
  }
  p.r.z = (rng.uniform() - 0.5) * 50.0; // small axial span within assembly
  p.u = sample_isotropic(rng);
  // Use the chi of UO2 for seeding; in practice fission sites overwrite this
  p.g = sample_energy_group(make_uo2().chi, rng);
  return p;
}

// --- Transport kernel ------------------------------------------------------

void transport_history(Particle& p, const Mgxs* materials, Tallies& tallies, Rng& rng,
  double& produced_neutrons, std::vector<Particle>& fission_bank)
{
  MgxsView xs_set[kMaterials];
  for (int i = 0; i < kMaterials; ++i) xs_set[i] = view_xs(materials[i]);

  while (true) {
    const int m = material_from_position(p.r);
    if (m < 0) return; // leaked (reflective box should prevent this)
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

      return; // capture
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

// --- CPU eigenvalue driver --------------------------------------------------

void run_eigenvalue_cpu(int batches, int inactive, int histories, uint64_t seed)
{
  Mgxs materials[kMaterials] {make_uo2(), make_mox43(), make_mox70(), make_mox87(), make_guide_tube(), make_fission_chamber(), make_water()};

  // Precompute fuel pin indices for sampling.
  std::vector<int> fuel_pins;
  fuel_pins.reserve(17 * 17);
  for (int j = 0; j < 17; ++j) {
    for (int i = 0; i < 17; ++i) {
      const int mat = kAssemblyMap[j * 17 + i];
      if (pin_is_fuel(mat)) fuel_pins.push_back(j * 17 + i);
    }
  }

  Tallies tallies(kMaterials, kGroups);

  std::vector<Particle> source_bank(histories);
  Rng seeder(seed);
  for (auto& p : source_bank) p = sample_source(seeder, fuel_pins);

  for (int b = 0; b < batches; ++b) {
    Tallies batch_tallies(kMaterials, kGroups);
    double produced = 0.0;
    std::vector<Particle> next_bank;
    next_bank.reserve(histories * 2);

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
  const char* labels[kMaterials] {"UO2", "MOX4.3", "MOX7.0", "MOX8.7", "Guide", "FC", "Water"};
  for (int m = 0; m < kMaterials; ++m) {
    std::cout << std::setw(7) << labels[m] << ": ";
    for (int g = 0; g < kGroups; ++g) {
      const double val = tallies.flux_tracklength[m * kGroups + g];
      std::cout << std::setw(10) << std::setprecision(4) << val;
    }
    std::cout << "\n";
  }
}

#ifdef __CUDACC__
// --- CUDA helpers -----------------------------------------------------------

inline void check_cuda(cudaError_t err, const char* msg)
{
  if (err != cudaSuccess) throw std::runtime_error(std::string(msg) + ": " + cudaGetErrorString(err));
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

D void atomic_accumulate(double* arr, int idx, double val) { atomicAdd(&arr[idx], val); }

D void append_fission(Particle* bank, int capacity, int* counter, const Particle& q, int copies)
{
  const int slot = atomicAdd(counter, copies);
  for (int n = 0; n < copies; ++n) {
    if (slot + n < capacity) bank[slot + n] = q;
  }
}

D double fast_log(double x) { return log(x); }

D Vec3 sample_isotropic_d(Rng& rng)
{
  const double mu = 2.0 * rng.uniform() - 1.0;
  const double phi = 2.0 * M_PI * rng.uniform();
  const double sin_t = sqrt(fmax(0.0, 1.0 - mu * mu));
  return {sin_t * cos(phi), sin_t * sin(phi), mu};
}

D int material_from_position_d(const Vec3& r)
{
  const double ix = (r.x + kAssemblyHalf) / kPitch;
  const double iy = (r.y + kAssemblyHalf) / kPitch;
  const int i = static_cast<int>(floor(ix));
  const int j = static_cast<int>(floor(iy));
  if (i < 0 || i >= 17 || j < 0 || j >= 17 || fabs(r.z) >= 200.0) return -1;
  const int mat = kAssemblyMap[j * 17 + i];
  if (mat == 5 || mat == 6) return mat;
  const double cx = -kAssemblyHalf + (i + 0.5) * kPitch;
  const double cy = -kAssemblyHalf + (j + 0.5) * kPitch;
  const double dx = r.x - cx;
  const double dy = r.y - cy;
  if (dx * dx + dy * dy <= kFuelRadius * kFuelRadius) return mat;
  return 6;
}

D SurfaceHit distance_to_surfaces_d(const Particle& p)
{
  const double bounds[3][2] {{-kAssemblyHalf, kAssemblyHalf}, {-kAssemblyHalf, kAssemblyHalf}, {-200.0, 200.0}};
  double min_d = 1e30;
  int axis = -1;
  const double dir[3] {p.u.x, p.u.y, p.u.z};
  const double pos[3] {p.r.x, p.r.y, p.r.z};
  for (int a = 0; a < 3; ++a) {
    if (fabs(dir[a]) < 1e-12) continue;
    for (int s = 0; s < 2; ++s) {
      const double d = (bounds[a][s] - pos[a]) / dir[a];
      if (d > 1e-10 && d < min_d) { min_d = d; axis = a; }
    }
  }
  return {min_d, axis};
}

D FuelHit distance_to_fuel_boundary_d(const Particle& p)
{
  const double ix = (p.r.x + kAssemblyHalf) / kPitch;
  const double iy = (p.r.y + kAssemblyHalf) / kPitch;
  const int i = static_cast<int>(floor(ix));
  const int j = static_cast<int>(floor(iy));
  if (i < 0 || i >= 17 || j < 0 || j >= 17) return {1e30, false};
  const int mat = kAssemblyMap[j * 17 + i];
  if (mat < 1 || mat > 3) return {1e30, false};
  const double cx = -kAssemblyHalf + (i + 0.5) * kPitch;
  const double cy = -kAssemblyHalf + (j + 0.5) * kPitch;
  const double dx = p.r.x - cx;
  const double dy = p.r.y - cy;
  const double b = dx * p.u.x + dy * p.u.y;
  const double c = dx * dx + dy * dy - kFuelRadius * kFuelRadius;
  const double disc = b * b - c;
  if (disc <= 0.0) return {1e30, false};
  const double sqrt_disc = sqrt(disc);
  const double d1 = -b - sqrt_disc;
  const double d2 = -b + sqrt_disc;
  double d_exit = 1e30;
  if (d1 > 1e-10) d_exit = d1; else if (d2 > 1e-10) d_exit = d2;
  const bool leaving = dx * dx + dy * dy <= kFuelRadius * kFuelRadius;
  return {d_exit, leaving};
}

D void reflect_d(Particle& p, int axis)
{
  if (axis == 0) p.u.x *= -1.0;
  if (axis == 1) p.u.y *= -1.0;
  if (axis == 2) p.u.z *= -1.0;
}

D void transport_kernel(Particle* source, int histories, MgxsView fuel, MgxsView mox43, MgxsView mox70,
  MgxsView mox87, MgxsView guide, MgxsView fc, MgxsView water, TalliesView tallies, Particle* bank,
  int max_bank, int* bank_counter, double* produced, uint64_t seed, int batch)
{
  MgxsView xs_set[kMaterials] {fuel, mox43, mox70, mox87, guide, fc, water};
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= histories) return;

  Particle p = source[tid];
  Rng rng(seed + static_cast<uint64_t>(batch) * histories + tid);

  while (true) {
    const int m = material_from_position_d(p.r);
    if (m < 0) return;
    const MgxsView xs = xs_set[m];
    const double sig_t = xs.total[p.g];
    const double free_path = -fast_log(rng.uniform()) / sig_t;
    SurfaceHit surf = distance_to_surfaces_d(p);
    FuelHit fuel_hit = distance_to_fuel_boundary_d(p);
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
        p.u = sample_isotropic_d(rng);
        continue;
      }

      if (xi < scatter_prob + fission_prob && xs.nu_fission[p.g] > 0.0) {
        const double nu = xs.nu_fission[p.g];
        const int num_new = static_cast<int>(floor(nu + rng.uniform()));
        atomic_accumulate(tallies.nu_fission, m * kGroups + p.g, nu * p.w);
        atomic_accumulate(produced, 0, nu * p.w);
        for (int n = 0; n < num_new; ++n) {
          Particle q;
          q.r = p.r;
          q.u = sample_isotropic_d(rng);
          q.g = sample_energy_group(xs.chi, rng);
          q.w = p.w;
          append_fission(bank, max_bank, bank_counter, q, 1);
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

    if (surf.distance < fuel_hit.distance) reflect_d(p, surf.axis);
  }
}

// --- CUDA host driver ------------------------------------------------------

void run_eigenvalue_cuda(int batches, int inactive, int histories, uint64_t seed)
{
  Mgxs fuel = make_uo2();
  Mgxs mox43 = make_mox43();
  Mgxs mox70 = make_mox70();
  Mgxs mox87 = make_mox87();
  Mgxs guide = make_guide_tube();
  Mgxs fc = make_fission_chamber();
  Mgxs water = make_water();

  DeviceXsBuffers fuel_dev {}, mox43_dev {}, mox70_dev {}, mox87_dev {}, guide_dev {}, fc_dev {}, water_dev {};
  MgxsView fuel_v {}, mox43_v {}, mox70_v {}, mox87_v {}, guide_v {}, fc_v {}, water_v {};
  allocate_device_xs(fuel, fuel_dev, fuel_v);
  allocate_device_xs(mox43, mox43_dev, mox43_v);
  allocate_device_xs(mox70, mox70_dev, mox70_v);
  allocate_device_xs(mox87, mox87_dev, mox87_v);
  allocate_device_xs(guide, guide_dev, guide_v);
  allocate_device_xs(fc, fc_dev, fc_v);
  allocate_device_xs(water, water_dev, water_v);

  const int tally_size = kMaterials * kGroups;
  DeviceBuffers buffers {};
  check_cuda(cudaMalloc(&buffers.source, sizeof(Particle) * histories), "alloc source");
  const int max_bank = histories * 8;
  check_cuda(cudaMalloc(&buffers.next_bank, sizeof(Particle) * max_bank), "alloc bank");
  check_cuda(cudaMalloc(&buffers.bank_count, sizeof(int)), "alloc bank count");
  check_cuda(cudaMalloc(&buffers.produced_sum, sizeof(double)), "alloc produced");
  check_cuda(cudaMalloc(&buffers.flux, sizeof(double) * tally_size), "alloc flux");
  check_cuda(cudaMalloc(&buffers.nu_fission, sizeof(double) * tally_size), "alloc nu fission");

  // CPU-side source bank.
  std::vector<int> fuel_pins;
  for (int j = 0; j < 17; ++j) {
    for (int i = 0; i < 17; ++i) if (pin_is_fuel(kAssemblyMap[j * 17 + i])) fuel_pins.push_back(j * 17 + i);
  }
  std::vector<Particle> source_bank(histories);
  Rng seeder(seed);
  for (auto& p : source_bank) p = sample_source(seeder, fuel_pins);

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
    transport_kernel<<<blocks, threads>>>(buffers.source, histories, fuel_v, mox43_v, mox70_v, mox87_v, guide_v, fc_v,
      water_v, tallies_v, buffers.next_bank, max_bank, buffers.bank_count, buffers.produced_sum, seed, b);
    check_cuda(cudaDeviceSynchronize(), "run kernel");

    int bank_size = 0;
    check_cuda(cudaMemcpy(&bank_size, buffers.bank_count, sizeof(int), cudaMemcpyDeviceToHost), "copy bank size");
    std::vector<Particle> next_bank(bank_size);
    if (bank_size > 0) check_cuda(cudaMemcpy(next_bank.data(), buffers.next_bank, sizeof(Particle) * bank_size, cudaMemcpyDeviceToHost), "copy bank");
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
  const char* labels[kMaterials] {"UO2", "MOX4.3", "MOX7.0", "MOX8.7", "Guide", "FC", "Water"};
  for (int m = 0; m < kMaterials; ++m) {
    std::cout << std::setw(7) << labels[m] << ": ";
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
  cudaFree(mox43_dev.total);
  cudaFree(mox43_dev.absorption);
  cudaFree(mox43_dev.fission);
  cudaFree(mox43_dev.nu_fission);
  cudaFree(mox43_dev.chi);
  cudaFree(mox43_dev.scatter);
  cudaFree(mox70_dev.total);
  cudaFree(mox70_dev.absorption);
  cudaFree(mox70_dev.fission);
  cudaFree(mox70_dev.nu_fission);
  cudaFree(mox70_dev.chi);
  cudaFree(mox70_dev.scatter);
  cudaFree(mox87_dev.total);
  cudaFree(mox87_dev.absorption);
  cudaFree(mox87_dev.fission);
  cudaFree(mox87_dev.nu_fission);
  cudaFree(mox87_dev.chi);
  cudaFree(mox87_dev.scatter);
  cudaFree(guide_dev.total);
  cudaFree(guide_dev.absorption);
  cudaFree(guide_dev.fission);
  cudaFree(guide_dev.nu_fission);
  cudaFree(guide_dev.chi);
  cudaFree(guide_dev.scatter);
  cudaFree(fc_dev.total);
  cudaFree(fc_dev.absorption);
  cudaFree(fc_dev.fission);
  cudaFree(fc_dev.nu_fission);
  cudaFree(fc_dev.chi);
  cudaFree(fc_dev.scatter);
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
    const uint64_t seed = 20240520;
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
