#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace simple_transport {

// ============================= 向量与工具 ====================================
// Vec3 为三维向量工具类，对应 OpenMC 中 Particle 的位置和方向向量。
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  Vec3() = default;
  Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

  Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  Vec3 operator-(const Vec3& other) const { return {x - other.x, y - other.y, z - other.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }

  Vec3& operator+=(const Vec3& other) {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalise(const Vec3& v) {
  double n = norm(v);
  if (n == 0.0) return {0.0, 0.0, 1.0};
  return {v.x / n, v.y / n, v.z / n};
}

constexpr double kPi = 3.14159265358979323846;

// ============================= 材料系统 =====================================
// Material 与 OpenMC 的 Material 类类似，封装多群宏观截面与裂变谱。
class Material {
 public:
  Material(std::string name, std::vector<double> sigma_a, std::vector<double> sigma_f,
           std::vector<double> nu, std::vector<std::vector<double>> sigma_s,
           std::vector<double> chi)
      : name_(std::move(name)), sigma_a_(std::move(sigma_a)),
        sigma_f_(std::move(sigma_f)), nu_(std::move(nu)),
        sigma_s_(std::move(sigma_s)), chi_(std::move(chi)) {
    std::size_t n = sigma_a_.size();
    if (sigma_f_.size() != n || nu_.size() != n || sigma_s_.size() != n || chi_.size() != n) {
      throw std::runtime_error("材料截面数组长度不一致");
    }
  }

  [[nodiscard]] std::size_t groups() const noexcept { return sigma_a_.size(); }
  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] double sigma_a(std::size_t g) const { return sigma_a_.at(g); }
  [[nodiscard]] double sigma_f(std::size_t g) const { return sigma_f_.at(g); }
  [[nodiscard]] double nu(std::size_t g) const { return nu_.at(g); }
  [[nodiscard]] const std::vector<double>& chi() const noexcept { return chi_; }

  [[nodiscard]] double sigma_t(std::size_t g) const {
    double scatter = std::accumulate(sigma_s_[g].begin(), sigma_s_[g].end(), 0.0);
    return sigma_a_[g] + sigma_f_[g] + scatter;
  }

  [[nodiscard]] const std::vector<double>& scatter_row(std::size_t g) const {
    return sigma_s_.at(g);
  }

 private:
  std::string name_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> sigma_s_;
  std::vector<double> chi_;
};

// ============================= 几何描述 =====================================
// 1. 边界条件枚举，完全对应 OpenMC 中的 boundary 类型。
enum class Boundary { Interface, Vacuum, Reflective };

// 2. 轴向枚举，用于描述平面法向量方向。
enum class Axis { X = 0, Y = 1, Z = 2 };

// 3. Surface 表示与坐标轴对齐的平面（x/y/z = 常数）。
class Surface {
 public:
  Surface(std::string name, Axis axis, double coord, Boundary boundary)
      : name_(std::move(name)), axis_(axis), coordinate_(coord), boundary_(boundary) {}

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] Axis axis() const noexcept { return axis_; }
  [[nodiscard]] double coordinate() const noexcept { return coordinate_; }
  [[nodiscard]] Boundary boundary() const noexcept { return boundary_; }

 private:
  std::string name_;
  Axis axis_;
  double coordinate_;
  Boundary boundary_;
};

// 4. Cell 保存材料、六个面以及邻接关系，模拟三维笛卡尔网格单元。
class Cell {
 public:
  Cell(std::string name, int id, const Surface* x_minus, const Surface* x_plus,
       const Surface* y_minus, const Surface* y_plus, const Surface* z_minus,
       const Surface* z_plus, const Material* material)
      : name_(std::move(name)), id_(id), material_(material) {
    surfaces_[static_cast<int>(Axis::X)][0] = x_minus;
    surfaces_[static_cast<int>(Axis::X)][1] = x_plus;
    surfaces_[static_cast<int>(Axis::Y)][0] = y_minus;
    surfaces_[static_cast<int>(Axis::Y)][1] = y_plus;
    surfaces_[static_cast<int>(Axis::Z)][0] = z_minus;
    surfaces_[static_cast<int>(Axis::Z)][1] = z_plus;
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] int id() const noexcept { return id_; }
  [[nodiscard]] const Material& material() const { return *material_; }

  [[nodiscard]] bool contains(const Vec3& r) const noexcept {
    return r.x >= surfaces_[0][0]->coordinate() && r.x < surfaces_[0][1]->coordinate() &&
           r.y >= surfaces_[1][0]->coordinate() && r.y < surfaces_[1][1]->coordinate() &&
           r.z >= surfaces_[2][0]->coordinate() && r.z < surfaces_[2][1]->coordinate();
  }

  void set_neighbour(Axis axis, bool positive, const Cell* neighbour) {
    neighbours_[static_cast<int>(axis)][positive ? 1 : 0] = neighbour;
  }

  [[nodiscard]] const Cell* neighbour(Axis axis, bool positive) const noexcept {
    return neighbours_[static_cast<int>(axis)][positive ? 1 : 0];
  }

  [[nodiscard]] const Surface& surface(Axis axis, bool positive) const noexcept {
    return *surfaces_[static_cast<int>(axis)][positive ? 1 : 0];
  }

  struct BoundaryHit {
    double distance = std::numeric_limits<double>::infinity();
    const Surface* surface = nullptr;
    Axis axis = Axis::X;
    bool positive = true;
  };

  [[nodiscard]] BoundaryHit distance_to_boundary(const Vec3& position,
                                                 const Vec3& direction) const {
    BoundaryHit hit;
    for (int ax = 0; ax < 3; ++ax) {
      double component = (ax == 0) ? direction.x : (ax == 1 ? direction.y : direction.z);
      if (std::abs(component) < 1e-14) continue;
      bool positive = component > 0.0;
      const Surface* s = surfaces_[ax][positive ? 1 : 0];
      double coord = s->coordinate();
      double pos = (ax == 0) ? position.x : (ax == 1 ? position.y : position.z);
      double distance = (coord - pos) / component;
      if (distance <= 1e-14) continue;
      if (distance < hit.distance) {
        hit.distance = distance;
        hit.surface = s;
        hit.axis = static_cast<Axis>(ax);
        hit.positive = positive;
      }
    }
    return hit;
  }

 private:
  std::string name_;
  int id_;
  const Material* material_;
  const Surface* surfaces_[3][2] = {};
  const Cell* neighbours_[3][2] = {};
};

// 5. Universe 与 Lattice 的封装沿袭 OpenMC 的层次组织思想。
class Universe {
 public:
  explicit Universe(std::string name) : name_(std::move(name)) {}

  Cell& add_cell(Cell cell) {
    cells_.push_back(std::move(cell));
    return cells_.back();
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }

  [[nodiscard]] const std::vector<Cell>& cells() const noexcept { return cells_; }

  [[nodiscard]] Cell& cell(std::size_t index) { return cells_.at(index); }

  [[nodiscard]] const Cell* find_cell(const Vec3& r) const noexcept {
    for (const auto& cell : cells_) {
      if (cell.contains(r)) return &cell;
    }
    return nullptr;
  }

 private:
  std::string name_;
  std::vector<Cell> cells_;
};

class Lattice3D {
 public:
  Lattice3D(std::string name, Vec3 pitch, std::array<std::size_t, 3> dims, Vec3 origin)
      : name_(std::move(name)), pitch_(pitch), dims_(dims), origin_(origin) {
    universes_.resize(dims_[0] * dims_[1] * dims_[2], nullptr);
  }

  void set_universe(std::size_t i, std::size_t j, std::size_t k, const Universe* universe) {
    universes_.at(index(i, j, k)) = universe;
  }

  [[nodiscard]] const Universe* universe_at(const Vec3& r) const noexcept {
    double lx = r.x - origin_.x;
    double ly = r.y - origin_.y;
    double lz = r.z - origin_.z;
    if (lx < 0.0 || ly < 0.0 || lz < 0.0) return nullptr;
    std::size_t i = static_cast<std::size_t>(lx / pitch_.x);
    std::size_t j = static_cast<std::size_t>(ly / pitch_.y);
    std::size_t k = static_cast<std::size_t>(lz / pitch_.z);
    if (i >= dims_[0] || j >= dims_[1] || k >= dims_[2]) return nullptr;
    return universes_[index(i, j, k)];
  }

  [[nodiscard]] const std::vector<const Universe*>& universes() const noexcept {
    return universes_;
  }

 private:
  [[nodiscard]] std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return (k * dims_[1] + j) * dims_[0] + i;
  }

  std::string name_;
  Vec3 pitch_;
  std::array<std::size_t, 3> dims_;
  Vec3 origin_;
  std::vector<const Universe*> universes_;
};

// 6. Geometry 负责根宇宙与晶格组合，提供粒子定位。
class Geometry {
 public:
  Geometry(const Universe& root, const Lattice3D& lattice) : root_(root), lattice_(lattice) {}

  [[nodiscard]] const Cell* locate(const Vec3& r) const noexcept {
    if (const Universe* u = lattice_.universe_at(r)) {
      if (const Cell* c = u->find_cell(r)) return c;
    }
    return root_.find_cell(r);
  }

  [[nodiscard]] std::vector<const Cell*> enumerate_cells() const {
    std::vector<const Cell*> result;
    for (const auto& cell : root_.cells()) result.push_back(&cell);
    for (const Universe* u : lattice_.universes()) {
      if (!u) continue;
      for (const auto& cell : u->cells()) result.push_back(&cell);
    }
    return result;
  }

 private:
  const Universe& root_;
  const Lattice3D& lattice_;
};

// ============================= 粒子与统计 ====================================
class Particle {
 public:
  void set_position(Vec3 r) { position_ = r; }
  void set_direction(Vec3 u) { direction_ = normalise(u); }
  void set_weight(double w) { weight_ = w; }
  void set_group(std::size_t g) { group_ = g; }
  void set_cell(const Cell* cell) { cell_ = cell; }
  void set_alive(bool alive) { alive_ = alive; }

  [[nodiscard]] const Vec3& position() const noexcept { return position_; }
  [[nodiscard]] const Vec3& direction() const noexcept { return direction_; }
  [[nodiscard]] double weight() const noexcept { return weight_; }
  [[nodiscard]] std::size_t group() const noexcept { return group_; }
  [[nodiscard]] bool alive() const noexcept { return alive_; }
  [[nodiscard]] const Cell* cell() const noexcept { return cell_; }

  void move(double distance) { position_ += direction_ * distance; }

 private:
  Vec3 position_;
  Vec3 direction_ = {0.0, 0.0, 1.0};
  double weight_ = 1.0;
  std::size_t group_ = 0;
  bool alive_ = true;
  const Cell* cell_ = nullptr;
};

class Bank {
 public:
  void clear() { particles_.clear(); }
  void add(const Particle& p) { particles_.push_back(p); }
  [[nodiscard]] std::size_t size() const noexcept { return particles_.size(); }

  [[nodiscard]] Particle& operator[](std::size_t i) { return particles_.at(i); }
  [[nodiscard]] const Particle& operator[](std::size_t i) const { return particles_.at(i); }

  [[nodiscard]] std::vector<Particle>::iterator begin() { return particles_.begin(); }
  [[nodiscard]] std::vector<Particle>::iterator end() { return particles_.end(); }

 private:
  std::vector<Particle> particles_;
};

class TallyManager {
 public:
  TallyManager(std::size_t cell_count, std::size_t group_count)
      : cell_pathlength_(cell_count, std::vector<double>(group_count, 0.0)) {}

  void score_flux(const Particle& p, double track_length) {
    if (!p.cell()) return;
    auto& row = cell_pathlength_.at(p.cell()->id());
    row.at(p.group()) += track_length * p.weight();
  }

  void score_leakage(const Particle& p) { leakage_ += p.weight(); }

  [[nodiscard]] const std::vector<std::vector<double>>& cell_pathlength() const {
    return cell_pathlength_;
  }

  [[nodiscard]] double leakage() const noexcept { return leakage_; }

 private:
  std::vector<std::vector<double>> cell_pathlength_;
  double leakage_ = 0.0;
};

// ============================= 随机与物理 ====================================
class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed) : generator_(seed), uniform_(0.0, 1.0) {}

  double sample_uniform() { return uniform_(generator_); }

  Vec3 sample_isotropic_direction() {
    double mu = 2.0 * sample_uniform() - 1.0;
    double phi = 2.0 * kPi * sample_uniform();
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
  }

  std::mt19937& generator() noexcept { return generator_; }

 private:
  std::mt19937 generator_;
  std::uniform_real_distribution<double> uniform_;
};

class PhysicsDriver {
 public:
  PhysicsDriver(const Geometry& geometry, unsigned seed)
      : geometry_(geometry), rng_(seed) {}

  void initialise_source(Bank& source, std::size_t particle_count, std::size_t group_count,
                         const Vec3& source_position) {
    source.clear();
    for (std::size_t n = 0; n < particle_count; ++n) {
      Particle p;
      p.set_position(source_position);
      p.set_direction(rng_.sample_isotropic_direction());
      p.set_group(n % group_count);
      p.set_weight(1.0);
      p.set_alive(true);
      const Cell* cell = geometry_.locate(p.position());
      if (!cell) throw std::runtime_error("源初始化位置超出几何范围");
      p.set_cell(cell);
      source.add(p);
    }
  }

  double transport_generation(Bank& source, Bank& fission_bank, TallyManager& tallies) {
    fission_bank.clear();
    double produced = 0.0;
    for (auto& particle : source) {
      transport_particle(particle, fission_bank, tallies, produced);
    }
    if (source.size() == 0) throw std::runtime_error("源粒子数量为零");
    return produced / static_cast<double>(source.size());
  }

  void resample_source(const Bank& fission_bank, Bank& source) {
    if (fission_bank.size() == 0) {
      throw std::runtime_error("裂变粒子库为空，系统处于次临界状态");
    }
    std::uniform_int_distribution<std::size_t> pick(0, fission_bank.size() - 1);
    for (auto& particle : source) {
      particle = fission_bank[pick(rng_.generator())];
      particle.set_weight(1.0);
      particle.set_alive(true);
      particle.set_cell(geometry_.locate(particle.position()));
    }
  }

 private:
  const Geometry& geometry_;
  RandomEngine rng_;

  void transport_particle(Particle& p, Bank& fission_bank, TallyManager& tallies,
                          double& produced) {
    while (p.alive() && p.cell() != nullptr) {
      const Material& mat = p.cell()->material();
      double sigma_t = mat.sigma_t(p.group());
      double xi = rng_.sample_uniform();
      double free_path = -std::log(1.0 - xi) / sigma_t;
      auto hit = p.cell()->distance_to_boundary(p.position(), p.direction());
      double distance = std::min(free_path, hit.distance);

      tallies.score_flux(p, distance);
      p.move(distance);

      if (hit.surface && hit.distance <= free_path) {
        handle_boundary(p, hit, tallies);
        continue;
      }

      collide(p, mat, fission_bank, produced);
    }
  }

  void handle_boundary(Particle& p, const Cell::BoundaryHit& hit, TallyManager& tallies) {
    if (!hit.surface) {
      tallies.score_leakage(p);
      p.set_alive(false);
      return;
    }
    const Surface& surface = *hit.surface;
    switch (surface.boundary()) {
      case Boundary::Vacuum:
        tallies.score_leakage(p);
        p.set_alive(false);
        break;
      case Boundary::Reflective: {
        Vec3 dir = p.direction();
        Vec3 pos = p.position();
        double eps = 1e-9;
        if (hit.axis == Axis::X) {
          dir.x = -dir.x;
          pos.x += hit.positive ? -eps : eps;
        } else if (hit.axis == Axis::Y) {
          dir.y = -dir.y;
          pos.y += hit.positive ? -eps : eps;
        } else {
          dir.z = -dir.z;
          pos.z += hit.positive ? -eps : eps;
        }
        p.set_direction(dir);
        p.set_position(pos);
        break;
      }
      case Boundary::Interface: {
        const Cell* next = p.cell()->neighbour(hit.axis, hit.positive);
        if (!next) {
          tallies.score_leakage(p);
          p.set_alive(false);
        } else {
          Vec3 pos = p.position();
          double eps = 1e-9;
          if (hit.axis == Axis::X) {
            pos.x += hit.positive ? eps : -eps;
          } else if (hit.axis == Axis::Y) {
            pos.y += hit.positive ? eps : -eps;
          } else {
            pos.z += hit.positive ? eps : -eps;
          }
          p.set_position(pos);
          p.set_cell(next);
        }
        break;
      }
    }
    if (p.alive()) {
      p.set_cell(geometry_.locate(p.position()));
    }
  }

  void collide(Particle& p, const Material& mat, Bank& fission_bank, double& produced) {
    double sigma_a = mat.sigma_a(p.group());
    double sigma_f = mat.sigma_f(p.group());
    const auto& row = mat.scatter_row(p.group());
    double scatter_sum = std::accumulate(row.begin(), row.end(), 0.0);
    double sigma_t = sigma_a + sigma_f + scatter_sum;
    double xi = rng_.sample_uniform();

    double cumulative = sigma_a / sigma_t;
    if (xi < cumulative) {
      p.set_alive(false);
      return;
    }

    xi -= cumulative;
    cumulative = scatter_sum / sigma_t;
    if (xi < cumulative) {
      scatter(p, row);
      return;
    }

    fission(p, mat, fission_bank, produced);
  }

  void scatter(Particle& p, const std::vector<double>& row) {
    double total = std::accumulate(row.begin(), row.end(), 0.0);
    double xi = rng_.sample_uniform() * total;
    double accum = 0.0;
    for (std::size_t g = 0; g < row.size(); ++g) {
      accum += row[g];
      if (xi < accum) {
        p.set_group(g);
        break;
      }
    }
    p.set_direction(rng_.sample_isotropic_direction());
  }

  void fission(const Particle& p, const Material& mat, Bank& fission_bank, double& produced) {
    double nu = mat.nu(p.group());
    int secondaries = static_cast<int>(std::floor(nu));
    if (rng_.sample_uniform() < nu - secondaries) ++secondaries;
    produced += static_cast<double>(secondaries);
    if (secondaries == 0) return;

    std::discrete_distribution<std::size_t> spectrum(mat.chi().begin(), mat.chi().end());
    for (int i = 0; i < secondaries; ++i) {
      Particle q;
      Vec3 dir = rng_.sample_isotropic_direction();
      q.set_direction(dir);
      q.set_position(p.position() + dir * 1e-7);
      q.set_group(spectrum(rng_.generator()));
      q.set_weight(1.0);
      q.set_alive(true);
      const Cell* cell = geometry_.locate(q.position());
      q.set_cell(cell);
      fission_bank.add(q);
    }
  }
};

// ============================= 模拟封装 =====================================
struct SimulationSettings {
  std::size_t groups;
  std::size_t histories;
  std::size_t generations;
  Vec3 source_position;
  unsigned seed;
};

class Simulation {
 public:
  Simulation(const Geometry& geometry, SimulationSettings settings)
      : geometry_(geometry), settings_(settings),
        tallies_(compute_cell_capacity(geometry), settings.groups),
        physics_(geometry_, settings.seed) {
    physics_.initialise_source(source_, settings.histories, settings.groups,
                               settings.source_position);
  }

  void run() {
    std::cout << std::fixed << std::setprecision(5);
    for (std::size_t g = 0; g < settings_.generations; ++g) {
      double keff = physics_.transport_generation(source_, fission_bank_, tallies_);
      physics_.resample_source(fission_bank_, source_);
      std::cout << "第 " << g + 1 << " 代有效增殖因子估计 k-effective = " << keff << '\n';
    }

    std::cout << "\n按单元与能群划分的路径长度计数 (cm)：\n";
    const auto cells = geometry_.enumerate_cells();
    for (const Cell* cell : cells) {
      if (!cell) continue;
      std::cout << "  单元 " << cell->name() << "：";
      for (std::size_t g = 0; g < settings_.groups; ++g) {
        std::cout << " 群" << g << '='
                  << tallies_.cell_pathlength().at(cell->id()).at(g) /
                         static_cast<double>(settings_.generations);
      }
      std::cout << '\n';
    }

    std::cout << "\n总泄漏权重："
              << tallies_.leakage() / static_cast<double>(settings_.generations) << '\n';
  }

 private:
  static std::size_t compute_cell_capacity(const Geometry& geometry) {
    std::size_t max_id = 0;
    for (const Cell* cell : geometry.enumerate_cells()) {
      if (cell && static_cast<std::size_t>(cell->id()) > max_id) {
        max_id = static_cast<std::size_t>(cell->id());
      }
    }
    return max_id + 1;
  }

  const Geometry& geometry_;
  SimulationSettings settings_;
  TallyManager tallies_;
  PhysicsDriver physics_;
  Bank source_;
  Bank fission_bank_;
};

}  // namespace simple_transport

// ============================= 驱动程序 =====================================
struct DemoModel {
  // 材料定义
  simple_transport::Material fuel_fresh;
  simple_transport::Material fuel_burned;
  simple_transport::Material moderator;
  simple_transport::Material reflector;

  // 平面定义：x/y/z 方向的最小、最大表面及内部界面。
  simple_transport::Surface x0;
  simple_transport::Surface x1;
  simple_transport::Surface x2;
  simple_transport::Surface y0;
  simple_transport::Surface y1;
  simple_transport::Surface y2;
  simple_transport::Surface z0;
  simple_transport::Surface z1;
  simple_transport::Surface z2;

  // 根宇宙与晶格
  simple_transport::Universe root_universe;
  simple_transport::Universe u000;
  simple_transport::Universe u100;
  simple_transport::Universe u010;
  simple_transport::Universe u110;
  simple_transport::Universe u001;
  simple_transport::Universe u101;
  simple_transport::Universe u011;
  simple_transport::Universe u111;
  simple_transport::Lattice3D lattice;
  simple_transport::Geometry geometry;

  DemoModel()
      : fuel_fresh("fuel_fresh", {0.010, 0.07, 0.18}, {0.022, 0.10, 0.14}, {2.6, 2.4, 2.1},
                   {{0.32, 0.04, 0.01}, {0.02, 0.38, 0.05}, {0.01, 0.09, 0.28}},
                   {0.75, 0.20, 0.05}),
        fuel_burned("fuel_burned", {0.013, 0.085, 0.22}, {0.018, 0.11, 0.12}, {2.3, 2.1, 2.0},
                    {{0.28, 0.05, 0.02}, {0.03, 0.36, 0.06}, {0.02, 0.08, 0.30}},
                    {0.65, 0.25, 0.10}),
        moderator("moderator", {0.001, 0.009, 0.05}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                  {{0.22, 0.05, 0.01}, {0.02, 0.34, 0.09}, {0.01, 0.06, 0.38}},
                  {0.92, 0.08, 0.0}),
        reflector("reflector", {0.0004, 0.004, 0.035}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                  {{0.27, 0.02, 0.01}, {0.01, 0.28, 0.10}, {0.01, 0.04, 0.38}}, {1.0, 0.0, 0.0}),
        x0("x0", simple_transport::Axis::X, 0.0, simple_transport::Boundary::Vacuum),
        x1("x1", simple_transport::Axis::X, 2.0, simple_transport::Boundary::Interface),
        x2("x2", simple_transport::Axis::X, 4.0, simple_transport::Boundary::Vacuum),
        y0("y0", simple_transport::Axis::Y, 0.0, simple_transport::Boundary::Vacuum),
        y1("y1", simple_transport::Axis::Y, 2.0, simple_transport::Boundary::Interface),
        y2("y2", simple_transport::Axis::Y, 4.0, simple_transport::Boundary::Vacuum),
        z0("z0", simple_transport::Axis::Z, 0.0, simple_transport::Boundary::Reflective),
        z1("z1", simple_transport::Axis::Z, 3.0, simple_transport::Boundary::Interface),
        z2("z2", simple_transport::Axis::Z, 6.0, simple_transport::Boundary::Vacuum),
        root_universe("root_universe"), u000("u000"), u100("u100"), u010("u010"),
        u110("u110"), u001("u001"), u101("u101"), u011("u011"), u111("u111"),
        lattice("assembly", {2.0, 2.0, 3.0}, {2, 2, 2}, {0.0, 0.0, 0.0}),
        geometry(root_universe, lattice) {
    // 创建每个晶格位置的单元。
    u000.add_cell({"fuel_lower_fresh", 0, &x0, &x1, &y0, &y1, &z0, &z1, &fuel_fresh});
    u100.add_cell({"fuel_lower_burned", 1, &x1, &x2, &y0, &y1, &z0, &z1, &fuel_burned});
    u010.add_cell({"reflector_lower_front", 2, &x0, &x1, &y1, &y2, &z0, &z1, &reflector});
    u110.add_cell({"reflector_lower_back", 3, &x1, &x2, &y1, &y2, &z0, &z1, &reflector});
    u001.add_cell({"moderator_upper_left", 4, &x0, &x1, &y0, &y1, &z1, &z2, &moderator});
    u101.add_cell({"moderator_upper_right", 5, &x1, &x2, &y0, &y1, &z1, &z2, &moderator});
    u011.add_cell({"reflector_upper_front", 6, &x0, &x1, &y1, &y2, &z1, &z2, &reflector});
    u111.add_cell({"reflector_upper_back", 7, &x1, &x2, &y1, &y2, &z1, &z2, &reflector});

    // 绑定晶格映射。
    lattice.set_universe(0, 0, 0, &u000);
    lattice.set_universe(1, 0, 0, &u100);
    lattice.set_universe(0, 1, 0, &u010);
    lattice.set_universe(1, 1, 0, &u110);
    lattice.set_universe(0, 0, 1, &u001);
    lattice.set_universe(1, 0, 1, &u101);
    lattice.set_universe(0, 1, 1, &u011);
    lattice.set_universe(1, 1, 1, &u111);

    // 将单元指针按照 (i, j, k) 存储，便于批量建立邻接关系。
    simple_transport::Cell* cells[2][2][2];
    cells[0][0][0] = &u000.cell(0);
    cells[1][0][0] = &u100.cell(0);
    cells[0][1][0] = &u010.cell(0);
    cells[1][1][0] = &u110.cell(0);
    cells[0][0][1] = &u001.cell(0);
    cells[1][0][1] = &u101.cell(0);
    cells[0][1][1] = &u011.cell(0);
    cells[1][1][1] = &u111.cell(0);

    for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
        for (std::size_t k = 0; k < 2; ++k) {
          simple_transport::Cell* cell = cells[i][j][k];
          if (i + 1 < 2) cell->set_neighbour(simple_transport::Axis::X, true, cells[i + 1][j][k]);
          if (i > 0) cell->set_neighbour(simple_transport::Axis::X, false, cells[i - 1][j][k]);
          if (j + 1 < 2) cell->set_neighbour(simple_transport::Axis::Y, true, cells[i][j + 1][k]);
          if (j > 0) cell->set_neighbour(simple_transport::Axis::Y, false, cells[i][j - 1][k]);
          if (k + 1 < 2) cell->set_neighbour(simple_transport::Axis::Z, true, cells[i][j][k + 1]);
          if (k > 0) cell->set_neighbour(simple_transport::Axis::Z, false, cells[i][j][k - 1]);
        }
      }
    }
  }
};

int main() {
  DemoModel model;
  simple_transport::SimulationSettings settings{3, 10000, 15, {1.2, 0.8, 1.5}, 202405u};
  simple_transport::Simulation simulation(model.geometry, settings);
  simulation.run();
}
