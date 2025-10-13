#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// 本示例整合材料、几何、粒子、随机数等对象，重现 OpenMC 中最核心的物理驱动循环。
// 相比完整示例，这里只追踪单一立方体内的粒子，便于逐行理解：
//   - 如何抽样自由程；
//   - 如何处理边界条件；
//   - 如何根据截面执行碰撞分支。

namespace lesson_physics {

// ============================= 工具与数据结构 =============================
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;

  Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalise(const Vec3& v) {
  double n = norm(v);
  if (n == 0.0) return {0.0, 0.0, 1.0};
  return {v.x / n, v.y / n, v.z / n};
}

class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed = 123u) : rng_(seed) {}

  double uniform() { return uniform_(rng_); }

  Vec3 isotropic_direction() {
    double mu = 2.0 * uniform() - 1.0;  // cos(theta)
    double phi = 2.0 * kPi * uniform();
    double sin_theta = std::sqrt(1.0 - mu * mu);
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
  }

  std::size_t pick(const std::vector<double>& pdf) {
    std::discrete_distribution<std::size_t> dist(pdf.begin(), pdf.end());
    return dist(rng_);
  }

 private:
  static constexpr double kPi = 3.14159265358979323846;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

class Material {
 public:
  Material(std::string name, std::vector<double> sigma_a, std::vector<double> sigma_f,
           std::vector<double> nu, std::vector<std::vector<double>> sigma_s,
           std::vector<double> chi)
      : name_(std::move(name)), sigma_a_(std::move(sigma_a)), sigma_f_(std::move(sigma_f)),
        nu_(std::move(nu)), sigma_s_(std::move(sigma_s)), chi_(std::move(chi)) {
    std::size_t n = sigma_a_.size();
    if (sigma_f_.size() != n || nu_.size() != n || sigma_s_.size() != n || chi_.size() != n) {
      throw std::runtime_error("材料截面维度不一致");
    }
    for (const auto& row : sigma_s_) {
      if (row.size() != n) throw std::runtime_error("散射矩阵不是方阵");
    }
  }

  [[nodiscard]] double sigma_t(std::size_t g) const {
    double scatter = std::accumulate(sigma_s_[g].begin(), sigma_s_[g].end(), 0.0);
    return sigma_a_[g] + sigma_f_[g] + scatter;
  }

  [[nodiscard]] double sigma_a(std::size_t g) const { return sigma_a_.at(g); }
  [[nodiscard]] double sigma_f(std::size_t g) const { return sigma_f_.at(g); }
  [[nodiscard]] double nu(std::size_t g) const { return nu_.at(g); }
  [[nodiscard]] const std::vector<double>& chi() const noexcept { return chi_; }
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

enum class Boundary { Interface, Vacuum, Reflective };
enum class Axis { X = 0, Y = 1, Z = 2 };

class Surface {
 public:
  Surface(std::string name, Axis axis, double coord, Boundary boundary)
      : name_(std::move(name)), axis_(axis), coordinate_(coord), boundary_(boundary) {}

  [[nodiscard]] double coordinate() const noexcept { return coordinate_; }
  [[nodiscard]] Boundary boundary() const noexcept { return boundary_; }
  [[nodiscard]] const std::string& name() const noexcept { return name_; }

 private:
  std::string name_;
  Axis axis_;
  double coordinate_;
  Boundary boundary_;
};

class Cell {
 public:
  Cell(std::string name, const Surface* x_minus, const Surface* x_plus, const Surface* y_minus,
       const Surface* y_plus, const Surface* z_minus, const Surface* z_plus, const Material* mat)
      : name_(std::move(name)), material_(mat) {
    surfaces_[0][0] = x_minus;
    surfaces_[0][1] = x_plus;
    surfaces_[1][0] = y_minus;
    surfaces_[1][1] = y_plus;
    surfaces_[2][0] = z_minus;
    surfaces_[2][1] = z_plus;
  }

  struct BoundaryHit {
    double distance = std::numeric_limits<double>::infinity();
    const Surface* surface = nullptr;
    Axis axis = Axis::X;
    bool positive = true;
  };

  BoundaryHit distance_to_boundary(const Vec3& position, const Vec3& direction) const {
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

  [[nodiscard]] const Material& material() const { return *material_; }
  [[nodiscard]] const std::string& name() const noexcept { return name_; }

 private:
  std::string name_;
  const Material* material_;
  const Surface* surfaces_[3][2] = {};
};

struct Particle {
  Vec3 position;
  Vec3 direction;
  std::size_t group = 0;
  double weight = 1.0;
  bool alive = true;
};

// ============================= 物理驱动 =============================
class PhysicsDriver {
 public:
  PhysicsDriver(const Material& material, const Cell& cell)
      : material_(material), cell_(cell) {}

  void transport(Particle& p, RandomEngine& rng) {
    std::cout << "--- 追踪开始（入射群=" << p.group << ") ---\n";
    for (int event = 0; event < 20 && p.alive; ++event) {
      double sigma_t = material_.sigma_t(p.group);
      double free_path = -std::log(1.0 - rng.uniform()) / sigma_t;
      Cell::BoundaryHit hit = cell_.distance_to_boundary(p.position, p.direction);

      if (hit.surface && hit.distance < free_path) {
        // 先走到边界。根据边界条件决定后续动作。
        p.position = p.position + p.direction * (hit.distance + 1e-9);
        handle_boundary(p, hit, rng);
      } else {
        // 在体内发生碰撞，执行反应抽样。
        p.position = p.position + p.direction * free_path;
        collide(p, rng);
      }
    }
    if (!p.alive) {
      std::cout << "粒子终止，当前位置 = (" << p.position.x << ", " << p.position.y << ", "
                << p.position.z << ")\n";
    } else {
      std::cout << "达到事件上限，停止追踪。\n";
    }
  }

 private:
  void handle_boundary(Particle& p, const Cell::BoundaryHit& hit, RandomEngine& rng) {
    switch (hit.surface->boundary()) {
      case Boundary::Vacuum:
        std::cout << "  击中真空边界 " << hit.surface->name() << "，粒子被吸收。\n";
        p.alive = false;
        break;
      case Boundary::Reflective:
        std::cout << "  击中反射边界 " << hit.surface->name() << "，方向翻转。\n";
        reflect(p.direction, hit.axis);
        break;
      case Boundary::Interface:
        std::cout << "  穿越界面 " << hit.surface->name()
                  << "（示例中仍在同一单元，保持状态）。\n";
        // 本教学场景只有一个单元，故不改变 Cell；可扩展到多单元切换。
        p.position = p.position + p.direction * 1e-7;  // 防止数值抖动
        break;
    }
  }

  void reflect(Vec3& dir, Axis axis) {
    if (axis == Axis::X) dir.x = -dir.x;
    if (axis == Axis::Y) dir.y = -dir.y;
    if (axis == Axis::Z) dir.z = -dir.z;
  }

  void collide(Particle& p, RandomEngine& rng) {
    double sigma_t = material_.sigma_t(p.group);
    double xi = rng.uniform() * sigma_t;
    double cumulative = material_.sigma_a(p.group);
    if (xi < cumulative) {
      std::cout << "  发生吸收，粒子死亡。\n";
      p.alive = false;
      return;
    }
    cumulative += material_.sigma_f(p.group);
    if (xi < cumulative) {
      std::size_t new_group = rng.pick(material_.chi());
      std::cout << "  发生裂变，产生二次粒子群 = " << new_group
                << "（本示例不追踪二次粒子）。\n";
      // 主粒子仍按散射处理，演示逻辑；真实 OpenMC 会将二次粒子写入银行。
    }
    // 默认进入散射：
    std::size_t scatter_group = rng.pick(material_.scatter_row(p.group));
    p.group = scatter_group;
    p.direction = normalise(rng.isotropic_direction());
    std::cout << "  发生散射，新群 = " << p.group << ", 新方向 = (" << p.direction.x << ", "
              << p.direction.y << ", " << p.direction.z << ")\n";
  }

  const Material& material_;
  const Cell& cell_;
};

}  // namespace lesson_physics

int main() {
  using namespace lesson_physics;

  // 1. 构造材料与几何：单一 8 cm 立方体，部分面反射部分真空。
  Material material("fuel", {0.02, 0.03}, {0.05, 0.02}, {2.5, 2.5},
                    {{0.03, 0.01}, {0.02, 0.04}}, {0.7, 0.3});

  Surface x0("x0", Axis::X, 0.0, Boundary::Reflective);
  Surface x1("x1", Axis::X, 8.0, Boundary::Interface);
  Surface y0("y0", Axis::Y, 0.0, Boundary::Reflective);
  Surface y1("y1", Axis::Y, 8.0, Boundary::Interface);
  Surface z0("z0", Axis::Z, 0.0, Boundary::Vacuum);
  Surface z1("z1", Axis::Z, 8.0, Boundary::Interface);

  Cell cell("fuel_cell", &x0, &x1, &y0, &y1, &z0, &z1, &material);

  // 2. 初始化粒子状态：位置在中心，方向倾向 z 轴。
  Particle particle;
  particle.position = {4.0, 4.0, 4.0};
  particle.direction = normalise({0.3, -0.1, 1.0});
  particle.group = 0;

  // 3. 执行物理驱动，逐事件打印处理流程。
  RandomEngine rng(2024u);
  PhysicsDriver driver(material, cell);
  driver.transport(particle, rng);
}
