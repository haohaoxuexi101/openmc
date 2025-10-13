#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// 本课程演示如何用“接口 + 继承 + 组合”复刻 OpenMC 的面向对象骨架。
// 核心思路：
//   * 用抽象基类 `Surface` 承担接口职责，派生类只补充几何细节；
//   * `Cell` 组合多个 `Surface` 并记录邻接关系，体现“组合优于继承”的建模方式；
//   * `Particle` 仅保存必要的输运状态，通过接口向几何/材料请求信息；
//   * `Geometry` 集中管理对象生命周期，外部算法只与公开接口交互；
//   * `TransportInspector` 作为策略类，说明如何将“推进算法”与数据对象解耦。

namespace lesson_oop_architecture {

// ============================= 基础向量工具 =============================
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  Vec3 operator+(const Vec3& other) const { return {x + other.x, y + other.y, z + other.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) {
  double n = norm(v);
  if (n == 0.0) return {0.0, 0.0, 1.0};
  return {v.x / n, v.y / n, v.z / n};
}

// ============================= 几何接口层 =============================
class Surface {
 public:
  explicit Surface(std::string name) : name_(std::move(name)) {}
  virtual ~Surface() = default;

  virtual double distance(const Vec3& pos, const Vec3& dir) const = 0;
  virtual Vec3 normal(const Vec3& pos) const = 0;

  const std::string& name() const { return name_; }

 private:
  std::string name_;
};

class AxisAlignedPlane : public Surface {
 public:
  enum class Axis { X, Y, Z };

  AxisAlignedPlane(std::string name, Axis axis, double position)
      : Surface(std::move(name)), axis_(axis), position_(position) {}

  double distance(const Vec3& pos, const Vec3& dir) const override {
    double numerator = position_ - component(pos);
    double denominator = component(dir);
    if (std::abs(denominator) < 1e-14) {
      return std::numeric_limits<double>::infinity();
    }
    double d = numerator / denominator;
    if (d <= 1e-12) {
      return std::numeric_limits<double>::infinity();
    }
    return d;
  }

  Vec3 normal(const Vec3&) const override { return axis_direction(); }

  Axis axis() const { return axis_; }
  double position() const { return position_; }

  double evaluate(const Vec3& pos) const { return component(pos) - position_; }
  Vec3 axis_direction() const {
    switch (axis_) {
      case Axis::X:
        return {1.0, 0.0, 0.0};
      case Axis::Y:
        return {0.0, 1.0, 0.0};
      case Axis::Z:
      default:
        return {0.0, 0.0, 1.0};
    }
  }

 private:
  double component(const Vec3& v) const {
    switch (axis_) {
      case Axis::X:
        return v.x;
      case Axis::Y:
        return v.y;
      case Axis::Z:
      default:
        return v.z;
    }
  }

  Axis axis_;
  double position_;
};

// ============================= 材料对象（组合示例） =============================
class Material {
 public:
  Material(std::string name, std::vector<double> sigma_t)
      : name_(std::move(name)), sigma_t_(std::move(sigma_t)) {
    if (sigma_t_.empty()) {
      throw std::invalid_argument("材料至少需要一个能群的总截面");
    }
  }

  const std::string& name() const { return name_; }
  double sigma_t(std::size_t g) const { return sigma_t_.at(g); }

 private:
  std::string name_;
  std::vector<double> sigma_t_;
};

// ============================= Cell：组合 + 邻接 =============================
class Cell {
 public:
  enum class BoundaryType { Interface, Vacuum, Reflective };

  struct Face {
    const AxisAlignedPlane* plane = nullptr;
    int sense = +1;  // +1 表示内部位于 plane 的正侧（evaluate >= 0），-1 表示内部位于负侧。
    BoundaryType boundary = BoundaryType::Interface;
    Cell* neighbor = nullptr;  // 界面时指向相邻单元。
  };

  Cell(std::string name, const Material* material) : name_(std::move(name)), material_(material) {}

  Face& add_face(const AxisAlignedPlane* plane, int sense, BoundaryType type, Cell* neighbor = nullptr) {
    faces_.push_back(Face{plane, sense, type, neighbor});
    return faces_.back();
  }

  const std::string& name() const { return name_; }
  const Material* material() const { return material_; }
  const std::vector<Face>& faces() const { return faces_; }
  Face& face_mut(std::size_t index) { return faces_.at(index); }

  bool contains(const Vec3& pos) const {
    for (const auto& face : faces_) {
      double signed_value = static_cast<double>(face.sense) * face.plane->evaluate(pos);
      if (signed_value < -1e-10) {
        return false;
      }
    }
    return true;
  }

  std::pair<double, const Face*> distance_to_boundary(const Vec3& pos, const Vec3& dir) const {
    double best = std::numeric_limits<double>::infinity();
    const Face* hit = nullptr;
    for (const auto& face : faces_) {
      Vec3 outward = outward_normal(face);
      double projection = dot(outward, dir);
      if (projection <= 1e-14) {
        continue;  // 粒子朝向内部或平行于该面。
      }
      double d = face.plane->distance(pos, dir);
      if (d < best) {
        best = d;
        hit = &face;
      }
    }
    return {best, hit};
  }

 private:
  Vec3 outward_normal(const Face& face) const {
    Vec3 axis_dir = face.plane->axis_direction();
    if (face.sense > 0) {
      // 内部在正侧，则外法向指向负侧。
      return {-axis_dir.x, -axis_dir.y, -axis_dir.z};
    }
    return axis_dir;
  }

  std::string name_;
  const Material* material_;
  std::vector<Face> faces_;
};

// ============================= 粒子对象 =============================
struct Particle {
  Vec3 position;
  Vec3 direction;
  std::size_t group = 0;
  const Cell* cell = nullptr;
};

// ============================= Geometry 管理器 =============================
class Geometry {
 public:
  template <typename SurfaceT, typename... Args>
  SurfaceT* make_surface(Args&&... args) {
    auto ptr = std::make_shared<SurfaceT>(std::forward<Args>(args)...);
    SurfaceT* raw = ptr.get();
    surfaces_.push_back(std::move(ptr));
    return raw;
  }

  Cell* make_cell(std::string name, const Material* material) {
    cells_.push_back(std::make_unique<Cell>(std::move(name), material));
    return cells_.back().get();
  }

  const Cell* locate(const Vec3& pos) const {
    for (const auto& cell : cells_) {
      if (cell->contains(pos)) {
        return cell.get();
      }
    }
    return nullptr;
  }

 private:
  std::vector<std::shared_ptr<Surface>> surfaces_;
  std::vector<std::unique_ptr<Cell>> cells_;
};

// ============================= 算法策略：粒子推进检查器 =============================
class TransportInspector {
 public:
  explicit TransportInspector(const Geometry& geometry) : geometry_(geometry) {}

  void move(Particle& particle, double requested_distance) const {
    auto [distance, face] = particle.cell->distance_to_boundary(particle.position, particle.direction);
    if (!face || distance > requested_distance) {
      particle.position = particle.position + particle.direction * requested_distance;
      std::cout << "粒子在单元内平移 " << requested_distance << " cm，当前位置 = (" << particle.position.x
                << ", " << particle.position.y << ", " << particle.position.z << ")\n";
      return;
    }

    std::cout << "粒子在 " << distance << " cm 处击中面 " << face->plane->name() << "，边界类型 = ";
    switch (face->boundary) {
      case Cell::BoundaryType::Vacuum:
        std::cout << "真空（终止历史）\n";
        break;
      case Cell::BoundaryType::Reflective:
        std::cout << "反射（翻转方向）\n";
        handle_reflection(particle, *face, distance);
        break;
      case Cell::BoundaryType::Interface:
        std::cout << "界面（跳入邻接单元）\n";
        handle_interface(particle, *face, distance);
        break;
    }
  }

 private:
  void handle_reflection(Particle& particle, const Cell::Face& face, double distance) const {
    particle.position = particle.position + particle.direction * (distance - 1e-9);
    Vec3 axis = face.plane->axis_direction();
    double component = dot(particle.direction, axis);
    particle.direction = particle.direction + axis * (-2.0 * component);
    particle.direction = normalize(particle.direction);
  }

  void handle_interface(Particle& particle, const Cell::Face& face, double distance) const {
    particle.position = particle.position + particle.direction * (distance + 1e-9);
    if (face.neighbor) {
      particle.cell = face.neighbor;
      std::cout << "粒子进入单元 " << particle.cell->name() << "，材料 = " << particle.cell->material()->name()
                << "\n";
    } else {
      std::cout << "未设置邻居，视为真空终止\n";
    }
  }

  const Geometry& geometry_;
};

// ============================= 示例模型构建 =============================
struct DemoModel {
  Geometry geometry;
  std::vector<Material> materials;
  Cell* fuel_cell = nullptr;
  Cell* moderator_cell = nullptr;

  DemoModel() {
    materials.emplace_back("fresh_fuel", std::vector<double>{1.2, 0.9, 0.6});
    materials.emplace_back("moderator", std::vector<double>{0.4, 0.3, 0.2});

    auto* x0 = geometry.make_surface<AxisAlignedPlane>("x0", AxisAlignedPlane::Axis::X, 0.0);
    auto* x_mid = geometry.make_surface<AxisAlignedPlane>("x_mid", AxisAlignedPlane::Axis::X, 0.5);
    auto* x1 = geometry.make_surface<AxisAlignedPlane>("x1", AxisAlignedPlane::Axis::X, 1.0);
    auto* y0 = geometry.make_surface<AxisAlignedPlane>("y0", AxisAlignedPlane::Axis::Y, 0.0);
    auto* y1 = geometry.make_surface<AxisAlignedPlane>("y1", AxisAlignedPlane::Axis::Y, 1.0);
    auto* z0 = geometry.make_surface<AxisAlignedPlane>("z0", AxisAlignedPlane::Axis::Z, 0.0);
    auto* z1 = geometry.make_surface<AxisAlignedPlane>("z1", AxisAlignedPlane::Axis::Z, 1.0);

    fuel_cell = geometry.make_cell("fuel_cell", &materials[0]);
    fuel_cell->add_face(x0, +1, Cell::BoundaryType::Reflective);
    fuel_cell->add_face(x_mid, -1, Cell::BoundaryType::Interface);
    fuel_cell->add_face(y0, +1, Cell::BoundaryType::Reflective);
    fuel_cell->add_face(y1, -1, Cell::BoundaryType::Reflective);
    fuel_cell->add_face(z0, +1, Cell::BoundaryType::Reflective);
    fuel_cell->add_face(z1, -1, Cell::BoundaryType::Reflective);

    moderator_cell = geometry.make_cell("moderator_cell", &materials[1]);
    moderator_cell->add_face(x_mid, +1, Cell::BoundaryType::Interface, fuel_cell);
    moderator_cell->add_face(x1, -1, Cell::BoundaryType::Vacuum);
    moderator_cell->add_face(y0, +1, Cell::BoundaryType::Reflective);
    moderator_cell->add_face(y1, -1, Cell::BoundaryType::Reflective);
    moderator_cell->add_face(z0, +1, Cell::BoundaryType::Reflective);
    moderator_cell->add_face(z1, -1, Cell::BoundaryType::Reflective);

    // 建立界面互指。
    fuel_cell->face_mut(1).neighbor = moderator_cell;
  }
};

// 上面 fuel->faces() 返回 const 引用，无法修改 neighbor。为了演示，我们为 DemoModel 定义
// 辅助函数，以安全方式写入邻接。更好的做法是让 Cell 暴露 face_mut() 接口。
inline void link_neighbors(Cell* a, std::size_t face_a, Cell* b, std::size_t face_b) {
  a->face_mut(face_a).neighbor = b;
  b->face_mut(face_b).neighbor = a;
}

void run_demo() {
  DemoModel model;
  Cell* fuel = model.fuel_cell;
  Cell* moderator = model.moderator_cell;

  // 绑定界面邻接，并验证 Geometry::locate 能够找到同一对象。
  link_neighbors(fuel, 1, moderator, 0);
  const Cell* located = model.geometry.locate({0.2, 0.3, 0.4});
  if (located != fuel) {
    std::cerr << "定位检查失败：期望 fuel_cell" << std::endl;
  }

  Particle particle;
  particle.position = {0.2, 0.5, 0.5};
  particle.direction = normalize({0.6, 0.1, 0.0});
  particle.group = 0;
  particle.cell = fuel;

  std::cout << "初始单元: " << particle.cell->name() << " 材料 = " << particle.cell->material()->name()
            << "\n";

  TransportInspector inspector(model.geometry);
  inspector.move(particle, 0.4);
  inspector.move(particle, 0.4);
}

}  // namespace lesson_oop_architecture

int main() {
  using namespace lesson_oop_architecture;
  std::cout << std::fixed << std::setprecision(5);
  run_demo();
}
