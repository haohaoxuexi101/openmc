#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>

// 课程 1：几何状态（GeometryState）的对象建模示例
// 借鉴说明：本节直接对照 OpenMC 源码 `src/geometry/surface.h` 与
// `src/geometry/cell.h` 中 Surface 基类及 Cell 边界判定逻辑的接口拆解，
// 通过更小型的 LessonGeometryState 组合器复刻 "面 -> 单元 -> 几何状态"
// 的责任链。为避免与正式源码命名冲突，示例类均使用 Lesson 前缀并保持
// 独立实现，不再引入额外命名空间。

// ------------------ 面定义层（抽象接口） ------------------
class Surface {
public:
  virtual ~Surface() = default;
  // 判定粒子位于表面法向的哪一侧：>0 表示正侧，<0 表示负侧
  virtual double evaluate(double x, double y, double z) const = 0;
  virtual std::string name() const = 0;
};

// 平面：OpenMC 源码中 Plane 族的极简版
class Plane : public Surface {
public:
  Plane(std::string name, double a, double b, double c, double d)
      : name_(std::move(name)), a_(a), b_(b), c_(c), d_(d) {}

  double evaluate(double x, double y, double z) const override {
    return a_ * x + b_ * y + c_ * z - d_;
  }

  std::string name() const override { return name_; }

private:
  std::string name_;
  double a_, b_, c_, d_;
};

// 球面：演示继承提供不同 evaluate 实现
class Sphere : public Surface {
public:
  Sphere(std::string name, double x0, double y0, double z0, double r)
      : name_(std::move(name)), x0_(x0), y0_(y0), z0_(z0), r_(r) {}

  double evaluate(double x, double y, double z) const override {
    double dx = x - x0_;
    double dy = y - y0_;
    double dz = z - z0_;
    return std::sqrt(dx * dx + dy * dy + dz * dz) - r_;
  }

  std::string name() const override { return name_; }

private:
  std::string name_;
  double x0_, y0_, z0_, r_;
};

// ------------------ 单元定义层（组合接口） ------------------
class Cell {
public:
  Cell(std::string name, std::vector<std::shared_ptr<Surface>> surfaces)
      : name_(std::move(name)), surfaces_(std::move(surfaces)) {}

  const std::string& name() const { return name_; }

  // 判断点是否位于该 Cell 内部：所有面值均小于等于 0
  bool contains(double x, double y, double z) const {
    for (const auto& surf : surfaces_) {
      if (surf->evaluate(x, y, z) > 0.0) {
        return false;
      }
    }
    return true;
  }

private:
  std::string name_;
  std::vector<std::shared_ptr<Surface>> surfaces_;
};

// ------------------ 几何状态：粒子所处的空间单元 ------------------
class LessonGeometryState {
public:
  explicit LessonGeometryState(std::vector<Cell> cells)
      : cells_(std::move(cells)) {}

  // 给定空间坐标，定位粒子所在的 Cell
  const Cell* locate(double x, double y, double z) const {
    for (const auto& cell : cells_) {
      if (cell.contains(x, y, z)) {
        return &cell;
      }
    }
    return nullptr;  // 未找到说明粒子逃逸
  }

private:
  std::vector<Cell> cells_;
};

int main() {
  // 准备几何：一个立方体和内嵌球体，体现继承 + 组合
  auto px = std::make_shared<Plane>("px", 1.0, 0.0, 0.0, 1.0);
  auto nx = std::make_shared<Plane>("nx", -1.0, 0.0, 0.0, 1.0);
  auto py = std::make_shared<Plane>("py", 0.0, 1.0, 0.0, 1.0);
  auto ny = std::make_shared<Plane>("ny", 0.0, -1.0, 0.0, 1.0);
  auto pz = std::make_shared<Plane>("pz", 0.0, 0.0, 1.0, 1.0);
  auto nz = std::make_shared<Plane>("nz", 0.0, 0.0, -1.0, 1.0);
  auto sphere = std::make_shared<Sphere>("fuel", 0.0, 0.0, 0.0, 0.6);

  Cell fuel("燃料球", {px, nx, py, ny, pz, nz, sphere});
  Cell moderator("慢化剂方块", {px, nx, py, ny, pz, nz});

  LessonGeometryState geometry({fuel, moderator});

  // 测试定位
  std::vector<std::tuple<double, double, double>> test_points = {
      {0.0, 0.0, 0.0}, {0.5, 0.5, 0.5}, {0.8, 0.0, 0.0}, {1.2, 0.0, 0.0}};

  for (const auto& [x, y, z] : test_points) {
    const Cell* cell = geometry.locate(x, y, z);
    if (cell) {
      std::cout << "粒子位于 Cell: " << cell->name() << "\n";
    } else {
      std::cout << "粒子已经逃逸几何\n";
    }
  }

  return 0;
}
