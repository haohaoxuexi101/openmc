#include <array>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// 本示例聚焦 OpenMC 的几何层级，将 Surface / Cell / Universe 的设计思路抽离成
// 独立教学代码。核心目标：
//   - 展示笛卡尔平面的封装方式；
//   - 说明 Cell 如何组合 6 个面并维护邻接关系；
//   - 演示 Universe 作为容器对 Cell 做统一管理。

namespace lesson_geometry {

enum class Boundary { Interface, Vacuum, Reflective };
enum class Axis { X = 0, Y = 1, Z = 2 };

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

class Cell {
 public:
  Cell(std::string name, const Surface* x_minus, const Surface* x_plus, const Surface* y_minus,
       const Surface* y_plus, const Surface* z_minus, const Surface* z_plus)
      : name_(std::move(name)) {
    surfaces_[0][0] = x_minus;
    surfaces_[0][1] = x_plus;
    surfaces_[1][0] = y_minus;
    surfaces_[1][1] = y_plus;
    surfaces_[2][0] = z_minus;
    surfaces_[2][1] = z_plus;
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }

  void set_neighbour(Axis axis, bool positive, const Cell* neighbour) {
    neighbours_[static_cast<int>(axis)][positive ? 1 : 0] = neighbour;
  }

  [[nodiscard]] const Cell* neighbour(Axis axis, bool positive) const noexcept {
    return neighbours_[static_cast<int>(axis)][positive ? 1 : 0];
  }

  struct BoundaryHit {
    double distance = std::numeric_limits<double>::infinity();
    const Surface* surface = nullptr;
    Axis axis = Axis::X;
    bool positive = true;
  };

  BoundaryHit distance_to_boundary(const std::array<double, 3>& position,
                                   const std::array<double, 3>& direction) const {
    BoundaryHit hit;
    for (int ax = 0; ax < 3; ++ax) {
      double component = direction[ax];
      if (std::abs(component) < 1e-14) continue;  // 防御式编程：避免除以零
      bool positive = component > 0.0;
      const Surface* s = surfaces_[ax][positive ? 1 : 0];
      double coord = s->coordinate();
      double pos = position[ax];
      double distance = (coord - pos) / component;
      if (distance <= 1e-14) continue;  // 保证距离为正数
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
  const Surface* surfaces_[3][2] = {};
  const Cell* neighbours_[3][2] = {};
};

class Universe {
 public:
  explicit Universe(std::string name) : name_(std::move(name)) {}

  Cell& add_cell(Cell cell) {
    cells_.push_back(std::move(cell));
    return cells_.back();
  }

  [[nodiscard]] const std::vector<Cell>& cells() const noexcept { return cells_; }

 private:
  std::string name_;
  std::vector<Cell> cells_;
};

}  // namespace lesson_geometry

int main() {
  using namespace lesson_geometry;

  // 1. 构造一个 10 cm 立方体：六个平面各自标记边界条件。
  Surface x0("x0", Axis::X, 0.0, Boundary::Reflective);
  Surface x1("x1", Axis::X, 10.0, Boundary::Interface);
  Surface y0("y0", Axis::Y, 0.0, Boundary::Reflective);
  Surface y1("y1", Axis::Y, 10.0, Boundary::Interface);
  Surface z0("z0", Axis::Z, 0.0, Boundary::Vacuum);
  Surface z1("z1", Axis::Z, 10.0, Boundary::Interface);

  // 2. 创建一个 Cell，并展示如何设置邻接（此例仅在 +X 方向连接一个虚拟邻居）。
  Universe universe("fuel_universe");
  Cell& cell = universe.add_cell(Cell("fuel_cell", &x0, &x1, &y0, &y1, &z0, &z1));
  cell.set_neighbour(Axis::X, true, nullptr);  // 无实际邻居，保持空指针表示边界

  // 3. 计算一个粒子到最近边界的距离，体现几何搜索的算法。
  std::array<double, 3> position{5.0, 5.0, 5.0};
  std::array<double, 3> direction{0.2, -0.4, 0.9};
  Cell::BoundaryHit hit = cell.distance_to_boundary(position, direction);

  std::cout << "粒子位于 Cell: " << cell.name() << "\n";
  std::cout << "最近边界 = " << hit.surface->name() << ", 距离 = " << hit.distance << " cm\n";
  std::cout << "该边界条件 = " << (hit.surface->boundary() == Boundary::Reflective ? "反射" :
                                      hit.surface->boundary() == Boundary::Vacuum ? "真空" : "界面")
            << "\n";
}
