#include <array>
#include <iostream>
#include <string>
#include <vector>

// 本示例专门拆解 OpenMC 在几何定位阶段的“宇宙 + 晶格 + 几何”组合模式。
// 目标是理解粒子如何从坐标映射到 Cell，以及为何需要 Universe/Lattice 的层次。

namespace lesson_lattice {

enum class Boundary { Interface, Vacuum, Reflective };
enum class Axis { X = 0, Y = 1, Z = 2 };

class Surface {
 public:
  Surface(std::string name, Axis axis, double coord, Boundary boundary)
      : name_(std::move(name)), axis_(axis), coordinate_(coord), boundary_(boundary) {}

  [[nodiscard]] double coordinate() const noexcept { return coordinate_; }
  [[nodiscard]] const std::string& name() const noexcept { return name_; }

 private:
  std::string name_;
  Axis axis_;
  double coordinate_;
  Boundary boundary_;
};

class Cell {
 public:
  Cell(std::string name, int id, const Surface* x_minus, const Surface* x_plus,
       const Surface* y_minus, const Surface* y_plus, const Surface* z_minus,
       const Surface* z_plus)
      : name_(std::move(name)), id_(id) {
    surfaces_[0][0] = x_minus;
    surfaces_[0][1] = x_plus;
    surfaces_[1][0] = y_minus;
    surfaces_[1][1] = y_plus;
    surfaces_[2][0] = z_minus;
    surfaces_[2][1] = z_plus;
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] int id() const noexcept { return id_; }

  [[nodiscard]] bool contains(const std::array<double, 3>& r) const noexcept {
    return r[0] >= surfaces_[0][0]->coordinate() && r[0] < surfaces_[0][1]->coordinate() &&
           r[1] >= surfaces_[1][0]->coordinate() && r[1] < surfaces_[1][1]->coordinate() &&
           r[2] >= surfaces_[2][0]->coordinate() && r[2] < surfaces_[2][1]->coordinate();
  }

 private:
  std::string name_;
  int id_;
  const Surface* surfaces_[3][2] = {};
};

class Universe {
 public:
  explicit Universe(std::string name) : name_(std::move(name)) {}

  Cell& add_cell(Cell cell) {
    cells_.push_back(std::move(cell));
    return cells_.back();
  }

  [[nodiscard]] const std::vector<Cell>& cells() const noexcept { return cells_; }

  [[nodiscard]] const Cell* find_cell(const std::array<double, 3>& r) const noexcept {
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
  Lattice3D(std::string name, std::array<double, 3> pitch, std::array<std::size_t, 3> dims,
            std::array<double, 3> origin)
      : name_(std::move(name)), pitch_(pitch), dims_(dims), origin_(origin) {
    universes_.resize(dims_[0] * dims_[1] * dims_[2], nullptr);
  }

  void set_universe(std::size_t i, std::size_t j, std::size_t k, const Universe* universe) {
    universes_.at(index(i, j, k)) = universe;
  }

  [[nodiscard]] const Universe* universe_at(const std::array<double, 3>& r) const noexcept {
    double lx = r[0] - origin_[0];
    double ly = r[1] - origin_[1];
    double lz = r[2] - origin_[2];
    if (lx < 0.0 || ly < 0.0 || lz < 0.0) return nullptr;
    std::size_t i = static_cast<std::size_t>(lx / pitch_[0]);
    std::size_t j = static_cast<std::size_t>(ly / pitch_[1]);
    std::size_t k = static_cast<std::size_t>(lz / pitch_[2]);
    if (i >= dims_[0] || j >= dims_[1] || k >= dims_[2]) return nullptr;
    return universes_[index(i, j, k)];
  }

 private:
  [[nodiscard]] std::size_t index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
    return (k * dims_[1] + j) * dims_[0] + i;
  }

  std::string name_;
  std::array<double, 3> pitch_;
  std::array<std::size_t, 3> dims_;
  std::array<double, 3> origin_;
  std::vector<const Universe*> universes_;
};

class Geometry {
 public:
  Geometry(const Universe& root, const Lattice3D& lattice)
      : root_(root), lattice_(lattice) {}

  [[nodiscard]] const Cell* locate(const std::array<double, 3>& r) const noexcept {
    if (const Universe* u = lattice_.universe_at(r)) {
      if (const Cell* c = u->find_cell(r)) return c;
    }
    return root_.find_cell(r);
  }

 private:
  const Universe& root_;
  const Lattice3D& lattice_;
};

}  // namespace lesson_lattice

int main() {
  using namespace lesson_lattice;

  // 1. 构建表面：每个 Cell 共享这些平面，体现 OpenMC 的“共享几何实体”思想。
  Surface x0("x0", Axis::X, 0.0, Boundary::Vacuum);
  Surface x1("x1", Axis::X, 2.0, Boundary::Interface);
  Surface x2("x2", Axis::X, 4.0, Boundary::Vacuum);
  Surface y0("y0", Axis::Y, 0.0, Boundary::Vacuum);
  Surface y1("y1", Axis::Y, 2.0, Boundary::Interface);
  Surface y2("y2", Axis::Y, 4.0, Boundary::Vacuum);
  Surface z0("z0", Axis::Z, 0.0, Boundary::Reflective);
  Surface z1("z1", Axis::Z, 3.0, Boundary::Interface);
  Surface z2("z2", Axis::Z, 6.0, Boundary::Vacuum);

  // 2. 创建 2×2×2 晶格，仿照主示例中的 DemoModel。
  Universe root("root");
  Universe u000("u000");
  Universe u100("u100");
  Universe u010("u010");
  Universe u110("u110");
  Universe u001("u001");
  Universe u101("u101");
  Universe u011("u011");
  Universe u111("u111");

  u000.add_cell({"fuel_lower_fresh", 0, &x0, &x1, &y0, &y1, &z0, &z1});
  u100.add_cell({"fuel_lower_burned", 1, &x1, &x2, &y0, &y1, &z0, &z1});
  u010.add_cell({"reflector_lower_front", 2, &x0, &x1, &y1, &y2, &z0, &z1});
  u110.add_cell({"reflector_lower_back", 3, &x1, &x2, &y1, &y2, &z0, &z1});
  u001.add_cell({"moderator_upper_left", 4, &x0, &x1, &y0, &y1, &z1, &z2});
  u101.add_cell({"moderator_upper_right", 5, &x1, &x2, &y0, &y1, &z1, &z2});
  u011.add_cell({"reflector_upper_front", 6, &x0, &x1, &y1, &y2, &z1, &z2});
  u111.add_cell({"reflector_upper_back", 7, &x1, &x2, &y1, &y2, &z1, &z2});

  Lattice3D lattice("assembly", {2.0, 2.0, 3.0}, {2, 2, 2}, {0.0, 0.0, 0.0});
  lattice.set_universe(0, 0, 0, &u000);
  lattice.set_universe(1, 0, 0, &u100);
  lattice.set_universe(0, 1, 0, &u010);
  lattice.set_universe(1, 1, 0, &u110);
  lattice.set_universe(0, 0, 1, &u001);
  lattice.set_universe(1, 0, 1, &u101);
  lattice.set_universe(0, 1, 1, &u011);
  lattice.set_universe(1, 1, 1, &u111);

  Geometry geometry(root, lattice);

  // 3. 通过 locate 演示粒子定位。坐标选自不同单元中点，便于人工验证。
  std::array<std::array<double, 3>, 3> samples{{{1.0, 1.0, 1.0}, {3.5, 0.5, 1.5}, {0.5, 3.5, 4.5}}};

  for (const auto& r : samples) {
    if (const Cell* cell = geometry.locate(r)) {
      std::cout << "位置 (" << r[0] << ", " << r[1] << ", " << r[2]
                << ") 位于单元 " << cell->name() << " (ID=" << cell->id() << ")\n";
    } else {
      std::cout << "位置 (" << r[0] << ", " << r[1] << ", " << r[2]
                << ") 未命中任何单元（落在根宇宙外）\n";
    }
  }
}
