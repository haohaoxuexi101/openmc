#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

// 课程 4：组合关系——Cell 与 Material 绑定
// 借鉴说明：Cell 的材料指针来自 OpenMC `src/geometry/cell.h` 中的 material_ 字段，
// MaterialRegistry 的查询模式则参考 `src/material.cpp` 里的材料表注册逻辑。
// 教学实现保持 Lesson 前缀类名即可避开与正式类冲突，因此不再使用多层命名空间。

class Material {
public:
  Material(std::string name, double xs_total, double xs_absorb)
      : name_(std::move(name)), xs_total_(xs_total), xs_absorb_(xs_absorb) {}

  const std::string& name() const { return name_; }
  double xs_total() const { return xs_total_; }
  double xs_absorb() const { return xs_absorb_; }

private:
  std::string name_;
  double xs_total_;
  double xs_absorb_;
};

class MaterialRegistry {
public:
  void add(std::shared_ptr<Material> mat) {
    materials_[mat->name()] = std::move(mat);
  }

  std::shared_ptr<Material> get(const std::string& name) const {
    return materials_.at(name);
  }

private:
  std::map<std::string, std::shared_ptr<Material>> materials_;
};

class Cell {
public:
  Cell(std::string name, std::string material_name)
      : name_(std::move(name)), material_name_(std::move(material_name)) {}

  const std::string& name() const { return name_; }
  const std::string& material_name() const { return material_name_; }

private:
  std::string name_;
  std::string material_name_;
};

class LessonGeometryStateBinding {
public:
  LessonGeometryStateBinding(MaterialRegistry registry, std::vector<Cell> cells)
      : registry_(std::move(registry)), cells_(std::move(cells)) {}

  void report() const {
    std::cout << "几何与材料绑定一览:\n";
    for (const auto& cell : cells_) {
      auto mat = registry_.get(cell.material_name());
      std::cout << "  Cell " << cell.name() << " -> Material " << mat->name()
                << " (Σt=" << mat->xs_total() << ", Σa=" << mat->xs_absorb()
                << ")\n";
    }
  }

private:
  MaterialRegistry registry_;
  std::vector<Cell> cells_;
};

int main() {
  MaterialRegistry registry;
  registry.add(std::make_shared<Material>("燃料", 1.2, 0.5));
  registry.add(std::make_shared<Material>("慢化剂", 0.7, 0.01));

  LessonGeometryStateBinding geometry(
      registry, {Cell("燃料芯块", "燃料"), Cell("慢化剂包壳", "慢化剂")});

  geometry.report();
  return 0;
}
