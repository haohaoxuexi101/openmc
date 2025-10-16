#include <iostream>
#include <memory>
#include <string>
#include <vector>

// Lesson 19 作为汇总案例，将前面 Region/Universe 组合成完整的层级：
// - Region 保持纯虚接口，强调扩展性；
// - Cell/Universe 负责组织与描述；
// - main 输出层级与查询结果，帮助理解真实几何树遍历。

struct Position {
  double x;
  double y;
  double z;
};

class Region {
public:
  virtual ~Region() = default;
  virtual bool contains(const Position& r) const = 0;
  virtual std::string description() const = 0;
};

class ZHalfspace final : public Region {
public:
  ZHalfspace(double plane_z, bool less) : plane_z_{plane_z}, less_{less} {}

  bool contains(const Position& r) const override
  {
    return less_ ? (r.z <= plane_z_) : (r.z >= plane_z_);
  }

  std::string description() const override
  {
    return (less_ ? "z <= " : "z >= ") + std::to_string(plane_z_);
  }

private:
  double plane_z_;
  bool less_;
};

class CylinderRegion final : public Region {
public:
  explicit CylinderRegion(double radius) : radius_{radius} {}

  bool contains(const Position& r) const override
  {
    return (r.x * r.x + r.y * r.y) <= radius_ * radius_;
  }

  std::string description() const override
  {
    return "r <= " + std::to_string(radius_);
  }

private:
  double radius_;
};

class IntersectionRegion final : public Region {
public:
  explicit IntersectionRegion(std::vector<std::shared_ptr<Region>> children)
    : children_{std::move(children)}
  {
  }

  bool contains(const Position& r) const override
  {
    for (const auto& child : children_) {
      if (!child->contains(r)) return false;
    }
    return true;
  }

  std::string description() const override
  {
    std::string desc = "(";
    for (std::size_t i = 0; i < children_.size(); ++i) {
      if (i > 0) desc += " AND ";
      desc += children_[i]->description();
    }
    desc += ")";
    return desc;
  }

private:
  std::vector<std::shared_ptr<Region>> children_;
};

class Cell {
public:
  Cell(std::string name, std::shared_ptr<Region> region) : name_{std::move(name)}, region_{std::move(region)} {}

  bool contains(const Position& r) const { return region_->contains(r); }
  const std::string& name() const { return name_; }
  std::string region_description() const { return region_->description(); }

private:
  std::string name_;
  std::shared_ptr<Region> region_;
};

class Universe {
public:
  void add_cell(Cell cell) { cells_.push_back(std::move(cell)); }

  const Cell* find_cell(const Position& r) const
  {
    for (const auto& cell : cells_) {
      if (cell.contains(r)) return &cell;
    }
    return nullptr;
  }

  void describe() const
  {
    for (const auto& cell : cells_) {
      std::cout << cell.name() << " region " << cell.region_description() << '\n';
    }
  }

private:
  std::vector<Cell> cells_;
};

int main()
{
  auto fuel_region = std::make_shared<IntersectionRegion>(
    std::vector<std::shared_ptr<Region>>{
      std::make_shared<CylinderRegion>(0.4),
      std::make_shared<ZHalfspace>(10.0, true),
      std::make_shared<ZHalfspace>(0.0, false),
    });

  auto moderator_region = std::make_shared<IntersectionRegion>(
    std::vector<std::shared_ptr<Region>>{
      std::make_shared<CylinderRegion>(0.6),
      std::make_shared<ZHalfspace>(10.0, true),
      std::make_shared<ZHalfspace>(0.0, false),
    });

  Universe pin;
  pin.add_cell(Cell{"fuel", fuel_region});
  pin.add_cell(Cell{"moderator", moderator_region});

  pin.describe();

  Position samples[] = {{0.1, 0.1, 5.0}, {0.55, 0.0, 3.0}, {0.7, 0.0, 5.0}};
  for (const auto& pos : samples) {
    if (const Cell* cell = pin.find_cell(pos)) {
      std::cout << "Point (" << pos.x << ',' << pos.y << ',' << pos.z << ") -> " << cell->name()
                << '\n';
    } else {
      std::cout << "Point (" << pos.x << ',' << pos.y << ',' << pos.z << ") -> outside" << '\n';
    }
  }
}
