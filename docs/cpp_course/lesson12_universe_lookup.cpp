#include <iostream>
#include <memory>
#include <string>
#include <vector>

// 本课延续上一节：
// - Cell 保存名称 + 区域指针，实现职责单一的数据结构。
// - Universe 只负责线性搜索并打印诊断，模拟 OpenMC 几何树顶层遍历。
// - main 演示在层次结构中查找点属于哪一个 cell。

struct Position {
  double x;
  double y;
  double z;
};

class Region {
public:
  virtual ~Region() = default;
  virtual bool contains(const Position& r) const = 0;
};

class CylinderRegion final : public Region {
public:
  explicit CylinderRegion(double radius) : radius_{radius} {}

  bool contains(const Position& r) const override
  {
    return (r.x * r.x + r.y * r.y) <= radius_ * radius_;
  }

private:
  double radius_;
};

class ZSlab final : public Region {
public:
  ZSlab(double lower, double upper) : lower_{lower}, upper_{upper} {}

  bool contains(const Position& r) const override
  {
    return r.z >= lower_ && r.z <= upper_;
  }

private:
  double lower_;
  double upper_;
};

class AndRegion final : public Region {
public:
  AndRegion(std::shared_ptr<Region> a, std::shared_ptr<Region> b)
    : a_{std::move(a)}
    , b_{std::move(b)}
  {
  }

  bool contains(const Position& r) const override
  {
    return a_->contains(r) && b_->contains(r);
  }

private:
  std::shared_ptr<Region> a_;
  std::shared_ptr<Region> b_;
};

class Cell {
public:
  Cell(std::string name, std::shared_ptr<Region> region)
    : name_{std::move(name)}
    , region_{std::move(region)}
  {
  }

  bool contains(const Position& r) const { return region_->contains(r); }
  const std::string& name() const { return name_; }

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

  void dump() const
  {
    std::cout << "Universe contains " << cells_.size() << " cells." << '\n';
  }

private:
  std::vector<Cell> cells_;
};

int main()
{
  auto axial = std::make_shared<ZSlab>(0.0, 10.0);
  Universe pin;
  pin.add_cell(Cell{"fuel", std::make_shared<AndRegion>(std::make_shared<CylinderRegion>(0.4), axial)});
  pin.add_cell(Cell{"moderator", std::make_shared<AndRegion>(std::make_shared<CylinderRegion>(0.6), axial)});

  pin.dump();

  std::vector<Position> tests{{0.2, 0.0, 5.0}, {0.5, 0.0, 2.0}, {0.7, 0.0, 5.0}};
  for (const auto& pos : tests) {
    if (const Cell* c = pin.find_cell(pos)) {
      std::cout << "Point -> " << c->name() << '\n';
    } else {
      std::cout << "Point -> outside" << '\n';
    }
  }
}
