#include <iostream>
#include <memory>
#include <string>
#include <vector>

// 本课展示“布尔组合 + 多态区域”的设计：
// - Region 定义 contains/describe 接口，强调语义清晰的 OO 抽象。
// - 组合类 IntersectionRegion/UnionRegion 像 OpenMC 的 CSG 布尔运算。
// - main 演示如何拼装燃料/包壳区域并打印判断结果。

struct Position {
  double x;
  double y;
  double z;
};

class Region {
public:
  virtual ~Region() = default;
  virtual bool contains(const Position& r) const = 0;
  virtual std::string describe() const = 0;
};

class ZHalfspace final : public Region {
public:
  ZHalfspace(double plane_z, bool less) : plane_z_{plane_z}, less_{less} {}

  bool contains(const Position& r) const override
  {
    return less_ ? (r.z <= plane_z_) : (r.z >= plane_z_);
  }

  std::string describe() const override
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

  std::string describe() const override
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

  std::string describe() const override
  {
    std::string text;
    for (std::size_t i = 0; i < children_.size(); ++i) {
      if (i) text += " AND ";
      text += children_[i]->describe();
    }
    return text;
  }

private:
  std::vector<std::shared_ptr<Region>> children_;
};

class UnionRegion final : public Region {
public:
  explicit UnionRegion(std::vector<std::shared_ptr<Region>> children)
    : children_{std::move(children)}
  {
  }

  bool contains(const Position& r) const override
  {
    for (const auto& child : children_) {
      if (child->contains(r)) return true;
    }
    return false;
  }

  std::string describe() const override
  {
    std::string text;
    for (std::size_t i = 0; i < children_.size(); ++i) {
      if (i) text += " OR ";
      text += children_[i]->describe();
    }
    return text;
  }

private:
  std::vector<std::shared_ptr<Region>> children_;
};

int main()
{
  auto axial_bounds = std::make_shared<IntersectionRegion>(
    std::vector<std::shared_ptr<Region>>{
      std::make_shared<ZHalfspace>(0.0, false),
      std::make_shared<ZHalfspace>(10.0, true),
    });

  auto fuel = std::make_shared<IntersectionRegion>(
    std::vector<std::shared_ptr<Region>>{
      std::make_shared<CylinderRegion>(0.4),
      axial_bounds,
    });

  auto cladding = std::make_shared<IntersectionRegion>(
    std::vector<std::shared_ptr<Region>>{
      std::make_shared<CylinderRegion>(0.45),
      axial_bounds,
    });

  UnionRegion pin_stack{std::vector<std::shared_ptr<Region>>{fuel, cladding}};

  std::cout << "Pin region description: " << pin_stack.describe() << '\n';

  std::vector<Position> samples{{0.3, 0.0, 5.0}, {0.43, 0.0, 5.0}, {0.5, 0.0, 5.0}};
  for (const auto& s : samples) {
    std::cout << "Point (" << s.x << ',' << s.y << ',' << s.z << ") -> "
              << (pin_stack.contains(s) ? "inside" : "outside") << '\n';
  }
}
