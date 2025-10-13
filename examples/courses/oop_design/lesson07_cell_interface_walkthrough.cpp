#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// 课程 7：Cell 布尔树与距离接口拆解
// -----------------------------------
// 借鉴说明：Region/Cell 结构对照 OpenMC `src/geometry/cell.h` 中 Region 派生类
// 与 distance 计算的入口逻辑；示例将布尔树精简为 Halfspace + Intersection，
// 用以说明 contains/distance 的职责划分，而材料索引与最近交距迭代仍保留
// 源码中的核心思想。

// ------------------ 几何元素定义 ------------------
class Surface {
public:
  virtual ~Surface() = default;
  virtual double evaluate(double x, double y, double z) const = 0;
  virtual double distance(double x, double y, double z, double u,
                          double v, double w) const = 0;
};

class Plane : public Surface {
public:
  Plane(std::string name, double nx, double ny, double nz, double d)
      : name_(std::move(name)), nx_(nx), ny_(ny), nz_(nz), d_(d) {}

  double evaluate(double x, double y, double z) const override {
    return nx_ * x + ny_ * y + nz_ * z - d_;
  }

  double distance(double x, double y, double z, double u, double v,
                  double w) const override {
    double denom = nx_ * u + ny_ * v + nz_ * w;
    if (std::abs(denom) < 1e-12) {
      return std::numeric_limits<double>::infinity();
    }
    double t = (d_ - (nx_ * x + ny_ * y + nz_ * z)) / denom;
    return t > 0.0 ? t : std::numeric_limits<double>::infinity();
  }

private:
  std::string name_;
  double nx_, ny_, nz_, d_;
};

// ------------------ Region 布尔节点 ------------------
struct Position {
  double x, y, z;
};

class Region {
public:
  virtual ~Region() = default;
  virtual bool contains(Position r) const = 0;
};

class Halfspace : public Region {
public:
  Halfspace(std::shared_ptr<Surface> surface, int sense)
      : surface_(std::move(surface)), sense_(sense) {}

  bool contains(Position r) const override {
    double value = surface_->evaluate(r.x, r.y, r.z);
    return sense_ > 0 ? value <= 0.0 : value >= 0.0;
  }

  const Surface& surface() const { return *surface_; }
  int sense() const { return sense_; }

private:
  std::shared_ptr<Surface> surface_;
  int sense_;
};

class RegionIntersection : public Region {
public:
  explicit RegionIntersection(std::vector<std::shared_ptr<Region>> regions)
      : regions_(std::move(regions)) {}

  bool contains(Position r) const override {
    for (const auto& region : regions_) {
      if (!region->contains(r)) {
        return false;
      }
    }
    return true;
  }

private:
  std::vector<std::shared_ptr<Region>> regions_;
};

// ------------------ Cell 封装 Region 与材料索引 ------------------
class Cell {
public:
  Cell(int id, std::shared_ptr<Region> region, int material)
      : id_(id), region_(std::move(region)), material_(material) {}

  bool contains(Position r) const { return region_->contains(r); }

  int material() const { return material_; }

  // 模拟 `cell.h` 的 distance：遍历所有构成 Region 的半空间，取最近交点
  double distance(Position r, double u, double v, double w) const {
    double min_t = std::numeric_limits<double>::infinity();
    for (const auto& hs : halfspaces_) {
      double t = hs->surface().distance(r.x, r.y, r.z, u, v, w);
      if (t < min_t && hs->contains({r.x + u * 1e-10, r.y + v * 1e-10, r.z + w * 1e-10}) == false) {
        min_t = t;
      }
    }
    return min_t;
  }

  void set_halfspaces(std::vector<std::shared_ptr<Halfspace>> hs) {
    halfspaces_ = std::move(hs);
  }

private:
  int id_;
  std::shared_ptr<Region> region_;
  int material_;
  std::vector<std::shared_ptr<Halfspace>> halfspaces_;
};

// ------------------ 演示粒子导航流程 ------------------
int main() {
  auto px = std::make_shared<Plane>("px", 1.0, 0.0, 0.0, 1.0);
  auto nx = std::make_shared<Plane>("nx", -1.0, 0.0, 0.0, 1.0);
  auto py = std::make_shared<Plane>("py", 0.0, 1.0, 0.0, 1.0);
  auto ny = std::make_shared<Plane>("ny", 0.0, -1.0, 0.0, 1.0);
  auto pz = std::make_shared<Plane>("pz", 0.0, 0.0, 1.0, 1.0);
  auto nz = std::make_shared<Plane>("nz", 0.0, 0.0, -1.0, 1.0);

  auto hs_px = std::make_shared<Halfspace>(px, +1);
  auto hs_nx = std::make_shared<Halfspace>(nx, +1);
  auto hs_py = std::make_shared<Halfspace>(py, +1);
  auto hs_ny = std::make_shared<Halfspace>(ny, +1);
  auto hs_pz = std::make_shared<Halfspace>(pz, +1);
  auto hs_nz = std::make_shared<Halfspace>(nz, +1);

  auto cube_region = std::make_shared<RegionIntersection>(std::vector<std::shared_ptr<Region>>{
      hs_px, hs_nx, hs_py, hs_ny, hs_pz, hs_nz});

  Cell cube_cell(0, cube_region, /*material*/ 0);
  cube_cell.set_halfspaces({hs_px, hs_nx, hs_py, hs_ny, hs_pz, hs_nz});

  Position r{0.2, -0.1, 0.3};
  double u = 0.6, v = 0.2, w = -0.3;

  std::cout << std::boolalpha;
  std::cout << "粒子是否在 Cell 内? " << cube_cell.contains(r) << "\n";
  std::cout << "沿方向 (" << u << ", " << v << ", " << w << ") 的离开距离 = "
            << cube_cell.distance(r, u, v, w) << " cm\n";

  return 0;
}

