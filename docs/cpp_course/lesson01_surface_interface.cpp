#include <cmath>
#include <iostream>
#include <string>
#include <limits>

// 本示例完全自包含，展示了 OpenMC 中曲面抽象接口的核心理念。
// - Surface 为粒子追踪提供 evaluate/distance/normal 等纯虚接口。
// - DummyHDF5Group 模拟 HDF5 序列化，证明接口不仅限于几何计算。
// - main 函数通过平面示例演示 evaluate、distance、法向量的用法。

struct Position {
  double x{};
  double y{};
  double z{};
};

struct Direction {
  double x{};
  double y{};
  double z{};
};

// 简易归一化工具，用于保持和 OpenMC 相同的几何推导前提。
Direction normalize(Direction d)
{
  const double mag = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
  return {d.x / mag, d.y / mag, d.z / mag};
}

struct DummyHDF5Group {
  void write(const std::string& key, double value) const
  {
    std::cout << "[HDF5] " << key << " = " << value << '\n';
  }
};

class Surface {
public:
  virtual ~Surface() = default;

  virtual double evaluate(const Position& r) const = 0;

  virtual double distance(const Position& r, const Direction& u) const = 0;

  virtual Direction normal(const Position& r) const = 0;

  virtual void to_hdf5(const DummyHDF5Group& grp) const = 0;

  int id() const noexcept { return id_; }
  const std::string& name() const noexcept { return name_; }

protected:
  Surface(int id, std::string name) : id_{id}, name_{std::move(name)} {}

private:
  int id_;
  std::string name_;
};

class PlaneSurface final : public Surface {
public:
  PlaneSurface(int id, std::string name, Direction normal, double d)
    : Surface{id, std::move(name)}
    , normal_{normalize(normal)}
    , offset_{d}
  {
  }

  double evaluate(const Position& r) const override
  {
    return normal_.x * r.x + normal_.y * r.y + normal_.z * r.z - offset_;
  }

  double distance(const Position& r, const Direction& u) const override
  {
    const double denom = normal_.x * u.x + normal_.y * u.y + normal_.z * u.z;
    if (std::abs(denom) < 1e-12) {
      return std::numeric_limits<double>::infinity();
    }
    const double numer = -evaluate(r);
    const double t = numer / denom;
    return t > 0.0 ? t : std::numeric_limits<double>::infinity();
  }

  Direction normal(const Position&) const override { return normal_; }

  void to_hdf5(const DummyHDF5Group& grp) const override
  {
    grp.write("offset", offset_);
  }

private:
  Direction normal_;
  double offset_;
};

int main()
{
  PlaneSurface surface{42, "fuel_top", {0.0, 0.0, 1.0}, 10.0};
  Position pos{0.0, 0.0, 7.5};
  Direction dir = normalize({0.0, 0.0, 1.0});

  std::cout << "Surface: " << surface.name() << " (id=" << surface.id() << ")\n";
  std::cout << "Signed distance (evaluate) = " << surface.evaluate(pos) << '\n';
  std::cout << "Ray distance = " << surface.distance(pos, dir) << '\n';
  Direction n = surface.normal(pos);
  std::cout << "Normal = (" << n.x << ", " << n.y << ", " << n.z << ")\n";

  DummyHDF5Group grp;
  surface.to_hdf5(grp);
}
