#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

// 本节展示三个典型派生曲面：平面、球面、圆柱。所有逻辑均可编译运行，
// 并输出粒子与曲面的交点距离，帮助理解 OpenMC 几何追踪的核心计算。

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

Direction normalize(Direction d)
{
  const double mag = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
  return {d.x / mag, d.y / mag, d.z / mag};
}

class Surface {
public:
  virtual ~Surface() = default;
  virtual double distance(const Position& r, const Direction& u) const = 0;
  virtual std::string description() const = 0;
};

class PlaneSurface final : public Surface {
public:
  PlaneSurface(Direction normal, double offset) : normal_{normalize(normal)}, offset_{offset} {}

  double distance(const Position& r, const Direction& u) const override
  {
    const double denom = normal_.x * u.x + normal_.y * u.y + normal_.z * u.z;
    if (std::abs(denom) < 1e-12) {
      return std::numeric_limits<double>::infinity();
    }
    const double numer = offset_ - (normal_.x * r.x + normal_.y * r.y + normal_.z * r.z);
    const double t = numer / denom;
    return t > 0.0 ? t : std::numeric_limits<double>::infinity();
  }

  std::string description() const override
  {
    return "Plane(z=" + std::to_string(offset_) + ")";
  }

private:
  Direction normal_;
  double offset_;
};

class SphereSurface final : public Surface {
public:
  SphereSurface(Position center, double radius) : center_{center}, radius_{radius} {}

  double distance(const Position& r, const Direction& u) const override
  {
    const double dx = r.x - center_.x;
    const double dy = r.y - center_.y;
    const double dz = r.z - center_.z;
    const double b = 2.0 * (dx * u.x + dy * u.y + dz * u.z);
    const double c = dx * dx + dy * dy + dz * dz - radius_ * radius_;
    const double discriminant = b * b - 4 * c;
    if (discriminant < 0) {
      return std::numeric_limits<double>::infinity();
    }
    const double sqrt_disc = std::sqrt(discriminant);
    const double t1 = (-b - sqrt_disc) / 2.0;
    const double t2 = (-b + sqrt_disc) / 2.0;
    if (t1 > 0.0) return t1;
    if (t2 > 0.0) return t2;
    return std::numeric_limits<double>::infinity();
  }

  std::string description() const override
  {
    return "Sphere(r=" + std::to_string(radius_) + ")";
  }

private:
  Position center_;
  double radius_;
};

class CylinderSurface final : public Surface {
public:
  CylinderSurface(double radius) : radius_{radius} {}

  double distance(const Position& r, const Direction& u) const override
  {
    // 仅考虑 z 轴圆柱，和 OpenMC 默认方向保持一致。
    const double a = u.x * u.x + u.y * u.y;
    const double b = 2.0 * (r.x * u.x + r.y * u.y);
    const double c = r.x * r.x + r.y * r.y - radius_ * radius_;
    if (std::abs(a) < 1e-12) {
      return std::numeric_limits<double>::infinity();
    }
    const double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
      return std::numeric_limits<double>::infinity();
    }
    const double sqrt_disc = std::sqrt(discriminant);
    const double t1 = (-b - sqrt_disc) / (2.0 * a);
    const double t2 = (-b + sqrt_disc) / (2.0 * a);
    if (t1 > 0.0) return t1;
    if (t2 > 0.0) return t2;
    return std::numeric_limits<double>::infinity();
  }

  std::string description() const override
  {
    return "Cylinder(r=" + std::to_string(radius_) + ")";
  }

private:
  double radius_;
};

int main()
{
  Position start{0.1, 0.1, -2.0};
  Direction dir = normalize({0.0, 0.0, 1.0});

  std::vector<std::unique_ptr<Surface>> surfaces;
  surfaces.emplace_back(std::make_unique<PlaneSurface>(Direction{0.0, 0.0, 1.0}, 0.0));
  surfaces.emplace_back(std::make_unique<SphereSurface>(Position{0.0, 0.0, 4.0}, 2.0));
  surfaces.emplace_back(std::make_unique<CylinderSurface>(0.5));

  std::cout << "Starting from (" << start.x << ", " << start.y << ", " << start.z << ")\n";
  for (const auto& surface : surfaces) {
    std::cout << "Distance to " << surface->description() << " = "
              << surface->distance(start, dir) << '\n';
  }
}
