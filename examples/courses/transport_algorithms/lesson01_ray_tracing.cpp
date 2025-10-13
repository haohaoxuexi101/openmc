#include <iostream>
#include <vector>
#include <cmath>
#include <array>
#include <limits>
#include <algorithm>

// 课程 1：射线追踪基础
// 演示如何沿射线寻找与平面交点，模拟 OpenMC 中的几何追踪。

struct Plane {
  double nx, ny, nz; // 法向量
  double d;          // 平面方程 nx*x + ny*y + nz*z = d
};

struct Ray {
  std::array<double, 3> r; // 起点
  std::array<double, 3> u; // 单位方向
};

// 计算射线到平面的距离；若平行返回无穷大
double distance_to_plane(const Ray& ray, const Plane& plane) {
  double denom = plane.nx * ray.u[0] + plane.ny * ray.u[1] + plane.nz * ray.u[2];
  if (std::abs(denom) < 1e-12) return std::numeric_limits<double>::infinity();
  double numer = plane.d - (plane.nx * ray.r[0] + plane.ny * ray.r[1] + plane.nz * ray.r[2]);
  double t = numer / denom;
  return t > 0.0 ? t : std::numeric_limits<double>::infinity();
}

int main() {
  Plane planes[] = {{1, 0, 0, 1}, {-1, 0, 0, 1}, {0, 1, 0, 1}, {0, -1, 0, 1},
                    {0, 0, 1, 1}, {0, 0, -1, 1}};

  Ray ray{{0.0, 0.0, 0.0}, {0.6, 0.3, 0.7}};
  double norm = std::sqrt(ray.u[0]*ray.u[0] + ray.u[1]*ray.u[1] + ray.u[2]*ray.u[2]);
  for (double& c : ray.u) c /= norm;

  double min_dist = std::numeric_limits<double>::infinity();
  for (const auto& plane : planes) {
    min_dist = std::min(min_dist, distance_to_plane(ray, plane));
  }

  std::cout << "最近交距 = " << min_dist << " cm\n";
  return 0;
}
