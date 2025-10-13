#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// 课程 2：随机射线与各向同性方向抽样
// ------------------------------------
// 本课展示如何从三维体源中均匀采样射线起点，并生成各向同性方向。
// random ray 求解器通过这种方式替代逐粒子跟踪，使得一个射线段即可代表
// 多个体元的平均通量贡献。

struct Vec3 {
  double x, y, z;
};

struct AxisAlignedBox {
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;
};

Vec3 sample_uniform_point(const AxisAlignedBox& box, std::mt19937& rng) {
  std::uniform_real_distribution<double> dist_x(box.xmin, box.xmax);
  std::uniform_real_distribution<double> dist_y(box.ymin, box.ymax);
  std::uniform_real_distribution<double> dist_z(box.zmin, box.zmax);
  return {dist_x(rng), dist_y(rng), dist_z(rng)};
}

Vec3 sample_isotropic_direction(std::mt19937& rng) {
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  double mu = 2.0 * dist01(rng) - 1.0;
  double phi = 2.0 * M_PI * dist01(rng);
  double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
}

int main() {
  AxisAlignedBox source_box{-2.0, 2.0, -2.0, 2.0, -2.0, 2.0};
  std::mt19937 rng(20240613u);

  const int samples = 5;
  std::cout << "随机射线样本 (起点 + 方向)" << '\n';
  for (int i = 0; i < samples; ++i) {
    Vec3 r = sample_uniform_point(source_box, rng);
    Vec3 u = sample_isotropic_direction(rng);
    std::cout << std::fixed << std::setprecision(3)
              << "起点(" << r.x << ", " << r.y << ", " << r.z << ")"
              << " 方向(" << u.x << ", " << u.y << ", " << u.z << ")" << '\n';
  }

  return 0;
}

