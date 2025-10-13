#include <iomanip>
#include <iostream>
#include <limits>
#include <cmath>

// 课程 3：区域穿越与距离计算
// ----------------------------
// 该示例实现 random ray 求解器的几何核心：给定起点与方向，
// 计算射线离开当前区域的距离，并判断撞击到哪个面。
// 我们使用轴对齐盒体，可直接通过逐面参数解得到距离。

struct Vec3 {
  double x, y, z;
};

struct AxisAlignedBox {
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;
};

struct ExitInfo {
  double distance;
  int face; // 0:+x,1:-x,2:+y,3:-y,4:+z,5:-z
};

ExitInfo distance_to_exit(const AxisAlignedBox& box, Vec3 r, Vec3 u) {
  double min_t = std::numeric_limits<double>::infinity();
  int face = -1;

  auto check = [&](double plane_pos, double coord, double direction, int face_id) {
    if (std::abs(direction) < 1e-12) return;
    double t = (plane_pos - coord) / direction;
    if (t > 1e-12 && t < min_t) {
      min_t = t;
      face = face_id;
    }
  };

  check(box.xmax, r.x, u.x, 0);
  check(box.xmin, r.x, u.x, 1);
  check(box.ymax, r.y, u.y, 2);
  check(box.ymin, r.y, u.y, 3);
  check(box.zmax, r.z, u.z, 4);
  check(box.zmin, r.z, u.z, 5);

  return {min_t, face};
}

int main() {
  AxisAlignedBox region{-2.0, 2.0, -3.0, 3.0, -1.0, 1.5};
  Vec3 start{0.0, -1.0, 0.5};
  Vec3 direction{0.7, 0.5, 0.1};

  ExitInfo info = distance_to_exit(region, start, direction);
  std::cout << std::fixed << std::setprecision(4)
            << "离开距离 = " << info.distance << " cm, 撞击面编号 = " << info.face
            << '\n';

  return 0;
}

