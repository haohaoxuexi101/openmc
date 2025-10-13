#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// 课程 6：随机射线求解器主循环
// --------------------------------
// 汇总前五节课的成果，构建一个最小 random ray 多群求解循环：
// - 构造材料与区域；
// - 随机抽样射线并计算在每个区域内的轨迹长度；
// - 依据 Σ_t 计算群内无碰撞存活概率，对曲线长度计分；
// - 累计裂变贡献并输出近似的 k-effective 估计。

struct Vec3 {
  double x, y, z;
};

struct AxisAlignedBox {
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;

  bool contains(Vec3 r) const {
    return r.x >= xmin && r.x <= xmax && r.y >= ymin && r.y <= ymax && r.z >= zmin && r.z <= zmax;
  }
};

struct Material {
  std::string name;
  std::vector<double> sigma_t;
  std::vector<double> nu_sigma_f;
};

struct Region {
  AxisAlignedBox box;
  int material_id;
};

struct Tallies {
  std::vector<double> track_length;
  std::vector<double> fission_source;
};

Vec3 sample_point(const AxisAlignedBox& box, std::mt19937& rng) {
  std::uniform_real_distribution<double> dx(box.xmin, box.xmax);
  std::uniform_real_distribution<double> dy(box.ymin, box.ymax);
  std::uniform_real_distribution<double> dz(box.zmin, box.zmax);
  return {dx(rng), dy(rng), dz(rng)};
}

Vec3 sample_direction(std::mt19937& rng) {
  std::uniform_real_distribution<double> dist01(0.0, 1.0);
  double mu = 2.0 * dist01(rng) - 1.0;
  double phi = 2.0 * M_PI * dist01(rng);
  double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
}

int locate_region(const std::vector<Region>& regions, Vec3 r) {
  for (size_t i = 0; i < regions.size(); ++i) {
    if (regions[i].box.contains(r)) return static_cast<int>(i);
  }
  return -1;
}

double distance_to_exit(const AxisAlignedBox& box, Vec3 r, Vec3 u) {
  double min_t = 1e30;
  auto check = [&](double plane, double coord, double dir) {
    if (std::abs(dir) < 1e-12) return;
    double t = (plane - coord) / dir;
    if (t > 1e-12) min_t = std::min(min_t, t);
  };
  check(box.xmax, r.x, u.x);
  check(box.xmin, r.x, u.x);
  check(box.ymax, r.y, u.y);
  check(box.ymin, r.y, u.y);
  check(box.zmax, r.z, u.z);
  check(box.zmin, r.z, u.z);
  return min_t;
}

int main() {
  std::vector<Material> materials = {
      {"燃料", {0.6, 0.5, 0.4}, {0.3, 0.25, 0.15}},
      {"慢化剂", {0.2, 0.15, 0.1}, {0.0, 0.0, 0.0}},
  };

  std::vector<Region> regions = {
      {{-3, 3, -3, 3, -3, 3}, 0},
      {{-5, 5, -5, 5, -5, 5}, 1},
  };

  Tallies tallies;
  tallies.track_length.assign(3, 0.0);
  tallies.fission_source.assign(3, 0.0);

  AxisAlignedBox source_box = regions[0].box;
  std::mt19937 rng(20240613u);

  const int histories = 1000;
  for (int h = 0; h < histories; ++h) {
    Vec3 r = sample_point(source_box, rng);
    Vec3 u = sample_direction(rng);
    int group = h % 3;
    double weight = 1.0;

    int region_id = locate_region(regions, r);
    if (region_id < 0) continue;

    const Region* region = &regions[region_id];
    const Material& mat = materials[region->material_id];

    double distance = distance_to_exit(region->box, r, u);
    double attenuation = std::exp(-mat.sigma_t[group] * distance);
    tallies.track_length[group] += weight * distance;
    tallies.fission_source[group] += weight * mat.nu_sigma_f[group] * attenuation;
  }

  double total_track = tallies.track_length[0] + tallies.track_length[1] + tallies.track_length[2];
  double total_source = tallies.fission_source[0] + tallies.fission_source[1] + tallies.fission_source[2];
  double k_estimate = (total_source > 0.0) ? total_source / static_cast<double>(histories) : 0.0;

  std::cout << std::fixed << std::setprecision(4);
  for (int g = 0; g < 3; ++g) {
    std::cout << "群" << g << " track-length = " << tallies.track_length[g]
              << ", 裂变源贡献 = " << tallies.fission_source[g] << '\n';
  }
  std::cout << "总轨迹长度 = " << total_track << '\n';
  std::cout << "k 有效近似 = " << k_estimate << '\n';

  return 0;
}

