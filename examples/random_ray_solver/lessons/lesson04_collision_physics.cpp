#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <tuple>
#include <vector>

// 课程 4：多群碰撞抽样
// ----------------------
// 本课演示 random ray 求解器如何根据多群截面抽样碰撞类型：
// 1. 依据 Σ_t 抽取自由程；
// 2. 使用 Σ_s、νΣ_f、Σ_a = Σ_t - Σ_s - νΣ_f 计算离散概率；
// 3. 输出抽中的事件及其对统计量的贡献。

struct GroupXS {
  double sigma_t;
  double sigma_s;
  double nu_sigma_f;
};

struct Material {
  std::string name;
  std::vector<GroupXS> groups;
};

struct CollisionResult {
  std::string event;
  double score;
};

CollisionResult sample_collision(const Material& mat, int group, std::mt19937& rng) {
  const GroupXS& xs = mat.groups[group];
  double absorption = std::max(0.0, xs.sigma_t - xs.sigma_s - xs.nu_sigma_f);

  std::vector<std::pair<std::string, double>> channels = {
      {"散射", xs.sigma_s},
      {"裂变", xs.nu_sigma_f},
      {"吸收", absorption},
  };

  double total = 0.0;
  for (const auto& ch : channels) total += ch.second;

  std::uniform_real_distribution<double> dist(0.0, total);
  double r = dist(rng);

  std::string chosen;
  double prefix = 0.0;
  for (const auto& ch : channels) {
    prefix += ch.second;
    if (r <= prefix) {
      chosen = ch.first;
      break;
    }
  }

  double contribution = (chosen == "散射") ? xs.sigma_s / xs.sigma_t : xs.nu_sigma_f / xs.sigma_t;
  return {chosen, contribution};
}

int main() {
  Material fuel{"燃料", {{0.6, 0.25, 0.3}, {0.5, 0.23, 0.22}, {0.4, 0.2, 0.15}}};
  std::mt19937 rng(20240613u);

  const int group = 0; // 快群
  std::cout << "材料: " << fuel.name << ", 群: " << group << '\n';
  for (int i = 0; i < 5; ++i) {
    CollisionResult res = sample_collision(fuel, group, rng);
    std::cout << "历史 " << i << " -> " << res.event << ", 贡献因子 = " << res.score << '\n';
  }

  return 0;
}

