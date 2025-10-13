#include <algorithm>
#include <deque>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// 课程 5：裂变银行与源项再生
// ----------------------------
// random ray 求解器在处理裂变时需要把次级中子放入“银行”中。
// 本课展示如何维护一个队列，将同一批历史产生的裂变粒子收集、归一化，
// 并为下一代求解生成新的源项。

struct BankEntry {
  double x, y, z;
  double energy_group;
  double weight;
};

int main() {
  std::deque<BankEntry> fission_bank;
  std::mt19937 rng(20240613u);
  std::uniform_real_distribution<double> dist_pos(-2.0, 2.0);
  std::discrete_distribution<int> dist_group({3, 4, 5});

  // 模拟 5 次裂变事件，每次可能产生 2~4 个次级粒子
  for (int history = 0; history < 5; ++history) {
    int produced = 2 + history % 3;
    for (int n = 0; n < produced; ++n) {
      fission_bank.push_back({dist_pos(rng), dist_pos(rng), dist_pos(rng),
                              static_cast<double>(dist_group(rng)), 1.0});
    }
  }

  std::cout << "裂变银行粒子总数 = " << fission_bank.size() << '\n';

  // 归一化：每个粒子权重除以总数
  double total = static_cast<double>(fission_bank.size());
  for (auto& entry : fission_bank) {
    entry.weight /= total;
  }

  // 生成下一代源项
  std::vector<BankEntry> next_generation;
  next_generation.insert(next_generation.end(), fission_bank.begin(), fission_bank.end());

  std::cout << "下一代源项示例 (前 3 个粒子)" << '\n';
  for (size_t i = 0; i < std::min<size_t>(3, next_generation.size()); ++i) {
    const auto& p = next_generation[i];
    std::cout << std::fixed << std::setprecision(3)
              << "r=(" << p.x << ", " << p.y << ", " << p.z << "), g=" << p.energy_group
              << ", w=" << p.weight << '\n';
  }

  return 0;
}

