#include <iostream>
#include <random>
#include <vector>
#include <string>

// 课程 2：碰撞类型抽样
// 模拟 OpenMC 中根据反应截面随机选择散射/吸收/裂变的逻辑。

struct ReactionChannel {
  std::string name;
  double xs; // 截面
};

int main() {
  std::vector<ReactionChannel> channels = {{"散射", 0.8}, {"吸收", 0.1}, {"裂变", 0.2}};
  double total = 0.0;
  for (auto& ch : channels) total += ch.xs;

  std::mt19937 rng(2024);
  std::uniform_real_distribution<double> uni(0.0, total);

  for (int i = 0; i < 5; ++i) {
    double xi = uni(rng);
    double accum = 0.0;
    for (const auto& ch : channels) {
      accum += ch.xs;
      if (xi <= accum) {
        std::cout << "样本 " << i << " -> " << ch.name << "\n";
        break;
      }
    }
  }
  return 0;
}
