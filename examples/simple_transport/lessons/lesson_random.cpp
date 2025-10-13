#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// 本示例提取 OpenMC 中隶属于 random 模块的关键接口，展示如何封装随机引擎
// 与“采样策略”解耦。在 OpenMC 源码里，C++ 层提供统一的 RNG 包装，Python
// 层和物理求解均通过接口调用，这里用中文注释完整复刻其设计思路。

namespace lesson_random {

// RandomEngine 仿照 OpenMC 的包裹器，内部使用 std::mt19937 保存随机状态，
// 对外暴露“采样方法”接口：均匀数、各向同性方向、离散分布选择等。
class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed) : generator_(seed), uniform_(0.0, 1.0) {}

  // 采样 [0, 1) 均匀分布。该接口在 OpenMC 中被广泛用于自由程、反应分支选择等。
  double sample_uniform() { return uniform_(generator_); }

  // 抽样各向同性方向。这里采用余弦和方位角分解，与 OpenMC 中 direction::sample_isotropic 一致。
  std::array<double, 3> sample_isotropic_direction() {
    double mu = 2.0 * sample_uniform() - 1.0;              // cos(theta)
    double phi = 2.0 * kPi * sample_uniform();             // 方位角
    double sin_theta = std::sqrt(1.0 - mu * mu);           // 正弦项
    return {sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu};
  }

  // 离散分布抽样：OpenMC 的 tally、材料裂变谱都会使用 std::discrete_distribution。
  std::size_t sample_discrete(const std::vector<double>& weights) {
    std::discrete_distribution<std::size_t> dist(weights.begin(), weights.end());
    return dist(generator_);
  }

  // 暴露底层引擎引用，便于像 OpenMC 那样传递给 std::uniform_int_distribution。
  std::mt19937& generator() noexcept { return generator_; }

 private:
  static constexpr double kPi = 3.14159265358979323846;
  std::mt19937 generator_;
  std::uniform_real_distribution<double> uniform_;
};

}  // namespace lesson_random

int main() {
  using namespace lesson_random;

  RandomEngine rng(202406u);
  std::cout << std::fixed << std::setprecision(6);

  // 1. 连续均匀分布：演示“接口层与算法层分离”。
  std::cout << "均匀随机数示例：";
  for (int i = 0; i < 5; ++i) {
    std::cout << ' ' << rng.sample_uniform();
  }
  std::cout << "\n\n";

  // 2. 各向同性方向抽样：输出方向余弦，说明接口直接返回值对象，调用方无需关心内部状态。
  auto direction = rng.sample_isotropic_direction();
  std::cout << "各向同性方向抽样 (ux, uy, uz) = (" << direction[0] << ", " << direction[1]
            << ", " << direction[2] << ")\n\n";

  // 3. 离散分布：模拟裂变谱或散射概率的选择。
  std::vector<double> chi{0.7, 0.2, 0.1};
  std::vector<int> histogram(chi.size(), 0);
  for (int n = 0; n < 1000; ++n) {
    std::size_t g = rng.sample_discrete(chi);
    histogram[g]++;
  }
  std::cout << "离散抽样直方图：";
  for (std::size_t g = 0; g < histogram.size(); ++g) {
    std::cout << " 群" << g << '=' << histogram[g];
  }
  std::cout << '\n';
}
