#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// 本课程补充 OpenMC 中“概率表 (Probability Table)”的数据结构与抽样流程。
// 在未解析共振区 (URR) 内，OpenMC 读取 ACE/HEATED_PROBABILITY_TABLE，
// 根据积分区间采样总截面、散射截面、裂变截面等。该机制是保证共振自遮蔽的关键。
// 本示例通过简化的五点概率表演示：
//   - ProbabilityTable 存储 P_i, Σ_t,i, Σ_s,i, νΣ_f,i；
//   - sample(double xi) 根据随机数返回一个条目；
//   - collision(double xi_reaction) 根据截面比例选择反应类型。

namespace lesson_probability {

struct ProbabilityTableEntry {
  double probability{};  // 累积概率间隔宽度
  double sigma_t{};
  double sigma_s{};
  double nu_sigma_f{};
};

class ProbabilityTable {
 public:
  explicit ProbabilityTable(std::vector<ProbabilityTableEntry> entries)
      : entries_(std::move(entries)) {
    normalize();
  }

  const ProbabilityTableEntry& sample(double xi) const {
    double cumulative = 0.0;
    for (const auto& entry : entries_) {
      cumulative += entry.probability;
      if (xi < cumulative) {
        return entry;
      }
    }
    return entries_.back();
  }

  std::string collision(double xi_reaction, const ProbabilityTableEntry& entry) const {
    double absorb = entry.sigma_t - entry.sigma_s - entry.nu_sigma_f;
    double cumulative = absorb;
    if (xi_reaction < cumulative) {
      return "absorption";
    }
    cumulative += entry.sigma_s;
    if (xi_reaction < cumulative) {
      return "scatter";
    }
    return "fission";
  }

 private:
  void normalize() {
    double sum = 0.0;
    for (const auto& entry : entries_) {
      sum += entry.probability;
    }
    for (auto& entry : entries_) {
      entry.probability /= sum;
    }
  }

  std::vector<ProbabilityTableEntry> entries_;
};

}  // namespace lesson_probability

int main() {
  using namespace lesson_probability;

  ProbabilityTable table({{0.1, 6.0, 3.0, 1.0},
                          {0.2, 5.0, 2.5, 1.2},
                          {0.3, 4.0, 2.0, 1.0},
                          {0.25, 3.5, 1.8, 0.9},
                          {0.15, 3.0, 1.5, 0.7}});

  std::mt19937 rng{777u};
  std::uniform_real_distribution<double> uniform(0.0, 1.0);

  for (int i = 0; i < 5; ++i) {
    double xi = uniform(rng);
    const auto& entry = table.sample(xi);
    double xi_reaction = uniform(rng) * entry.sigma_t;
    std::string reaction = table.collision(xi_reaction, entry);
    std::cout << "样本 " << i << ": Σ_t = " << std::setw(6) << entry.sigma_t
              << ", Σ_s = " << std::setw(6) << entry.sigma_s << ", νΣ_f = " << std::setw(6)
              << entry.nu_sigma_f << ", 反应 = " << reaction << "\n";
  }
}
