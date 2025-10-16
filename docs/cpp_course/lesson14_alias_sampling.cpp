#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

// AliasTable 体现“预处理 + O(1) 采样”的算法思想：
// - build 函数构造概率与别名列表，仿照 OpenMC 的源项能量采样。
// - sample 使用一个随机源对象，依赖组合而非继承，便于替换。
// - main 对比目标分布与采样频率，输出直观统计。

class AliasTable {
public:
  explicit AliasTable(std::vector<double> probabilities)
  {
    const std::size_t n = probabilities.size();
    table_.resize(n);

    std::vector<double> scaled(n);
    std::vector<std::size_t> small;
    std::vector<std::size_t> large;

    for (std::size_t i = 0; i < n; ++i) {
      scaled[i] = probabilities[i] * n;
      if (scaled[i] < 1.0) {
        small.push_back(i);
      } else {
        large.push_back(i);
      }
    }

    while (!small.empty() && !large.empty()) {
      std::size_t s = small.back();
      small.pop_back();
      std::size_t l = large.back();

      table_[s] = {scaled[s], l};
      scaled[l] = (scaled[l] + scaled[s]) - 1.0;

      if (scaled[l] < 1.0) {
        small.push_back(l);
        large.pop_back();
      }
    }

    for (std::size_t g : large) {
      table_[g] = {1.0, g};
    }
    for (std::size_t g : small) {
      table_[g] = {1.0, g};
    }
  }

  std::size_t sample(std::mt19937& rng) const
  {
    std::uniform_int_distribution<std::size_t> pick(0, table_.size() - 1);
    std::uniform_real_distribution<double> split(0.0, 1.0);

    std::size_t column = pick(rng);
    const Entry& entry = table_[column];
    double draw = split(rng);
    return draw < entry.prob ? column : entry.alias;
  }

private:
  struct Entry {
    double prob;
    std::size_t alias;
  };
  std::vector<Entry> table_;
};

int main()
{
  std::vector<double> distribution{0.05, 0.15, 0.5, 0.2, 0.1};
  AliasTable table{distribution};

  std::mt19937 rng{42};
  std::vector<int> counts(distribution.size());
  const int samples = 10000;
  for (int i = 0; i < samples; ++i) {
    ++counts[table.sample(rng)];
  }

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Index  target  sampled\n";
  for (std::size_t i = 0; i < distribution.size(); ++i) {
    double target = distribution[i];
    double freq = static_cast<double>(counts[i]) / samples;
    std::cout << std::setw(5) << i << std::setw(8) << target << std::setw(9) << freq << '\n';
  }
}
