#include <algorithm>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>

// BankCompactor 模拟 OpenMC 次级粒子库的“优先级合并”算法：
// - Secondary 粒子记录能量与重要性；
// - Compactor 使用最小堆保持固定容量，体现贪心策略；
// - main 显示被保留的粒子，解释负载均衡的意义。

struct Secondary {
  double energy_eV;
  double importance;
};

struct ImportanceCompare {
  bool operator()(const Secondary& a, const Secondary& b) const
  {
    return a.importance > b.importance; // priority_queue 默认大顶堆
  }
};

class BankCompactor {
public:
  explicit BankCompactor(std::size_t capacity) : capacity_{capacity} {}

  void push(Secondary s)
  {
    if (queue_.size() < capacity_) {
      queue_.push(std::move(s));
    } else if (s.importance > queue_.top().importance) {
      queue_.pop();
      queue_.push(std::move(s));
    }
  }

  std::vector<Secondary> extract_sorted()
  {
    std::vector<Secondary> result;
    while (!queue_.empty()) {
      result.push_back(queue_.top());
      queue_.pop();
    }
    std::sort(result.begin(), result.end(), [](const Secondary& a, const Secondary& b) {
      return a.importance > b.importance;
    });
    return result;
  }

private:
  std::size_t capacity_;
  std::priority_queue<Secondary, std::vector<Secondary>, ImportanceCompare> queue_;
};

int main()
{
  BankCompactor compactor{3};
  std::vector<Secondary> candidates{
    {2e5, 0.2},
    {1e5, 0.4},
    {5e5, 0.1},
    {3e5, 0.8},
    {4e5, 0.6},
  };

  for (const auto& s : candidates) {
    compactor.push(s);
  }

  std::cout << "Retained particles (importance desc):\n";
  for (const auto& s : compactor.extract_sorted()) {
    std::cout << "E = " << s.energy_eV << " eV, importance = " << s.importance << '\n';
  }
}
