#include <iostream>
#include <vector>
#include <numeric>

// 课程 5：负载划分策略
// 演示如何依据历史数量与权重划分任务，模拟 Bank 分发。

std::vector<int> partition(int tasks, int ranks) {
  std::vector<int> counts(ranks, tasks / ranks);
  int remainder = tasks % ranks;
  for (int i = 0; i < remainder; ++i) counts[i] += 1;
  return counts;
}

int main() {
  int histories = 17;
  int ranks = 4;
  auto counts = partition(histories, ranks);

  int offset = 0;
  for (int r = 0; r < ranks; ++r) {
    std::cout << "进程 " << r << " 处理 " << counts[r] << " 条历史, 起始编号 "
              << offset << "\n";
    offset += counts[r];
  }
  return 0;
}
