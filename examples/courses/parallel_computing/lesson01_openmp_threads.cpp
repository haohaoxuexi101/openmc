#include <iostream>
#include <omp.h>
#include <vector>

// 课程 1：OpenMP 粒子并行
// 展示如何使用 OpenMP 将粒子历史切分为线程并行处理。

int main() {
  std::vector<int> histories(8);
  for (size_t i = 0; i < histories.size(); ++i) histories[i] = static_cast<int>(i);

  #pragma omp parallel for schedule(static)
  for (int i = 0; i < static_cast<int>(histories.size()); ++i) {
    int tid = omp_get_thread_num();
    std::cout << "线程 " << tid << " 处理历史 " << histories[i] << "\n";
  }
  return 0;
}
