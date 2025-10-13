#include <iostream>
#include <omp.h>
#include <random>

// 课程 3：线程私有随机数流
// 展示如何确保每个 OpenMP 线程拥有独立 RNG，避免相关性。

int main() {
  int histories = 8;

  #pragma omp parallel for
  for (int i = 0; i < histories; ++i) {
    int tid = omp_get_thread_num();
    std::seed_seq seed{tid, i};
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    double xi = uni(rng);
    #pragma omp critical
    {
      std::cout << "线程 " << tid << " 历史 " << i << " 随机数 = " << xi << "\n";
    }
  }
  return 0;
}
