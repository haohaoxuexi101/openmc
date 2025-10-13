#include <omp.h>
#include <chrono>
#include <iostream>
#include <random>
#include <thread>
#include <vector>

// 课程 7：OpenMP 动态任务 + 裂变银行工作窃取
// ------------------------------------------
// 本课模拟代际求解中不同粒子历史的计算时间差异，演示如何使用 OpenMP 任务与
// 工作窃取策略平衡负载。我们创建一个粒子队列，每个粒子拥有不同的“追踪时长”。

struct ParticleWork {
  int id;
  int steps; // 模拟该粒子需要的追踪步数
};

int main() {
  std::vector<ParticleWork> bank;
  for (int i = 0; i < 16; ++i) {
    bank.push_back({i, 5 + (i * 3) % 11});
  }

  std::cout << "OpenMP 任务分发示例\n";

#pragma omp parallel
  {
#pragma omp single
    {
      for (const auto& work : bank) {
#pragma omp task firstprivate(work)
        {
          // 通过 sleep 模拟不同粒子的追踪耗时
          std::this_thread::sleep_for(std::chrono::milliseconds(5 * work.steps));
#pragma omp critical
          {
            std::cout << "线程 " << omp_get_thread_num() << " 完成粒子 " << work.id
                      << " (步数=" << work.steps << ")\n";
          }
        }
      }
    }
  }

  std::cout << "全部粒子完成。\n";
  return 0;
}

