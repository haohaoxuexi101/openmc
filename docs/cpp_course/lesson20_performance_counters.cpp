#include <chrono>
#include <iostream>
#include <map>
#include <string>
#include <thread>

// Lesson 20 总结性能度量：ScopedTimer 以 RAII 自动累计阶段耗时，
// 模拟 OpenMC 对热点函数的统计；main 通过三次循环输出时间总和。

class ScopedTimer {
public:
  ScopedTimer(std::string label, std::map<std::string, double>& accum)
    : label_{std::move(label)}
    , accum_{accum}
    , start_{std::chrono::high_resolution_clock::now()}
  {
  }

  ~ScopedTimer()
  {
    auto end = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(end - start_).count();
    accum_[label_] += elapsed;
  }

private:
  std::string label_;
  std::map<std::string, double>& accum_;
  std::chrono::high_resolution_clock::time_point start_;
};

void simulate_material_lookup()
{
  std::this_thread::sleep_for(std::chrono::milliseconds(10));
}

void simulate_collision()
{
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
}

int main()
{
  std::map<std::string, double> timers;

  for (int i = 0; i < 3; ++i) {
    {
      ScopedTimer timer{"material_lookup", timers};
      simulate_material_lookup();
    }
    {
      ScopedTimer timer{"collision", timers};
      simulate_collision();
    }
  }

  for (const auto& [label, seconds] : timers) {
    std::cout << label << ": " << seconds << " s" << '\n';
  }
}
