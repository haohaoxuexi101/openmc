#include <cstdint>
#include <iomanip>
#include <iostream>
#include <random>

// RNGStream 模拟 OpenMC 对随机数子流的管理：
// - 每个流保存主种子与步进间隔，实现线程安全的分段使用。
// - spawn_substream 通过复制状态 + 偏移生成子流，避免共享状态。
// - main 展示主流与子流的采样互不干扰。

class RNGStream {
public:
  RNGStream(std::uint64_t seed, std::uint64_t stride)
    : engine_{seed}
    , stride_{stride}
  {
  }

  double uniform()
  {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(engine_);
  }

  RNGStream spawn_substream(std::uint64_t skip) const
  {
    RNGStream copy = *this;
    copy.advance(skip * copy.stride_);
    return copy;
  }

  void advance(std::uint64_t steps) const
  {
    for (std::uint64_t i = 0; i < steps; ++i) {
      engine_.discard(1);
    }
  }

private:
  mutable std::mt19937_64 engine_;
  std::uint64_t stride_;
};

int main()
{
  RNGStream master{12345, 4};
  RNGStream worker_a = master.spawn_substream(0);
  RNGStream worker_b = master.spawn_substream(1);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Master draws: " << master.uniform() << ", " << master.uniform() << '\n';
  std::cout << "Worker A draws: " << worker_a.uniform() << ", " << worker_a.uniform() << '\n';
  std::cout << "Worker B draws: " << worker_b.uniform() << ", " << worker_b.uniform() << '\n';
}
