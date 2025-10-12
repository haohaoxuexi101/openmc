#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

// 本示例提炼 OpenMC 中“粒子状态 + 粒子银行”组合模式。核心看点：
//   - 粒子状态结构体如何保持最少数据但又便于扩展；
//   - Bank 使用 std::vector 完成 push / sample / swap，强调 STL 与物理算法的结合；
//   - 演示重采样（resample）逻辑，与 OpenMC power iteration 中的 source_bank ↔ fission_bank
//     交换完全对应。

namespace lesson_particle_bank {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;
};

struct Particle {
  Vec3 position;   // 位置向量
  Vec3 direction;  // 方向余弦
  double weight = 1.0;
  int group = 0;   // 能群编号
};

class Bank {
 public:
  void push(const Particle& p) { particles_.push_back(p); }

  [[nodiscard]] std::size_t size() const noexcept { return particles_.size(); }

  Particle& operator[](std::size_t i) { return particles_.at(i); }
  const Particle& operator[](std::size_t i) const { return particles_.at(i); }

  // 重采样接口：给定随机引擎，从当前银行中等概率抽样，生成新的银行。
  template <typename RNG>
  void resample_to(std::size_t target_size, RNG& rng) {
    if (particles_.empty()) {
      throw std::runtime_error("无法从空银行重采样");
    }
    std::uniform_int_distribution<std::size_t> dist(0, particles_.size() - 1);
    std::vector<Particle> new_bank;
    new_bank.reserve(target_size);
    for (std::size_t i = 0; i < target_size; ++i) {
      new_bank.push_back(particles_[dist(rng)]);
    }
    particles_.swap(new_bank);
  }

  auto begin() { return particles_.begin(); }
  auto end() { return particles_.end(); }
  auto begin() const { return particles_.begin(); }
  auto end() const { return particles_.end(); }

 private:
  std::vector<Particle> particles_;
};

}  // namespace lesson_particle_bank

int main() {
  using namespace lesson_particle_bank;

  // 1. 构造一个裂变银行并写入若干粒子，展示 push 接口的简单性。
  Bank fission_bank;
  for (int i = 0; i < 5; ++i) {
    Particle p;
    p.position = {0.1 * i, 0.2 * i, 0.3 * i};
    p.direction = {0.0, 0.0, 1.0};
    p.weight = 1.0;
    p.group = i % 3;
    fission_bank.push(p);
  }

  std::cout << "初始裂变银行粒子数 = " << fission_bank.size() << "\n";
  for (const auto& p : fission_bank) {
    std::cout << "  group=" << p.group << ", pos_z=" << p.position.z << "\n";
  }

  // 2. 通过 resample_to 将银行尺寸缩减为 3，模拟 k-eigenvalue 迭代中的源重采样。
  std::mt19937 rng(777u);
  fission_bank.resample_to(3, rng);

  std::cout << "重采样后粒子数 = " << fission_bank.size() << "\n";
  for (std::size_t i = 0; i < fission_bank.size(); ++i) {
    const Particle& p = fission_bank[i];
    std::cout << std::setw(2) << i << ": group=" << p.group << ", pos_z=" << p.position.z
              << ", weight=" << p.weight << "\n";
  }
}
