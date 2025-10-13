#include <iostream>
#include <array>
#include <string>
#include <random>
#include <cmath>

// 课程 2：粒子对象的状态生命周期
// 借鉴说明：核心接口与字段命名参考 OpenMC `src/particle_data.h` 内 Particle
// 结构体的状态字段，特别是位置 (r)、方向 (u)、能群 (g) 以及 alive 标志。
// 示例在精简实现的同时保留 move / scatter / absorb 的接口节奏，方便读者
// 对照源码理解状态迁移的边界控制。

class Particle {
public:
  Particle(std::string tag, std::array<double, 3> r, std::array<double, 3> u,
           int g)
      : tag_(std::move(tag)), r_(r), u_(u), group_(g), alive_(true) {}

  void move(double distance) {
    // 位移接口：限定只有活粒子才能移动
    if (!alive_) return;
    for (int i = 0; i < 3; ++i) {
      r_[i] += distance * u_[i];
    }
  }

  void scatter(const std::array<double, 3>& new_dir, int new_group) {
    // 散射接口：改变方向与能群
    if (!alive_) return;
    u_ = new_dir;
    group_ = new_group;
  }

  void absorb() { alive_ = false; }

  bool alive() const { return alive_; }
  const std::string& tag() const { return tag_; }
  int group() const { return group_; }

  void print_state() const {
    std::cout << "粒子" << tag_ << " 位置:(" << r_[0] << ", " << r_[1] << ", "
              << r_[2] << ") 能群:" << group_ << " 状态:" << (alive_ ? "存活" : "死亡")
              << "\n";
  }

private:
  std::string tag_;
  std::array<double, 3> r_;
  std::array<double, 3> u_;
  int group_;
  bool alive_;
};

int main() {
  // 随机方向演示，模拟 OpenMC 中使用随机引擎更新粒子状态
  std::mt19937 rng(42);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);

  Particle neutron("n0", {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}, 1);
  neutron.print_state();

  // 1) 传播 1 cm
  neutron.move(1.0);
  neutron.print_state();

  // 2) 发生散射，随机生成新方向与能群
  std::array<double, 3> new_dir{uniform(rng), uniform(rng), uniform(rng)};
  double norm = std::sqrt(new_dir[0] * new_dir[0] + new_dir[1] * new_dir[1] +
                          new_dir[2] * new_dir[2]);
  for (double& c : new_dir) c /= norm;
  neutron.scatter(new_dir, 2);
  neutron.print_state();

  // 3) 被吸收
  neutron.absorb();
  neutron.print_state();

  return 0;
}
