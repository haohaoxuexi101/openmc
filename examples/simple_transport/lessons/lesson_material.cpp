#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// 本示例拆解 OpenMC 中 Material 类的核心职责：
//   1. 存储多群宏观截面；
//   2. 提供派生属性（总截面）；
//   3. 将物理意义映射到面向对象接口。
// 通过中文注释解释每个函数的设计哲学，并附带一次 Monte Carlo 采样演示。

namespace lesson_material {

class Material {
 public:
  Material(std::string name, std::vector<double> sigma_a, std::vector<double> sigma_f,
           std::vector<double> nu, std::vector<std::vector<double>> sigma_s,
           std::vector<double> chi)
      : name_(std::move(name)), sigma_a_(std::move(sigma_a)), sigma_f_(std::move(sigma_f)),
        nu_(std::move(nu)), sigma_s_(std::move(sigma_s)), chi_(std::move(chi)) {
    // 防御式编程：在构造阶段即验证数据一致性，避免运行期才发现问题。
    std::size_t groups = sigma_a_.size();
    if (sigma_f_.size() != groups || nu_.size() != groups || sigma_s_.size() != groups ||
        chi_.size() != groups) {
      throw std::runtime_error("多群截面维度不一致");
    }
    for (const auto& row : sigma_s_) {
      if (row.size() != groups) {
        throw std::runtime_error("散射矩阵非方阵");
      }
    }
  }

  [[nodiscard]] const std::string& name() const noexcept { return name_; }
  [[nodiscard]] std::size_t groups() const noexcept { return sigma_a_.size(); }
  [[nodiscard]] double sigma_a(std::size_t g) const { return sigma_a_.at(g); }
  [[nodiscard]] double sigma_f(std::size_t g) const { return sigma_f_.at(g); }
  [[nodiscard]] double nu(std::size_t g) const { return nu_.at(g); }
  [[nodiscard]] const std::vector<double>& chi() const noexcept { return chi_; }

  // σ_t = σ_a + σ_f + Σ σ_s，直接在成员函数中实现，保证外部使用时不必重复手工求和。
  [[nodiscard]] double sigma_t(std::size_t g) const {
    double scatter = std::accumulate(sigma_s_[g].begin(), sigma_s_[g].end(), 0.0);
    return sigma_a_[g] + sigma_f_[g] + scatter;
  }

  [[nodiscard]] const std::vector<double>& scatter_row(std::size_t g) const {
    return sigma_s_.at(g);
  }

 private:
  std::string name_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> sigma_s_;
  std::vector<double> chi_;
};

// 该结构体模拟 OpenMC 中“碰撞结果”，方便将 Material 的接口映射到物理分支决策。
struct CollisionResult {
  std::string channel;  // 反应类型标签：scatter / absorption / fission
  std::size_t outgoing_group = 0;  // 若发生散射或裂变，给出新能群编号
};

// 随机引擎类模仿 OpenMC 的 Random 类封装，演示组合的力量。
class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed = 42u) : rng_(seed) {}

  double uniform() { return uniform_(rng_); }

  std::size_t pick_group(const std::vector<double>& pdf) {
    std::discrete_distribution<std::size_t> dist(pdf.begin(), pdf.end());
    return dist(rng_);
  }

 private:
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_{0.0, 1.0};
};

// 利用 Material 提供的接口执行一次碰撞抽样，完整模拟吸收/散射/裂变流程。
inline CollisionResult sample_collision(const Material& mat, std::size_t incoming_group,
                                        RandomEngine& rng) {
  double sigma_t = mat.sigma_t(incoming_group);
  double xi = rng.uniform() * sigma_t;
  double cumulative = mat.sigma_a(incoming_group);
  if (xi < cumulative) {
    return {"absorption", incoming_group};
  }
  cumulative += mat.sigma_f(incoming_group);
  if (xi < cumulative) {
    // 裂变产生的能群由裂变中子能谱 chi 控制。
    return {"fission", rng.pick_group(mat.chi())};
  }
  // 剩余概率用于散射，根据散射矩阵在当前行中的分布抽样去向能群。
  const auto& row = mat.scatter_row(incoming_group);
  return {"scatter", rng.pick_group(row)};
}

}  // namespace lesson_material

int main() {
  using lesson_material::CollisionResult;
  using lesson_material::Material;
  using lesson_material::RandomEngine;
  using lesson_material::sample_collision;

  // 构造一个三群材料，数据取自简化版燃料，展示多群接口完整性。
  Material fuel("fuel", {0.01, 0.02, 0.03}, {0.04, 0.03, 0.02}, {2.5, 2.5, 2.5},
                {{0.02, 0.00, 0.00}, {0.01, 0.03, 0.00}, {0.00, 0.02, 0.04}}, {0.6, 0.3, 0.1});

  std::cout << "材料名称: " << fuel.name() << "\n";
  for (std::size_t g = 0; g < fuel.groups(); ++g) {
    std::cout << "  群 " << g << ": σ_t = " << fuel.sigma_t(g)
              << ", νσ_f = " << fuel.nu(g) * fuel.sigma_f(g) << "\n";
  }

  // 执行 5 次碰撞抽样，输出每次反应类型与能群变化，帮助读者把握接口节奏。
  RandomEngine rng(2023u);
  std::size_t group = 0;
  for (int i = 0; i < 5; ++i) {
    CollisionResult result = sample_collision(fuel, group, rng);
    std::cout << std::setw(2) << i << ": 事件 = " << result.channel
              << ", 入射群 = " << group << ", 出射群 = " << result.outgoing_group << "\n";
    group = result.outgoing_group;
  }
}
