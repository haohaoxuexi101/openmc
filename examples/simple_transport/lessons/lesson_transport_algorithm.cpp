#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// 本课程专注解释 OpenMC 式中子输运算法的关键步骤：
//   1. 根据总截面抽样自由程；
//   2. 根据反应截面比例选择吸收 / 散射 / 裂变分支；
//   3. 散射更新能群与方向，裂变产生二次粒子并累积 k-effective；
//   4. 使用粒子银行重采样完成代际循环。
// 所有接口均保持高度模块化，方便直接拓展到复杂几何或并行框架中。

namespace lesson_transport_algorithm {

// ============================= 基础向量与随机工具 =============================
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;
};

class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed = 20240717u)
      : rng_(seed), uniform_(0.0, 1.0) {}

  double uniform() { return uniform_(rng_); }

  double sample_exponential(double sigma_t) {
    double xi = std::max(uniform(), 1e-12);
    return -std::log(xi) / sigma_t;
  }

  Vec3 sample_isotropic_direction() {
    double u = 2.0 * uniform() - 1.0;
    double phi = 2.0 * 3.14159265358979323846 * uniform();
    double r = std::sqrt(std::max(0.0, 1.0 - u * u));
    return {r * std::cos(phi), r * std::sin(phi), u};
  }

  std::size_t sample_from_cdf(const std::vector<double>& cdf) {
    double xi = uniform();
    for (std::size_t i = 0; i < cdf.size(); ++i) {
      if (xi <= cdf[i]) return i;
    }
    return cdf.size() - 1;
  }

  std::size_t sample_uniform_index(std::size_t size) {
    if (size == 0) throw std::runtime_error("无法从空集合抽样");
    std::uniform_int_distribution<std::size_t> dist(0, size - 1);
    return dist(rng_);
  }

 private:
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_;
};

// ============================= 多群材料：截面与散射矩阵 =============================
class MultiGroupMaterial {
 public:
  struct CollisionBranch {
    enum class Type { Absorption, Scatter, Fission };
    Type type;
    std::size_t target_group = 0;
  };

  MultiGroupMaterial(std::string name,
                      std::vector<double> sigma_a,
                      std::vector<double> sigma_s,
                      std::vector<double> sigma_f,
                      std::vector<double> nu,
                      std::vector<std::vector<double>> scatter_matrix,
                      std::vector<double> chi)
      : name_(std::move(name)),
        sigma_a_(std::move(sigma_a)),
        sigma_s_(std::move(sigma_s)),
        sigma_f_(std::move(sigma_f)),
        nu_(std::move(nu)),
        scatter_matrix_(std::move(scatter_matrix)),
        chi_(std::move(chi)) {
    validate_sizes();
    build_collision_tables();
  }

  const std::string& name() const { return name_; }
  std::size_t groups() const { return sigma_a_.size(); }
  double sigma_t(std::size_t g) const { return sigma_t_.at(g); }
  double nu(std::size_t g) const { return nu_.at(g); }

  const std::vector<CollisionBranch>& branches(std::size_t g) const { return branches_.at(g); }
  std::size_t sample_branch(std::size_t g, RandomEngine& rng) const {
    return rng.sample_from_cdf(branch_cdf_.at(g));
  }
  std::size_t sample_fission_group(RandomEngine& rng) const { return rng.sample_from_cdf(chi_cdf_); }

 private:
  void validate_sizes() const {
    std::size_t g = sigma_a_.size();
    auto check = [&](const std::vector<double>& vec, const char* label) {
      if (vec.size() != g) {
        throw std::invalid_argument(std::string(label) + " 的长度必须等于能群数");
      }
    };
    check(sigma_s_, "sigma_s");
    check(sigma_f_, "sigma_f");
    check(nu_, "nu");
    if (chi_.size() != g) {
      throw std::invalid_argument("裂变谱 chi 的长度必须等于能群数");
    }
    if (scatter_matrix_.size() != g) {
      throw std::invalid_argument("散射矩阵行数必须等于能群数");
    }
    for (const auto& row : scatter_matrix_) {
      if (row.size() != g) {
        throw std::invalid_argument("散射矩阵必须是 g×g");
      }
    }
  }

  void build_collision_tables() {
    std::size_t g = groups();
    sigma_t_.assign(g, 0.0);
    branches_.assign(g, {});
    branch_cdf_.assign(g, {});

    for (std::size_t i = 0; i < g; ++i) {
      sigma_t_[i] = sigma_a_[i] + sigma_s_[i] + sigma_f_[i];
      if (sigma_t_[i] <= 0.0) {
        throw std::invalid_argument("总截面不能为零");
      }

      std::vector<CollisionBranch> branches;
      std::vector<double> probabilities;

      if (sigma_a_[i] > 0.0) {
        branches.push_back({CollisionBranch::Type::Absorption, i});
        probabilities.push_back(sigma_a_[i]);
      }
      if (sigma_s_[i] > 0.0) {
        for (std::size_t j = 0; j < g; ++j) {
          double rate = scatter_matrix_[i][j];
          if (rate <= 0.0) continue;
          branches.push_back({CollisionBranch::Type::Scatter, j});
          probabilities.push_back(rate);
        }
      }
      if (sigma_f_[i] > 0.0) {
        branches.push_back({CollisionBranch::Type::Fission, i});
        probabilities.push_back(sigma_f_[i]);
      }

      double total_rate = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
      for (double& p : probabilities) {
        p /= total_rate;
      }

      std::vector<double> cdf;
      cdf.reserve(probabilities.size());
      double cumulative = 0.0;
      for (double p : probabilities) {
        cumulative += p;
        cdf.push_back(cumulative);
      }

      branches_[i] = std::move(branches);
      branch_cdf_[i] = std::move(cdf);
    }

    double chi_sum = std::accumulate(chi_.begin(), chi_.end(), 0.0);
    if (chi_sum <= 0.0) {
      throw std::invalid_argument("裂变谱总和必须大于零");
    }
    double cumulative = 0.0;
    chi_cdf_.clear();
    for (double value : chi_) {
      cumulative += value / chi_sum;
      chi_cdf_.push_back(cumulative);
    }
  }

  std::string name_;
  std::vector<double> sigma_a_;
  std::vector<double> sigma_s_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> scatter_matrix_;
  std::vector<double> chi_;

  std::vector<double> sigma_t_;
  std::vector<std::vector<CollisionBranch>> branches_;
  std::vector<std::vector<double>> branch_cdf_;
  std::vector<double> chi_cdf_;
};

// ============================= 粒子与历史结果 =============================
struct Particle {
  Vec3 position{0.0, 0.0, 0.0};
  Vec3 direction{0.0, 0.0, 1.0};
  std::size_t group = 0;
  double weight = 1.0;
};

struct HistoryTallies {
  double track_length = 0.0;
  double produced_neutrons = 0.0;
  std::vector<Particle> fission_bank;
};

// ============================= 物理驱动：单历史追踪 =============================
class TransportKernel {
 public:
  explicit TransportKernel(const MultiGroupMaterial& material) : material_(material) {}

  HistoryTallies transport(Particle particle, RandomEngine& rng) const {
    HistoryTallies tallies;
    constexpr int kMaxCollisions = 256;

    for (int step = 0; step < kMaxCollisions; ++step) {
      double sigma_t = material_.sigma_t(particle.group);
      double free_path = rng.sample_exponential(sigma_t);
      tallies.track_length += particle.weight * free_path;

      std::size_t branch_index = material_.sample_branch(particle.group, rng);
      const auto& branch = material_.branches(particle.group).at(branch_index);

      if (branch.type == MultiGroupMaterial::CollisionBranch::Type::Absorption) {
        break;  // 粒子被吸收，历史结束。
      } else if (branch.type == MultiGroupMaterial::CollisionBranch::Type::Scatter) {
        particle.group = branch.target_group;
        particle.direction = rng.sample_isotropic_direction();
        continue;  // 继续追踪后续碰撞。
      } else {
        // 裂变：记录产生的中子数，并生成二次粒子写入裂变银行。
        double nu = material_.nu(particle.group);
        tallies.produced_neutrons += particle.weight * nu;

        int children = std::max(1, static_cast<int>(std::round(nu)));
        double child_weight = particle.weight * nu / static_cast<double>(children);
        for (int i = 0; i < children; ++i) {
          Particle newborn;
          newborn.position = particle.position;
          newborn.direction = rng.sample_isotropic_direction();
          newborn.group = material_.sample_fission_group(rng);
          newborn.weight = child_weight;
          tallies.fission_bank.push_back(newborn);
        }
        break;
      }
    }

    return tallies;
  }

 private:
  const MultiGroupMaterial& material_;
};

// ============================= 代际驱动：k-effective 估计 =============================
class EigenvalueDriver {
 public:
  EigenvalueDriver(const MultiGroupMaterial& material, int histories_per_generation, int generations)
      : material_(material),
        kernel_(material),
        histories_(histories_per_generation),
        generations_(generations),
        rng_(12345u) {
    source_bank_.resize(histories_);
    for (int i = 0; i < histories_; ++i) {
      source_bank_[i].group = i % material_.groups();
      source_bank_[i].direction = rng_.sample_isotropic_direction();
    }
  }

  void run() {
    std::cout << "===== 中子输运算法教学：多群单区代际计算 =====\n";
    for (int gen = 0; gen < generations_; ++gen) {
      double produced = 0.0;
      double track_sum = 0.0;
      std::vector<Particle> fission_bank;

      for (const auto& particle : source_bank_) {
        auto tallies = kernel_.transport(particle, rng_);
        produced += tallies.produced_neutrons;
        track_sum += tallies.track_length;
        fission_bank.insert(fission_bank.end(), tallies.fission_bank.begin(), tallies.fission_bank.end());
      }

      double k_eff = produced / static_cast<double>(source_bank_.size());
      std::cout << "第 " << (gen + 1) << " 代 k-effective = " << std::setprecision(5) << k_eff
                << "，平均 track-length = " << track_sum / static_cast<double>(source_bank_.size()) << "\n";

      if (fission_bank.empty()) {
        std::cout << "裂变银行为空，计算提前结束。\n";
        break;
      }

      resample_source(fission_bank);
    }
  }

 private:
  void resample_source(const std::vector<Particle>& fission_bank) {
    std::vector<Particle> next(histories_);
    for (int i = 0; i < histories_; ++i) {
      std::size_t index = rng_.sample_uniform_index(fission_bank.size());
      next[i] = fission_bank[index];
    }
    source_bank_ = next;
  }

  const MultiGroupMaterial& material_;
  TransportKernel kernel_;
  int histories_;
  int generations_;
  RandomEngine rng_;
  std::vector<Particle> source_bank_;
};

// ============================= 教学示例入口 =============================
void run_demo() {
  // 定义 3 群材料参数，散射矩阵满足行和 = sigma_s。
  MultiGroupMaterial material(
      "tutorial_material",
      /*sigma_a=*/{0.12, 0.08, 0.05},
      /*sigma_s=*/{0.30, 0.45, 0.20},
      /*sigma_f=*/{0.40, 0.35, 0.50},
      /*nu=*/{2.45, 2.50, 2.60},
      /*scatter_matrix=*/{
          {0.05, 0.18, 0.07},
          {0.02, 0.32, 0.11},
          {0.01, 0.09, 0.10},
      },
      /*chi=*/{0.65, 0.25, 0.10});

  EigenvalueDriver driver(material, /*histories_per_generation=*/64, /*generations=*/4);
  driver.run();
}

}  // namespace lesson_transport_algorithm

int main() {
  using namespace lesson_transport_algorithm;
  run_demo();
}
