#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// 本课程展示如何按照 OpenMC 的接口风格设计“并行代际驱动器”：
//   * 每个线程维护独立的随机数流与 `TransportKernel`，保证无锁计算；
//   * 主线程负责划分粒子区间、聚合计分与裂变银行；
//   * 通过 RAII 封装线程结果，避免数据竞争并保持面向对象结构清晰。

namespace lesson_parallel_driver {

// ============================= 基础设施：向量与随机引擎 =============================
struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 1.0;
};

class RandomEngine {
 public:
  explicit RandomEngine(unsigned seed) : rng_(seed), uniform_(0.0, 1.0) {}

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
    std::uniform_int_distribution<std::size_t> dist(0, size - 1);
    return dist(rng_);
  }

 private:
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_;
};

// ============================= 多群截面容器 =============================
class MultiGroupXS {
 public:
  struct Branch {
    enum class Type { Absorption, Scatter, Fission };
    Type type;
    std::size_t target_group = 0;
  };

  MultiGroupXS(std::vector<double> sigma_a,
               std::vector<double> sigma_s,
               std::vector<double> sigma_f,
               std::vector<double> nu,
               std::vector<std::vector<double>> scatter_matrix,
               std::vector<double> chi)
      : sigma_a_(std::move(sigma_a)),
        sigma_s_(std::move(sigma_s)),
        sigma_f_(std::move(sigma_f)),
        nu_(std::move(nu)),
        scatter_matrix_(std::move(scatter_matrix)),
        chi_(std::move(chi)) {
    validate();
    build_tables();
  }

  std::size_t groups() const { return sigma_a_.size(); }
  double sigma_t(std::size_t g) const { return sigma_t_.at(g); }
  double nu(std::size_t g) const { return nu_.at(g); }

  const std::vector<Branch>& branches(std::size_t g) const { return branches_.at(g); }
  std::size_t sample_branch(std::size_t g, RandomEngine& rng) const { return rng.sample_from_cdf(branch_cdf_.at(g)); }
  std::size_t sample_fission_group(RandomEngine& rng) const { return rng.sample_from_cdf(chi_cdf_); }

 private:
  void validate() const {
    std::size_t g = sigma_a_.size();
    auto check = [&](const std::vector<double>& vec, const char* label) {
      if (vec.size() != g) {
        throw std::invalid_argument(std::string(label) + " 长度不匹配能群数");
      }
    };
    check(sigma_s_, "sigma_s");
    check(sigma_f_, "sigma_f");
    check(nu_, "nu");
    if (scatter_matrix_.size() != g) throw std::invalid_argument("散射矩阵行数错误");
    for (const auto& row : scatter_matrix_) {
      if (row.size() != g) throw std::invalid_argument("散射矩阵列数错误");
    }
    if (chi_.size() != g) throw std::invalid_argument("裂变谱长度应等于能群数");
  }

  void build_tables() {
    std::size_t g = groups();
    sigma_t_.assign(g, 0.0);
    branches_.assign(g, {});
    branch_cdf_.assign(g, {});

    for (std::size_t i = 0; i < g; ++i) {
      sigma_t_[i] = sigma_a_[i] + sigma_s_[i] + sigma_f_[i];
      if (sigma_t_[i] <= 0.0) throw std::invalid_argument("总截面不能为零");

      std::vector<Branch> table;
      std::vector<double> probs;
      if (sigma_a_[i] > 0.0) {
        table.push_back({Branch::Type::Absorption, i});
        probs.push_back(sigma_a_[i]);
      }
      if (sigma_s_[i] > 0.0) {
        for (std::size_t j = 0; j < g; ++j) {
          double rate = scatter_matrix_[i][j];
          if (rate <= 0.0) continue;
          table.push_back({Branch::Type::Scatter, j});
          probs.push_back(rate);
        }
      }
      if (sigma_f_[i] > 0.0) {
        table.push_back({Branch::Type::Fission, i});
        probs.push_back(sigma_f_[i]);
      }

      double total = std::accumulate(probs.begin(), probs.end(), 0.0);
      std::vector<double> cdf;
      double accum = 0.0;
      for (double p : probs) {
        accum += p / total;
        cdf.push_back(accum);
      }
      branches_[i] = std::move(table);
      branch_cdf_[i] = std::move(cdf);
    }

    double chi_total = std::accumulate(chi_.begin(), chi_.end(), 0.0);
    double accum = 0.0;
    chi_cdf_.clear();
    for (double value : chi_) {
      accum += value / chi_total;
      chi_cdf_.push_back(accum);
    }
  }

  std::vector<double> sigma_a_;
  std::vector<double> sigma_s_;
  std::vector<double> sigma_f_;
  std::vector<double> nu_;
  std::vector<std::vector<double>> scatter_matrix_;
  std::vector<double> chi_;

  std::vector<double> sigma_t_;
  std::vector<std::vector<Branch>> branches_;
  std::vector<std::vector<double>> branch_cdf_;
  std::vector<double> chi_cdf_;
};

// ============================= 粒子状态与线程结果 =============================
struct Particle {
  Vec3 position{0.0, 0.0, 0.0};
  Vec3 direction{0.0, 0.0, 1.0};
  std::size_t group = 0;
  double weight = 1.0;
};

struct ThreadResult {
  double produced = 0.0;
  double track = 0.0;
  std::vector<Particle> fission_bank;
};

// ============================= 无共享状态的输运核 =============================
class TransportKernel {
 public:
  explicit TransportKernel(const MultiGroupXS& xs) : xs_(xs) {}

  ThreadResult run_batch(const std::vector<Particle>& batch, RandomEngine& rng) const {
    ThreadResult result;
    for (const auto& particle : batch) {
      result = accumulate(result, transport_single(particle, rng));
    }
    return result;
  }

 private:
  ThreadResult transport_single(Particle particle, RandomEngine& rng) const {
    ThreadResult result;
    constexpr int kMaxCollisions = 256;

    for (int step = 0; step < kMaxCollisions; ++step) {
      double free_path = rng.sample_exponential(xs_.sigma_t(particle.group));
      result.track += particle.weight * free_path;

      std::size_t branch_idx = xs_.sample_branch(particle.group, rng);
      const auto& branch = xs_.branches(particle.group).at(branch_idx);

      if (branch.type == MultiGroupXS::Branch::Type::Absorption) {
        break;
      } else if (branch.type == MultiGroupXS::Branch::Type::Scatter) {
        particle.group = branch.target_group;
        particle.direction = rng.sample_isotropic_direction();
        continue;
      } else {
        double nu = xs_.nu(particle.group);
        result.produced += particle.weight * nu;
        int children = std::max(1, static_cast<int>(std::round(nu)));
        double child_weight = particle.weight * nu / static_cast<double>(children);
        for (int n = 0; n < children; ++n) {
          Particle newborn;
          newborn.group = xs_.sample_fission_group(rng);
          newborn.direction = rng.sample_isotropic_direction();
          newborn.weight = child_weight;
          result.fission_bank.push_back(newborn);
        }
        break;
      }
    }
    return result;
  }

  ThreadResult accumulate(ThreadResult lhs, const ThreadResult& rhs) const {
    lhs.produced += rhs.produced;
    lhs.track += rhs.track;
    lhs.fission_bank.insert(lhs.fission_bank.end(), rhs.fission_bank.begin(), rhs.fission_bank.end());
    return lhs;
  }

  const MultiGroupXS& xs_;
};

// ============================= 并行驱动器 =============================
class ParallelEigenvalueDriver {
 public:
  ParallelEigenvalueDriver(const MultiGroupXS& xs, int histories_per_gen, int generations, int threads)
      : xs_(xs),
        histories_(histories_per_gen),
        generations_(generations),
        threads_(threads),
        kernel_(xs) {
    if (threads_ <= 0) throw std::invalid_argument("线程数必须为正");
    source_bank_.resize(histories_);
    RandomEngine rng(98765u);
    for (auto& p : source_bank_) {
      p.group = rng.sample_uniform_index(xs_.groups());
      p.direction = rng.sample_isotropic_direction();
    }
  }

  void run() {
    std::cout << "===== 并行代际驱动教学 =====\n";
    for (int gen = 0; gen < generations_; ++gen) {
      auto results = launch_threads();

      double produced = 0.0;
      double track = 0.0;
      std::vector<Particle> fission_bank;
      for (const auto& res : results) {
        produced += res.produced;
        track += res.track;
        fission_bank.insert(fission_bank.end(), res.fission_bank.begin(), res.fission_bank.end());
      }

      double k_eff = produced / static_cast<double>(source_bank_.size());
      std::cout << "第 " << (gen + 1) << " 代 (线程数=" << threads_ << ") k-effective = " << std::setprecision(5)
                << k_eff << "，平均 track-length = " << track / source_bank_.size() << "\n";

      if (fission_bank.empty()) {
        std::cout << "裂变银行为空，终止计算。\n";
        break;
      }

      resample_source(fission_bank);
    }
  }

 private:
  std::vector<ThreadResult> launch_threads() {
    std::vector<ThreadResult> results(threads_);
    std::vector<std::thread> workers;
    workers.reserve(threads_);

    int chunk = histories_ / threads_;
    int remainder = histories_ % threads_;
    int offset = 0;

    for (int t = 0; t < threads_; ++t) {
      int count = chunk + (t < remainder ? 1 : 0);
      std::vector<Particle> slice(source_bank_.begin() + offset, source_bank_.begin() + offset + count);
      offset += count;

      workers.emplace_back([this, t, slice = std::move(slice), &results]() {
        RandomEngine local_rng(1234u + 7919u * static_cast<unsigned>(t));
        ThreadResult local = kernel_.run_batch(slice, local_rng);
        results[t] = std::move(local);
      });
    }

    for (auto& worker : workers) {
      worker.join();
    }
    return results;
  }

  void resample_source(const std::vector<Particle>& fission_bank) {
    RandomEngine rng(45678u);
    std::vector<Particle> next(histories_);
    for (auto& particle : next) {
      std::size_t idx = rng.sample_uniform_index(fission_bank.size());
      particle = fission_bank[idx];
    }
    source_bank_ = std::move(next);
  }

  const MultiGroupXS& xs_;
  int histories_;
  int generations_;
  int threads_;
  TransportKernel kernel_;
  std::vector<Particle> source_bank_;
};

// ============================= 示例入口 =============================
void run_demo() {
  MultiGroupXS xs(/*sigma_a=*/{0.1, 0.05},
                  /*sigma_s=*/{0.35, 0.45},
                  /*sigma_f=*/{0.45, 0.40},
                  /*nu=*/{2.45, 2.60},
                  /*scatter_matrix=*/{{0.10, 0.25}, {0.05, 0.40}},
                  /*chi=*/{0.70, 0.30});

  ParallelEigenvalueDriver driver(xs, /*histories_per_gen=*/120, /*generations=*/4, /*threads=*/4);
  driver.run();
}

}  // namespace lesson_parallel_driver

int main() {
  using namespace lesson_parallel_driver;
  run_demo();
}
