#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

// 本课程拆解 OpenMC 中“分布与抽样”模块的设计理念：
//   - 抽象基类 Distribution，暴露统一的 sample() 接口；
//   - PiecewiseLinearDistribution 对应 ENDF TAB1 数据，用面积归一化实现采样；
//   - AliasMethod 展示如何将离散分布转换为常数时间抽样结构，与 OpenMC 的
//     `DiscreteDistribution::build_alias_table` 相呼应；
//   - MixtureDistribution 演示组合多个分布，实现复合抽样。
// 所有类均配套中文注释说明接口选择与数值细节，帮助理解 OpenMC 的抽样体系。

namespace lesson_distribution {

class Distribution {
 public:
  virtual ~Distribution() = default;
  virtual double sample(std::mt19937& rng) const = 0;
};

// -------- 分段线性分布：对应连续能量散射、能量分布常见格式 --------
class PiecewiseLinearDistribution : public Distribution {
 public:
  PiecewiseLinearDistribution(std::vector<double> grid, std::vector<double> pdf)
      : grid_(std::move(grid)), pdf_(std::move(pdf)) {
    if (grid_.size() != pdf_.size()) {
      throw std::runtime_error("网格与 pdf 长度不一致");
    }
    normalize();
  }

  double sample(std::mt19937& rng) const override {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double xi = uniform(rng);
    std::size_t upper = std::lower_bound(cdf_.begin(), cdf_.end(), xi) - cdf_.begin();
    if (upper == 0) {
      upper = 1;
    }
    std::size_t lower = upper - 1;
    double x0 = grid_[lower];
    double x1 = grid_[upper];
    double pdf0 = pdf_[lower];
    double pdf1 = pdf_[upper];
    double slope = (pdf1 - pdf0) / (x1 - x0);
    double cdf0 = cdf_[lower];
    double cdf1 = cdf_[upper];
    double target = xi - cdf0;
    if (std::fabs(slope) < 1e-14) {
      return x0 + target / pdf0;
    }
    double a = 0.5 * slope;
    double b = pdf0;
    double c = -target;
    double discriminant = std::max(0.0, b * b - 4.0 * a * c);
    double root = (-b + std::sqrt(discriminant)) / (2.0 * a);
    return x0 + root;
  }

 private:
  void normalize() {
    double integral = 0.0;
    for (std::size_t i = 0; i + 1 < grid_.size(); ++i) {
      double width = grid_[i + 1] - grid_[i];
      integral += 0.5 * (pdf_[i] + pdf_[i + 1]) * width;
    }
    for (double& p : pdf_) {
      p /= integral;
    }
    cdf_.resize(grid_.size());
    cdf_[0] = 0.0;
    for (std::size_t i = 0; i + 1 < grid_.size(); ++i) {
      double width = grid_[i + 1] - grid_[i];
      cdf_[i + 1] = cdf_[i] + 0.5 * (pdf_[i] + pdf_[i + 1]) * width;
    }
  }

  std::vector<double> grid_;
  mutable std::vector<double> pdf_;
  std::vector<double> cdf_;
};

// -------- Alias Method：离散分布的常数时间抽样 --------
class AliasTable {
 public:
  explicit AliasTable(std::vector<double> probabilities)
      : prob_(probabilities.size()), alias_(probabilities.size()) {
    std::size_t n = probabilities.size();
    std::vector<double> scaled(probabilities);
    std::vector<std::size_t> small;
    std::vector<std::size_t> large;

    for (std::size_t i = 0; i < n; ++i) {
      scaled[i] *= n;
      if (scaled[i] < 1.0) {
        small.push_back(i);
      } else {
        large.push_back(i);
      }
    }

    while (!small.empty() && !large.empty()) {
      std::size_t s = small.back();
      small.pop_back();
      std::size_t l = large.back();
      alias_[s] = l;
      prob_[s] = scaled[s];
      scaled[l] = scaled[l] - (1.0 - scaled[s]);
      if (scaled[l] < 1.0) {
        large.pop_back();
        small.push_back(l);
      }
    }
    for (std::size_t idx : large) {
      prob_[idx] = 1.0;
    }
    for (std::size_t idx : small) {
      prob_[idx] = 1.0;
    }
  }

  std::size_t sample(std::mt19937& rng) const {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double xi = uniform(rng);
    double y = uniform(rng);
    std::size_t column = static_cast<std::size_t>(xi * prob_.size());
    double threshold = prob_[column];
    return y < threshold ? column : alias_[column];
  }

 private:
  std::vector<double> prob_;
  std::vector<std::size_t> alias_;
};

class AliasDistribution : public Distribution {
 public:
  AliasDistribution(std::vector<double> values, std::vector<double> probabilities)
      : values_(std::move(values)), alias_(std::move(probabilities)) {}

  double sample(std::mt19937& rng) const override {
    std::size_t idx = alias_.sample(rng);
    return values_[idx];
  }

 private:
  std::vector<double> values_;
  AliasTable alias_;
};

// -------- MixtureDistribution：组合多个分布 --------
class MixtureDistribution : public Distribution {
 public:
  MixtureDistribution(std::vector<double> weights, std::vector<std::unique_ptr<Distribution>> parts)
      : weights_(std::move(weights)), parts_(std::move(parts)) {
    double sum = std::accumulate(weights_.begin(), weights_.end(), 0.0);
    for (double& w : weights_) {
      w /= sum;
    }
  }

  double sample(std::mt19937& rng) const override {
    std::discrete_distribution<std::size_t> select(weights_.begin(), weights_.end());
    std::size_t idx = select(rng);
    return parts_[idx]->sample(rng);
  }

 private:
  std::vector<double> weights_;
  std::vector<std::unique_ptr<Distribution>> parts_;
};

}  // namespace lesson_distribution

int main() {
  using namespace lesson_distribution;

  std::mt19937 rng{2025u};

  PiecewiseLinearDistribution slowing_down({0.0, 0.5, 1.0, 2.0}, {0.2, 0.5, 1.0, 1.4});
  AliasDistribution discrete({0.0, 0.5, 1.0}, {0.2, 0.3, 0.5});

  std::vector<std::unique_ptr<Distribution>> parts;
  parts.emplace_back(std::make_unique<PiecewiseLinearDistribution>(slowing_down));
  parts.emplace_back(std::make_unique<AliasDistribution>(discrete));
  MixtureDistribution mix({0.7, 0.3}, std::move(parts));

  std::cout << "分段线性分布样本:\n";
  for (int i = 0; i < 5; ++i) {
    std::cout << std::setw(10) << slowing_down.sample(rng) << "\n";
  }

  std::cout << "\nAlias 抽样样本:\n";
  for (int i = 0; i < 5; ++i) {
    std::cout << std::setw(10) << discrete.sample(rng) << "\n";
  }

  std::cout << "\n复合分布样本:\n";
  for (int i = 0; i < 5; ++i) {
    std::cout << std::setw(10) << mix.sample(rng) << "\n";
  }
}
