#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// 本课程展示“各向异性散射处理”的两种典型路径：
//   1. 多群模式：通过 Legendre 多项式展开（P0/P1）构造角度概率密度；
//   2. 连续能量模式：利用分段线性角度分布 (tabulated distribution) 采样余弦。
// 代码直接映射 OpenMC 中 `anisotropic_scatter.cpp`、`angle_distribution.cpp`
// 的设计方式：统一的接口，内部根据模式采用不同的抽样策略。

namespace lesson_scatter {

// ----------- 多群各向异性散射：Legendre 展开 + 拒绝采样 -----------

// P_l(mu) 计算器，只实现到 l=2 以演示结构，实际 OpenMC 支持更高阶。
inline double legendre(int l, double mu) {
  switch (l) {
    case 0:
      return 1.0;
    case 1:
      return mu;
    case 2:
      return 0.5 * (3.0 * mu * mu - 1.0);
    default:
      throw std::runtime_error("未实现的阶数");
  }
}

// 多群散射截面结构：每个阶数 (l) 一个矩阵 Σ_s,l^{g→g'}。
struct MultigroupScatterData {
  std::size_t groups{};
  std::vector<std::vector<double>> sigma_l;  // sigma_l[l][g*groups + g']

  double pdf(int l, std::size_t g, std::size_t gp) const {
    return sigma_l[l][g * groups + gp];
  }
};

// 使用 Legendre 展开构造角度概率密度函数。
class MultigroupAnisotropicSampler {
 public:
  MultigroupAnisotropicSampler(const MultigroupScatterData& data, std::mt19937 rng)
      : data_(data), rng_(std::move(rng)), uniform_(-1.0, 1.0), choice_(0.0, 1.0) {}

  // 返回 (出射能群, cosθ)，结构与 OpenMC `sample_mg_angle` 相似。
  std::pair<std::size_t, double> sample(std::size_t incoming_group) {
    // 1. 根据各向同性 (P0) 分布选择出射能群。
    std::vector<double> pdf(data_.groups);
    for (std::size_t gp = 0; gp < data_.groups; ++gp) {
      pdf[gp] = data_.pdf(0, incoming_group, gp);
    }
    std::discrete_distribution<std::size_t> select(pdf.begin(), pdf.end());
    std::size_t outgoing_group = select(rng_);

    // 2. 使用拒绝采样构建角度分布 f(mu) = 0.5 * [1 + a1 * 3mu + a2 * (5mu^2-1)]。
    double sigma_p0 = data_.pdf(0, incoming_group, outgoing_group);
    double sigma_p1 = data_.sigma_l.size() > 1 ? data_.pdf(1, incoming_group, outgoing_group) : 0.0;
    double sigma_p2 = data_.sigma_l.size() > 2 ? data_.pdf(2, incoming_group, outgoing_group) : 0.0;

    // 系数归一化：a_l = Σ_s,l / Σ_s,0。
    double a1 = sigma_p0 > 0.0 ? sigma_p1 / sigma_p0 : 0.0;
    double a2 = sigma_p0 > 0.0 ? sigma_p2 / sigma_p0 : 0.0;

    double mu = 0.0;
    while (true) {
      double candidate = uniform_(rng_);
      double pdf_val = 0.5 * (legendre(0, candidate) + 3.0 * a1 * legendre(1, candidate) +
                              5.0 * a2 * legendre(2, candidate));
      double bound = 0.5 * (1.0 + 3.0 * std::fabs(a1) + 5.0 * std::fabs(a2));
      if (choice_(rng_) * bound <= pdf_val) {
        mu = candidate;
        break;
      }
    }
    return {outgoing_group, mu};
  }

 private:
  const MultigroupScatterData& data_;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> uniform_;
  std::uniform_real_distribution<double> choice_;
};

// ----------- 连续能量各向异性散射：分段线性角度分布 -----------

// 模拟 OpenMC `AngleDistribution` 的 tabular 格式：
// mu_grid: 角余弦节点；cdf: 累积分布；pdf: 用于倾斜重构。
struct TabularAngularDistribution {
  std::vector<double> mu_grid;
  std::vector<double> pdf;
  std::vector<double> cdf;

  // 预处理：将 pdf 归一化并构建 cdf。
  void normalize() {
    double integral = 0.0;
    for (std::size_t i = 0; i + 1 < mu_grid.size(); ++i) {
      double mu0 = mu_grid[i];
      double mu1 = mu_grid[i + 1];
      double area = 0.5 * (pdf[i] + pdf[i + 1]) * (mu1 - mu0);
      integral += area;
    }
    for (double& value : pdf) {
      value /= integral;
    }
    cdf.resize(mu_grid.size());
    cdf[0] = 0.0;
    for (std::size_t i = 0; i + 1 < mu_grid.size(); ++i) {
      double mu0 = mu_grid[i];
      double mu1 = mu_grid[i + 1];
      double area = 0.5 * (pdf[i] + pdf[i + 1]) * (mu1 - mu0);
      cdf[i + 1] = cdf[i] + area;
    }
  }

  double sample_mu(std::mt19937& rng) const {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    double xi = uniform(rng);
    // 二分查找定位段。
    std::size_t upper = std::lower_bound(cdf.begin(), cdf.end(), xi) - cdf.begin();
    if (upper == 0) {
      upper = 1;
    }
    std::size_t lower = upper - 1;
    double cdf0 = cdf[lower];
    double cdf1 = cdf[upper];
    double mu0 = mu_grid[lower];
    double mu1 = mu_grid[upper];
    double pdf0 = pdf[lower];
    double pdf1 = pdf[upper];

    double width = mu1 - mu0;
    double slope = (pdf1 - pdf0) / width;
    double a = 0.5 * slope;
    double b = pdf0;
    double c = cdf0 - xi;
    double mu_rel;
    if (std::fabs(a) < 1e-12) {
      mu_rel = -c / b;
    } else {
      double discriminant = std::max(0.0, b * b - 4.0 * a * c);
      mu_rel = (-b + std::sqrt(discriminant)) / (2.0 * a);
    }
    return mu0 + mu_rel;
  }
};

}  // namespace lesson_scatter

int main() {
  using namespace lesson_scatter;

  // 多群示例：两群系统，给出 P0/P1/P2 三个矩阵（为方便仅设置对角元素）。
  MultigroupScatterData mg_data;
  mg_data.groups = 2;
  mg_data.sigma_l.resize(3);
  for (auto& row : mg_data.sigma_l) {
    row.assign(mg_data.groups * mg_data.groups, 0.0);
  }
  mg_data.sigma_l[0][0] = 0.4;  // Σ_s,0^{0→0}
  mg_data.sigma_l[0][3] = 0.3;  // Σ_s,0^{1→1}
  mg_data.sigma_l[1][0] = 0.05;
  mg_data.sigma_l[1][3] = -0.02;
  mg_data.sigma_l[2][0] = 0.01;
  mg_data.sigma_l[2][3] = 0.005;

  MultigroupAnisotropicSampler mg_sampler(mg_data, std::mt19937{2024u});
  for (int i = 0; i < 5; ++i) {
    auto [group, mu] = mg_sampler.sample(0);
    std::cout << "多群散射样本 " << i << ": 出射群 = " << group << ", cosθ = " << std::setw(8)
              << mu << "\n";
  }

  std::cout << "\n";

  // 连续能量示例：构造一个前向偏置的角度分布。
  TabularAngularDistribution dist;
  dist.mu_grid = {-1.0, -0.5, 0.0, 0.5, 1.0};
  dist.pdf = {0.2, 0.3, 0.6, 0.9, 1.4};
  dist.normalize();

  std::mt19937 rng{1234u};
  for (int i = 0; i < 5; ++i) {
    double mu = dist.sample_mu(rng);
    std::cout << "连续能量散射样本 " << i << ": cosθ = " << std::setw(8) << mu << "\n";
  }
}
