#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// 本课程聚焦“多普勒展宽 (Doppler broadening)”的基本算法。
// OpenMC 在加载核数据时，会根据温度对共振截面进行高斯卷积，从而反映原子热运动。
// 真实实现位于 `openmc/endf` 与 `openmc/thermal` 模块，包含复杂的表格插值。
// 这里通过一个简化模型演示关键思想：
//   1. 原始截面 sigma_0(E) 由一组离散点给出；
//   2. 通过高斯核 K(E, E', T) = 1/(sqrt(pi)*alpha) * exp(-(E'-E)^2/alpha^2) 做卷积；
//   3. alpha 与温度成正比，公式与自由气体模型一致；
//   4. 提供接口 evaluate(E, T) 与 broaden(T) 生成整张表格。

namespace lesson_doppler {

struct ResonanceCurve {
  std::vector<double> energy;
  std::vector<double> sigma0;

  void validate() const {
    if (energy.size() != sigma0.size()) {
      throw std::runtime_error("共振曲线维度不一致");
    }
  }
};

class DopplerBroaden {
 public:
  explicit DopplerBroaden(ResonanceCurve curve, double awr)
      : curve_(std::move(curve)), awr_(awr) {
    curve_.validate();
  }

  double evaluate(double E, double temperature) const {
    double alpha = alpha_width(E, temperature);
    double prefactor = 1.0 / (std::sqrt(M_PI) * alpha);
    double sum = 0.0;
    for (std::size_t i = 0; i < curve_.energy.size(); ++i) {
      double Ei = curve_.energy[i];
      double sigma = curve_.sigma0[i];
      double width = interval_width(i);
      double weight = std::exp(-std::pow(Ei - E, 2) / (alpha * alpha));
      sum += sigma * weight * width;
    }
    return prefactor * sum;
  }

  std::vector<double> broaden(double temperature) const {
    std::vector<double> result(curve_.energy.size());
    for (std::size_t i = 0; i < curve_.energy.size(); ++i) {
      result[i] = evaluate(curve_.energy[i], temperature);
    }
    return result;
  }

 private:
  double alpha_width(double E, double temperature) const {
    const double kB = 8.617333262e-5;  // eV/K
    double beta = std::sqrt(2.0 * kB * temperature / awr_);
    return beta * std::sqrt(E);
  }

  double interval_width(std::size_t i) const {
    if (curve_.energy.size() == 1) {
      return 1.0;
    }
    if (i == 0) {
      return curve_.energy[1] - curve_.energy[0];
    }
    if (i + 1 == curve_.energy.size()) {
      return curve_.energy[i] - curve_.energy[i - 1];
    }
    double left = curve_.energy[i] - curve_.energy[i - 1];
    double right = curve_.energy[i + 1] - curve_.energy[i];
    return 0.5 * (left + right);
  }

  ResonanceCurve curve_;
  double awr_;  // 相对原子质量比 AWR
};

}  // namespace lesson_doppler

int main() {
  using namespace lesson_doppler;

  ResonanceCurve curve{{0.1, 0.2, 0.3, 0.4, 0.5}, {1000.0, 4000.0, 8000.0, 4000.0, 1000.0}};
  DopplerBroaden broaden(curve, 238.0);

  double T_cold = 300.0;
  double T_hot = 1200.0;
  auto sigma_cold = broaden.broaden(T_cold);
  auto sigma_hot = broaden.broaden(T_hot);

  std::cout << "能量(eV)  σ(T=" << T_cold << "K)   σ(T=" << T_hot << "K)\n";
  for (std::size_t i = 0; i < curve.energy.size(); ++i) {
    std::cout << std::setw(8) << curve.energy[i] << std::setw(14) << sigma_cold[i]
              << std::setw(14) << sigma_hot[i] << "\n";
  }
}
