#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// 本课程专注于“多群截面制作（MGXS builder）”的核心概念。
// OpenMC 在读取连续能量核数据后，会根据能群结构与指定权重函数
// 对微观截面做加权平均，生成多群宏观截面。本示例抽象出以下步骤：
//   1. NuclideXS：存储核素在细能格上的截面（含各 Legendre 阶数系数）。
//   2. Spectrum：定义折合权重（如无限介质光谱）。
//   3. MGXSEntry：对每个能群求加权平均，输出宏观 σ_t、σ_s(l=0,1)、νσ_f、χ。
//   4. MaterialMixer：根据材料组成与原子数密度组合多个核素。
// 所有接口均模仿 OpenMC `openmc/mgxs.h` 中的数据结构，并通过中文注释详解设计思路。

namespace lesson_mgxs_builder {

// 细能格截面表，包含总截面、吸收、裂变以及各阶 Legendre 散射系数。
struct NuclideXS {
  std::string name;                       // 核素名称，仅用于输出。
  std::vector<double> energy;             // 细能格节点（升序排列）。
  std::vector<double> sigma_t;            // 每个节点的总截面。
  std::vector<double> sigma_a;            // 吸收截面。
  std::vector<double> nu_sigma_f;         // νΣ_f，直接存储便于折合。
  std::vector<double> p0;                 // 零阶 Legendre 系数（各向同性分量）。
  std::vector<double> p1;                 // 一阶 Legendre 系数（各向异性分量）。

  // 通过断言确保所有数组长度一致，模仿 OpenMC 在构造时的防御式检查。
  void validate() const {
    auto node_count = energy.size();
    if (node_count < 2) {
      throw std::runtime_error("能量节点至少需要两个");
    }
    std::size_t segment_count = node_count - 1;
    if (sigma_t.size() != segment_count || sigma_a.size() != segment_count ||
        nu_sigma_f.size() != segment_count || p0.size() != segment_count ||
        p1.size() != segment_count) {
      throw std::runtime_error("核素截面表维度不一致");
    }
  }
};

// 折合权重，示例使用一份“伪无限介质”能谱，数组长度与 energy 匹配。
struct Spectrum {
  std::vector<double> energy;
  std::vector<double> flux;  // 每个细能格节点的权重（相对值即可）。

  void normalize() {
    double sum = std::accumulate(flux.begin(), flux.end(), 0.0);
    for (double& value : flux) {
      value /= sum;
    }
  }
};

// 单一核素折合到多群后的结果。
struct MGXSEntry {
  std::vector<double> sigma_t;        // 多群总截面。
  std::vector<double> sigma_a;        // 多群吸收截面。
  std::vector<double> nu_sigma_f;     // 多群 νΣ_f。
  std::vector<double> scatter_p0;     // P0 散射矩阵的行向量（下标 [g*G + g']）。
  std::vector<double> scatter_p1;     // P1 散射矩阵。
  std::vector<double> chi;            // 裂变谱（本示例与能群数相同）。
};

// MGXSBuilder 负责执行“能群折合”。接口设计仿照 OpenMC `MGXS::collapse_to_group`。
class MGXSBuilder {
 public:
  MGXSBuilder(std::vector<double> group_bounds, Spectrum spectrum)
      : bounds_(std::move(group_bounds)), spectrum_(std::move(spectrum)) {
    spectrum_.normalize();
  }

  // 对单个核素执行折合，返回 MGXSEntry。
  MGXSEntry collapse(const NuclideXS& xs) const {
    xs.validate();
    MGXSEntry result;
    std::size_t groups = bounds_.size() - 1;
    result.sigma_t.resize(groups, 0.0);
    result.sigma_a.resize(groups, 0.0);
    result.nu_sigma_f.resize(groups, 0.0);
    result.scatter_p0.assign(groups * groups, 0.0);
    result.scatter_p1.assign(groups * groups, 0.0);
    result.chi.resize(groups, 0.0);

    // 遍历细能格，将每个能段累加到对应多群。
    for (std::size_t g = 0; g < groups; ++g) {
      double e_low = bounds_[g];
      double e_high = bounds_[g + 1];
      // 查找细能格区间。
      for (std::size_t i = 0; i + 1 < xs.energy.size(); ++i) {
        double seg_low = xs.energy[i];
        double seg_high = xs.energy[i + 1];
        if (seg_high <= e_low || seg_low >= e_high) {
          continue;  // 能段不与群范围重叠。
        }
        double overlap = std::min(seg_high, e_high) - std::max(seg_low, e_low);
        if (overlap <= 0.0) {
          continue;
        }
        double weight = segment_weight(i, overlap);
        result.sigma_t[g] += weight * xs.sigma_t[i];
        result.sigma_a[g] += weight * xs.sigma_a[i];
        result.nu_sigma_f[g] += weight * xs.nu_sigma_f[i];
        result.scatter_p0[g * groups + g] += weight * xs.p0[i];
        result.scatter_p1[g * groups + g] += weight * xs.p1[i];
        result.chi[g] += weight * xs.nu_sigma_f[i];
      }
    }

    // χ 需要按 νΣ_f 权重归一化。
    double chi_norm = std::accumulate(result.chi.begin(), result.chi.end(), 0.0);
    if (chi_norm > 0.0) {
      for (double& val : result.chi) {
        val /= chi_norm;
      }
    }
    return result;
  }

 private:
  // segment_weight 根据谱权重对能段求平均。
  double segment_weight(std::size_t index, double overlap) const {
    double w = spectrum_.flux[index] + spectrum_.flux[index + 1];
    return 0.5 * w * overlap;
  }

  std::vector<double> bounds_;
  Spectrum spectrum_;
};

// 材料配方：记录核素名称与原子数密度。
struct MaterialRecipe {
  std::string name;
  std::vector<std::pair<std::string, double>> nuclides;  // (核素, 原子密度)
};

// MaterialMixer 将核素多群截面乘以原子密度后叠加，得到宏观截面。
class MaterialMixer {
 public:
  explicit MaterialMixer(std::map<std::string, MGXSEntry> library)
      : library_(std::move(library)) {}

  MGXSEntry mix(const MaterialRecipe& recipe) const {
    auto groups = library_.begin()->second.sigma_t.size();
    MGXSEntry macro;
    macro.sigma_t.assign(groups, 0.0);
    macro.sigma_a.assign(groups, 0.0);
    macro.nu_sigma_f.assign(groups, 0.0);
    macro.scatter_p0.assign(groups * groups, 0.0);
    macro.scatter_p1.assign(groups * groups, 0.0);
    macro.chi.assign(groups, 0.0);

    for (const auto& [name, density] : recipe.nuclides) {
      const auto& xs = library_.at(name);
      for (std::size_t g = 0; g < groups; ++g) {
        macro.sigma_t[g] += density * xs.sigma_t[g];
        macro.sigma_a[g] += density * xs.sigma_a[g];
        macro.nu_sigma_f[g] += density * xs.nu_sigma_f[g];
        macro.chi[g] += density * xs.chi[g];
      }
      for (std::size_t idx = 0; idx < groups * groups; ++idx) {
        macro.scatter_p0[idx] += density * xs.scatter_p0[idx];
        macro.scatter_p1[idx] += density * xs.scatter_p1[idx];
      }
    }

    double chi_norm = std::accumulate(macro.chi.begin(), macro.chi.end(), 0.0);
    if (chi_norm > 0.0) {
      for (double& val : macro.chi) {
        val /= chi_norm;
      }
    }
    return macro;
  }

 private:
  std::map<std::string, MGXSEntry> library_;
};

}  // namespace lesson_mgxs_builder

int main() {
  using namespace lesson_mgxs_builder;

  // Step 1: 准备两种核素的细能格截面数据（数据为演示用）。
  NuclideXS u235{"U-235",
                 {0.0, 0.5, 1.0, 2.0},
                 {6.0, 5.5, 4.0},
                 {1.5, 1.2, 1.0},
                 {1.8, 1.6, 1.4},
                 {4.0, 3.5, 2.5},
                 {0.5, 0.4, 0.3}};
  NuclideXS u238{"U-238",
                 {0.0, 0.5, 1.0, 2.0},
                 {7.0, 6.0, 5.0},
                 {0.5, 0.4, 0.3},
                 {0.0, 0.0, 0.0},
                 {5.5, 4.8, 3.6},
                 {0.2, 0.15, 0.1}};

  Spectrum spectrum{{0.0, 0.5, 1.0, 2.0}, {1.0, 0.8, 0.6, 0.4}};

  MGXSBuilder builder({0.0, 0.6, 2.0}, spectrum);
  MGXSEntry u235_group = builder.collapse(u235);
  MGXSEntry u238_group = builder.collapse(u238);

  std::map<std::string, MGXSEntry> library{{u235.name, u235_group}, {u238.name, u238_group}};
  MaterialMixer mixer(library);

  MaterialRecipe recipe{"燃料", {{"U-235", 0.05}, {"U-238", 0.95}}};
  MGXSEntry macro = mixer.mix(recipe);

  std::cout << "=== 多群宏观截面 ===\n";
  for (std::size_t g = 0; g < macro.sigma_t.size(); ++g) {
    std::cout << "群 " << g << ": Σ_t = " << std::setw(8) << macro.sigma_t[g]
              << ", Σ_a = " << std::setw(8) << macro.sigma_a[g]
              << ", νΣ_f = " << std::setw(8) << macro.nu_sigma_f[g]
              << ", χ = " << std::setw(6) << macro.chi[g] << "\n";
  }

  std::cout << "\n=== P0 散射矩阵（对角元素示例）===\n";
  std::size_t groups = macro.sigma_t.size();
  for (std::size_t g = 0; g < groups; ++g) {
    std::cout << "Σ_s^{(" << g << "→" << g << ")} = "
              << macro.scatter_p0[g * groups + g] << "\n";
  }

  std::cout << "\n=== P1 各向异性系数（对角元素）===\n";
  for (std::size_t g = 0; g < groups; ++g) {
    std::cout << "Σ_{s,1}^{(" << g << "→" << g << ")} = "
              << macro.scatter_p1[g * groups + g] << "\n";
  }
}
