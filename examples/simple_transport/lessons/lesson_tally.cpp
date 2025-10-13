#include <iomanip>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

// 本示例拆解 OpenMC 计分系统中最常见的 track-length tally。通过一个极简
// 的粒子轨迹集，展示 TallyManager 如何按照 (cell, group) 维度累积数据。

namespace lesson_tally {

struct CellInfo {
  std::string name;
  int id;
};

struct ParticleState {
  const CellInfo* cell = nullptr;
  std::size_t group = 0;
  double weight = 1.0;
};

class TallyManager {
 public:
  TallyManager(std::size_t cell_count, std::size_t group_count)
      : cell_pathlength_(cell_count, std::vector<double>(group_count, 0.0)) {}

  // track-length 计分的核心接口：以“路径长度 × 权重”累积。
  void score_flux(const ParticleState& p, double track_length) {
    if (!p.cell) return;  // 空指针用于表示粒子已经泄漏或尚未定位。
    auto& row = cell_pathlength_.at(p.cell->id);
    row.at(p.group) += track_length * p.weight;
  }

  void score_leakage(const ParticleState& p) { leakage_ += p.weight; }

  [[nodiscard]] const std::vector<std::vector<double>>& cell_pathlength() const {
    return cell_pathlength_;
  }

  [[nodiscard]] double leakage() const noexcept { return leakage_; }

 private:
  std::vector<std::vector<double>> cell_pathlength_;
  double leakage_ = 0.0;
};

}  // namespace lesson_tally

int main() {
  using namespace lesson_tally;
  std::cout << std::fixed << std::setprecision(4);

  // 1. 准备两个 Cell：模拟燃料与慢化剂。ID 即数组索引，保持与 OpenMC 一致的映射策略。
  CellInfo fuel{"fuel", 0};
  CellInfo moderator{"moderator", 1};

  // 2. 构造粒子轨迹：每个元素代表一段路径长度及其所属单元、能群、权重。
  std::vector<std::pair<ParticleState, double>> tracks = {
      {{&fuel, 0, 1.0}, 2.5},   // 燃料群0，长度 2.5 cm
      {{&fuel, 1, 0.8}, 1.2},   // 燃料群1，权重 0.8
      {{&moderator, 0, 1.0}, 3.7},
      {{&moderator, 2, 1.0}, 0.9},
      {{nullptr, 0, 1.0}, 0.0}  // 泄漏粒子，用于演示 score_leakage
  };

  TallyManager tallies(2, 3);

  for (const auto& [state, length] : tracks) {
    if (state.cell == nullptr) {
      tallies.score_leakage(state);
    } else {
      tallies.score_flux(state, length);
    }
  }

  std::cout << "按单元/能群累积的路径长度：\n";
  for (const auto& cell : {fuel, moderator}) {
    std::cout << "  " << cell.name << ':';
    const auto& row = tallies.cell_pathlength().at(cell.id);
    for (std::size_t g = 0; g < row.size(); ++g) {
      std::cout << " 群" << g << '=' << row[g];
    }
    std::cout << '\n';
  }

  std::cout << "总泄漏权重 = " << tallies.leakage() << '\n';
}
