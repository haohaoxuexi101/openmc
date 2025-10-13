#include <array>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// 本课程演示 OpenMC 几何加速结构中的“近邻表 (neighbor list)”思想。
// 在复杂几何中，粒子穿越单元边界时需要迅速找到下一个可能的 Cell。
// OpenMC 通过在初始化阶段构建邻接映射（`Cell::neighbors_` 等）来避免运行期搜索。
// 这里我们在 3×2×2 的笛卡尔晶格上构建邻接表，并通过接口展示如何查询。

namespace lesson_neighbor {

struct CellId {
  int i{};
  int j{};
  int k{};
};

// 简化的晶格描述：记录尺寸与单元索引。
struct Lattice3D {
  int nx{};
  int ny{};
  int nz{};

  int flat_index(int i, int j, int k) const { return (k * ny + j) * nx + i; }
};

// 邻接信息：六个方向的邻居索引，-1 表示边界（可绑定边界条件）。
struct NeighborEntry {
  std::array<int, 6> neighbors{};  // +/-x, +/-y, +/-z
};

class NeighborTable {
 public:
  NeighborTable(const Lattice3D& lattice, std::vector<std::string> universe)
      : lattice_(lattice), universe_(std::move(universe)) {
    build();
  }

  const NeighborEntry& query(int flat) const { return table_.at(flat); }

  const std::string& universe_name(int flat) const { return universe_.at(flat); }

 private:
  void build() {
    table_.resize(universe_.size());
    for (int k = 0; k < lattice_.nz; ++k) {
      for (int j = 0; j < lattice_.ny; ++j) {
        for (int i = 0; i < lattice_.nx; ++i) {
          int flat = lattice_.flat_index(i, j, k);
          NeighborEntry entry{};
          entry.neighbors[0] = (i > 0) ? lattice_.flat_index(i - 1, j, k) : -1;
          entry.neighbors[1] = (i + 1 < lattice_.nx) ? lattice_.flat_index(i + 1, j, k) : -1;
          entry.neighbors[2] = (j > 0) ? lattice_.flat_index(i, j - 1, k) : -1;
          entry.neighbors[3] = (j + 1 < lattice_.ny) ? lattice_.flat_index(i, j + 1, k) : -1;
          entry.neighbors[4] = (k > 0) ? lattice_.flat_index(i, j, k - 1) : -1;
          entry.neighbors[5] = (k + 1 < lattice_.nz) ? lattice_.flat_index(i, j, k + 1) : -1;
          table_[flat] = entry;
        }
      }
    }
  }

  Lattice3D lattice_;
  std::vector<std::string> universe_;
  std::vector<NeighborEntry> table_;
};

}  // namespace lesson_neighbor

int main() {
  using namespace lesson_neighbor;

  Lattice3D lattice{3, 2, 2};
  std::vector<std::string> universe_names(lattice.nx * lattice.ny * lattice.nz);
  for (int k = 0; k < lattice.nz; ++k) {
    for (int j = 0; j < lattice.ny; ++j) {
      for (int i = 0; i < lattice.nx; ++i) {
        int flat = lattice.flat_index(i, j, k);
        universe_names[flat] = "U_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
      }
    }
  }

  NeighborTable table(lattice, universe_names);

  int target = lattice.flat_index(1, 0, 1);
  const NeighborEntry& entry = table.query(target);

  std::cout << "查询单元: " << table.universe_name(target) << "\n";
  static const char* labels[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};
  for (int dir = 0; dir < 6; ++dir) {
    int neighbor_index = entry.neighbors[dir];
    if (neighbor_index == -1) {
      std::cout << std::setw(3) << labels[dir] << " -> 边界 (真空/反射可在此挂接)\n";
    } else {
      std::cout << std::setw(3) << labels[dir] << " -> " << table.universe_name(neighbor_index)
                << " (flat index = " << neighbor_index << ")\n";
    }
  }
}
