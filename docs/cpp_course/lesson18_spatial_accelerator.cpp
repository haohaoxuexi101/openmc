#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <vector>

// SpatialAccelerator 体现加速结构如何封装：
// - CellBounds 保存包围盒；
// - UniformGrid 负责建立索引列表，强调“构建期预处理 + 查询期常量时间”；
// - main 展示查询流程并打印候选 cell，说明减少逐个测试的优势。

struct AABB {
  std::array<double, 3> min;
  std::array<double, 3> max;

  bool contains(double x, double y, double z) const
  {
    return x >= min[0] && x <= max[0] && y >= min[1] && y <= max[1] && z >= min[2] && z <= max[2];
  }
};

struct CellBounds {
  std::string name;
  AABB bounds;
};

class UniformGrid {
public:
  UniformGrid(std::array<double, 3> low, std::array<double, 3> high, std::array<int, 3> div)
    : low_{low}
    , high_{high}
    , div_{div}
  {
    std::size_t total = static_cast<std::size_t>(div[0] * div[1] * div[2]);
    buckets_.resize(total);
  }

  void insert(std::size_t cell_index, const AABB& box)
  {
    auto to_index = [&](double value, int axis) {
      double normalized = (value - low_[axis]) / (high_[axis] - low_[axis]);
      int idx = static_cast<int>(normalized * div_[axis]);
      return std::clamp(idx, 0, div_[axis] - 1);
    };

    std::array<int, 3> lo{
      to_index(box.min[0], 0),
      to_index(box.min[1], 1),
      to_index(box.min[2], 2),
    };
    std::array<int, 3> hi{
      to_index(box.max[0], 0),
      to_index(box.max[1], 1),
      to_index(box.max[2], 2),
    };

    for (int k = lo[2]; k <= hi[2]; ++k) {
      for (int j = lo[1]; j <= hi[1]; ++j) {
        for (int i = lo[0]; i <= hi[0]; ++i) {
          buckets_[flat_index(i, j, k)].push_back(cell_index);
        }
      }
    }
  }

  const std::vector<std::size_t>& bucket(double x, double y, double z) const
  {
    auto index = [&](double value, int axis) {
      double normalized = (value - low_[axis]) / (high_[axis] - low_[axis]);
      int idx = static_cast<int>(normalized * div_[axis]);
      idx = std::clamp(idx, 0, div_[axis] - 1);
      return idx;
    };
    int i = index(x, 0);
    int j = index(y, 1);
    int k = index(z, 2);
    return buckets_[flat_index(i, j, k)];
  }

private:
  std::size_t flat_index(int i, int j, int k) const
  {
    return static_cast<std::size_t>((k * div_[1] + j) * div_[0] + i);
  }

  std::array<double, 3> low_;
  std::array<double, 3> high_;
  std::array<int, 3> div_;
  std::vector<std::vector<std::size_t>> buckets_;
};

int main()
{
  std::vector<CellBounds> cells{
    {"fuel", {{-0.4, -0.4, 0.0}, {0.4, 0.4, 10.0}}},
    {"moderator", {{-0.6, -0.6, 0.0}, {0.6, 0.6, 10.0}}},
    {"guide_tube", {{-0.1, -0.1, 0.0}, {0.1, 0.1, 10.0}}},
  };

  UniformGrid grid{{-1.0, -1.0, 0.0}, {1.0, 1.0, 10.0}, {4, 4, 2}};
  for (std::size_t i = 0; i < cells.size(); ++i) {
    grid.insert(i, cells[i].bounds);
  }

  std::vector<std::array<double, 3>> queries{{0.05, 0.05, 5.0}, {0.5, 0.5, 5.0}};
  for (const auto& q : queries) {
    const auto& bucket = grid.bucket(q[0], q[1], q[2]);
    std::cout << "Query (" << q[0] << ',' << q[1] << ',' << q[2] << ") candidates:";
    for (std::size_t idx : bucket) {
      if (cells[idx].bounds.contains(q[0], q[1], q[2])) {
        std::cout << ' ' << cells[idx].name;
      }
    }
    std::cout << '\n';
  }
}
