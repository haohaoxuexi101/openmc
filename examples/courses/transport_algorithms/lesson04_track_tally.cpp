#include <iostream>
#include <vector>

// 课程 4：轨迹长度计数
// 展示如何积累 track-length tally，并与单元路径关联。

struct Segment {
  int cell_id;
  double length;
};

int main() {
  std::vector<Segment> history = {{0, 0.5}, {1, 1.2}, {0, 0.3}};
  std::vector<double> tally(2, 0.0);

  for (const auto& seg : history) {
    tally[seg.cell_id] += seg.length;
  }

  std::cout << "Cell0 轨迹长度 = " << tally[0] << "\n";
  std::cout << "Cell1 轨迹长度 = " << tally[1] << "\n";
  return 0;
}
