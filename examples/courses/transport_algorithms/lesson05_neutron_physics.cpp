#include <iostream>
#include <vector>

// 课程 5：多群中子物理简表
// 展示如何把散射矩阵、裂变产额等核数据组织成面向对象结构，
// 便于在求解器中调用。

struct GroupPhysics {
  double sigma_t;
  double sigma_s;
  double nu_sigma_f;
};

class MultiGroupSet {
public:
  explicit MultiGroupSet(std::vector<GroupPhysics> groups)
      : groups_(std::move(groups)) {}

  double absorption(size_t g) const { return groups_[g].sigma_t - groups_[g].sigma_s; }
  double reproduction(size_t g) const { return groups_[g].nu_sigma_f; }

  void dump() const {
    for (size_t g = 0; g < groups_.size(); ++g) {
      std::cout << "能群 " << g << " 吸收 = " << absorption(g)
                << " 裂变产额 = " << reproduction(g) << "\n";
    }
  }

private:
  std::vector<GroupPhysics> groups_;
};

int main() {
  MultiGroupSet mg({{1.0, 0.6, 0.4}, {0.8, 0.7, 0.1}});
  mg.dump();
  return 0;
}
