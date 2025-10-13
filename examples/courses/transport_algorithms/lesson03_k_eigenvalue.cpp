#include <iostream>
#include <vector>
#include <numeric>

// 课程 3：k-eigenvalue 幂迭代
// 本示例通过简单的两能群模型展示特征值迭代的主要步骤。

struct GroupData {
  double nu_sigma_f;
  double sigma_a;
};

int main() {
  std::vector<GroupData> groups = {{0.4, 0.2}, {0.1, 0.05}};

  double k = 1.0;
  std::vector<double> flux = {1.0, 1.0};

  for (int iter = 0; iter < 5; ++iter) {
    // 1) 生成裂变源
    double production = 0.0;
    for (size_t g = 0; g < groups.size(); ++g) {
      production += groups[g].nu_sigma_f * flux[g];
    }

    // 2) 归一化得到新的通量
    for (size_t g = 0; g < flux.size(); ++g) {
      flux[g] = groups[g].nu_sigma_f / groups[g].sigma_a;
    }

    double flux_sum = std::accumulate(flux.begin(), flux.end(), 0.0);
    for (double& val : flux) val /= flux_sum;

    double new_k = production;
    std::cout << "迭代 " << iter << " k = " << new_k << " 通量 = [" << flux[0]
              << ", " << flux[1] << "]\n";
    k = new_k;
  }

  std::cout << "最终 k ≈ " << k << "\n";
  return 0;
}
