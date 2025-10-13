#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// 课程 1：域与材料配置解析
// --------------------------
// 本课展示 random ray 求解器如何读取输入配置并构造“材料-区域”映射。
// 为了聚焦结构，我们手工构建一个 2×2×1 的盒体装配：燃料、慢化剂、反射层。
// 代码模仿 OpenMC 的做法，使用轻量结构体存储多群截面与几何包围盒。

struct MaterialXS {
  std::string name;
  std::vector<double> sigma_t; // 总截面
  std::vector<double> sigma_s; // 散射截面（各群自散射）
  std::vector<double> nu_sigma_f; // 裂变截面
};

struct AxisAlignedBox {
  double xmin, xmax;
  double ymin, ymax;
  double zmin, zmax;
};

struct RegionConfig {
  std::string label;
  int material_id;
  AxisAlignedBox bounds;
};

int main() {
  std::vector<MaterialXS> materials = {
      {"燃料", {0.6, 0.5, 0.4}, {0.2, 0.25, 0.3}, {0.3, 0.25, 0.15}},
      {"慢化剂", {0.2, 0.15, 0.1}, {0.18, 0.13, 0.08}, {0.0, 0.0, 0.0}},
      {"反射层", {0.3, 0.28, 0.22}, {0.25, 0.23, 0.2}, {0.0, 0.0, 0.0}},
  };

  std::vector<RegionConfig> regions = {
      {"中部燃料", 0, {-5.0, 5.0, -5.0, 5.0, -2.0, 2.0}},
      {"上方慢化剂", 1, {-5.0, 5.0, -5.0, 5.0, 2.0, 6.0}},
      {"下方反射层", 2, {-5.0, 5.0, -5.0, 5.0, -6.0, -2.0}},
  };

  std::cout << "==== 材料表 (三群) ====" << '\n';
  for (size_t i = 0; i < materials.size(); ++i) {
    const auto& m = materials[i];
    std::cout << i << ": " << m.name << '\n';
    for (size_t g = 0; g < m.sigma_t.size(); ++g) {
      std::cout << "  群" << g << ": Σ_t=" << m.sigma_t[g]
                << ", Σ_s=" << m.sigma_s[g]
                << ", νΣ_f=" << m.nu_sigma_f[g] << '\n';
    }
  }

  std::cout << "\n==== 几何区域 ====" << '\n';
  for (const auto& r : regions) {
    std::cout << r.label << " 使用材料 #" << r.material_id << "，包围盒："
              << "x=[" << r.bounds.xmin << ", " << r.bounds.xmax << "] cm"
              << " y=[" << r.bounds.ymin << ", " << r.bounds.ymax << "] cm"
              << " z=[" << r.bounds.zmin << ", " << r.bounds.zmax << "] cm" << '\n';
  }

  return 0;
}

