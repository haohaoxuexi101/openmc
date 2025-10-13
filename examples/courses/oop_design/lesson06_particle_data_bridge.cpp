#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

// 课程 6：ParticleData 接口桥接示例
// -----------------------------------
// 本课聚焦 OpenMC `particle_data.h` 的继承关系：ParticleData 继承自 GeometryState。
// 我们用一个精简示例说明为何要让几何状态作为基类，粒子状态作为派生类，
// 并演示如何封装 getter/setter 以保证传输环节的性能。

// ------------------ 几何状态基类（对应 geometry.h） ------------------
class GeometryState {
public:
  struct Coord {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
  };

  GeometryState() = default;
  explicit GeometryState(Coord r) : r_(r) {}

  const Coord& coord() const { return r_; }
  void set_coord(double x, double y, double z) { r_ = {x, y, z}; }

  int cell() const { return cell_; }
  void set_cell(int c) { cell_ = c; }

  int surface() const { return surface_; }
  void set_surface(int s) { surface_ = s; }

protected:
  Coord r_{};        // 位置坐标
  int cell_ = -1;    // 当前所在 cell ID
  int surface_ = -1; // 最近穿越的 surface ID
};

// ------------------ 粒子状态派生类（对应 particle_data.h） ------------------
class ParticleDataLite : public GeometryState {
public:
  ParticleDataLite() = default;

  ParticleDataLite(uint64_t id, Coord r, std::array<uint64_t, 4> seeds)
      : GeometryState(r), id_(id), seeds_(seeds) {}

  uint64_t id() const { return id_; }
  void set_id(uint64_t id) { id_ = id; }

  double energy() const { return E_; }
  void set_energy(double E) { E_ = E; }

  double weight() const { return wgt_; }
  void set_weight(double w) { wgt_ = w; }

  std::array<uint64_t, 4>& seeds() { return seeds_; }
  const std::array<uint64_t, 4>& seeds() const { return seeds_; }

  void push_path(double segment) { path_length_ += segment; }
  double path_length() const { return path_length_; }

  void score_absorption(double macro_xs) { absorption_score_ += wgt_ * macro_xs; }
  double absorption_score() const { return absorption_score_; }

private:
  uint64_t id_ = 0;                        // 粒子唯一 ID
  double E_ = 1.0;                         // 能量 (MeV)
  double wgt_ = 1.0;                       // 权重
  std::array<uint64_t, 4> seeds_{{1, 2, 3, 4}}; // 多流随机种子
  double path_length_ = 0.0;               // 轨迹长度累积
  double absorption_score_ = 0.0;          // 吸收计分缓存
};

// ------------------ 示例流程 ------------------
int main() {
  // 创建粒子并设置初始几何状态
  ParticleDataLite particle(42, {0.0, 0.0, 0.0}, {2024, 6, 1, 7});
  particle.set_cell(10);
  particle.set_surface(3);
  particle.set_energy(2.0); // MeV

  // 模拟一次自由程与吸收计分
  particle.push_path(1.25);      // cm
  particle.score_absorption(0.4);

  // 输出粒子状态，观察 GeometryState 字段与 ParticleData 字段的结合
  std::cout << "粒子ID: " << particle.id() << "\n";
  std::cout << "位置: (" << particle.coord().x << ", " << particle.coord().y << ", "
            << particle.coord().z << ")\n";
  std::cout << "当前Cell ID: " << particle.cell() << ", 最近表面ID: " << particle.surface() << "\n";
  std::cout << "能量: " << particle.energy() << " MeV, 权重: " << particle.weight() << "\n";
  std::cout << "轨迹长度累计: " << particle.path_length() << " cm\n";
  std::cout << "吸收计分累计: " << particle.absorption_score() << "\n";

  // 演示随机种子与多线程友好性
  std::cout << "随机种子流: ";
  for (auto seed : particle.seeds()) {
    std::cout << seed << ' ';
  }
  std::cout << "\n";

  return 0;
}

