#include <array>
#include <iostream>
#include <string>
#include <vector>

// lesson21_geometry_state.cpp
// ---------------------------------
// 以可执行示例还原 OpenMC::GeometryState 的职责：
// 1. 为几何遍历缓存一条“坐标栈”，避免在不同模块间重复计算局部坐标。
// 2. 提供事件记录（上一碰撞点、当前边界）与异常标记接口。
// 3. 与粒子状态解耦，允许绘图、体素剖分等工具直接操控几何层级。
// 下文的实现保持结构紧凑，让读者理解组合 + 轻量面向对象的优势。

using Position = std::array<double, 3>;

struct LevelInfo {
  int cell_id;
  int universe_id;
  std::array<double, 3> local_pos;
};

class GeometryState {
public:
  GeometryState() = default;

  void reset(Position new_pos)
  {
    levels_.clear();
    position_world_ = new_pos;
    last_collision_ = new_pos;
    on_surface_ = -1;
  }

  void push_level(int cell, int universe, std::array<double, 3> local)
  {
    levels_.push_back(LevelInfo{cell, universe, local});
  }

  void pop_level()
  {
    if (!levels_.empty()) levels_.pop_back();
  }

  void record_collision()
  {
    last_collision_ = position_world_;
    on_surface_ = -1;
  }

  void advance(double distance, std::array<double, 3> direction)
  {
    for (int i = 0; i < 3; ++i) {
      position_world_[i] += distance * direction[i];
    }
  }

  void mark_on_surface(int surface_id)
  {
    on_surface_ = surface_id;
  }

  [[nodiscard]] bool lost() const { return lost_; }

  void mark_lost(std::string reason)
  {
    lost_ = true;
    message_ = std::move(reason);
  }

  void print_state() const
  {
    std::cout << "world=(" << position_world_[0] << ',' << position_world_[1] << ','
              << position_world_[2] << ") levels=" << levels_.size();
    if (on_surface_ >= 0) {
      std::cout << " on_surface=" << on_surface_;
    }
    std::cout << '\n';
    for (size_t i = 0; i < levels_.size(); ++i) {
      const auto& lvl = levels_[i];
      std::cout << "  level " << i << " cell=" << lvl.cell_id << " universe="
                << lvl.universe_id << " local=(" << lvl.local_pos[0] << ','
                << lvl.local_pos[1] << ',' << lvl.local_pos[2] << ")\n";
    }
    if (lost_) {
      std::cout << "  LOST: " << message_ << '\n';
    }
  }

private:
  Position position_world_ {0.0, 0.0, 0.0};
  Position last_collision_ {0.0, 0.0, 0.0};
  int on_surface_ {-1};
  bool lost_ {false};
  std::string message_;
  std::vector<LevelInfo> levels_;
};

int main()
{
  GeometryState state;
  state.reset({0.0, 0.0, 0.0});
  state.push_level(10, 1, {0.0, 0.0, 0.0});
  state.advance(2.0, {0.0, 0.0, 1.0});
  state.push_level(42, 7, {0.0, 0.0, 0.2});
  state.mark_on_surface(99);
  state.print_state();
  state.record_collision();
  state.advance(1.0, {1.0, 0.0, 0.0});
  state.pop_level();
  state.mark_lost("cell lookup failed after crossing surface 99");
  state.print_state();
}
