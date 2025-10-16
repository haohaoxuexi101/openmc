#include <array>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

// lesson22_particle_data.cpp
// ---------------------------------
// 模拟 OpenMC::ParticleData 对 GeometryState 的继承：
// * 粒子数据继续复用坐标栈、失效标记等工具；
// * 在其上叠加截面缓存、材料索引、事件统计等字段；
// * 利用 std::optional 表达“需要懒加载”的昂贵数据。
// 演示粒子进入不同材料时如何更新缓存，体现组合 + 缓存友好布局。

class GeometryState {
public:
  using Position = std::array<double, 3>;

  void reset(Position p)
  {
    position_world_ = p;
    on_surface_ = -1;
    lost_ = false;
  }

  void move(double d, std::array<double, 3> u)
  {
    for (int i = 0; i < 3; ++i) {
      position_world_[i] += d * u[i];
    }
  }

  void mark_on_surface(int surface) { on_surface_ = surface; }

  void mark_lost(std::string msg)
  {
    lost_ = true;
    lost_msg_ = std::move(msg);
  }

  void print_geometry() const
  {
    std::cout << "pos=(" << position_world_[0] << ',' << position_world_[1] << ','
              << position_world_[2] << ")";
    if (on_surface_ >= 0) std::cout << " surface=" << on_surface_;
    if (lost_) std::cout << " LOST:" << lost_msg_;
    std::cout << '\n';
  }

private:
  Position position_world_ {0.0, 0.0, 0.0};
  int on_surface_ {-1};
  bool lost_ {false};
  std::string lost_msg_;
};

struct CrossSectionCache {
  int material_id {-1};
  std::optional<double> total_xs;
  std::optional<double> absorption_xs;
};

class ParticleData : public GeometryState {
public:
  void set_material(int material_id)
  {
    if (material_id != cache_.material_id) {
      cache_.material_id = material_id;
      cache_.total_xs.reset();
      cache_.absorption_xs.reset();
      events_.push_back("material->" + std::to_string(material_id));
    }
  }

  void set_energy(double e)
  {
    if (energy_ != e) {
      energy_ = e;
      cache_.total_xs.reset();
      cache_.absorption_xs.reset();
      events_.push_back("energy->" + std::to_string(e));
    }
  }

  double total_xs()
  {
    if (!cache_.total_xs) {
      cache_.total_xs = 0.01 * (cache_.material_id + 1) * energy_;
      events_.push_back("compute total_xs");
    }
    return *cache_.total_xs;
  }

  double absorption_xs()
  {
    if (!cache_.absorption_xs) {
      cache_.absorption_xs = total_xs() * 0.2;
      events_.push_back("compute absorption_xs");
    }
    return *cache_.absorption_xs;
  }

  void print_state() const
  {
    print_geometry();
    std::cout << "  energy=" << energy_ << " material=" << cache_.material_id
              << " cache(total=" << (cache_.total_xs ? "Y" : "N")
              << ",abs=" << (cache_.absorption_xs ? "Y" : "N") << ")\n";
    std::cout << "  events:";
    for (const auto& e : events_) std::cout << ' ' << e;
    std::cout << '\n';
  }

private:
  double energy_ {1.0};
  CrossSectionCache cache_ {};
  std::vector<std::string> events_;
};

int main()
{
  ParticleData p;
  p.reset({0.0, 0.0, 0.0});
  p.set_material(3);
  p.set_energy(2.0);
  p.total_xs();
  p.absorption_xs();
  p.print_state();
  p.move(1.0, {0.0, 0.0, 1.0});
  p.set_material(5);
  p.total_xs();
  p.print_state();
}
