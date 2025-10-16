#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <vector>

// lesson23_particle_lifecycle.cpp
// ---------------------------------
// 将 GeometryState 与 ParticleData 组合成可推进的 Particle：
// * GeometryState 负责几何（位置、表面、坐标栈）。
// * ParticleData 负责物理参数（能量、材料、截面缓存）。
// * Particle 提供事件驱动接口，封装推进、碰撞与次级粒子生成统计。
// 这样的三层设计正是 OpenMC 中避免“巨型 God object”又保持缓存亲和的核心技巧。

using Vec3 = std::array<double, 3>;

class GeometryState {
public:
  void init(Vec3 r, Vec3 u)
  {
    r_ = r;
    u_ = u;
    surface_ = -1;
    history_.clear();
  }

  void move(double distance)
  {
    for (int i = 0; i < 3; ++i) {
      r_[i] += distance * u_[i];
    }
    history_.push_back({r_, surface_});
  }

  void set_direction(Vec3 dir) { u_ = dir; }

  void mark_surface(int surface) { surface_ = surface; }

  const Vec3& position() const { return r_; }
  const Vec3& direction() const { return u_; }

  void print_geometry() const
  {
    std::cout << "pos=(" << r_[0] << ',' << r_[1] << ',' << r_[2] << ")"
              << " dir=(" << u_[0] << ',' << u_[1] << ',' << u_[2] << ")";
    if (surface_ >= 0) std::cout << " surface=" << surface_;
    std::cout << '\n';
    std::cout << "  history:";
    for (const auto& entry : history_) {
      std::cout << " [" << entry.surface << ":" << entry.pos[2] << ']';
    }
    std::cout << '\n';
  }

private:
  struct HistoryEntry {
    Vec3 pos;
    int surface;
  };

  Vec3 r_ {0.0, 0.0, 0.0};
  Vec3 u_ {0.0, 0.0, 1.0};
  int surface_ {-1};
  std::vector<HistoryEntry> history_;
};

struct Material {
  std::string name;
  double total_coeff;
  double scatter_fraction;
};

class ParticleData : public GeometryState {
public:
  void set_material(const Material* material)
  {
    if (material_ != material) {
      material_ = material;
      cache_total_.reset();
      cache_absorption_.reset();
      log_.push_back("set material -> " + material->name);
    }
  }

  void set_energy(double e)
  {
    if (energy_ != e) {
      energy_ = e;
      cache_total_.reset();
      cache_absorption_.reset();
      log_.push_back("set energy -> " + std::to_string(e));
    }
  }

protected:
  double total_xs()
  {
    if (!cache_total_) {
      cache_total_ = material_->total_coeff * std::sqrt(energy_);
      log_.push_back("compute total_xs");
    }
    return *cache_total_;
  }

  double absorption_xs()
  {
    if (!cache_absorption_) {
      cache_absorption_ = total_xs() * (1.0 - material_->scatter_fraction);
      log_.push_back("compute absorption_xs");
    }
    return *cache_absorption_;
  }

  void dump_log() const
  {
    std::cout << "  cache log:";
    for (const auto& item : log_) std::cout << ' ' << item;
    std::cout << '\n';
  }

  const Material* material_ {nullptr};
  double energy_ {1.0};
  std::optional<double> cache_total_;
  std::optional<double> cache_absorption_;
  std::vector<std::string> log_;
};

class Particle : public ParticleData {
public:
  Particle(int64_t id, const Material* material, Vec3 r, Vec3 u)
    : id_{id}
  {
    init(r, u);
    set_material(material);
  }

  void event_calculate_xs()
  {
    macro_total_ = total_xs();
    macro_absorption_ = absorption_xs();
  }

  void event_advance(double distance)
  {
    move(distance);
    track_length_ += distance;
  }

  void event_collision(std::mt19937& rng)
  {
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    const double xi = uni(rng);
    ++collisions_;
    if (xi < macro_absorption_ / macro_total_) {
      alive_ = false;
      log_.push_back("absorbed");
    } else {
      scatter_angle_ = 2.0 * xi - 1.0;
      set_direction({std::sqrt(1 - scatter_angle_ * scatter_angle_), 0.0, scatter_angle_});
      log_.push_back("scattered mu=" + std::to_string(scatter_angle_));
    }
  }

  void print_report() const
  {
    std::cout << "Particle " << id_ << (alive_ ? " alive" : " dead")
              << " collisions=" << collisions_ << " track=" << track_length_
              << '\n';
    print_geometry();
    dump_log();
    std::cout << "  physics log:";
    for (const auto& entry : log_) std::cout << ' ' << entry;
    std::cout << '\n';
  }

private:
  int64_t id_;
  bool alive_ {true};
  double macro_total_ {0.0};
  double macro_absorption_ {0.0};
  double track_length_ {0.0};
  int collisions_ {0};
  double scatter_angle_ {1.0};
  std::vector<std::string> log_;
};

int main()
{
  const Material fuel{"fuel", 1.2, 0.85};
  const Material moderator{"moderator", 0.5, 0.95};

  std::mt19937 rng(1337);

  Particle particle{1, &fuel, {0.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
  particle.set_energy(2.0); // 激发缓存无效化
  particle.event_calculate_xs();
  particle.event_advance(1.0);
  particle.event_collision(rng);
  particle.print_report();

  particle.set_material(&moderator); // 重置缓存并记录日志
  particle.set_energy(0.5);
  particle.event_calculate_xs();
  particle.event_advance(0.3);
  particle.event_collision(rng);
  particle.print_report();
}
