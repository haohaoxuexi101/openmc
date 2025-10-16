#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

// 粒子状态包含位置、方向与存活标记。
struct Particle {
  double x{};
  double y{};
  double z{};
  double mu{}; // cos(theta) 相对于法向
  bool alive{true};
};

class BoundaryCondition {
public:
  virtual ~BoundaryCondition() = default;
  virtual void handle(Particle& p) const = 0;
  virtual std::string name() const = 0;
};

class VacuumBC final : public BoundaryCondition {
public:
  void handle(Particle& p) const override { p.alive = false; }
  std::string name() const override { return "Vacuum"; }
};

class ReflectiveBC final : public BoundaryCondition {
public:
  void handle(Particle& p) const override
  {
    p.mu = -p.mu;
    p.z -= 1e-6; // 推回边界内侧，避免无限循环。
  }

  std::string name() const override { return "Reflective"; }
};

class WhiteBC final : public BoundaryCondition {
public:
  explicit WhiteBC(unsigned seed) : rng_{seed}, dist_{0.0, 1.0} {}

  void handle(Particle& p) const override
  {
    // 模拟各向同性漫反射，采样新方向余弦。
    double xi = dist_(rng_);
    p.mu = std::sqrt(xi);
    p.z -= 1e-6;
  }

  std::string name() const override { return "White"; }

private:
  mutable std::mt19937 rng_;
  mutable std::uniform_real_distribution<double> dist_;
};

class PeriodicBC final : public BoundaryCondition {
public:
  PeriodicBC(double pitch) : pitch_{pitch} {}

  void handle(Particle& p) const override
  {
    p.z -= pitch_;
  }

  std::string name() const override { return "Periodic"; }

private:
  double pitch_;
};

void simulate(const BoundaryCondition& bc)
{
  Particle particle{0.0, 0.0, 1.0, -0.8, true};
  std::cout << "Boundary: " << bc.name() << '\n';
  bc.handle(particle);
  std::cout << "  alive=" << particle.alive << ", mu=" << particle.mu
            << ", z=" << particle.z << '\n';
}

int main()
{
  std::vector<std::unique_ptr<BoundaryCondition>> bcs;
  bcs.emplace_back(std::make_unique<VacuumBC>());
  bcs.emplace_back(std::make_unique<ReflectiveBC>());
  bcs.emplace_back(std::make_unique<WhiteBC>(1234));
  bcs.emplace_back(std::make_unique<PeriodicBC>(2.0));

  for (const auto& bc : bcs) {
    simulate(*bc);
  }
}
