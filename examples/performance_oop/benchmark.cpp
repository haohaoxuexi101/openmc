#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace detail {

constexpr int kParticleCount = 300'000;
constexpr int kSteps = 150;
constexpr double kStepLength = 1e-3;

struct GeometryState {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double ux{1.0};
  double uy{0.0};
  double uz{0.0};

  void move(double distance) noexcept {
    x += ux * distance;
    y += uy * distance;
    z += uz * distance;
  }
};

struct ParticleData : public GeometryState {
  double energy{1.0};
  double weight{1.0};

  void collide(double microscopic_cs) noexcept {
    // 简化的“碰撞”：能量衰减，权重略有调整
    energy *= std::exp(-microscopic_cs);
    weight *= 1.0 - 0.25 * microscopic_cs;
  }
};

struct Particle : public ParticleData {
  void transport(double step, double cs) noexcept {
    move(step);
    collide(cs);
  }
};

struct AggregateGeometryState {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double ux{1.0};
  double uy{0.0};
  double uz{0.0};
};

struct AggregateParticleData {
  double energy{1.0};
  double weight{1.0};
};

struct AggregateParticle {
  AggregateGeometryState* geom{nullptr};
  AggregateParticleData* data{nullptr};

  void transport(double step, double cs) noexcept {
    geom->x += geom->ux * step;
    geom->y += geom->uy * step;
    geom->z += geom->uz * step;
    data->energy *= std::exp(-cs);
    data->weight *= 1.0 - 0.25 * cs;
  }
};

struct SoAParticleView {
  double* x;
  double* y;
  double* z;
  double* ux;
  double* uy;
  double* uz;
  double* energy;
  double* weight;

  void transport(std::size_t i, double step, double cs) const noexcept {
    x[i] += ux[i] * step;
    y[i] += uy[i] * step;
    z[i] += uz[i] * step;
    energy[i] *= std::exp(-cs);
    weight[i] *= 1.0 - 0.25 * cs;
  }
};

struct RunResult {
  std::string name;
  double milliseconds;
  double energy_sum;
};

inline double compute_energy_sum(const std::vector<Particle>& particles) {
  return std::accumulate(particles.begin(), particles.end(), 0.0,
                         [](double acc, const Particle& p) {
                           return acc + p.energy * p.weight;
                         });
}

inline double compute_energy_sum(const std::vector<std::unique_ptr<AggregateParticleData>>& data) {
  return std::accumulate(data.begin(), data.end(), 0.0,
                         [](double acc, const std::unique_ptr<AggregateParticleData>& d) {
                           return acc + d->energy * d->weight;
                         });
}

inline double compute_energy_sum(const std::vector<double>& energy,
                                 const std::vector<double>& weight) {
  double acc = 0.0;
  for (std::size_t i = 0; i < energy.size(); ++i) {
    acc += energy[i] * weight[i];
  }
  return acc;
}

template <typename F>
RunResult time_scenario(const std::string& name, F&& f) {
  const auto start = std::chrono::steady_clock::now();
  const double energy = f();
  const auto end = std::chrono::steady_clock::now();
  const double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
  return {name, elapsed, energy};
}

inline std::vector<double> build_cross_sections() {
  std::vector<double> xs(kSteps);
  std::mt19937 rng(1337);
  std::uniform_real_distribution<double> dist(0.2, 0.6);
  for (double& v : xs) {
    v = dist(rng);
  }
  return xs;
}

}  // namespace detail

int main() {
  using namespace detail;

  const auto cross_sections = build_cross_sections();

  auto polymorphic = time_scenario("polymorphic", [&]() {
    std::vector<Particle> particles(kParticleCount);
    for (int step = 0; step < kSteps; ++step) {
      double cs = cross_sections[step];
      for (Particle& p : particles) {
        p.transport(kStepLength, cs);
      }
    }
    return compute_energy_sum(particles);
  });

  auto aggregate = time_scenario("aggregate", [&]() {
    std::vector<std::unique_ptr<AggregateGeometryState>> geometry(kParticleCount);
    std::vector<std::unique_ptr<AggregateParticleData>> data(kParticleCount);
    std::vector<AggregateParticle> particles(kParticleCount);

    for (int i = 0; i < kParticleCount; ++i) {
      geometry[i] = std::make_unique<AggregateGeometryState>();
      data[i] = std::make_unique<AggregateParticleData>();
      particles[i].geom = geometry[i].get();
      particles[i].data = data[i].get();
    }

    for (int step = 0; step < kSteps; ++step) {
      double cs = cross_sections[step];
      for (AggregateParticle& p : particles) {
        p.transport(kStepLength, cs);
      }
    }
    return compute_energy_sum(data);
  });

  auto soa = time_scenario("soa", [&]() {
    std::vector<double> x(kParticleCount, 0.0);
    std::vector<double> y(kParticleCount, 0.0);
    std::vector<double> z(kParticleCount, 0.0);
    std::vector<double> ux(kParticleCount, 1.0);
    std::vector<double> uy(kParticleCount, 0.0);
    std::vector<double> uz(kParticleCount, 0.0);
    std::vector<double> energy(kParticleCount, 1.0);
    std::vector<double> weight(kParticleCount, 1.0);

    SoAParticleView view{x.data(), y.data(), z.data(), ux.data(), uy.data(),
                         uz.data(), energy.data(), weight.data()};

    for (int step = 0; step < kSteps; ++step) {
      double cs = cross_sections[step];
      for (int i = 0; i < kParticleCount; ++i) {
        view.transport(i, kStepLength, cs);
      }
    }
    return compute_energy_sum(energy, weight);
  });

  const std::array<RunResult, 3> runs{polymorphic, aggregate, soa};

  for (const auto& run : runs) {
    std::cout << "Scenario: " << std::left << std::setw(12) << run.name
              << " | Time: " << std::setw(8) << std::fixed << std::setprecision(2)
              << run.milliseconds << " ms | Final energy sum: " << std::scientific
              << std::setprecision(6) << run.energy_sum << '\n';
  }

  const double reference = polymorphic.energy_sum;
  const bool matches = std::all_of(runs.begin(), runs.end(), [&](const RunResult& r) {
    return std::abs(r.energy_sum - reference) < 1e-6 * std::abs(reference);
  });

  if (matches) {
    std::cout << "All strategies produce equivalent physics results." << std::endl;
    return 0;
  }

  std::cerr << "Mismatch detected between strategies." << std::endl;
  return 1;
}
