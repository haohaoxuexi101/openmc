#include <chrono>
#include <iostream>
#include <optional>
#include <random>

struct ParticleState {
  double energy{1.0e5};
  double material_temperature{600.0};
  double total_xs{0.0};
  double collision_time{0.0};
};

class CrossSectionCache {
public:
  bool try_get(double energy, double temperature, double& xs) const
  {
    if (cached_energy_ && cached_temperature_ && *cached_energy_ == energy &&
        *cached_temperature_ == temperature) {
      xs = cached_xs_;
      return true;
    }
    return false;
  }

  void store(double energy, double temperature, double xs)
  {
    cached_energy_ = energy;
    cached_temperature_ = temperature;
    cached_xs_ = xs;
  }

private:
  std::optional<double> cached_energy_;
  std::optional<double> cached_temperature_;
  double cached_xs_{0.0};
};

class ParticlePipeline {
public:
  explicit ParticlePipeline(unsigned seed) : rng_{seed}, exp_(1.0) {}

  void event_calculate_xs(ParticleState& state)
  {
    double xs;
    if (cache_.try_get(state.energy, state.material_temperature, xs)) {
      std::cout << "cache hit\n";
    } else {
      xs = 0.02 * state.energy / 1.0e5 + 0.001 * (state.material_temperature - 600.0);
      cache_.store(state.energy, state.material_temperature, xs);
      std::cout << "cache miss -> xs=" << xs << '\n';
    }
    state.total_xs = xs;
  }

  void event_sample_free_flight(ParticleState& state)
  {
    const double distance = exp_(rng_) / state.total_xs;
    state.collision_time += distance / 3.0e7; // 模拟光速的一部分
    std::cout << "sampled distance=" << distance << '\n';
  }

  void event_scatter(ParticleState& state)
  {
    std::uniform_real_distribution<double> uniform(0.5, 1.0);
    state.energy *= uniform(rng_);
    std::cout << "new energy=" << state.energy << '\n';
  }

private:
  CrossSectionCache cache_;
  std::mt19937 rng_;
  std::exponential_distribution<double> exp_;
};

int main()
{
  ParticleState state;
  ParticlePipeline pipeline{2024};

  for (int i = 0; i < 2; ++i) {
    pipeline.event_calculate_xs(state);
    pipeline.event_sample_free_flight(state);
    pipeline.event_scatter(state);
  }
}
