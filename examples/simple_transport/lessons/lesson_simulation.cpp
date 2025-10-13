#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// 本示例聚焦 Simulation 层的职责：协调 PhysicsDriver、Bank、Tally 三者。
// 为降低物理复杂度，我们构造一个“虚拟”物理驱动 StubPhysicsDriver，
// 但接口严格遵循主程序，从而展示 OpenMC 风格的顶层控制循环。

namespace lesson_simulation {

struct Particle {
  double weight = 1.0;
};

class Bank {
 public:
  void clear() { particles_.clear(); }
  void add(const Particle& p) { particles_.push_back(p); }
  [[nodiscard]] std::size_t size() const noexcept { return particles_.size(); }

  Particle& operator[](std::size_t i) { return particles_.at(i); }
  const Particle& operator[](std::size_t i) const { return particles_.at(i); }

  std::vector<Particle>::iterator begin() { return particles_.begin(); }
  std::vector<Particle>::iterator end() { return particles_.end(); }

 private:
  std::vector<Particle> particles_;
};

// 用于展示统计管理的最小实现：记录每代产生的“裂变粒子总权重”。
class GenerationTally {
 public:
  void record(double produced) { production_history_.push_back(produced); }
  const std::vector<double>& history() const noexcept { return production_history_; }

 private:
  std::vector<double> production_history_;
};

class StubPhysicsDriver {
 public:
  explicit StubPhysicsDriver(unsigned seed) : rng_(seed), fluctuation_(0.9, 1.1) {}

  void initialise_source(Bank& source, std::size_t histories) {
    source.clear();
    for (std::size_t i = 0; i < histories; ++i) {
      source.add(Particle{1.0});  // 每个源粒子初始权重为 1.0
    }
  }

  double transport_generation(Bank& source, Bank& fission_bank, GenerationTally& tally) {
    fission_bank.clear();
    double produced = 0.0;

    for (auto& particle : source) {
      double local_k = base_keff_ + 0.02 * static_cast<double>(generation_index_);
      double new_weight = particle.weight * local_k * fluctuation_(rng_);
      produced += new_weight;
      tally.record(new_weight);
      fission_bank.add(Particle{new_weight});
    }

    ++generation_index_;
    if (source.size() == 0) throw std::runtime_error("源粒子数量为 0");
    return produced / static_cast<double>(source.size());
  }

  void resample_source(const Bank& fission_bank, Bank& source) {
    if (fission_bank.size() == 0) {
      throw std::runtime_error("裂变库为空，无法重采样");
    }
    std::uniform_int_distribution<std::size_t> pick(0, fission_bank.size() - 1);
    for (auto& particle : source) {
      particle = fission_bank[pick(rng_)];
      particle.weight = 1.0;  // 与 OpenMC 一样，将源粒子重新归一化
    }
  }

 private:
  double base_keff_ = 1.02;  // 基准增殖因子，模拟一个略超临界系统
  std::size_t generation_index_ = 0;
  std::mt19937 rng_;
  std::uniform_real_distribution<double> fluctuation_;
};

struct SimulationSettings {
  std::size_t histories = 1000;
  std::size_t generations = 5;
  unsigned seed = 202406u;
};

class Simulation {
 public:
  Simulation(SimulationSettings settings)
      : settings_(settings), physics_(settings.seed) {
    physics_.initialise_source(source_, settings_.histories);
  }

  void run() {
    std::cout << std::fixed << std::setprecision(5);
    for (std::size_t g = 0; g < settings_.generations; ++g) {
      double keff = physics_.transport_generation(source_, fission_bank_, tally_);
      physics_.resample_source(fission_bank_, source_);
      std::cout << "第 " << g + 1 << " 代 k-effective 估计 = " << keff << '\n';
    }

    std::cout << "\n每代产生的总权重记录：\n";
    const auto& history = tally_.history();
    for (std::size_t i = 0; i < history.size(); ++i) {
      std::cout << "  event " << i + 1 << " -> " << history[i] << '\n';
    }
  }

 private:
  SimulationSettings settings_;
  StubPhysicsDriver physics_;
  GenerationTally tally_;
  Bank source_;
  Bank fission_bank_;
};

}  // namespace lesson_simulation

int main() {
  using namespace lesson_simulation;

  SimulationSettings settings;
  settings.histories = 4;     // 减少粒子数量，便于输出观察
  settings.generations = 3;   // 演示 power iteration 流程
  settings.seed = 13579u;

  Simulation simulation(settings);
  simulation.run();
}
