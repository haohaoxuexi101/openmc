#include <cstdint>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// This example is inspired by the OpenMC Monte Carlo code base.  The goal is to
// demonstrate how the GeometryState, ParticleData, and Particle classes relate
// to each other in a simplified, self-contained setting.  The real OpenMC code
// contains many more data members and algorithms, but the class hierarchy is
// the same: GeometryState manages spatial information, ParticleData adds
// particle-specific attributes, and Particle supplies the high-level behaviour
// for transporting a particle through space.
// -----------------------------------------------------------------------------

// GeometryState ---------------------------------------------------------------
// Keeps track of coordinates and allows the particle to be moved around.  The
// real OpenMC class stores arrays of nested coordinates, but here we only keep a
// single position and direction so the example stays compact.
class GeometryState {
public:
  GeometryState() = default;

  void init_from_position_direction(const std::vector<double>& position,
    const std::vector<double>& direction)
  {
    if (position.size() != 3 || direction.size() != 3) {
      throw std::invalid_argument("Position and direction must have size 3");
    }
    position_ = position;
    direction_ = direction;
    history_.clear();
    history_.push_back(position_);
  }

  virtual void mark_as_lost(const std::string& message)
  {
    std::ostringstream out;
    out << "[GeometryState] Particle " << id_ << " lost: " << message;
    std::cout << out.str() << '\n';
  }

  void move_distance(double distance)
  {
    for (std::size_t i = 0; i < 3; ++i) {
      position_[i] += direction_[i] * distance;
    }
    history_.push_back(position_);
  }

  void set_direction(const std::vector<double>& new_direction)
  {
    if (new_direction.size() != 3) {
      throw std::invalid_argument("Direction must have size 3");
    }
    direction_ = new_direction;
  }

  const std::vector<double>& position() const { return position_; }
  const std::vector<double>& direction() const { return direction_; }

  void set_id(std::int64_t id) { id_ = id; }
  std::int64_t id() const { return id_; }

  const std::vector<std::vector<double>>& history() const { return history_; }

protected:
  std::vector<double> position_ {0.0, 0.0, 0.0};
  std::vector<double> direction_ {0.0, 0.0, 1.0};
  std::int64_t id_ {0};
  std::vector<std::vector<double>> history_;
};

// ParticleData ----------------------------------------------------------------
// Adds physical quantities such as energy, weight, and simulation time.  In the
// real code this class also stores large caches of cross sections and random
// number seeds.  Here we keep only a few essential members to highlight the
// idea of extending a base class with additional responsibilities.
class ParticleData : public GeometryState {
public:
  ParticleData() = default;

  void set_energy(double energy_eV) { energy_eV_ = energy_eV; }
  double energy() const { return energy_eV_; }

  void set_weight(double weight) { weight_ = weight; }
  double weight() const { return weight_; }

  void set_time(double time_s) { time_ = time_s; }
  double time() const { return time_; }

  void record_collision()
  {
    ++collision_count_;
    last_collision_energy_eV_ = energy_eV_;
  }

  int collisions() const { return collision_count_; }
  double last_collision_energy() const { return last_collision_energy_eV_; }

private:
  double energy_eV_ {1.0e6};
  double weight_ {1.0};
  double time_ {0.0};
  int collision_count_ {0};
  double last_collision_energy_eV_ {0.0};
};

// Particle --------------------------------------------------------------------
// Provides behaviour using the state stored in ParticleData and GeometryState.
// We implement a few simple actions: streaming (moving forward), colliding, and
// printing diagnostic information.
class Particle : public ParticleData {
public:
  void stream(double distance)
  {
    std::cout << "[Particle] Streaming for " << distance << " cm\n";
    move_distance(distance);
    set_time(time() + distance / speed_of_light_cm_per_s());
  }

  void collide(double energy_loss_fraction)
  {
    if (energy_loss_fraction < 0.0 || energy_loss_fraction > 1.0) {
      throw std::invalid_argument("Energy loss fraction must be between 0 and 1");
    }

    double new_energy = energy() * (1.0 - energy_loss_fraction);
    std::cout << "[Particle] Collision: energy " << energy() << " -> "
              << new_energy << " eV\n";
    set_energy(new_energy);
    record_collision();
  }

  void print_summary() const
  {
    std::cout << "\n=== Particle summary ===\n";
    std::cout << "ID: " << id() << '\n';
    std::cout << "Position: (" << format_vector(position()) << ")\n";
    std::cout << "Direction: (" << format_vector(direction()) << ")\n";
    std::cout << "Energy: " << energy() << " eV\n";
    std::cout << "Weight: " << weight() << '\n';
    std::cout << "Simulation time: " << time() << " s\n";
    std::cout << "Collisions: " << collisions() << '\n';
    std::cout << "Last collision energy: " << last_collision_energy() << " eV\n";
    std::cout << "Trajectory points:\n";
    const auto& hist = history();
    for (std::size_t i = 0; i < hist.size(); ++i) {
      std::cout << "  Step " << std::setw(2) << i << ": (" << format_vector(hist[i])
                << ")\n";
    }
  }

private:
  static double speed_of_light_cm_per_s() { return 2.99792458e10; }

  static std::string format_vector(const std::vector<double>& v)
  {
    std::ostringstream out;
    out << std::fixed << std::setprecision(2)
        << v[0] << ", " << v[1] << ", " << v[2];
    return out.str();
  }
};

// Demonstration ---------------------------------------------------------------
int main()
{
  try {
    Particle p;
    p.set_id(42);
    p.init_from_position_direction({0.0, 0.0, 0.0}, {0.0, 0.0, 1.0});
    p.set_energy(2.0e6);
    p.set_weight(0.8);

    p.stream(10.0);               // Move forward by 10 cm
    p.collide(0.25);              // Lose 25 % of the energy
    p.set_direction({0.0, 1.0, 0.0});
    p.stream(5.0);                // Turn and move sideways
    p.collide(0.10);              // Another collision

    p.print_summary();
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << '\n';
    return 1;
  }

  return 0;
}
