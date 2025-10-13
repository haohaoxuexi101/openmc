#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/position.h"
#include "openmc/surface.h"
#include "openmc/universe.h"

namespace tutorial {

// Simple helper describing how a surface bounds a cell.
struct SurfaceSense {
  int surface_id{};
  bool positive{}; // true means inside corresponds to evaluate <= 0
};

//------------------------------------------------------------------------------
// Lightweight planar surface that delegates to openmc::Surface
//------------------------------------------------------------------------------

class PlaneSurface : public openmc::Surface {
public:
  PlaneSurface(std::string name, openmc::Position normal, double d);

  double evaluate(openmc::Position r) const override;
  double distance(openmc::Position r, openmc::Direction u, bool coincident) const override;
  openmc::Direction normal(openmc::Position /*r*/) const override;

protected:
  void to_hdf5_inner(hid_t) const override {}

private:
  openmc::Position normal_;
  double d_{};
};

class SphereSurface : public openmc::Surface {
public:
  SphereSurface(std::string name, openmc::Position center, double radius);

  double evaluate(openmc::Position r) const override;
  double distance(openmc::Position r, openmc::Direction u, bool coincident) const override;
  openmc::Direction normal(openmc::Position r) const override;

protected:
  void to_hdf5_inner(hid_t) const override {}

private:
  openmc::Position center_;
  double radius_{};
};

//------------------------------------------------------------------------------
// Material probabilities mirrored onto a minimal openmc::Material handle
//------------------------------------------------------------------------------

class Material {
public:
  enum class Interaction { Scatter, Absorption, Fission };

  Material(std::string name, double total_xs, double scatter_ratio, double absorption_ratio,
    double fission_ratio);
  Material(Material&&) noexcept = default;
  Material& operator=(Material&&) noexcept = default;
  Material(const Material&) = delete;
  Material& operator=(const Material&) = delete;

  const std::string& name() const { return name_; }
  double total_xs() const { return total_xs_; }

  Interaction sample_interaction(double xi) const;
  double mean_free_path(double speed) const;

  openmc::Material& openmc_handle() { return *material_; }
  const openmc::Material& openmc_handle() const { return *material_; }

private:
  std::string name_;
  double total_xs_{};
  double scatter_weight_{};
  double absorption_weight_{};
  double fission_weight_{};
  std::unique_ptr<openmc::Material> material_;
};

class Cell {
public:
  Cell(std::string name, int material_id, std::vector<SurfaceSense> surfaces);

  const std::string& name() const { return name_; }
  int material_id() const { return material_id_; }
  const std::vector<SurfaceSense>& surfaces() const { return surfaces_; }

private:
  std::string name_;
  int material_id_{};
  std::vector<SurfaceSense> surfaces_;
};

class Universe {
public:
  Universe(std::string name, std::vector<int> cell_ids);

  const std::string& name() const { return name_; }
  const std::vector<int>& cell_ids() const { return cell_ids_; }

private:
  std::string name_;
  std::vector<int> cell_ids_;
};

class Lattice {
public:
  using Index3 = std::array<int, 3>;

  Lattice(std::string name, Index3 dims, openmc::Position pitch, std::vector<int> universe_ids);

  const std::string& name() const { return name_; }
  const Index3& dimensions() const { return dims_; }
  const openmc::Position& pitch() const { return pitch_; }
  int universe_at(int i, int j, int k) const;

private:
  std::string name_;
  Index3 dims_{};
  openmc::Position pitch_;
  std::vector<int> universe_ids_;
};

class GeometryRegistry {
public:
  GeometryRegistry();

  int add_material(Material material);
  int add_surface(std::shared_ptr<openmc::Surface> surface);
  int add_cell(Cell cell);
  int add_universe(Universe universe);
  int add_lattice(Lattice lattice);

  const Material& material(int id) const;
  const openmc::Surface& surface(int id) const;
  const Cell& cell(int id) const;
  const Universe& universe(int id) const;
  const Lattice& lattice(int id) const;

private:
  std::vector<Material> materials_;
  std::vector<std::shared_ptr<openmc::Surface>> surfaces_;
  std::vector<Cell> cells_;
  std::vector<Universe> universes_;
  std::vector<Lattice> lattices_;
};

struct CoordinateLevel {
  int universe_id{};
  int cell_id{};
  openmc::Position position{};
  openmc::Direction direction{};
};

class GeometryState {
public:
  explicit GeometryState(const GeometryRegistry* registry);

  void attach(openmc::GeometryState* state);
  void push(const CoordinateLevel& level);
  void pop();
  CoordinateLevel& current();
  const CoordinateLevel& current() const;
  const GeometryRegistry& registry() const;
  const Cell& current_cell() const;
  const Material& current_material() const;
  int depth() const;

  openmc::GeometryState& openmc_state();
  const openmc::GeometryState& openmc_state() const;

private:
  friend class Particle;
  void sync_to_openmc();

  const GeometryRegistry* registry_;
  openmc::GeometryState owned_state_{};
  openmc::GeometryState* state_{nullptr};
  std::vector<CoordinateLevel> stack_;
};

struct ParticleData {
  std::uint64_t id{};
  double weight{1.0};
  double energy{1.0};
  double age{0.0};
  std::mt19937_64 rng;
  openmc::ParticleData* openmc{nullptr};

  explicit ParticleData(std::uint64_t seed);
  double sample_unit();
};

struct TallyResult {
  double tracks{0.0};
  double collisions{0.0};
  double absorption{0.0};
  double fission{0.0};
};

class Particle {
public:
  Particle(GeometryState geometry, ParticleData data);
  Particle(Particle&&) noexcept = default;
  Particle& operator=(Particle&&) noexcept = default;
  Particle(const Particle&) = delete;
  Particle& operator=(const Particle&) = delete;

  GeometryState& geometry_state() { return geometry_; }
  const GeometryState& geometry_state() const { return geometry_; }

  ParticleData& data() { return data_; }
  const ParticleData& data() const { return data_; }

  openmc::Direction direction() const { return direction_; }
  void set_direction(openmc::Direction dir);

  bool alive() const { return alive_; }
  void kill() { alive_ = false; }

  double distance_to_boundary() const;
  void move(double distance);
  Material::Interaction collide();

private:
  void sync_particle_geometry();

  GeometryState geometry_;
  ParticleData data_;
  openmc::Particle particle_;
  openmc::Direction direction_{0.0, 0.0, 1.0};
  bool alive_{true};
};

GeometryRegistry build_demo_registry();
Particle make_demo_particle(const GeometryRegistry& registry, std::uint64_t seed = 1u);
void print_state_summary(const Particle& particle);
void run_chapter_simulation(const std::string& title, int histories,
  std::function<void(Particle&, TallyResult&)> stepper);

} // namespace tutorial

