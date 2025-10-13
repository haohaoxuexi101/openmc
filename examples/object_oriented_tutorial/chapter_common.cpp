#include "chapter_common.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

namespace tutorial {

namespace {
constexpr double kEps = 1e-12;

openmc::Direction normalize(openmc::Direction v)
{
  double n = v.norm();
  if (n <= 0.0) {
    return {0.0, 0.0, 1.0};
  }
  return v / n;
}

openmc::Position add_scaled(const openmc::Position& r, const openmc::Direction& u, double s)
{
  return {r.x + s * u.x, r.y + s * u.y, r.z + s * u.z};
}

} // namespace

//------------------------------------------------------------------------------
// PlaneSurface
//------------------------------------------------------------------------------

PlaneSurface::PlaneSurface(std::string name, openmc::Position normal, double d)
  : normal_(normalize(normal)), d_(d)
{
  name_ = std::move(name);
}

double PlaneSurface::evaluate(openmc::Position r) const
{
  return normal_.dot(r) - d_;
}

double PlaneSurface::distance(
  openmc::Position r, openmc::Direction u, bool /*coincident*/) const
{
  double denom = normal_.dot(u);
  if (std::fabs(denom) < kEps) {
    return std::numeric_limits<double>::infinity();
  }
  double numer = d_ - normal_.dot(r);
  double dist = numer / denom;
  if (dist <= kEps) {
    return std::numeric_limits<double>::infinity();
  }
  return dist;
}

openmc::Direction PlaneSurface::normal(openmc::Position /*r*/) const
{
  return normal_;
}

//------------------------------------------------------------------------------
// SphereSurface
//------------------------------------------------------------------------------

SphereSurface::SphereSurface(std::string name, openmc::Position center, double radius)
  : center_(center), radius_(radius)
{
  name_ = std::move(name);
}

double SphereSurface::evaluate(openmc::Position r) const
{
  auto diff = r - center_;
  return diff.dot(diff) - radius_ * radius_;
}

double SphereSurface::distance(
  openmc::Position r, openmc::Direction u, bool /*coincident*/) const
{
  auto offset = r - center_;
  double b = 2.0 * offset.dot(u);
  double c = offset.dot(offset) - radius_ * radius_;
  double disc = b * b - 4.0 * c;
  if (disc < 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  double sqrt_disc = std::sqrt(disc);
  double t1 = (-b - sqrt_disc) / 2.0;
  double t2 = (-b + sqrt_disc) / 2.0;
  double candidate = std::numeric_limits<double>::infinity();
  if (t1 > kEps) candidate = t1;
  if (t2 > kEps) candidate = std::min(candidate, t2);
  return candidate;
}

openmc::Direction SphereSurface::normal(openmc::Position r) const
{
  return normalize(r - center_);
}

//------------------------------------------------------------------------------
// Material
//------------------------------------------------------------------------------

Material::Material(std::string name, double total_xs, double scatter_ratio,
  double absorption_ratio, double fission_ratio)
  : name_(std::move(name)), total_xs_(total_xs),
    material_(std::make_unique<openmc::Material>())
{
  double sum = scatter_ratio + absorption_ratio + fission_ratio;
  if (sum <= 0.0) {
    throw std::invalid_argument("Material interaction ratios must be positive");
  }
  scatter_weight_ = scatter_ratio / sum;
  absorption_weight_ = absorption_ratio / sum;
  fission_weight_ = fission_ratio / sum;
  material_->set_name(name_);
  material_->set_density(total_xs_, "atom/b-cm");
}

Material::Interaction Material::sample_interaction(double xi) const
{
  if (xi < scatter_weight_) return Interaction::Scatter;
  xi -= scatter_weight_;
  if (xi < absorption_weight_) return Interaction::Absorption;
  return Interaction::Fission;
}

double Material::mean_free_path(double speed) const
{
  if (total_xs_ <= 0.0 || speed <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  return 1.0 / (total_xs_ * speed);
}

//------------------------------------------------------------------------------
// Cell / Universe / Lattice
//------------------------------------------------------------------------------

Cell::Cell(std::string name, int material_id, std::vector<SurfaceSense> surfaces)
  : name_(std::move(name)), material_id_(material_id), surfaces_(std::move(surfaces))
{}

Universe::Universe(std::string name, std::vector<int> cell_ids)
  : name_(std::move(name)), cell_ids_(std::move(cell_ids))
{}

Lattice::Lattice(std::string name, Index3 dims, openmc::Position pitch,
  std::vector<int> universe_ids)
  : name_(std::move(name)), dims_(dims), pitch_(pitch),
    universe_ids_(std::move(universe_ids))
{
  if (static_cast<std::size_t>(dims_[0] * dims_[1] * dims_[2]) != universe_ids_.size()) {
    throw std::invalid_argument("Universe list does not match lattice dimensions");
  }
}

int Lattice::universe_at(int i, int j, int k) const
{
  if (i < 0 || i >= dims_[0] || j < 0 || j >= dims_[1] || k < 0 || k >= dims_[2]) {
    throw std::out_of_range("Lattice index out of range");
  }
  std::size_t idx = static_cast<std::size_t>(i + dims_[0] * (j + dims_[1] * k));
  return universe_ids_.at(idx);
}

//------------------------------------------------------------------------------
// GeometryRegistry
//------------------------------------------------------------------------------

GeometryRegistry::GeometryRegistry()
{
  openmc::model::n_coord_levels = 3;
  openmc::model::root_universe = 0;
}

int GeometryRegistry::add_material(Material material)
{
  int id = static_cast<int>(materials_.size());
  material.openmc_handle().set_id(id + 1);
  materials_.push_back(std::move(material));
  return id;
}

int GeometryRegistry::add_surface(std::shared_ptr<openmc::Surface> surface)
{
  surfaces_.push_back(std::move(surface));
  return static_cast<int>(surfaces_.size() - 1);
}

int GeometryRegistry::add_cell(Cell cell)
{
  cells_.push_back(std::move(cell));
  return static_cast<int>(cells_.size() - 1);
}

int GeometryRegistry::add_universe(Universe universe)
{
  universes_.push_back(std::move(universe));
  return static_cast<int>(universes_.size() - 1);
}

int GeometryRegistry::add_lattice(Lattice lattice)
{
  lattices_.push_back(std::move(lattice));
  return static_cast<int>(lattices_.size() - 1);
}

const Material& GeometryRegistry::material(int id) const
{
  return materials_.at(static_cast<std::size_t>(id));
}

const openmc::Surface& GeometryRegistry::surface(int id) const
{
  return *surfaces_.at(static_cast<std::size_t>(id));
}

const Cell& GeometryRegistry::cell(int id) const
{
  return cells_.at(static_cast<std::size_t>(id));
}

const Universe& GeometryRegistry::universe(int id) const
{
  return universes_.at(static_cast<std::size_t>(id));
}

const Lattice& GeometryRegistry::lattice(int id) const
{
  return lattices_.at(static_cast<std::size_t>(id));
}

//------------------------------------------------------------------------------
// GeometryState
//------------------------------------------------------------------------------

GeometryState::GeometryState(const GeometryRegistry* registry)
  : registry_(registry)
{
  if (!registry_) {
    throw std::invalid_argument("GeometryState requires a valid registry");
  }
  state_ = &owned_state_;
  sync_to_openmc();
}

void GeometryState::attach(openmc::GeometryState* state)
{
  state_ = state ? state : &owned_state_;
  sync_to_openmc();
}

void GeometryState::push(const CoordinateLevel& level)
{
  if (stack_.size() >= static_cast<std::size_t>(openmc::model::n_coord_levels)) {
    throw std::runtime_error("Too many coordinate levels for demo model");
  }
  stack_.push_back(level);
  sync_to_openmc();
}

void GeometryState::pop()
{
  if (stack_.empty()) {
    throw std::runtime_error("Cannot pop from empty geometry stack");
  }
  stack_.pop_back();
  sync_to_openmc();
}

CoordinateLevel& GeometryState::current()
{
  if (stack_.empty()) {
    throw std::runtime_error("Geometry stack is empty");
  }
  return stack_.back();
}

const CoordinateLevel& GeometryState::current() const
{
  if (stack_.empty()) {
    throw std::runtime_error("Geometry stack is empty");
  }
  return stack_.back();
}

const GeometryRegistry& GeometryState::registry() const
{
  return *registry_;
}

const Cell& GeometryState::current_cell() const
{
  return registry_->cell(current().cell_id);
}

const Material& GeometryState::current_material() const
{
  return registry_->material(current_cell().material_id());
}

int GeometryState::depth() const
{
  return static_cast<int>(stack_.size());
}

openmc::GeometryState& GeometryState::openmc_state()
{
  return *(state_ ? state_ : &owned_state_);
}

const openmc::GeometryState& GeometryState::openmc_state() const
{
  return *(state_ ? state_ : &owned_state_);
}

void GeometryState::sync_to_openmc()
{
  openmc::GeometryState& state_ref = *(state_ ? state_ : &owned_state_);
  state_ref.clear();
  for (std::size_t i = 0; i < stack_.size(); ++i) {
    auto& coord = state_ref.coord(static_cast<int>(i));
    coord.r = stack_[i].position;
    coord.u = normalize(stack_[i].direction);
    coord.cell = stack_[i].cell_id;
    coord.universe = stack_[i].universe_id;
  }
  state_ref.n_coord() = static_cast<int>(stack_.size() > 0 ? stack_.size() : 1);
  if (!stack_.empty()) {
    state_ref.surface() = openmc::SURFACE_NONE;
    state_ref.material() = current_cell().material_id();
    state_ref.u() = normalize(stack_.front().direction);
    state_ref.r() = stack_.front().position;
  }
}

//------------------------------------------------------------------------------
// ParticleData
//------------------------------------------------------------------------------

ParticleData::ParticleData(std::uint64_t seed) : id(seed), rng(seed) {}

double ParticleData::sample_unit()
{
  std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
  return unit_dist(rng);
}

//------------------------------------------------------------------------------
// Particle
//------------------------------------------------------------------------------

Particle::Particle(GeometryState geometry, ParticleData data)
  : geometry_(std::move(geometry)), data_(std::move(data))
{
  geometry_.attach(&particle_);
  geometry_.sync_to_openmc();
  data_.openmc = &particle_;

  particle_.type() = openmc::ParticleType::neutron;
  particle_.E() = data_.energy;
  particle_.wgt() = data_.weight;
  particle_.time() = data_.age;
  if (geometry_.depth() > 0) {
    direction_ = normalize(geometry_.current().direction);
    geometry_.current().direction = direction_;
    particle_.u() = direction_;
  }
}

void Particle::sync_particle_geometry()
{
  geometry_.sync_to_openmc();
  particle_.u() = direction_;
  particle_.r() = geometry_.openmc_state().r();
  particle_.material() = geometry_.openmc_state().material();
}

void Particle::set_direction(openmc::Direction dir)
{
  direction_ = normalize(dir);
  geometry_.current().direction = direction_;
  sync_particle_geometry();
}

double Particle::distance_to_boundary() const
{
  const auto& cell = geometry_.current_cell();
  const auto& registry = geometry_.registry();
  double min_distance = std::numeric_limits<double>::infinity();
  for (const auto& sense : cell.surfaces()) {
    const auto& surface = registry.surface(sense.surface_id);
    double distance = surface.distance(geometry_.current().position, direction_, false);
    if (distance < min_distance) {
      min_distance = distance;
    }
  }
  return min_distance;
}

void Particle::move(double distance)
{
  geometry_.current().position = add_scaled(geometry_.current().position, direction_, distance);
  data_.age += distance;
  particle_.time() = data_.age;
  sync_particle_geometry();
}

Material::Interaction Particle::collide()
{
  auto interaction = geometry_.current_material().sample_interaction(data_.sample_unit());
  switch (interaction) {
  case Material::Interaction::Scatter: {
    double mu = 2.0 * data_.sample_unit() - 1.0;
    constexpr double pi = 3.14159265358979323846;
    double phi = 2.0 * pi * data_.sample_unit();
    double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
    set_direction({sin_theta * std::cos(phi), sin_theta * std::sin(phi), mu});
    break;
  }
  case Material::Interaction::Absorption:
    alive_ = false;
    data_.weight = 0.0;
    particle_.wgt() = 0.0;
    break;
  case Material::Interaction::Fission:
    data_.weight *= 2.0;
    particle_.wgt() = data_.weight;
    break;
  }
  return interaction;
}

//------------------------------------------------------------------------------
// Demo helpers
//------------------------------------------------------------------------------

static GeometryRegistry build_registry_impl()
{
  GeometryRegistry registry;

  int fuel = registry.add_material(Material{"Fuel", 0.8, 0.6, 0.3, 0.1});
  int moderator = registry.add_material(Material{"Moderator", 0.2, 0.9, 0.1, 0.0});

  auto left = registry.add_surface(
    std::make_shared<PlaneSurface>("Left", openmc::Position{-1.0, 0.0, 0.0}, -0.5));
  auto right = registry.add_surface(
    std::make_shared<PlaneSurface>("Right", openmc::Position{1.0, 0.0, 0.0}, 0.5));
  auto bottom = registry.add_surface(
    std::make_shared<PlaneSurface>("Bottom", openmc::Position{0.0, -1.0, 0.0}, -0.5));
  auto top = registry.add_surface(
    std::make_shared<PlaneSurface>("Top", openmc::Position{0.0, 1.0, 0.0}, 0.5));
  auto sphere = registry.add_surface(
    std::make_shared<SphereSurface>("FuelSphere", openmc::Position{0.0, 0.0, 0.0}, 0.3));

  int fuel_cell =
    registry.add_cell(Cell{"FuelCell", fuel, {{sphere, false}}});
  int moderator_cell = registry.add_cell(Cell{"ModeratorCell", moderator,
    {{left, true}, {right, true}, {bottom, true}, {top, true}}});

  int core_universe = registry.add_universe(Universe{"Core", {fuel_cell, moderator_cell}});

  Lattice::Index3 dims{2, 2, 1};
  std::vector<int> lattice_universes = {core_universe, core_universe, core_universe, core_universe};
  registry.add_lattice(
    Lattice{"QuarterCore", dims, openmc::Position{0.5, 0.5, 1.0}, lattice_universes});

  return registry;
}

GeometryRegistry build_demo_registry()
{
  return build_registry_impl();
}

Particle make_demo_particle(const GeometryRegistry& registry, std::uint64_t seed)
{
  GeometryState geometry{&registry};
  CoordinateLevel level;
  level.universe_id = 0;
  level.cell_id = 0;
  level.position = {0.0, 0.0, 0.0};
  level.direction = {0.0, 0.0, 1.0};
  geometry.push(level);

  ParticleData data{seed};
  data.energy = 1.0;
  data.weight = 1.0;
  return Particle{std::move(geometry), std::move(data)};
}

void print_state_summary(const Particle& particle)
{
  const auto& state = particle.geometry_state().current();
  const auto& cell = particle.geometry_state().current_cell();
  const auto& material = particle.geometry_state().current_material();

  std::cout << "Particle at (" << state.position.x << ", " << state.position.y << ", "
            << state.position.z << ")\n";
  std::cout << "  Cell: " << cell.name() << " (material " << material.name() << ")\n";
  auto dir = particle.direction();
  std::cout << "  Direction: (" << dir.x << ", " << dir.y << ", " << dir.z << ")\n";
}

void run_chapter_simulation(const std::string& title, int histories,
  std::function<void(Particle&, TallyResult&)> stepper)
{
  std::cout << "=== " << title << " ===\n";
  GeometryRegistry registry = build_demo_registry();
  TallyResult totals;

  for (int n = 0; n < histories; ++n) {
    Particle particle = make_demo_particle(registry, static_cast<std::uint64_t>(n + 1));
    TallyResult local;
    while (particle.alive()) {
      stepper(particle, local);
    }
    totals.tracks += local.tracks;
    totals.collisions += local.collisions;
    totals.absorption += local.absorption;
    totals.fission += local.fission;
  }

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "Tracks: " << totals.tracks / histories
            << "  Collisions: " << totals.collisions / histories
            << "  Absorption: " << totals.absorption / histories
            << "  Fission weight: " << totals.fission / histories << "\n\n";
}

} // namespace tutorial

