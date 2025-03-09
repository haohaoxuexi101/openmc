#include "openmc/simulation.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/container_util.h"
#include "openmc/cross_sections.h"
// #include "openmc/distribution.h"
// #include "openmc/distribution_spatial.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/initialize.h"
#include "openmc/material.h"
#include "openmc/mcpl_interface.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
#include "openmc/particle_data.h"
#include "openmc/photon.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/state_point.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/trigger.h"
#include "openmc/timer.h"
#include "openmc/track_output.h"
#include "openmc/weight_windows.h"
// #include "openmc/xml_interface.h"

#include "lattice.h"
#include "openmc/physics.h"
#include "openmc/physics_mg.h"
#include "openmc/tallies/tally_scoring.h"
// particle.h中没有#include "openmc/mgxs_interface.h"
// 因此这里必须包含
#include "openmc/cell.h"
#include "openmc/mgxs_interface.h"

namespace openmc {
class DeltaParticle : public openmc::Particle {
public:
  DeltaParticle() = default;
  void delta_calculate_xs();
  void delta_advance();
  void delta_cross_surface();
  void delta_collide();
  void get_all_materials();
  std::unordered_map<int, double> macro_xs_t;
  // Macroscopic cross sections
  MacroXS& max_macro_xs() { return max_macro_xs_; }
  const MacroXS& max_macro_xs() const { return max_macro_xs_; }
  bool vir_collision_ {true};

  // private:
  MacroXS max_macro_xs_;
};

bool find_cell_inner(
  GeometryState& p, const NeighborList* neighbor_list, bool verbose)
{
  // Find which cell of this universe the particle is in.  Use the neighbor list
  // to shorten the search if one was provided.
  bool found = false;
  int32_t i_cell = C_NONE;
  if (neighbor_list) {
    for (auto it = neighbor_list->cbegin(); it != neighbor_list->cend(); ++it) {
      i_cell = *it;

      // Make sure the search cell is in the same universe.
      int i_universe = p.lowest_coord().universe;
      if (model::cells[i_cell]->universe_ != i_universe)
        continue;

      // Check if this cell contains the particle.
      Position r {p.r_local()};
      Direction u {p.u_local()};
      auto surf = p.surface();
      if (model::cells[i_cell]->contains(r, u, surf)) {
        p.lowest_coord().cell = i_cell;
        found = true;
        break;
      }
    }

    // If we're attempting a neighbor list search and fail, we
    // now know we should return false. This will trigger an
    // exhaustive search from neighbor_list_find_cell and make
    // the result from that be appended to the neighbor list.
    if (!found) {
      return found;
    }
  }

  // Check successively lower coordinate levels until finding material fill
  for (;; ++p.n_coord()) {
    // If we did not attempt to use neighbor lists, i_cell is still C_NONE.  In
    // that case, we should now do an exhaustive search to find the right value
    // of i_cell.
    //
    // Alternatively, neighbor list searches could have succeeded, but we found
    // that the fill of the neighbor cell was another universe. As such, in the
    // code below this conditional, we set i_cell back to C_NONE to indicate
    // that.
    if (i_cell == C_NONE) {
      int i_universe = p.lowest_coord().universe;
      const auto& univ {model::universes[i_universe]};
      found = univ->find_cell(p);
    }

    if (!found) {
      return found;
    }
    i_cell = p.lowest_coord().cell;

    // Announce the cell that the particle is entering.
    if (found && verbose) {
      auto msg = fmt::format("    Entering cell {}", model::cells[i_cell]->id_);
      write_message(msg, 1);
    }

    Cell& c {*model::cells[i_cell]};
    if (c.type_ == Fill::MATERIAL) {
      // Found a material cell which means this is the lowest coord level.

      p.cell_instance() = 0;
      // Find the distribcell instance number.
      if (c.distribcell_index_ >= 0) {
        p.cell_instance() = cell_instance_at_level(p, p.n_coord() - 1);
      }

      // Set the material and temperature.
      p.material_last() = p.material();
      p.material() = c.material(p.cell_instance());
      p.sqrtkT_last() = p.sqrtkT();
      p.sqrtkT() = c.sqrtkT(p.cell_instance());

      return true;

    } else if (c.type_ == Fill::UNIVERSE) {
      //========================================================================
      //! Found a lower universe, update this coord level then search the next.

      // Set the lower coordinate level universe.
      auto& coord {p.coord(p.n_coord())};
      coord.universe = c.fill_;

      // Set the position and direction.
      coord.r = p.r_local();
      coord.u = p.u_local();

      // Apply translation.
      coord.r -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

    } else if (c.type_ == Fill::LATTICE) {
      //========================================================================
      //! Found a lower lattice, update this coord level then search the next.

      Lattice& lat {*model::lattices[c.fill_]};

      // Set the position and direction.
      auto& coord {p.coord(p.n_coord())};
      coord.r = p.r_local();
      coord.u = p.u_local();

      // Apply translation.
      coord.r -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

      // Determine lattice indices.
      auto& i_xyz {coord.lattice_i};
      lat.get_indices(coord.r, coord.u, i_xyz);

      // Get local position in appropriate lattice cell
      coord.r = lat.get_local_position(coord.r, i_xyz);

      // Set lattice indices.
      coord.lattice = c.fill_;

      // Set the lower coordinate level universe.
      if (lat.are_valid_indices(i_xyz)) {
        coord.universe = lat[i_xyz];
      } else {
        if (lat.outer_ != NO_OUTER_UNIVERSE) {
          coord.universe = lat.outer_;
        } else {
          p.mark_as_lost(fmt::format(
            "Particle {} left lattice {}, but it has no outer definition.",
            p.id(), lat.id_));
        }
      }
    }
    i_cell = C_NONE; // trip non-neighbor cell search at next iteration
    found = false;
  }

  return found;
}

bool delta_exhaustive_find_cell(GeometryState& p, bool verbose)
{
  // auto& coord {p.lowest_coord()};
  // auto& lat {*model::lattices[coord.lattice]};
  // // Set the new coordinate position.
  // const auto& upper_coord {p.coord(p.n_coord() - 2)};
  // const auto& cell {model::cells[upper_coord.cell]};
  // Position r = upper_coord.r;
  // int i_universe = p.lowest_coord().universe;
  // if (i_universe == C_NONE) {
  p.coord(0).universe = model::root_universe;
  p.n_coord() = 1;
  // p.r_local() = lat.get_local_position(r, coord.lattice_i);
  int i_universe = model::root_universe;
  // }
  // Reset all the deeper coordinate levels.
  for (int i = p.n_coord(); i < model::n_coord_levels; i++) {
    p.coord(i).reset();
  }
  return find_cell_inner(p, nullptr, verbose);
}

void delta_cross_lattice(
  GeometryState& p, const BoundaryInfo& boundary, bool verbose)
{
  auto& coord {p.lowest_coord()};
  auto& lat {*model::lattices[coord.lattice]};

  if (verbose) {
    write_message(
      fmt::format("    Crossing lattice {}. Current position ({},{},{}). r={}",
        lat.id_, coord.lattice_i[0], coord.lattice_i[1], coord.lattice_i[2],
        p.r()),
      1);
  }

  // Set the lattice indices.
  coord.lattice_i[0] += boundary.lattice_translation[0];
  coord.lattice_i[1] += boundary.lattice_translation[1];
  coord.lattice_i[2] += boundary.lattice_translation[2];

  // Set the new coordinate position.
  const auto& upper_coord {p.coord(p.n_coord() - 2)};
  const auto& cell {model::cells[upper_coord.cell]};
  Position r = upper_coord.r;
  r -= cell->translation_;
  if (!cell->rotation_.empty()) {
    r = r.rotate(cell->rotation_);
  }
  p.r_local() = lat.get_local_position(r, coord.lattice_i);

  if (!lat.are_valid_indices(coord.lattice_i)) {
    // The particle is outside the lattice.  Search for it from the base coords.
    p.n_coord() = 1;
    bool found = exhaustive_find_cell(p);

    if (!found) {
      p.mark_as_lost(fmt::format("Particle {} could not be located after "
                                 "crossing a boundary of lattice {}",
        p.id(), lat.id_));
    }

  } else {
    // Find cell in next lattice element.
    p.lowest_coord().universe = lat[coord.lattice_i];
    bool found = exhaustive_find_cell(p);

    if (!found) {
      // A particle crossing the corner of a lattice tile may not be found.  In
      // this case, search for it from the base coords.
      p.n_coord() = 1;
      bool found = exhaustive_find_cell(p);
      if (!found) {
        p.mark_as_lost(fmt::format("Particle {} could not be located after "
                                   "crossing a boundary of lattice {}",
          p.id(), lat.id_));
      }
    }
  }
}

void DeltaParticle::delta_cross_surface()
{
  // Saving previous cell data
  for (int j = 0; j < n_coord(); ++j) {
    cell_last(j) = coord(j).cell;
  }
  n_coord_last() = n_coord();

  // Set surface that particle is on and adjust coordinate levels
  surface() = boundary().surface;
  n_coord() = boundary().coord_level;

  if (boundary().lattice_translation[0] != 0 ||
      boundary().lattice_translation[1] != 0 ||
      boundary().lattice_translation[2] != 0) {
    // Particle crosses lattice boundary

    bool verbose = settings::verbosity >= 10 || trace();
    delta_cross_lattice(*this, boundary(), verbose);
    event() = TallyEvent::LATTICE;
  } else {
    // Particle crosses surface
    // TODO: off-by-one
    const auto& surf {model::surfaces[surface_index()].get()};
    // If BC, add particle to surface source before crossing surface
    if (surf->surf_source_ && surf->bc_) {
      add_surf_source_to_bank(*this, *surf);
    }
    this->cross_surface(*surf);
    // If no BC, add particle to surface source after crossing surface
    if (surf->surf_source_ && !surf->bc_) {
      add_surf_source_to_bank(*this, *surf);
    }
    if (settings::weight_window_checkpoint_surface) {
      apply_weight_windows(*this);
    }
    event() = TallyEvent::SURFACE;
  }
  // Score cell to cell partial currents
  if (!model::active_surface_tallies.empty()) {
    score_surface_tally(*this, model::active_surface_tallies);
  }
}

void DeltaParticle::delta_calculate_xs()
{
  // Set the random number stream
  stream() = STREAM_TRACKING;

  // Store pre-collision particle properties
  wgt_last() = wgt();
  E_last() = E();
  u_last() = u();
  r_last() = r();
  time_last() = time();

  // Reset event variables
  event() = TallyEvent::KILL;
  event_nuclide() = NUCLIDE_NONE;
  event_mt() = REACTION_NONE;

  // If the cell hasn't been determined based on the particle's location,
  // initiate a search for the current cell. This generally happens at the
  // beginning of the history and again for any secondary particles
  if (lowest_coord().cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
      return;
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = lowest_coord().cell;

    // Initialize last cells from current cell
    for (int j = 0; j < n_coord(); ++j) {
      cell_last(j) = coord(j).cell;
    }
    n_coord_last() = n_coord();
  }

  // Write particle track.
  if (write_track())
    write_particle_track(*this);

  if (settings::check_overlaps)
    check_cell_overlap(*this);

  for (auto& c : model::cells) {
    // Ignore non-material cells
    if (c->material_.size() == 0) {
      // std::cout << "the cell id without material: " << c->id_ << std::endl;
      continue;
    }
    // std::cout << "the cell id: " << c->id_ << std::endl;
    for (auto i_mat : c->material_) {
      if (i_mat != MATERIAL_VOID) {
        const auto& mat {model::materials[i_mat]};
        std::cout << "the material id: " << mat->id_ << std::endl;
        model::materials[i_mat]->calculate_xs(*this);
        std::cout << "attempt to calculate the xs" << std::endl;
        std::cout << "the total xs: " << this->macro_xs().total << std::endl;
        std::cout
          << "----------------------------------------------------------------"
          << std::endl;

        // using an unordered map to store the total macroscopic xs with the
        // cell id
        this->macro_xs_t[c->id_] = this->macro_xs().total;
      }
    }
  }

  if (lowest_coord().cell == C_NONE) {
    if (!exhaustive_find_cell(*this)) {
      mark_as_lost(
        "Could not find the cell containing particle " + std::to_string(id()));
      return;
    }

    // Set birth cell attribute
    if (cell_born() == C_NONE)
      cell_born() = lowest_coord().cell;

    // Initialize last cells from current cell
    for (int j = 0; j < n_coord(); ++j) {
      cell_last(j) = coord(j).cell;
    }
    n_coord_last() = n_coord();
  }

  if (material() != MATERIAL_VOID) {
    if (settings::run_CE) {
      if (material() != material_last() || sqrtkT() != sqrtkT_last()) {
        // If the material is the same as the last material and the
        // temperature hasn't changed, we don't need to lookup cross
        // sections again.
        model::materials[material()]->calculate_xs(*this);
      }
    }
  }

  // get the maximum total macroscopic xs from the unordered map
  double max_sigma_t = 0.0;
  int max_cell_id = 0;
  for (auto& i : this->macro_xs_t) {
    if (i.second > max_sigma_t) {
      max_sigma_t = i.second;
      max_cell_id = i.first;
    }
  }

  this->max_macro_xs().total = max_sigma_t;
  std::cout << "the max total macroscopic xs: " << max_sigma_t << std::endl;
  // std::cout << "the cell id with the max total macroscopic xs: " <<
  // max_cell_id
  //           << std::endl;
}

void DeltaParticle::delta_advance()
{

  // Find the distance to the nearest boundary
  // boundary() = distance_to_boundary(*this);

  // Sample a distance to collision

  max_macro_xs().total = 10;

  if (type() == ParticleType::electron || type() == ParticleType::positron) {
    collision_distance() = 0.0;
  } else if (macro_xs().total == 0.0) {
    collision_distance() = INFINITY;
  } else {
    collision_distance() =
      -std::log(prn(current_seed())) / max_macro_xs().total;
  }

  // Select smaller of the two distances
  // double distance = std::min(boundary().distance, collision_distance());
  double distance = collision_distance();

  // Advance particle in space and time
  // Short-term solution until the surface source is revised and we can use
  // this->move_distance(distance)

  // Saving previous cell data
  for (int j = 0; j < n_coord(); ++j) {
    cell_last(j) = coord(j).cell;
  }
  n_coord_last() = n_coord();

  for (int j = 0; j < n_coord(); ++j) {
    coord(j).r += distance * coord(j).u;
  }

  if (delta_exhaustive_find_cell(*this, false)) {
    int cell_index = coord(n_coord() - 1).cell;
    std::cout << "find the cell: " << coord(n_coord() - 1).cell << std::endl;
  }

  if (this->macro_xs().total >
      prn(current_seed()) * this->max_macro_xs().total) {
    this->vir_collision_ = false;
  } else {
    std::cout << "the virtual collision distance: " << distance << std::endl;
  }
}

void DeltaParticle::delta_collide()
{
  // Score collision estimate of keff
  if (settings::run_mode == RunMode::EIGENVALUE &&
      type() == ParticleType::neutron) {
    keff_tally_collision() += wgt() * macro_xs().nu_fission / macro_xs().total;
  }

  // Score surface current tallies -- this has to be done before the
  // collision since the direction of the particle will change and we need
  // to use the pre-collision direction to figure out what mesh surfaces
  // were crossed

  if (!model::active_meshsurf_tallies.empty())
    score_surface_tally(*this, model::active_meshsurf_tallies);

  // Clear surface component
  surface() = SURFACE_NONE;

  if (settings::run_CE) {
    collision(*this);
  } else {
    collision_mg(*this);
  }

  // Score collision estimator tallies -- this is done after a collision
  // has occurred rather than before because we need information on the
  // outgoing energy for any tallies with an outgoing energy filter
  if (!model::active_collision_tallies.empty())
    score_collision_tally(*this);
  if (!model::active_analog_tallies.empty()) {
    if (settings::run_CE) {
      score_analog_tally_ce(*this);
    } else {
      score_analog_tally_mg(*this);
    }
  }

  if (!model::active_pulse_height_tallies.empty() &&
      type() == ParticleType::photon) {
    pht_collision_energy();
  }

  // Reset banked weight during collision
  n_bank() = 0;
  bank_second_E() = 0.0;
  wgt_bank() = 0.0;
  zero_delayed_bank();

  // Reset fission logical
  fission() = false;

  // Save coordinates for tallying purposes
  r_last_current() = r();

  // Set last material to none since cross sections will need to be
  // re-evaluated
  material_last() = C_NONE;

  // Set all directions to base level -- right now, after a collision, only
  // the base level directions are changed
  for (int j = 0; j < n_coord() - 1; ++j) {
    if (coord(j + 1).rotated) {
      // If next level is rotated, apply rotation matrix
      const auto& m {model::cells[coord(j).cell]->rotation_};
      const auto& u {coord(j).u};
      coord(j + 1).u = u.rotate(m);
    } else {
      // Otherwise, copy this level's direction
      coord(j + 1).u = coord(j).u;
    }
  }

  // Score flux derivative accumulators for differential tallies.
  if (!model::active_tallies.empty())
    score_collision_derivative(*this);

#ifdef DAGMC
  history().reset();
#endif
}

void DeltaParticle::get_all_materials()
{
  for (auto& c : model::cells) {
    // Ignore non-material cells
    if (c->material_.size() == 0) {
      std::cout << "the cell id without material: " << c->id_ << std::endl;
      continue;
    }
    std::cout << "the cell id: " << c->id_ << std::endl;
    for (auto i_mat : c->material_) {
      if (i_mat != MATERIAL_VOID) {
        const auto& mat {model::materials[i_mat]};
        std::cout << "the material id: " << mat->id_ << std::endl;
        // model::materials[i_mat]->calculate_xs(*this);
        std::cout << "attempt to calculate the xs" << std::endl;
      }
    }
  }
}

} // namespace openmc

int main()
{

  openmc::settings::path_input = std::string(
    "/home/ssn/ssn_mc/openmc/OpenMC_CPP_TESTS/openmc_develop_test/");

  // openmc::read_model_xml();

  openmc::read_separate_xml_files();

  openmc_simulation_init();

  // // 局部类必须内联定义！！！
  // class DeltaParticle : public openmc::Particle {
  // public:
  //   void delta_calculate_xs()
  //   {
  //     std::cout << "DeltaParticle::delta_calculate_xs()" << std::endl;
  //   }
  // };

  openmc::DeltaParticle p;

  openmc::initialize_history(p, 1);

  // p.get_all_materials();
  while (p.vir_collision_ == true) {
    // p.delta_calculate_xs();
    p.event_calculate_xs();
    // p.event_advance();
    p.delta_advance();
    // p.delta_cross_surface();
    // p.event_cross_surface();
    // event_cross_surface()之后材料截面会变化
  }
  //   std::cout << p.r().x << std::endl;

  return 0;
}