#include "chapter_common.hpp"

#include <algorithm>
#include <cmath>

using namespace tutorial;

int main() {
  run_chapter_simulation("Collision statistics", 50, [](Particle& particle, TallyResult& tally) {
    double boundary = particle.distance_to_boundary();
    double collision = -std::log(particle.data().sample_unit()) / particle.geometry_state().current_material().total_xs();
    double step = std::min(boundary, collision);
    particle.move(step);
    tally.tracks += step;
    if (step >= collision - 1e-12) {
      auto interaction = particle.collide();
      ++tally.collisions;
      if (interaction == Material::Interaction::Absorption) {
        tally.absorption += particle.data().weight;
      } else if (interaction == Material::Interaction::Fission) {
        tally.fission += particle.data().weight;
      }
    }
    if (step >= boundary - 1e-12) {
      particle.kill();
    }
  });
  return 0;
}
