#include "chapter_common.hpp"

#include <iostream>

using namespace tutorial;

int main() {
  GeometryRegistry registry = build_demo_registry();
  Particle particle = make_demo_particle(registry, 7u);
  particle.set_direction({0.2, 0.0, 1.0});

  double distance = particle.distance_to_boundary();
  std::cout << "Initial distance to boundary: " << distance << "\n";
  particle.move(distance);
  print_state_summary(particle);
  return 0;
}
