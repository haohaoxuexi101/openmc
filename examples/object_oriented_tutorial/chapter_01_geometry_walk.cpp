#include "chapter_common.hpp"

#include <iostream>

using namespace tutorial;

int main() {
  GeometryRegistry registry = build_demo_registry();
  Particle particle = make_demo_particle(registry, 1u);
  print_state_summary(particle);

  std::cout << "Geometry stack depth: " << particle.geometry_state().depth() << "\n";
  std::cout << "Distance to next boundary: " << particle.distance_to_boundary() << "\n";
  return 0;
}
