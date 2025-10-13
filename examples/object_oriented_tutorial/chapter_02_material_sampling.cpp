#include "chapter_common.hpp"

#include <iostream>
#include <map>

using namespace tutorial;

int main() {
  GeometryRegistry registry = build_demo_registry();
  Particle particle = make_demo_particle(registry, 42u);

  std::map<Material::Interaction, int> counts;
  for (int n = 0; n < 1000; ++n) {
    auto interaction = particle.collide();
    counts[interaction]++;
    if (!particle.alive()) {
      particle = make_demo_particle(registry, static_cast<std::uint64_t>(n + 43));
    }
  }

  std::cout << "Sampling distribution after 1000 draws:\n";
  std::cout << "  Scatter: " << counts[Material::Interaction::Scatter] << "\n";
  std::cout << "  Absorption: " << counts[Material::Interaction::Absorption] << "\n";
  std::cout << "  Fission: " << counts[Material::Interaction::Fission] << "\n";
  return 0;
}
