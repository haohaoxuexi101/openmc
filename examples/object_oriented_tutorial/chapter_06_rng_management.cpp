#include "chapter_common.hpp"

#include <iomanip>
#include <iostream>

using namespace tutorial;

int main() {
  GeometryRegistry registry = build_demo_registry();
  Particle a = make_demo_particle(registry, 1234u);
  Particle b = make_demo_particle(registry, 1234u);

  std::cout << std::fixed << std::setprecision(6);
  for (int i = 0; i < 5; ++i) {
    std::cout << "Draw " << i << ": " << a.data().sample_unit() << " vs " << b.data().sample_unit() << '\n';
  }
  return 0;
}
