#include "chapter_common.hpp"

#include <iostream>

using namespace tutorial;

int main() {
  GeometryRegistry registry = build_demo_registry();
  const auto& lattice = registry.lattice(0);
  std::cout << "Lattice '" << lattice.name() << "' dimensions: "
            << lattice.dimensions()[0] << "x" << lattice.dimensions()[1] << "x" << lattice.dimensions()[2] << "\n";
  for (int j = 0; j < lattice.dimensions()[1]; ++j) {
    for (int i = 0; i < lattice.dimensions()[0]; ++i) {
      std::cout << lattice.universe_at(i, j, 0) << ' ';
    }
    std::cout << '\n';
  }
  return 0;
}
