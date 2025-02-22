// geometry.cpp
#include "model.h"
#include <iostream>

void print_cells_info()
{
  for (const auto& pair : model::cell_map) {
    std::cout << "Cell ID: " << pair.first << " has data: " << pair.second
              << std::endl;
  }

  for (const auto& cell : model::cells) {
    if (cell) {
      std::cout << "Cell ID: " << cell->id << std::endl;
    }
  }
}

int main()
{
  model::initialize_cells();
  //   initialize_cells();
  print_cells_info();
  return 0;
}
