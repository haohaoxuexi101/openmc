// model.cpp
#include "model.h"

namespace model {
std::vector<std::unique_ptr<Cell>> cells;
std::unordered_map<int32_t, int32_t> cell_map;
void initialize_cells()
{
  cells.push_back(std::make_unique<Cell>());
  cells[0]->id = 1;
  cell_map[1] = 100;
}
} // namespace model
// void initialize_cells()
// {
//   model::cells.push_back(std::make_unique<Cell>());
//   model::cells[0]->id = 1;
//   model::cell_map[1] = 100;
// }