// model.h
#ifndef MODEL_H
#define MODEL_H

#include <memory>
#include <unordered_map>
#include <vector>

class Cell {
public:
  int id;
};

namespace model {
extern std::vector<std::unique_ptr<Cell>> cells;
extern std::unordered_map<int32_t, int32_t> cell_map;
void initialize_cells();
} // namespace model
// void initialize_cells();
#endif // MODEL_H
