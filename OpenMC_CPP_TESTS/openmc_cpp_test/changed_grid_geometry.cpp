#include "openmc/geometry.h"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <unordered_map>
#include <vector>
#include <cmath>

#include "openmc/array.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/lattice.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/string_utils.h"
#include "openmc/surface.h"

namespace openmc
{

    //==============================================================================
    // Global variables
    //==============================================================================

    namespace model
    {

        int root_universe{-1};
        int n_coord_levels;

        vector<int64_t> overlap_check_count;

    } // namespace model

    //==============================================================================
    // Grid Structure Implementation for Cell Lookup
    //==============================================================================

    struct GridCell
    {
        std::vector<int> cell_indices; // Stores cell indices in this grid cell
    };

    std::unordered_map<int, GridCell> grid;
    int grid_size = 10;        // Example number of grid cells per axis
    double grid_spacing = 1.0; // Example spacing between grid cells

    // Function to get grid index based on position
    int get_grid_index(const Position &pos)
    {
        int x_idx = static_cast<int>(pos.x / grid_spacing);
        int y_idx = static_cast<int>(pos.y / grid_spacing);
        int z_idx = static_cast<int>(pos.z / grid_spacing);
        return x_idx + y_idx * grid_size + z_idx * grid_size * grid_size;
    }

    // Populate grid with cell indices and positions
    void populate_grid(const std::vector<int> &cell_indices, const std::vector<Position> &positions)
    {
        for (size_t i = 0; i < cell_indices.size(); ++i)
        {
            int grid_index = get_grid_index(positions[i]);
            grid[grid_index].cell_indices.push_back(cell_indices[i]);
        }
    }

    // Search grid for a cell at a given position
    bool search_grid(const Position &pos, int &found_index)
    {
        int grid_index = get_grid_index(pos);
        if (grid.find(grid_index) != grid.end())
        {
            found_index = grid[grid_index].cell_indices[0]; // Return first matching cell for simplicity
            return true;
        }
        return false;
    }

    //==============================================================================
    // Non-member functions
    //==============================================================================

    bool check_cell_overlap(GeometryState &p, bool error)
    {
        int n_coord = p.n_coord();
        bool overlap_detected = false;

        std::vector<int> cell_indices;
        std::vector<Position> cell_positions;
        for (int j = 0; j < n_coord; j++)
        {
            Universe &univ = *model::universes[p.coord(j).universe];
            for (auto index_cell : univ.cells_)
            {
                cell_indices.push_back(index_cell);
                cell_positions.push_back(p.coord(j).r);
            }
        }

        populate_grid(cell_indices, cell_positions);

        for (int j = 0; j < n_coord; j++)
        {
            int found_index;
            if (search_grid(p.coord(j).r, found_index))
            {
                Cell &c = *model::cells[found_index];
                if (c.contains(p.coord(j).r, p.coord(j).u, p.surface()))
                {
                    if (found_index != p.coord(j).cell)
                    {
                        if (error)
                        {
                            fatal_error(fmt::format("Overlapping cells detected: {}, {} on universe {}",
                                                    c.id_, model::cells[p.coord(j).cell]->id_, p.coord(j).universe));
                        }
                        overlap_detected = true;
                        model::overlap_check_count[found_index]++;
                    }
                }
            }
        }

        return overlap_detected;
    }

    // 其余代码保持不变...

} // namespace openmc
