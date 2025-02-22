#include "openmc/geometry.h"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <algorithm>

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
    // KDTree Implementation
    //==============================================================================

    class KDTree
    {
    public:
        struct Node
        {
            Position position;
            int cell_index;
            Node *left;
            Node *right;
            Node(const Position &pos, int idx) : position(pos), cell_index(idx), left(nullptr), right(nullptr) {}
        };

        Node *root;

        KDTree() : root(nullptr) {}

        void build(std::vector<int> &cells, const std::vector<Position> &positions)
        {
            root = build_recursive(cells, positions, 0, cells.size(), 0);
        }

        Node *build_recursive(std::vector<int> &cells, const std::vector<Position> &positions, int start, int end, int depth)
        {
            if (start >= end)
                return nullptr;
            int mid = (start + end) / 2;
            auto comp = [depth, &positions](int a, int b)
            {
                return positions[a][depth % 3] < positions[b][depth % 3];
            };
            std::nth_element(cells.begin() + start, cells.begin() + mid, cells.begin() + end, comp);
            Node *node = new Node(positions[cells[mid]], cells[mid]);
            node->left = build_recursive(cells, positions, start, mid, depth + 1);
            node->right = build_recursive(cells, positions, mid + 1, end, depth + 1);
            return node;
        }

        int search(const Position &pos, Node *node, int depth)
        {
            if (!node)
                return -1;
            if (pos == node->position)
                return node->cell_index;
            if (pos[depth % 3] < node->position[depth % 3])
            {
                return search(pos, node->left, depth + 1);
            }
            else
            {
                return search(pos, node->right, depth + 1);
            }
        }

        int search(const Position &pos)
        {
            return search(pos, root, 0);
        }
    };

    //==============================================================================
    // Non-member functions
    //==============================================================================

    bool check_cell_overlap(GeometryState &p, bool error)
    {
        int n_coord = p.n_coord();
        bool overlap_detected = false;

        // 创建 KDTree 并填充数据
        std::vector<int> cell_indices;
        std::vector<Position> cell_positions;
        for (int j = 0; j < n_coord; j++)
        {
            Universe &univ = *model::universes[p.coord(j).universe];
            for (auto index_cell : univ.cells_)
            {
                cell_indices.push_back(index_cell);
                cell_positions.push_back(p.coord(j).r); // 假设每个单元的参考点
            }
        }

        KDTree cell_tree;
        cell_tree.build(cell_indices, cell_positions);

        // 使用 KDTree 查询
        for (int j = 0; j < n_coord; j++)
        {
            Universe &univ = *model::universes[p.coord(j).universe];
            int index_cell = cell_tree.search(p.coord(j).r);
            if (index_cell != -1)
            {
                Cell &c = *model::cells[index_cell];
                if (c.contains(p.coord(j).r, p.coord(j).u, p.surface()))
                {
                    if (index_cell != p.coord(j).cell)
                    {
                        if (error)
                        {
                            fatal_error(fmt::format("Overlapping cells detected: {}, {} on universe {}",
                                                    c.id_, model::cells[p.coord(j).cell]->id_, univ.id_));
                        }
                        overlap_detected = true;
                        model::overlap_check_count[index_cell]++;
                    }
                }
            }
        }

        return overlap_detected;
    }

    // 其余的代码保持不变...

} // namespace openmc
