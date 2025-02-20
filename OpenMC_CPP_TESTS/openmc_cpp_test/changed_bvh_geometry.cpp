#include "openmc/geometry.h"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <algorithm>
#include <vector>
#include <memory>
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
    // BVH (Bounding Volume Hierarchy) Implementation
    //==============================================================================

    // Structure representing an axis-aligned bounding box (AABB)
    struct AABB
    {
        double xmin;
        double xmax;
        double ymin;
        double ymax;
        double zmin;
        double zmax;

        // Check if a given position is within this bounding box
        bool contains(const Position &pos) const
        {
            return (pos.x >= xmin && pos.x <= xmax &&
                    pos.y >= ymin && pos.y <= ymax &&
                    pos.z >= zmin && pos.z <= zmax);
        }
    };

    // Class representing a node in the Bounding Volume Hierarchy (BVH)
    class BVHNode
    {
    public:
        AABB bounding_box;              // The bounding box that encapsulates this node
        std::unique_ptr<BVHNode> left;  // Pointer to the left child node
        std::unique_ptr<BVHNode> right; // Pointer to the right child node
        int cell_index = -1;            // Index of the cell for leaf nodes; -1 if not a leaf

        BVHNode() = default;
    };

    // Function to build a BVH tree from a list of cell indices and their bounding boxes
    std::unique_ptr<BVHNode> build_bvh(std::vector<int> &cell_indices, const std::vector<AABB> &cell_bounding_boxes, int start, int end)
    {
        if (start >= end)
            return nullptr;

        auto node = std::make_unique<BVHNode>();

        // Compute the bounding box for this node
        AABB combined_box = cell_bounding_boxes[cell_indices[start]];
        for (int i = start + 1; i < end; ++i)
        {
            const auto &box = cell_bounding_boxes[cell_indices[i]];
            combined_box.xmin = std::min(combined_box.xmin, box.xmin);
            combined_box.xmax = std::max(combined_box.xmax, box.xmax);
            combined_box.ymin = std::min(combined_box.ymin, box.ymin);
            combined_box.ymax = std::max(combined_box.ymax, box.ymax);
            combined_box.zmin = std::min(combined_box.zmin, box.zmin);
            combined_box.zmax = std::max(combined_box.zmax, box.zmax);
        }
        node->bounding_box = combined_box;

        // If there is only one cell, this is a leaf node
        if (end - start == 1)
        {
            node->cell_index = cell_indices[start];
            return node;
        }

        // Split the cell indices by the median along the x-axis for simplicity
        int mid = (start + end) / 2;
        std::nth_element(cell_indices.begin() + start, cell_indices.begin() + mid, cell_indices.begin() + end,
                         [&cell_bounding_boxes](int a, int b)
                         {
                             return cell_bounding_boxes[a].xmin < cell_bounding_boxes[b].xmin;
                         });

        // Recursively build the left and right subtrees
        node->left = build_bvh(cell_indices, cell_bounding_boxes, start, mid);
        node->right = build_bvh(cell_indices, cell_bounding_boxes, mid, end);

        return node;
    }

    // Function to query the BVH tree for a cell containing the given position
    bool query_bvh(const BVHNode *node, const Position &pos, int &found_cell)
    {
        if (!node || !node->bounding_box.contains(pos))
            return false;

        // If this is a leaf node, check if the position matches the cell
        if (node->cell_index != -1)
        {
            found_cell = node->cell_index;
            return true;
        }

        // Recursively query the left and right subtrees
        return query_bvh(node->left.get(), pos, found_cell) || query_bvh(node->right.get(), pos, found_cell);
    }

    //==============================================================================
    // Non-member functions
    //==============================================================================

    // Function to check if there is any cell overlap using the BVH structure
    bool check_cell_overlap(GeometryState &p, bool error)
    {
        int n_coord = p.n_coord();
        bool overlap_detected = false;

        // Create BVH and populate data
        std::vector<int> cell_indices;
        std::vector<AABB> cell_bounding_boxes;
        for (int j = 0; j < n_coord; j++)
        {
            Universe &univ = *model::universes[p.coord(j).universe];
            for (auto index_cell : univ.cells_)
            {
                cell_indices.push_back(index_cell);
                Cell &c = *model::cells[index_cell];
                AABB cell_box;
                cell_box.xmin = c.bounding_box().xmin;
                cell_box.xmax = c.bounding_box().xmax;
                cell_box.ymin = c.bounding_box().ymin;
                cell_box.ymax = c.bounding_box().ymax;
                cell_box.zmin = c.bounding_box().zmin;
                cell_box.zmax = c.bounding_box().zmax;
                cell_bounding_boxes.push_back(cell_box);
            }
        }

        // Build the BVH tree
        auto bvh_root = build_bvh(cell_indices, cell_bounding_boxes, 0, cell_indices.size());

        // Use BVH to query for cell overlaps
        for (int j = 0; j < n_coord; j++)
        {
            int found_index;
            if (query_bvh(bvh_root.get(), p.coord(j).r, found_index))
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

    // The rest of the code remains unchanged...

} // namespace openmc
