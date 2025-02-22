#include "openmc/geometry_aux.h"

#include <algorithm> // for std::max
#include <sstream>
#include <unordered_set>

#include <fmt/core.h>
#include <pugixml.hpp>

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/dagmc.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/lattice.h"
#include "openmc/material.h"
#include "openmc/settings.h"
#include "openmc/surface.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/filter_cell_instance.h"
#include "openmc/tallies/filter_distribcell.h"

namespace openmc {
void ssn_read_cells(pugi::xml_node node)
{
  // Count the number of cells.
  int n_cells = 0;
  for (pugi::xml_node cell_node : node.children("cell")) {
    n_cells++;
  }

  // Loop over XML cell elements and populate the array.
  model::cells.reserve(n_cells);
  for (pugi::xml_node cell_node : node.children("cell")) {
    model::cells.push_back(make_unique<CSGCell>(cell_node));
  }

  // Fill the cell map.
  for (int i = 0; i < model::cells.size(); i++) {
    int32_t id = model::cells[i]->id_;
    auto search = model::cell_map.find(id);
    if (search == model::cell_map.end()) {
      model::cell_map[id] = i;
    } else {
      fatal_error(
        fmt::format("Two or more cells use the same unique ID: {}", id));
    }
  }

  read_dagmc_universes(node);

  populate_universes();

  // Allocate the cell overlap count if necessary.
  if (settings::check_overlaps) {
    model::overlap_check_count.resize(model::cells.size(), 0);
  }

  if (model::cells.size() == 0) {
    fatal_error("No cells were found in the geometry.xml file");
  }
}

void ssn_read_geometry_xml()
{
  // Display output message
  write_message("Reading geometry XML file...", 5);

  // Check if geometry.xml exists
  std::string filename = settings::path_input + "geometry.xml";
  if (!file_exists(filename)) {
    fatal_error("Geometry XML file '" + filename + "' does not exist!");
  }

  // Parse settings.xml file
  pugi::xml_document doc;
  auto result = doc.load_file(filename.c_str());
  if (!result) {
    fatal_error("Error processing geometry.xml file.");
  }

  // Get root element
  pugi::xml_node root = doc.document_element();

  //   read_geometry_xml(root);
  read_surfaces(root);

  ssn_read_cells(root);
}
} // namespace openmc

int main()
{
  openmc::settings::path_input = std::string(
    "/home/ssn/ssn_mc/openmc/OpenMC_CPP_TESTS/openmc_develop_test/");
  openmc::ssn_read_geometry_xml();
  return 0;
}