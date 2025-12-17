#include <chrono>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>
#include "openmc/capi.h"

namespace {

constexpr double kBoltzmann = 8.617333262145e-5; // eV/K

void write_file(const std::string& path, const std::string& content)
{
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Unable to open file: " + path);
  }
  out << content;
}

void write_materials_xml(const std::string& xs_path)
{
  std::string xml = R"XML(<?xml version="1.0"?>
<materials cross_sections=")XML" + xs_path + R"XML(">
  <material id="1" name="UO2 fuel">
    <density units="macro" value="1.0"/>
    <macroscopic name="UO2"/>
  </material>
  <material id="2" name="Light water moderator">
    <density units="macro" value="1.0"/>
    <macroscopic name="LWTR"/>
  </material>
</materials>
)XML";
  write_file("materials.xml", xml);
}

void write_geometry_xml()
{
  std::string xml = R"XML(<?xml version="1.0"?>
<geometry>
  <surface id="1" type="xplane" x0="-1.89" boundary="reflective"/>
  <surface id="2" type="xplane" x0="1.89" boundary="reflective"/>
  <surface id="3" type="yplane" y0="-1.89" boundary="reflective"/>
  <surface id="4" type="yplane" y0="1.89" boundary="reflective"/>
  <surface id="5" type="zplane" z0="-150" boundary="reflective"/>
  <surface id="6" type="zplane" z0="150" boundary="reflective"/>
  <surface id="7" type="zcylinder" r="0.54"/>

  <cell id="1" material="1" region="-7" name="fuel pellet"/>
  <cell id="2" material="2" region="+7" name="moderator"/>
  <cell id="3" fill="10" region="+1 -2 +3 -4 +5 -6" name="root"/>

  <lattice id="10" name="3D multi-region core">
    <dimension>3 3 3</dimension>
    <lower_left>-1.89 -1.89 -150.0</lower_left>
    <width>1.26 1.26 100.0</width>
    <universes>
      <zplane>
        2 2 2
        2 2 2
        2 2 2
      </zplane>
      <zplane>
        1 1 1
        1 2 1
        1 1 1
      </zplane>
      <zplane>
        2 2 2
        2 2 2
        2 2 2
      </zplane>
    </universes>
  </lattice>
</geometry>
)XML";
  write_file("geometry.xml", xml);
}

void write_settings_xml(int seed)
{
  std::string xml = R"XML(<?xml version="1.0"?>
<settings>
  <energy_mode>multi-group</energy_mode>
  <seed>)XML" + std::to_string(seed) + R"XML(</seed>
  <batches>40</batches>
  <inactive>5</inactive>
  <particles>10000</particles>
  <run_mode>eigenvalue</run_mode>
  <source>
    <space type="box">
      <parameters>-1.89 -1.89 -150 1.89 1.89 150</parameters>
    </space>
    <constraints>
      <fissionable>true</fissionable>
    </constraints>
  </source>
</settings>
)XML";
  write_file("settings.xml", xml);
}

void write_tallies_xml()
{
  std::string xml = R"XML(<?xml version="1.0"?>
<tallies>
  <tally id="1" name="region-averaged flux">
    <filters>
      <cellfilter>1 2</cellfilter>
    </filters>
    <scores>flux fission nu-fission</scores>
  </tally>
</tallies>
)XML";
  write_file("tallies.xml", xml);
}

hid_t create_group(hid_t parent, const std::string& name)
{
  hid_t group = H5Gcreate(parent, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (group < 0) {
    throw std::runtime_error("Failed to create group " + name);
  }
  return group;
}

void write_string_attribute(hid_t obj, const std::string& name, const std::string& value)
{
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, value.size() + 1);
  H5Tset_strpad(type, H5T_STR_NULLTERM);
  hid_t attr = H5Acreate(obj, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, type, value.c_str());
  H5Aclose(attr);
  H5Sclose(space);
  H5Tclose(type);
}

void write_int_attribute(hid_t obj, const std::string& name, int value)
{
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t attr = H5Acreate(obj, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, &value);
  H5Aclose(attr);
  H5Sclose(space);
}

void write_double_attribute(hid_t obj, const std::string& name, double value)
{
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t attr = H5Acreate(obj, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_DOUBLE, &value);
  H5Aclose(attr);
  H5Sclose(space);
}

void write_bool_attribute(hid_t obj, const std::string& name, bool value)
{
  int v = value ? 1 : 0;
  write_int_attribute(obj, name, v);
}

void write_int_array_attribute(hid_t obj, const std::string& name, const std::vector<int>& values)
{
  hsize_t dims[1] {values.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t attr = H5Acreate(obj, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT, values.data());
  H5Aclose(attr);
  H5Sclose(space);
}

void write_double_array_attribute(hid_t obj, const std::string& name, const std::vector<double>& values)
{
  hsize_t dims[1] {values.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t attr = H5Acreate(obj, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_DOUBLE, values.data());
  H5Aclose(attr);
  H5Sclose(space);
}

void write_dataset_1d(hid_t parent, const std::string& name, const std::vector<double>& data)
{
  hsize_t dims[1] {data.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t dset = H5Dcreate(parent, name.c_str(), H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  H5Dclose(dset);
  H5Sclose(space);
}

void write_dataset_1d_int(hid_t parent, const std::string& name, const std::vector<int>& data)
{
  hsize_t dims[1] {data.size()};
  hid_t space = H5Screate_simple(1, dims, nullptr);
  hid_t dset = H5Dcreate(parent, name.c_str(), H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
  H5Dclose(dset);
  H5Sclose(space);
}

void write_mgxs_block(hid_t parent, const std::string& name, bool fissionable,
  const std::vector<double>& total, const std::vector<double>& absorption,
  const std::vector<double>& fission, const std::vector<double>& nu_fission,
  const std::vector<double>& chi, const std::vector<double>& scatter)
{
  hid_t group = create_group(parent, name);

  write_double_attribute(group, "atomic_weight_ratio", 1.0);
  write_bool_attribute(group, "fissionable", fissionable);
  write_string_attribute(group, "representation", "isotropic");
  write_string_attribute(group, "scatter_format", "legendre");
  write_int_attribute(group, "order", 0);
  write_string_attribute(group, "scatter_shape", "[Order][G][G']");

  hid_t kts = create_group(group, "kTs");
  const double kT = 293.6 * kBoltzmann;
  write_dataset_1d(kts, "293K", {kT});
  H5Gclose(kts);

  hid_t temperature = create_group(group, "293K");
  write_dataset_1d(temperature, "total", total);
  write_dataset_1d(temperature, "absorption", absorption);
  if (fissionable) {
    write_dataset_1d(temperature, "fission", fission);
    write_dataset_1d(temperature, "nu-fission", nu_fission);
    write_dataset_1d(temperature, "chi", chi);
  }

  hid_t scatter_group = create_group(temperature, "scatter_data");
  const int G = static_cast<int>(total.size());
  std::vector<int> gmin(G, 1);
  std::vector<int> gmax(G, G);
  write_dataset_1d_int(scatter_group, "g_min", gmin);
  write_dataset_1d_int(scatter_group, "g_max", gmax);
  write_dataset_1d(scatter_group, "scatter_matrix", scatter);
  H5Gclose(scatter_group);

  H5Gclose(temperature);
  H5Gclose(group);
}

void write_mgxs_library(const std::string& path)
{
  hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    throw std::runtime_error("Could not create mgxs library at " + path);
  }

  write_string_attribute(file, "filetype", "mgxs");
  write_int_array_attribute(file, "version", {1, 0});

  const std::vector<double> group_edges {
    1e-5, 0.0635, 10.0, 1.0e2, 1.0e3, 0.5e6, 1.0e6, 20.0e6};
  write_int_attribute(file, "energy_groups", static_cast<int>(group_edges.size() - 1));
  write_int_attribute(file, "delayed_groups", 0);
  write_double_array_attribute(file, "group structure", group_edges);

  const std::vector<double> uo2_total {0.1779492, 0.3298048, 0.4803882, 0.5543674, 0.3118013, 0.3951678, 0.5644058};
  const std::vector<double> uo2_abs {8.0248E-03, 3.7174E-03, 2.6769E-02, 9.6236E-02, 3.0020E-02, 1.1126E-01, 2.8278E-01};
  const std::vector<double> uo2_fiss {7.21206E-03, 8.19301E-04, 6.45320E-03, 1.85648E-02, 1.78084E-02, 8.30348E-02, 2.16004E-01};
  const std::vector<double> uo2_nu_fiss {2.005998E-02, 2.027303E-03, 1.570599E-02, 4.518301E-02, 4.334208E-02, 2.020901E-01, 5.257105E-01};
  const std::vector<double> uo2_chi {5.8791E-01, 4.1176E-01, 3.3906E-04, 1.1761E-07, 0.0, 0.0, 0.0};

  const std::vector<double> uo2_scatter {
    0.1275370, 0.0423780, 0.0000094, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.3244560, 0.0016314, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.4509400, 0.0026792, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.4525650, 0.0055664, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0001253, 0.2714010, 0.0102550, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0012968, 0.2658020, 0.0168090,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0085458, 0.2730800};

  const std::vector<double> h2o_total {0.15920605, 0.412969593, 0.59030986, 0.58435, 0.718, 1.2544497, 2.650379};
  const std::vector<double> h2o_abs {6.0105E-04, 1.5793E-05, 3.3716E-04, 1.9406E-03, 5.7416E-03, 1.5001E-02, 3.7239E-02};
  const std::vector<double> h2o_scatter {
    0.0444777, 0.1134000, 0.0007235, 0.0000037, 0.0000001, 0.0, 0.0,
    0.0, 0.2823340, 0.1299400, 0.0006234, 0.0000480, 0.0000074, 0.0000010,
    0.0, 0.0, 0.3452560, 0.2245700, 0.0169990, 0.0026443, 0.0005034,
    0.0, 0.0, 0.0, 0.0910284, 0.4155100, 0.0637320, 0.0121390,
    0.0, 0.0, 0.0, 0.0000714, 0.1391380, 0.5118200, 0.0612290,
    0.0, 0.0, 0.0, 0.0, 0.0022157, 0.6999130, 0.5373200,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.1324400, 2.4807000};

  write_mgxs_block(file, "UO2", true, uo2_total, uo2_abs, uo2_fiss, uo2_nu_fiss, uo2_chi, uo2_scatter);
  write_mgxs_block(file, "LWTR", false, h2o_total, h2o_abs, {}, {}, {}, h2o_scatter);

  H5Fclose(file);
}

} // namespace

int main()
{
  try {
    const std::string xs_path = "mgxs.h5";
    write_mgxs_library(xs_path);
    write_materials_xml(xs_path);
    write_geometry_xml();
    const auto seed = static_cast<int>(std::chrono::system_clock::now().time_since_epoch().count() % 1000000);
    write_settings_xml(seed);
    write_tallies_xml();

    openmc_set_seed(seed);
    int err = openmc_init(0, nullptr, nullptr);
    if (err != 0) {
      throw std::runtime_error("openmc_init failed with error code " + std::to_string(err));
    }

    err = openmc_run();
    if (err != 0) {
      throw std::runtime_error("openmc_run failed with error code " + std::to_string(err));
    }

    err = openmc_finalize();
    if (err != 0) {
      throw std::runtime_error("openmc_finalize failed with error code " + std::to_string(err));
    }

    std::cout << "Simulation complete. Statepoints and tallies are now available." << std::endl;
  } catch (const std::exception& ex) {
    std::cerr << "Fatal error: " << ex.what() << std::endl;
    return 1;
  }

  return 0;
}
