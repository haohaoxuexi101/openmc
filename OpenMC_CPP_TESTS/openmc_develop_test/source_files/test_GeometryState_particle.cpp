#include "openmc/simulation.h"

#include "openmc/bank.h"
#include "openmc/capi.h"
#include "openmc/container_util.h"
#include "openmc/cross_sections.h"
// #include "openmc/distribution.h"
// #include "openmc/distribution_spatial.h"
#include "openmc/eigenvalue.h"
#include "openmc/error.h"
#include "openmc/event.h"
#include "openmc/geometry.h"
#include "openmc/geometry_aux.h"
#include "openmc/initialize.h"
#include "openmc/material.h"
#include "openmc/mcpl_interface.h"
#include "openmc/message_passing.h"
#include "openmc/nuclide.h"
#include "openmc/output.h"
#include "openmc/particle.h"
// #include "openmc/particle_data.h"
#include "openmc/photon.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/source.h"
#include "openmc/state_point.h"
#include "openmc/tallies/derivative.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/tallies/trigger.h"
#include "openmc/timer.h"
#include "openmc/track_output.h"
#include "openmc/weight_windows.h"
// #include "openmc/xml_interface.h"

int main()
{
  // openmc::model::n_coord_levels = 2;
  // 如果需要完整的初始化变量，那么就要使用openmc::read_separate_xml_files()
  // 或者将其中的函数抽离出来，单独调用
  // 但是单独定义变量是不行的
  // 材料几何截面参数无法独立存在？
  openmc::settings::path_input = std::string(
    "/home/ssn/ssn_mc/openmc/OpenMC_CPP_TESTS/openmc_develop_test/");

  class MyParticle : public openmc::Particle {
  public:
    void get_xyz(double* xyz) const
    {
      xyz[0] = r().x;
      xyz[1] = r().y;
      xyz[2] = r().z;
      std::cout << "xyz: " << xyz[0] << " " << xyz[1] << " " << xyz[2]
                << std::endl;
    }
  };
  // 源抽样错误的首要原因在于，没有提供正确的model.xml
  // 提供的是seperate的xml文件（特征值模式），而不是model.xml（固定源模式）
  // 优于xml文件的误用，导致用了固定源的函数，而没有使用特征值模式的函数

  // fixed source模式下，从外部源中抽样一个粒子
  // if (openmc::read_model_xml()) {
  //   std::cout << "read the model xml" << std::endl;
  // }
  // openmc::Particle p;
  // MyParticle p;
  // uint64_t seed = 1;
  // auto site = openmc::sample_external_source(&seed);
  // p.from_source(&site);

  // eigenvalue模式下，从source_bank中取出一个粒子
  // 为什么需要read_settings_xml()？
  // 必须read并且finalize，否则有问题？
  openmc::read_separate_xml_files();
  // openmc::read_settings_xml();
  // openmc::read_cross_sections_xml();
  // openmc::read_materials_xml();
  // openmc::read_geometry_xml();

  // openmc::finalize_geometry();
  // openmc::finalize_cross_sections();

  openmc_simulation_init();

  // openmc::Particle p;
  MyParticle p;

  // openmc_simulation_init();

  // openmc::initialize_source();
  openmc::initialize_history(p, 1);

  // p.from_source(&openmc::simulation::source_bank[0]);

  // double xyz[3];
  // p.get_xyz(xyz);
  std::cout << p.r().x << " " << p.r().y << " " << p.r().z << std::endl;

  // bool ifoverlap = openmc::check_cell_overlap(p, true);
  // if (ifoverlap) {
  //   std::cout << "overlap" << std::endl;
  // } else {
  //   std::cout << "not overlap" << std::endl;
  // }
  bool foundif = openmc::exhaustive_find_cell(p);
  if (!foundif) {
    std::cout << "not found" << std::endl;
  } else {
    std::cout << "found" << std::endl;
  }

  openmc::BoundaryInfo bd;

  bd = openmc::distance_to_boundary(p);

  std::cout << bd.distance << std::endl;

  return 0;
}