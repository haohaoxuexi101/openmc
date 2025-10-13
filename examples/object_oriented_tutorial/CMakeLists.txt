set(tutorial_sources
  chapter_common.cpp
)

add_library(tutorial_common STATIC ${tutorial_sources})
target_include_directories(tutorial_common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

set(chapters
  01_geometry_walk
  02_material_sampling
  03_surface_tracking
  04_collision_statistics
  05_lattice_navigation
  06_rng_management
  07_tally_coupling
  08_branching_histories
  09_directional_bias
  10_full_demo
)

foreach(chapter IN LISTS chapters)
  set(target_name tutorial_chapter_${chapter})
  add_executable(${target_name} chapter_${chapter}.cpp)
  target_link_libraries(${target_name} PRIVATE tutorial_common)
endforeach()
