# Object-Oriented Particle Mini-Toolkit

This teaching example condenses the documentation chapters on `GeometryState`,
`ParticleData`, and `Particle` into ten self-contained C++ programs. Each
chapter builds on a shared miniature toolkit implemented in
`chapter_common.hpp`/`chapter_common.cpp`. The toolkit now wraps the actual
OpenMC C++ types — `openmc::GeometryState`, `openmc::ParticleData`,
`openmc::Particle`, and the geometry primitives from `openmc::Surface` through
`openmc::Universe` — so every executable exercises the real data structures
used in production transport. The programs illustrate how geometry objects,
particle state, materials, and tallies cooperate during a Monte Carlo transport
walk.

## Building

```bash
cmake -S . -B build -DOPENMC_BUILD_TESTS=OFF
cmake --build build --target tutorial_chapter_01
```

Replace the target with any of the chapter executables listed below. Binaries
are placed under `build/bin`.

## Chapters

| Chapter | Source file | Focus |
| --- | --- | --- |
| 1 | `chapter_01_geometry_walk.cpp` | Inspecting geometry stacks and boundary distance |
| 2 | `chapter_02_material_sampling.cpp` | Material interaction sampling |
| 3 | `chapter_03_surface_tracking.cpp` | Streaming toward surfaces |
| 4 | `chapter_04_collision_statistics.cpp` | Coupled collision/boundary logic |
| 5 | `chapter_05_lattice_navigation.cpp` | Lattice addressing |
| 6 | `chapter_06_rng_management.cpp` | Reproducible particle RNG |
| 7 | `chapter_07_tally_coupling.cpp` | Tally updates tied to events |
| 8 | `chapter_08_branching_histories.cpp` | Fission weight tracking |
| 9 | `chapter_09_directional_bias.cpp` | Direction management and biasing |
| 10 | `chapter_10_full_demo.cpp` | Integrated transport mini-loop |

Each program can be executed independently and prints its own instructional
narrative to standard output.

