# 章节 8：Lattice、Universe 与递归坐标管理

本章节解释 `Lattice` 如何与 `Universe`、`GeometryState` 的层级栈协作，实现多级几何的粒子定位。

## Lattice 成员复述

- `LatticeType type_`：区分正交与六角晶格。
- `std::array<int, 3> shape_`：晶格尺寸。
- `std::vector<UniverseIndex> universes_`：每个格点关联的宇宙索引。
- `Position origin_`、`std::array<double, 3> pitch_`：定义晶格的空间布局。
- `std::vector<double> lower_`, `upper_`：六角晶格额外的坐标范围定义。
- `UniverseIndex outer_`：粒子出界时的外围宇宙。【F:include/openmc/lattice.h†L32-L214】

## Universe 与 Lattice 的协同

- `Universe::find_cell(LocalCoord& coord)` 遍历 `cells_`，根据 `coord.r` 与单元区域判定粒子所在单元，必要时更新 `coord.cell`。【F:include/openmc/universe.h†L112-L215】
- 当 `Cell` 被晶格填充时，`GeometryState::current_coord` 的 `lattice` 与 `lattice_index` 成员会更新，指示下一层级为晶格坐标。【F:include/openmc/particle_data.h†L80-L175】
- `GeometryState::level_down()` 在进入晶格或子宇宙时复制 `LocalCoord`，保持父层级信息，方便粒子逆向退出时使用。

## 粒子推进流程

1. `Particle::event_advance` 计算到晶格边界的距离，必要时调整 `LocalCoord::lattice_index`。
2. 若跨越晶格单元，`GeometryState::update_local_coord` 根据 `Lattice::universes_` 设置新的 `UniverseIndex` 并推入层级栈。
3. `Universe::find_cell` 在新宇宙中继续查找具体 `Cell`，形成递归定位链条。【F:src/geometry.cpp†L340-L587】

## 设计优势

- 层级栈避免递归函数调用，改用数组复制，提升性能。
- 晶格索引成员与 `Lattice` 的 `shape_`、`pitch_` 一一对应，便于在 tallies 中统计格点信息。
