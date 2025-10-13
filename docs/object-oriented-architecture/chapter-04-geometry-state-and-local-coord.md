# 章节 4：GeometryState 与 LocalCoord 层级跟踪

本章节详解几何上下文结构 `GeometryState` 与 `LocalCoord`，说明其成员字段如何协同跟踪粒子在嵌套宇宙中的位置。

## LocalCoord 成员

- `Position r`：当前位置的局部坐标。
- `Direction u`：当前方向向量，适用于晶格内局部坐标系。
- `CellIndex cell`、`UniverseIndex universe`、`LatticeIndex lattice`：标识粒子所在实体。
- `std::array<int, 3> lattice_index`：记录在晶格中的 (`i`, `j`, `k`) 索引。
- `int level`：当前层级深度。【F:include/openmc/particle_data.h†L80-L175】

## GeometryState 成员

- `LocalCoord current_coord`：指向 `coord` 数组中的当前层级。
- `std::array<LocalCoord, MAX_COORD> coord`：固定长度的层级栈，避免动态分配。
- `int n_coord`：层级深度计数。
- `SurfaceIndex surface_last`：记录上一次穿越的表面。
- `CellIndex cell_last`：记录上一个所在单元。
- `int level_last`：辅助回溯时的层级信息。
- `bool on_boundary`：标识粒子是否位于几何边界，供边界条件逻辑使用。【F:include/openmc/geometry_state.h†L30-L118】

## 协同流程

- `GeometryState::level_down()`/`level_up()` 修改 `n_coord` 并更新 `current_coord` 指针，用于进入或退出子宇宙。【F:include/openmc/geometry_state.h†L86-L118】
- `Particle::cross_surface` 调用 `GeometryState::update_local_coord()` 刷新 `LocalCoord` 中的 `cell` 与 `universe` 成员，确保与几何层级同步。【F:src/particle.cpp†L414-L611】
- `GeometryState` 与 `ParticleData::cell_instance_id` 等字段一同使用，避免在统计与 tallies 中丢失层级信息。【F:include/openmc/particle_data.h†L180-L257】

## 设计要点

- 栈式 `coord` 数组通过常量 `MAX_COORD` 控制最大嵌套深度，保证缓存友好。
- 成员 `surface_last`、`cell_last` 为 tallies 和边界条件提供上下文，可在粒子推进时快速判断重复穿越或新进入单元。
