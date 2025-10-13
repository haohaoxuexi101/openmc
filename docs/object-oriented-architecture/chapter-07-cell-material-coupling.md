# 章节 7：Cell、Material 与 Particle 的耦合路径

本章节描述粒子进入单元后如何通过几何和材料模块获取核数据，列出关键成员之间的关联。

## 数据流步骤

1. `GeometryState::current_coord.cell` 指出粒子所在 `Cell` 索引。【F:include/openmc/geometry_state.h†L30-L118】
2. 全局容器 `model::cells`（在 `src/geometry.cpp` 初始化）根据索引返回 `Cell` 对象，读取其 `fill_` 成员判定填充类型。【F:src/geometry.cpp†L210-L411】【F:include/openmc/cell.h†L124-L305】
3. 若 `fill_` 指向 `Material`，通过 `Cell::material_` 及温度映射数组获取材料索引；若指向 `Universe` 或 `Lattice`，则更新 `GeometryState` 进入下一层级。
4. `ParticleData::material` 存储材料句柄，供 `event_calculate_xs` 调用 `Material::calculate_xs`。【F:include/openmc/particle_data.h†L120-L204】【F:include/openmc/material.h†L170-L265】

## 关键成员说明

- `Cell::fill_`：`struct Fill { FillType type; int32_t id; };` 用于区分材料、宇宙或晶格填充。【F:include/openmc/cell.h†L180-L245】
- `Cell::material_`：`gsl::span<MaterialIndex>` 存储可用材料索引，结合 `sqrtkT_` 匹配温度。【F:include/openmc/cell.h†L217-L305】
- `Material::nuclides_`：`std::vector<int>` 列出核素索引；`Material::atom_density_` 保存核素原子密度。【F:include/openmc/material.h†L120-L265】
- `ParticleData::macro_xs` 与 `ParticleData::micro_xs`：承载由 `Material` 与 `Nuclide` 计算的截面数据，用于后续碰撞采样。【F:include/openmc/particle_data.h†L167-L257】

## 耦合特性

- `Particle` 并不直接持有 `Cell` 或 `Material` 对象，而是通过索引与全局模型解耦，便于重建几何或材料数据库。
- 材料温度和密度变化仅需更新 `Material` 成员即可；`ParticleData` 通过 `last_material` 与 `last_E` 判断缓存是否有效，避免重复计算。
