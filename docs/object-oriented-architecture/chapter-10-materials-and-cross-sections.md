# 章节 10：材料、核素与截面缓存

本章节梳理 `Material`、`Nuclide`、`Reaction` 与 `ParticleData` 缓存结构之间的成员关系。

## Material 成员

- `int32_t id_`、`std::string name_`：材料标识。
- `std::vector<int> nuclides_`：核素索引列表。
- `std::vector<double> atom_density_`：对应核素的原子密度。
- `double density_`、`double volume_`：宏观密度与体积。
- `TemperatureRanges temps_`：温度区间与插值信息。
- `bool fissionable_`：是否可裂变。
- 方法 `calculate_xs(ParticleData&, const CrossSections&)` 根据粒子状态计算宏观截面。【F:include/openmc/material.h†L42-L265】

## Nuclide 成员

- `int32_t id_`、`std::string name_`：核素标识。
- `EnergyGrid unionized_grid_`：对数联合能量网格。
- `Reaction xs_[MAX_REACTION]`：反应数据数组。
- `ThermalScattering* thermal_`：热化截面接口指针。
- 函数 `calculate_micro_xs`、`sample_scatter` 供粒子事件调用。【F:include/openmc/nuclide.h†L40-L322】

## Reaction 成员

- `ReactionType type_`：反应类型。
- `double q_value_`：能量释放。
- `Tabular data_`、`AngleDistribution angle_dist_`：截面与角度分布数据。
- `EnergyDistribution energy_dist_`：次级能量分布。【F:include/openmc/reaction.h†L34-L296】

## ParticleData 缓存字段

- `MacroXS macro_xs` 包含 `total`, `absorption`, `nu_fission`, `scatter`, `fission` 等成员。
- `MicroXS micro_xs` 存储每个核素的微观截面数组与散射矩阵索引。
- `NuclideMicroXS last_nuclide_xs` 保存最近一次核素截面及能量。
- `double last_E`、`MaterialHandle last_material`、`bool need_depletion_rx` 组合判定缓存失效条件。【F:include/openmc/particle_data.h†L167-L257】

## 数据流

1. `Particle::event_calculate_xs` 根据 `ParticleData::material` 获取 `Material`。
2. `Material::calculate_xs` 遍历 `nuclides_`，调用 `Nuclide::calculate_micro_xs` 并写入 `ParticleData::micro_xs`。
3. 累积结果写入 `ParticleData::macro_xs`，供 `event_collide` 采样反应类型。【F:src/material.cpp†L140-L391】【F:src/physics_common.cpp†L78-L287】

## 设计亮点

- 将截面缓存放置于 `ParticleData` 内部，减少跨模块内存访问。
- `Material` 的核素列表使用向量，便于燃耗更新时动态调整。
- `Nuclide` 将概率表与热化模型封装为成员指针，允许在运行时切换实现。
