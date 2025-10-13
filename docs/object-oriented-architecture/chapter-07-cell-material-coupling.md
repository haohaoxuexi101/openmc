# 章节 7：Cell、Material 与粒子状态的多层耦合

本章节在既有几何栈的基础上，深入分析 `Cell`、`Universe`、`Lattice` 填充策略如何驱动材料选择、温度映射与截面计算。通过五个单元展开，涵盖成员说明、数据流与扩展建议。

---

## 单元 7.1：Cell 对象的结构化描述

1. **基础字段**：`Cell::id_`、`name_`、`universe_` 定义唯一身份；`fill_` 描述填充类型与目标 ID，可指向 `Material`、`Universe` 或 `Lattice`。【F:include/openmc/cell.h†L150-L320】
2. **材料表**：`vector<int32_t> material_` 与 `vector<double> sqrtkT_` 存储分布式材料及温度实例，支持每个 cell instance 拥有独立核数据。【F:include/openmc/cell.h†L250-L283】
3. **几何接口**：`contains`、`distance`、`surfaces` 提供几何判定；`find_parent_cells`、`get_contained_cells` 支持层级遍历与 distribcell 构建。【F:include/openmc/cell.h†L162-L308】

---

## 单元 7.2：粒子进入 Cell 的数据流

1. **定位**：`GeometryState::lowest_coord().cell` 保存当前 cell 索引，由 `geometry::find_cell` 计算，必要时遍历 `Universe` 与 `Lattice` 树。【F:include/openmc/particle_data.h†L270-L353】【F:src/geometry.cpp†L210-L411】
2. **填充判定**：
   - 若 `fill_.type == Fill::MATERIAL`，直接返回材料索引；
   - 若填充 `Universe`/`Lattice`，调用 `GeometryState::level_down()` 推入新层级，并继续递归搜索，直至找到真实材料。【F:include/openmc/cell.h†L180-L305】【F:src/geometry.cpp†L299-L358】
3. **材料/温度映射**：
   - `cell_instance()` 提供 distribcell 偏移，结合 `Cell::material(instance)` 选择材料索引；
   - `Cell::sqrtkT(instance)` 提供温度，写入 `GeometryState::sqrtkT()`，驱动材料截面插值。【F:include/openmc/particle_data.h†L266-L374】【F:include/openmc/cell.h†L250-L283】
4. **缓存同步**：`ParticleData::material()` 与 `material_last()` 保留当前/上一次材料索引，用于判断截面缓存是否有效，从而避免重复计算。【F:include/openmc/particle_data.h†L336-L373】

---

## 单元 7.3：Material 对象与核数据准备

1. **核素与密度**：`Material::nuclide_`、`atom_density_`、`density_`、`density_gpcc_` 定义核素组成与宏观密度；`mat_nuclide_index_` 提供快速索引，以 O(1) 速度查找目标核素。【F:include/openmc/material.h†L38-L200】
2. **热化与截面计算**：
   - `init_thermal()` 与 `ThermalTable` 结构关联热中子散射表；
   - `calculate_xs(Particle& p)` 读取 `ParticleData` 中的能量、方向、温度，填充宏观截面并触发微观缓存更新。【F:include/openmc/material.h†L38-L200】
3. **几何依赖**：材料对象本身不感知几何，但依赖 `Particle` 提供的 `material()`、`sqrtkT()`。通过这种反向依赖，实现几何-材料解耦同时保持高性能。

---

## 单元 7.4：抽样、统计与材料的交互

1. **反应抽样**：`physics::sample_nuclide`、`sample_fission` 等函数使用 `ParticleData::macro_xs()`、`neutron_xs()`，这些数据源自当前 `Material` 的计算结果。【F:include/openmc/physics.h†L24-L104】【F:src/physics.cpp†L91-L569】
2. **统计计数**：tally 模块通过 `GeometryState::cell_last()`、`material_last()` 区分不同 cell/material 组合，实现分布式材料统计；`Particle::event_collide` 利用 `macro_xs().nu_fission` 贡献碰撞估计。【F:include/openmc/particle_data.h†L279-L373】【F:src/particle.cpp†L316-L493】
3. **权重调整**：在材料吸收强烈的区域，权重窗口使用 `ww_factor()`、`split_particle` 根据当前材料的泄漏特性调整粒子数量，依赖材料提供的宏观截面上下界。【F:include/openmc/particle_data.h†L633-L640】【F:include/openmc/physics.h†L99-L104】

---

## 单元 7.5：扩展与调试建议

1. **材料切换**：若引入温度依赖材料或燃耗模拟，应确保在更新材料数据后调用 `ParticleData::invalidate_neutron_xs()`，以强制刷新缓存并避免使用过期截面。【F:include/openmc/particle_data.h†L649-L655】
2. **几何诊断**：在复杂嵌套中，可利用 `GeometryState::mark_as_lost()` 打印粒子 ID、坐标与层级，结合 `Cell::name()`、`Material::name()` 快速定位问题。【F:include/openmc/particle_data.h†L212-L239】【F:include/openmc/cell.h†L233-L240】【F:include/openmc/material.h†L110-L124】
3. **多物理耦合**：扩展燃耗或热工耦合时，建议通过材料 ID/instance 与外部求解器同步更新密度、温度，再触发 `Material::finalize()`，确保粒子侧读取的是最新数据。【F:include/openmc/material.h†L48-L200】

> **课程总结**：`Cell`/`Material` 之间的关系通过 `GeometryState` 与 `ParticleData` 的索引桥接，实现了几何、材料、统计三者的松耦合而高效率的协作，为复杂物理模拟提供了可扩展框架。
