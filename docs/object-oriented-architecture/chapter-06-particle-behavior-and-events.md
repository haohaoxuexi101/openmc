# 章节 6：Particle 事件流、抽样策略与计数体系（深度课程）

本章节面向 `Particle` 类的行为模型，详细说明事件循环如何串联几何、材料、截面与统计模块，并给出抽样与计数的完整流程。按七个课程单元展开，涵盖从初始化、推进、碰撞到死亡的全生命周期。

---

## 单元 6.1：事件驱动框架概览

1. **核心方法**：`event_calculate_xs`、`event_advance`、`event_cross_surface`、`event_collide`、`event_revive_from_secondary`、`event_death` 构成主循环，按顺序处理物理和几何事件。【F:include/openmc/particle.h†L66-L122】【F:src/particle.cpp†L200-L493】
2. **状态来源**：所有事件读写的字段均定义在 `ParticleData`/`GeometryState`，保证接口统一，便于 CPU/GPU 共用逻辑。
3. **流程控制**：外部调度器根据粒子存活状态 (`alive()`)、权重 (`wgt()`)、次级银行等决定循环终止或复活。

---

## 单元 6.2：初始化与 Source 绑定

1. **`from_source`**：将 `SourceSite` 中的坐标、方向、能量、时间、权重写入 `ParticleData` 成员，并调用 `GeometryState::init_from_r_u` 重置层级栈，确保坐标一致性。【F:include/openmc/particle.h†L58-L65】【F:include/openmc/particle_data.h†L223-L254】
2. **出生单元识别**：若源记录未指定 cell，`exhaustive_find_cell` 使用 `GeometryState` 栈搜索，写入 `cell_born()` 以供后续脉冲高度与 progeny 统计。【F:src/particle.cpp†L397-L439】【F:src/geometry.cpp†L210-L358】
3. **随机种子设置**：`ParticleData::seeds_` 通过全局工作索引或源权重初始化，确保多个进程的粒子轨迹可复现。

---

## 单元 6.3：截面计算与材料耦合

1. **材料索引**：`GeometryState::material()` 与 `material_last()` 在粒子进入新 cell 时更新，驱动 `event_calculate_xs` 选择正确材料对象。【F:include/openmc/particle_data.h†L336-L370】【F:src/particle.cpp†L200-L212】
2. **宏/微观截面更新**：
   - 对连续能量模式：调用 `Material::calculate_xs` 填充 `MacroXS` 并刷新 `NuclideMicroXS` 缓存；
   - 对多群模式：访问 `data::mg.macro_xs_` 并更新群索引 (`g()`, `g_last()`)。【F:include/openmc/material.h†L38-L200】【F:src/particle.cpp†L200-L212】
3. **温度缓存**：`GeometryState::sqrtkT()` 提供材料温度，必要时触发热化散射或概率表逻辑，确保与 `ParticleData` 截面缓存一致。

---

## 单元 6.4：推进、抽样与边界处理

1. **距离抽样**：
   - 通过 `distance_to_boundary(*this)` 获取最近几何边界，并将结果写入 `boundary()`；
   - 使用 `macro_xs().total` 与 `-log(prn(current_seed()))` 抽样碰撞距离，保存到 `collision_distance()`。【F:src/particle.cpp†L213-L270】
2. **运动更新**：逐层更新 `coord(j).r` 与全局时间 `time()`，如果超出时间截断则回退距离并置零权重。
3. **统计计数**：轨迹长度与 `keff` track-length 项根据移动距离立即累积，减少事后遍历成本。【F:src/particle.cpp†L250-L269】
4. **边界条件**：当 `boundary().surface != SURFACE_NONE` 时，`event_cross_surface` 决定粒子穿越表面或跨晶格：
   - **表面**：调用 `Surface::reflect`/`cross_*_bc` 应用真空、反射、周期边界；
   - **晶格**：`cross_lattice` 更新 `LocalCoord::lattice_i`、`coord_level`，保持宇宙栈正确。【F:src/particle.cpp†L272-L314】【F:src/geometry.cpp†L299-L358】【F:include/openmc/surface.h†L36-L119】【F:include/openmc/lattice.h†L88-L138】
5. **Cell 重新定位**：更新完层级后，再次调用 `geometry::find_cell` 确认粒子位于新单元，写入 `cell_last()` 与 `material()` 等字段，保证后续截面计算连续。【F:src/geometry.cpp†L210-L411】

---

## 单元 6.5：碰撞与反应抽样

1. **反应选择**：根据粒子类型调用 `physics::sample_*_reaction`，流程包括：
   - `sample_nuclide` 基于材料中各 nuclide 的总截面抽样；
   - `sample_fission` 或其他反应根据微观截面 CDF 采样具体反应道；
   - `sample_secondary_photons`、`elastic_scatter` 等函数生成次级粒子或更新方向、能量。【F:include/openmc/physics.h†L24-L104】【F:src/physics.cpp†L91-L1163】
2. **能量/方向更新**：`ParticleData` 的 `E()`、`mu_`、`u()` 被写入新值，同时 `E_last()`、`wgt_last()` 保存碰撞前状态，供 tally 与脉冲高度使用。【F:include/openmc/particle_data.h†L519-L588】
3. **统计贡献**：
   - `keff_tally_collision()` 用宏观裂变截面与权重计算碰撞估计；
   - `score_collision_tally`、`score_pulse_height_tally` 等函数读取 `GeometryState` 的 cell/surface 信息记录统计。【F:src/particle.cpp†L316-L493】
4. **次级粒子生成**：通过 `create_secondary` 或 `secondary_bank()` 缓存新粒子，当主粒子死亡时 `event_revive_from_secondary` 依次取出，继承层级与统计信息。【F:include/openmc/particle.h†L42-L57】【F:src/particle.cpp†L397-L452】

---

## 单元 6.6：复活、分裂与权重窗口

1. **事件计数**：`n_event()` 自增，超过阈值触发警告并置零权重，防止无限循环。【F:include/openmc/particle_data.h†L631-L637】【F:src/particle.cpp†L397-L405】
2. **次级复活**：当主粒子死亡且银行非空时，`from_source` 加载银行尾部粒子，重置 `bank_second_E()` 与计数器，保持统计独立性。【F:src/particle.cpp†L397-L452】
3. **分裂/轮盘**：`split` 与 `physics::split_particle` 使用 `wgt()`、`ww_factor()` 控制粒子数量，在 `event_cross_surface` 或 `event_collide` 中调用，依赖 `ParticleData` 的 RNG 流保障公平抽样。【F:include/openmc/particle.h†L52-L57】【F:include/openmc/particle_data.h†L633-L640】【F:include/openmc/physics.h†L99-L104】

---

## 单元 6.7：死亡阶段与全局归约

1. **脉冲高度与轨迹输出**：若启用跟踪，`event_death` 写入最终轨迹、调用 `score_pulse_height_tally`，利用 `GeometryState` 的 cell 层级映射统计位置。【F:src/particle.cpp†L455-L520】
2. **全局计数归约**：通过 OpenMP 原子将 `keff_tally_*` 累积到全局变量，并重置本地计数器，确保后续粒子不会混淆贡献。【F:src/particle.cpp†L455-L492】
3. **Progeny 记录**：在临界计算中，根据 `id()` 计算偏移写入 `simulation::progeny_per_particle`，为裂变银行排序提供输入。【F:src/particle.cpp†L486-L492】
4. **DAGMC 支持**：若启用 CAD 几何，死亡时重置 `history()`，避免下一粒子继承错误的射线历史。【F:src/particle.cpp†L455-L465】

> **课程总结**：`Particle` 的事件机制以 `GeometryState`/`ParticleData` 为数据骨架，将抽样、计数与几何遍历紧密耦合。掌握该流程，可以在不破坏整体结构的前提下扩展新物理或优化抽样算法。
