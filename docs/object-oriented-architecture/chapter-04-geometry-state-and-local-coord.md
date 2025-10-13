# 章节 4：GeometryState 与 LocalCoord 层级跟踪（深度课程）

本章节聚焦 `GeometryState` 与 `LocalCoord` 在多层几何中的协同机制，系统说明它们的成员设计、与 `Particle`/`ParticleData` 的组合策略，以及与 `Cell`、`Universe`、`Lattice`、`Surface`、`Material` 的双向耦合。课程划分为六个学习单元，可单独阅读也可串联成体系化培训。

---

## 单元 4.1：层级坐标体系概览

1. **设计目标**：在嵌套宇宙和晶格内追踪粒子位置，确保不依赖 RTTI/虚函数即可在 CPU 与 GPU 之间共享数据布局。
2. **核心策略**：
   - `GeometryState` 以组合保存全局与局部坐标栈、上一层级信息以及与材料访问相关的索引。【F:include/openmc/particle_data.h†L208-L375】
   - `LocalCoord` 负责单层局部坐标描述，包括坐标、方向、隶属单元、宇宙、晶格索引，以及旋转标记，支持晶格旋转/反射等几何变换。【F:include/openmc/particle_data.h†L81-L95】
3. **运行时表现**：粒子推进时 `GeometryState` 像一座“多层记忆宫殿”，使粒子能够在宇宙/晶格树中下潜或回溯，而无需额外堆分配或昂贵的 map 查询。

---

## 单元 4.2：LocalCoord 成员级解析

| 成员 | 类型/含义 | 与其他对象的联系 |
|------|-----------|--------------------|
| `r` | `Position`：该层的空间位置 | 在 `GeometryState::coord_` 栈中复制，随粒子推进同步更新，用于 `Cell::contains` 与 `Surface::distance` 判定。【F:include/openmc/particle_data.h†L88-L95】【F:include/openmc/cell.h†L162-L186】 |
| `u` | `Direction`：该层方向 | 提供局部坐标系内方向，结合晶格旋转用于重新计算方向余弦。 |
| `cell`/`universe`/`lattice` | `int`：索引 | 直接映射到 `model::cells`、`model::universes`、`model::lattices` 的向量索引，保障材料和几何访问一致。【F:include/openmc/cell.h†L150-L320】【F:include/openmc/universe.h†L28-L48】【F:include/openmc/lattice.h†L46-L144】 |
| `lattice_i` | `array<int,3>`：晶格索引 | 与 `Lattice::get_indices` 互相更新，实现晶格跨越与周期边界处理。【F:include/openmc/lattice.h†L88-L114】 |
| `rotated` | `bool`：是否旋转 | 供 `LocalCoord::rotate` 与晶格旋转矩阵互操作，记录本层使用的旋转是否已应用。 |

*重置逻辑*：`LocalCoord::reset()` 会清空上述索引，确保粒子复活或回溯时不会残留旧层数据。【F:include/openmc/particle_data.h†L85-L95】

---

## 单元 4.3：GeometryState 成员深度剖析

1. **坐标栈**：
   - `vector<LocalCoord> coord_` 保存每个层级的 `LocalCoord`，结合 `n_coord_` 标记当前深度，既能应对深度递归又保证内存连续性。【F:include/openmc/particle_data.h†L270-L353】
   - `lowest_coord()`、`r_local()`、`u_local()` 等访问器用于获取最深层局部信息，在晶格和嵌套宇宙中至关重要。
2. **历史快照**：
   - `cell_last_` 与 `n_coord_last_` 记录穿越表面前的层级，用于 tally 和边界条件判断；`surface_`/`surface_index()` 表征当前边界。【F:include/openmc/particle_data.h†L279-L366】
   - `r_last_current_`、`r_last_`、`u_last_` 记录碰撞前后空间状态，供轨迹、碰撞核和加速算法使用。
3. **材料上下文**：
   - `material_`、`material_last_`、`sqrtkT_` 等字段提供材料索引与温度缓存，使 `Particle::event_calculate_xs` 可以直接调用 `Material::calculate_xs` 更新微/宏观截面。【F:include/openmc/particle_data.h†L336-L373】【F:include/openmc/material.h†L38-L200】
4. **边界预测**：
   - `BoundaryInfo boundary_` 缓存最近边界距离及跨越后层级，为 `Particle::event_advance` 与 `geometry::distance_to_boundary` 的联动提供接口。【F:include/openmc/particle_data.h†L186-L203】【F:src/particle.cpp†L213-L270】【F:src/geometry.cpp†L359-L411】
5. **唯一标识**：`id_` 与 `cell_instance_` 协助并行时的错误报告、分布式材料索引偏移和 progeny 统计。

---

## 单元 4.4：与几何实体的联动流程

1. **定位 Cell**：
   - 几何搜索通过 `geometry::find_cell`/`Universe::find_cell` 将粒子嵌入正确的 `Cell`，将索引写入 `lowest_coord().cell` 并更新材料索引。【F:include/openmc/universe.h†L28-L48】【F:src/geometry.cpp†L103-L411】
   - `Cell::distance` 使用 `GeometryState` 的位置与方向计算下一次表面交点，同时可利用 `GeometryState*` 指针获取栈信息。【F:include/openmc/cell.h†L183-L199】
2. **穿越 Surface**：
   - 当 `boundary().surface` 为真实表面时，`Particle::event_cross_surface` 会调用 `Surface::distance`、`Surface::reflect` 等函数，并更新 `GeometryState::surface_`、`cell_last_`、`n_coord_` 等字段。【F:include/openmc/surface.h†L36-L119】【F:src/particle.cpp†L272-L314】
3. **进入/退出 Lattice**：
   - `GeometryState::boundary().lattice_translation` 指示晶格跨越方向，`cross_lattice` 根据 `Lattice::distance` 与 `LocalCoord::lattice_i` 更新层级与局部坐标。【F:include/openmc/lattice.h†L88-L144】【F:src/geometry.cpp†L299-L358】
4. **Universe 栈管理**：
   - 当 `Cell::fill` 指向 `Universe`/`Lattice` 时，通过 `level_down()` 推入新的 `LocalCoord`，复制父层位置方向；退出时 `level_up()` 回溯并恢复父层信息，确保跨层 tally 与碰撞序列连续。【F:include/openmc/particle_data.h†L223-L378】【F:src/geometry.cpp†L210-L282】

---

## 单元 4.5：驱动物理模拟的数据流

1. **截面计算**：`Particle::event_calculate_xs` 基于 `GeometryState::material()` 查找材料，并通过 `Material::calculate_xs` 填充 `ParticleData` 中的交叉截面缓存，用于后续抽样。【F:src/particle.cpp†L200-L212】【F:include/openmc/material.h†L38-L200】
2. **距离抽样**：`Particle::event_advance` 调用 `distance_to_boundary(*this)` 计算几何距离，同时根据 `macro_xs().total` 与 `current_seed()` 采样碰撞距离，两者较小者决定下一事件类型。【F:src/particle.cpp†L213-L270】
3. **计数与统计**：
   - 轨迹长度、表面电流等 tally 直接读取 `GeometryState` 的 `n_coord_last_`、`surface_` 等字段，避免重新定位几何。【F:src/particle.cpp†L250-L314】【F:docs/object-oriented-architecture/chapter-11-tallies-and-sharing.md†L15-L33】
   - 当粒子死亡或穿越边界时，`GeometryState` 提供足够上下文以记录 `keff`、脉冲高度与 progeny 数据。
4. **权重与事件循环**：`GeometryState::id()` 与 `ParticleData::n_event()` 确保在极端情况下（如无限循环）能打印可追踪的错误信息并终止粒子。

---

## 单元 4.6：面向对象整合与扩展注意事项

1. **组合优先**：`Particle` 通过继承 `ParticleData`（进而继承 `GeometryState`）重用坐标与材料缓存，但不引入额外虚函数，从而保持 POD 特性，利于 GPU/CPU 共用数据布局。【F:include/openmc/particle.h†L23-L122】【F:include/openmc/particle_data.h†L416-L672】
2. **服务分离**：几何、材料、截面抽样分别通过非成员函数（例如 `geometry::distance_to_boundary`、`physics::collision`）访问 `GeometryState`，保证模块之间只通过轻量引用和索引交互。【F:src/geometry.cpp†L359-L411】【F:include/openmc/physics.h†L24-L104】
3. **扩展指南**：
   - 增加新的坐标层级字段时，需要同步更新 `clear()`、`level_down()` 等栈管理逻辑，并评估 `MAX_COORD` 常量所需的上限。
   - 若引入新的几何实体（例如自适应网格），应复用 `LocalCoord` 的 `rotated` 与 `lattice_i` 扩展机制，避免额外状态分散。
4. **调试策略**：利用 `GeometryState::mark_as_lost()` 的错误报告能力，可以在几何搜索失败时打印粒子 ID、层级与坐标，有助于复杂几何调试。【F:include/openmc/particle_data.h†L212-L239】【F:src/particle_data.cpp†L18-L87】

> **课程总结**：理解 `GeometryState`/`LocalCoord` 的栈式管理与跨模块接口，是掌握 OpenMC 几何与物理内核协同的首要步骤。后续章节将沿此数据流继续剖析 `ParticleData`、`Particle` 以及材料与反应抽样的关系。
