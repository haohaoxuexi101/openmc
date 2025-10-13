# 章节 5：ParticleData 内存布局与多物理耦合（深度课程）

本章节围绕 `ParticleData` 的成员体系、缓存策略、随机抽样与统计计数的接口设计展开，重点说明与 `GeometryState`、`Material`、`Surface`、`Cell`、`Universe`、`Lattice` 等模块的耦合路径，并给出高性能实现要点。章节划分为七个独立单元。

---

## 单元 5.1：类继承关系与角色定位

1. **继承链**：`Particle` ← `ParticleData` ← `GeometryState`。其中 `ParticleData` 负责物理属性、统计缓存与随机数流，而 `GeometryState` 管理几何栈。【F:include/openmc/particle.h†L23-L122】【F:include/openmc/particle_data.h†L208-L672】
2. **组合对象**：
   - 通过 `MacroXS`、`NuclideMicroXS`、`ElementMicroXS` 缓存材料依赖的截面数据；
   - 通过 `SourceSite`、`NuBank` 等结构承载二次粒子与裂变银行信息。
3. **OOP 策略**：保持 POD 兼容，便于 GPU 分支采用结构数组布局，同时统一接口访问函数（`E()`、`r()`、`macro_xs()` 等）。【F:include/openmc/particle_data.h†L400-L672】

---

## 单元 5.2：能量、群组与截面缓存

| 成员/方法 | 功能 | 耦合路径 |
|-----------|------|----------|
| `double& E()` / `double& E_last()` | 当前与上一次碰撞能量 | 供 `physics::collision`、`Material::calculate_xs` 决定散射或裂变产物分布。【F:include/openmc/particle_data.h†L519-L552】【F:include/openmc/material.h†L38-L200】 |
| `int& g()` / `int& g_last()` | 多群索引 | 与多群截面 `data::mg` 交互，支撑 GPU 事件模式。【F:include/openmc/particle_data.h†L519-L535】【F:src/particle.cpp†L200-L212】 |
| `NuclideMicroXS` & `ElementMicroXS` | 微观截面缓存 | 在 `Particle::event_calculate_xs` 中刷新，避免重复插值。【F:include/openmc/particle_data.h†L498-L513】【F:src/particle.cpp†L200-L212】 |
| `MacroXS macro_xs_` | 宏观截面 | 与 `Material::calculate_xs`、`Particle::event_advance` 联动决定碰撞距离与反应抽样。【F:include/openmc/particle_data.h†L508-L510】【F:src/particle.cpp†L213-L270】 |
| `CacheDataMG mg_xs_cache_` | 多群缓存 | 存储温度、角度索引，减少重复求解。【F:include/openmc/particle_data.h†L511-L514】 |

*高性能提示*：截面缓存放置在 `ParticleData` 内部，保证与粒子数据同处缓存线，减少访存延迟；更新函数统一通过引用返回，避免复制。

---

## 单元 5.3：随机数流与抽样机制

1. **随机种子**：`uint64_t seeds_[N_STREAMS]` 与 `int& stream()` 允许根据事件类型切换随机流，确保碰撞、几何、次级采样互不干扰。【F:include/openmc/particle_data.h†L456-L647】
2. **抽样接口**：
   - `double& mu_`、`double& time()` 等字段在碰撞抽样后写入结果，为脉冲高度、光子次级产额提供输入；
   - `uint64_t* current_seed()` 提供原始种子指针，供 `physics::sample_*` 系列函数使用，确保 Monte Carlo 可再现性。【F:include/openmc/particle_data.h†L515-L647】【F:include/openmc/physics.h†L24-L104】
3. **距离与事件选择**：`Particle::event_advance` 通过 `macro_xs().total` 与 `prn(current_seed())` 抽样自由程，比较 `GeometryState::boundary().distance` 判定下一事件是碰撞还是跨越表面。【F:src/particle.cpp†L213-L270】
4. **反应抽样**：`physics::collision` 读取 `ParticleData` 的能量、材料索引，分别调用 `sample_nuclide`、`sample_fission`、`sample_secondary_photons` 等函数，返回的次级粒子写入 `secondary_bank_`。【F:include/openmc/physics.h†L24-L104】【F:src/physics.cpp†L91-L1163】

---

## 单元 5.4：统计计数与 Tallies

1. **内部累计量**：`keff_tally_*` 字段累积轨迹、碰撞、吸收、泄漏贡献，`Particle::event_death` 将其通过原子操作汇总到全局变量，并在死亡前清零。【F:include/openmc/particle_data.h†L619-L644】【F:src/particle.cpp†L455-L492】
2. **碰撞/事件计数**：`n_collision()`、`n_event()` 记录历史次数，驱动分裂/轮盘算法与极端事件保护（如 `max_particle_events`）。【F:include/openmc/particle_data.h†L577-L637】【F:src/particle.cpp†L397-L405】
3. **Filter 匹配缓存**：`vector<FilterMatch> filter_matches_` 与 `flux_derivs_` 减少每次 tally 重复搜索，结合 `GeometryState::cell_last_` 形成几何-统计高速缓存。【F:include/openmc/particle_data.h†L463-L607】
4. **Pulse-height 支持**：`pht_storage_`、`tracks_` 与 `write_track_` 配合 `Particle::pht_*` 函数，记录能量沉积与轨迹。穿越表面时 `GeometryState` 的层级数据用于识别所属 cell。【F:include/openmc/particle_data.h†L467-L618】【F:src/particle.cpp†L494-L520】

---

## 单元 5.5：粒子银行与裂变后代管理

1. **Secondary Bank**：`vector<SourceSite> secondary_bank_` 存储次级粒子，`Particle::event_revive_from_secondary` 将银行尾部取出重新初始化粒子状态，并借助 `GeometryState` 重新定位出生单元。【F:include/openmc/particle_data.h†L459-L616】【F:src/particle.cpp†L397-L452】
2. **Fission Bank**：`NuBank` 队列和 `n_delayed_bank_` 数组用于延迟中子统计算法，与 `create_fission_sites`、`sample_fission_neutron` 耦合，实现平均与抽样两阶段处理。【F:include/openmc/particle_data.h†L444-L575】【F:src/physics.cpp†L91-L1157】
3. **Progeny 统计**：`int64_t& n_progeny()` 在粒子死亡时写入 `simulation::progeny_per_particle`，提供后续银行排序与并行归约依据。【F:include/openmc/particle_data.h†L638-L644】【F:src/particle.cpp†L486-L492】

---

## 单元 5.6：与几何实体的耦合与一致性维护

1. **Cell/Material 映射**：`cell_born_`、`material()`、`material_last()` 等成员保持粒子所处材料一致性，当 `Cell` 变换填充（Material/Universe/Lattice）时，`GeometryState` 更新 `coord_`，`ParticleData` 则更新材料与温度缓存。【F:include/openmc/particle_data.h†L336-L580】【F:include/openmc/cell.h†L150-L320】
2. **Surface 边界条件**：`surface()` 与 `boundary()` 中的 `coord_level`、`lattice_translation` 为 `Particle::cross_surface`、`Surface::reflect` 提供上下文，确保脉冲高度与表面 tally 识别正确的几何层级。【F:include/openmc/particle_data.h†L316-L333】【F:include/openmc/surface.h†L36-L119】【F:src/particle.cpp†L272-L314】
3. **Lattice 分布式实例**：`cell_instance()` 和 `lattice_i` 索引结合 `Lattice::offset`，用于 distribcell tallies，在晶格复制几何中准确区分每个 cell 实例。【F:include/openmc/particle_data.h†L266-L368】【F:include/openmc/lattice.h†L68-L138】

---

## 单元 5.7：性能优化与扩展建议

1. **内存局部性**：将常用标量（能量、权重、时间）紧凑存放在类前部，配合连续存储的向量容器（截面缓存、银行），最大化顺序访问效率。
2. **懒初始化**：`ParticleData()` 构造函数延迟分配银行与轨迹存储，只有在 tallies 或次级粒子功能启用时才扩大容量，以减少空跑成本。【F:include/openmc/particle_data.h†L490-L672】
3. **扩展流程**：添加新的统计量时，应同时更新 `event_death()` 归约、`zero_*` 重置函数，并评估 `Particle::event_revive_from_secondary` 是否需要处理相应缓存。
4. **多线程安全**：类本身不含锁，所有并行同步通过外部 OpenMP 原子或 MPI 归约完成，开发扩展时应避免在 `ParticleData` 内部引入共享状态的全局静态变量。

> **课程总结**：`ParticleData` 将几何状态、物理缓存、随机抽样与统计计数统一封装，为 `Particle` 的事件驱动流程提供基础。理解每个成员的耦合路径，是扩展物理模型或优化性能的关键前提。
