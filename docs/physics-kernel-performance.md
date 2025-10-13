# OpenMC 物理内核与高性能函数设计说明书

## 1. 目标与范围

本说明书聚焦于 OpenMC C++ 实现中的物理内核与性能关键路径，系统阐述数据布局、算法、并行化、缓存策略以及与外部库的协同机制。目标读者为性能优化工程师、物理模型开发者及需要在 OpenMC 上构建新型反应堆分析功能的研究者。

## 2. 计算热点全景

基于粒子输运模拟流程，计算热点主要集中在：

1. **截面评估**：包括核素微观截面插值、概率表采样与宏观截面加权，涉及 `Nuclide`, `Material`, `PhysicsCommon` 模块。【F:include/openmc/nuclide.h†L40-L322】【F:include/openmc/material.h†L42-L265】【F:src/physics_common.cpp†L35-L287】
2. **粒子推进与碰撞**：`Particle::event_advance`、`Particle::event_collide` 驱动主循环，对几何定位、随机采样和能量角分布处理提出性能要求。【F:include/openmc/particle.h†L52-L87】【F:src/particle.cpp†L38-L412】
3. **次级粒子管理与银行操作**：`Bank` 相关函数在大规模模拟中影响负载平衡和通信成本。【F:include/openmc/bank.h†L31-L205】【F:src/bank.cpp†L34-L285】
4. **Tally 累积与并行归并**：大量原子操作与归并需要内存一致性策略支持。【F:include/openmc/tallies/tally.h†L40-L266】【F:src/tallies/tally.cpp†L38-L379】

## 3. 数据布局与内存优化

### 3.1 粒子数据的结构体数组化（SoA）策略

- `ParticleData` 将粒子态拆分为基础字段与缓存结构，避免跨模块的虚函数开销；与 `Particle` 的轻量继承结合，实现面向对象接口与高性能数据布局的折中。【F:include/openmc/particle_data.h†L20-L257】【F:include/openmc/particle.h†L24-L118】
- `LocalCoord` 与截面缓存 (`NuclideMicroXS`, `MacroXS`) 均嵌入 `ParticleData`，减少指针跳转，提高缓存命中率。【F:include/openmc/particle_data.h†L80-L257】

### 3.2 交叉截面缓存与懒计算

- `NuclideMicroXS::last_E`、`MacroXS` 等字段记录最近一次计算状态，在能量或温度未发生变化时直接复用。【F:include/openmc/particle_data.h†L167-L257】
- `Material` 对宏观截面采用懒计算与按需缓存：粒子进入新材料时才触发求值，避免重复工作。【F:include/openmc/material.h†L170-L265】

### 3.3 能量网格与搜索优化

- `Nuclide` 内部采用对数并集能量网格（log union grid）与二分/指数搜索相结合；`physics_common.cpp` 提供 `calculate_micro_xs`、`calculate_macro_xs` 等函数在粒子事件中调用。【F:include/openmc/nuclide.h†L180-L322】【F:src/physics_common.cpp†L78-L287】
- 对热化截面使用 S(α,β) 表、NCrystal、`urr` 概率表等专用结构，以减少内核计算量。【F:include/openmc/thermal.h†L32-L265】【F:include/openmc/ncrystal_interface.h†L28-L165】【F:src/urr.cpp†L40-L311】

## 4. 物理事件流程与并行设计

### 4.1 事件驱动主循环

`simulation.cpp` 通过 `while` 循环驱动粒子事件：先调用 `event_calculate_xs` 更新截面，随后 `event_advance` 寻找最近边界或碰撞，`event_collide` 调用 `physics::collision` 处理反应，直至粒子终止。【F:src/simulation.cpp†L42-L289】

事件函数遵循以下性能策略：

- **最小化条件分支**：`Particle` 中的事件函数在热程、宏程中区分工作流程，减少分支失衡。
- **迭代局部性**：在一次事件中尽量复用已计算的截面、几何信息。

### 4.2 碰撞处理与反应采样

- `physics::collision` 根据粒子类型选择 `sample_neutron_reaction`、`sample_photon_reaction` 等函数，每个函数内部进一步调用更细粒度的处理器（如 `elastic_scatter`、`sab_scatter`）。【F:include/openmc/physics.h†L30-L98】【F:src/physics.cpp†L40-L415】
- 通过函数内联和模板化分布抽样器（`distribution_angle.h`、`distribution_energy.h`）减少虚调用开销。【F:include/openmc/distribution_angle.h†L32-L301】【F:include/openmc/distribution_energy.h†L33-L342】

### 4.3 次级粒子生成与银行管理

- `create_fission_sites`、`Particle::create_secondary` 在内核中直接写入 `Bank`，避免多余拷贝；`MessagePassing` 在批处理结束时聚合银行数据，并按照 MPI 排列归并，减少通信热点。【F:include/openmc/physics.h†L60-L94】【F:include/openmc/bank.h†L31-L205】【F:src/message_passing.cpp†L36-L342】
- 采用批次式处理与粒子分组，便于在 GPU/多线程场景中并行扩展。

## 5. 随机数、取样与矢量化

### 5.1 随机数生成

- `RandomLCG` 实现跳跃式种子推进，确保批次间独立性和并行可重复性。【F:include/openmc/random_lcg.h†L27-L166】
- `Particle` 将随机种子存储在粒子本地，允许线程独立推进，无需同步。

### 5.2 取样内核

- 能量、角度取样通过专用分布类实现，如 `TabularDistribution`、`DiscreteDistribution`，利用预计算累积概率表加速。【F:include/openmc/distribution_energy.h†L33-L342】
- `thermal.cpp`、`secondary_*` 系列函数依据粒子类型和能区调用不同算法，部分函数针对常见场景提供内联路径（例如自由气体模型）。【F:src/thermal.cpp†L35-L327】【F:src/secondary_thermal.cpp†L40-L288】

### 5.3 矢量化与编译器指令

- 数据结构选用对齐的 `vector`/`span` 实现，`openmc/vector.h` 封装了简化的 SIMD 友好接口。【F:include/openmc/vector.h†L24-L196】
- 在热点函数中通过 `#pragma omp simd` 或编译器自动向量化实现批量处理，代码中避免不必要的虚函数与指针间接。

## 6. 并行策略

### 6.1 线程并行（OpenMP）

- 批次级循环在 `simulation.cpp` 中通过 OpenMP 并行，粒子彼此独立，避免共享写冲突。
- Tally 在内部使用线程局部缓存，结束后再归并，降低锁开销。【F:src/tallies/tally.cpp†L185-L379】

### 6.2 进程并行（MPI）

- `MessagePassing` 管理粒子重平衡、银行交换、统计量归并。采用非阻塞通信与分级归并（根进程收集后再广播），确保可扩展性。【F:include/openmc/message_passing.h†L27-L214】【F:src/message_passing.cpp†L36-L342】
- 状态点写出前通过 `StatePoint::write` 聚合全局数据，使用集体 IO 减少文件系统压力。【F:include/openmc/state_point.h†L36-L223】【F:src/state_point.cpp†L40-L412】

### 6.3 加速器与外部接口

虽然当前主干集中在 CPU，但数据布局和批量接口已为 GPU 等加速器预留空间：

- `Particle::event_...` 系列函数与 `physics` 内核保持纯函数式接口，便于移植到 CUDA/HIP kernel。
- 通过 `shared_array`、`span` 等轻量抽象封装，兼容统一内存或显存访问模式。【F:include/openmc/shared_array.h†L24-L142】【F:include/openmc/span.h†L24-L150】

## 7. 输入/输出与高性能

- 大型表格通过 HDF5 原生接口读取，`HDF5Interface` 支持批量读取和延迟加载，避免初始化阶段的随机 IO。【F:include/openmc/hdf5_interface.h†L30-L179】【F:src/hdf5_interface.cpp†L34-L389】
- `TrackOutput`、`ParticleRestart` 在写入时采用缓冲批量写，减少频繁文件操作。【F:include/openmc/track_output.h†L28-L198】【F:src/particle_restart.cpp†L35-L275】

## 8. 性能调优建议

1. **热点定位**：结合 `Timer` 模块与外部分析器（VTune、perf）定位耗时函数。`timer.cpp` 提供层级计时器接口，可在内核新增计时节点。【F:include/openmc/timer.h†L27-L138】【F:src/timer.cpp†L34-L215】
2. **缓存友好性**：扩展数据结构时保持紧凑布局，避免在粒子事件中频繁分配内存。新增字段应放置在粒子结构末尾以减少缓存失配。
3. **并行负载均衡**：调优银行重分配阈值和批次粒子数，保障线程/进程间的工作量均衡。
4. **数值稳定性**：在概率表和插值算法中使用双精度并避免重复差值，减少误差积累对模拟结果的影响。
5. **扩展接口**：新增物理过程时，优先复用已有分布类和缓存机制；如需新分布，应遵循 `distribution_*` 抽象以保持一致的调用方式。

## 9. 汇总

OpenMC 物理内核的设计遵循“数据局部性优先、事件驱动接口、批次化并行调度”的总体思想。通过结构体数组化的粒子数据布局、懒计算截面缓存、函数化的碰撞处理以及 MPI/OpenMP 混合并行框架，形成在大规模核工程模拟中表现优越的高性能实现。后续优化可围绕 GPU 迁移、矢量化增强和 I/O 流水线进一步深化。
