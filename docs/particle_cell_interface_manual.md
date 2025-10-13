# OpenMC 粒子与几何对象接口深度解析

> 本说明书聚焦 `include/openmc/particle_data.h` 与 `include/openmc/cell.h` 两个核心头文件，逐行拆解其类层级、组合关系与接口设计意图。
> 目标是帮助读者理解 OpenMC 如何通过继承、组合与策略模式，将粒子的状态管理与几何导航紧密耦合，同时保持可扩展性与并行友好性。
>
> 文档结构概览：
> 1. **总览**——理解 `ParticleData` 继承自 `GeometryState` 的动机与整体职责。
> 2. **几何接口**——解析 `Surface`/`Region`/`Cell` 的布尔树与距离计算设计。
> 3. **状态字段**——拆解几何坐标栈、粒子属性与核数据缓存的组合方式。
> 4. **协作流程**——串联定位、距离计算、碰撞处理到跨界更新的接口调用链。
> 5. **关联速查**——以表格总结各类/结构体之间的依赖关系与设计要点。
> 6. **总结建议**——提炼继承与组合的设计经验并指向进一步实践。
> 7. **课程对照**——列出配套教学代码，帮助读者在工程中复现上述模式。

## 1. `ParticleData` 总览

`ParticleData` 是 OpenMC 中每个传输粒子的“胖状态对象”。其类定义位于 `include/openmc/particle_data.h`，继承自 `GeometryState`，后者定义在 `include/openmc/geometry.h`，提供与几何定位相关的成员（坐标、方向、当前单元/宇宙/晶格索引等）。【F:include/openmc/particle_data.h†L416-L474】【F:include/openmc/geometry.h†L28-L125】

继承的设计关键点：

- **运行时连续性**：粒子追踪时频繁需要同时读取空间坐标与几何层级（如当前 `cell`、`surface`、`lattice`），继承让这些字段与粒子状态驻留在同一结构体中，避免在热循环中进行额外的指针跳转。
- **接口复用**：`GeometryState` 暴露了 `coord()`、`surface()`、`cell()` 等访问器；`ParticleData` 直接继承后即可在传输阶段复用这些接口，而无需二次封装。【F:include/openmc/geometry.h†L72-L108】
- **责任分割**：`GeometryState` 不关心能量、权重等物理属性；`ParticleData` 则通过组合其它结构（如 `Bank`、`SourceSite`）承载这些信息，实现关注点分离。

`ParticleData` 的构造函数在源文件 `src/particle_data.cpp` 中实现，初始化几何状态并设置核数据缓存、随机种子、统计量等成员。【F:src/particle_data.cpp†L24-L113】

## 2. 几何对象接口

`include/openmc/cell.h` 描述了几何组合的核心：`Surface` 定义在 `surface.h` 中，`Region` 与 `Cell` 的布尔组合接口集中在该头文件。【F:include/openmc/cell.h†L20-L208】

### 2.1 Region 树

- `Region` 是抽象基类，定义 `contains(Position r)`、`distance(...)`、`bounding_box()` 等纯虚函数，代表布尔集合。派生类如 `Halfspace`、`Intersection`、`Union` 分别对应基本半空间与布尔运算节点。【F:include/openmc/cell.h†L56-L146】
- 每个 `Cell` 在构造时解析 XML，生成 `Region` 后缀表达式并持有一个 `std::unique_ptr<Region>`，表示该几何体由布尔树描述；树节点存储在 `expression_` 数组中，通过栈式遍历执行短路判定。【F:include/openmc/cell.h†L61-L146】
- `Region::contains` 通过“简单单元/复杂单元”两个分支分别使用短路和逆波兰栈算法，保证既有逻辑正确性又能在常见场景下压缩算术量。【F:include/openmc/cell.h†L81-L145】

### 2.2 Cell 的接口层

`Cell` 提供如下关键成员函数：【F:include/openmc/cell.h†L150-L240】

- `contains(Position r, Direction u, int32_t on_surface)`：委托给 `Region` 布尔树判定点是否落在单元内，同时允许通过 `on_surface` 快速处理跨界退化情况。
- `distance(Position r, Direction u, int32_t on_surface, GeometryState* p)`：计算粒子沿运动方向离开单元的距离，并将 `GeometryState` 传入以记录潜在的嵌套层级或邻接信息。
- `temperature(instance)`、`set_temperature`、`set_rotation`：管理每个实例的物性与坐标变换，接口保持 `const`/非 `const` 对称，方便多线程读取。【F:include/openmc/cell.h†L213-L232】
- `fill()`、`material()` 族函数控制单元是材料、宇宙还是晶格填充，通过 `Fill` 枚举将三种情形统一在单一接口下。【F:include/openmc/cell.h†L29-L35】【F:include/openmc/cell.h†L340-L383】

接口设计亮点：

- **布尔树缓存**：`Cell` 保存 `Region` 表达式、预处理的邻接表，并在距离计算中重用 `Surface::distance()` 结果，减少重复解析。【F:include/openmc/cell.h†L374-L424】
- **延迟解析策略**：许多 setter（如 `set_fill`, `set_temperature`）在 XML 加载时完成配置，运行阶段保持只读，确保多线程读取安全并利于向量化。

## 3. 几何状态堆栈与粒子字段

`GeometryState` 在 `particle_data.h` 中定义，负责跟踪粒子穿越嵌套宇宙/晶格时的“坐标栈”。【F:include/openmc/particle_data.h†L208-L375】

### 3.1 LocalCoord 层级

- `LocalCoord` 保存单个层级的局部坐标（位置 `r`、方向 `u`、所属 `cell/universe/lattice` 及晶格索引），并提供 `reset()`、`rotate()` 等函数管理局部旋转与清空操作。【F:include/openmc/particle_data.h†L81-L95】
- `GeometryState::coord_` 是 `std::vector<LocalCoord>`，`n_coord_` 标识当前有效层数；`lowest_coord()` 始终指向最内层局部坐标，方便在材料/晶格切换时定位。【F:include/openmc/particle_data.h†L262-L307】
- `surface_`、`BoundaryInfo` 字段记录粒子下一次离开界面的距离、索引及跨晶格偏移，供 `distance_to_boundary()` 与 `cross_lattice()` 使用。【F:include/openmc/particle_data.h†L316-L375】【F:include/openmc/geometry.h†L68-L79】

### 3.2 运行时状态缓存

- `material_`、`material_last_`、`sqrtkT_` 等字段在几何和核数据之间搭桥，避免每次碰撞都查表；`clear()` 会重置所有层级并将材料索引设为 `C_NONE`，用于生成新粒子或处理几何丢失。【F:include/openmc/particle_data.h†L223-L374】
- `r_last_current_`、`r_last_`、`u_last_` 记录粒子上一次碰撞或反射事件的空间数据，方便 tally 和事件驱动模式回放。【F:include/openmc/particle_data.h†L286-L305】

### 3.3 粒子属性

`ParticleData` 在 `GeometryState` 基础上组合核数据缓存、物理量、统计与随机数接口，是连接几何、核数据、物理三大模块的“胖对象”。【F:include/openmc/particle_data.h†L416-L672】

1. **核数据缓存**：`neutron_xs_`、`photon_xs_`、`macro_xs_`、`mg_xs_cache_` 等结构使得多群与连续能量路径可以复用同一接口，通过 `invalidate_neutron_xs()` 强制刷新缓存。【F:include/openmc/particle_data.h†L421-L655】
2. **生命周期与统计**：`type_`、`E_`、`g_`、`wgt_`、`n_collision_`、`keff_tally_*` 等字段服务于物理调度与 tallies。【F:include/openmc/particle_data.h†L426-L624】
3. **随机与次级管理**：`seeds_`、`stream_`、`secondary_bank_`、`nu_bank_`、`n_delayed_bank_[]` 等成员控制随机数流和裂变二次粒子，确保与 MPI/OpenMP 协调一致。【F:include/openmc/particle_data.h†L456-L575】【F:include/openmc/particle_data.h†L589-L615】
4. **调试与性能**：`trace_`、`tracks_`、`current_work_` 字段仅在调试或事件驱动模式启用，默认不分配内存开销；接口采用 `inline` getter/setter 减少粒子循环成本。【F:include/openmc/particle_data.h†L467-L672】

`ParticleData` 还提供 `zero_delayed_bank()`、`zero_flux_derivs()` 等批量工具函数，确保批次间状态复位不遗漏任何字段。【F:include/openmc/particle_data.h†L656-L672】

## 4. 几何-粒子接口协作流程

粒子追踪时，几何与粒子状态通过以下步骤协作：

1. **初始化**：`GeometryState::init_from_r_u()` 设置全局/局部坐标、清空表面标识；`ParticleData` 构造函数进一步初始化能量、随机数种子与核数据缓存。【F:include/openmc/particle_data.h†L236-L253】【F:src/particle_data.cpp†L24-L113】
2. **定位**：`exhaustive_find_cell()` 或 `neighbor_list_find_cell()` 根据 `GeometryState` 的 `coord_` 栈递归搜索宇宙/晶格，必要时调用 `Cell::contains()` 判定包含关系。【F:include/openmc/geometry.h†L56-L79】【F:include/openmc/cell.h†L150-L208】
3. **距离计算**：`distance_to_boundary()` 调用当前 `Cell::distance()`，返回 `BoundaryInfo`，并更新 `GeometryState::surface()` 及 `boundary()`。【F:include/openmc/geometry.h†L68-L79】【F:include/openmc/particle_data.h†L327-L339】
4. **碰撞处理**：物理模块根据 `ParticleData::macro_xs()` 或 `neutron_xs()` 采样自由程与反应通道，碰撞后通过 `material_last_`、`r_last_` 等缓存进行统计或俄罗斯轮盘调整。【F:include/openmc/particle_data.h†L421-L575】
5. **跨界更新**：当 `BoundaryInfo::surface != SURFACE_NONE` 时，`cross_lattice()` 更新晶格索引并调用 `GeometryState::move_distance()` 推进粒子位置，同时刷新 `coord_` 栈中的嵌套层级。【F:include/openmc/geometry.h†L68-L79】【F:include/openmc/particle_data.h†L236-L307】

该流程通过“继承 + 组合”的接口层搭建出松耦合协作：几何模块只依赖 `GeometryState`，物理模块只感知 `ParticleData` 的 getter/setter，而粒子对象则以内联函数连接两端。

## 5. 类与结构体关联速查表

| 模块 | 核心类型 | 关联字段/接口 | 设计要点 |
| ---- | -------- | ------------- | -------- |
| 几何状态 | `GeometryState`、`LocalCoord`、`BoundaryInfo` | `coord()`, `surface()`, `boundary()`, `move_distance()` | 坐标栈描述宇宙/晶格嵌套，`BoundaryInfo` 统一传递界面事件。【F:include/openmc/particle_data.h†L81-L375】【F:include/openmc/geometry.h†L68-L79】 |
| 几何体 | `Region`、`Cell` | `contains()`, `distance()`, `bounding_box()` | 布尔树 + 虚函数提供统一判定接口，支持材料/宇宙/晶格填充。【F:include/openmc/cell.h†L56-L383】 |
| 粒子状态 | `ParticleData`、`SourceSite`、`NuBank` | `neutron_xs()`, `E()`, `secondary_bank()` | 继承几何状态，组合核数据缓存与统计量，通过内联 getter 保持性能。【F:include/openmc/particle_data.h†L33-L672】 |
| 物理调度 | `Particle`（实现于 `particle.cpp`） | `event_calculate_xs()`, `create_secondary()` | 物理流程调用 `ParticleData` 接口，并在必要时更新几何状态以保持一致性。【F:src/particle.cpp†L1-L126】 |

## 6. 文档级总结与实践建议

- **继承封装共享几何信息**：将所有空间/几何字段集中在 `GeometryState`，让 `ParticleData`、`Ray` 等派生类共享一致接口，避免跨模块同步成本。【F:include/openmc/particle_data.h†L208-L375】【F:include/openmc/plot.h†L139-L223】
- **组合缓存加速物理核**：核数据缓存与统计量以 `vector`/结构体组合方式嵌入 `ParticleData`，可在不同物理模式间复用，且通过内联访问避免虚函数开销。【F:include/openmc/particle_data.h†L421-L672】
- **接口面向流程设计**：几何接口负责 `contains`/`distance`，几何状态提供 `move_distance`/`cross_lattice`，粒子接口暴露 `macro_xs`/`secondary_bank`。各模块通过明确职责层级实现高 cohesion、低 coupling。
- **教学示例联动**：面向对象课程第 6、7 课及随机射线系列示例将上述接口理念拆解为可编译程序，读者可逐步练习并回到正式源码进行对照巩固。【F:examples/courses/oop_design/lesson06_particle_data_bridge.cpp†L1-L120】【F:examples/random_ray_solver/lessons/lesson06_solver_driver.cpp†L1-L146】

通过本文档与示例的对照学习，可以全面掌握 `particle_data.h` 与 `cell.h` 的类/结构体关联、面向对象特性及接口设计理念，为自研求解器或扩展 OpenMC 奠定扎实基础。

## 7. 课程与示例

为帮助读者将上述理念应用到自定义求解器，本说明书配套以下课程示例：

- `examples/courses/oop_design/lesson01_geometry_state.cpp`：使用 `lesson::oop::geometry_state::LessonGeometryState` 精炼演示几何状态的职责划分。命名上的 `Lesson` 前缀提醒读者该类型仅为教学版实现，与 `include/openmc/particle_data.h` 中的正式 `GeometryState` 不会冲突。
- `examples/courses/oop_design/lesson06_particle_data_bridge.cpp`：展示如何通过继承教学版 `GeometryStateBridge` 构建轻量版 `ParticleDataLite`，并说明哪些字段适合放在基类/派生类中。
- `examples/courses/oop_design/lesson07_cell_interface_walkthrough.cpp`：拆解 `Cell` 的布尔树与距离计算接口，演示如何以组合方式拼装 Region 节点。
- `examples/random_ray_solver/lessons/lesson04_collision_physics.cpp` 与 `lesson05_bank_management.cpp`：说明随机射线求解器在碰撞抽样、裂变银行中的状态管理如何映射到 `ParticleData`/`Cell` 的接口。

配合本文档阅读实际头文件与课程代码，可以逐步理解 OpenMC 的对象模型，并在自研代码中复用其设计模式。

