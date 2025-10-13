# OpenMC 粒子与几何对象接口深度解析

> 本说明书聚焦 `include/openmc/particle_data.h` 与 `include/openmc/cell.h` 两个核心头文件，逐行拆解其类层级、组合关系与接口设计意图。
> 目标是帮助读者理解 OpenMC 如何通过继承、组合与策略模式，将粒子的状态管理与几何导航紧密耦合，同时保持可扩展性与并行友好性。
>
> 文档结构分为四部分：
> 1. 总览：理解 `ParticleData` 继承自 `GeometryState` 的动机与整体职责。
> 2. 几何层：解析 `Surface`/`Region`/`Cell` 的接口、布尔运算及缓存策略。
> 3. 粒子层：说明 `ParticleData` 如何组合核数据、随机数种子、碰撞历史等信息，并通过接口暴露给物理模块。
> 4. 设计实例：提供三个精炼的 C++ 课程文件，演示如何在自定义求解器中复用这些接口理念。

## 1. `ParticleData` 总览

`ParticleData` 是 OpenMC 中每个传输粒子的“胖状态对象”。其类定义位于 `include/openmc/particle_data.h`，继承自 `GeometryState`，后者定义在 `include/openmc/geometry.h`，提供与几何定位相关的成员（坐标、方向、当前单元/宇宙/晶格索引等）。【F:include/openmc/particle_data.h†L416-L474】【F:include/openmc/geometry.h†L28-L125】

继承的设计关键点：

- **运行时连续性**：粒子追踪时频繁需要同时读取空间坐标与几何层级（如当前 `cell`、`surface`、`lattice`），继承让这些字段与粒子状态驻留在同一结构体中，避免在热循环中进行额外的指针跳转。
- **接口复用**：`GeometryState` 暴露了 `coord()`、`surface()`、`cell()` 等访问器；`ParticleData` 直接继承后即可在传输阶段复用这些接口，而无需二次封装。【F:include/openmc/geometry.h†L72-L108】
- **责任分割**：`GeometryState` 不关心能量、权重等物理属性；`ParticleData` 则通过组合其它结构（如 `Bank`、`SourceSite`）承载这些信息，实现关注点分离。

`ParticleData` 的构造函数在源文件 `src/particle_data.cpp` 中实现，初始化几何状态并设置核数据缓存、随机种子、统计量等成员。【F:src/particle_data.cpp†L24-L113】

## 2. 几何对象接口

`include/openmc/cell.h` 描述了几何组合的核心：`Surface` 定义在 `surface.h` 中，`Region` 与 `Cell` 的布尔组合接口集中在该头文件。【F:include/openmc/cell.h†L20-L198】

### 2.1 Region 树

- `Region` 是抽象基类，定义 `contains(Position r)` 和 `bounding_box()` 等纯虚函数，代表布尔集合。派生类如 `Halfspace`、`Intersection`、`Union` 分别对应基本半空间与布尔运算节点。【F:include/openmc/cell.h†L38-L135】
- 每个 `Cell` 持有一个 `std::unique_ptr<Region>`，表示该几何体由布尔树描述。布尔树的存在允许在几何复杂时仍保持统一的接口：只要实现 `contains` 与 `distance`，上层即可通过多态调用处理不同形状。

### 2.2 Cell 的接口层

`Cell` 提供如下关键成员函数：【F:include/openmc/cell.h†L142-L198】

- `contains(Position r)`：调用 `Region` 布尔树判定点是否落在单元内。
- `distance(const Particle& p)`：计算粒子沿运动方向离开单元的距离。该函数内部会遍历单元边界，调用每个 `Surface::distance()`，选择最近的候选，并处理周期性/格子重复等情况。
- `material_`、`temperature_` 等属性：指向材料表、温度因子；通过接口 `material()`、`temperature()` 暴露给粒子物理模块。
- `fill_`：用于 Universe 嵌套，当 `cell` 被其它宇宙填充时，`fill()` 接口返回嵌套的几何索引。

接口设计亮点：

- **布尔树缓存**：`Cell` 持有 `Universe` 层级索引，使得粒子穿越不同 Universe 时可快速切换布尔树指针，避免重复解析 XML。
- **延迟解析策略**：许多 setter（如 `set_fill`, `set_temperature`）在 XML 加载时填充，运行阶段则保持只读，确保多线程读取安全。

## 3. 粒子状态字段

`ParticleData` 将粒子运行时需要的信息分成四类字段：【F:include/openmc/particle_data.h†L475-L655】

1. **标识与生命周期**：包括 `id_`、`generation_`、`alive_` 等，用于控制主循环、输出历史。
2. **物理量**：`E_`（能量）、`wgt_`（权重）、`mu_`/`phi_`（角度）、`uvw_`（方向向量），以及 `cell_last_`、`material_last_` 等缓存，可减少重复查询。
3. **统计累积**：`path_length_`、`num_collisions_`、`rxn_` 等字段为 tally 提供即时数据。
4. **随机与次级粒子管理**：`uint64_t seeds_[4]` 储存多个独立随机流；`std::vector<Bank>` 作为裂变银行缓存；`trace_` 记录调试用的空间轨迹。

接口层面，`ParticleData` 提供大量内联 getter/setter，例如 `energy()`, `set_energy(double)`, `neutron_xs()` 等。这些函数大多定义为 `inline` 以避免函数调用开销，同时确保粒子传输内循环中仍能保持高性能。【F:include/openmc/particle_data.h†L520-L655】

此外，`ParticleData` 暴露 `create_secondary(ParticleBank&)`、`from_source(const SourceSite&, uint64_t id)` 等成员，用于与源项、次级粒子数据结构交互，体现组合模式的应用。【F:include/openmc/particle_data.h†L656-L795】

## 4. 课程与示例

为帮助读者将上述理念应用到自定义求解器，本说明书配套以下课程示例：

- `examples/courses/oop_design/lesson01_geometry_state.cpp`：使用 `lesson::oop::geometry_state::LessonGeometryState` 精炼演示几何状态的职责划分。命名上的 `Lesson` 前缀提醒读者该类型仅为教学版实现，与 `include/openmc/particle_data.h` 中的正式 `GeometryState` 不会冲突。
- `examples/courses/oop_design/lesson06_particle_data_bridge.cpp`：展示如何通过继承教学版 `GeometryStateBridge` 构建轻量版 `ParticleDataLite`，并说明哪些字段适合放在基类/派生类中。
- `examples/courses/oop_design/lesson07_cell_interface_walkthrough.cpp`：拆解 `Cell` 的布尔树与距离计算接口，演示如何以组合方式拼装 Region 节点。
- `examples/random_ray_solver/lessons/lesson04_collision_physics.cpp` 与 `lesson05_bank_management.cpp`：说明随机射线求解器在碰撞抽样、裂变银行中的状态管理如何映射到 `ParticleData`/`Cell` 的接口。

配合本文档阅读实际头文件与课程代码，可以逐步理解 OpenMC 的对象模型，并在自研代码中复用其设计模式。

