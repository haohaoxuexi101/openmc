# OpenMC C++ 代码设计手册

> 本手册旨在帮助读者深入理解 OpenMC 的 C++ 核心实现。内容覆盖工程目录、类与模块职责、运行时数据流、并行策略、输入输出机制以及常用编程技巧。所有说明均基于 `src/` 与 `include/openmc/` 目录下的实际源码，并辅以与 Python 接口、示例程序的关联说明。

## 1. 总览与设计哲学

### 1.1 设计目标

OpenMC 的 C++ 层面追求以下目标：

- **物理正确性优先**：确保连续能量、多群、中子与光子等复杂物理过程的准确模拟。
- **可扩展性与可维护性**：通过模块化命名空间、头源文件分离、清晰的数据结构支持新物理模型的扩展。
- **跨语言互操作**：核心 C++ 负责重型计算，Python 接口通过 C API (`openmc/capi.h`) 驱动建模与后处理。
- **高性能并行**：广泛使用 OpenMP、MPI、向量化以及缓存友好的内存布局支持大规模计算。

### 1.2 核心命名空间与全局结构

- `namespace openmc`：所有公共类、函数、全局状态的顶级命名空间。
- `namespace model`：在多个源文件中共享的只读/可变全局模型数据（如材料、几何、计数器）。例如 `model::materials`、`model::universes` 均在对应模块中定义并在全局使用。
- `namespace data` 与 `namespace simulation`：分别负责核数据缓存与运行时统计（见 `src/particle_data.cpp`、`src/simulation.cpp`）。

这种分区使得不同功能块之间可以通过轻量级全局引用共享状态，同时依靠头文件声明确保编译期一致性。

### 1.3 构建系统

- 使用顶层 `CMakeLists.txt` 驱动构建，细节可参见 `cmake/` 子目录配置。
- 外部依赖（`fmt`, `xtensor`, `pugixml`, `HDF5` 等）通过 `vendor/` 或 `cmake/Find*.cmake` 管理。
- `include/openmc/` 暴露公共 API，`src/` 提供实现；例如 `include/openmc/material.h` 定义 `Material` 类接口，对应实现位于 `src/material.cpp`。

## 2. 工程目录与模块职责

### 2.1 顶层目录

- `src/`：C++ 核心实现，后续章节逐模块解读。
- `include/openmc/`：与 `src/` 一一对应的头文件，提供类声明、常量、模板工具等。
- `docs/`：开发者文档、理论手册（本文件即新增于此）。
- `examples/`：示例与回归测试输入。`examples/simple_transport/main.cpp` 是仿 OpenMC 面向对象风格的教学代码。
  自 3D 版本起，该目录额外提供 `lessons/` 子目录，将 Vec3、Material、几何、粒子银行、PhysicsDriver、随机引擎、晶格定位、tally 管理与 Simulation 顶层控制等特性拆解为可独立编译的小程序，方便逐模块对照源码学习。【F:examples/simple_transport/README.md†L27-L58】

  `lessons/` 当前包含以下教学代码，每个都附带中文注释、控制台演示：

  | 文件 | 目的 |
  | ---- | ---- |
  | `lesson_vec3.cpp` | 掌握三维向量工具类的值语义与运算封装。 |
  | `lesson_material.cpp` | 理解多群截面数据结构与碰撞分支概率。 |
  | `lesson_oop_architecture.cpp` | 聚焦 Surface/Cell/Geometry/Particle 的接口划分与组合实践。 |
  | `lesson_geometry.cpp` | 分析平面、单元、宇宙的组合关系与距离算法。 |
  | `lesson_particle_bank.cpp` | 学习粒子状态与 source/fission bank 的管理接口。 |
  | `lesson_physics.cpp` | 追踪粒子自由程、边界处理与散射/裂变抽样。 |
  | `lesson_random.cpp` | 解析随机引擎封装、均匀与各向同性抽样、离散采样接口。 |
  | `lesson_lattice.cpp` | 演示 Universe + Lattice3D + Geometry 的三维定位流程。 |
  | `lesson_tally.cpp` | 拆解 track-length tally 如何按 (cell, group) 累积统计量。 |
  | `lesson_simulation.cpp` | 展示 Simulation 设置、物理驱动、重采样与统计管理的顶层协作。 |
  | `lesson_mgxs_builder.cpp` | 复刻多群截面折合、材料混配与 Legendre 系数管理流程。 |
  | `lesson_scatter_anisotropic.cpp` | 对照多群与连续能量各向异性散射的抽样接口。 |
  | `lesson_distribution_sampling.cpp` | 系统说明分布基类、分段线性抽样与 Alias 表。 |
  | `lesson_neighbor_table.cpp` | 展示笛卡尔晶格邻接表构建及边界条件挂接点。 |
  | `lesson_doppler_broadening.cpp` | 解释多普勒展宽的高斯卷积与温度依赖实现。 |
  | `lesson_probability_table.cpp` | 演示未解析共振区概率表的采样逻辑。 |
  | `lesson_transport_algorithm.cpp` | 复刻多群输运主循环、裂变银行与代际重采样流程。 |
  | `lesson_parallel_driver.cpp` | 演示线程粒子分片、局部随机流与结果归并的并行驱动。 |
  > **提示**：完整的三维随机射线多群求解器已独立为 `examples/random_ray_solver/main.cpp`，可用于综合实践材料/几何/粒子/统计的协同设计。

### 2.2 编译入口

- `src/main.cpp` 是 C++ 可执行文件的入口，实现命令行解析、MPI 初始化、调用 `openmc_run()` 驱动主循环。
- `src/initialize.cpp` 与 `src/finalize.cpp` 分别管理运行前后资源准备与清理。

## 3. 公用基础设施

### 3.1 错误与日志

- `src/error.cpp`：实现 `fatal_error`, `warning`, `write_message` 等，依靠 `fmt` 格式化并在 MPI 环境中收敛信息。
- `include/openmc/error.h`：声明异常类与错误接口，广泛被各模块 `#include`，如 `src/material.cpp` 在解析 XML 失败时调用 `fatal_error()`。【F:src/material.cpp†L41-L66】

### 3.2 实用工具

- `src/file_utils.cpp` 与 `include/openmc/file_utils.h`：封装文件系统操作，如路径组合、存在性检查，被核数据加载与输出模块使用。
- `src/string_utils.cpp`：提供字符串拆分、大小写转换等辅助函数。
- `src/math_functions.cpp`：实现向量数学、矩阵运算与特殊函数，供几何与物理模块调用。
- `src/timer.cpp`：基于 RAII 的计时器，`simulation::time_total` 等计时变量在 `src/simulation.cpp` 中使用以统计运行时间。【F:src/simulation.cpp†L1-L64】

### 3.3 随机数设施

- `src/random_lcg.cpp`：实现可重现的线性同余随机数生成器及流管理，粒子对象通过 `Particle::seeds` 维护多个流。
- `src/random_dist.cpp`：提供通用随机采样工具（如多项式、指数分布采样），在物理模块内被调用。
- `lessons/lesson_distribution_sampling.cpp` 展示了如何将 Piecewise Linear、Alias 表与复合分布组合，模拟 `openmc/random_dist.h` 中多态分布接口的设计。代码通过基类 `Distribution` 与派生实现说明“接口 + 策略”模式。【F:examples/simple_transport/lessons/lesson_distribution_sampling.cpp†L1-L173】

## 4. 核数据与材料层

### 4.1 核数据加载

- `src/cross_sections.cpp`：解析 `cross_sections.xml`，构建 HDF5 数据路径并缓存到 `data::libraries`。
- `src/xsdata.cpp` 与 `src/nuclide.cpp`：加载连续能量核数据表、预处理散射核、构建温度依赖插值。
- `src/mgxs.cpp`、`src/mgxs_interface.cpp`：实现多群截面的稠密/稀疏存储、插值接口。
- `src/thermal.cpp`：处理热中子散射数据 S(α,β)。

### 4.2 材料对象

- `include/openmc/material.h` 定义 `Material` 成员（密度、成分、温度等）。
- `src/material.cpp` 实现 XML 解析、密度单位转换、核素索引映射、与 NCrystal 的对接等逻辑。在构造函数中读取 `<density>`、`<nuclide>` 节点，将宏观截面转换为内部核素表示，并通过 `init_nuclide_index()` 建立材料内核素到全局核素数组的索引映射。【F:src/material.cpp†L33-L133】【F:src/simulation.cpp†L78-L95】
- `model::materials` 是 `std::vector<unique_ptr<Material>>`，供粒子追踪与截面查询使用。

### 4.3 反应与截面对象

- `src/reaction.cpp`, `src/reaction_product.cpp`：封装 ENDF 反应数据、二次粒子采样器。
- `src/scattdata.cpp`：散射角度、能量分布处理。
- `src/secondary_*` 系列实现不同类型的二次粒子生成（Kalbach-Mann、热散射等）。

#### 4.4 多群折合、概率表与多普勒展宽

- `lessons/lesson_mgxs_builder.cpp` 手把手演练如何从细能格核素截面折合出多群 Σ_t/Σ_s/νΣ_f，并与材料配比组合，对应 OpenMC `MGXS` 类的构造流程。【F:examples/simple_transport/lessons/lesson_mgxs_builder.cpp†L1-L207】
- `lessons/lesson_probability_table.cpp` 提炼 URR 概率表的采样逻辑，模拟 `openmc/probability_table.h` 中根据随机数选择截面条目的过程。【F:examples/simple_transport/lessons/lesson_probability_table.cpp†L1-L85】
- `lessons/lesson_doppler_broadening.cpp` 通过高斯卷积演示温度展宽，帮助理解 `openmc/endf.cpp` 中的多普勒处理接口及其对截面平滑的作用。【F:examples/simple_transport/lessons/lesson_doppler_broadening.cpp†L1-L102】

## 5. 几何建模层

### 5.1 空间坐标结构

- `include/openmc/geometry.h` 定义 `GeometryState`（保存坐标层级、当前单元、所在表面等）。
- `src/geometry.cpp` 提供几何搜索、邻接列表、重叠检测等核心算法。例如 `find_cell_inner()` 利用邻居列表或穷举遍历确定粒子所在 `Cell`，并递归进入嵌套 `Universe`。【F:src/geometry.cpp†L1-L118】
- `lessons/lesson_random_ray_solver.cpp` 展示了如何在三维盒体集合中通过逆序遍历确定区域优先级，复刻 OpenMC `Cell` 树的定位逻辑。【F:examples/simple_transport/lessons/lesson_random_ray_solver.cpp†L103-L156】
- `lessons/lesson_oop_architecture.cpp` 以接口层 (`Surface`)、组合层 (`Cell`/`Geometry`) 与状态层 (`Particle`) 分别实现类职责，演示继承与组合协同时如何保持线程安全的值语义。【F:examples/simple_transport/lessons/lesson_oop_architecture.cpp†L1-L260】

### 5.2 基本构件

- `src/surface.cpp`：定义 `Surface` 基类及平面、圆柱、球面等派生体，支持解析 XML 参数和计算表面距离。
- `src/cell.cpp`：将布尔几何、填充（材料/宇宙/晶格）绑定到面片集合，提供 `contains()` 判断与路径长度限制。
- `src/universe.cpp`：管理 `Universe` 与其包含的 `Cell` 集合。
- `src/lattice.cpp`：实现规则与三角晶格，支持多维索引转换、周期边界；几何搜索时通过 `lat.offset()` 更新分布式单元实例号。【F:src/geometry.cpp†L65-L112】
- `src/mesh.cpp`：用于空间网格划分，支持体积分布 tally、加速器权窗等。

### 5.3 几何预处理与加速

- `src/geometry_aux.cpp`：提供几何初始化、包围盒计算、几何体积估计等辅助工具。
- `src/dagmc.cpp`：与 DAGMC 几何引擎对接，实现 CAD 网格追踪。
- `src/random_ray/`：实现 Random Ray 求解器的几何剖分、矩阵装配。
- `lessons/lesson_neighbor_table.cpp` 给出了 3×2×2 晶格的邻接表构建示例，说明如何在初始化阶段预填充 `Cell::neighbors_`，减少追踪时的几何搜索成本。【F:examples/simple_transport/lessons/lesson_neighbor_table.cpp†L1-L94】

## 6. 粒子状态与生命周期

### 6.1 粒子类

- `include/openmc/particle.h` 定义 `Particle` 的状态属性（位置、方向、能量、权重、历史计数）。
- `src/particle.cpp` 实现粒子速度计算、从源点构造、二次粒子缓存 (`secondary_bank()`)、碰撞后状态更新等。`Particle::event_calculate_xs()` 会设置随机数流、缓存碰撞前属性、触发局部几何搜索与截面计算，为随后的物理处理做准备。【F:src/particle.cpp†L1-L126】

### 6.2 源项与粒子银行

- `src/source.cpp`：根据输入配置或外部文件采样初始粒子；支持固定源、裂变源、自定义分布。
- `src/bank.cpp`：实现 source/fission/secondary bank 的结构体与 I/O。
- `src/eigenvalue.cpp`：在 k-eigenvalue 模式下处理裂变银行的归一化、广义 Batcher 法。

### 6.3 事件与跟踪

- `src/event.cpp`：在事件驱动（event-based）模式下管理粒子队列，支持大规模并行。
- `src/track_output.cpp`：将粒子轨迹写入 VTK/VTU 文件以便后处理。

## 7. 物理过程实现

### 7.1 碰撞核心

- `src/physics.cpp`：实现碰撞调度与多粒子类型支持。`collision()` 根据粒子类型调用 `sample_neutron_reaction()`、`sample_photon_reaction()` 等，并在需要时创建裂变源、执行生存偏倚、触发俄罗斯轮盘等。【F:src/physics.cpp†L1-L101】
- `src/physics_common.cpp`：放置散射角采样、移动性修正等共享工具。
- `src/physics_mg.cpp`：多群模式下的截面插值与群间迁移。
- `src/photon.cpp`、`src/bremsstrahlung.cpp`：处理光子与次级电子物理。
- `lessons/lesson_scatter_anisotropic.cpp` 将 Legendre 展开、多群散射角拒绝采样与连续能量分段线性分布并列展示，帮助读者对照 `physics_common.cpp` 与 `angle_distribution.cpp` 的策略接口。【F:examples/simple_transport/lessons/lesson_scatter_anisotropic.cpp†L1-L154】
- `lessons/lesson_transport_algorithm.cpp` 以教学示例复刻多群自由程抽样、碰撞分支选择、裂变库累积与代际重采样流程，可与 `physics.cpp`/`simulation.cpp` 对照理解事件流向。【F:examples/simple_transport/lessons/lesson_transport_algorithm.cpp†L1-L350】

### 7.2 自由程采样与几何交互

- `src/particle.cpp` 与 `src/geometry.cpp` 协同：`Particle::event_advance()`（后续代码段）根据当前材料的总截面抽取碰撞距离，并调用几何模块定位下一个表面或碰撞点。
- `src/boundary_condition.cpp`：定义反射、真空、周期等边界条件行为。

### 7.3 统计与迭代流程

- `src/simulation.cpp`：驱动批次循环、世代循环、历史循环；调用 `openmc_next_batch`、`openmc_reset` 管理统计量、tally 与权重窗口。初始化阶段会构造银行、设置核素索引并在必要时恢复状态点。【F:src/simulation.cpp†L1-L115】
- `src/settings.cpp`：解析 `settings.xml`，将 XML 节点映射为 `settings::` 全局变量（如能量截断、并行模式、加速器开关）。

## 8. 计分系统（Tallies）

- `src/tallies/tally.cpp`：`Tally` 基类，维护滤波器组合、核函数、结果数组以及失配检查。
- `src/tallies/filter*.cpp`：每个文件实现一个 `Filter` 派生类，例如 `filter_cell.cpp` 根据单元 ID 过滤，`filter_energy.cpp` 处理能量区间。
- `src/tallies/tally_scoring.cpp`：调度计分流程，根据粒子事件调用相应 Filter 与 Score。
- `src/tallies/trigger.cpp`：实现基于不确定度的自动停止准则。
- Tallies 在 `simulation.cpp` 初始化阶段调用 `set_strides()`、`init_results()`，确保 OpenMP 线程安全。【F:src/simulation.cpp†L87-L101】
- `examples/random_ray_solver/main.cpp` 提供完整的三维随机射线多群求解器，涵盖材料/几何配置解析、粒子银行裂变再生以及曲线长度计分，为理解 OpenMC tally 架构与碰撞处理提供综合实战代码。【F:examples/random_ray_solver/main.cpp†L1-L786】

## 9. 输入、输出与状态点

- `src/xml_interface.cpp`：封装 PugiXML 操作，提供 `get_node_value`、`check_for_node` 等辅助方法。
- `src/settings.cpp`、`src/material.cpp`、`src/geometry.cpp` 等模块均借助上述接口解析 XML。
- `src/state_point.cpp`：负责状态点 HDF5 文件输出，包含 tallies、k-effective、几何体积等信息。
- `src/summary.cpp`：生成 `summary.h5`，记录输入模型概述供 Python API 快速加载。
- `src/output.cpp`：构造屏幕输出、写日志文件。

## 10. 并行与加速特性

- `src/message_passing.cpp`：MPI 包装与通信工具，提供 `bcast`, `reduce` 等函数。
- `src/weight_windows.cpp`：实现自适应权窗算法，配合 MPI 分布式收敛。
- `src/cmfd_solver.cpp`：关联协同多群扩散加速（CMFD），通过松弛迭代与 Monte Carlo 模块互换通量。
- `src/event.cpp` 与 `settings::event_based`：提供事件驱动模拟，允许大量在途粒子跨线程调度，提高大规模并行效率。
- `lessons/lesson_parallel_driver.cpp` 以线程安全的粒子分片、局部 RNG、结果归并示例展示如何实现“多线程 power iteration”，帮助读者在阅读 `simulation.cpp` 与 OpenMP/MPI 代码前先掌握并行化要点。【F:examples/simple_transport/lessons/lesson_parallel_driver.cpp†L1-L356】
- OpenMP：许多循环使用 `#pragma omp` 并行化，例如 `src/geometry.cpp` 的重叠检查累加使用 `#pragma omp atomic`。【F:src/geometry.cpp†L29-L57】

## 11. 扩展接口

- `src/dagmc.cpp`：与 CAD 中的 MOAB / DAGMC 集成，实现基于三角网格的几何追踪。
- `src/mcpl_interface.cpp`、`src/ncrystal_interface.cpp`：对接外部粒子文件（MCPL）与晶体散射库（NCrystal）。
- C API 函数集中在 `include/openmc/capi.h` 声明，并在 `src/simulation.cpp`、`src/state_point.cpp` 等文件内以 `extern "C"` 形式实现，供 Python API 调用。
- `src/mgxs_interface.cpp`：允许 Python 端加载用户自定义多群截面并注入 C++ 核心。

## 12. 示例程序的关联

- `examples/simple_transport/main.cpp` 提供教学版三维输运器，模仿 OpenMC 分层结构：
  - **材料类**：封装多群截面表，调用方式参考 `Material::init_nuclide_index()`。
  - **几何类**：构建 `Cell`、`Lattice` 树形结构，并在行走时递归查询邻居。
  - **粒子类**：维护位置、方向、能量群，与真实 `Particle` 接口类似。
  - **模拟驱动**：按批次、世代循环，并输出类似 `k-effective` 的统计。
- 该示例可直接编译运行，用于理解真实源码的职责分配与交互方式。

## 13. 编程技巧与风格

### 13.1 头源分离与前向声明

- 所有类在 `include/openmc/*.h` 中声明，源文件通过 `#include` 对应头文件实现成员函数。通过减少包含依赖降低编译时间。
- 在头文件中尽量使用前向声明以减少互相依赖，例如 `include/openmc/material.h` 中前向声明 `class Nuclide;`，源文件再包含 `openmc/nuclide.h`。

### 13.2 RAII 与资源管理

- 文件句柄、HDF5 对象、随机数流等均通过 RAII 包装自动释放。
- `unique_ptr` 广泛用于持有模型对象 (`model::materials`, `model::universes`)；引用时采用裸指针以避免额外开销。

### 13.3 数据局部性与缓存友好

- 粒子碰撞流程中尽量保持结构体连续存储（例如 `Particle::secondary_bank()` 使用 `std::vector<Bank>`）。
- 截面数据 (`XSData`) 与核素数据 (`Nuclide`) 提前展开至平坦数组，减少运行时查找。

### 13.4 并行安全

- 使用 OpenMP 原子操作或线程局部缓存避免竞争，如几何重叠计数的 `#pragma omp atomic`。
- MPI 相关函数统一封装在 `message_passing.cpp`，确保不同模块共享一致的通信语义。

### 13.5 错误防御

- 解析输入时大量使用 `fatal_error` 提供可读错误信息，保证用户模型错误不会导致未定义行为。
- 采用 `assert` 与条件检查确保边界情况（如 `Material` 中密度单位合法性）。

## 14. 继承、组合与接口设计详解

### 14.1 几何层的继承树

OpenMC 将几何抽象为“曲面—单元—宇宙—晶格”的层级：

- `Surface` 作为抽象基类，规定 `evaluate`、`distance`、`normal` 等纯虚接口，所有具体曲面（轴对齐平面、圆柱、球面、一般平面等）继承并覆写这些成员，以多态方式服务于几何搜索和边界条件。【F:include/openmc/surface.h†L36-L195】
- `Cell` 自身同样是抽象类，将布尔区域判断与距离计算接口化，具体的 `CSGCell`、`DAGCell` 等派生类负责实现。派生体通过组合 `Surface` 指针来表达布尔体积，再通过 `Fill` 枚举将材料、子宇宙或晶格装配起来，从而在“继承表达行为、组合承载结构”之间取得平衡。【F:include/openmc/cell.h†L29-L199】
- `GeometryState` 则是复合类型：它持有坐标栈、上一次碰撞位置、当前所在表面等多层数据，为粒子行走算法提供一个“组合所有需要信息的接口对象”。这也是为何几何模块大量传递 `GeometryState&` ——通过组合将粒子状态与几何工具解耦。【F:include/openmc/particle_data.h†L206-L318】

### 14.2 组合驱动的粒子—材料耦合

- `Material` 类通过组合 `vector<int>`、`xt::xtensor` 等容器保存密度、核素索引、热散射表，同时暴露 `calculate_xs(Particle&)` 等接口，从而在不暴露内部细节的情况下与粒子交互。粒子对象只需调用接口，即可在碰撞前得到当前材料的横截面缓存，这体现了“接口分离 + 组合存储”带来的高内聚。【F:include/openmc/material.h†L35-L195】
- 粒子追踪时，`Particle` 持有 `GeometryState`（位置/方向/坐标层级）并引用全局材料表；材料与几何并非通过继承关联，而是以组合方式共享 `model::materials`、`model::cells` 等向量，既方便在 MPI 进程间广播，也便于在运行过程中按索引访问。【F:include/openmc/particle_data.h†L206-L314】【F:include/openmc/material.h†L25-L75】

### 14.3 接口与工厂模式

- Tallies 子系统以 `Tally` 为核心，接口提供 `set_filters`、`set_scores` 等高层方法，而计分细节交由过滤器组合决定。`Tally::create` 充当工厂函数，确保在构造时自动注册到全局容器并分配 ID，使得调用者完全通过接口而非构造函数直接管理生命周期。【F:include/openmc/tallies/tally.h†L24-L196】
- `Filter` 同样以抽象基类形式存在，纯虚函数 `from_xml`、`get_all_bins`、`text_label` 等要求派生类完整定义“如何筛选事件”。其静态模板 `Filter::create` 与字符串工厂实现了统一构建入口，内部则通过组合 `model::tally_filters`（`vector<unique_ptr<Filter>>`）集中管理所有派生对象，实现接口层的透明扩展。【F:include/openmc/tallies/filter.h†L19-L179】
- 边界条件模块使用策略模式：`BoundaryCondition` 提供统一接口 `handle_particle`，具体的真空、镜面、白噪声、周期边界均以派生类实现，从而允许在解析 XML 时按类型动态构造合适的策略对象，运行期通过多态回调处理粒子状态。【F:include/openmc/boundary_condition.h†L20-L149】

### 14.4 教学案例：在示例代码中练习继承与组合

教学示例 `examples/simple_transport/main.cpp` 将上述理念以更易读的方式呈现：

1. **组合搭建几何骨架**：`Surface`、`Cell`、`Universe`、`Lattice3D` 与 `Geometry` 五个类互相组合，借助面指针、邻接表与晶格索引构建出 2×2×2 的体元装配；读者可尝试新增一个 `Surface` 并在 `Cell` 构造中传入，从而模拟“添加局部反射面”的扩展。【F:examples/simple_transport/main.cpp†L96-L280】
2. **策略化边界处理**：`Cell::BoundaryHit` 负责记录最邻近表面，`PhysicsDriver::handle_boundary` 再根据 `Surface::boundary()` 的枚举值执行真空吸收、镜面反射或跨界跳转，与 OpenMC 的策略模式对应。若想加入“白反射”，可参考 `Boundary::Reflective` 的分支扩展新的处理函数。【F:examples/simple_transport/main.cpp†L151-L479】
3. **接口化物理流程**：`PhysicsDriver::transport_generation`、`transport_particle`、`collide` 将源项初始化、自由程抽样、裂变银行管理解耦，模拟了真实代码中“接口驱动”的层次；通过继承或组合独立的散射核对象，即可进一步练习如何替换物理模型而不破坏外部接口。【F:examples/simple_transport/main.cpp†L371-L620】
4. **组合统计模块**：`TallyManager` 聚合二维向量存储 `(cell, group)` 路径长度，`Simulation` 则组合 `Geometry`、`PhysicsDriver` 与 `Bank` 对象完成批次驱动。尝试在 `TallyManager` 中添加新的成员（例如角向矩），即可体会到组合带来的扩展性。【F:examples/simple_transport/main.cpp†L326-L733】

通过将 OpenMC 源码与教学示例对照，可以系统体会“抽象接口表达策略、组合容器承载数据”的核心哲学，并在日常练习中不断巩固。

## 15. 附录：C++ 源码索引

下表对 `src/` 及子目录的 `.cpp` 文件进行逐一说明，便于交叉参考：

| 文件 | 摘要 |
| --- | --- |
| `src/bank.cpp` | 定义源银行、裂变银行的数据结构与序列化。 |
| `src/boundary_condition.cpp` | 实现真空、反射、周期等边界条件。 |
| `src/bremsstrahlung.cpp` | 描述制动辐射产生的次级粒子。 |
| `src/cell.cpp` | Cell 几何布尔运算、材料/宇宙填充逻辑。 |
| `src/cmfd_solver.cpp` | CMFD 加速器求解。 |
| `src/cross_sections.cpp` | 交叉截面 XML 解析与库索引。 |
| `src/dagmc.cpp` | DAGMC 几何接口。 |
| `src/distribution.cpp` | 通用分布基类实现。 |
| `src/distribution_angle.cpp` | 角度分布采样。 |
| `src/distribution_energy.cpp` | 能量分布采样。 |
| `src/distribution_multi.cpp` | 多维复合分布。 |
| `src/distribution_spatial.cpp` | 空间分布采样。 |
| `src/eigenvalue.cpp` | k 特征值迭代与裂变源管理。 |
| `src/endf.cpp` | ENDF 数据结构解析工具。 |
| `src/error.cpp` | 错误处理与日志。 |
| `src/event.cpp` | 事件驱动模拟框架。 |
| `src/external/quartic_solver.cpp` | 四次方程求解器（外部依赖）。 |
| `src/file_utils.cpp` | 文件系统工具。 |
| `src/finalize.cpp` | 模拟结束时释放资源。 |
| `src/geometry.cpp` | 几何搜索、重叠检测、邻居加速。 |
| `src/geometry_aux.cpp` | 几何辅助预处理。 |
| `src/hdf5_interface.cpp` | HDF5 读写封装。 |
| `src/initialize.cpp` | 模拟初始化流程。 |
| `src/lattice.cpp` | 规则/三角晶格实现。 |
| `src/main.cpp` | 程序入口、命令行解析。 |
| `src/material.cpp` | 材料类实现与核素装填。 |
| `src/math_functions.cpp` | 数值算法与特殊函数。 |
| `src/mcpl_interface.cpp` | MCPL 粒子文件读写。 |
| `src/mesh.cpp` | 网格体积计算与 tally 支撑。 |
| `src/message_passing.cpp` | MPI 通信封装。 |
| `src/mgxs.cpp` | 多群截面存储与组合。 |
| `src/mgxs_interface.cpp` | Python 端多群截面对接。 |
| `src/ncrystal_interface.cpp` | 晶体散射数据接口。 |
| `src/nuclide.cpp` | 核素数据加载与截面计算。 |
| `src/output.cpp` | 屏幕输出与日志格式化。 |
| `src/particle.cpp` | 粒子状态管理、碰撞预处理。 |
| `src/particle_data.cpp` | 全局粒子相关常量与缓存。 |
| `src/particle_restart.cpp` | 粒子重启文件读写。 |
| `src/photon.cpp` | 光子物理过程。 |
| `src/physics.cpp` | 碰撞调度与反应采样。 |
| `src/physics_common.cpp` | 碰撞共用算法。 |
| `src/physics_mg.cpp` | 多群碰撞算法。 |
| `src/plot.cpp` | 剖分输出（2D/3D Plot）。 |
| `src/position.cpp` | 位置与方向向量工具。 |
| `src/progress_bar.cpp` | 进度条显示。 |
| `src/random_dist.cpp` | 随机抽样工具。 |
| `src/random_lcg.cpp` | 线性同余随机数生成器。 |
| `src/random_ray/*.cpp` | Random Ray 求解器组件（源域、矩阵、仿真等）。 |
| `src/reaction.cpp` | 反应截面与概率计算。 |
| `src/reaction_product.cpp` | 二次粒子抽样。 |
| `src/scattdata.cpp` | 散射数据处理。 |
| `src/secondary_correlated.cpp` | 相关二次粒子产生。 |
| `src/secondary_kalbach.cpp` | Kalbach-Mann 分布采样。 |
| `src/secondary_nbody.cpp` | 多体二次粒子。 |
| `src/secondary_thermal.cpp` | 热散射二次粒子。 |
| `src/secondary_uncorrelated.cpp` | 非相关二次粒子。 |
| `src/settings.cpp` | 输入设置解析。 |
| `src/simulation.cpp` | 主模拟驱动与批次循环。 |
| `src/source.cpp` | 粒子源采样。 |
| `src/state_point.cpp` | 状态点 HDF5 输出。 |
| `src/string_utils.cpp` | 字符串工具。 |
| `src/summary.cpp` | 输入模型概述输出。 |
| `src/surface.cpp` | 面片定义与距离计算。 |
| `src/tallies/derivative.cpp` | 导数计分。 |
| `src/tallies/filter_*.cpp` | 各类过滤器实现。 |
| `src/tallies/tally.cpp` | 计分器主逻辑。 |
| `src/tallies/tally_scoring.cpp` | 计分调度。 |
| `src/tallies/trigger.cpp` | 触发器控制。 |
| `src/thermal.cpp` | 热散射核。 |
| `src/timer.cpp` | 计时器。 |
| `src/track_output.cpp` | 轨迹输出。 |
| `src/universe.cpp` | Universe 与 Cell 集合。 |
| `src/urr.cpp` | 未解析共振区采样。 |
| `src/volume_calc.cpp` | 几何体积蒙特卡洛估计。 |
| `src/weight_windows.cpp` | 权窗自适应控制。 |
| `src/wmp.cpp` | Windowed Multipole 截面计算。 |
| `src/xml_interface.cpp` | XML 解析辅助。 |
| `examples/simple_transport/main.cpp` | 教学版 3D 多群输运器。 |

> 如需逐文件深入阅读，建议配合本手册的章节指引、头文件声明以及 IDE 的跳转功能定位类与函数实现。

## 16. 结语

OpenMC 的 C++ 代码通过清晰的分层结构将复杂的核物理、几何建模与并行计算有机结合。本手册汇总了关键组件的职责与交互方式，结合 `examples/simple_transport/` 中的教学示例，可帮助读者迅速建立整体认知，并在阅读源码或贡献新特性时作为参照。
