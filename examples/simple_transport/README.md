# 三维多群中子输运示例（OpenMC 风格版）

本示例给出一个可直接编译运行的 **三维多群 k-eigenvalue** 程序。代码严格对齐 OpenMC
的面向对象层次：材料（`Material`）、表面/单元/宇宙/晶格（`Surface`、`Cell`、`Universe`、
`Lattice3D`、`Geometry`）、粒子状态与随机引擎（`Particle`、`Bank`、`RandomEngine`）、
物理驱动（`PhysicsDriver`）以及顶层模拟封装（`Simulation`）。与上一版一维示例相比，
此版本实现了**三维立方晶格的体元定位、真实的方向余弦输运以及基于粒子银行的 k_eff
估计**，可直接映射到 OpenMC 源码的核心组织方式。

## 编译与运行

```console
cd examples/simple_transport
g++ -std=c++17 -O2 -o simple_transport main.cpp
./simple_transport
```

典型输出（随机数不同会有微小差异）：

```
第 1 代有效增殖因子估计 k-effective = 1.16742
...
按单元与能群划分的路径长度计数 (cm)：
  单元 fuel_lower_fresh： 群0=... 群1=... 群2=...
  单元 fuel_lower_burned： 群0=... 群1=... 群2=...
  ...
总泄漏权重：...
```

输出中可观察到三维网格中不同体元、不同能群的路径长度估计以及代内的 k_eff 收敛。

### 特性拆解学习代码

> 为了更好地理解 OpenMC 各个模块的职责，本目录新增 `lessons/` 子目录，将主要特性拆解
> 成可独立编译运行的小程序。每个示例都配套中文注释，展示接口设计、组合与继承的实
> 践方式。

| 文件 | 学习目标 | 编译示例 |
| ---- | -------- | -------- |
| `lessons/lesson_vec3.cpp` | 解析三维向量工具类 `Vec3` 的值语义与算子设计。 | `g++ -std=c++17 -O2 -o lesson_vec3 lessons/lesson_vec3.cpp` |
| `lessons/lesson_material.cpp` | 拆解 `Material` 的多群截面存储、总截面计算与碰撞抽样。 | `g++ -std=c++17 -O2 -o lesson_material lessons/lesson_material.cpp` |
| `lessons/lesson_geometry.cpp` | 演示 `Surface`/`Cell`/`Universe` 的组合关系与边界距离计算。 | `g++ -std=c++17 -O2 -o lesson_geometry lessons/lesson_geometry.cpp` |
| `lessons/lesson_particle_bank.cpp` | 聚焦 `Particle` 与 `Bank` 组合、粒子银行重采样逻辑。 | `g++ -std=c++17 -O2 -o lesson_bank lessons/lesson_particle_bank.cpp` |
| `lessons/lesson_physics.cpp` | 将材料、几何、粒子、随机引擎整合成迷你 `PhysicsDriver`。 | `g++ -std=c++17 -O2 -o lesson_physics lessons/lesson_physics.cpp` |
| `lessons/lesson_random.cpp` | 展示随机引擎封装、均匀采样、各向同性方向与离散分布。 | `g++ -std=c++17 -O2 -o lesson_random lessons/lesson_random.cpp` |
| `lessons/lesson_lattice.cpp` | 讲解 Universe + Lattice3D + Geometry 的三维定位路径。 | `g++ -std=c++17 -O2 -o lesson_lattice lessons/lesson_lattice.cpp` |
| `lessons/lesson_tally.cpp` | 细化 track-length tally 的路径长度积累与泄漏统计。 | `g++ -std=c++17 -O2 -o lesson_tally lessons/lesson_tally.cpp` |
| `lessons/lesson_simulation.cpp` | 演示 Simulation 设置、物理驱动、重采样与统计的顶层协同。 | `g++ -std=c++17 -O2 -o lesson_simulation lessons/lesson_simulation.cpp` |
| `lessons/lesson_mgxs_builder.cpp` | 复刻多群截面制作流程：核素细能格折合与材料混配。 | `g++ -std=c++17 -O2 -o lesson_mgxs lessons/lesson_mgxs_builder.cpp` |
| `lessons/lesson_scatter_anisotropic.cpp` | 对照 OpenMC 的多群与连续能量各向异性散射抽样策略。 | `g++ -std=c++17 -O2 -o lesson_scatter lessons/lesson_scatter_anisotropic.cpp` |
| `lessons/lesson_distribution_sampling.cpp` | 拆解分布接口、分段线性抽样、Alias 表与复合抽样。 | `g++ -std=c++17 -O2 -o lesson_dist lessons/lesson_distribution_sampling.cpp` |
| `lessons/lesson_neighbor_table.cpp` | 展示笛卡尔晶格的邻接表构建与边界挂接策略。 | `g++ -std=c++17 -O2 -o lesson_neighbor lessons/lesson_neighbor_table.cpp` |
| `lessons/lesson_doppler_broadening.cpp` | 解释多普勒展宽的高斯卷积思想与温度依赖。 | `g++ -std=c++17 -O2 -o lesson_doppler lessons/lesson_doppler_broadening.cpp` |
| `lessons/lesson_probability_table.cpp` | 演示未解析共振区的概率表抽样机制。 | `g++ -std=c++17 -O2 -o lesson_prob lessons/lesson_probability_table.cpp` |

所有小程序均可直接运行，并在控制台输出关键步骤，帮助读者由浅入深掌握 OpenMC 风格的
面向对象实现。

## 理论说明：如何将 Monte Carlo 物理映射到代码

> 下面每一节都按“物理理论 ➜ 代码落地”展开，帮助读者理解 OpenMC 为何采用这种
> 对象划分，并在示例中通过中文注释完整复现。

### 1. 多群输运与宏观截面

- **理论背景**：多群方法将连续能谱离散为有限群，输运方程中的吸收、裂变、散射截面
  在每个能群上取宏观平均。总截面满足
  \( \Sigma^g_t = \Sigma^g_a + \Sigma^g_f + \sum_{g'} \Sigma^{g \to g'}_s \)。
- **代码映射**：`Material` 在构造时接收 `sigma_a/sigma_f/nu/sigma_s/chi`，并通过
  [`sigma_t`](main.cpp) 内的求和实现上述公式；碰撞处理时，`PhysicsDriver::collide`
  按照总截面比例将吸收、散射、裂变分支概率映射到随机抽样，完全复刻 OpenMC 中
  `collision()` 的写法。
- **设计要点**：所有数据使用 `std::vector<double>`，既直观又支持与真实截面库接口，
  也是 OpenMC 在 C++ 核心中常用的存储形式。

### 2. 三维几何与体元定位

- **理论背景**：在 OpenMC 中，几何由一系列边界曲面（Surface）通过布尔组合定义。
  对于笛卡尔晶格，核心操作是对粒子当前位置与方向求解到相邻平面的距离，并判断是
  否穿越界面。三维情形下，需同时对 `x/y/z` 三个方向计算自由程。
- **代码映射**：
  - `Surface` 枚举了 `Axis`、`Boundary` 两种属性，对应 OpenMC 中平面定义和边界类型。
  - `Cell::distance_to_boundary` 按照方向余弦分别计算到正、负方向平面的距离，取最小
    值即得到下一个交点；这与 OpenMC `surface_t::distance()` 的做法一致。
  - `Lattice3D` 提供 `universe_at()`，先将粒子坐标投影到晶格索引，再返回对应宇宙；
    `Geometry::locate` 则先查询晶格，再回退到根宇宙，完整模拟了 OpenMC 的定位流程。
  - 在示例几何中，`DemoModel` 构建了 2×2×2 的立方晶格，每个单元对应一个 Universe，
    并通过循环自动建立六个方向的邻接指针，等价于 OpenMC 中的网格邻接链接。

### 3. 粒子输运、散射与裂变抽样

- **理论背景**：中子自由程满足指数分布 \( s = -\ln(1-\xi)/\Sigma_t \)。当粒子到达界面
  时根据边界条件处理：真空吸收、反射改变方向、界面则进入相邻体元。碰撞处按截面
  权重决定吸收/散射/裂变。
- **代码映射**：
  - `PhysicsDriver::transport_particle` 中首先抽样自由程 `free_path`，再与最近界面距离
    `hit.distance` 取最小值；若界面更近，调用 `handle_boundary`，否则进入 `collide`。
  - `handle_boundary` 对真空、反射、界面三种情况分别处理：反射通过翻转对应方向余弦
    并做 `eps` 级位移，确保不会反复击中同一平面；界面则切换到相邻 `Cell`。
  - `collide` 依据截面比例实现多分支抽样；`scatter` 重采样能群与方向，`fission` 通过
    `std::discrete_distribution` 根据裂变谱抽样能群，并将二次粒子写入裂变库，完全对标
    OpenMC 在 `physics/fission.cpp` 中的策略。

### 4. k-eigenvalue 迭代与粒子银行

- **理论背景**：K-eigenvalue 问题常用 power iteration。每一代追踪固定数量的源粒子，
  将裂变产生的粒子等概率重采样作为下一代源粒子，k_eff 估计为“本代产生的裂变中子数
  / 初始粒子数”。
- **代码映射**：`PhysicsDriver::transport_generation` 返回 `produced/particles` 即代内
  k_eff；`resample_source` 则从裂变库中均匀抽样替换源库粒子，保持总粒子数不变。该
  逻辑对应 OpenMC `power_iteration()` 中的 `source_bank`/`fission_bank` 交换过程。

### 5. 统计量与 tallies

- **理论背景**：OpenMC 中的 track-length tally 利用 \( \phi \approx \sum w \cdot l / V \)。
- **代码映射**：`TallyManager` 存储 `(cell, group)` 的路径长度累积；`Simulation::run`
  在输出时除以代数，即得到平均轨迹长度。若需真实通量，还可再除体积，本示例保留
  了进一步扩展的接口。

### 6. 3D 几何的数值稳定技巧

- 在界面处理时对位置进行 `1e-9` 的微调，避免浮点误差导致粒子卡在界面；这是 OpenMC
  中常见的“epsilon 推离”技巧。
- `distance_to_boundary` 中忽略绝对值小于 `1e-14` 的方向余弦，避免除零；若粒子严格平行
  某轴，则忽略对应平面，保证算法稳定。
- 新产生的裂变粒子位置加上 `dir * 1e-7` 的偏移，避免立即与原界面重合。

## 编程技巧与类设计指南

| 模块 | 示例代码中的类 | 说明 |
| ---- | --------------- | ---- |
| 向量工具 | `Vec3` | 提供向量加减与归一化，简化方向余弦运算。 |
| 材料 | `Material` | 构造函数断言截面长度一致，暴露 `sigma_t()` 等辅助函数。 |
| 几何 | `Surface` / `Cell` / `Universe` / `Lattice3D` / `Geometry` | 组合成完整的三维几何定位链路。`Cell` 内部保存六个方向的邻居指针，便于在界面切换。 |
| 粒子与银行 | `Particle` / `Bank` | 轻量化封装当前状态，银行直接使用 `std::vector` 支持范围 for。 |
| 随机数 | `RandomEngine` | 内部使用 `std::mt19937`，提供均匀分布与各向同性方向抽样接口。 |
| 物理驱动 | `PhysicsDriver` | 划分为 `transport_particle`、`handle_boundary`、`collide` 等函数，代码结构与 OpenMC 的模块划分一致。 |
| 统计与主循环 | `TallyManager` / `Simulation` | 构造时初始化 tallies，`run()` 中循环迭代、重采样、输出统计。 |

## 继承、组合与接口教学案例

1. **仿真实例中的组合关系**：`Surface`、`Cell`、`Universe`、`Lattice3D`、`Geometry` 五个类彼此组合构成体元树，每个 `Cell` 同时持有六个面指针和邻接指针，体现“组合承载结构、接口暴露行为”的思路。【F:examples/simple_transport/main.cpp†L96-L280】
2. **策略化边界接口**：`Cell::BoundaryHit` 负责存储最近碰撞信息，`PhysicsDriver::handle_boundary` 根据 `Surface::boundary()` 的不同枚举实现真空、反射或跨界策略，可据此练习新增“白反射”或“周期”分支。【F:examples/simple_transport/main.cpp†L151-L479】
3. **接口驱动物理循环**：`PhysicsDriver::transport_generation` 将单代追踪封装为接口，内部再调用 `transport_particle` 和 `collide`，方便替换散射模型或计分逻辑，同时 `TallyManager` 被组合到 `Simulation` 中集中管理统计量，完全贴合 OpenMC 的接口划分方式。【F:examples/simple_transport/main.cpp†L371-L620】
4. **进一步的继承练习**：可在本示例基础上新增一个继承自 `Surface` 的“斜切平面”类，或继承 `PhysicsDriver` 改写 `collide` 函数，实现自定义的散射核；通过替换派生类即可验证多态接口的威力。【F:examples/simple_transport/main.cpp†L96-L620】

**面向对象划分的好处**：

1. **职责单一**：每个类只关注自身职责，例如 `Cell` 只负责几何信息，碰撞逻辑完全由
   `PhysicsDriver` 决定，便于测试与替换。
2. **易于扩展**：想要增加圆柱表面、能谱展开或并行化时，无需改动其他模块；只要保持
   相同接口即可。
3. **贴合 OpenMC 源码**：阅读 OpenMC 时，可直接将示例中的中文注释与 OpenMC 对应文件
   对照，降低入门门槛。

## 示例几何：2×2×2 立方晶格

- 下层四个单元：两个燃料块（新鲜燃料、燃耗燃料）和两块反射层。
- 上层四个单元：前半部为慢化剂，后半部为反射层，实现轴向与方位异质性。
- `z=0` 平面设为反射边界，模拟对称面；`x=0/4`、`y=0/4`、`z=6` 为真空边界，允许粒子逸出。
- 源位置位于燃料区内部 `{1.2, 0.8, 1.5}`，对应 `DemoModel` 最终在 `main()` 中设置。

邻接关系通过双层循环自动建立，保证每个体元的 X/Y/Z 三个方向都有正确的邻居指针；
这与 OpenMC 在构建笛卡尔晶格时自动填充邻接映射的逻辑一致。

## 可继续深入的练习

1. **加入体积计算与通量归一**：在 `Cell` 中补充体积属性，输出时除以体积即可得到真实
   通量估计。
2. **实现吸收截面温度依赖**：在 `Material` 中存储温度表，使用插值生成 `sigma_*`。
3. **扩展 tallies**：仿照 OpenMC，在 `TallyManager` 中增加能谱、泄漏角分布等统计量。
4. **并行化**：为 `PhysicsDriver::transport_generation` 添加 OpenMP 并行循环或任务队列，
   体验大规模 Monte Carlo 的性能优化。

通过阅读与运行本示例，读者可以看到 OpenMC 面向对象风格在三维情形下如何组织：从
材料数据、几何定位到粒子输运与统计量输出，每一步都有对应类负责，且全部配有详尽
中文注释，便于对照官方源码进一步深入学习。
