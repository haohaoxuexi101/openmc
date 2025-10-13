# 三维随机射线多群求解器

本示例目录提供一个可直接编译运行的 **Random Ray Monte Carlo** 求解器，通过纯 C++
复刻 OpenMC 的面向对象划分，将材料、几何、粒子、计分等模块解耦并支持多群多区域输运。
所有源码均附带详细中文注释，便于学习 OpenMC 的设计理念。

## 功能特性

- 任意群数的宏观截面输入，支持吸收/散射/裂变与平均中子产生数 `nu`。
- 三维轴对齐盒体几何，可按嵌套顺序堆叠构成燃料芯块、慢化剂、反射层等多区域体系。
- 体积分布的随机射线采样，使用曲线长度计分器统计区域-能群通量、吸收率、散射率、裂变率。
- 原位抽样裂变二次粒子，采用粒子银行深度优先追踪所有二次历史。
- 可通过输入文件配置历史数、随机种子、最大追踪步数以及源项能谱。

## 编译与运行

```bash
cd examples/random_ray_solver
g++ -std=c++17 -O2 -o random_ray_solver main.cpp
./random_ray_solver --input config/example_reactor.rr --histories 50000
```

可选参数：

- `--input <path>` 指定配置文件路径（默认为 `config/example_reactor.rr`）。
- `--histories <N>` 覆盖输入文件中的历史数设置。
- `--seed <S>` 自定义随机数种子，便于重复实验。

## 输入文件格式

输入文件采用分节语法：

- `[global]`：总群数、历史数、随机种子、最大步数等全局参数。
- `[material <name>]`：材料截面，使用 `sigma_s_rowX` 表示散射矩阵的第 `X` 行。
- `[region <name>]`：定义区域的包围盒与所用材料，区域按照出现顺序自动设置优先级。
- `[source <name>]`：可选源项设置，包含粒子抽样区域及能群概率分布。

参考 `config/example_reactor.rr` 可快速创建新的系统组合。

## 与 OpenMC 的设计映射

| 本示例模块 | OpenMC 对应模块 | 说明 |
| --- | --- | --- |
| `Material` | `openmc::Material` | 保存多群宏观截面并提供 `sigma_t`、`nu`、`chi` 等接口 |
| `Geometry` + `Region` | `openmc::Geometry` / `openmc::Cell` | 采用嵌套盒体描述体包含关系，并提供 `locate`/`distance_to_boundary` |
| `Particle` / `ParticleBank` | `openmc::Particle` / `Bank` | 粒子状态封装及裂变二次粒子堆栈 |
| `Tallies` | `openmc::Tally` + `Filter` | 记录区域-能群二维计分表，实现曲线长度估计 |
| `RandomRaySolver` | `openmc::Simulation` + `openmc::physics::transport_particle` | 封装求解流程、碰撞处理与粒子源抽样 |

通过阅读 `main.cpp` 可以直观理解 OpenMC 如何在 C++ 中运用组合优先的对象模型组织复杂的核输运仿真。
