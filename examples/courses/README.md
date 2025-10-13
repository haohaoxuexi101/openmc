# OpenMC 风格专题课程总览

本目录将之前零散的 Lesson 拓展为三套系统课程，每套不少于五节课，
通过可直接编译运行的 C++ 示例提炼 OpenMC 源码的核心设计思想。
若希望对 `particle_data.h` 与 `cell.h` 进行逐字段研读，请结合
《OpenMC 粒子与几何对象接口深度解析》文档与本目录的面向对象课程学习。

## 课程一：面向对象架构（oop_design）
| Lesson | 主题 | 亮点 |
| --- | --- | --- |
| lesson01_geometry_state.cpp | 几何状态建模 | Surface/Cell 继承 + GeometryState 组合 |
| lesson02_particle_lifecycle.cpp | 粒子状态生命周期 | move/scatter/absorb 接口守护状态一致性 |
| lesson03_navigation_context.cpp | 定位上下文解耦 | Context 持有 Navigator，演示接口分离 |
| lesson04_material_binding.cpp | 几何-材料组合 | MaterialRegistry 绑定 Cell 材料 |
| lesson05_interface_strategy.cpp | 策略模式 | 计分策略抽象，映射 Tally 设计 |
| lesson06_particle_data_bridge.cpp | ParticleData 接口桥接 | 解释 GeometryState 继承 + 状态字段布局 |
| lesson07_cell_interface_walkthrough.cpp | Cell 布尔树拆解 | Region 多态与 distance 接口详解 |

编译示例：
```bash
g++ -std=c++17 -O2 -o lesson01 examples/courses/oop_design/lesson01_geometry_state.cpp
```

## 课程二：输运算法实践（transport_algorithms）
| Lesson | 主题 | 亮点 |
| --- | --- | --- |
| lesson01_ray_tracing.cpp | 射线追踪 | 多平面最近交距，解释几何寻道 |
| lesson02_collision_sampling.cpp | 碰撞抽样 | Σ 表驱动随机分支 |
| lesson03_k_eigenvalue.cpp | k-eigenvalue 幂迭代 | 多群裂变-吸收反馈 |
| lesson04_track_tally.cpp | 轨迹长度计数 | Cell 级路径积分 |
| lesson05_neutron_physics.cpp | 多群核数据 | GroupPhysics 组织吸收/裂变数据 |

编译示例：
```bash
g++ -std=c++17 -O2 -o lesson_transport examples/courses/transport_algorithms/lesson01_ray_tracing.cpp
```

## 课程三：并行与高性能（parallel_computing）
| Lesson | 主题 | 亮点 |
| --- | --- | --- |
| lesson01_openmp_threads.cpp | OpenMP 粒子切分 | 静态调度演示 |
| lesson02_mpi_bank_reduce.cpp | MPI 裂变银行归并 | Allreduce 汇总粒子数 |
| lesson03_thread_local_rng.cpp | 线程私有 RNG | seed_seq 避免随机相关 |
| lesson04_hybrid_exchange.cpp | 混合并行数据交换 | OpenMP 聚合 + MPI 汇总 |
| lesson05_work_partition.cpp | 负载划分策略 | 均匀划分 Bank 的工具函数 |
| lesson06_async_pipeline.cpp | 非阻塞裂变流水线 | MPI Isend/Irecv 与计算重叠 |
| lesson07_dynamic_tasking.cpp | OpenMP 动态任务 | 粒子工作窃取与负载均衡 |

编译示例：
```bash
g++ -std=c++17 -O2 -fopenmp -o lesson_parallel examples/courses/parallel_computing/lesson01_openmp_threads.cpp
mpicxx -std=c++17 -O2 -o lesson_mpi examples/courses/parallel_computing/lesson02_mpi_bank_reduce.cpp
```

> 提示：MPI 课程需在安装 MPI 的环境中编译运行，OpenMP 课程需使用支持 `-fopenmp` 的编译器。
