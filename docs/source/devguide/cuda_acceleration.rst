.. _cuda-acceleration:

OpenMC CUDA 加速方案蓝图
==========================

本节面向希望在 OpenMC 中实现 GPU（以 CUDA 为代表）加速的开发者，给出一套兼顾
物理正确性与工程可落地性的总体方案。内容覆盖总体目标、内存布局、核函数划分、
随机数与统计量处理、构建流程以及验证策略，可作为实现 GPU 版本中子输运的路线
图。

设计目标与约束
----------------

* **兼容性**：保持与现有 CPU 版本一致的物理模型、输入格式和结果精度；允许在
  单机 CPU、单 GPU 与多 GPU 之间自由切换。
* **可维护性**：避免为 GPU 复制一份完全独立的代码路径，尽可能在数据结构与算
  法层面抽象出可在不同后端复用的组件。
* **扩展性**：支持后续增加 ROCm、SYCL 等其他后端；为多 GPU、分布式运行预留
  通信接口。

总体架构
--------

1. **分层调度**：在 C++ 层新增 ``transport::Driver`` 抽象，负责批次/世代调度；其
   中定义纯虚接口 ``launch_history_batch``，由 CPU 与 CUDA 子类分别实现。
2. **数据准备**：利用已有的 ``Particle``、``Material``、``Tally`` 数据结构构建
   GPU 友好的结构化数组（SoA），并通过一次性拷贝传输至显存。
3. **执行流程**：

   #. CPU 端完成源采样，将粒子初始状态写入页锁定内存；
   #. CUDA 驱动将粒子批量传入 GPU，执行历史并行的追踪核函数；
   #. 追踪完成后，返回宏观统计量（计数器、吸收/泄漏 tallies 等）。

数据布局策略
--------------

* **几何与截面**：

  * 将面向对象的几何 DAG 展平为 ``DeviceGeometry``，其中包含包围盒、曲面参数等
    数组，可由光线追踪算法直接索引。
  * 对连续能群截面使用 ``windowed multipole`` 或 ``log-log`` 插值表，将必要系数
    预编码为结构化数组并上传至常量/纹理内存，提高缓存命中。

* **粒子状态**：采用 SoA 布局存储 ``position``、``direction``、``energy``、
  ``material`` 等字段，配合压缩标记（bitmask）标识粒子存活状态，便于在 kernel
  内部进行流压缩（stream compaction）。

CUDA 核函数划分
----------------

* ``transport_history_kernel``：每个线程处理一条历史，核心步骤为几何追踪、选取反
  应类型、取样次级粒子与 tallies。为减少 warp 发散，可在 kernel 内将弹性散射、
  吸收等路径拆分为小状态机。
* ``surface_crossing_kernel``：对在一步步长内可能跨越多曲面的粒子使用子核函数配
  合共享内存存储候选曲面，实现批量求最小距离。
* ``tally_reduce_kernel``：对 tallies 的局部累加使用分层规约（warp-level shuffle +
  block-level shared memory），最后再与 CPU 侧累加。

核函数示例
-----------

.. code-block:: cuda

   __global__ void transport_history_kernel(DeviceState state,
                                            DeviceTallies tallies,
                                            curandStatePhilox4_32_10_t* rng) {
     const int tid = blockIdx.x * blockDim.x + threadIdx.x;
     if (tid >= state.num_active) return;

     Particle p = load_particle(state, tid);
     curandStatePhilox4_32_10_t local_rng = rng[tid];

     while (p.alive) {
       InteractionData xs = sample_macro_xs(p, state.data, local_rng);
       double distance = sample_distance(xs.total, local_rng);
       SurfaceHit hit = distance_to_boundary(p, distance, state.geometry);

       if (hit.crossed) {
         apply_surface_bc(p, hit);
         continue;
       }

       process_collision(p, xs, local_rng, tallies);
     }

     rng[tid] = local_rng;
     store_particle(state, tid, p);
   }

随机数与再现性
--------------

* 选用 ``curandStatePhilox4_32_10_t`` 作为基础 RNG，保证跨 GPU、跨批次的可重复性。
* 使用历史 ID（批次号、世代号、线程号）作为子序列种子，并与 CPU 版本的 XORWOW
  序列映射关系记录在案，便于交叉验证。

统计量与 Tallies
-----------------

* 在 GPU 端维护 ``DeviceTallies``，对区域剂量、通量、反应率等进行原子加。
* 对需要高精度的 tally（例如关联矩阵）采用 ``Kahan`` 补偿或 ``double-double``
  技巧，并在 kernel 末尾执行块级规约，降低原子操作开销。
* 支持 ``event-based`` 与 ``history-based`` tally 两种模式：前者在每次碰撞后
  写入，适合高并发；后者在历史结束时一次写入，便于保持与 CPU 版本一致。

内存管理与流控制
------------------

* 通过 ``cudaMallocAsync`` 与 ``cudaMemPool`` 管理长期驻留的几何、截面数据；粒子缓
  冲使用 ``cudaHostAlloc`` 分配的页锁定内存以提升传输效率。
* 利用多个 CUDA 流实现 ``copy-in → kernel → copy-out`` 的流水线，隐藏 PCIe 延迟。
* 若启用多 GPU，基于 MPI 的域分解与 ``cudaMemcpyPeerAsync`` 实现粒子交换。

构建与配置
------------

1. 在 ``CMakeLists.txt`` 中新增 ``OPENMC_USE_CUDA`` 选项，调用 ``find_package(CUDAToolkit)``
   并设置 ``target_compile_options`` 使 ``openmc`` 可同时生成 CPU 与 GPU 版本。
2. 提供 Python API 入口：``openmc.config['cuda'] = True``，在运行时通过 ``libopenmc``
   的 C API 切换后端。
3. 在 Dockerfile 中添加基于 ``nvidia/cuda`` 的构建 stage，自动安装驱动依赖。

验证与性能评估
----------------

* **单元测试**：对几何追踪、宏观截面取样、tally 累加等核心函数提供 GPU 版测试，
  与 CPU 结果逐项比对。
* **回归测试**：选取标准基准（C5G7、BEAVRS、JENDL pin cell 等），比较 keff、反
  应率向量及统计误差。
* **性能分析**：使用 ``Nsight Compute``、``Nsight Systems`` 量化占用率、内存带宽
  与热点 kernel；结合 Roofline 模型评估优化空间。

分阶段路线图
--------------

1. **阶段一（验证原型）**：完成单 GPU、均匀材料、简化几何的历史并行原型，跑通
   小规模基准并验证精度。
2. **阶段二（功能完善）**：加入多材料、S(a,β) 热散射、概率表等复杂物理，完善
   tally 与输出。
3. **阶段三（工程化）**：实现多 GPU 支持、与现有 MPI 并行融合、完成 CI 自动化
   与文档更新。
4. **阶段四（社区推广）**：发布性能数据、撰写用户指南，并针对关键硬件平台（如
   NVIDIA H100、A100）提供优化建议。

进一步参考
------------

* Liu et al., *Accelerating Monte Carlo Neutron Transport on GPUs*, Annals of
  Nuclear Energy, 2020.
* Hamilton et al., *Event-Based Algorithms for Monte Carlo Neutron Transport*,
  Nuclear Science and Engineering, 2016.

