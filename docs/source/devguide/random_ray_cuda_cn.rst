.. _random_ray_cuda_impl_cn:

============================================
随机射线 CUDA 技术规范（中文版）
============================================

概述
====

CUDA 路径通过启动每线程处理单个能群的设备内核，加速
:func:`openmc::RandomRay::attenuate_flux_flat_source` 中的平源区衰减例程。
构建时可使用 :cmake:`OPENMC_USE_CUDA` 缓存选项来启用 CUDA 支持，该选项
会将 CUDA 设为一等语言、链接 ``CUDA::cudart``，并把
``src/random_ray/random_ray_cuda.cu`` 编译进核心共享库（参见
:file:`CMakeLists.txt`）。

公开接口
========

主机与设备的桥接接口定义在 ``include/openmc/random_ray/random_ray_cuda.h`` 中，
提供两个入口函数：

``bool openmc::random_ray::cuda::available()``
    若运行时检测到支持 CUDA 的设备且编译时启用了 ``OPENMC_USE_CUDA``，
    则返回 ``true``。该检查在函数首次调用时惰性执行（参见
    :file:`include/openmc/random_ray/random_ray_cuda.h` 与
    :file:`src/random_ray/random_ray_cuda.cpp`）。

``void openmc::random_ray::cuda::attenuate_flat_source(...)``
    在把宏观截面、源项以及角通量向量复制到设备后执行 GPU 内核。若 CUDA
    支持不可用，封装器会抛出致命错误以防止混用不一致的构建（参见
    :file:`include/openmc/random_ray/random_ray_cuda.h` 与
    :file:`src/random_ray/random_ray_cuda.cpp`）。

运行时行为
==========

CUDA 翻译单元维护线程局部的缓冲区缓存，因此重复调用时无需重新分配设备内存。
只有当所需能群数超过已有容量时才会调用 :func:`cudaMalloc`。
所有缓冲区会在线程结束或容量扩张时自动释放（参见
:file:`src/random_ray/random_ray_cuda.cu`）。

``available()`` 会惰性地查询 ``cudaGetDeviceCount`` 并记录是否存在至少一块设备。
若 CUDA 运行时返回除 ``cudaErrorNoDevice`` 之外的错误，例程会通过
``cudaGetLastError`` 清除状态，以便后续 API 调用继续进行。因此即使在缺少驱动的
系统上，主机代码也可以安全地轮询该函数（参见
:file:`src/random_ray/random_ray_cuda.cu`）。

内核语义
========

``flat_source_kernel`` 为每个线程分配一个能群，利用共享的指数逼近计算更新后的角通量
与 Delta-psi。主机封装器会复制输入数据到设备、以 128 个线程为一块发射内核，并在
CPU 获得标量通量累积互斥锁之前将结果复制回主机（参见
:file:`src/random_ray/random_ray_cuda.cu` 与
:file:`src/random_ray/random_ray.cpp`）。

与随机射线求解器的集成
========================

``RandomRay::attenuate_flux_flat_source`` 会在每次调用时检查
``cuda::available()``。当设备路径成功时，该函数跳过 CPU 标量循环，直接执行源区簿记；
若不存在可用设备或材料为空腔，则继续沿用原有的 CPU 算法（参见
:file:`src/random_ray/random_ray.cpp`）。

错误处理
========

所有 CUDA 运行时调用均被 ``check_cuda`` 包裹，任何失败都会携带 CUDA 诊断字符串触发致命错误。
如果未启用 ``OPENMC_USE_CUDA`` 却调用加速例程，同样会抛出致命错误，确保用户不会在物理结果过时的
混合构建上运行计算（参见 :file:`src/random_ray/random_ray_cuda.cu` 与
:file:`src/random_ray/random_ray_cuda.cpp`）。
