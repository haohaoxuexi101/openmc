.. _random_ray_cuda_theory_cn:

============================================
随机射线 CUDA 理论补充（中文版）
============================================

.. note::
   本文档为 :ref:`random_ray_cuda_theory` 的中文翻译，
   介绍在启用 ``OPENMC_USE_CUDA`` 选项后 Random Ray 求解器
   对平源区（flat source region）执行 GPU 加速衰减计算的
   物理背景与数值实现。

动机
====

Random Ray 求解器在每个平源区内需要求解一维输运方程来衰减角通量：

.. math::

   \psi_{\text{out}} = \psi_{\text{in}} e^{-\Sigma_t s}
   + \frac{q}{\Sigma_t} \left(1 - e^{-\Sigma_t s}\right),

其中 :math:`s` 为射线段长度，:math:`\Sigma_t` 是对应能群的总宏观截面，
而 :math:`q` 为各向同性体源项。在 CPU 实现中，用于统计量的贡献为

.. math::

   \Delta \psi = \left(\psi_{\text{in}} - \frac{q}{\Sigma_t}\right)
   \left(1 - e^{-\Sigma_t s}\right),

该项被累加到平源区的新标量通量估计中。当每批次需要追踪成千上万条射线
并覆盖上百个能群时，衰减循环会成为运行时间瓶颈。通过将该循环卸载到
CUDA，可以在保持相同数学形式的同时，充分利用 GPU 的 SIMT 执行模型
来摊薄计算代价。

:math:`1 - e^{-\tau}` 的有理逼近
=================================

CUDA 内核若直接调用 ``<cmath>`` 中的 :math:`\exp` 会带来额外的设备端
代码开销与精度损失。因此求解器采用 C.J. Josey 等人提出的平源射线追踪
有理逼近。辅助函数 :func:`openmc::cjosey_exponential` 使用七阶分子和
分母多项式来评估 :math:`1 - e^{-\tau}`，主机与设备端共享完全一致的
系数（参见 :file:`include/openmc/random_ray/exponential.h`）。GPU 内核同样
调用该函数，从而保证 CPU 与 GPU 衰减路径的数值一致性（参见
:file:`src/random_ray/random_ray_cuda.cu`）。

单精度工作数组
================

随机射线求解器将角通量与 Delta-psi 工作数组存放在 ``std::vector<float>``
中，以匹配多群数据的布局。CUDA 内核以单精度缓冲区接收源项、角通量与
Delta-psi，但对宏观截面仍保持双精度，以避免在计算
:math:`\tau = \Sigma_t s` 时出现下溢。这与 CPU 实现保持一致：CPU 同样先
用双精度截面乘以段长，然后在执行多项式评估前转换为 ``float``（参见
:file:`src/random_ray/random_ray.cpp` 与
:file:`src/random_ray/random_ray_cuda.cu`）。

激活段与非激活段
==================

随机射线轨迹通过“死亡”（非激活）与“存活”（激活）阶段来避免源迭代偏差。
CUDA 加速仅在对统计量有贡献的激活段中生效。一旦射线仍处于非激活阶段，
或处于空腔材料中，CPU 例程会继续处理衰减并保护免受虚空区域影响。由于
两条代码路径都写入同一角通量与 Delta-psi 数组，GPU 内核可以透明地为
:func:`openmc::RandomRay::attenuate_flux_flat_source` 后续的标量通量累积
提供数据。若射线处于空腔材料或能群数为零，则完全绕过 GPU，回退到标量
CPU 循环（参见 :file:`src/random_ray/random_ray.cpp` 与
:file:`src/random_ray/random_ray_cuda.cu`）。

故障处理
========

若运行环境缺少支持 CUDA 的设备，或 OpenMC 在构建时未启用
``OPENMC_USE_CUDA``，求解器会在接触任何设备缓冲区前报告致命错误（参见
:file:`src/random_ray/random_ray_cuda.cpp` 与
:file:`src/random_ray/random_ray_cuda.cu`）。这能确保生产计算不会在物理
实现不匹配时悄然混用 CPU 与 GPU 结果。
