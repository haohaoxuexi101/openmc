.. _random_ray_cuda_theory:

=====================================
Random Ray CUDA Theory Supplement
=====================================

.. note::
   本节补充说明 Random Ray 求解器中平源区（flat source region）在 GPU 上的
   传输模型。它与 :ref:`random_ray` 主章节的推导一致，但聚焦于
   ``OPENMC_USE_CUDA`` 选项启用时的指数衰减近似与数据流。

Motivation
==========

The random ray solver attenuates the angular flux in each flat source region by
solving the 1-D transport equation

.. math::

   \psi_{\text{out}} = \psi_{\text{in}} e^{-\Sigma_t s}
   + \frac{q}{\Sigma_t} \left(1 - e^{-\Sigma_t s}\right),

where :math:`s` is the ray segment length, :math:`\Sigma_t` is the total
macroscopic cross section for the energy group, and :math:`q` is the isotropic
volume source.  In the CPU implementation the contribution recorded for tallies
is

.. math::

   \Delta \psi = \left(\psi_{\text{in}} - \frac{q}{\Sigma_t}\right)
   \left(1 - e^{-\Sigma_t s}\right),

which is accumulated into the new scalar flux estimator for the flat source
region.  Transporting thousands of rays per batch with hundreds of energy groups
causes this attenuation loop to dominate runtime.  Offloading the loop to CUDA
preserves the same mathematics while amortizing the work over the GPU's SIMT
execution model.

Rational Approximation of :math:`1 - e^{-\tau}`
===============================================

The CUDA kernel cannot call :math:`\exp` from ``<cmath>`` directly without
incurring extra device code and precision penalties.  Instead the solver relies
on the rational approximation introduced by C.J. Josey et al. for flat source
ray tracing.  The helper function :func:`openmc::cjosey_exponential` evaluates
this approximation for :math:`1 - e^{-\tau}` using a seventh-order numerator and
denominator polynomial, with identical coefficients shared between host and
device builds (see :file:`include/openmc/random_ray/exponential.h`).  The same
function is invoked inside the GPU kernel so that CPU and GPU attenuation paths
stay numerically consistent (see :file:`src/random_ray/random_ray_cuda.cu`).

Single-precision Working Arrays
===============================

The random ray solver stores its angular flux and delta-psi work arrays in
``std::vector<float>`` to match the multigroup data layout.  The CUDA kernel
accepts single-precision buffers for the source, angular flux, and delta-psi, but
keeps macroscopic cross sections in double precision to avoid underflow when
computing :math:`\tau = \Sigma_t s`.  This mirrors the CPU implementation, which
also multiplies double-precision cross sections by the segment length before
casting to ``float`` for the polynomial evaluation (compare
:file:`src/random_ray/random_ray.cpp` and :file:`src/random_ray/random_ray_cuda.cu`).

Active vs. Inactive Path Lengths
================================

Random ray tracks maintain "dead" (inactive) and "live" (active) phases to avoid
biasing the source iteration.  CUDA acceleration only applies during segments
that contribute to tallies, i.e., once a ray becomes active.  The CPU routine
still governs dead-zone attenuation and guards against void regions.  Because
both code paths write through the same angular flux and delta-psi arrays, the GPU
kernel transparently feeds the subsequent scalar flux accumulation in
:func:`openmc::RandomRay::attenuate_flux_flat_source`.  Rays in void material or
with zero energy groups bypass the device entirely and fall back to the scalar
CPU loop (see :file:`src/random_ray/random_ray.cpp` and
:file:`src/random_ray/random_ray_cuda.cu`).

Failure Handling
================

If no CUDA-capable device is present, or if OpenMC was compiled without
``OPENMC_USE_CUDA``, the solver reports a fatal error before any device buffers
are touched (see :file:`src/random_ray/random_ray_cuda.cpp` and
:file:`src/random_ray/random_ray_cuda.cu`).  This ensures that production
simulations never silently mix mismatched physics implementations.
