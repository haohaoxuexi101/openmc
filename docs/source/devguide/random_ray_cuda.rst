.. _random_ray_cuda_impl:

====================================
Random Ray CUDA Technical Specification
====================================

Overview
========

The CUDA path accelerates the flat-source attenuation routine inside
:func:`openmc::RandomRay::attenuate_flux_flat_source` by launching a device kernel
that processes one energy group per thread.  Compilation of the CUDA sources is
controlled by the :cmake:`OPENMC_USE_CUDA` cache option, which enables CUDA as a
first-class language, links ``CUDA::cudart``, and compiles
``src/random_ray/random_ray_cuda.cu`` into the core shared library (see
:file:`CMakeLists.txt`).

Public Interface
================

The host/device bridge is declared in
``include/openmc/random_ray/random_ray_cuda.h``.  Two entry points exist:

``bool openmc::random_ray::cuda::available()``
    Returns ``true`` when a CUDA-capable device is detected at runtime and the
    code was built with ``OPENMC_USE_CUDA``.  The check happens lazily the first
    time the function is called (see :file:`include/openmc/random_ray/random_ray_cuda.h`
    and :file:`src/random_ray/random_ray_cuda.cpp`).

``void openmc::random_ray::cuda::attenuate_flat_source(...)``
    Executes the GPU kernel after copying the macroscopic cross sections, source
    term, and angular flux vectors to the device.  The wrapper throws a fatal
    error when CUDA support is absent to guard against mismatched builds (see
    :file:`include/openmc/random_ray/random_ray_cuda.h` and
    :file:`src/random_ray/random_ray_cuda.cpp`).

Runtime Behavior
================

The CUDA translation unit maintains a thread-local buffer cache so repeated
invocations avoid re-allocating device memory.  Calls to :func:`cudaMalloc` only
occur when the requested number of energy groups exceeds the cached capacity.
All buffers are released automatically when the thread exits or when the
capacity grows (see :file:`src/random_ray/random_ray_cuda.cu`).

``available()`` lazily queries ``cudaGetDeviceCount`` and records whether at
least one device exists.  If the CUDA runtime reports an error other than
``cudaErrorNoDevice``, the status is cleared via ``cudaGetLastError`` so that
subsequent API calls can proceed.  Host code can therefore poll the function
safely even on systems without drivers (see :file:`src/random_ray/random_ray_cuda.cu`).

Kernel Semantics
================

``flat_source_kernel`` assigns one energy group to each thread and computes the
updated angular flux and delta-psi using the shared exponential approximation.
The host wrapper copies inputs to the device, launches blocks of 128 threads, and
copies results back before the scalar flux accumulation lock is acquired on the
CPU (see :file:`src/random_ray/random_ray_cuda.cu` and
:file:`src/random_ray/random_ray.cpp`).

Integration into the Random Ray Solver
======================================

``RandomRay::attenuate_flux_flat_source`` checks ``cuda::available()`` on every
call.  When the device path succeeds the function skips the scalar CPU loop and
continues directly to the bookkeeping guarded by the source-region mutex.  If
no device is present, or if the material is void, the original CPU algorithm is
used unchanged (see :file:`src/random_ray/random_ray.cpp`).

Error Handling
==============

All CUDA runtime calls are wrapped in ``check_cuda`` so that any failure emits a
fatal error with the CUDA diagnostic string.  Invoking the acceleration routines
without enabling ``OPENMC_USE_CUDA`` likewise triggers a fatal error, ensuring
users never inadvertently run a hybrid build with stale physics results (see
:file:`src/random_ray/random_ray_cuda.cu` and
:file:`src/random_ray/random_ray_cuda.cpp`).
