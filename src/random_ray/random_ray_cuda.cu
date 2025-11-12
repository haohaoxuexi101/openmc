#ifdef OPENMC_USE_CUDA

#include "openmc/random_ray/random_ray_cuda.h"

#include "openmc/error.h"
#include "openmc/random_ray/exponential.h"

#include <cuda_runtime.h>

#include <string>

namespace openmc {
namespace random_ray {
namespace cuda {
namespace detail {

namespace {

struct DeviceState {
  bool initialized {false};
  bool available {false};
};

inline DeviceState& state()
{
  static DeviceState s;
  return s;
}

inline void check_cuda(cudaError_t result, const char* msg)
{
  if (result != cudaSuccess) {
    fatal_error(std::string{"CUDA error during "} + msg + ": " +
      cudaGetErrorString(result));
  }
}

inline void initialize_state()
{
  auto& s = state();
  if (s.initialized)
    return;

  int device_count = 0;
  auto status = cudaGetDeviceCount(&device_count);
  s.available = (status == cudaSuccess && device_count > 0);
  if (status != cudaSuccess && status != cudaErrorNoDevice) {
    // Clear the error so subsequent CUDA calls do not fail immediately
    cudaGetLastError();
  }
  s.initialized = true;
}

struct DeviceBuffers {
  double* sigma_t {nullptr};
  float* source {nullptr};
  float* angular_flux {nullptr};
  float* delta_psi {nullptr};
  int capacity {0};

  void resize(int n)
  {
    if (n <= capacity)
      return;

    free();
    check_cuda(cudaMalloc(&sigma_t, n * sizeof(double)), "cudaMalloc(sigma_t)");
    check_cuda(cudaMalloc(&source, n * sizeof(float)), "cudaMalloc(source)");
    check_cuda(cudaMalloc(&angular_flux, n * sizeof(float)),
      "cudaMalloc(angular_flux)");
    check_cuda(cudaMalloc(&delta_psi, n * sizeof(float)), "cudaMalloc(delta_psi)");
    capacity = n;
  }

  void free()
  {
    if (sigma_t) {
      cudaFree(sigma_t);
      sigma_t = nullptr;
    }
    if (source) {
      cudaFree(source);
      source = nullptr;
    }
    if (angular_flux) {
      cudaFree(angular_flux);
      angular_flux = nullptr;
    }
    if (delta_psi) {
      cudaFree(delta_psi);
      delta_psi = nullptr;
    }
    capacity = 0;
  }

  ~DeviceBuffers() { free(); }
};

thread_local DeviceBuffers buffers;

__global__ void flat_source_kernel(int n, double distance, const double* sigma_t,
  const float* source, float* angular_flux, float* delta_psi)
{
  int g = blockIdx.x * blockDim.x + threadIdx.x;
  if (g >= n)
    return;

  float tau = static_cast<float>(sigma_t[g] * distance);
  float exponential = cjosey_exponential(tau);
  float new_delta_psi = (angular_flux[g] - source[g]) * exponential;
  delta_psi[g] = new_delta_psi;
  angular_flux[g] -= new_delta_psi;
}

} // namespace

bool available()
{
  initialize_state();
  return state().available;
}

void attenuate_flat_source(int negroups, double distance, const double* sigma_t,
  const float* source, float* angular_flux, float* delta_psi)
{
  initialize_state();
  if (!state().available) {
    fatal_error("CUDA random ray attenuation requested but no CUDA-capable device "
                "is available.");
  }

  if (negroups <= 0)
    return;

  buffers.resize(negroups);

  check_cuda(cudaMemcpy(buffers.sigma_t, sigma_t, negroups * sizeof(double),
              cudaMemcpyHostToDevice),
    "cudaMemcpy(sigma_t)");
  check_cuda(cudaMemcpy(buffers.source, source, negroups * sizeof(float),
              cudaMemcpyHostToDevice),
    "cudaMemcpy(source)");
  check_cuda(cudaMemcpy(buffers.angular_flux, angular_flux,
              negroups * sizeof(float), cudaMemcpyHostToDevice),
    "cudaMemcpy(angular_flux)");

  constexpr int block_size = 128;
  int grid_size = (negroups + block_size - 1) / block_size;
  flat_source_kernel<<<grid_size, block_size>>>(negroups, distance, buffers.sigma_t,
    buffers.source, buffers.angular_flux, buffers.delta_psi);
  check_cuda(cudaGetLastError(), "flat_source_kernel launch");

  check_cuda(cudaMemcpy(delta_psi, buffers.delta_psi, negroups * sizeof(float),
              cudaMemcpyDeviceToHost),
    "cudaMemcpy(delta_psi)");
  check_cuda(cudaMemcpy(angular_flux, buffers.angular_flux,
              negroups * sizeof(float), cudaMemcpyDeviceToHost),
    "cudaMemcpy(angular_flux back)");
}

} // namespace detail
} // namespace cuda
} // namespace random_ray
} // namespace openmc

#endif // OPENMC_USE_CUDA
