#include "openmc/gpu/driver.h"

#ifdef OPENMC_USE_CUDA

#include <cuda_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/gpu/device_constants.h"
#include "openmc/particle.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"

namespace openmc {
namespace cuda {
namespace {

#define CUDA_CHECK(call)                                                       \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fatal_error(std::string{"CUDA error: "} + cudaGetErrorString(err));     \
    }                                                                          \
  } while (0)

// Simple RAII wrapper for device memory allocations.
template<class T>
class DeviceBuffer {
public:
  DeviceBuffer(size_t count) : size_{count}
  {
    if (size_ == 0) {
      data_ = nullptr;
      return;
    }
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void**>(&data_), size_ * sizeof(T)));
  }

  ~DeviceBuffer()
  {
    if (data_)
      cudaFree(data_);
  }

  T* data() { return data_; }
  const T* data() const { return data_; }
  size_t size() const { return size_; }

private:
  T* data_ {nullptr};
  size_t size_ {0};
};

__device__ __forceinline__ double prn_device(uint64_t& seed)
{
  constexpr uint64_t prn_mult {6364136223846793005ULL};
  constexpr uint64_t prn_add {1442695040888963407ULL};

  seed = prn_mult * seed + prn_add;
  uint64_t word = ((seed >> ((seed >> 59u) + 5u)) ^ seed) *
    12605985483714917081ull;
  uint64_t result = (word >> 43u) ^ word;
  return ldexp(static_cast<double>(result), -64);
}

__global__ void sample_distance_kernel(int n, const int* types,
  const double* macro_totals, const double* boundary_distances, uint64_t* seeds,
  double* collision_distances, double* travel_distances, int* collision_flags,
  double inf)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n)
    return;

  int type = types[tid];
  double macro = macro_totals[tid];
  double boundary = boundary_distances[tid];
  uint64_t seed = seeds[tid];

  double collision_distance;
  if (type == PARTICLE_ELECTRON || type == PARTICLE_POSITRON) {
    collision_distance = 0.0;
  } else if (macro <= 0.0) {
    collision_distance = inf;
  } else {
    double xi = prn_device(seed);
    collision_distance = -log(xi) / macro;
  }

  double distance = fmin(boundary, collision_distance);

  collision_distances[tid] = collision_distance;
  travel_distances[tid] = distance;
  collision_flags[tid] = (collision_distance <= boundary) ? 1 : 0;
  seeds[tid] = seed;
}

void sample_distances_gpu_impl(int n, const int* types,
  const double* macro_totals, const double* boundary_distances, uint64_t* seeds,
  double* collision_distances, double* travel_distances, int* collision_flags)
{
  if (n == 0)
    return;

  DeviceBuffer<int> d_types(n);
  DeviceBuffer<double> d_macro_totals(n);
  DeviceBuffer<double> d_boundary_distances(n);
  DeviceBuffer<uint64_t> d_seeds(n);
  DeviceBuffer<double> d_collision_distances(n);
  DeviceBuffer<double> d_travel_distances(n);
  DeviceBuffer<int> d_collision_flags(n);

  CUDA_CHECK(cudaMemcpy(d_types.data(), types, n * sizeof(int),
    cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_macro_totals.data(), macro_totals, n * sizeof(double),
    cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_boundary_distances.data(), boundary_distances,
    n * sizeof(double), cudaMemcpyHostToDevice));
  CUDA_CHECK(cudaMemcpy(d_seeds.data(), seeds, n * sizeof(uint64_t),
    cudaMemcpyHostToDevice));

  const double inf = std::numeric_limits<double>::infinity();
  int block_size = std::max(1, settings::cuda_block_size);
  int grid_size = (n + block_size - 1) / block_size;

  sample_distance_kernel<<<grid_size, block_size>>>(n, d_types.data(),
    d_macro_totals.data(), d_boundary_distances.data(), d_seeds.data(),
    d_collision_distances.data(), d_travel_distances.data(),
    d_collision_flags.data(), inf);
  CUDA_CHECK(cudaGetLastError());
  CUDA_CHECK(cudaDeviceSynchronize());

  CUDA_CHECK(cudaMemcpy(seeds, d_seeds.data(), n * sizeof(uint64_t),
    cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(collision_distances, d_collision_distances.data(),
    n * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(travel_distances, d_travel_distances.data(),
    n * sizeof(double), cudaMemcpyDeviceToHost));
  CUDA_CHECK(cudaMemcpy(collision_flags, d_collision_flags.data(),
    n * sizeof(int), cudaMemcpyDeviceToHost));
}

bool runtime_initialized {false};

} // namespace

void initialize_runtime()
{
  if (runtime_initialized)
    return;

  CUDA_CHECK(cudaFree(nullptr));
  int device_count = 0;
  CUDA_CHECK(cudaGetDeviceCount(&device_count));
  if (device_count == 0) {
    fatal_error("OPENMC_USE_CUDA is enabled, but no CUDA-capable device was "
                "found.");
  }
  runtime_initialized = true;
}

void finalize_runtime()
{
  if (!runtime_initialized)
    return;

  CUDA_CHECK(cudaDeviceSynchronize());
  runtime_initialized = false;
}

bool process_advance_event_queue(SharedArray<EventQueueItem>& queue)
{
  if (!settings::cuda_enabled || !settings::cuda_accelerate_advance)
    return false;

  initialize_runtime();

  const size_t total = queue.size();
  if (total == 0) {
    queue.resize(0);
    return true;
  }

  static_assert(static_cast<int>(ParticleType::neutron) == PARTICLE_NEUTRON,
    "ParticleType::neutron enum mismatch with CUDA constants");
  static_assert(static_cast<int>(ParticleType::photon) == PARTICLE_PHOTON,
    "ParticleType::photon enum mismatch with CUDA constants");
  static_assert(static_cast<int>(ParticleType::electron) == PARTICLE_ELECTRON,
    "ParticleType::electron enum mismatch with CUDA constants");
  static_assert(static_cast<int>(ParticleType::positron) == PARTICLE_POSITRON,
    "ParticleType::positron enum mismatch with CUDA constants");

  const size_t max_batch =
    static_cast<size_t>(std::max(1, settings::cuda_max_batch));
  const size_t chunk_capacity = std::min(max_batch, total);

  std::vector<int> types(chunk_capacity);
  std::vector<double> macro_totals(chunk_capacity);
  std::vector<double> boundary_distances(chunk_capacity);
  std::vector<uint64_t> seeds(chunk_capacity);
  std::vector<double> collision_distances(chunk_capacity);
  std::vector<double> travel_distances(chunk_capacity);
  std::vector<int> collision_flags(chunk_capacity);
  std::vector<int64_t> indices(chunk_capacity);

  for (size_t offset = 0; offset < total; offset += chunk_capacity) {
    size_t count = std::min(chunk_capacity, total - offset);

    for (size_t i = 0; i < count; ++i) {
      auto& item = queue[offset + i];
      Particle& p = simulation::particles[item.idx];
      p.boundary() = distance_to_boundary(p);

      types[i] = static_cast<int>(p.type());
      macro_totals[i] = p.macro_xs().total;
      boundary_distances[i] = p.boundary().distance;
      seeds[i] = *p.current_seed();
      indices[i] = item.idx;
    }

    sample_distances_gpu_impl(static_cast<int>(count), types.data(),
      macro_totals.data(), boundary_distances.data(), seeds.data(),
      collision_distances.data(), travel_distances.data(),
      collision_flags.data());

    for (size_t i = 0; i < count; ++i) {
      Particle& p = simulation::particles[indices[i]];
      *p.current_seed() = seeds[i];
      p.collision_distance() = collision_distances[i];
      double distance = travel_distances[i];
      p.advance_along_distance(distance);

      if (!p.alive())
        continue;

      EventQueueItem next_item {p, indices[i]};
      if (collision_flags[i]) {
        simulation::collision_queue.thread_safe_append(next_item);
      } else {
        simulation::surface_crossing_queue.thread_safe_append(next_item);
      }
    }
  }

  queue.resize(0);
  return true;
}

} // namespace cuda
} // namespace openmc

#endif // OPENMC_USE_CUDA
