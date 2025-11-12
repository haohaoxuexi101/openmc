#ifndef OPENMC_RANDOM_RAY_RANDOM_RAY_CUDA_H
#define OPENMC_RANDOM_RAY_RANDOM_RAY_CUDA_H

#include <cstddef>

namespace openmc {
namespace random_ray {
namespace cuda {

bool available();

void attenuate_flat_source(int negroups, double distance, const double* sigma_t,
  const float* source, float* angular_flux, float* delta_psi);

} // namespace cuda
} // namespace random_ray
} // namespace openmc

#endif // OPENMC_RANDOM_RAY_RANDOM_RAY_CUDA_H
