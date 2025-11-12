#include "openmc/random_ray/random_ray_cuda.h"

#include "openmc/error.h"

namespace openmc {
namespace random_ray {
namespace cuda {

#ifdef OPENMC_USE_CUDA
namespace detail {
bool available();
void attenuate_flat_source(int negroups, double distance, const double* sigma_t,
  const float* source, float* angular_flux, float* delta_psi);
} // namespace detail
#endif

bool available()
{
#ifdef OPENMC_USE_CUDA
  return detail::available();
#else
  return false;
#endif
}

void attenuate_flat_source(int negroups, double distance, const double* sigma_t,
  const float* source, float* angular_flux, float* delta_psi)
{
#ifdef OPENMC_USE_CUDA
  detail::attenuate_flat_source(negroups, distance, sigma_t, source, angular_flux,
    delta_psi);
#else
  fatal_error("Random ray CUDA attenuation requested but CUDA support was not "
              "enabled at build time.");
#endif
}

} // namespace cuda
} // namespace random_ray
} // namespace openmc
