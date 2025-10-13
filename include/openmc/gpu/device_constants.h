#ifndef OPENMC_GPU_DEVICE_CONSTANTS_H
#define OPENMC_GPU_DEVICE_CONSTANTS_H

namespace openmc {
namespace cuda {

constexpr int PARTICLE_NEUTRON = 0;
constexpr int PARTICLE_PHOTON = 1;
constexpr int PARTICLE_ELECTRON = 2;
constexpr int PARTICLE_POSITRON = 3;

} // namespace cuda
} // namespace openmc

#endif // OPENMC_GPU_DEVICE_CONSTANTS_H
