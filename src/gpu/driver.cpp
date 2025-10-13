#include "openmc/gpu/driver.h"

#ifndef OPENMC_USE_CUDA

namespace openmc {
namespace cuda {

void initialize_runtime() {}

void finalize_runtime() {}

bool process_advance_event_queue(SharedArray<EventQueueItem>&)
{
  return false;
}

} // namespace cuda
} // namespace openmc

#endif // OPENMC_USE_CUDA
