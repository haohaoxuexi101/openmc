#ifndef OPENMC_GPU_DRIVER_H
#define OPENMC_GPU_DRIVER_H

#include "openmc/event.h"
#include "openmc/shared_array.h"

namespace openmc {
namespace cuda {

//! Initialize CUDA runtime resources when GPU acceleration is enabled.
void initialize_runtime();

//! Release CUDA runtime resources.
void finalize_runtime();

//! Attempt to process the advance-particle event queue on the GPU.
//!
//! \param queue Event queue containing particles that require advance events.
//! \return True if the queue was processed on the GPU, false otherwise.
bool process_advance_event_queue(SharedArray<EventQueueItem>& queue);

} // namespace cuda
} // namespace openmc

#endif // OPENMC_GPU_DRIVER_H
