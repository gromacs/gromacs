/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2010, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \libinternal \file
 *  \brief Declare functions for detection and initialization for GPU devices.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_GPU_UTILS_H
#define GMX_GPU_UTILS_GPU_UTILS_H

#include <cstdio>

#include <string>
#include <vector>

#include "gromacs/gpu_utils/gpu_macros.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_device_info_t;
struct gmx_gpu_info_t;

namespace gmx
{
class MDLogger;
}

//! Enum which is only used to describe transfer calls at the moment
enum class GpuApiCallBehavior
{
    Sync,
    Async
};

//! Types of actions associated to waiting or checking the completion of GPU tasks
enum class GpuTaskCompletion
{
    Wait, /*<< Issue a blocking wait for the task */
    Check /*<< Only check whether the task has completed */
};

/*! \brief Return whether GPUs can be detected
 *
 * Returns true when this is a build of \Gromacs configured to support
 * GPU usage, GPU detection is not disabled by an environment variable
 * and a valid device driver, ICD, and/or runtime was detected.
 * Does not throw. */
bool canPerformGpuDetection();

/*! \brief Return whether GPU detection is functioning correctly
 *
 * Returns true when this is a build of \Gromacs configured to support
 * GPU usage, and a valid device driver, ICD, and/or runtime was detected.
 *
 * This function is not intended to be called from build
 * configurations that do not support GPUs, and there will be no
 * descriptive message in that case.
 *
 * \param[out] errorMessage  When returning false on a build configured with
 *                           GPU support and non-nullptr was passed,
 *                           the string contains a descriptive message about
 *                           why GPUs cannot be detected.
 *
 * Does not throw. */
GPU_FUNC_QUALIFIER
bool isGpuDetectionFunctional(std::string* GPU_FUNC_ARGUMENT(errorMessage))
        GPU_FUNC_TERM_WITH_RETURN(false);

/*! \brief Find all GPUs in the system.
 *
 *  Will detect every GPU supported by the device driver in use.
 *  Must only be called if canPerformGpuDetection() has returned true.
 *  This routine also checks for the compatibility of each and fill the
 *  gpu_info->gpu_dev array with the required information on each the
 *  device: ID, device properties, status.
 *
 *  Note that this function leaves the GPU runtime API error state clean;
 *  this is implemented ATM in the CUDA flavor.
 *  TODO: check if errors do propagate in OpenCL as they do in CUDA and
 *  whether there is a mechanism to "clear" them.
 *
 *  \param[in] gpu_info    pointer to structure holding GPU information.
 *
 *  \throws                InternalError if a GPU API returns an unexpected failure (because
 *                         the call to canDetectGpus() should always prevent this occuring)
 */
GPU_FUNC_QUALIFIER
void findGpus(gmx_gpu_info_t* GPU_FUNC_ARGUMENT(gpu_info)) GPU_FUNC_TERM;

/*! \brief Return a container of the detected GPUs that are compatible.
 *
 * This function filters the result of the detection for compatible
 * GPUs, based on the previously run compatibility tests.
 *
 * \param[in]     gpu_info    Information detected about GPUs, including compatibility.
 * \return                    vector of IDs of GPUs already recorded as compatible */
std::vector<int> getCompatibleGpus(const gmx_gpu_info_t& gpu_info);

/*! \brief Return a string describing how compatible the GPU with given \c index is.
 *
 * \param[in]   gpu_info    Information about detected GPUs
 * \param[in]   index       index of GPU to ask about
 * \returns                 A null-terminated C string describing the compatibility status, useful for error messages.
 */
const char* getGpuCompatibilityDescription(const gmx_gpu_info_t& gpu_info, int index);

/*! \brief Frees the gpu_dev and dev_use array fields of \p gpu_info.
 *
 * \param[in]    gpu_info    pointer to structure holding GPU information
 */
void free_gpu_info(const gmx_gpu_info_t* gpu_info);

/*! \brief Initializes the GPU described by \c deviceInfo.
 *
 * TODO Doxygen complains about these - probably a Doxygen bug, since
 * the patterns here are the same as elsewhere in this header.
 *
 * \param[in]    deviceInfo   device info of the GPU to initialize
 *
 * Issues a fatal error for any critical errors that occur during
 * initialization.
 */
GPU_FUNC_QUALIFIER
void init_gpu(const gmx_device_info_t* GPU_FUNC_ARGUMENT(deviceInfo)) GPU_FUNC_TERM;

/*! \brief Frees up the CUDA GPU used by the active context at the time of calling.
 *
 * If \c deviceInfo is nullptr, then it is understood that no device
 * was selected so no context is active to be freed. Otherwise, the
 * context is explicitly destroyed and therefore all data uploaded to
 * the GPU is lost. This must only be called when none of this data is
 * required anymore, because subsequent attempts to free memory
 * associated with the context will otherwise fail.
 *
 * Calls gmx_warning upon errors.
 *
 * \param[in]  deviceInfo   device info of the GPU to clean up for
 *
 * \returns                 true if no error occurs during the freeing.
 */
CUDA_FUNC_QUALIFIER
void free_gpu(const gmx_device_info_t* CUDA_FUNC_ARGUMENT(deviceInfo)) CUDA_FUNC_TERM;

/*! \brief Return a pointer to the device info for \c deviceId
 *
 * \param[in] gpu_info      GPU info of all detected devices in the system.
 * \param[in] deviceId      ID for the GPU device requested.
 *
 * \returns                 Pointer to the device info for \c deviceId.
 */
GPU_FUNC_QUALIFIER
gmx_device_info_t* getDeviceInfo(const gmx_gpu_info_t& GPU_FUNC_ARGUMENT(gpu_info),
                                 int GPU_FUNC_ARGUMENT(deviceId)) GPU_FUNC_TERM_WITH_RETURN(nullptr);

/*! \brief Returns the device ID of the CUDA GPU currently in use.
 *
 * The GPU used is the one that is active at the time of the call in the active context.
 *
 * \returns                 device ID of the GPU in use at the time of the call
 */
CUDA_FUNC_QUALIFIER
int get_current_cuda_gpu_device_id() CUDA_FUNC_TERM_WITH_RETURN(-1);

/*! \brief Formats and returns a device information string for a given GPU.
 *
 * Given an index *directly* into the array of available GPUs (gpu_dev)
 * returns a formatted info string for the respective GPU which includes
 * ID, name, compute capability, and detection status.
 *
 * \param[out]  s           pointer to output string (has to be allocated externally)
 * \param[in]   gpu_info    Information about detected GPUs
 * \param[in]   index       an index *directly* into the array of available GPUs
 */
GPU_FUNC_QUALIFIER
void get_gpu_device_info_string(char*                 GPU_FUNC_ARGUMENT(s),
                                const gmx_gpu_info_t& GPU_FUNC_ARGUMENT(gpu_info),
                                int                   GPU_FUNC_ARGUMENT(index)) GPU_FUNC_TERM;


/*! \brief Returns the size of the gpu_dev_info struct.
 *
 * The size of gpu_dev_info can be used for allocation and communication.
 *
 * \returns                 size in bytes of gpu_dev_info
 */
GPU_FUNC_QUALIFIER
size_t sizeof_gpu_dev_info() GPU_FUNC_TERM_WITH_RETURN(0);

//! Get status of device with specified index
int gpu_info_get_stat(const gmx_gpu_info_t& info, int index);

/*! \brief Check if GROMACS has been built with GPU support.
 *
 * \param[in] error Pointer to error string or nullptr.
 * \todo Move this to NB module once it exists.
 */
bool buildSupportsNonbondedOnGpu(std::string* error);

/*! \brief Starts the GPU profiler if mdrun is being profiled.
 *
 *  When a profiler run is in progress (based on the presence of the NVPROF_ID
 *  env. var.), the profiler is started to begin collecting data during the
 *  rest of the run (or until stopGpuProfiler is called).
 *
 *  Note that this is implemented only for the CUDA API.
 */
CUDA_FUNC_QUALIFIER
void startGpuProfiler() CUDA_FUNC_TERM;


/*! \brief Resets the GPU profiler if mdrun is being profiled.
 *
 * When a profiler run is in progress (based on the presence of the NVPROF_ID
 * env. var.), the profiler data is restet in order to eliminate the data collected
 * from the preceding part fo the run.
 *
 * This function should typically be called at the mdrun counter reset time.
 *
 * Note that this is implemented only for the CUDA API.
 */
CUDA_FUNC_QUALIFIER
void resetGpuProfiler() CUDA_FUNC_TERM;


/*! \brief Stops the CUDA profiler if mdrun is being profiled.
 *
 *  This function can be called at cleanup when skipping recording
 *  recording subsequent API calls from being traces/profiled is desired,
 *  e.g. before uninitialization.
 *
 *  Note that this is implemented only for the CUDA API.
 */
CUDA_FUNC_QUALIFIER
void stopGpuProfiler() CUDA_FUNC_TERM;

//! Tells whether the host buffer was pinned for non-blocking transfers. Only implemented for CUDA.
CUDA_FUNC_QUALIFIER
bool isHostMemoryPinned(const void* CUDA_FUNC_ARGUMENT(h_ptr)) CUDA_FUNC_TERM_WITH_RETURN(false);

/*! \brief Enable peer access between GPUs where supported
 * \param[in] gpuIdsToUse   List of GPU IDs in use
 * \param[in] mdlog         Logger object
 */
CUDA_FUNC_QUALIFIER
void setupGpuDevicePeerAccess(const std::vector<int>& CUDA_FUNC_ARGUMENT(gpuIdsToUse),
                              const gmx::MDLogger&    CUDA_FUNC_ARGUMENT(mdlog)) CUDA_FUNC_TERM;

#endif
