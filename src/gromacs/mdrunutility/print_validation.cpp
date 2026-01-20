/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 *
 * \brief Defines functions that write a message about the validation state of features
 *
 * \ingroup module_mdrunutility
 */

#include "gmxpre.h"

#include "print_validation.h"

#include "config.h"

#include "gromacs/hardware/device_information.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/message_string_collector.h"

namespace gmx
{

void logValidationMessages(const MDLogger&            mdlog,
                           const SimulationWorkload&  simulationWorkload,
                           const bool                 haveFillerParticlesInLocalState,
                           const bool                 useH5mdOutputFile,
                           const bool                 useModularSimulator,
                           const IntegrationAlgorithm integrationAlgorithm,
                           const DeviceInformation*   deviceInfo)
{
    // NOTE: when adding/removing features here, update docs/user-guide/validation.rst
    MessageStringCollector experimentalFeatures;
    experimentalFeatures.startContext(
            "You are using the following feature(s) which are experimental. We do not "
            "recommend this for production simulations. Check your results. If you find "
            "issues, report them on the GROMACS forum or at "
            "https://gitlab.com/gromacs/gromacs/-/issues");
    experimentalFeatures.appendIf(simulationWorkload.useMdGpuGraph && GMX_GPU_CUDA, "CUDA graphs");
    experimentalFeatures.appendIf(simulationWorkload.useMdGpuGraph && GMX_GPU_SYCL, "SYCL graphs");
    experimentalFeatures.appendIf(simulationWorkload.useNvshmem, "NVSHMEM");
    experimentalFeatures.appendIf(haveFillerParticlesInLocalState, "direct halo communication");
    experimentalFeatures.appendIf(GMX_ENABLE_NBNXM_CPU_VECTORIZATION,
                                  "non-bonded kernel vectorization");
    experimentalFeatures.appendIf(GMX_SYCL_ENABLE_HANDLER_FREE_SUBMISSION,
                                  "oneAPI experimental handler-free submission");
    experimentalFeatures.appendIf(GMX_GPU_FFT_VKFFT && deviceInfo != nullptr
                                          && deviceInfo->deviceVendor != DeviceVendor::Amd
                                          && deviceInfo->deviceVendor != DeviceVendor::Apple,
                                  "VkFFT library on a non-AMD and non-Apple platform");
    experimentalFeatures.appendIf(GMX_GPU_FFT_ONEMATH, "OneMath FFT library on any platform");
    experimentalFeatures.appendIf(GMX_GPU_FFT_BBFFT, "double-batched FFT library");
    experimentalFeatures.appendIf(GMX_GPU_HIP && simulationWorkload.useGpuPmeDecomposition,
                                  "PME GPU decomposition with the AMD HIP GPU backend");
    experimentalFeatures.appendIf(useH5mdOutputFile, "H5md trajectory output");
    experimentalFeatures.finishContext();
    if (!experimentalFeatures.isEmpty())
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText(experimentalFeatures.toString());
    }

    // NOTE: when adding/removing features here, update docs/user-guide/validation.rst
    MessageStringCollector validationPendingFeatures;
    validationPendingFeatures.startContext(
            "You are using the following feature(s) which are partially, but not fully "
            "validated. Please check your results. If you find issues, report them on the "
            "GROMACS forum or at https://gitlab.com/gromacs/gromacs/-/issues");
    validationPendingFeatures.appendIf(
            useModularSimulator && integrationAlgorithm != IntegrationAlgorithm::VV,
            "modular simulator with an integrator other than VV");
    validationPendingFeatures.appendIf(GMX_SYCL_DPCPP && deviceInfo != nullptr
                                               && deviceInfo->deviceVendor != DeviceVendor::Intel,
                                       "oneAPI SYCL GPU backend on a non-Intel platform");
    validationPendingFeatures.appendIf(
            GMX_SYCL_ACPP && deviceInfo != nullptr && deviceInfo->deviceVendor != DeviceVendor::Amd,
            "AdaptiveCpp SYCL GPU backend on a non-AMD platform");
    validationPendingFeatures.appendIf(GMX_GPU_HIP && simulationWorkload.useGpuPmeDecomposition,
                                       "PME GPU decomposition with the AMD HIP GPU backend");
    validationPendingFeatures.appendIf(simulationWorkload.useGpuNonbondedFE,
                                       "Non-bonded free-energy calculations on the GPU");
    // TODO MdModules doesn't have a way to report which
    // ForceProviders are active (nor have they even been built at
    // this point of mdrunner())
    validationPendingFeatures.appendIf(GMX_USE_EXT_FMM, "an external FMM library");
    // TODO report if NNPot force provider is active
    // TODO report if Colvars plug-in is active
    // TODO report if Plumed plug-in is active
    validationPendingFeatures.finishContext();
    if (!validationPendingFeatures.isEmpty())
    {
        GMX_LOG(mdlog.warning).asParagraph().appendText(validationPendingFeatures.toString());
    }
}

} // namespace gmx
