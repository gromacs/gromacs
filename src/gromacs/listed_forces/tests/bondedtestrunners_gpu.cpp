/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief GPU implementation of the bonded force test runner.
 *
 * Sets up minimal idef and NBAtomDataGpu, runs the bonded kernel on the device,
 * and copies forces and energy back for comparison with the CPU reference.
 *
 * \ingroup module_listed_forces
 */
#include "gmxpre.h"

#include "config.h"

#include <algorithm>
#include <array>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/capabilities.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/vectypes.h"

#include "testutils/test_device.h"

#include "bondedtestrunners.h"

#if GMX_GPU
#    include "gromacs/gpu_utils/devicebuffer.h"
#    include "gromacs/gpu_utils/gputraits.h"
#    include "gromacs/mdtypes/enerdata.h"
#    include "gromacs/mdtypes/simulation_workload.h"
#    include "gromacs/nbnxm/gpu_types_common.h"
#    include "gromacs/pbcutil/ishift.h"
#    include "gromacs/topology/forcefieldparameters.h"
#    include "gromacs/topology/idef.h"
#    include "gromacs/topology/ifunc.h"
#endif

namespace gmx
{
namespace test
{

bool BondedDeviceTestRunner::supportsFlavor(BondedKernelFlavor flavor, const iListInput& input, real lambda) const
{
    if constexpr (!GpuConfigurationCapabilities::Bonded)
    {
        GMX_UNUSED_VALUE(flavor);
        GMX_UNUSED_VALUE(input);
        GMX_UNUSED_VALUE(lambda);
        return false;
    }
    // GPU bonded supports only the full flavor (forces + virial + energy), no FEP
    if (flavor != BondedKernelFlavor::ForcesAndVirialAndEnergy || input.fep)
    {
        return false;
    }
    if (!input.ftype.has_value())
    {
        return false;
    }
    return std::find(fTypesOnGpu.begin(), fTypesOnGpu.end(), input.ftype.value()) != fTypesOnGpu.end();
}

#if GMX_GPU && !GMX_GPU_OPENCL

void BondedDeviceTestRunner::run(const iListInput&           input,
                                 const PaddedVector<RVec>&   x,
                                 const t_pbc&                pbc,
                                 real                        lambda,
                                 const std::vector<t_iatom>& iatoms,
                                 BondedKernelFlavor          flavor,
                                 OutputQuantities*           output)
{
    GMX_UNUSED_VALUE(lambda);
    GMX_RELEASE_ASSERT(flavor == BondedKernelFlavor::ForcesAndVirialAndEnergy,
                       "GPU runner only supports ForcesAndVirialAndEnergy flavor");

    const DeviceContext& deviceContext = testDevice_.deviceContext();
    const DeviceStream&  deviceStream  = testDevice_.deviceStream();
    deviceContext.activate();

    const int  ftype    = static_cast<int>(input.ftype.value());
    const int  numAtoms = c_numAtomsBondedTest;
    const real scaleFac = 1.0f;

    // Build minimal force-field parameters (single type)
    gmx_ffparams_t ffparams;
    ffparams.functype = { input.ftype.value() };
    ffparams.iparams  = { input.iparams };

    // Build minimal interaction definitions
    InteractionDefinitions idef(ffparams);
    idef.il[ftype].iatoms                   = iatoms;
    idef.ilsort                             = ilsortNO_FE;
    idef.numNonperturbedInteractions[ftype] = static_cast<int>(iatoms.size());

    // Identity atom order (same as input coordinates)
    std::vector<int> nbnxnAtomOrder(numAtoms);
    for (int i = 0; i < numAtoms; i++)
    {
        nbnxnAtomOrder[i] = i;
    }

    // Host-side xq (x, y, z, q) for upload
    std::vector<real> h_xq;
    h_xq.reserve(numAtoms * 4);
    for (int i = 0; i < numAtoms; i++)
    {
        h_xq.emplace_back(x[i][XX]);
        h_xq.emplace_back(x[i][YY]);
        h_xq.emplace_back(x[i][ZZ]);
        h_xq.emplace_back(c_bondedTestCharges[i]);
    }

    // Allocate device buffers for NBAtomDataGpu (only xq, f, fShift used by bonded)
    NBAtomDataGpu nbAtomDataGpu = {};
    nbAtomDataGpu.numAtoms      = numAtoms;
    nbAtomDataGpu.numAtomsLocal = numAtoms;
    nbAtomDataGpu.numAtomsAlloc = numAtoms;

    allocateDeviceBuffer(&nbAtomDataGpu.xq, numAtoms, deviceContext);
    allocateDeviceBuffer(&nbAtomDataGpu.f, numAtoms, deviceContext);
    allocateDeviceBuffer(&nbAtomDataGpu.fShift, c_numShiftVectors, deviceContext);

    copyToDeviceBuffer(&nbAtomDataGpu.xq,
                       reinterpret_cast<Float4*>(h_xq.data()),
                       0,
                       numAtoms,
                       deviceStream,
                       GpuApiCallBehavior::Sync,
                       nullptr);
    clearDeviceBufferAsync(&nbAtomDataGpu.f, 0, numAtoms, deviceStream);
    clearDeviceBufferAsync(&nbAtomDataGpu.fShift, 0, c_numShiftVectors, deviceStream);

    // ListedForcesGpu (wallcycle can be nullptr – start/stop check for null)
    ListedForcesGpu listedForcesGpu(ffparams, scaleFac, 1, deviceContext, deviceStream, nullptr);

    listedForcesGpu.updateHaveInteractions(idef);
    listedForcesGpu.updateInteractionListsAndDeviceBuffers(nbnxnAtomOrder, idef, &nbAtomDataGpu);

    if (!listedForcesGpu.haveInteractions())
    {
        freeDeviceBuffer(&nbAtomDataGpu.xq);
        freeDeviceBuffer(&nbAtomDataGpu.f);
        freeDeviceBuffer(&nbAtomDataGpu.fShift);
        FAIL() << "ListedForcesGpu reported no interactions for type " << ftype;
    }

    listedForcesGpu.setPbc(pbc.pbcType, pbc.box, false);

    StepWorkload stepWork;
    stepWork.computeVirial = true;
    stepWork.computeEnergy = true;
    listedForcesGpu.launchKernel(stepWork);
    listedForcesGpu.launchEnergyTransfer();

    gmx_enerdata_t enerd(1, nullptr);
    listedForcesGpu.waitAccumulateEnergyTerms(&enerd);

    // Copy forces and shift forces back
    std::vector<Float3> h_f(numAtoms);
    std::vector<Float3> h_fShift(c_numShiftVectors);
    copyFromDeviceBuffer(
            h_f.data(), &nbAtomDataGpu.f, 0, numAtoms, deviceStream, GpuApiCallBehavior::Sync, nullptr);
    copyFromDeviceBuffer(h_fShift.data(),
                         &nbAtomDataGpu.fShift,
                         0,
                         c_numShiftVectors,
                         deviceStream,
                         GpuApiCallBehavior::Sync,
                         nullptr);

    // Fill test output (Float3 is RVec on HIP, uses [0],[1],[2])
    output->energy    = enerd.term[ftype];
    output->dvdlambda = 0.0;
    for (int i = 0; i < numAtoms; i++)
    {
        output->f[i][XX] = h_f[i][XX];
        output->f[i][YY] = h_f[i][YY];
        output->f[i][ZZ] = h_f[i][ZZ];
        output->f[i][3]  = 0.0;
    }
    for (int s = 0; s < c_numShiftVectors; s++)
    {
        output->fshift[s][XX] = h_fShift[s][XX];
        output->fshift[s][YY] = h_fShift[s][YY];
        output->fshift[s][ZZ] = h_fShift[s][ZZ];
    }

    freeDeviceBuffer(&nbAtomDataGpu.xq);
    freeDeviceBuffer(&nbAtomDataGpu.f);
    freeDeviceBuffer(&nbAtomDataGpu.fShift);
}

#else

void BondedDeviceTestRunner::run(const iListInput& /* input */,
                                 const PaddedVector<RVec>& /* x */,
                                 const t_pbc& /* pbc */,
                                 real /* lambda */,
                                 const std::vector<t_iatom>& /* iatoms */,
                                 BondedKernelFlavor /* flavor */,
                                 OutputQuantities* /* output */)
{
    GMX_UNUSED_VALUE(testDevice_);
    FAIL() << "Dummy bonded GPU runner was called instead of the real one.";
}

#endif

} // namespace test
} // namespace gmx
