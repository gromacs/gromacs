/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020,2021, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief Implements SETTLE using GPU
 *
 * This file contains implementation for the data management of GPU version of SETTLE constraints algorithm.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "settle_gpu.h"

#include <assert.h>
#include <stdio.h>

#include <cmath>

#include <algorithm>

#include "gromacs/gpu_utils/devicebuffer.h"
#include "gromacs/gpu_utils/gputraits.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/settle_gpu_internal.h"
#include "gromacs/pbcutil/pbc.h"

namespace gmx
{

void SettleGpu::apply(const DeviceBuffer<Float3>& d_x,
                      DeviceBuffer<Float3>        d_xp,
                      const bool                  updateVelocities,
                      DeviceBuffer<Float3>        d_v,
                      const real                  invdt,
                      const bool                  computeVirial,
                      tensor                      virialScaled,
                      const PbcAiuc&              pbcAiuc)
{

    // Early exit if no settles
    if (numSettles_ == 0)
    {
        return;
    }

    if (computeVirial)
    {
        // Fill with zeros so the values can be reduced to it
        // Only 6 values are needed because virial is symmetrical
        clearDeviceBufferAsync(&d_virialScaled_, 0, 6, deviceStream_);
    }

    launchSettleGpuKernel(numSettles_,
                          d_atomIds_,
                          settleParameters_,
                          d_x,
                          d_xp,
                          updateVelocities,
                          d_v,
                          invdt,
                          computeVirial,
                          d_virialScaled_,
                          pbcAiuc,
                          deviceStream_);


    if (computeVirial)
    {
        copyFromDeviceBuffer(
                h_virialScaled_.data(), &d_virialScaled_, 0, 6, deviceStream_, GpuApiCallBehavior::Sync, nullptr);

        // Mapping [XX, XY, XZ, YY, YZ, ZZ] internal format to a tensor object
        virialScaled[XX][XX] += h_virialScaled_[0];
        virialScaled[XX][YY] += h_virialScaled_[1];
        virialScaled[XX][ZZ] += h_virialScaled_[2];

        virialScaled[YY][XX] += h_virialScaled_[1];
        virialScaled[YY][YY] += h_virialScaled_[3];
        virialScaled[YY][ZZ] += h_virialScaled_[4];

        virialScaled[ZZ][XX] += h_virialScaled_[2];
        virialScaled[ZZ][YY] += h_virialScaled_[4];
        virialScaled[ZZ][ZZ] += h_virialScaled_[5];
    }
}

SettleGpu::SettleGpu(const gmx_mtop_t& mtop, const DeviceContext& deviceContext, const DeviceStream& deviceStream) :
    deviceContext_(deviceContext), deviceStream_(deviceStream)
{
    static_assert(sizeof(real) == sizeof(float),
                  "Real numbers should be in single precision in GPU code.");

    // This is to prevent the assertion failure for the systems without water
    int totalSettles = 0;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        totalSettles += iatoms.size() / nral1;
    }
    if (totalSettles == 0)
    {
        return;
    }

    // TODO This should be lifted to a separate subroutine that gets the values of Oxygen and
    // Hydrogen masses, checks if they are consistent across the topology and if there is no more
    // than two values for each mass if the free energy perturbation is enabled. In later case,
    // masses may need to be updated on a regular basis (i.e. in set(...) method).
    // TODO Do the checks for FEP
    real mO = -1.0;
    real mH = -1.0;

    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int        nral1           = 1 + NRAL(F_SETTLE);
        InteractionList  interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        std::vector<int> iatoms          = interactionList.iatoms;
        for (unsigned i = 0; i < iatoms.size() / nral1; i++)
        {
            WaterMolecule settler;
            settler.ow1 = iatoms[i * nral1 + 1]; // Oxygen index
            settler.hw2 = iatoms[i * nral1 + 2]; // First hydrogen index
            settler.hw3 = iatoms[i * nral1 + 3]; // Second hydrogen index
            t_atom ow1  = mtop.moltype[mt].atoms.atom[settler.ow1];
            t_atom hw2  = mtop.moltype[mt].atoms.atom[settler.hw2];
            t_atom hw3  = mtop.moltype[mt].atoms.atom[settler.hw3];

            if (mO < 0)
            {
                mO = ow1.m;
            }
            if (mH < 0)
            {
                mH = hw2.m;
            }
            GMX_RELEASE_ASSERT(mO == ow1.m,
                               "Topology has different values for oxygen mass. Should be identical "
                               "in order to use SETTLE.");
            GMX_RELEASE_ASSERT(hw2.m == hw3.m && hw2.m == mH,
                               "Topology has different values for hydrogen mass. Should be "
                               "identical in order to use SETTLE.");
        }
    }

    GMX_RELEASE_ASSERT(mO > 0, "Could not find oxygen mass in the topology. Needed in SETTLE.");
    GMX_RELEASE_ASSERT(mH > 0, "Could not find hydrogen mass in the topology. Needed in SETTLE.");

    // TODO Very similar to SETTLE initialization on CPU. Should be lifted to a separate method
    // (one that gets dOH and dHH values and checks them for consistency)
    int settle_type = -1;
    for (unsigned mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const int       nral1           = 1 + NRAL(F_SETTLE);
        InteractionList interactionList = mtop.moltype[mt].ilist[F_SETTLE];
        for (int i = 0; i < interactionList.size(); i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = interactionList.iatoms[i];
            }
            else if (interactionList.iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block "
                          "of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong "
                          "approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }

    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    real dOH = mtop.ffparams.iparams[settle_type].settle.doh;
    real dHH = mtop.ffparams.iparams[settle_type].settle.dhh;

    settleParameters_ = settleParameters(mO, mH, 1.0 / mO, 1.0 / mH, dOH, dHH);

    allocateDeviceBuffer(&d_virialScaled_, 6, deviceContext_);
    h_virialScaled_.resize(6);
}

SettleGpu::~SettleGpu()
{
    // Early exit if there is no settles
    if (numSettles_ == 0)
    {
        return;
    }
    freeDeviceBuffer(&d_virialScaled_);
    if (numAtomIdsAlloc_ > 0)
    {
        freeDeviceBuffer(&d_atomIds_);
    }
}

void SettleGpu::set(const InteractionDefinitions& idef)
{
    const int              nral1     = 1 + NRAL(F_SETTLE);
    const InteractionList& il_settle = idef.il[F_SETTLE];
    ArrayRef<const int>    iatoms    = il_settle.iatoms;
    numSettles_                      = il_settle.size() / nral1;

    reallocateDeviceBuffer(&d_atomIds_, numSettles_, &numAtomIds_, &numAtomIdsAlloc_, deviceContext_);
    h_atomIds_.resize(numSettles_);
    for (int i = 0; i < numSettles_; i++)
    {
        WaterMolecule settler;
        settler.ow1   = iatoms[i * nral1 + 1]; // Oxygen index
        settler.hw2   = iatoms[i * nral1 + 2]; // First hydrogen index
        settler.hw3   = iatoms[i * nral1 + 3]; // Second hydrogen index
        h_atomIds_[i] = settler;
    }
    copyToDeviceBuffer(
            &d_atomIds_, h_atomIds_.data(), 0, numSettles_, deviceStream_, GpuApiCallBehavior::Sync, nullptr);
}

} // namespace gmx
