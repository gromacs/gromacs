/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * \brief Defines the class to accumulate the data needed for the Leap-Frog integrator tests
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "leapfrogtestdata.h"

#include <cassert>
#include <cmath>

#include <algorithm>
#include <array>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

LeapFrogTestData::LeapFrogTestData(int        numAtoms,
                                   real       timestep,
                                   const rvec v0,
                                   const rvec f0,
                                   int        numTCoupleGroups,
                                   int        nstpcouple) :
    numAtoms_(numAtoms),
    timestep_(timestep),
    x0_(numAtoms),
    x_(numAtoms),
    xPrime_(numAtoms),
    v0_(numAtoms),
    v_(numAtoms),
    f_(numAtoms),
    inverseMasses_(numAtoms),
    inverseMassesPerDim_(numAtoms),
    kineticEnergyData_(std::vector<real>(numTCoupleGroups == 0 ? 1 : numTCoupleGroups, 0),
                       EnsembleTemperatureSetting::NotAvailable,
                       0.0,
                       false,
                       0.0,
                       1),
    numTCoupleGroups_(numTCoupleGroups)
{
    mdAtoms_.nr = numAtoms_;

    for (int i = 0; i < numAtoms_; i++)
    {
        // Typical PBC box size is tens of nanometers
        x_[i][XX] = (i % 21) * 1.0;
        x_[i][YY] = 6.5 + (i % 13) * (-1.0);
        x_[i][ZZ] = (i % 32) * (0.0);

        for (int d = 0; d < DIM; d++)
        {
            xPrime_[i][d] = 0.0;
            // Thermal velocity is ~1 nm/ps (|v0| = 1-2 nm/ps)
            v_[i][d] = v0[d];
            // TODO Check what value typical MD forces have (now ~ 1 kJ/mol/nm)
            f_[i][d] = f0[d];

            x0_[i][d] = x_[i][d];
            v0_[i][d] = v_[i][d];
        }
        // Atom masses are ~1-100 g/mol
        inverseMasses_[i] = 1.0 / (1.0 + i % 100);
        for (int d = 0; d < DIM; d++)
        {
            inverseMassesPerDim_[i][d] = inverseMasses_[i];
        }
    }
    mdAtoms_.invmass       = inverseMasses_;
    mdAtoms_.invMassPerDim = inverseMassesPerDim_;

    // Temperature coupling
    mdAtoms_.cTC.resize(numAtoms_);

    // To do temperature coupling at each step
    inputRecord_.nsttcouple = 1;

    if (numTCoupleGroups_ == 0)
    {
        inputRecord_.etc = TemperatureCoupling::No;
        for (int i = 0; i < numAtoms_; i++)
        {
            mdAtoms_.cTC[i] = 0;
        }
        t_grp_tcstat temperatureCouplingGroupData;
        temperatureCouplingGroupData.lambda = 1.0;
        kineticEnergyData_.tcstat[0]        = temperatureCouplingGroupData;
    }
    else
    {
        inputRecord_.etc = TemperatureCoupling::Yes;
        for (int i = 0; i < numAtoms_; i++)
        {
            mdAtoms_.cTC[i] = i % numTCoupleGroups_;
        }
        for (int i = 0; i < numTCoupleGroups_; i++)
        {
            real         tCoupleLambda = 1.0 - (i + 1.0) / 10.0;
            t_grp_tcstat temperatureCouplingGroupData;
            temperatureCouplingGroupData.lambda = tCoupleLambda;
            kineticEnergyData_.tcstat[i]        = temperatureCouplingGroupData;
        }
    }

    inputRecord_.eI      = IntegrationAlgorithm::MD;
    inputRecord_.delta_t = timestep_;

    state_.box[XX][XX] = 10.0;
    state_.box[XX][YY] = 0.0;
    state_.box[XX][ZZ] = 0.0;

    state_.box[YY][XX] = 0.0;
    state_.box[YY][YY] = 10.0;
    state_.box[YY][ZZ] = 0.0;

    state_.box[ZZ][XX] = 0.0;
    state_.box[ZZ][YY] = 0.0;
    state_.box[ZZ][ZZ] = 10.0;

    mdAtoms_.homenr                   = numAtoms_;
    mdAtoms_.haveVsites               = false;
    mdAtoms_.havePartiallyFrozenAtoms = false;

    update_ = std::make_unique<Update>(inputRecord_, kineticEnergyData_, nullptr);
    update_->updateAfterPartition(numAtoms,
                                  gmx::ArrayRef<const unsigned short>(),
                                  mdAtoms_.cTC,
                                  gmx::ArrayRef<const unsigned short>());

    doPressureCouple_ = (nstpcouple != 0);

    if (doPressureCouple_)
    {
        inputRecord_.pressureCouplingOptions.epc        = PressureCoupling::ParrinelloRahman;
        inputRecord_.pressureCouplingOptions.nstpcouple = nstpcouple;
        dtPressureCouple_ = inputRecord_.pressureCouplingOptions.nstpcouple * inputRecord_.delta_t;

        velocityScalingMatrix_(XX, XX) = 1.2;
        velocityScalingMatrix_(XX, YY) = 0.0;
        velocityScalingMatrix_(XX, ZZ) = 0.0;

        velocityScalingMatrix_(YY, XX) = 0.0;
        velocityScalingMatrix_(YY, YY) = 0.8;
        velocityScalingMatrix_(YY, ZZ) = 0.0;

        velocityScalingMatrix_(ZZ, XX) = 0.0;
        velocityScalingMatrix_(ZZ, YY) = 0.0;
        velocityScalingMatrix_(ZZ, ZZ) = 0.9;
    }
    else
    {
        inputRecord_.pressureCouplingOptions.epc = PressureCoupling::No;
        velocityScalingMatrix_(XX, XX)           = 1.0;
        velocityScalingMatrix_(XX, YY)           = 0.0;
        velocityScalingMatrix_(XX, ZZ)           = 0.0;

        velocityScalingMatrix_(YY, XX) = 0.0;
        velocityScalingMatrix_(YY, YY) = 1.0;
        velocityScalingMatrix_(YY, ZZ) = 0.0;

        velocityScalingMatrix_(ZZ, XX) = 0.0;
        velocityScalingMatrix_(ZZ, YY) = 0.0;
        velocityScalingMatrix_(ZZ, ZZ) = 1.0;
    }
}

} // namespace test
} // namespace gmx
