/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2024- The GROMACS Authors
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
 * \brief Defines the class to accumulate the data needed for the Langevin integrator tests
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */
#include "gmxpre.h"

#include "langevintestdata.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{
namespace test
{

LangevinTestData::LangevinTestData(int        numAtoms,
                                   real       timestep,
                                   const RVec v0,
                                   const RVec f0,
                                   int        numTCoupleGroups,
                                   real       temperature,
                                   real       tauT,
                                   int        seed) :
    numAtoms_(numAtoms),
    timestep_(timestep),
    temperature_(temperature),
    tauT_(tauT),
    seed_(seed),
    x0_(numAtoms),
    x_(numAtoms),
    xPrime_(numAtoms),
    v0_(numAtoms),
    v_(numAtoms),
    f_(numAtoms),
    inverseMasses_(numAtoms),
    inverseMassesPerDim_(numAtoms),
    kineticEnergyData_(std::vector<real>(numTCoupleGroups == 0 ? 1 : numTCoupleGroups, temperature),
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

    mdAtoms_.ptype.resize(numAtoms_);

    // To do temperature coupling at each step
    inputRecord_.nsttcouple = 1;
    inputRecord_.opts.ngtc  = numTCoupleGroups_;
    snew(inputRecord_.opts.tau_t, numTCoupleGroups_);
    snew(inputRecord_.opts.anneal_time, numTCoupleGroups_);
    snew(inputRecord_.opts.anneal_temp, numTCoupleGroups_);

    GMX_ASSERT(numTCoupleGroups_ > 0,
               "The Langevin (SD) integrator needs temperature coupling groups.");

    inputRecord_.etc = TemperatureCoupling::Yes;
    for (int i = 0; i < numAtoms_; i++)
    {
        mdAtoms_.cTC[i]   = i % numTCoupleGroups_;
        mdAtoms_.ptype[i] = ParticleType::Atom;
    }
    for (int i = 0; i < numTCoupleGroups_; i++)
    {
        real         tCoupleLambda = 1.0 - (i + 1.0) / 10.0;
        t_grp_tcstat temperatureCouplingGroupData;
        temperatureCouplingGroupData.lambda = tCoupleLambda;
        kineticEnergyData_.tcstat[i]        = temperatureCouplingGroupData;
        inputRecord_.opts.tau_t[i]          = tauT;
    }
    snew(inputRecord_.opts.nFreeze, 1);
    snew(inputRecord_.opts.acceleration, 1);
    for (int d = 0; d < DIM; d++)
    {
        inputRecord_.opts.nFreeze[0][d]      = 0;
        inputRecord_.opts.acceleration[0][d] = 0;
    }

    inputRecord_.eI      = IntegrationAlgorithm::SD1;
    inputRecord_.delta_t = timestep_;
    inputRecord_.ld_seed = seed;

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
}

} // namespace test
} // namespace gmx
