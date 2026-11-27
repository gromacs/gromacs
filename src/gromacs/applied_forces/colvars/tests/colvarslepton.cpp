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
 * \brief Tests Lepton custom functions through the Colvars force provider.
 *
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include <map>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/applied_forces/colvars/colvarsforceprovider.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{
namespace test
{
namespace
{

class ColvarsLeptonTest : public ::testing::TestWithParam<real>
{
protected:
    void SetUp() override
    {
        init_t_atoms(&atoms_, 3, false);
        for (int i = 0; i < atoms_.nr; ++i)
        {
            atoms_.atom[i].m = 1;
            atoms_.atom[i].q = 0;
        }
    }

    void TearDown() override { done_atom(&atoms_); }

    t_atoms atoms_ = {};
};

TEST_P(ColvarsLeptonTest, CalculatesCustomFunctionEnergyAndForces)
{
    const std::string colvarsConfigString = R"(units gromacs
        colvarsTrajFrequency 0
        colvar {
            name squared_distance
            customFunction r^2
            distance {
                name r
                group1 {
                    atomNumbers 1
                }
                group2 {
                    atomNumbers 2
                }
            }
        }
        harmonic {
            colvars squared_distance
            forceConstant 2000
            centers 0.3
        })";

    const real              distance           = GetParam();
    const std::vector<RVec> coordinates        = { RVec{ 0, 0, 0 },
                                                   RVec{ real(0.6) * distance, real(0.8) * distance, 0 },
                                                   RVec{ 1, 1, 1 } };
    const std::vector<RVec> colvarsCoordinates = { coordinates[0], coordinates[1] };
    // box is not actually used, given PbcType::No
    matrix                                   box = { { 3, 0, 0 }, { 0, 3, 0 }, { 0, 0, 3 } };
    const MpiComm                            mpiComm(MpiComm::SingleRank{});
    LocalAtomSetManager                      atomSetManager;
    MDLogger                                 logger;
    const std::map<std::string, std::string> inputStrings;
    const ColvarsForceProviderState          state;

    ColvarsForceProvider forceProvider(colvarsConfigString,
                                       atoms_,
                                       PbcType::No,
                                       logger,
                                       inputStrings,
                                       300,
                                       123456,
                                       &atomSetManager,
                                       mpiComm,
                                       nullptr,
                                       0.002,
                                       colvarsCoordinates,
                                       "",
                                       state);

    ForceProviderInput forceProviderInput(coordinates, atoms_.nr, {}, {}, 0.0, 0, box, mpiComm, nullptr);
    std::vector<RVec>   forces(atoms_.nr, RVec{ 0, 0, 0 });
    ForceWithVirial     forceWithVirial(forces, true);
    gmx_enerdata_t      energies(1, nullptr);
    ForceProviderOutput forceProviderOutput(&forceWithVirial, &energies);

    forceProvider.calculateForces(forceProviderInput, &forceProviderOutput);

    // E = k/2 * (r^2 - c)^2; the chain rule contributes a factor of 2*r.
    const double displacement   = double(distance) * double(distance) - 0.3;
    const double expectedEnergy = 0.5 * 2000 * displacement * displacement;
    const double forceScale     = 2 * 2000 * displacement;
    EXPECT_NEAR(expectedEnergy, energies.term[InteractionFunction::CenterOfMassPullingEnergy], 1e-3);
    for (int dim = 0; dim < DIM; ++dim)
    {
        SCOPED_TRACE(dim);
        const double expectedForce = forceScale * coordinates[1][dim];
        EXPECT_NEAR(expectedForce, forces[0][dim], 1e-3);
        EXPECT_NEAR(-expectedForce, forces[1][dim], 1e-3);
        EXPECT_NEAR(0.0, forces[2][dim], 1e-3);
    }
}

INSTANTIATE_TEST_SUITE_P(Distances, ColvarsLeptonTest, ::testing::Values(real(0.4), real(0.8)));

} // namespace
} // namespace test
} // namespace gmx
