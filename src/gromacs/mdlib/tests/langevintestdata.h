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
 * \brief Tests for the Langevin integrator
 *
 * \todo Add PBC handling test.
 * \todo Reference values tests.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \author Magnus Lundborg <magnus.lundborg@scilifelab.se>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_LANGEVINTESTDATA_H
#define GMX_MDLIB_TESTS_LANGEVINTESTDATA_H

#include <cstdint>

#include <memory>
#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

/* \brief Declares the class to accumulate the data needed for the Langevin integrator tests
 *
 */
class LangevinTestData
{
public:
    //! Number of atoms in the system
    int numAtoms_;
    //! Integration timestep
    real timestep_;
    //! Current step number
    int64_t step_;
    //! Temperature
    real temperature_;
    //! Tau-t (inverse friction constant)
    real tauT_;
    //! Random seed
    int seed_;

    //! Initial coordinates
    PaddedVector<RVec> x0_;
    //! Current coordinates
    PaddedVector<RVec> x_;
    //! Coordinates after integrator update
    PaddedVector<RVec> xPrime_;
    //! Initial velocities
    PaddedVector<RVec> v0_;
    //! Current velocities
    PaddedVector<RVec> v_;
    //! External forces
    PaddedVector<RVec> f_;
    //! Inverse masses of the particles
    PaddedVector<real> inverseMasses_;
    //! Inverse masses of the particles per dimension
    std::vector<RVec> inverseMassesPerDim_;

    //! MD atoms structure in which inverse masses will be passed to the integrator
    t_mdatoms mdAtoms_;
    //! Input record (to get integrator type, temperature and pressure coupling)
    t_inputrec inputRecord_;
    //! System state
    t_state state_;
    //! Force calculation data
    t_fcdata forceCalculationData_;
    //! Kinetic energy data (to disable non-equilibrium MD integration)
    gmx_ekindata_t kineticEnergyData_;
    //! Update data
    std::unique_ptr<Update> update_;

    //! Number of temperature coupling groups
    int numTCoupleGroups_;

    LangevinTestData(int        numAtoms,
                     real       timestep,
                     const RVec v0,
                     const RVec f0,
                     int        numTCoupleGroups,
                     real       temperature,
                     real       tauT,
                     int        seed);
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_LANGEVINTESTDATA_H
