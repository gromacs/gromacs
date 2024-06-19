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
 * \brief Tests for the Leap-Frog integrator
 *
 * \todo Add anisotropic Parrinello-Rahman and other pressure coupling schemes
 * \todo Add PBC handling test.
 * \todo Reference values tests.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_LEAPFROGTESTDATA_H
#define GMX_MDLIB_TESTS_LEAPFROGTESTDATA_H

#include <memory>
#include <vector>

#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/math/matrix.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/update.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{

/* \brief Declares the class to accumulate the data needed for the Leap-Frog integrator tests
 *
 */
class LeapFrogTestData
{
public:
    //! Number of atoms in the system
    int numAtoms_;
    //! Integration timestep
    real timestep_;

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

    //! If the pressure coupling is enabled
    bool doPressureCouple_;
    //! Period between pressure coupling steps
    float dtPressureCouple_;
    //! Matrix for Parrinello-Rahman velocity scaling
    Matrix3x3 velocityScalingMatrix_;

    /*! \brief Constructor.
     *
     * \param[in]  numAtoms          Number of atoms in the system
     * \param[in]  timestep          Integration timestep
     * \param[in]  v0                Initial velocity (same for all particles)
     * \param[in]  f0                External constant force, acting on all particles
     * \param[in]  numTCoupleGroups  Number of temperature coupling groups (zero for no temperature coupling)
     * \param[in]  nstpcouple        Number of steps between pressure coupling steps (zero for no pressure coupling)
     */
    LeapFrogTestData(int numAtoms, real timestep, const rvec v0, const rvec f0, int numTCoupleGroups, int nstpcouple);
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_LEAPFROGTESTDATA_H
