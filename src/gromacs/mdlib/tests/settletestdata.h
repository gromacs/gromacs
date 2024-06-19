/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 * \brief SETTLE tests header.
 *
 * Declares the class that accumulates SETTLE test data.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_TESTS_SETTLETESTDATA_H
#define GMX_MDLIB_TESTS_SETTLETESTDATA_H

#include <memory>
#include <vector>

#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

/* \brief SETTLE test data object.
 *
 * Initializes and stores data necessary to run SETTLE constraints, including
 * atom coordinates and velocities, virial, system topology and some parameters.
 */
class SettleTestData
{
public:
    //! Initial (undisturbed) positions
    PaddedVector<gmx::RVec> x_;
    //! Updated water atom positions to constrain
    PaddedVector<gmx::RVec> xPrime_;
    //! Water atom velocities to constrain
    PaddedVector<gmx::RVec> v_;
    //! SETTLE virial
    tensor virial_ = { { 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0 } };

    //! Global topology
    gmx_mtop_t mtop_;
    //! Number of atoms
    int numAtoms_ = 0;
    //! Atom masses
    std::vector<real> masses_;
    //! Reciprocal masses
    std::vector<real> inverseMasses_;
    //! Local topology
    std::unique_ptr<InteractionDefinitions> idef_;

    //! Inverse timestep
    const real reciprocalTimeStep_ = 1.0 / 0.002;
    //! Target distance between oxygen and hydrogens
    const real dOH_ = 0.09572;
    //! Target distance between hydrogens
    const real dHH_ = 0.15139;
    //! Mass of oxygen atom
    const real oxygenMass_ = 15.9994;
    //! Mass of hydrogen atom
    const real hydrogenMass_ = 1.008;

    //! Stride for array with atom indexes
    const int atomsPerSettle_ = NRAL(F_SETTLE);

    /*! \brief Construct the object and initialize the data structures.
     *
     * \param[in] numSettles   Number of SETTLE constraints in the system.
     *
     */
    SettleTestData(int numSettles);

    ~SettleTestData();
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_SETTLETESTDATA_H
