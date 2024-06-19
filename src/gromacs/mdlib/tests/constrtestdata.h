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
 * \brief SHAKE and LINCS tests header.
 *
 * Contains description and constructor for the test data accumulating object,
 * declares CPU- and GPU-based functions used to apply SHAKE or LINCS on the
 * test data.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_CONSTRTESTDATA_H
#define GMX_MDLIB_TESTS_CONSTRTESTDATA_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/hardware/device_management.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/lincs.h"
#include "gromacs/mdlib/shake.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{
namespace test
{

/* \brief
 * Constraints test data structure.
 *
 * Structure to collect all the necessary data, including system coordinates and topology,
 * constraints information, etc. The structure can be reset and reused.
 */
class ConstraintsTestData
{
public:
    //! Human-friendly name for a system
    std::string title_;
    //! Number of atoms
    int numAtoms_;
    //! Topology
    gmx_mtop_t mtop_;
    //! Masses
    std::vector<real> masses_;
    //! Inverse masses
    std::vector<real> invmass_;
    //! Input record (info that usually in .mdp file)
    t_inputrec ir_;
    //! Local topology
    std::unique_ptr<InteractionDefinitions> idef_;
    //! Computational time array (normally used to benchmark performance)
    t_nrnb nrnb_;

    //! Inverse timestep
    real invdt_;
    //! Number of flexible constraints
    int nflexcon_ = 0;
    //! Whether the virial should be computed
    bool computeVirial_;
    //! Scaled virial
    tensor virialScaled_;
    //! If the free energy is computed
    bool compute_dHdLambda_;
    //! If there are atoms with perturbed mass
    bool hasMassPerturbed_ = false;
    //! Lambda value
    real lambda_ = 0.0;
    //! For free energy computation
    real dHdLambda_;

    //! Coordinates before the timestep
    PaddedVector<RVec> x_;
    //! Coordinates after timestep, output for the constraints
    PaddedVector<RVec> xPrime_;
    //! Backup for coordinates (for reset)
    PaddedVector<RVec> xPrime0_;
    //! Intermediate set of coordinates (normally used for projection correction)
    PaddedVector<RVec> xPrime2_;
    //! Velocities
    PaddedVector<RVec> v_;
    //! Backup for velocities (for reset)
    PaddedVector<RVec> v0_;

    //! Constraints data (type1-i1-j1-type2-i2-j2-...)
    std::vector<int> constraints_;
    //! Target lengths for all constraint types
    std::vector<real> constraintsR0_;

    /*! \brief
     * Constructor for the object with all parameters and variables needed by constraints algorithms.
     *
     * This constructor assembles stubs for all the data structures, required to initialize
     * and apply LINCS and SHAKE constraints. The coordinates and velocities before constraining
     * are saved to allow for reset. The constraints data are stored for testing after constraints
     * were applied.
     *
     * \param[in]  title                Human-friendly name of the system.
     * \param[in]  numAtoms             Number of atoms in the system.
     * \param[in]  masses               Atom masses. Size of this vector should be equal to numAtoms.
     * \param[in]  constraints          List of constraints, organized in triples of integers.
     *                                  First integer is the index of type for a constraint, second
     *                                  and third are the indices of constrained atoms. The types
     *                                  of constraints should be sequential but not necessarily
     *                                  start from zero (which is the way they normally are in
     *                                  GROMACS).
     * \param[in]  constraintsR0        Target values for bond lengths for bonds of each type. The
     *                                  size of this vector should be equal to the total number of
     *                                  unique types in constraints vector.
     * \param[in]  computeVirial        Whether the virial should be computed.
     * \param[in]  compute_dHdLambda    Whether free energy should be computed.
     * \param[in]  initialTime          Initial time.
     * \param[in]  timestep             Timestep.
     * \param[in]  x                    Coordinates before integration step.
     * \param[in]  xPrime               Coordinates after integration step, but before constraining.
     * \param[in]  v                    Velocities before constraining.
     * \param[in]  shakeTolerance       Target tolerance for SHAKE.
     * \param[in]  shakeUseSOR          Use successive over-relaxation method for SHAKE iterations.
     *                                  The general formula is:
     *                                     x_n+1 = (1-omega)*x_n + omega*f(x_n),
     *                                  where omega = 1 if SOR is off and may be < 1 if SOR is on.
     * \param[in]  lincsNumIterations   Number of iterations used to compute the inverse matrix.
     * \param[in]  lincsExpansionOrder  The order for algorithm that adjusts the direction of the
     *                                  bond after constraints are applied.
     * \param[in]  lincsWarnAngle       The threshold value for the change in bond angle. When
     *                                  exceeded the program will issue a warning.
     *
     */
    ConstraintsTestData(const std::string&       title,
                        int                      numAtoms,
                        std::vector<real>        masses,
                        std::vector<int>         constraints,
                        std::vector<real>        constraintsR0,
                        bool                     computeVirial,
                        bool                     compute_dHdLambda,
                        real                     initialTime,
                        real                     timestep,
                        const std::vector<RVec>& x,
                        const std::vector<RVec>& xPrime,
                        const std::vector<RVec>& v,
                        real                     shakeTolerance,
                        gmx_bool                 shakeUseSOR,
                        int                      lincsNumIterations,
                        int                      lincsExpansionOrder,
                        real                     lincsWarnAngle);

    /*! \brief
     * Reset the data structure so it can be reused.
     *
     * Set the coordinates and velocities back to their values before
     * constraining. The scaled virial tensor and dHdLambda are zeroed.
     *
     */
    void reset();
};

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_CONSTRTESTDATA_H
