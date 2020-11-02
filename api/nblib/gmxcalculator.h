/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Implements a force calculator based on GROMACS data structures.
 *
 * Intended for internal use inside the ForceCalculator.
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#ifndef NBLIB_GMXCALCULATOR_H
#define NBLIB_GMXCALCULATOR_H

#include <memory>

#include "nblib/vector.h"

struct nonbonded_verlet_t;
struct t_forcerec;
struct t_nrnb;
struct interaction_const_t;
struct gmx_enerdata_t;

namespace gmx
{
template<typename T>
class ArrayRef;
class StepWorkload;
} // namespace gmx

namespace nblib
{
class Box;
class NbvSetupUtil;
class SimulationState;
struct NBKernelOptions;

/*! \brief GROMACS non-bonded force calculation backend
 *
 * This class encapsulates the various GROMACS data structures and their interplay
 * from the NBLIB user. The class is a private member of the ForceCalculator and
 * is not intended for the public interface.
 *
 * Handles the task of storing the simulation problem description using the internal
 * representation used within GROMACS. It currently supports short range non-bonded
 * interactions (PP) on a single node CPU.
 *
 */

class GmxForceCalculator final
{
public:
    GmxForceCalculator();

    ~GmxForceCalculator();

    //! Compute forces and return
    void compute(gmx::ArrayRef<const gmx::RVec> coordinateInput, gmx::ArrayRef<gmx::RVec> forceOutput);

    //! Puts particles on a grid based on bounds specified by the box (for every NS step)
    void setParticlesOnGrid(gmx::ArrayRef<const int>       particleInfoAllVdw,
                            gmx::ArrayRef<const gmx::RVec> coordinates,
                            const Box&                     box);

private:
    //! Friend to allow setting up private members in this class
    friend class NbvSetupUtil;

    //! Non-Bonded Verlet object for force calculation
    std::unique_ptr<nonbonded_verlet_t> nbv_;

    //! Only nbfp and shift_vec are used
    std::unique_ptr<t_forcerec> forcerec_;

    //! Parameters for various interactions in the system
    std::unique_ptr<interaction_const_t> interactionConst_;

    //! Tasks to perform in an MD Step
    std::unique_ptr<gmx::StepWorkload> stepWork_;

    //! Energies of different interaction types; currently only needed as an argument for dispatchNonbondedKernel
    std::unique_ptr<gmx_enerdata_t> enerd_;

    //! Non-bonded flop counter; currently only needed as an argument for dispatchNonbondedKernel
    std::unique_ptr<t_nrnb> nrnb_;

    //! Legacy matrix for box
    matrix box_{ { 0 } };
};

} // namespace nblib

#endif // NBLIB_GMXCALCULATOR_H
