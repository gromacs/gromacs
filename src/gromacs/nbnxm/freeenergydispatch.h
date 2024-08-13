/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *
 * \brief
 * Declares the free-energy kernel dispatch class
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */
#ifndef GMX_NBNXM_FREEENERGYDISPATCH_H
#define GMX_NBNXM_FREEENERGYDISPATCH_H

#include <memory>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/threaded_force_buffer.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct gmx_enerdata_t;
struct gmx_wallcycle;
struct interaction_const_t;
struct t_lambda;
struct t_nrnb;

namespace gmx
{
class PairlistSets;
template<typename>
class ArrayRefWithPadding;
class ForceWithShiftForces;
class StepWorkload;

/*! \internal
 *  \brief Temporary data and methods for handling dispatching of the nbnxm free-energy kernels
 */
class FreeEnergyDispatch
{
public:
    //! Constructor
    FreeEnergyDispatch(int numEnergyGroups);

    //! Sets up the threaded force buffer and the reduction, should be called after constructing the pair lists
    void setupFepThreadedForceBuffer(int numAtomsForce, const PairlistSets& pairlistSets);

    //! Dispatches the non-bonded free-energy kernels, thread parallel and reduces the output
    void dispatchFreeEnergyKernels(const PairlistSets&                    pairlistSets,
                                   const ArrayRefWithPadding<const RVec>& coords,
                                   ForceWithShiftForces*                  forceWithShiftForces,
                                   bool                                   useSimd,
                                   int                                    ntype,
                                   const interaction_const_t&             ic,
                                   ArrayRef<const RVec>                   shiftvec,
                                   ArrayRef<const real>                   nbfp,
                                   ArrayRef<const real>                   nbfp_grid,
                                   ArrayRef<const real>                   chargeA,
                                   ArrayRef<const real>                   chargeB,
                                   ArrayRef<const int>                    typeA,
                                   ArrayRef<const int>                    typeB,
                                   ArrayRef<const real>                   lambda,
                                   gmx_enerdata_t*                        enerd,
                                   const StepWorkload&                    stepWork,
                                   t_nrnb*                                nrnb,
                                   gmx_wallcycle*                         wcycle);

private:
    //! Temporary array for storing foreign lambda group pair energies
    gmx_grppairener_t foreignGroupPairEnergies_;

    //! Threaded force buffer for nonbonded FEP
    ThreadedForceBuffer<RVec> threadedForceBuffer_;
    //! Threaded buffer for nonbonded FEP foreign energies and dVdl, no forces, so numAtoms = 0
    ThreadedForceBuffer<RVec> threadedForeignEnergyBuffer_;
};

} // namespace gmx

#endif // GMX_NBNXM_FREEENERGYDISPATCH_H
