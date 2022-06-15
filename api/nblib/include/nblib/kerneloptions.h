/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
/*! \inpublicapi \file
 * \brief
 * Implements nblib kernel setup options
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef NBLIB_KERNELOPTIONS_H
#define NBLIB_KERNELOPTIONS_H

#include <memory>

#include "nblib/basicdefinitions.h"

namespace nblib
{

//! Enum for selecting the SIMD kernel type
enum class SimdKernels : int
{
    SimdAuto,
    SimdNo,
    Simd4XM,
    Simd2XMM,
    Count
};

//! Enum for selecting the combination rule
enum class CombinationRule : int
{
    Geometric,
    LorentzBerthelot,
    None,
    Count
};

//! Enum for selecting coulomb type
enum class CoulombType : int
{
    Pme,
    Cutoff,
    ReactionField,
    Count
};

/*! \internal \brief
 * The options for the nonbonded kernel caller
 */
struct NBKernelOptions final
{
    //! Whether to use a GPU, currently GPUs are not supported
    bool useGpu = false;
    //! The number of OpenMP threads to use
    int numOpenMPThreads = 1;
    //! The SIMD type for the kernel
    SimdKernels nbnxmSimd = SimdKernels::SimdAuto;
    //! The LJ combination rule
    CombinationRule ljCombinationRule = CombinationRule::Geometric;
    //! The pairlist and interaction cut-off
    real pairlistCutoff = 1.0;
    //! The Coulomb interaction function
    CoulombType coulombType = CoulombType::Pme;
    //! Whether to use tabulated PME grid correction instead of analytical, not applicable with simd=no
    bool useTabulatedEwaldCorr = false;
    //! The number of iterations for each kernel
    int numIterations = 100;
    //! The time step
    real timestep = 0.001;
};

} // namespace nblib

#endif // NBLIB_KERNELOPTIONS_H
