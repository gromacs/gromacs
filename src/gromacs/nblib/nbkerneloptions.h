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
/*! \internal \file
 * \brief
 * Implements nblib kernel setup options
 *
 * \author Berk Hess <hess@kth.se>
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */
#ifndef GROMACS_SETUP_H
#define GROMACS_SETUP_H

#include <memory>

#include "gromacs/math/vectypes.h"

#include "nbkerneldef.h"

struct interaction_const_t;
struct nonbonded_verlet_t;

namespace nblib
{

/*! \internal \brief
 * The options for the nonbonded kernel caller
 */
struct NBKernelOptions
{
    //! Whether to use a GPU, currently GPUs are not supported
    bool useGpu = false;
    //! The number of OpenMP threads to use
    int numThreads = 1;
    //! The SIMD type for the kernel
    BenchMarkKernels nbnxmSimd = BenchMarkKernels::SimdAuto;
    //! The LJ combination rule
    CombinationRule ljCombinationRule = CombinationRule::Geometric;
    //! Use i-cluster half-LJ optimization for clusters with <= half LJ
    bool useHalfLJOptimization = false;
    //! The pairlist and interaction cut-off
    real pairlistCutoff = 1.0;
    //! Whether to compute energies (shift forces for virial are always computed on CPU)
    bool computeVirialAndEnergy = false;
    //! The Coulomb interaction function
    BenchMarkCoulomb coulombType = BenchMarkCoulomb::Pme;
    //! Whether to use tabulated PME grid correction instead of analytical, not applicable with simd=no
    bool useTabulatedEwaldCorr = false;
    //! The number of iterations for each kernel
    int numIterations = 100;
    //! Print cycles/pair instead of pairs/cycle
    bool cyclesPerPair = false;
    //! The time step
    real timestep = 0.001;
};

/*! \brief
 * Sets up and runs nonbonded kernel calls
 *
 * The simulated system is a box of 12 Argon molecules scaled
 * by the factor \p sizeFactor, which has to be a power of 2.
 * Timings can be printed to stdout.
 *
 * \param[in,out] system The particle system to compute nonbonded forces for
 * \param[in] options How the benchmark will be run.
 * \param[in] printTimings Whether to print cycle counters
 */
// void nbKernel(NBKernelSystem& system, const NBKernelOptions& options, const bool& printTimings);


} // namespace nblib

#endif // GROMACS_SETUP_H
