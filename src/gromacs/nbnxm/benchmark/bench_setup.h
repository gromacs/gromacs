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

/*! \libinternal \file
 * \brief
 * This file declares functions for setting up kernel benchmarks
 *
 * \author Berk Hess <hess@kth.se>
 *
 * \inlibraryapi
 * \ingroup module_nbnxm
 */

#ifndef GMX_NBNXN_BENCH_SETUP_H
#define GMX_NBNXN_BENCH_SETUP_H

#include <string>

#include "gromacs/utility/real.h"

namespace gmx
{

//! Enum for selecting the SIMD kernel type for benchmarks
enum class NbnxmBenchMarkKernels : int
{
    SimdAuto,
    SimdNo,
    Simd4XM,
    Simd2XMM,
    Count
};

//! Enum for selecting the combination rule for kernel benchmarks
enum class NbnxmBenchMarkCombRule : int
{
    RuleGeom,
    RuleLB,
    RuleNone,
    Count
};

//! Enum for selecting coulomb type for kernel benchmarks
enum class NbnxmBenchMarkCoulomb : int
{
    Pme,
    ReactionField,
    Count
};

/*! \internal \brief
 * The options for the kernel benchmarks
 */
struct NbnxmKernelBenchOptions
{
    //! Whether to use a GPU, currently GPUs are not supported
    bool useGpu = false;
    //! The number of OpenMP threads to use
    int numThreads = 1;
    //! The SIMD type for the kernel
    NbnxmBenchMarkKernels nbnxmSimd = NbnxmBenchMarkKernels::SimdAuto;
    //! The LJ combination rule
    NbnxmBenchMarkCombRule ljCombinationRule = NbnxmBenchMarkCombRule::RuleGeom;
    //! Use i-cluster half-LJ optimization for clusters with <= half LJ
    bool useHalfLJOptimization = false;
    //! The pairlist and interaction cut-off
    real pairlistCutoff = 1.0;
    //! The Coulomb Ewald coefficient
    real ewaldcoeff_q = 0;
    //! Whether to compute energies (shift forces for virial are always computed on CPU)
    bool computeVirialAndEnergy = false;
    //! The Coulomb interaction function
    NbnxmBenchMarkCoulomb coulombType = NbnxmBenchMarkCoulomb::Pme;
    //! Whether to use tabulated PME grid correction instead of analytical, not applicable with simd=no
    bool useTabulatedEwaldCorr = false;
    //! Whether to run all combinations of Coulomb type, combination rule and SIMD
    bool doAll = false;
    //! Number of iterations to run before running each kernel benchmark, currently always 1
    int numPreIterations = 1;
    //! The number of iterations for each kernel
    int numIterations = 100;
    //! The number of (untimed) iterations to run at startup to warm up the CPU
    int numWarmupIterations = 0;
    //! Print cycles/pair instead of pairs/cycle
    bool cyclesPerPair = false;
    //! Report in micro seconds instead of cycles
    bool reportTime = false;
    //! Also report into a csv file
    std::string outputFile;
};

/*! \brief
 * Sets up and runs one or more Nbnxm kernel benchmarks
 *
 * The simulated system is a box of 1000 SPC/E water molecules scaled
 * by the factor \p sizeFactor, which has to be a power of 2.
 * One or more benchmarks are run, as specified by \p options.
 * Benchmark settings and timings are printed to stdout.
 *
 * \param[in] sizeFactor How much should the system size be increased.
 * \param[in] options How the benchmark will be run.
 */
void bench(int sizeFactor, const NbnxmKernelBenchOptions& options);

} // namespace gmx

#endif
