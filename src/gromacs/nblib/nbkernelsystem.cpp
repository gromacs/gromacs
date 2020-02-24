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
#include "gmxpre.h"

#include "nbkerneloptions.h"

#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/utility/logger.h"

namespace nblib
{

//! Add the options instance to the list for all requested kernel SIMD types
//! TODO This should be refactored so that if SimdAuto is set only one kernel
//!      layout is chosen.
//! TODO This should be refactored to only return the desired kernel layout
    static void expandSimdOptionAndPushBack(const NBKernelOptions& options, std::vector<NBKernelOptions>* optionsList)
    {
        if (options.nbnxmSimd == BenchMarkKernels::SimdAuto)
        {
            bool addedInstance = false;
#ifdef GMX_NBNXN_SIMD_4XN
            optionsList->push_back(options);
            optionsList->back().nbnxmSimd = BenchMarkKernels::Simd4XM;
            addedInstance                 = true;
#endif
#ifdef GMX_NBNXN_SIMD_2XNN
            optionsList->push_back(options);
            optionsList->back().nbnxmSimd = BenchMarkKernels::Simd2XMM;
            addedInstance                 = true;
#endif
            if (!addedInstance)
            {
                optionsList->push_back(options);
                optionsList->back().nbnxmSimd = BenchMarkKernels::SimdNo;
            }
        }
        else
        {
            optionsList->push_back(options);
        }
    }


// void nbKernel(NBKernelSystem        &system,
//              const NBKernelOptions &options,
//              const bool            &printTimings)
//{
//    // We don't want to call gmx_omp_nthreads_init(), so we init what we need
//    gmx_omp_nthreads_set(emntPairsearch, options.numThreads);
//    gmx_omp_nthreads_set(emntNonbonded, options.numThreads);
//
//    real                       minBoxSize = norm(system.box[XX]);
//    for (int dim = YY; dim < DIM; dim++)
//    {
//        minBoxSize = std::min(minBoxSize, norm(system.box[dim]));
//    }
//    if (options.pairlistCutoff > 0.5*minBoxSize)
//    {
//        gmx_fatal(FARGS, "The cut-off should be shorter than half the box size");
//    }
//
//    std::vector<NBKernelOptions> optionsList;
//    expandSimdOptionAndPushBack(options, &optionsList);
//    GMX_RELEASE_ASSERT(!optionsList.empty(), "Expect at least one benchmark setup");
//
//    // setupAndRunInstance(system, optionsList[0], printTimings);
//}

} // namespace nblib