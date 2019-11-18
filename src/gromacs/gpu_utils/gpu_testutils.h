/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 *  \brief Declare functions for detection of GPU devices, specific for tests.
 *
 *  \todo This should eventually go to src/testutils
 *
 *  \author Artem Zhmurov <zhmurov@gmail.com>
 *
 *  \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GPU_TESTUTILS_H
#define GMX_GPU_UTILS_GPU_TESTUTILS_H

/*! \brief Checks if there is a compatible GPU to run the computations on
 *
 * There are several reasons why code can not rune on the GPU:
 * 1. The GPU can not be detected, because there is none in the system.
 * 2. GPU detection is disabled by GMX_DISABLE_GPU_DETECTION environmental variable.
 * 3. GPUs are detected, but none of them is compatible.
 * This function checks all these conditions and returns true only if there at least
 * one GPU that can be used for computations.
 *
 * \returns True, if there a GPU that can be used for computations
 */
bool canComputeOnGpu();

#endif // GMX_GPU_UTILS_GPU_TESTUTILS_H
