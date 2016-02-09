/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016, by the GROMACS development team, led by
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
#ifndef CUDA_ARCH_UTILS_CUH_
#define CUDA_ARCH_UTILS_CUH_

/*! \file
 *  \brief CUDA arch dependent definitions.
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */

/* GMX_PTX_ARCH is set to the virtual arch (PTX) version targeted by
 * the current compiler pass or zero for the host pass and it is
 * intended to be used instead of __CUDA_ARCH__.
 */
#ifndef __CUDA_ARCH__
    #define GMX_PTX_ARCH 0
#else
    #define GMX_PTX_ARCH __CUDA_ARCH__
#endif

/* Until CC 5.2 and likely for the near future all NVIDIA architectures
   have a warp size of 32, but this could change later. If it does, the
   following constants should depend on the value of GMX_PTX_ARCH.
 */
static const int warp_size      = 32;
static const int warp_size_log2 = 5;

#endif /* CUDA_ARCH_UTILS_CUH_ */
