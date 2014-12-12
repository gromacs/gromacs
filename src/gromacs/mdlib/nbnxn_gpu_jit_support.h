/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
 *  \brief Declares functions that support JIT compilation (e.g. for OpenCL)
 *
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 */

#ifndef GMX_MDLIB_NBNXN_GPU_JIT_SUPPORT_H
#define GMX_MDLIB_NBNXN_GPU_JIT_SUPPORT_H

#include "gromacs/gmxlib/gpu_utils/gpu_macros.h"
#include "gromacs/legacyheaders/types/hw_info.h"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/simple.h"

struct gmx_gpu_info_t;

/*! \brief Handles any JIT compilation of nbnxn kernels for the GPU given by \p mygpu */
GPU_FUNC_QUALIFIER void
nbnxn_gpu_compile_kernels(int                       gmx_unused  mygpu,
                          int                       gmx_unused  rank,
                          const gmx_gpu_info_t      gmx_unused *gpu_info,
                          const gmx_gpu_opt_t       gmx_unused *gpu_opt,
                          const interaction_const_t gmx_unused *ic) GPU_FUNC_TERM

#endif
