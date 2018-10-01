/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * \brief Declares helper functionality for device transfers for tests
 * for GPU host allocator.
 *
 * Undefined symbols in Google Test, GROMACS use of -Wundef, and the
 * implementation of FindCUDA.cmake and/or nvcc mean that no
 * compilation unit should include a gtest header while being compiled
 * by nvcc. None of -isystem, -Wno-undef, nor the pragma GCC
 * diagnostic work.
 *
 * Thus, this header isolates CUDA-specific functionality to its own
 * translation unit. The OpenCL and no-GPU implementations do not
 * require this separation, but do so for consistency.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 */
#ifndef GMX_GPU_UTILS_TESTS_DEVICETRANSFERS_H
#define GMX_GPU_UTILS_TESTS_DEVICETRANSFERS_H

#include "gromacs/utility/arrayref.h"

struct gmx_gpu_info_t;

namespace gmx
{

/*! \brief Helper function for GPU test code to be platform agnostic.
 *
 * Transfers \c input to device 0, if present, and transfers it back
 * into \c output. Both sizes must match. If no devices are present,
 * do a simple host-side buffer copy instead.
 *
 * \throws InternalError  Upon any GPU API error condition. */
void doDeviceTransfers(const gmx_gpu_info_t &gpuInfo,
                       ArrayRef<const char>  input,
                       ArrayRef<char>        output);

} // namespace gmx

#endif
