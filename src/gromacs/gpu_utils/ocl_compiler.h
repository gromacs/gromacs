/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
 *  \brief Declare infrastructure for OpenCL JIT compilation
 *
 *  \author Dimitrios Karkoulis <dimitris.karkoulis@gmail.com>
 *  \author Anca Hamuraru <anca@streamcomputing.eu>
 *  \author Teemu Virolainen <teemu@streamcomputing.eu>
 *  \author Mark Abraham <mark.j.abraham@gmail.com>
 *  \inlibraryapi
 */
#ifndef GMX_GPU_UTILS_OCL_COMPILER_H
#define GMX_GPU_UTILS_OCL_COMPILER_H

#include <string>

#include "gromacs/gpu_utils/oclutils.h"
#include "gromacs/hardware/gpu_hw_info.h"

namespace gmx
{
namespace ocl
{

/*! \brief Compile the specified kernel for the context and device.
 *
 * \param[out] fplog                 Open file pointer for log output
 * \param[in]  kernelBaseFilename    The name of the kernel source file to compile, e.g. "nbnxn_ocl_kernels.cl"
 * \param[in]  extraDefines          Preprocessor defines required by the calling code, e.g. for configuring the kernels
 * \param[in]  context               OpenCL context on the device to compile for
 * \param[in]  deviceId              OpenCL device id of the device to compile for
 * \param[in]  deviceVendorId        Enumerator of the device vendor to compile for
 *
 * \returns The compiled OpenCL program
 *
 * \todo Consider whether we can parallelize the compilation of all
 * the kernels by compiling them in separate programs - but since the
 * resulting programs can't refer to each other, that might lead to
 * bloat of util code?
 *
 * \throws std::bad_alloc  if out of memory.
 *         FileIOError     if a file I/O error prevents returning a valid compiled program.
 *         InternalError   if an OpenCL API error prevents returning a valid compiled program. */
cl_program
compileProgram(FILE              *fplog,
               const std::string &kernelBaseFilename,
               const std::string &extraDefines,
               cl_context         context,
               cl_device_id       deviceId,
               ocl_vendor_id_t    deviceVendorId);

} // namespace
} // namespace

#endif
