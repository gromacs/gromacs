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
/*! \libinternal \file
 * \brief
 * Wraps the complexity of including OpenCL in Gromacs.
 *
 * Because OpenCL 2.0 is not officially supported widely, \Gromacs
 * uses earlier interfaces. Some of those have been deprecated in 2.0,
 * and generate warnings, which we need to suppress.
 *
 * Additionally, this code wraps they way that things work differently
 * on Apple platforms.
 *
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GMXOPENCL_H
#define GMX_GPU_UTILS_GMXOPENCL_H

/*! \brief Declare to OpenCL SDKs that we intend to use OpenCL API
   features that were deprecated in 1.2 or 2.0, so that they don't
   warn about it. */
///@{
#  define CL_USE_DEPRECATED_OPENCL_1_1_APIS
#  define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#  define CL_USE_DEPRECATED_OPENCL_2_0_APIS
///@}
#  ifdef __APPLE__
#    include <OpenCL/opencl.h>
#  else
#    include <CL/opencl.h>
#  endif

#endif
