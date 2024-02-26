/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
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
#ifndef GMX_GPU_UTILS_MACROS_H
#define GMX_GPU_UTILS_MACROS_H

#include "config.h"

#include "gromacs/utility/basedefinitions.h" // for gmx_unused

/* These macros that let us define inlineable null implementations so
   that non-GPU Gromacs can run with no overhead without conditionality
   everywhere a GPU function is called. */
#define REAL_FUNC_QUALIFIER
#define REAL_FUNC_ARGUMENT(arg) arg
#define REAL_FUNC_TERM
#define REAL_FUNC_TERM_WITH_RETURN(arg)

#define NULL_FUNC_QUALIFIER gmx_unused static
#define NULL_FUNC_ARGUMENT(arg) arg gmx_unused
#define NULL_FUNC_TERM \
    {                  \
    }
#define NULL_FUNC_TERM_WITH_RETURN(arg) \
    {                                   \
        return (arg);                   \
    }

#ifdef DOXYGEN

/* Doxygen build appreciates always having argument names, and doesn't
 * care about duplicate function definitions. */
#    define GPU_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#    define GPU_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#    define GPU_FUNC_TERM REAL_FUNC_TERM
#    define GPU_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    define CUDA_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#    define CUDA_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#    define CUDA_FUNC_TERM REAL_FUNC_TERM
#    define CUDA_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    define OPENCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#    define OPENCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#    define OPENCL_FUNC_TERM REAL_FUNC_TERM
#    define OPENCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    define SYCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#    define SYCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#    define SYCL_FUNC_TERM REAL_FUNC_TERM
#    define SYCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)

#else // Not DOXYGEN

/* GPU support is enabled, so these functions will have real code defined somewhere */
#    if GMX_GPU && !GMX_GPU_HIP
#        define GPU_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#        define GPU_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#        define GPU_FUNC_TERM REAL_FUNC_TERM
#        define GPU_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    else
#        define GPU_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#        define GPU_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#        define GPU_FUNC_TERM NULL_FUNC_TERM
#        define GPU_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#    endif

/* Enable and disable platform-specific function implementations */
#    if GMX_GPU_OPENCL
#        define OPENCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#        define OPENCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#        define OPENCL_FUNC_TERM REAL_FUNC_TERM
#        define OPENCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    else
#        define OPENCL_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#        define OPENCL_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#        define OPENCL_FUNC_TERM NULL_FUNC_TERM
#        define OPENCL_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#    endif

#    if GMX_GPU_CUDA
#        define CUDA_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#        define CUDA_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#        define CUDA_FUNC_TERM REAL_FUNC_TERM
#        define CUDA_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    else
#        define CUDA_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#        define CUDA_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#        define CUDA_FUNC_TERM NULL_FUNC_TERM
#        define CUDA_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#    endif

#    if GMX_GPU_SYCL
#        define SYCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#        define SYCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#        define SYCL_FUNC_TERM REAL_FUNC_TERM
#        define SYCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#    else
#        define SYCL_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#        define SYCL_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#        define SYCL_FUNC_TERM NULL_FUNC_TERM
#        define SYCL_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#    endif

#endif // ifdef DOXYGEN

#endif // GMX_GPU_UTILS_MACROS_H
