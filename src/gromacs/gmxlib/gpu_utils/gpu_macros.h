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
#ifndef GMX_GMXLIB_GPU_UTILS_MACROS_H
#define GMX_GMXLIB_GPU_UTILS_MACROS_H

#include "config.h"

/* These macros that let us define inlineable null implementations so
   that non-GPU Gromacs can run with no overhead without conditionality
   everywhere a GPU function is called. */
#define REAL_FUNC_QUALIFIER
#define REAL_FUNC_ARGUMENT(arg) arg
#define REAL_FUNC_TERM ;
#define REAL_FUNC_TERM_WITH_RETURN(arg) ;

#define NULL_FUNC_QUALIFIER static
#define NULL_FUNC_ARGUMENT(arg) /*arg*/
#define NULL_FUNC_TERM {}
#define NULL_FUNC_TERM_WITH_RETURN(arg) { return (arg); }

#ifdef DOXYGEN

/* Doxygen build appreciates always having argument names, and doesn't
 * care about duplicate function definitions. */
#define GPU_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define GPU_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define GPU_FUNC_TERM REAL_FUNC_TERM
#define GPU_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#define CUDA_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define CUDA_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define CUDA_FUNC_TERM REAL_FUNC_TERM
#define CUDA_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#define OPENCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define OPENCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define OPENCL_FUNC_TERM REAL_FUNC_TERM
#define OPENCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)

#elif defined GMX_GPU

/* GPU support is enabled, so these functions will have real code
 * defined somewhere */
#define GPU_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define GPU_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define GPU_FUNC_TERM REAL_FUNC_TERM
#define GPU_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)

#  if defined GMX_USE_OPENCL

/* OpenCL support is enabled, so CUDA-specific functions need empty
 * implementations, while OpenCL-specific functions will have real
 * code defined somewhere. */
#define CUDA_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#define CUDA_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#define CUDA_FUNC_TERM NULL_FUNC_TERM
#define CUDA_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#define OPENCL_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define OPENCL_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define OPENCL_FUNC_TERM REAL_FUNC_TERM
#define OPENCL_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)

#  else /* !(defined GMX_USE_OPENCL) */

/* CUDA support is enabled, so OpenCL-specific functions need empty
 * implementations, while CUDA-specific functions will have real
 * code defined somewhere. */
#define CUDA_FUNC_QUALIFIER REAL_FUNC_QUALIFIER
#define CUDA_FUNC_ARGUMENT REAL_FUNC_ARGUMENT
#define CUDA_FUNC_TERM REAL_FUNC_TERM
#define CUDA_FUNC_TERM_WITH_RETURN(arg) REAL_FUNC_TERM_WITH_RETURN(arg)
#define OPENCL_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#define OPENCL_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#define OPENCL_FUNC_TERM NULL_FUNC_TERM
#define OPENCL_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)

#  endif

#else /* !(defined DOXYGEN) && !(defined GMX_GPU) */

/* No GPU support is configured, so none of these functions will have
 * real definitions. */
#define GPU_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#define GPU_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#define GPU_FUNC_TERM NULL_FUNC_TERM
#define GPU_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#define CUDA_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#define CUDA_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#define CUDA_FUNC_TERM NULL_FUNC_TERM
#define CUDA_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)
#define OPENCL_FUNC_QUALIFIER NULL_FUNC_QUALIFIER
#define OPENCL_FUNC_ARGUMENT NULL_FUNC_ARGUMENT
#define OPENCL_FUNC_TERM NULL_FUNC_TERM
#define OPENCL_FUNC_TERM_WITH_RETURN(arg) NULL_FUNC_TERM_WITH_RETURN(arg)

#endif

#endif
