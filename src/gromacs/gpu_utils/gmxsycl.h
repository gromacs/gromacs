/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2020- The GROMACS Authors
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
 * Wraps the complexity of including SYCL in GROMACS.
 *
 * The __SYCL_COMPILER_VERSION macro is used to identify Intel DPCPP compiler.
 * See https://github.com/intel/llvm/pull/2998 for better proposal.
 *
 * Intel SYCL headers use symbol DIM as a template parameter, which gets broken by macro DIM defined
 * in gromacs/math/vectypes.h. Here, we include the SYCL header while temporary undefining this macro.
 * See https://github.com/intel/llvm/issues/2981.
 *
 * Different compilers, at the time of writing, have different names for some of the proposed features
 * of the SYCL2020 standard. For uniformity, they are all aliased in our custom sycl_2020 namespace.
 *
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GMXSYCL_H
#define GMX_GPU_UTILS_GMXSYCL_H

#include "config.h"

// For hipSYCL, we need to activate floating-point atomics
#if GMX_SYCL_HIPSYCL
#    define HIPSYCL_EXT_FP_ATOMICS
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wunused-variable"
#    pragma clang diagnostic ignored "-Wunused-parameter"
#    pragma clang diagnostic ignored "-Wmissing-noreturn"
#    pragma clang diagnostic ignored "-Wshadow-field"
#    pragma clang diagnostic ignored "-Wctad-maybe-unsupported"
#    pragma clang diagnostic ignored "-Wdeprecated-copy-dtor"
#    pragma clang diagnostic ignored "-Winconsistent-missing-destructor-override"
#    pragma clang diagnostic ignored "-Wunused-template"
#    pragma clang diagnostic ignored "-Wsign-compare"
#    pragma clang diagnostic ignored "-Wundefined-reinterpret-cast"
#    pragma clang diagnostic ignored "-Wdeprecated-copy"
#    pragma clang diagnostic ignored "-Wnewline-eof"
#    pragma clang diagnostic ignored "-Wextra-semi"
#    pragma clang diagnostic ignored "-Wsuggest-override"
#    pragma clang diagnostic ignored "-Wsuggest-destructor-override"
#    pragma clang diagnostic ignored "-Wgcc-compat"
#    include <SYCL/sycl.hpp>
#    pragma clang diagnostic pop
#else // DPC++
// Needed for CUDA targets https://github.com/intel/llvm/issues/5936, enabled for SPIR automatically
#    if defined(__SYCL_DEVICE_ONLY__) && defined(__NVPTX__)
#        define SYCL_USE_NATIVE_FP_ATOMICS 1
#    endif
// DPC++ has issues with DIM macro (https://github.com/intel/llvm/issues/2981)
// and has no SYCL/sycl.hpp up to oneAPI 2022.0
#    ifdef DIM
#        if DIM != 3
#            error "The workaround here assumes we use DIM=3."
#        else
#            undef DIM
#            include <CL/sycl.hpp>
#            define DIM 3
#        endif
#    else
#        include <CL/sycl.hpp>
#    endif
#endif

/* Exposing Intel-specific extensions in a manner compatible with SYCL2020 provisional spec.
 * Despite ICPX (up to 2021.3.0 at the least) having SYCL_LANGUAGE_VERSION=202001,
 * some parts of the spec are still in custom sycl::ONEAPI namespace (sycl::ext::oneapi in beta versions),
 * and some functions have different names. To make things easier to upgrade
 * in the future, this thin layer is added.
 * */
namespace sycl_2020
{
namespace detail
{
#if GMX_SYCL_DPCPP && defined(__INTEL_LLVM_COMPILER) && (__INTEL_LLVM_COMPILER < 20220100)
namespace origin = sycl::ext::oneapi;
#elif GMX_SYCL_HIPSYCL || GMX_SYCL_DPCPP
namespace origin = ::sycl;
#else
#    error "Unsupported version of SYCL compiler"
#endif
} // namespace detail

using detail::origin::atomic_ref;
using detail::origin::memory_order;
using detail::origin::memory_scope;

#if GMX_SYCL_DPCPP
template<typename dataT, int dimensions = 1>
using local_accessor =
        sycl::accessor<dataT, dimensions, sycl::access_mode::read_write, sycl::target::local>;
#elif GMX_SYCL_HIPSYCL
template<typename dataT, int dimensions = 1>
using local_accessor = sycl::local_accessor<dataT, dimensions>;
#endif

} // namespace sycl_2020

#endif
