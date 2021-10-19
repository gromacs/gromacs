/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
#endif


#ifdef DIM
#    if DIM != 3
#        error "The workaround here assumes we use DIM=3."
#    else
#        undef DIM
#        include <CL/sycl.hpp>
#        define DIM 3
#    endif
#else
#    include <CL/sycl.hpp>
#endif

#if GMX_SYCL_HIPSYCL
#    pragma clang diagnostic pop
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
#if GMX_SYCL_DPCPP
// Confirmed to work for 2021.1-beta10 (20201005) to 2021.3.0 (20210619).
// Deprecated in favor of sycl::ext::oneapi on 20210717 in https://github.com/intel/llvm/commit/d703f578.
// Removed on 20210927 with https://github.com/intel/llvm/pull/4488
#    if __clang_major__ >= 14
namespace origin = sycl::ext::oneapi;
#    else
namespace origin = cl::sycl::ONEAPI;
#    endif
#elif GMX_SYCL_HIPSYCL
namespace origin = cl::sycl;
#else
#    error "Unsupported version of SYCL compiler"
#endif
} // namespace detail

using detail::origin::memory_order;
using detail::origin::memory_scope;
using detail::origin::plus;
using detail::origin::sub_group;

#if GMX_SYCL_DPCPP
using detail::origin::atomic_ref;
template<typename... Args>
bool group_any_of(Args&&... args)
{
    return detail::origin::any_of(std::forward<Args>(args)...);
}
template<typename... Args>
auto group_reduce(Args&&... args) -> decltype(detail::origin::reduce(std::forward<Args>(args)...))
{
    return detail::origin::reduce(std::forward<Args>(args)...);
}
#elif GMX_SYCL_HIPSYCL
using detail::origin::atomic_ref;
using detail::origin::group_any_of;
using detail::origin::group_reduce;
#else
#    error "Unsupported SYCL compiler"
#endif

} // namespace sycl_2020

#endif
