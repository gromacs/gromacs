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
 * \inlibraryapi
 */

#ifndef GMX_GPU_UTILS_GMXSYCL_H
#define GMX_GPU_UTILS_GMXSYCL_H

#include "config.h"

#include <sycl/sycl.hpp>

/* Macro to optimize runtime performance by not recording unnecessary events.
 *
 * It relies on the availability of HIPSYCL_EXT_CG_PROPERTY_* extension, and is no-op for
 * other SYCL implementations. Macro can be used as follows (note the lack of comma after it):
 * `queue.submit(GMX_SYCL_DISCARD_EVENT [=](....))`.
 *
 * When this macro is added to `queue.submit`, the returned event should not be used!
 * As a consequence, patterns like `queue.submit(GMX_SYCL_DISCARD_EVENT [=](....)).wait()`
 * must be avoided. If you intend to use the returned event in any way, do not add this macro.
 *
 * The use of the returned event will not necessarily cause run-time errors, but can cause
 * performance degradation (specifically, in hipSYCL the synchronization will be sub-optimal).
 */
#if GMX_SYCL_HIPSYCL
namespace gmx::internal
{
static const sycl::property_list sc_syclDiscardEventProperty_list{
    sycl::property::command_group::hipSYCL_coarse_grained_events()
};
}
#    define GMX_SYCL_DISCARD_EVENT gmx::internal::sc_syclDiscardEventProperty_list,
#else // IntelLLVM does not support command-group properties
#    define GMX_SYCL_DISCARD_EVENT
#endif

#endif
