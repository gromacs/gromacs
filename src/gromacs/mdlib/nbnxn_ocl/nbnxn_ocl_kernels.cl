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
#define __IN_OPENCL_KERNEL__

/* Auxiliary kernels */
__kernel void
memset_f3(__global float3 *buf,const float value,const unsigned int Nbuf)
{
    unsigned int tidx = get_global_id(0);
    if(tidx < Nbuf)
        buf[tidx] = value;
}

__kernel void
memset_f2(__global float2 *buf,const float value,const unsigned int Nbuf)
{
    unsigned int tidx = get_global_id(0);
    if(tidx < Nbuf)
        buf[tidx] = value;
}

__kernel void
memset_f(__global float *buf,const float value,const unsigned int Nbuf)
{
    unsigned int tidx = get_global_id(0);
    if(tidx < Nbuf)
        buf[tidx] = value;
}

/* Very few data */
__kernel void
zero_e_fshift(__global float *fshift,__global float *e_lj,__global float *e_el,const unsigned int Nbuf)
{
    unsigned int tidx = get_global_id(0);
    if(tidx < Nbuf)
        fshift[tidx] = 0.0f;
    if(tidx==0)
    {
        *e_lj     = 0.0f;
        *e_el     = 0.0f;
    }
}

/* Generate pruning kernels. */
#define HAVE_FRESH_LIST 1
#include "nbnxn_ocl_kernel_pruneonly.clh"
#undef HAVE_FRESH_LIST
#include "nbnxn_ocl_kernel_pruneonly.clh"

#if defined GMX_OCL_FASTGEN
    #define FLAVOR_LEVEL_GENERATOR "nbnxn_ocl_kernels_fastgen.clh"
#elif defined GMX_OCL_FASTGEN_ADD_TWINCUT
    #define FLAVOR_LEVEL_GENERATOR "nbnxn_ocl_kernels_fastgen_add_twincut.clh"
#else
    #define FLAVOR_LEVEL_GENERATOR "nbnxn_ocl_kernels.clh"
#endif

/* Top-level kernel generation: will generate through multiple inclusion the
 * following flavors for all kernels:
 * - force-only output;
 * - force and energy output;
 * - force-only with pair list pruning;
 * - force and energy output with pair list pruning.
 */

/** Force only **/
#include FLAVOR_LEVEL_GENERATOR
/** Force & energy **/
#define CALC_ENERGIES
#include FLAVOR_LEVEL_GENERATOR
#undef CALC_ENERGIES

/*** Pair-list pruning kernels ***/
/** Force only **/
#define PRUNE_NBL
#include FLAVOR_LEVEL_GENERATOR
/** Force & energy **/
#define CALC_ENERGIES
#include FLAVOR_LEVEL_GENERATOR
#undef CALC_ENERGIES
#undef PRUNE_NBL
