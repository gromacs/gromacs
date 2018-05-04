/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
//This is a top-level file to generate all PME openCL kernels

/* SPREAD/SPLINE */


#define c_spreadMaxWarpsPerBlock 8
#define c_spreadMaxThreadsPerBlock (c_spreadMaxWarpsPerBlock * warp_size)
#define atomsPerBlock  (c_spreadMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM)


// splineAndSpread
#define computeSplines 1
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernel
#include "pme-spread-kernel.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spline
#define computeSplines 1
#define spreadCharges 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernel
#include "pme-spread-kernel.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spread
#define computeSplines 0
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernel
#include "pme-spread-kernel.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME


#undef atomsPerBlock


/* GATHER */

// FIXME these are duplicates of host-side consts
// moreover, c_gatherMaxWarpsPerBlock should be defien through c_gatherMaxThreadsPerBlock, probably
// same for spread
#define c_gatherMaxWarpsPerBlock 4
#define c_gatherMaxThreadsPerBlock (c_gatherMaxWarpsPerBlock * warp_size)
#define atomsPerBlock (c_gatherMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM)


// gather
#define overwriteForces 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherKernel
#include "pme-gather-kernel.cl"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME

// gather with reduction
#define overwriteForces 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherReduceWithInputKernel
#include "pme-gather-kernel.cl"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME


/* SOLVE */
//test
//#undef warp_size
#include "pme-ocl-definitely-common.h"

//GRID orderings
#define YZX 1
#define XYZ 2

// solve, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXKernel
#include "pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernel
#include "pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernel
#include "pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernel
#include "pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

#undef atomsPerBlock
