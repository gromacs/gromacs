/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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

/*! \internal \file
 *  \brief Top-level file for generating PME OpenCL kernels.
 * This file includes all the PME OpenCL kernel files multiple times, with additional defines.
 * Plenty of common constants/macros are also expected to be already defined by the compiler
 * (such as "order" which is PME interpolation order).
 * For details, please see how pme_program.cl is compiled in pme_gpu_program_impl_ocl.cpp.
 *
 * This file specifically expects the following work size definitions:
 *
 * - c_spreadWorkGroupSize is a spline/spread kernels work group size.
 * - c_solveMaxWorkGroupSize is a solve kernels maximum work group size.
 * - c_gatherWorkGroupSize is a gather kernels work group size.
 *
 *  \author Aleksei Iupinov <a.yupinov@gmail.com>
 */

// Assert placeholders, to not rip them out from OpenCL implementation - hopefully they come in handy some day with OpenCL 2
#define assert(a)

#define PmeOpenCLKernelParams PmeGpuKernelParamsBase

/* SPREAD/SPLINE */

#define atomsPerBlock (c_spreadWorkGroupSize / threadsPerAtom)
#define atomsPerWarp (warp_size / threadsPerAtom)

// spline/spread fused single grid (without FEP or FEP without energy and virial calculation at this step)
#define computeSplines 1
#define spreadCharges 1
#define numGrids 1 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernelSingle
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// spline on single grid (without FEP or FEP without energy and virial calculation at this step)
#define computeSplines 1
#define spreadCharges 0
#define numGrids 1 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernelSingle
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// spread on single grid (without FEP or FEP without energy and virial calculation at this step)
#define computeSplines 0
#define spreadCharges 1
#define numGrids 1 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernelSingle
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// spline/spread fused on two grids (FEP with energy and virial calculation at this step)
#define computeSplines 1
#define spreadCharges 1
#define numGrids 2 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernelDual
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// spline on two grids (FEP with energy and virial calculation at this step)
#define computeSplines 1
#define spreadCharges 0
#define numGrids 2 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernelDual
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// spread on two grids (FEP with energy and virial calculation at this step)
#define computeSplines 0
#define spreadCharges 1
#define numGrids 2 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernelDual
#include "pme_spread.clh"
#undef computeSplines
#undef spreadCharges
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME


#undef atomsPerBlock


/* GATHER */

#define atomsPerBlock (c_gatherWorkGroupSize / threadsPerAtom)

// gather using single grid
#define numGrids 1 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherKernelSingle
#include "pme_gather.clh"
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

// gather using two grids
#define numGrids 2 /* Using define to avoid a conditional inside the function. */
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherKernelDual
#include "pme_gather.clh"
#undef numGrids
#undef CUSTOMIZED_KERNEL_NAME

#undef atomsPerBlock

/* SOLVE */

// GridOrdering enum replacement
#define YZX 1
#define XYZ 2

// solve, YZX dimension order state A
#define gridOrdering YZX
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXKernelA
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order state A
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernelA
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order state A
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernelA
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order state A
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernelA
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, YZX dimension order state B
#define gridOrdering YZX
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXKernelB
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order state B
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernelB
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order state B
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernelB
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order state B
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernelB
#include "pme_solve.clh"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME
