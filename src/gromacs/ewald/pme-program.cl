//This is a top-level file to generate all PME openCL kernels

//FIXME #define CUSTOMIZED_KERNEL_NAME(x) x for CUDA

/* SPREAD/SPLINE */


#define c_spreadMaxWarpsPerBlock 8
#define c_spreadMaxThreadsPerBlock (c_spreadMaxWarpsPerBlock * warp_size)
#define atomsPerBlock  (c_spreadMaxThreadsPerBlock / PME_SPREADGATHER_THREADS_PER_ATOM)


// splineAndSpread
#define computeSplines 1
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineAndSpreadKernel
#include "../../ewald/pme-spread-kernel.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spline
#define computeSplines 1
#define spreadCharges 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSplineKernel
#include "../../ewald/pme-spread-kernel.cl"
#undef computeSplines
#undef spreadCharges
#undef CUSTOMIZED_KERNEL_NAME

// spread
#define computeSplines 0
#define spreadCharges 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSpreadKernel
#include "../../ewald/pme-spread-kernel.cl"
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
#include "../../ewald/pme-gather-kernel.cl"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME

// gather with reduction
#define overwriteForces 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeGatherReduceWithInputKernel
#include "../../ewald/pme-gather-kernel.cl"
#undef overwriteForces
#undef CUSTOMIZED_KERNEL_NAME


/* SOLVE */
//test
//#undef warp_size
#include "../../ewald/pme-ocl-definitely-common.h"

//GRID orderings
#define YZX 1
#define XYZ 2

// solve, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXKernel
#include "../../ewald/pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, YZX dimension order
#define gridOrdering YZX
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveYZXEnergyKernel
#include "../../ewald/pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 0
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZKernel
#include "../../ewald/pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

// solve with reduction, XYZ dimension order
#define gridOrdering XYZ
#define computeEnergyAndVirial 1
#define CUSTOMIZED_KERNEL_NAME(x) pmeSolveXYZEnergyKernel
#include "../../ewald/pme-solve-kernel.cl"
#undef gridOrdering
#undef computeEnergyAndVirial
#undef CUSTOMIZED_KERNEL_NAME

#undef atomsPerBlock
