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
