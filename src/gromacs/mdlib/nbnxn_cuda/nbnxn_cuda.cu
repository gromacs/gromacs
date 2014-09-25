/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <assert.h>

#if defined(_MSVC)
#include <limits>
#endif

#include <cuda.h>

#include "types/simple.h"
#include "types/nbnxn_pairlist.h"
#include "types/nb_verlet.h"
#include "types/ishift.h"
#include "types/force_flags.h"
#include "../nbnxn_consts.h"

#ifdef TMPI_ATOMICS
#include "thread_mpi/atomic.h"
#endif

#include "nbnxn_cuda_types.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"
#include "nbnxn_cuda.h"
#include "nbnxn_cuda_data_mgmt.h"

#if defined TEXOBJ_SUPPORTED && __CUDA_ARCH__ >= 300
#define USE_TEXOBJ
#endif

/*! Texture reference for LJ C6/C12 parameters; bound to cu_nbparam_t.nbfp */
texture<float, 1, cudaReadModeElementType> nbfp_texref;

/*! Texture reference for LJ-PME parameters; bound to cu_nbparam_t.nbfp_comb */
texture<float, 1, cudaReadModeElementType> nbfp_comb_texref;

/*! Texture reference for Ewald coulomb force table; bound to cu_nbparam_t.coulomb_tab */
texture<float, 1, cudaReadModeElementType> coulomb_tab_texref;

/* Convenience defines */
#define NCL_PER_SUPERCL         (NBNXN_GPU_NCLUSTER_PER_SUPERCLUSTER)
#define CL_SIZE                 (NBNXN_GPU_CLUSTER_SIZE)

/***** The kernels come here *****/
#include "nbnxn_cuda_kernel_utils.cuh"

/* Top-level kernel generation: will generate through multiple inclusion the
 * following flavors for all kernels:
 * - force-only output;
 * - force and energy output;
 * - force-only with pair list pruning;
 * - force and energy output with pair list pruning.
 */
/** Force only **/
#include "nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES

/*** Pair-list pruning kernels ***/
/** Force only **/
#define PRUNE_NBL
#include "nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL

/*! Nonbonded kernel function pointer type */
typedef void (*nbnxn_cu_kfunc_ptr_t)(const cu_atomdata_t,
                                     const cu_nbparam_t,
                                     const cu_plist_t,
                                     bool);

/*********************************/

/* XXX always/never run the energy/pruning kernels -- only for benchmarking purposes */
static bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);


/* Bit-pattern used for polling-based GPU synchronization. It is used as a float
 * and corresponds to having the exponent set to the maximum (127 -- single
 * precision) and the mantissa to 0.
 */
static unsigned int poll_wait_pattern = (0x7FU << 23);

/*! Returns the number of blocks to be used for the nonbonded GPU kernel. */
static inline int calc_nb_kernel_nblock(int nwork_units, cuda_dev_info_t *dinfo)
{
    int max_grid_x_size;

    assert(dinfo);

    max_grid_x_size = dinfo->prop.maxGridSize[0];

    /* do we exceed the grid x dimension limit? */
    if (nwork_units > max_grid_x_size)
    {
        gmx_fatal(FARGS, "Watch out, the input system is too large to simulate!\n"
                  "The number of nonbonded work units (=number of super-clusters) exceeds the"
                  "maximum grid size in x dimension (%d > %d)!", nwork_units, max_grid_x_size);
    }

    return nwork_units;
}


/* Constant arrays listing all kernel function pointers and enabling selection
   of a kernel in an elegant manner. */

/*! Pointers to the non-bonded kernels organized in 2-dim arrays by:
 *  electrostatics and VDW type.
 *
 *  Note that the row- and column-order of function pointers has to match the
 *  order of corresponding enumerated electrostatics and vdw types, resp.,
 *  defined in nbnxn_cuda_types.h.
 */

/*! Force-only kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_noener_noprune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_F_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_F_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_F_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_F_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_F_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_F_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_F_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_F_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_F_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_cuda      }
};

/*! Force + energy kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_noprune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_cuda              },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_cuda               },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_cuda          },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_cuda     },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_cuda               },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_cuda        }
};

/*! Force + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_noener_prune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_F_prune_cuda,             nbnxn_kernel_ElecCut_VdwLJFsw_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_prune_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_F_prune_cuda,              nbnxn_kernel_ElecRF_VdwLJFsw_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_prune_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_cuda,         nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_prune_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_prune_cuda,  nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_prune_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_F_prune_cuda,              nbnxn_kernel_ElecEw_VdwLJFsw_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_prune_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_cuda,       nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_prune_cuda      }
};

/*! Force + energy + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_prune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_prune_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_prune_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_prune_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_prune_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_prune_cuda      }
};

/*! Return a pointer to the kernel version to be executed at the current step. */
static inline nbnxn_cu_kfunc_ptr_t select_nbnxn_kernel(int  eeltype,
                                                       int  evdwtype,
                                                       bool bDoEne,
                                                       bool bDoPrune)
{
    nbnxn_cu_kfunc_ptr_t res;

    assert(eeltype < eelCuNR);
    assert(evdwtype < evdwCuNR);

    if (bDoEne)
    {
        if (bDoPrune)
        {
            res = nb_kfunc_ener_prune_ptr[eeltype][evdwtype];
        }
        else
        {
            res = nb_kfunc_ener_noprune_ptr[eeltype][evdwtype];
        }
    }
    else
    {
        if (bDoPrune)
        {
            res = nb_kfunc_noener_prune_ptr[eeltype][evdwtype];
        }
        else
        {
            res = nb_kfunc_noener_noprune_ptr[eeltype][evdwtype];
        }
    }

    return res;
}

/*! Calculates the amount of shared memory required by the CUDA kernel in use. */
static inline int calc_shmem_required()
{
    int shmem;

    /* size of shmem (force-buffers/xq/atom type preloading) */
    /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
    /* i-atom x+q in shared memory */
    shmem  = NCL_PER_SUPERCL * CL_SIZE * sizeof(float4);
    /* cj in shared memory, for both warps separately */
    shmem += 2 * NBNXN_GPU_JGROUP_SIZE * sizeof(int);
#ifdef IATYPE_SHMEM
    /* i-atom types in shared memory */
    shmem += NCL_PER_SUPERCL * CL_SIZE * sizeof(int);
#endif
#if __CUDA_ARCH__ < 300
    /* force reduction buffers in shared memory */
    shmem += CL_SIZE * CL_SIZE * 3 * sizeof(float);
#endif

    return shmem;
}

/*! As we execute nonbonded workload in separate streams, before launching
   the kernel we need to make sure that he following operations have completed:
   - atomdata allocation and related H2D transfers (every nstlist step);
   - pair list H2D transfer (every nstlist step);
   - shift vector H2D transfer (every nstlist step);
   - force (+shift force and energy) output clearing (every step).

   These operations are issued in the local stream at the beginning of the step
   and therefore always complete before the local kernel launch. The non-local
   kernel is launched after the local on the same device/context, so this is
   inherently scheduled after the operations in the local stream (including the
   above "misc_ops").
   However, for the sake of having a future-proof implementation, we use the
   misc_ops_done event to record the point in time when the above  operations
   are finished and synchronize with this event in the non-local stream.
 */
void nbnxn_cuda_launch_kernel(nbnxn_cuda_ptr_t        cu_nb,
                              const nbnxn_atomdata_t *nbatom,
                              int                     flags,
                              int                     iloc)
{
    cudaError_t          stat;
    int                  adat_begin, adat_len; /* local/nonlocal offset and length used for xq and f */
    /* CUDA kernel launch-related stuff */
    int                  shmem, nblock;
    dim3                 dim_block, dim_grid;
    nbnxn_cu_kfunc_ptr_t nb_kernel = NULL; /* fn pointer to the nonbonded kernel */

    cu_atomdata_t       *adat    = cu_nb->atdat;
    cu_nbparam_t        *nbp     = cu_nb->nbparam;
    cu_plist_t          *plist   = cu_nb->plist[iloc];
    cu_timers_t         *t       = cu_nb->timers;
    cudaStream_t         stream  = cu_nb->stream[iloc];

    bool                 bCalcEner   = flags & GMX_FORCE_VIRIAL;
    bool                 bCalcFshift = flags & GMX_FORCE_VIRIAL;
    bool                 bDoTime     = cu_nb->bDoTime;

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

    /* don't launch the kernel if there is no work to do */
    if (plist->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_I(iloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    /* When we get here all misc operations issues in the local stream are done,
       so we record that in the local stream and wait for it in the nonlocal one. */
    if (cu_nb->bUseTwoStreams)
    {
        if (iloc == eintLocal)
        {
            stat = cudaEventRecord(cu_nb->misc_ops_done, stream);
            CU_RET_ERR(stat, "cudaEventRecord on misc_ops_done failed");
        }
        else
        {
            stat = cudaStreamWaitEvent(stream, cu_nb->misc_ops_done, 0);
            CU_RET_ERR(stat, "cudaStreamWaitEvent on misc_ops_done failed");
        }
    }

    /* beginning of timed HtoD section */
    if (bDoTime)
    {
        stat = cudaEventRecord(t->start_nb_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* HtoD x, q */
    cu_copy_H2D_async(adat->xq + adat_begin, nbatom->x + adat_begin * 4,
                      adat_len * sizeof(*adat->xq), stream);

    if (bDoTime)
    {
        stat = cudaEventRecord(t->stop_nb_h2d[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* beginning of timed nonbonded calculation section */
    if (bDoTime)
    {
        stat = cudaEventRecord(t->start_nb_k[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* get the pointer to the kernel flavor we need to use */
    nb_kernel = select_nbnxn_kernel(nbp->eeltype,
                                    nbp->vdwtype,
                                    bCalcEner,
                                    plist->bDoPrune || always_prune);

    /* kernel launch config */
    nblock    = calc_nb_kernel_nblock(plist->nsci, cu_nb->dev_info);
    dim_block = dim3(CL_SIZE, CL_SIZE, 1);
    dim_grid  = dim3(nblock, 1, 1);
    shmem     = calc_shmem_required();

    if (debug)
    {
        fprintf(debug, "GPU launch configuration:\n\tThread block: %dx%dx%d\n\t"
                "Grid: %dx%d\n\t#Super-clusters/clusters: %d/%d (%d)\n",
                dim_block.x, dim_block.y, dim_block.z,
                dim_grid.x, dim_grid.y, plist->nsci*NCL_PER_SUPERCL,
                NCL_PER_SUPERCL, plist->na_c);
    }

    nb_kernel<<< dim_grid, dim_block, shmem, stream>>> (*adat, *nbp, *plist, bCalcFshift);
    CU_LAUNCH_ERR("k_calc_nb");

    if (bDoTime)
    {
        stat = cudaEventRecord(t->stop_nb_k[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

void nbnxn_cuda_launch_cpyback(nbnxn_cuda_ptr_t        cu_nb,
                               const nbnxn_atomdata_t *nbatom,
                               int                     flags,
                               int                     aloc)
{
    cudaError_t stat;
    int         adat_begin, adat_len, adat_end; /* local/nonlocal offset and length used for xq and f */
    int         iloc = -1;

    /* determine interaction locality from atom locality */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else
    {
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);
        gmx_incons(stmp);
    }

    cu_atomdata_t   *adat    = cu_nb->atdat;
    cu_timers_t     *t       = cu_nb->timers;
    bool             bDoTime = cu_nb->bDoTime;
    cudaStream_t     stream  = cu_nb->stream[iloc];

    bool             bCalcEner   = flags & GMX_FORCE_VIRIAL;
    bool             bCalcFshift = flags & GMX_FORCE_VIRIAL;

    /* don't launch copy-back if there was no work to do */
    if (cu_nb->plist[iloc]->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
        adat_end    = cu_nb->atdat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
        adat_end    = cu_nb->atdat->natoms;
    }

    /* beginning of timed D2H section */
    if (bDoTime)
    {
        stat = cudaEventRecord(t->start_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    if (!cu_nb->bUseStreamSync)
    {
        /* For safety reasons set a few (5%) forces to NaN. This way even if the
           polling "hack" fails with some future NVIDIA driver we'll get a crash. */
        for (int i = adat_begin; i < 3*adat_end + 2; i += adat_len/20)
        {
#ifdef NAN
            nbatom->out[0].f[i] = NAN;
#else
#  ifdef _MSVC
            if (numeric_limits<float>::has_quiet_NaN)
            {
                nbatom->out[0].f[i] = numeric_limits<float>::quiet_NaN();
            }
            else
#  endif
            {
                nbatom->out[0].f[i] = GMX_REAL_MAX;
            }
#endif
        }

        /* Set the last four bytes of the force array to a bit pattern
           which can't be the result of the force calculation:
           max exponent (127) and zero mantissa. */
        *(unsigned int*)&nbatom->out[0].f[adat_end*3 - 1] = poll_wait_pattern;
    }

    /* With DD the local D2H transfer can only start after the non-local
       has been launched. */
    if (iloc == eintLocal && cu_nb->bUseTwoStreams)
    {
        stat = cudaStreamWaitEvent(stream, cu_nb->nonlocal_done, 0);
        CU_RET_ERR(stat, "cudaStreamWaitEvent on nonlocal_done failed");
    }

    /* DtoH f */
    cu_copy_D2H_async(nbatom->out[0].f + adat_begin * 3, adat->f + adat_begin,
                      (adat_len)*sizeof(*adat->f), stream);

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we first need the non-local
       data back first. */
    if (iloc == eintNonlocal)
    {
        stat = cudaEventRecord(cu_nb->nonlocal_done, stream);
        CU_RET_ERR(stat, "cudaEventRecord on nonlocal_done failed");
    }

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH fshift */
        if (bCalcFshift)
        {
            cu_copy_D2H_async(cu_nb->nbst.fshift, adat->fshift,
                              SHIFTS * sizeof(*cu_nb->nbst.fshift), stream);
        }

        /* DtoH energies */
        if (bCalcEner)
        {
            cu_copy_D2H_async(cu_nb->nbst.e_lj, adat->e_lj,
                              sizeof(*cu_nb->nbst.e_lj), stream);
            cu_copy_D2H_async(cu_nb->nbst.e_el, adat->e_el,
                              sizeof(*cu_nb->nbst.e_el), stream);
        }
    }

    if (bDoTime)
    {
        stat = cudaEventRecord(t->stop_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

/* Atomic compare-exchange operation on unsigned values. It is used in
 * polling wait for the GPU.
 */
static inline bool atomic_cas(volatile unsigned int *ptr,
                              unsigned int           oldval,
                              unsigned int           newval)
{
    assert(ptr);

#ifdef TMPI_ATOMICS
    return tMPI_Atomic_cas((tMPI_Atomic_t *)ptr, oldval, newval);
#else
    gmx_incons("Atomic operations not available, atomic_cas() should not have been called!");
    return true;
#endif
}

void nbnxn_cuda_wait_gpu(nbnxn_cuda_ptr_t cu_nb,
                         const nbnxn_atomdata_t *nbatom,
                         int flags, int aloc,
                         real *e_lj, real *e_el, rvec *fshift)
{
    /* NOTE:  only implemented for single-precision at this time */
    cudaError_t            stat;
    int                    i, adat_end, iloc = -1;
    volatile unsigned int *poll_word;

    /* determine interaction locality from atom locality */
    if (LOCAL_A(aloc))
    {
        iloc = eintLocal;
    }
    else if (NONLOCAL_A(aloc))
    {
        iloc = eintNonlocal;
    }
    else
    {
        char stmp[STRLEN];
        sprintf(stmp, "Invalid atom locality passed (%d); valid here is only "
                "local (%d) or nonlocal (%d)", aloc, eatLocal, eatNonlocal);
        gmx_incons(stmp);
    }

    cu_plist_t      *plist    = cu_nb->plist[iloc];
    cu_timers_t     *timers   = cu_nb->timers;
    wallclock_gpu_t *timings  = cu_nb->timings;
    nb_staging       nbst     = cu_nb->nbst;

    bool             bCalcEner   = flags & GMX_FORCE_VIRIAL;
    bool             bCalcFshift = flags & GMX_FORCE_VIRIAL;

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

    /* don't launch wait/update timers & counters if there was no work to do

       NOTE: if timing with multiple GPUs (streams) becomes possible, the
       counters could end up being inconsistent due to not being incremented
       on some of the nodes! */
    if (cu_nb->plist[iloc]->nsci == 0)
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_end = cu_nb->atdat->natoms_local;
    }
    else
    {
        adat_end = cu_nb->atdat->natoms;
    }

    if (cu_nb->bUseStreamSync)
    {
        stat = cudaStreamSynchronize(cu_nb->stream[iloc]);
        CU_RET_ERR(stat, "cudaStreamSynchronize failed in cu_blockwait_nb");
    }
    else
    {
        /* Busy-wait until we get the signal pattern set in last byte
         * of the l/nl float vector. This pattern corresponds to a floating
         * point number which can't be the result of the force calculation
         * (maximum, 127 exponent and 0 mantissa).
         * The polling uses atomic compare-exchange.
         */
        poll_word = (volatile unsigned int*)&nbatom->out[0].f[adat_end*3 - 1];
        while (atomic_cas(poll_word, poll_wait_pattern, poll_wait_pattern))
        {
        }
    }

    /* timing data accumulation */
    if (cu_nb->bDoTime)
    {
        /* only increase counter once (at local F wait) */
        if (LOCAL_I(iloc))
        {
            timings->nb_c++;
            timings->ktime[plist->bDoPrune ? 1 : 0][bCalcEner ? 1 : 0].c += 1;
        }

        /* kernel timings */
        timings->ktime[plist->bDoPrune ? 1 : 0][bCalcEner ? 1 : 0].t +=
            cu_event_elapsed(timers->start_nb_k[iloc], timers->stop_nb_k[iloc]);

        /* X/q H2D and F D2H timings */
        timings->nb_h2d_t += cu_event_elapsed(timers->start_nb_h2d[iloc],
                                              timers->stop_nb_h2d[iloc]);
        timings->nb_d2h_t += cu_event_elapsed(timers->start_nb_d2h[iloc],
                                              timers->stop_nb_d2h[iloc]);

        /* only count atdat and pair-list H2D at pair-search step */
        if (plist->bDoPrune)
        {
            /* atdat transfer timing (add only once, at local F wait) */
            if (LOCAL_A(aloc))
            {
                timings->pl_h2d_c++;
                timings->pl_h2d_t += cu_event_elapsed(timers->start_atdat,
                                                      timers->stop_atdat);
            }

            timings->pl_h2d_t += cu_event_elapsed(timers->start_pl_h2d[iloc],
                                                  timers->stop_pl_h2d[iloc]);
        }
    }

    /* add up energies and shift forces (only once at local F wait) */
    if (LOCAL_I(iloc))
    {
        if (bCalcEner)
        {
            *e_lj += *nbst.e_lj;
            *e_el += *nbst.e_el;
        }

        if (bCalcFshift)
        {
            for (i = 0; i < SHIFTS; i++)
            {
                fshift[i][0] += nbst.fshift[i].x;
                fshift[i][1] += nbst.fshift[i].y;
                fshift[i][2] += nbst.fshift[i].z;
            }
        }
    }

    /* turn off pruning (doesn't matter if this is pair-search step or not) */
    plist->bDoPrune = false;
}

/*! Return the reference to the nbfp texture. */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_texref()
{
    return nbfp_texref;
}

/*! Return the reference to the nbfp_comb texture. */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_comb_texref()
{
    return nbfp_comb_texref;
}

/*! Return the reference to the coulomb_tab. */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_coulomb_tab_texref()
{
    return coulomb_tab_texref;
}

/*! Set up the cache configuration for the non-bonded kernels,
 */
void nbnxn_cuda_set_cacheconfig(cuda_dev_info_t *devinfo)
{
    cudaError_t stat;

    for (int i = 0; i < eelCuNR; i++)
    {
        for (int j = 0; j < evdwCuNR; j++)
        {
            if (devinfo->prop.major >= 3)
            {
                /* Default kernel on sm 3.x 48/16 kB Shared/L1 */
                stat = cudaFuncSetCacheConfig(nb_kfunc_ener_prune_ptr[i][j], cudaFuncCachePreferShared);
                stat = cudaFuncSetCacheConfig(nb_kfunc_ener_noprune_ptr[i][j], cudaFuncCachePreferShared);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_prune_ptr[i][j], cudaFuncCachePreferShared);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_noprune_ptr[i][j], cudaFuncCachePreferShared);
            }
            else
            {
                /* On Fermi prefer L1 gives 2% higher performance */
                /* Default kernel on sm_2.x 16/48 kB Shared/L1 */
                stat = cudaFuncSetCacheConfig(nb_kfunc_ener_prune_ptr[i][j], cudaFuncCachePreferL1);
                stat = cudaFuncSetCacheConfig(nb_kfunc_ener_noprune_ptr[i][j], cudaFuncCachePreferL1);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_prune_ptr[i][j], cudaFuncCachePreferL1);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_noprune_ptr[i][j], cudaFuncCachePreferL1);
            }
            CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
        }
    }
}
