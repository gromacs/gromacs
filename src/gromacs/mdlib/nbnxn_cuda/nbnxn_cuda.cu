/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \file
 *  \brief Define CUDA implementation of nbnxn_gpu.h
 *
 *  \author Szilard Pall <pall.szilard@gmail.com>
 */
#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <stdlib.h>

#include "gromacs/mdlib/nbnxn_gpu.h"

#if defined(_MSVC)
#include <limits>
#endif


#include "nbnxn_cuda.h"

#include "gromacs/gpu_utils/cudautils.cuh"
#include "gromacs/mdlib/force_flags.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_gpu_common.h"
#include "gromacs/mdlib/nbnxn_gpu_data_mgmt.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"
#include "gromacs/timing/gpu_timing.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/gmxassert.h"

#include "nbnxn_cuda_types.h"

/*
 * Texture references are created at compile-time and need to be declared
 * at file scope as global variables (see http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html#texture-reference-api).
 * The texture references below are used in two translation units;
 * we declare them here along the kernels that use them (when compiling legacy Fermi kernels),
 * and provide getters (see below) used by the data_mgmt module where the
 * textures are bound/unbound.
 * (In principle we could do it the other way arond, but that would likely require
 * device linking and we'd rather avoid technical hurdles.)
 */
/*! Texture reference for LJ C6/C12 parameters; bound to cu_nbparam_t.nbfp */
texture<float, 1, cudaReadModeElementType> nbfp_texref;

/*! Texture reference for LJ-PME parameters; bound to cu_nbparam_t.nbfp_comb */
texture<float, 1, cudaReadModeElementType> nbfp_comb_texref;

/*! Texture reference for Ewald coulomb force table; bound to cu_nbparam_t.coulomb_tab */
texture<float, 1, cudaReadModeElementType> coulomb_tab_texref;


/***** The kernel declarations/definitions come here *****/

/* Top-level kernel declaration generation: will generate through multiple
 * inclusion the following flavors for all kernel declarations:
 * - force-only output;
 * - force and energy output;
 * - force-only with pair list pruning;
 * - force and energy output with pair list pruning.
 */
#define FUNCTION_DECLARATION_ONLY
/** Force only **/
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES

/*** Pair-list pruning kernels ***/
/** Force only **/
#define PRUNE_NBL
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernels.cuh"
/** Force & energy **/
#define CALC_ENERGIES
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernels.cuh"
#undef CALC_ENERGIES
#undef PRUNE_NBL

/* Prune-only kernels */
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_pruneonly.cuh"
#undef FUNCTION_DECLARATION_ONLY

/* Now generate the function definitions if we are using a single compilation unit. */
#if GMX_CUDA_NB_SINGLE_COMPILATION_UNIT
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_F_noprune.cu"
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_F_prune.cu"
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_VF_noprune.cu"
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_VF_prune.cu"
#include "gromacs/mdlib/nbnxn_cuda/nbnxn_cuda_kernel_pruneonly.cu"
#else
/* Prevent compilation in multiple compilation unit mode for CC 2.x. Although we have
 * build-time checks to prevent this, the user could manually tweaks nvcc flags
 * which would lead to buggy kernels getting compiled.
 */
#if GMX_PTX_ARCH > 0 && GMX_PTX_ARCH <= 210
#error Due to an CUDA compiler bug, the CUDA non-bonded module can not be compiled with multiple compilation units for CC 2.x devices. If you have changed the nvcc flags manually, either use the GMX_CUDA_TARGET_* variables instead or set GMX_CUDA_NB_SINGLE_COMPILATION_UNIT=ON CMake option.
#endif
#endif /* GMX_CUDA_NB_SINGLE_COMPILATION_UNIT */


/*! Nonbonded kernel function pointer type */
typedef void (*nbnxn_cu_kfunc_ptr_t)(const cu_atomdata_t,
                                     const cu_nbparam_t,
                                     const cu_plist_t,
                                     bool);

/*********************************/

/* XXX switch between chevron and cudaLaunch (supported only in CUDA >=7.0)
   -- only for benchmarking purposes */
static const bool bUseCudaLaunchKernel =
    (GMX_CUDA_VERSION >= 7000) && (getenv("GMX_DISABLE_CUDALAUNCH") == NULL);

/* XXX always/never run the energy/pruning kernels -- only for benchmarking purposes */
static bool always_ener  = (getenv("GMX_GPU_ALWAYS_ENER") != NULL);
static bool never_ener   = (getenv("GMX_GPU_NEVER_ENER") != NULL);
static bool always_prune = (getenv("GMX_GPU_ALWAYS_PRUNE") != NULL);


/*! Returns the number of blocks to be used for the nonbonded GPU kernel. */
static inline int calc_nb_kernel_nblock(int nwork_units, const gmx_device_info_t *dinfo)
{
    int max_grid_x_size;

    assert(dinfo);
    /* CUDA does not accept grid dimension of 0 (which can happen e.g. with an
       empty domain) and that case should be handled before this point. */
    assert(nwork_units > 0);

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
    { nbnxn_kernel_ElecCut_VdwLJ_F_cuda,            nbnxn_kernel_ElecCut_VdwLJCombGeom_F_cuda,            nbnxn_kernel_ElecCut_VdwLJCombLB_F_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_F_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_F_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_F_cuda,             nbnxn_kernel_ElecRF_VdwLJCombGeom_F_cuda,             nbnxn_kernel_ElecRF_VdwLJCombLB_F_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_F_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_F_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_F_cuda,             nbnxn_kernel_ElecEw_VdwLJCombGeom_F_cuda,             nbnxn_kernel_ElecEw_VdwLJCombLB_F_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_F_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_F_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_cuda      }
};

/*! Force + energy kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_noprune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJCombLB_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJCombLB_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJCombLB_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_cuda      }
};

/*! Force + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_noener_prune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJCombGeom_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJCombLB_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_F_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_F_prune_cuda             },
    { nbnxn_kernel_ElecRF_VdwLJ_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJCombGeom_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJCombLB_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_F_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_F_prune_cuda              },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_F_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_F_prune_cuda         },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_F_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_F_prune_cuda  },
    { nbnxn_kernel_ElecEw_VdwLJ_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJCombGeom_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJCombLB_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_F_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_F_prune_cuda              },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_F_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_F_prune_cuda       }
};

/*! Force + energy + pruning kernel function pointers. */
static const nbnxn_cu_kfunc_ptr_t nb_kfunc_ener_prune_ptr[eelCuNR][evdwCuNR] =
{
    { nbnxn_kernel_ElecCut_VdwLJ_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJCombGeom_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJCombLB_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJFsw_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJPsw_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombGeom_VF_prune_cuda,            nbnxn_kernel_ElecCut_VdwLJEwCombLB_VF_prune_cuda            },
    { nbnxn_kernel_ElecRF_VdwLJ_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJCombLB_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJFsw_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJPsw_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecRF_VdwLJEwCombLB_VF_prune_cuda             },
    { nbnxn_kernel_ElecEwQSTab_VdwLJ_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombGeom_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJCombLB_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJFsw_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJPsw_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombGeom_VF_prune_cuda,        nbnxn_kernel_ElecEwQSTab_VdwLJEwCombLB_VF_prune_cuda        },
    { nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJ_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombGeom_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJCombLB_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJFsw_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJPsw_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombGeom_VF_prune_cuda, nbnxn_kernel_ElecEwQSTabTwinCut_VdwLJEwCombLB_VF_prune_cuda },
    { nbnxn_kernel_ElecEw_VdwLJ_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJCombLB_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJFsw_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJPsw_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombGeom_VF_prune_cuda,             nbnxn_kernel_ElecEw_VdwLJEwCombLB_VF_prune_cuda             },
    { nbnxn_kernel_ElecEwTwinCut_VdwLJ_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombGeom_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJCombLB_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJFsw_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJPsw_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombGeom_VF_prune_cuda,      nbnxn_kernel_ElecEwTwinCut_VdwLJEwCombLB_VF_prune_cuda      }
};

/*! Return a pointer to the kernel version to be executed at the current step. */
static inline nbnxn_cu_kfunc_ptr_t select_nbnxn_kernel(int                                  eeltype,
                                                       int                                  evdwtype,
                                                       bool                                 bDoEne,
                                                       bool                                 bDoPrune,
                                                       const gmx_device_info_t gmx_unused  *devInfo)
{
    nbnxn_cu_kfunc_ptr_t res;

    GMX_ASSERT(eeltype < eelCuNR,
               "The electrostatics type requested is not implemented in the CUDA kernels.");
    GMX_ASSERT(evdwtype < evdwCuNR,
               "The VdW type requested is not implemented in the CUDA kernels.");

    /* assert assumptions made by the kernels */
    GMX_ASSERT(c_nbnxnGpuClusterSize*c_nbnxnGpuClusterSize/c_nbnxnGpuClusterpairSplit == devInfo->prop.warpSize,
               "The CUDA kernels require the cluster_size_i*cluster_size_j/nbnxn_gpu_clusterpair_split to match the warp size of the architecture targeted.");

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

/*! \brief Calculates the amount of shared memory required by the nonbonded kernel in use. */
static inline int calc_shmem_required_nonbonded(const int num_threads_z, const gmx_device_info_t gmx_unused *dinfo, const cu_nbparam_t *nbp)
{
    int shmem;

    assert(dinfo);

    /* size of shmem (force-buffers/xq/atom type preloading) */
    /* NOTE: with the default kernel on sm3.0 we need shmem only for pre-loading */
    /* i-atom x+q in shared memory */
    shmem  = c_numClPerSupercl * c_clSize * sizeof(float4);
    /* cj in shared memory, for each warp separately */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);
    if (dinfo->prop.major >= 3)
    {
        if (nbp->vdwtype == evdwCuCUTCOMBGEOM ||
            nbp->vdwtype == evdwCuCUTCOMBLB)
        {
            /* i-atom LJ combination parameters in shared memory */
            shmem += c_numClPerSupercl * c_clSize * sizeof(float2);
        }
        else
        {
            /* i-atom types in shared memory */
            shmem += c_numClPerSupercl * c_clSize * sizeof(int);
        }
    }
    if (dinfo->prop.major < 3)
    {
        /* force reduction buffers in shared memory */
        shmem += c_clSize * c_clSize * 3 * sizeof(float);
    }
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
   kernel is launched after the local on the same device/context hence it is
   inherently scheduled after the operations in the local stream (including the
   above "misc_ops") on pre-GK110 devices with single hardware queue, but on later
   devices with multiple hardware queues the dependency needs to be enforced.
   We use the misc_ops_and_local_H2D_done event to record the point where
   the local x+q H2D (and all preceding) tasks are complete and synchronize
   with this event in the non-local stream before launching the non-bonded kernel.
 */
void nbnxn_gpu_launch_kernel(gmx_nbnxn_cuda_t       *nb,
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

    cu_atomdata_t       *adat    = nb->atdat;
    cu_nbparam_t        *nbp     = nb->nbparam;
    cu_plist_t          *plist   = nb->plist[iloc];
    cu_timers_t         *t       = nb->timers;
    cudaStream_t         stream  = nb->stream[iloc];

    bool                 bCalcEner   = flags & GMX_FORCE_ENERGY;
    bool                 bCalcFshift = flags & GMX_FORCE_VIRIAL;
    bool                 bDoTime     = nb->bDoTime;

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

    /* Don't launch the non-local kernel if there is no work to do.
       Doing the same for the local kernel is more complicated, since the
       local part of the force array also depends on the non-local kernel.
       So to avoid complicating the code and to reduce the risk of bugs,
       we always call the local kernel, the local x+q copy and later (not in
       this function) the stream wait, local f copyback and the f buffer
       clearing. All these operations, except for the local interaction kernel,
       are needed for the non-local interactions. The skip of the local kernel
       call is taken care of later in this function. */
    if (canSkipWork(nb, iloc))
    {
        plist->haveFreshList = false;

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

    /* When we get here all misc operations issues in the local stream as well as
       the local xq H2D are done,
       so we record that in the local stream and wait for it in the nonlocal one. */
    if (nb->bUseTwoStreams)
    {
        if (iloc == eintLocal)
        {
            stat = cudaEventRecord(nb->misc_ops_and_local_H2D_done, stream);
            CU_RET_ERR(stat, "cudaEventRecord on misc_ops_and_local_H2D_done failed");
        }
        else
        {
            stat = cudaStreamWaitEvent(stream, nb->misc_ops_and_local_H2D_done, 0);
            CU_RET_ERR(stat, "cudaStreamWaitEvent on misc_ops_and_local_H2D_done failed");
        }
    }

    if (nbp->useDynamicPruning && plist->haveFreshList)
    {
        /* Prunes for rlistOuter and rlistInner, sets plist->haveFreshList=false
           (TODO: ATM that's the way the timing accounting can distinguish between
           separate prune kernel and combined force+prune, maybe we need a better way?).
         */
        nbnxn_gpu_launch_kernel_pruneonly(nb, iloc, 1);
    }

    if (plist->nsci == 0)
    {
        /* Don't launch an empty local kernel (not allowed with CUDA) */
        return;
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
                                    (plist->haveFreshList && !nb->timers->didPrune[iloc]) || always_prune,
                                    nb->dev_info);

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int num_threads_z = 1;
    if (nb->dev_info->prop.major == 3 && nb->dev_info->prop.minor == 7)
    {
        num_threads_z = 2;
    }
    nblock    = calc_nb_kernel_nblock(plist->nsci, nb->dev_info);
    dim_block = dim3(c_clSize, c_clSize, num_threads_z);
    dim_grid  = dim3(nblock, 1, 1);
    shmem     = calc_shmem_required_nonbonded(num_threads_z, nb->dev_info, nbp);

    if (debug)
    {
        fprintf(debug, "Non-bonded GPU launch configuration:\n\tThread block: %ux%ux%u\n\t"
                "\tGrid: %ux%u\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %d\n",
                dim_block.x, dim_block.y, dim_block.z,
                dim_grid.x, dim_grid.y, plist->nsci*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c,
                shmem);
    }

    if (bUseCudaLaunchKernel)
    {
        gmx_unused void* kernel_args[4];
        kernel_args[0] = adat;
        kernel_args[1] = nbp;
        kernel_args[2] = plist;
        kernel_args[3] = &bCalcFshift;

#if GMX_CUDA_VERSION >= 7000
        cudaLaunchKernel((void *)nb_kernel, dim_grid, dim_block, kernel_args, shmem, stream);
#endif
    }
    else
    {
        nb_kernel<<< dim_grid, dim_block, shmem, stream>>> (*adat, *nbp, *plist, bCalcFshift);
    }
    CU_LAUNCH_ERR("k_calc_nb");

    if (bDoTime)
    {
        stat = cudaEventRecord(t->stop_nb_k[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

#if (defined(WIN32) || defined( _WIN32 ))
    /* Windows: force flushing WDDM queue */
    stat = cudaStreamQuery(stream);
#endif
}

/*! Calculates the amount of shared memory required by the CUDA kernel in use. */
static inline int calc_shmem_required_prune(const int num_threads_z)
{
    int shmem;

    /* i-atom x in shared memory */
    shmem  = c_numClPerSupercl * c_clSize * sizeof(float4);
    /* cj in shared memory, for each warp separately */
    shmem += num_threads_z * c_nbnxnGpuClusterpairSplit * c_nbnxnGpuJgroupSize * sizeof(int);

    return shmem;
}

void nbnxn_gpu_launch_kernel_pruneonly(gmx_nbnxn_cuda_t       *nb,
                                       int                     iloc,
                                       int                     numParts)
{
    cudaError_t          stat;

    cu_atomdata_t       *adat    = nb->atdat;
    cu_nbparam_t        *nbp     = nb->nbparam;
    cu_plist_t          *plist   = nb->plist[iloc];
    cu_timers_t         *t       = nb->timers;
    cudaStream_t         stream  = nb->stream[iloc];

    bool                 bDoTime = nb->bDoTime;

    if (plist->haveFreshList)
    {
        GMX_ASSERT(numParts == 1, "With first pruning we expect 1 part");

        /* Set rollingPruningNumParts to signal that it is not set */
        plist->rollingPruningNumParts = 0;
        plist->rollingPruningPart     = 0;
    }
    else
    {
        if (plist->rollingPruningNumParts == 0)
        {
            plist->rollingPruningNumParts = numParts;
        }
        else
        {
            GMX_ASSERT(numParts == plist->rollingPruningNumParts, "It is not allowed to change numParts in between list generation steps");
        }
    }

    /* Use a local variable for part and update in plist, so we can return here
     * without duplicating the part increment code.
     */
    int part = plist->rollingPruningPart;

    plist->rollingPruningPart++;
    if (plist->rollingPruningPart >= plist->rollingPruningNumParts)
    {
        plist->rollingPruningPart = 0;
    }

    /* Compute the number of list entries to prune in this pass */
    int numSciInPart = (plist->nsci - part)/numParts;

    /* Don't launch the kernel if there is no work to do (not allowed with CUDA) */
    if (numSciInPart <= 0)
    {
        plist->haveFreshList = false;

        return;
    }

    cudaEvent_t startEvent, stopEvent;
    if (bDoTime)
    {
        startEvent = (plist->haveFreshList ? t->start_prune_k[iloc] : t->start_rollingPrune_k[iloc]);
        stopEvent  = (plist->haveFreshList ? t->stop_prune_k[iloc]  : t->stop_rollingPrune_k[iloc]);
    }

    /* beginning of timed prune calculation section */
    if (bDoTime)
    {
        stat = cudaEventRecord(startEvent, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* Kernel launch config:
     * - The thread block dimensions match the size of i-clusters, j-clusters,
     *   and j-cluster concurrency, in x, y, and z, respectively.
     * - The 1D block-grid contains as many blocks as super-clusters.
     */
    int  num_threads_z  = c_cudaPruneKernelJ4Concurrency;
    int  nblock         = calc_nb_kernel_nblock(numSciInPart, nb->dev_info);
    dim3 dim_block      = dim3(c_clSize, c_clSize, num_threads_z);
    dim3 dim_grid       = dim3(nblock, 1, 1);
    int  shmem          = calc_shmem_required_prune(num_threads_z);

    if (debug)
    {
        fprintf(debug, "Pruning GPU kernel launch configuration:\n\tThread block: %dx%dx%d\n\t"
                "\tGrid: %dx%d\n\t#Super-clusters/clusters: %d/%d (%d)\n"
                "\tShMem: %d\n",
                dim_block.x, dim_block.y, dim_block.z,
                dim_grid.x, dim_grid.y, numSciInPart*c_numClPerSupercl,
                c_numClPerSupercl, plist->na_c,
                shmem);
    }

    if (bUseCudaLaunchKernel)
    {
        gmx_unused void* kernel_args[5];
        kernel_args[0] = adat;
        kernel_args[1] = nbp;
        kernel_args[2] = plist;
        kernel_args[3] = &numParts;
        kernel_args[4] = &part;

#if GMX_CUDA_VERSION >= 7000
        if (plist->haveFreshList)
        {
            cudaLaunchKernel((void *)nbnxn_kernel_prune_cuda<true>, dim_grid, dim_block, kernel_args, shmem, stream);
        }
        else
        {
            cudaLaunchKernel((void *)nbnxn_kernel_prune_cuda<false>, dim_grid, dim_block, kernel_args, shmem, stream);
        }
#endif
    }
    else
    {
        if (plist->haveFreshList)
        {
            nbnxn_kernel_prune_cuda<true><<< dim_grid, dim_block, shmem, stream>>> (*adat, *nbp, *plist, numParts, part);
        }
        else
        {
            nbnxn_kernel_prune_cuda<false><<< dim_grid, dim_block, shmem, stream>>> (*adat, *nbp, *plist, numParts, part);
        }
    }
    CU_LAUNCH_ERR("k_pruneonly");

    /* TODO: consider a more elegant way to track which kernel has been called
       (combined or separate 1st pass prune, rolling prune). */
    if (plist->haveFreshList)
    {
        plist->haveFreshList         = false;
        /* Mark that pruning has been done */
        nb->timers->didPrune[iloc] = true;
    }
    else
    {
        /* Mark that rolling pruning has been done */
        nb->timers->didRollingPrune[iloc] = true;
    }

    if (bDoTime)
    {
        stat = cudaEventRecord(stopEvent, stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

#if (defined(WIN32) || defined( _WIN32 ))
    /* Windows: force flushing WDDM queue */
    stat = cudaStreamQuery(stream);
#endif
}

void nbnxn_gpu_launch_cpyback(gmx_nbnxn_cuda_t       *nb,
                              const nbnxn_atomdata_t *nbatom,
                              int                     flags,
                              int                     aloc)
{
    cudaError_t stat;
    int         adat_begin, adat_len; /* local/nonlocal offset and length used for xq and f */
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

    cu_atomdata_t   *adat    = nb->atdat;
    cu_timers_t     *t       = nb->timers;
    bool             bDoTime = nb->bDoTime;
    cudaStream_t     stream  = nb->stream[iloc];

    bool             bCalcEner   = flags & GMX_FORCE_ENERGY;
    bool             bCalcFshift = flags & GMX_FORCE_VIRIAL;

    /* don't launch non-local copy-back if there was no non-local work to do */
    if (canSkipWork(nb, iloc))
    {
        return;
    }

    /* calculate the atom data index range based on locality */
    if (LOCAL_A(aloc))
    {
        adat_begin  = 0;
        adat_len    = adat->natoms_local;
    }
    else
    {
        adat_begin  = adat->natoms_local;
        adat_len    = adat->natoms - adat->natoms_local;
    }

    /* beginning of timed D2H section */
    if (bDoTime)
    {
        stat = cudaEventRecord(t->start_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }

    /* With DD the local D2H transfer can only start after the non-local
       kernel has finished. */
    if (iloc == eintLocal && nb->bUseTwoStreams)
    {
        stat = cudaStreamWaitEvent(stream, nb->nonlocal_done, 0);
        CU_RET_ERR(stat, "cudaStreamWaitEvent on nonlocal_done failed");
    }

    /* DtoH f */
    cu_copy_D2H_async(nbatom->out[0].f + adat_begin * 3, adat->f + adat_begin,
                      (adat_len)*sizeof(*adat->f), stream);

    /* After the non-local D2H is launched the nonlocal_done event can be
       recorded which signals that the local D2H can proceed. This event is not
       placed after the non-local kernel because we want the non-local data
       back first. */
    if (iloc == eintNonlocal)
    {
        stat = cudaEventRecord(nb->nonlocal_done, stream);
        CU_RET_ERR(stat, "cudaEventRecord on nonlocal_done failed");
    }

    /* only transfer energies in the local stream */
    if (LOCAL_I(iloc))
    {
        /* DtoH fshift */
        if (bCalcFshift)
        {
            cu_copy_D2H_async(nb->nbst.fshift, adat->fshift,
                              SHIFTS * sizeof(*nb->nbst.fshift), stream);
        }

        /* DtoH energies */
        if (bCalcEner)
        {
            cu_copy_D2H_async(nb->nbst.e_lj, adat->e_lj,
                              sizeof(*nb->nbst.e_lj), stream);
            cu_copy_D2H_async(nb->nbst.e_el, adat->e_el,
                              sizeof(*nb->nbst.e_el), stream);
        }
    }

    if (bDoTime)
    {
        stat = cudaEventRecord(t->stop_nb_d2h[iloc], stream);
        CU_RET_ERR(stat, "cudaEventRecord failed");
    }
}

/*! \brief Count pruning kernel time if either kernel has been triggered
 *
 *  We do the accounting for either of the two pruning kernel flavors:
 *   - 1st pass prune: ran during the current step (prior to the force kernel);
 *   - rolling prune:  ran at the end of the previous step (prior to the current step H2D xq);
 *
 * Note that the resetting of cu_timers_t::didPrune and cu_timers_t::didRollingPrune should happen
 * after calling this function.
 *
 * \param[in] timers      structs with CUDA timer objects
 * \param[inout] timings  GPU task timing data
 * \param[in] iloc        interaction locality
 */
static void countPruneKernelTime(const cu_timers_t   *timers,
                                 gmx_wallclock_gpu_t *timings,
                                 const int            iloc)
{
    // We might have not done any pruning (e.g. if we skipped with empty domains).
    if (!timers->didPrune[iloc] && !timers->didRollingPrune[iloc])
    {
        return;
    }

    if (timers->didPrune[iloc])
    {
        timings->pruneTime.c++;
        timings->pruneTime.t += cu_event_elapsed(timers->start_prune_k[iloc],
                                                 timers->stop_prune_k[iloc]);
    }
    if (timers->didRollingPrune[iloc])
    {
        timings->dynamicPruneTime.c++;
        timings->dynamicPruneTime.t += cu_event_elapsed(timers->start_rollingPrune_k[iloc],
                                                        timers->stop_rollingPrune_k[iloc]);
    }
}

void nbnxn_gpu_wait_for_gpu(gmx_nbnxn_cuda_t *nb,
                            int flags, int aloc,
                            real *e_lj, real *e_el, rvec *fshift)
{
    /* NOTE:  only implemented for single-precision at this time */
    cudaError_t stat;
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

    cu_plist_t                 *plist    = nb->plist[iloc];
    cu_timers_t                *timers   = nb->timers;
    struct gmx_wallclock_gpu_t *timings  = nb->timings;
    nb_staging                  nbst     = nb->nbst;

    bool                        bCalcEner   = flags & GMX_FORCE_ENERGY;
    bool                        bCalcFshift = flags & GMX_FORCE_VIRIAL;

    /* turn energy calculation always on/off (for debugging/testing only) */
    bCalcEner = (bCalcEner || always_ener) && !never_ener;

    /* Launch wait/update timers & counters and do reduction into staging buffers
       BUT skip it when during the non-local phase there was actually no work to do.
       This is consistent with nbnxn_gpu_launch_kernel.

       NOTE: if timing with multiple GPUs (streams) becomes possible, the
       counters could end up being inconsistent due to not being incremented
       on some of the nodes! */
    if (!canSkipWork(nb, iloc))
    {
        stat = cudaStreamSynchronize(nb->stream[iloc]);
        CU_RET_ERR(stat, "cudaStreamSynchronize failed in cu_blockwait_nb");

        /* timing data accumulation */
        if (nb->bDoTime)
        {
            /* only increase counter once (at local F wait) */
            if (LOCAL_I(iloc))
            {
                timings->nb_c++;
                timings->ktime[plist->haveFreshList ? 1 : 0][bCalcEner ? 1 : 0].c += 1;
            }

            /* kernel timings */
            timings->ktime[plist->haveFreshList ? 1 : 0][bCalcEner ? 1 : 0].t +=
                cu_event_elapsed(timers->start_nb_k[iloc], timers->stop_nb_k[iloc]);

            /* X/q H2D and F D2H timings */
            timings->nb_h2d_t += cu_event_elapsed(timers->start_nb_h2d[iloc],
                                                  timers->stop_nb_h2d[iloc]);
            timings->nb_d2h_t += cu_event_elapsed(timers->start_nb_d2h[iloc],
                                                  timers->stop_nb_d2h[iloc]);

            /* Count the pruning kernel times for both cases:1st pass (at search step)
               and rolling pruning (if called at the previous step).
               We do the accounting here as this is the only sync point where we
               know (without checking or additional sync-ing) that prune tasks in
               in the current stream have completed (having just blocking-waited
               for the force D2H). */
            countPruneKernelTime(timers, timings, iloc);

            /* only count atdat and pair-list H2D at pair-search step */
            if (timers->didPairlistH2D[iloc])
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

                /* Clear the timing flag for the next step */
                timers->didPairlistH2D[iloc] = false;
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
                for (int i = 0; i < SHIFTS; i++)
                {
                    fshift[i][0] += nbst.fshift[i].x;
                    fshift[i][1] += nbst.fshift[i].y;
                    fshift[i][2] += nbst.fshift[i].z;
                }
            }
        }
    }

    /* Always reset both pruning flags (doesn't hurt doing it even when timing is off). */
    timers->didPrune[iloc] = timers->didRollingPrune[iloc] = false;

    /* Turn off initial list pruning (doesn't hurt if this is not pair-search step). */
    plist->haveFreshList = false;
}

const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_texref()
{
    assert(!c_disableCudaTextures);
    return nbfp_texref;
}

const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_comb_texref()
{
    assert(!c_disableCudaTextures);
    return nbfp_comb_texref;
}

const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_coulomb_tab_texref()
{
    assert(!c_disableCudaTextures);
    return coulomb_tab_texref;
}

void nbnxn_cuda_set_cacheconfig(const gmx_device_info_t *devinfo)
{
    cudaError_t stat;

    for (int i = 0; i < eelCuNR; i++)
    {
        for (int j = 0; j < evdwCuNR; j++)
        {
            if (devinfo->prop.major >= 3)
            {
                /* Default kernel on sm 3.x and later 32/32 kB Shared/L1 */
                cudaFuncSetCacheConfig(nb_kfunc_ener_prune_ptr[i][j], cudaFuncCachePreferEqual);
                cudaFuncSetCacheConfig(nb_kfunc_ener_noprune_ptr[i][j], cudaFuncCachePreferEqual);
                cudaFuncSetCacheConfig(nb_kfunc_noener_prune_ptr[i][j], cudaFuncCachePreferEqual);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_noprune_ptr[i][j], cudaFuncCachePreferEqual);
            }
            else
            {
                /* On Fermi prefer L1 gives 2% higher performance */
                /* Default kernel on sm_2.x 16/48 kB Shared/L1 */
                cudaFuncSetCacheConfig(nb_kfunc_ener_prune_ptr[i][j], cudaFuncCachePreferL1);
                cudaFuncSetCacheConfig(nb_kfunc_ener_noprune_ptr[i][j], cudaFuncCachePreferL1);
                cudaFuncSetCacheConfig(nb_kfunc_noener_prune_ptr[i][j], cudaFuncCachePreferL1);
                stat = cudaFuncSetCacheConfig(nb_kfunc_noener_noprune_ptr[i][j], cudaFuncCachePreferL1);
            }
            CU_RET_ERR(stat, "cudaFuncSetCacheConfig failed");
        }
    }
}
