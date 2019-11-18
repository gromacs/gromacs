/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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

#ifndef GMX_NBNXM_PAIRLIST_H
#define GMX_NBNXM_PAIRLIST_H

#include <cstddef>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/locality.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/defaultinitializationallocator.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

// This file with constants is separate from this file to be able
// to include it during OpenCL jitting without including config.h
#include "constants.h"
#include "pairlistparams.h"

struct NbnxnPairlistCpuWork;
struct NbnxnPairlistGpuWork;


/* Convenience type for vector with aligned memory */
template<typename T>
using AlignedVector = std::vector<T, gmx::AlignedAllocator<T>>;

/* Convenience type for vector that avoids initialization at resize() */
template<typename T>
using FastVector = std::vector<T, gmx::DefaultInitializationAllocator<T>>;

/* A buffer data structure of 64 bytes
 * to be placed at the beginning and end of structs
 * to avoid cache invalidation of the real contents
 * of the struct by writes to neighboring memory.
 */
typedef struct
{
    int dummy[16];
} gmx_cache_protect_t;

/* This is the actual cluster-pair list j-entry.
 * cj is the j-cluster.
 * The interaction bits in excl are indexed i-major, j-minor.
 * The cj entries are sorted such that ones with exclusions come first.
 * This means that once a full mask (=NBNXN_INTERACTION_MASK_ALL)
 * is found, all subsequent j-entries in the i-entry also have full masks.
 */
struct nbnxn_cj_t
{
    int          cj;   /* The j-cluster                    */
    unsigned int excl; /* The exclusion (interaction) bits */
};

/* In nbnxn_ci_t the integer shift contains the shift in the lower 7 bits.
 * The upper bits contain information for non-bonded kernel optimization.
 * Simply calculating LJ and Coulomb for all pairs in a cluster pair is fine.
 * But three flags can be used to skip interactions, currently only for subc=0
 * !(shift & NBNXN_CI_DO_LJ(subc))   => we can skip LJ for all pairs
 * shift & NBNXN_CI_HALF_LJ(subc)    => we can skip LJ for the second half of i
 * !(shift & NBNXN_CI_DO_COUL(subc)) => we can skip Coulomb for all pairs
 */
#define NBNXN_CI_SHIFT 127
#define NBNXN_CI_DO_LJ(subc) (1 << (7 + 3 * (subc)))
#define NBNXN_CI_HALF_LJ(subc) (1 << (8 + 3 * (subc)))
#define NBNXN_CI_DO_COUL(subc) (1 << (9 + 3 * (subc)))

/* Cluster-pair Interaction masks
 * Bit i*j-cluster-size + j tells if atom i and j interact.
 */
// TODO: Rename according to convention when moving into Nbnxn namespace
/* All interaction mask is the same for all kernels */
constexpr unsigned int NBNXN_INTERACTION_MASK_ALL = 0xffffffffU;
/* 4x4 kernel diagonal mask */
constexpr unsigned int NBNXN_INTERACTION_MASK_DIAG = 0x08ceU;
/* 4x2 kernel diagonal masks */
constexpr unsigned int NBNXN_INTERACTION_MASK_DIAG_J2_0 = 0x0002U;
constexpr unsigned int NBNXN_INTERACTION_MASK_DIAG_J2_1 = 0x002fU;
/* 4x8 kernel diagonal masks */
constexpr unsigned int NBNXN_INTERACTION_MASK_DIAG_J8_0 = 0xf0f8fcfeU;
constexpr unsigned int NBNXN_INTERACTION_MASK_DIAG_J8_1 = 0x0080c0e0U;

/* Simple pair-list i-unit */
struct nbnxn_ci_t
{
    int ci;           /* i-cluster             */
    int shift;        /* Shift vector index plus possible flags, see above */
    int cj_ind_start; /* Start index into cj   */
    int cj_ind_end;   /* End index into cj     */
};

/* Grouped pair-list i-unit */
typedef struct
{
    /* Returns the number of j-cluster groups in this entry */
    int numJClusterGroups() const { return cj4_ind_end - cj4_ind_start; }

    int sci;           /* i-super-cluster       */
    int shift;         /* Shift vector index plus possible flags */
    int cj4_ind_start; /* Start index into cj4  */
    int cj4_ind_end;   /* End index into cj4    */
} nbnxn_sci_t;

/* Interaction data for a j-group for one warp */
struct nbnxn_im_ei_t
{
    // The i-cluster interactions mask for 1 warp
    unsigned int imask = 0U;
    // Index into the exclusion array for 1 warp, default index 0 which means no exclusions
    int excl_ind = 0;
};

typedef struct
{
    int           cj[c_nbnxnGpuJgroupSize];         /* The 4 j-clusters */
    nbnxn_im_ei_t imei[c_nbnxnGpuClusterpairSplit]; /* The i-cluster mask data       for 2 warps */
} nbnxn_cj4_t;

/* Struct for storing the atom-pair interaction bits for a cluster pair in a GPU pairlist */
struct nbnxn_excl_t
{
    /* Constructor, sets no exclusions, so all atom pairs interacting */
    nbnxn_excl_t()
    {
        for (unsigned int& pairEntry : pair)
        {
            pairEntry = NBNXN_INTERACTION_MASK_ALL;
        }
    }

    /* Topology exclusion interaction bits per warp */
    unsigned int pair[c_nbnxnGpuExclSize];
};

/* Cluster pairlist type for use on CPUs */
struct NbnxnPairlistCpu
{
    NbnxnPairlistCpu();

    gmx_cache_protect_t cp0;

    int                    na_ci;   /* The number of atoms per i-cluster        */
    int                    na_cj;   /* The number of atoms per j-cluster        */
    real                   rlist;   /* The radius for constructing the list     */
    FastVector<nbnxn_ci_t> ci;      /* The i-cluster list                       */
    FastVector<nbnxn_ci_t> ciOuter; /* The outer, unpruned i-cluster list       */

    FastVector<nbnxn_cj_t> cj;      /* The j-cluster list, size ncj             */
    FastVector<nbnxn_cj_t> cjOuter; /* The outer, unpruned j-cluster list       */
    int                    ncjInUse; /* The number of j-clusters that are used by ci entries in this list, will be <= cj.size() */

    int nci_tot; /* The total number of i clusters           */

    /* Working data storage for list construction */
    std::unique_ptr<NbnxnPairlistCpuWork> work;

    gmx_cache_protect_t cp1;
};

/* Cluster pairlist type, with extra hierarchies, for on the GPU
 *
 * NOTE: for better performance when combining lists over threads,
 *       all vectors should use default initialization. But when
 *       changing this, excl should be intialized when adding entries.
 */
struct NbnxnPairlistGpu
{
    /* Constructor
     *
     * \param[in] pinningPolicy  Sets the pinning policy for all buffers used on the GPU
     */
    NbnxnPairlistGpu(gmx::PinningPolicy pinningPolicy);

    gmx_cache_protect_t cp0;

    int  na_ci; /* The number of atoms per i-cluster        */
    int  na_cj; /* The number of atoms per j-cluster        */
    int  na_sc; /* The number of atoms per super cluster    */
    real rlist; /* The radius for constructing the list     */
    // The i-super-cluster list, indexes into cj4;
    gmx::HostVector<nbnxn_sci_t> sci;
    // The list of 4*j-cluster groups
    gmx::HostVector<nbnxn_cj4_t> cj4;
    // Atom interaction bits (non-exclusions)
    gmx::HostVector<nbnxn_excl_t> excl;
    // The total number of i-clusters
    int nci_tot;

    /* Working data storage for list construction */
    std::unique_ptr<NbnxnPairlistGpuWork> work;

    gmx_cache_protect_t cp1;
};

//! Initializes a free-energy pair-list
void nbnxn_init_pairlist_fep(t_nblist* nl);

#endif
