/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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

/*! \internal \file
 *  \brief
 *  Data types used internally in the nbnxn_cuda module.
 *
 *  \author Szilárd Páll <pall.szilard@gmail.com>
 *  \ingroup module_mdlib
 */

#ifndef NBNXN_CUDA_TYPES_H
#define NBNXN_CUDA_TYPES_H

#include "config.h"

#include "gromacs/gmxlib/cuda_tools/cudautils.cuh"
#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/mdlib/nbnxn_pairlist.h"

#ifndef HAVE_CUDA_TEXOBJ_SUPPORT
/** This typedef allows us to define only one version of struct cu_nbparam */
typedef int cudaTextureObject_t;
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Electrostatic CUDA kernel flavors.
 *
 *  Types of electrostatics implementations available in the CUDA non-bonded
 *  force kernels. These represent both the electrostatics types implemented
 *  by the kernels (cut-off, RF, and Ewald - a subset of what's defined in
 *  enums.h) as well as encode implementation details analytical/tabulated
 *  and single or twin cut-off (for Ewald kernels).
 *  Note that the cut-off and RF kernels have only analytical flavor and unlike
 *  in the CPU kernels, the tabulated kernels are ATM Ewald-only.
 *
 *  The row-order of pointers to different electrostatic kernels defined in
 *  nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 *  should match the order of enumerated types below.
 */
enum eelCu {
    eelCuCUT, eelCuRF, eelCuEWALD_TAB, eelCuEWALD_TAB_TWIN, eelCuEWALD_ANA, eelCuEWALD_ANA_TWIN, eelCuNR
};

/*! \brief VdW CUDA kernel flavors.
 *
 * The enumerates values correspond to the LJ implementations in the CUDA non-bonded
 * kernels.
 *
 * The column-order of pointers to different electrostatic kernels defined in
 * nbnxn_cuda.cu by the nb_*_kfunc_ptr function pointer table
 * should match the order of enumerated types below.
 */
enum evdwCu {
    evdwCuCUT, evdwCuFSWITCH, evdwCuPSWITCH, evdwCuEWALDGEOM, evdwCuEWALDLB, evdwCuNR
};

/* All structs prefixed with "cu_" hold data used in GPU calculations and
 * are passed to the kernels, except cu_timers_t. */
/*! \cond */
typedef struct cu_plist     cu_plist_t;
typedef struct cu_atomdata  cu_atomdata_t;
typedef struct cu_nbparam   cu_nbparam_t;
typedef struct cu_timers    cu_timers_t;
typedef struct nb_staging   nb_staging_t;
/*! \endcond */


/** \internal
 * \brief Staging area for temporary data downloaded from the GPU.
 *
 *  The energies/shift forces get downloaded here first, before getting added
 *  to the CPU-side aggregate values.
 */
struct nb_staging
{
    float   *e_lj;      /**< LJ energy            */
    float   *e_el;      /**< electrostatic energy */
    float3  *fshift;    /**< shift forces         */
};

/** \internal
 * \brief Nonbonded atom data - both inputs and outputs.
 */
struct cu_atomdata
{
    int      natoms;            /**< number of atoms                              */
    int      natoms_local;      /**< number of local atoms                        */
    int      nalloc;            /**< allocation size for the atom data (xq, f)    */

    float4  *xq;                /**< atom coordinates + charges, size natoms      */
    float3  *f;                 /**< force output array, size natoms              */

    float   *e_lj;              /**< LJ energy output, size 1                     */
    float   *e_el;              /**< Electrostatics energy input, size 1          */

    float3  *fshift;            /**< shift forces                                 */

    int      ntypes;            /**< number of atom types                         */
    int     *atom_types;        /**< atom type indices, size natoms               */

    float3  *shift_vec;         /**< shifts                                       */
    bool     bShiftVecUploaded; /**< true if the shift vector has been uploaded   */
};

/** \internal
 * \brief Parameters required for the CUDA nonbonded calculations.
 */
struct cu_nbparam
{

    int             eeltype;          /**< type of electrostatics, takes values from #eelCu */
    int             vdwtype;          /**< type of VdW impl., takes values from #evdwCu     */

    float           epsfac;           /**< charge multiplication factor                      */
    float           c_rf;             /**< Reaction-field/plain cutoff electrostatics const. */
    float           two_k_rf;         /**< Reaction-field electrostatics constant            */
    float           ewald_beta;       /**< Ewald/PME parameter                               */
    float           sh_ewald;         /**< Ewald/PME correction term substracted from the direct-space potential */
    float           sh_lj_ewald;      /**< LJ-Ewald/PME correction term added to the correction potential        */
    float           ewaldcoeff_lj;    /**< LJ-Ewald/PME coefficient                          */

    float           rcoulomb_sq;      /**< Coulomb cut-off squared                           */

    float           rvdw_sq;          /**< VdW cut-off squared                               */
    float           rvdw_switch;      /**< VdW switched cut-off                              */
    float           rlist_sq;         /**< pair-list cut-off squared                         */

    shift_consts_t  dispersion_shift; /**< VdW shift dispersion constants           */
    shift_consts_t  repulsion_shift;  /**< VdW shift repulsion constants            */
    switch_consts_t vdw_switch;       /**< VdW switch constants                     */

    /* LJ non-bonded parameters - accessed through texture memory */
    float               *nbfp;             /**< nonbonded parameter table with C6/C12 pairs per atom type-pair, 2*ntype^2 elements */
    cudaTextureObject_t  nbfp_texobj;      /**< texture object bound to nbfp                                                       */
    float               *nbfp_comb;        /**< nonbonded parameter table per atom type, 2*ntype elements                          */
    cudaTextureObject_t  nbfp_comb_texobj; /**< texture object bound to nbfp_texobj                                                */

    /* Ewald Coulomb force table data - accessed through texture memory */
    int                  coulomb_tab_size;   /**< table size (s.t. it fits in texture cache) */
    float                coulomb_tab_scale;  /**< table scale/spacing                        */
    float               *coulomb_tab;        /**< pointer to the table in the device memory  */
    cudaTextureObject_t  coulomb_tab_texobj; /**< texture object bound to coulomb_tab        */
};

/** \internal
 * \brief Pair list data.
 */
struct cu_plist
{
    int              na_c;        /**< number of atoms per cluster                  */

    int              nsci;        /**< size of sci, # of i clusters in the list     */
    int              sci_nalloc;  /**< allocation size of sci                       */
    nbnxn_sci_t     *sci;         /**< list of i-cluster ("super-clusters")         */

    int              ncj4;        /**< total # of 4*j clusters                      */
    int              cj4_nalloc;  /**< allocation size of cj4                       */
    nbnxn_cj4_t     *cj4;         /**< 4*j cluster list, contains j cluster number
                                       and index into the i cluster list            */
    nbnxn_excl_t    *excl;        /**< atom interaction bits                        */
    int              nexcl;       /**< count for excl                               */
    int              excl_nalloc; /**< allocation size of excl                      */

    bool             bDoPrune;    /**< true if pair-list pruning needs to be
                                       done during the  current step                */
};

/** \internal
 * \brief CUDA events used for timing GPU kernels and H2D/D2H transfers.
 *
 * The two-sized arrays hold the local and non-local values and should always
 * be indexed with eintLocal/eintNonlocal.
 */
struct cu_timers
{
    cudaEvent_t start_atdat;     /**< start event for atom data transfer (every PS step)             */
    cudaEvent_t stop_atdat;      /**< stop event for atom data transfer (every PS step)              */
    cudaEvent_t start_nb_h2d[2]; /**< start events for x/q H2D transfers (l/nl, every step)          */
    cudaEvent_t stop_nb_h2d[2];  /**< stop events for x/q H2D transfers (l/nl, every step)           */
    cudaEvent_t start_nb_d2h[2]; /**< start events for f D2H transfer (l/nl, every step)             */
    cudaEvent_t stop_nb_d2h[2];  /**< stop events for f D2H transfer (l/nl, every step)              */
    cudaEvent_t start_pl_h2d[2]; /**< start events for pair-list H2D transfers (l/nl, every PS step) */
    cudaEvent_t stop_pl_h2d[2];  /**< start events for pair-list H2D transfers (l/nl, every PS step) */
    cudaEvent_t start_nb_k[2];   /**< start event for non-bonded kernels (l/nl, every step)          */
    cudaEvent_t stop_nb_k[2];    /**< stop event non-bonded kernels (l/nl, every step)               */
};

/** \internal
 * \brief Main data structure for CUDA nonbonded force calculations.
 */
struct gmx_nbnxn_cuda_t
{
    struct gmx_device_info_t *dev_info;       /**< CUDA device information                              */
    bool                      bUseTwoStreams; /**< true if doing both local/non-local NB work on GPU    */
    bool                      bUseStreamSync; /**< true if the standard cudaStreamSynchronize is used
                                                   and not memory polling-based waiting                 */
    cu_atomdata_t            *atdat;          /**< atom data                                            */
    cu_nbparam_t             *nbparam;        /**< parameters required for the non-bonded calc.         */
    cu_plist_t               *plist[2];       /**< pair-list data structures (local and non-local)      */
    nb_staging_t              nbst;           /**< staging area where fshift/energies get downloaded    */

    cudaStream_t              stream[2];      /**< local and non-local GPU streams                      */

    /** events used for synchronization */
    cudaEvent_t    nonlocal_done;               /**< event triggered when the non-local non-bonded kernel
                                                   is done (and the local transfer can proceed)           */
    cudaEvent_t    misc_ops_and_local_H2D_done; /**< event triggered when the tasks issued in
                                                   the local stream that need to precede the
                                                   non-local force calculations are done
                                                   (e.g. f buffer 0-ing, local x/q H2D) */

    /* NOTE: With current CUDA versions (<=5.0) timing doesn't work with multiple
     * concurrent streams, so we won't time if both l/nl work is done on GPUs.
     * Timer init/uninit is still done even with timing off so only the condition
     * setting bDoTime needs to be change if this CUDA "feature" gets fixed. */
    bool                 bDoTime;   /**< True if event-based timing is enabled.               */
    cu_timers_t         *timers;    /**< CUDA event-based timers.                             */
    gmx_wallclock_gpu_t *timings;   /**< Timing data.                                         */
};

#ifdef __cplusplus
}
#endif

#endif  /* NBNXN_CUDA_TYPES_H */
