/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef NBNXN_CUDA_TYPES_H
#define NBNXN_CUDA_TYPES_H

#include "types/nbnxn_pairlist.h"
#include "types/nbnxn_cuda_types_ext.h"
#include "../../gmxlib/cuda_tools/cudautils.cuh"

#ifdef __cplusplus
extern "C" {
#endif

/*! Types of electrostatics implementations available in the CUDA non-bonded
 *  force kernels. These represent both the electrostatics types implemented
 *  by the kernels (cut-off, RF, and Ewald - a subset of what's defined in
 *  enums.h) as well as encode implementation details analytical/tabulated
 *  and single or twin cut-off (for Ewald kernels).
 *  Note that the cut-off and RF kernels have only analytical flavor and unlike
 *  in the CPU kernels, the tabulated kernels are ATM Ewald-only.
 *
 *  The order of pointers to different electrostatic kernels defined in
 *  nbnxn_cuda.cu by the nb_default_kfunc_ptr and nb_legacy_kfunc_ptr arrays
 *  should match the order of enumerated types below. */
enum {
    eelCuCUT, eelCuRF, eelCuEWALD_TAB, eelCuEWALD_TAB_TWIN, eelCuEWALD_ANA, eelCuEWALD_ANA_TWIN, eelCuNR
};

/*! Kernel flavors with different set of optimizations: default for CUDA <=v4.1
 *  compilers and legacy for earlier, 3.2 and 4.0 CUDA compilers. */
enum {
    eNbnxnCuKDefault, eNbnxnCuKLegacy, eNbnxnCuKNR
};

#define NBNXN_KVER_OLD(k)      (k == eNbnxnCuKOld)
#define NBNXN_KVER_LEGACY(k)   (k == eNbnxnCuKLegacy)
#define NBNXN_KVER_DEFAULT(k)  (k == eNbnxnCuKDefault)

/*! Non-bonded kernel versions. */

/* All structs prefixed with "cu_" hold data used in GPU calculations and
 * are passed to the kernels, except cu_timers_t. */
typedef struct cu_plist     cu_plist_t;
typedef struct cu_atomdata  cu_atomdata_t;
typedef struct cu_nbparam   cu_nbparam_t;
typedef struct cu_timers    cu_timers_t;
typedef struct nb_staging   nb_staging_t;


/*! Staging area for temporary data. The energies get downloaded here first,
 *  before getting added to the CPU-side aggregate values.
 */
struct nb_staging
{
    float   *e_lj;      /* LJ energy            */
    float   *e_el;      /* electrostatic energy */
    float3  *fshift;    /* shift forces         */
};

/*! Nonbonded atom data -- both inputs and outputs. */
struct cu_atomdata
{
    int      natoms;            /* number of atoms                              */
    int      natoms_local;      /* number of local atoms                        */
    int      nalloc;            /* allocation size for the atom data (xq, f)    */

    float4  *xq;                /* atom coordinates + charges, size natoms      */
    float3  *f;                 /* force output array, size natoms              */
    /* TODO: try float2 for the energies */
    float   *e_lj,              /* LJ energy output, size 1                     */
            *e_el;              /* Electrostatics energy input, size 1          */

    float3  *fshift;            /* shift forces                                 */

    int      ntypes;            /* number of atom types                         */
    int     *atom_types;        /* atom type indices, size natoms               */

    float3  *shift_vec;         /* shifts                                       */
    bool     bShiftVecUploaded; /* true if the shift vector has been uploaded   */
};

/*! Parameters required for the CUDA nonbonded calculations. */
struct cu_nbparam
{
    int      eeltype;        /* type of electrostatics                       */

    float    epsfac;         /* charge multiplication factor                 */
    float    c_rf, two_k_rf; /* Reaction-Field constants                     */
    float    ewald_beta;     /* Ewald/PME parameter                          */
    float    sh_ewald;       /* Ewald/PME  correction term                   */
    float    rvdw_sq;        /* VdW cut-off                                  */
    float    rcoulomb_sq;    /* Coulomb cut-off                              */
    float    rlist_sq;       /* pair-list cut-off                            */
    float    sh_invrc6;      /* LJ potential correction term                 */

    float   *nbfp;           /* nonbonded parameter table with C6/C12 pairs  */

    /* Ewald Coulomb force table */
    int      coulomb_tab_size;
    float    coulomb_tab_scale;
    float   *coulomb_tab;
};

/*! Pair list data */
struct cu_plist
{
    int              na_c;        /* number of atoms per cluster                  */

    int              nsci;        /* size of sci, # of i clusters in the list     */
    int              sci_nalloc;  /* allocation size of sci                       */
    nbnxn_sci_t     *sci;         /* list of i-cluster ("super-clusters")         */

    int              ncj4;        /* total # of 4*j clusters                      */
    int              cj4_nalloc;  /* allocation size of cj4                       */
    nbnxn_cj4_t     *cj4;         /* 4*j cluster list, contains j cluster number
                                     and index into the i cluster list            */
    nbnxn_excl_t    *excl;        /* atom interaction bits                        */
    int              nexcl;       /* count for excl                               */
    int              excl_nalloc; /* allocation size of excl                      */

    bool             bDoPrune;    /* true if pair-list pruning needs to be
                                     done during the  current step                */
};

/* CUDA events used for timing GPU kernels and H2D/D2H transfers.
 * The two-sized arrays hold the local and non-local values and should always
 * be indexed with eintLocal/eintNonlocal.
 */
struct cu_timers
{
    cudaEvent_t start_atdat, stop_atdat;         /* atom data transfer (every PS step)      */
    cudaEvent_t start_nb_h2d[2], stop_nb_h2d[2]; /* x/q H2D transfer (every step)           */
    cudaEvent_t start_nb_d2h[2], stop_nb_d2h[2]; /* f D2H transfer (every step)             */
    cudaEvent_t start_pl_h2d[2], stop_pl_h2d[2]; /* pair-list H2D transfer (every PS step)  */
    cudaEvent_t start_nb_k[2], stop_nb_k[2];     /* non-bonded kernels (every step)         */
};

/* Main data structure for CUDA nonbonded force calculations. */
struct nbnxn_cuda
{
    cuda_dev_info_t *dev_info;       /* CUDA device information                              */
    int              kernel_ver;     /* The version of the kernel to be executed on the
                                        device in use, possible values: eNbnxnCuK*           */
    bool             bUseTwoStreams; /* true if doing both local/non-local NB work on GPU    */
    bool             bUseStreamSync; /* true if the standard cudaStreamSynchronize is used
                                        and not memory polling-based waiting                 */
    cu_atomdata_t   *atdat;          /* atom data                                            */
    cu_nbparam_t    *nbparam;        /* parameters required for the non-bonded calc.         */
    cu_plist_t      *plist[2];       /* pair-list data structures (local and non-local)      */
    nb_staging_t     nbst;           /* staging area where fshift/energies get downloaded    */

    cudaStream_t     stream[2];      /* local and non-local GPU streams                      */

    /* events used for synchronization */
    cudaEvent_t    nonlocal_done, misc_ops_done;

    /* NOTE: With current CUDA versions (<=5.0) timing doesn't work with multiple
     * concurrent streams, so we won't time if both l/nl work is done on GPUs.
     * Timer init/uninit is still done even with timing off so only the condition
     * setting bDoTime needs to be change if this CUDA "feature" gets fixed. */
    bool             bDoTime;       /* True if event-based timing is enabled.               */
    cu_timers_t     *timers;        /* CUDA event-based timers.                             */
    wallclock_gpu_t *timings;       /* Timing data.                                         */
};

#ifdef __cplusplus
}
#endif

#endif  /* NBNXN_CUDA_TYPES_H */
