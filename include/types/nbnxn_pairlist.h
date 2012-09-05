/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2012, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 *
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */

#ifndef _nbnxn_pairlist_h
#define _nbnxn_pairlist_h

#ifdef __cplusplus
extern "C" {
#endif

/* A buffer data structure of 64 bytes
 * to be placed at the beginning and end of structs
 * to avoid cache invalidation of the real contents
 * of the struct by writes to neighboring memory.
 */
typedef struct {
    int dummy[16];
} gmx_cache_protect_t;

/* Abstract type for pair searching data */
typedef struct nbnxn_search * nbnxn_search_t;

/* Function that should return a pointer *ptr to memory
 * of size nbytes.
 * Error handling should be done within this function.
 */
typedef void gmx_nbat_alloc_t(void **ptr,size_t nbytes);

/* Function that should free the memory pointed to by *ptr.
 * NULL should not be passed to this function.
 */
typedef void gmx_nbat_free_t(void *ptr);

typedef struct {
    int      cj;    /* The j-cluster                    */
    unsigned excl;  /* The exclusion (interaction) bits */
} nbnxn_cj_t;

#define NBNXN_CI_SHIFT          127
#define NBNXN_CI_DO_LJ(subc)    (1<<(7+3*(subc)))
#define NBNXN_CI_HALF_LJ(subc)  (1<<(8+3*(subc)))
#define NBNXN_CI_DO_COUL(subc)  (1<<(9+3*(subc)))

/* Simple pair-list i-unit */
typedef struct {
    int ci;             /* i-cluster             */
    int shift;          /* Shift vector index plus possible flags */
    int cj_ind_start;   /* Start index into cj   */
    int cj_ind_end;     /* End index into cj     */
} nbnxn_ci_t;

/* Grouped pair-list i-unit */
typedef struct {
    int sci;            /* i-super-cluster       */
    int shift;          /* Shift vector index plus possible flags */
    int cj4_ind_start;  /* Start index into cj4  */
    int cj4_ind_end;    /* End index into cj4    */
} nbnxn_sci_t;

typedef struct {
    unsigned imask;        /* The i-cluster interactions mask for 1 warp  */
    int excl_ind;          /* Index into the exclusion array for 1 warp   */
} nbnxn_im_ei_t;

typedef struct {
    int cj[4];             /* The 4 j-clusters                            */
    nbnxn_im_ei_t imei[2]; /* The i-cluster mask data       for 2 warps   */
} nbnxn_cj4_t;

typedef struct {
    unsigned pair[32];     /* Exclusion bits for one warp,                *
                            * each unsigned has bit for 4*8 i clusters    */
} nbnxn_excl_t;

typedef struct {
    gmx_cache_protect_t cp0;

    gmx_nbat_alloc_t *alloc;
    gmx_nbat_free_t  *free;

    gmx_bool bSimple;      /* Simple list has na_sc=na_s and uses cj   *
                            * Complex list uses cj4                    */

    int      na_ci;        /* The number of atoms per i-cluster        */
    int      na_cj;        /* The number of atoms per j-cluster        */
    int      na_sc;        /* The number of atoms per super cluster    */
    real     rlist;        /* The radius for constructing the list     */
    int      nci;          /* The number of i-clusters in the list     */
    nbnxn_ci_t *ci;        /* The i-cluster list, size nci             */
    int      ci_nalloc;    /* The allocation size of ci                */
    int      nsci;         /* The number of i-super-clusters in the list */
    nbnxn_sci_t *sci;      /* The i-super-cluster list                 */
    int      sci_nalloc;   /* The allocation size of sci               */

    int      ncj;          /* The number of j-clusters in the list     */
    nbnxn_cj_t *cj;        /* The j-cluster list, size ncj             */
    int      cj_nalloc;    /* The allocation size of cj                */

    int      ncj4;         /* The total number of 4*j clusters         */
    nbnxn_cj4_t *cj4;      /* The 4*j cluster list, size ncj4          */
    int      cj4_nalloc;   /* The allocation size of cj4               */
    int      nexcl;        /* The count for excl                       */
    nbnxn_excl_t *excl;    /* Atom interaction bits (non-exclusions)   */
    int      excl_nalloc;  /* The allocation size for excl             */
    int      nci_tot;      /* The total number of i clusters           */

    struct nbnxn_list_work *work;

    gmx_cache_protect_t cp1;
} nbnxn_pairlist_t;

typedef struct {
    int          nnbl;      /* number of lists */
    nbnxn_pairlist_t **nbl; /* lists */
    gmx_bool     bCombined; /* TRUE if lists get combined into one (the 1st) */
    gmx_bool     bSimple;   /* TRUE if the list of of type "simple"
                               (na_sc=na_s, no super-clusters used) */
    int          natpair_ljq; /* Total number of atom pairs for LJ+Q kernel */
    int          natpair_lj;  /* Total number of atom pairs for LJ kernel   */
    int          natpair_q;   /* Total number of atom pairs for Q kernel    */
} nbnxn_pairlist_set_t;

enum { nbatXYZ, nbatXYZQ, nbatX4, nbatX8 };

typedef struct {
    real *f;      /* f, size natoms*fstride                             */
    real *fshift; /* Shift force array, size SHIFTS*DIM                 */
    int  nV;      /* The size of *Vvdw and *Vc                          */
    real *Vvdw;   /* Temporary Van der Waals group energy storage       */
    real *Vc;     /* Temporary Coulomb group energy storage             */
    int  nVS;     /* The size of *VSvdw and *VSc                        */
    real *VSvdw;  /* Temporary SIMD Van der Waals group energy storage  */
    real *VSc;    /* Temporary SIMD Coulomb group energy storage        */
} nbnxn_atomdata_output_t;

/* LJ combination rules: geometric, Lorentz-Berthelot, none */
enum { ljcrGEOM, ljcrLB, ljcrNONE, ljcrNR };

typedef struct {
    gmx_nbat_alloc_t *alloc;
    gmx_nbat_free_t  *free;
    int  ntype;      /* The number of different atom types                 */
    real *nbfp;      /* Lennard-Jones 6*C6 and 12*C12 params, size ntype^2*2 */
    int  comb_rule;  /* Combination rule, see enum above                   */
    real *nbfp_comb; /* LJ parameter per atom type, size ntype*2           */
    real *nbfp_s4;   /* As nbfp, but with stride 4, size ntype^2*4         */
    int  natoms;     /* Number of atoms                                    */
    int  natoms_local;  /* Number of local atoms                           */
    int  *type;      /* Atom types                                         */
    real *lj_comb;   /* LJ parameters per atom for combining for pairs     */
    int  XFormat;    /* The format of x (and q), enum                      */
    int  FFormat;    /* The format of f, enum                              */
    real *q;         /* Charges, can be NULL if incorporated in x          */
    int  na_c;       /* The number of atoms per cluster                    */
    int  nenergrp;   /* The number of energy groups                        */
    int  neg_2log;   /* Log2 of nenergrp                                   */
    int  *energrp;   /* The energy groups per cluster, can be NULL         */
    gmx_bool bDynamicBox; /* Do we need to update shift_vec every step?    */
    rvec *shift_vec; /* Shift vectors, copied from t_forcerec              */
    int  xstride;    /* stride for a coordinate in x (usually 3 or 4)      */
    int  fstride;    /* stride for a coordinate in f (usually 3 or 4)      */
    real *x;         /* x and possibly q, size natoms*xstride              */
    int  nout;       /* The number of force arrays                         */
    nbnxn_atomdata_output_t *out;  /* Output data structures               */
    int  nalloc;     /* Allocation size of all arrays (for x/f *x/fstride) */
} nbnxn_atomdata_t;

#ifdef __cplusplus
}
#endif

#endif
