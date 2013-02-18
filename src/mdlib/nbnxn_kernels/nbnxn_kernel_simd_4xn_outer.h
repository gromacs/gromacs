/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
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

#if !(GMX_NBNXN_SIMD_BITWIDTH == 128 || GMX_NBNXN_SIMD_BITWIDTH == 256)
#error "unsupported GMX_NBNXN_SIMD_BITWIDTH"
#endif

#ifdef GMX_NBNXN_HALF_WIDTH_SIMD
#define GMX_USE_HALF_WIDTH_SIMD_HERE
#endif
#include "gmx_simd_macros.h"

#define SUM_SIMD4(x) (x[0]+x[1]+x[2]+x[3])

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    GMX_SIMD_WIDTH_HERE

/* The stride of all the atom data arrays is max(UNROLLI,UNROLLJ) */
#if GMX_SIMD_WIDTH_HERE >= UNROLLI
#define STRIDE     GMX_SIMD_WIDTH_HERE
#else
#define STRIDE     UNROLLI
#endif

#if GMX_SIMD_WIDTH_HERE == 2
#define SUM_SIMD(x)  (x[0]+x[1])
#else
#if GMX_SIMD_WIDTH_HERE == 4
#define SUM_SIMD(x)  SUM_SIMD4(x)
#else
#if GMX_SIMD_WIDTH_HERE == 8
#define SUM_SIMD(x)  (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])
#else
#error "unsupported kernel configuration"
#endif
#endif
#endif


/* Decide if we should use the FDV0 table layout */
#if defined GMX_X86_AVX_256 && !defined GMX_USE_HALF_WIDTH_SIMD_HERE
/* With full AVX-256 SIMD, half SIMD-width table loads are optimal */
#if GMX_SIMD_WIDTH_HERE/2 == 4
#define TAB_FDV0
#endif
#else
/* We use the FDV0 table layout when we can use aligned table loads */
#if GMX_SIMD_WIDTH_HERE == 4
#define TAB_FDV0
#endif
#endif


#define SIMD_MASK_ALL   0xffffffff

#include "nbnxn_kernel_simd_utils.h"

/* All functionality defines are set here, except for:
 * CALC_ENERGIES, ENERGY_GROUPS which are defined before.
 * CHECK_EXCLS, which is set just before including the inner loop contents.
 * The combination rule defines, LJ_COMB_GEOM or LJ_COMB_LB are currently
 * set before calling the kernel function. We might want to move that
 * to inside the n-loop and have a different combination rule for different
 * ci's, as no combination rule gives a 50% performance hit for LJ.
 */

/* We always calculate shift forces, because it's cheap anyhow */
#define CALC_SHIFTFORCES

/* Assumes all LJ parameters are identical */
/* #define FIX_LJ_C */

/* The NBK_FUNC_NAME... macros below generate the whole zoo of kernels names
 * with all combinations off electrostatics (coul), LJ combination rules (ljc)
 * and energy calculations (ene), depending on the defines set.
 */

#define NBK_FUNC_NAME_C_LJC(base, coul, ljc, ene) base ## _ ## coul ## _comb_ ## ljc ## _ ## ene

#if defined LJ_COMB_GEOM
#define NBK_FUNC_NAME_C(base, coul, ene) NBK_FUNC_NAME_C_LJC(base, coul, geom, ene)
#else
#if defined LJ_COMB_LB
#define NBK_FUNC_NAME_C(base, coul, ene) NBK_FUNC_NAME_C_LJC(base, coul, lb, ene)
#else
#define NBK_FUNC_NAME_C(base, coul, ene) NBK_FUNC_NAME_C_LJC(base, coul, none, ene)
#endif
#endif

#ifdef CALC_COUL_RF
#define NBK_FUNC_NAME(base, ene) NBK_FUNC_NAME_C(base, rf, ene)
#endif
#ifdef CALC_COUL_TAB
#ifndef VDW_CUTOFF_CHECK
#define NBK_FUNC_NAME(base, ene) NBK_FUNC_NAME_C(base, tab, ene)
#else
#define NBK_FUNC_NAME(base, ene) NBK_FUNC_NAME_C(base, tab_twin, ene)
#endif
#endif
#ifdef CALC_COUL_EWALD
#ifndef VDW_CUTOFF_CHECK
#define NBK_FUNC_NAME(base, ene) NBK_FUNC_NAME_C(base, ewald, ene)
#else
#define NBK_FUNC_NAME(base, ene) NBK_FUNC_NAME_C(base, ewald_twin, ene)
#endif
#endif

static void
#ifndef CALC_ENERGIES
NBK_FUNC_NAME(nbnxn_kernel_simd_4xn, noener)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME(nbnxn_kernel_simd_4xn, ener)
#else
NBK_FUNC_NAME(nbnxn_kernel_simd_4xn, energrp)
#endif
#endif
#undef NBK_FUNC_NAME
#undef NBK_FUNC_NAME_C
#undef NBK_FUNC_NAME_C_LJC
(const nbnxn_pairlist_t     *nbl,
 const nbnxn_atomdata_t     *nbat,
 const interaction_const_t  *ic,
 rvec                       *shift_vec,
 real                       *f
#ifdef CALC_SHIFTFORCES
 ,
 real                       *fshift
#endif
#ifdef CALC_ENERGIES
 ,
 real                       *Vvdw,
 real                       *Vc
#endif
)
{
    const nbnxn_ci_t   *nbln;
    const nbnxn_cj_t   *l_cj;
    const int          *type;
    const real         *q;
    const real         *shiftvec;
    const real         *x;
    const real         *nbfp0, *nbfp1, *nbfp2 = NULL, *nbfp3 = NULL;
    real                facel;
    real               *nbfp_ptr;
    int                 nbfp_stride;
    int                 n, ci, ci_sh;
    int                 ish, ish3;
    gmx_bool            do_LJ, half_LJ, do_coul;
    int                 sci, scix, sciy, sciz, sci2;
    int                 cjind0, cjind1, cjind;
    int                 ip, jp;

#ifdef ENERGY_GROUPS
    int         Vstride_i;
    int         egps_ishift, egps_imask;
    int         egps_jshift, egps_jmask, egps_jstride;
    int         egps_i;
    real       *vvdwtp[UNROLLI];
    real       *vctp[UNROLLI];
#endif

    gmx_mm_pr  shX_S;
    gmx_mm_pr  shY_S;
    gmx_mm_pr  shZ_S;
    gmx_mm_pr  ix_S0, iy_S0, iz_S0;
    gmx_mm_pr  ix_S1, iy_S1, iz_S1;
    gmx_mm_pr  ix_S2, iy_S2, iz_S2;
    gmx_mm_pr  ix_S3, iy_S3, iz_S3;
    gmx_mm_pr  fix_S0, fiy_S0, fiz_S0;
    gmx_mm_pr  fix_S1, fiy_S1, fiz_S1;
    gmx_mm_pr  fix_S2, fiy_S2, fiz_S2;
    gmx_mm_pr  fix_S3, fiy_S3, fiz_S3;
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
    __m128     fix_S, fiy_S, fiz_S;
#else
    __m256d    fix_S, fiy_S, fiz_S;
#endif
#else
    __m128d    fix0_S, fiy0_S, fiz0_S;
    __m128d    fix2_S, fiy2_S, fiz2_S;
#endif

    gmx_mm_pr diag_jmi_S;
#if UNROLLI == UNROLLJ
    gmx_mm_pr diag_S0, diag_S1, diag_S2, diag_S3;
#else
    gmx_mm_pr diag0_S0, diag0_S1, diag0_S2, diag0_S3;
    gmx_mm_pr diag1_S0, diag1_S1, diag1_S2, diag1_S3;
#endif

#ifdef gmx_checkbitmask_epi32
    gmx_epi32 mask_S0, mask_S1, mask_S2, mask_S3;
#else
    gmx_mm_pr mask_S0, mask_S1, mask_S2, mask_S3;
#endif

    gmx_mm_pr  zero_S = gmx_set1_pr(0);

    gmx_mm_pr  one_S  = gmx_set1_pr(1.0);
    gmx_mm_pr  iq_S0  = gmx_setzero_pr();
    gmx_mm_pr  iq_S1  = gmx_setzero_pr();
    gmx_mm_pr  iq_S2  = gmx_setzero_pr();
    gmx_mm_pr  iq_S3  = gmx_setzero_pr();
    gmx_mm_pr  mrc_3_S;
#ifdef CALC_ENERGIES
    gmx_mm_pr  hrc_3_S, moh_rc_S;
#endif

#ifdef CALC_COUL_TAB
    /* Coulomb table variables */
    gmx_mm_pr   invtsp_S;
    const real *tab_coul_F;
#ifndef TAB_FDV0
    const real *tab_coul_V;
#endif
#if GMX_NBNXN_SIMD_BITWIDTH == 256
    int        ti0_array[2*GMX_SIMD_WIDTH_HERE-1], *ti0;
    int        ti1_array[2*GMX_SIMD_WIDTH_HERE-1], *ti1;
    int        ti2_array[2*GMX_SIMD_WIDTH_HERE-1], *ti2;
    int        ti3_array[2*GMX_SIMD_WIDTH_HERE-1], *ti3;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  mhalfsp_S;
#endif
#endif

#ifdef CALC_COUL_EWALD
    gmx_mm_pr beta2_S, beta_S;
#endif

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    gmx_mm_pr  sh_ewald_S;
#endif

#ifdef LJ_COMB_LB
    const real *ljc;

    gmx_mm_pr   hsig_i_S0, seps_i_S0;
    gmx_mm_pr   hsig_i_S1, seps_i_S1;
    gmx_mm_pr   hsig_i_S2, seps_i_S2;
    gmx_mm_pr   hsig_i_S3, seps_i_S3;
#else
#ifdef FIX_LJ_C
    real        pvdw_array[2*UNROLLI*UNROLLJ+3];
    real       *pvdw_c6, *pvdw_c12;
    gmx_mm_pr   c6_S0, c12_S0;
    gmx_mm_pr   c6_S1, c12_S1;
    gmx_mm_pr   c6_S2, c12_S2;
    gmx_mm_pr   c6_S3, c12_S3;
#endif

#ifdef LJ_COMB_GEOM
    const real *ljc;

    gmx_mm_pr   c6s_S0, c12s_S0;
    gmx_mm_pr   c6s_S1, c12s_S1;
    gmx_mm_pr   c6s_S2 = gmx_setzero_pr(), c12s_S2 = gmx_setzero_pr();
    gmx_mm_pr   c6s_S3 = gmx_setzero_pr(), c12s_S3 = gmx_setzero_pr();
#endif
#endif /* LJ_COMB_LB */

    gmx_mm_pr  vctot_S, Vvdwtot_S;
    gmx_mm_pr  sixth_S, twelveth_S;

    gmx_mm_pr  avoid_sing_S;
    gmx_mm_pr  rc2_S;
#ifdef VDW_CUTOFF_CHECK
    gmx_mm_pr  rcvdw2_S;
#endif

#ifdef CALC_ENERGIES
    gmx_mm_pr  sh_invrc6_S, sh_invrc12_S;

    /* cppcheck-suppress unassignedVariable */
    real       tmpsum_array[15], *tmpsum;
#endif
#ifdef CALC_SHIFTFORCES
    /* cppcheck-suppress unassignedVariable */
    real       shf_array[15], *shf;
#endif

    int ninner;

#ifdef COUNT_PAIRS
    int npair = 0;
#endif

#if defined LJ_COMB_GEOM || defined LJ_COMB_LB
    ljc = nbat->lj_comb;
#else
    /* No combination rule used */
#ifndef GMX_DOUBLE
    nbfp_ptr    = nbat->nbfp_s4;
#define NBFP_STRIDE  4
#else
    nbfp_ptr    = nbat->nbfp;
#define NBFP_STRIDE  2
#endif
    nbfp_stride = NBFP_STRIDE;
#endif

    /* Load j-i for the first i */
    diag_jmi_S = gmx_load_pr(nbat->simd_4xn_diag);
    /* Generate all the diagonal masks as comparison results */
#if UNROLLI == UNROLLJ
    diag_S0    = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag_S1    = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag_S2    = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag_S3    = gmx_cmplt_pr(zero_S, diag_jmi_S);
#else
#if UNROLLI == 2*UNROLLJ || 2*UNROLLI == UNROLLJ
    diag0_S0   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag0_S1   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag0_S2   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag0_S3   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);

#if UNROLLI == 2*UNROLLJ
    /* Load j-i for the second half of the j-cluster */
    diag_jmi_S = gmx_load_pr(nbat->simd_4xn_diag+UNROLLJ);
#endif

    diag1_S0   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag1_S1   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag1_S2   = gmx_cmplt_pr(zero_S, diag_jmi_S);
    diag_jmi_S = gmx_sub_pr(diag_jmi_S, one_S);
    diag1_S3   = gmx_cmplt_pr(zero_S, diag_jmi_S);
#endif
#endif

    /* Load masks for topology exclusion masking */
#ifdef gmx_checkbitmask_epi32
    mask_S0    = gmx_load_si(nbat->simd_excl_mask + 0*GMX_NBNXN_SIMD_BITWIDTH/32);
    mask_S1    = gmx_load_si(nbat->simd_excl_mask + 1*GMX_NBNXN_SIMD_BITWIDTH/32);
    mask_S2    = gmx_load_si(nbat->simd_excl_mask + 2*GMX_NBNXN_SIMD_BITWIDTH/32);
    mask_S3    = gmx_load_si(nbat->simd_excl_mask + 3*GMX_NBNXN_SIMD_BITWIDTH/32);
#else
    mask_S0    = gmx_load_pr((real *)nbat->simd_excl_mask + 0*UNROLLJ);
    mask_S1    = gmx_load_pr((real *)nbat->simd_excl_mask + 1*UNROLLJ);
    mask_S2    = gmx_load_pr((real *)nbat->simd_excl_mask + 2*UNROLLJ);
    mask_S3    = gmx_load_pr((real *)nbat->simd_excl_mask + 3*UNROLLJ);
#endif

#ifdef CALC_COUL_TAB
#if GMX_NBNXN_SIMD_BITWIDTH == 256
    /* Generate aligned table index pointers */
    ti0 = gmx_simd_align_int(ti0_array);
    ti1 = gmx_simd_align_int(ti1_array);
    ti2 = gmx_simd_align_int(ti2_array);
    ti3 = gmx_simd_align_int(ti3_array);
#endif

    invtsp_S  = gmx_set1_pr(ic->tabq_scale);
#ifdef CALC_ENERGIES
    mhalfsp_S = gmx_set1_pr(-0.5/ic->tabq_scale);
#endif

#ifdef TAB_FDV0
    tab_coul_F = ic->tabq_coul_FDV0;
#else
    tab_coul_F = ic->tabq_coul_F;
    tab_coul_V = ic->tabq_coul_V;
#endif
#endif /* CALC_COUL_TAB */

#ifdef CALC_COUL_EWALD
    beta2_S = gmx_set1_pr(ic->ewaldcoeff*ic->ewaldcoeff);
    beta_S  = gmx_set1_pr(ic->ewaldcoeff);
#endif

#if (defined CALC_COUL_TAB || defined CALC_COUL_EWALD) && defined CALC_ENERGIES
    sh_ewald_S = gmx_set1_pr(ic->sh_ewald);
#endif

    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

    avoid_sing_S = gmx_set1_pr(NBNXN_AVOID_SING_R2_INC);

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    rc2_S    = gmx_set1_pr(ic->rcoulomb*ic->rcoulomb);
#ifdef VDW_CUTOFF_CHECK
    rcvdw2_S = gmx_set1_pr(ic->rvdw*ic->rvdw);
#endif

#ifdef CALC_ENERGIES
    sixth_S      = gmx_set1_pr(1.0/6.0);
    twelveth_S   = gmx_set1_pr(1.0/12.0);

    sh_invrc6_S  = gmx_set1_pr(ic->sh_invrc6);
    sh_invrc12_S = gmx_set1_pr(ic->sh_invrc6*ic->sh_invrc6);
#endif

    mrc_3_S  = gmx_set1_pr(-2*ic->k_rf);

#ifdef CALC_ENERGIES
    hrc_3_S  = gmx_set1_pr(ic->k_rf);

    moh_rc_S = gmx_set1_pr(-ic->c_rf);
#endif

#ifdef CALC_ENERGIES
    tmpsum   = gmx_simd_align_real(tmpsum_array);
#endif
#ifdef CALC_SHIFTFORCES
    shf      = gmx_simd_align_real(shf_array);
#endif

#ifdef FIX_LJ_C
    pvdw_c6  = gmx_simd_align_real(pvdw_array+3);
    pvdw_c12 = pvdw_c6 + UNROLLI*UNROLLJ;

    for (jp = 0; jp < UNROLLJ; jp++)
    {
        pvdw_c6 [0*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [1*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [2*UNROLLJ+jp] = nbat->nbfp[0*2];
        pvdw_c6 [3*UNROLLJ+jp] = nbat->nbfp[0*2];

        pvdw_c12[0*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[1*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[2*UNROLLJ+jp] = nbat->nbfp[0*2+1];
        pvdw_c12[3*UNROLLJ+jp] = nbat->nbfp[0*2+1];
    }
    c6_S0            = gmx_load_pr(pvdw_c6 +0*UNROLLJ);
    c6_S1            = gmx_load_pr(pvdw_c6 +1*UNROLLJ);
    c6_S2            = gmx_load_pr(pvdw_c6 +2*UNROLLJ);
    c6_S3            = gmx_load_pr(pvdw_c6 +3*UNROLLJ);

    c12_S0           = gmx_load_pr(pvdw_c12+0*UNROLLJ);
    c12_S1           = gmx_load_pr(pvdw_c12+1*UNROLLJ);
    c12_S2           = gmx_load_pr(pvdw_c12+2*UNROLLJ);
    c12_S3           = gmx_load_pr(pvdw_c12+3*UNROLLJ);
#endif /* FIX_LJ_C */

#ifdef ENERGY_GROUPS
    egps_ishift  = nbat->neg_2log;
    egps_imask   = (1<<egps_ishift) - 1;
    egps_jshift  = 2*nbat->neg_2log;
    egps_jmask   = (1<<egps_jshift) - 1;
    egps_jstride = (UNROLLJ>>1)*UNROLLJ;
    /* Major division is over i-particle energy groups, determine the stride */
    Vstride_i    = nbat->nenergrp*(1<<nbat->neg_2log)*egps_jstride;
#endif

    l_cj = nbl->cj;

    ninner = 0;
    for (n = 0; n < nbl->nci; n++)
    {
        nbln = &nbl->ci[n];

        ish              = (nbln->shift & NBNXN_CI_SHIFT);
        ish3             = ish*3;
        cjind0           = nbln->cj_ind_start;
        cjind1           = nbln->cj_ind_end;
        ci               = nbln->ci;
        ci_sh            = (ish == CENTRAL ? ci : -1);

        shX_S = gmx_load1_pr(shiftvec+ish3);
        shY_S = gmx_load1_pr(shiftvec+ish3+1);
        shZ_S = gmx_load1_pr(shiftvec+ish3+2);

#if UNROLLJ <= 4
        sci              = ci*STRIDE;
        scix             = sci*DIM;
        sci2             = sci*2;
#else
        sci              = (ci>>1)*STRIDE;
        scix             = sci*DIM + (ci & 1)*(STRIDE>>1);
        sci2             = sci*2 + (ci & 1)*(STRIDE>>1);
        sci             += (ci & 1)*(STRIDE>>1);
#endif

        /* We have 5 LJ/C combinations, but use only three inner loops,
         * as the other combinations are unlikely and/or not much faster:
         * inner half-LJ + C for half-LJ + C / no-LJ + C
         * inner LJ + C      for full-LJ + C
         * inner LJ          for full-LJ + no-C / half-LJ + no-C
         */
        do_LJ   = (nbln->shift & NBNXN_CI_DO_LJ(0));
        do_coul = (nbln->shift & NBNXN_CI_DO_COUL(0));
        half_LJ = ((nbln->shift & NBNXN_CI_HALF_LJ(0)) || !do_LJ) && do_coul;

#ifdef ENERGY_GROUPS
        egps_i = nbat->energrp[ci];
        {
            int ia, egp_ia;

            for (ia = 0; ia < UNROLLI; ia++)
            {
                egp_ia     = (egps_i >> (ia*egps_ishift)) & egps_imask;
                vvdwtp[ia] = Vvdw + egp_ia*Vstride_i;
                vctp[ia]   = Vc   + egp_ia*Vstride_i;
            }
        }
#endif
#if defined CALC_ENERGIES
#if UNROLLJ == 4
        if (do_coul && l_cj[nbln->cj_ind_start].cj == ci_sh)
#endif
#if UNROLLJ == 2
        if (do_coul && l_cj[nbln->cj_ind_start].cj == (ci_sh<<1))
#endif
#if UNROLLJ == 8
        if (do_coul && l_cj[nbln->cj_ind_start].cj == (ci_sh>>1))
#endif
        {
            int  ia;
            real Vc_sub_self;

#ifdef CALC_COUL_RF
            Vc_sub_self = 0.5*ic->c_rf;
#endif
#ifdef CALC_COUL_TAB
#ifdef TAB_FDV0
            Vc_sub_self = 0.5*tab_coul_F[2];
#else
            Vc_sub_self = 0.5*tab_coul_V[0];
#endif
#endif
#ifdef CALC_COUL_EWALD
            /* beta/sqrt(pi) */
            Vc_sub_self = 0.5*ic->ewaldcoeff*M_2_SQRTPI;
#endif

            for (ia = 0; ia < UNROLLI; ia++)
            {
                real qi;

                qi = q[sci+ia];
#ifdef ENERGY_GROUPS
                vctp[ia][((egps_i>>(ia*egps_ishift)) & egps_imask)*egps_jstride]
#else
                Vc[0]
#endif
                    -= facel*qi*qi*Vc_sub_self;
            }
        }
#endif

        /* Load i atom data */
        sciy             = scix + STRIDE;
        sciz             = sciy + STRIDE;
        ix_S0          = gmx_add_pr(gmx_load1_pr(x+scix), shX_S);
        ix_S1          = gmx_add_pr(gmx_load1_pr(x+scix+1), shX_S);
        ix_S2          = gmx_add_pr(gmx_load1_pr(x+scix+2), shX_S);
        ix_S3          = gmx_add_pr(gmx_load1_pr(x+scix+3), shX_S);
        iy_S0          = gmx_add_pr(gmx_load1_pr(x+sciy), shY_S);
        iy_S1          = gmx_add_pr(gmx_load1_pr(x+sciy+1), shY_S);
        iy_S2          = gmx_add_pr(gmx_load1_pr(x+sciy+2), shY_S);
        iy_S3          = gmx_add_pr(gmx_load1_pr(x+sciy+3), shY_S);
        iz_S0          = gmx_add_pr(gmx_load1_pr(x+sciz), shZ_S);
        iz_S1          = gmx_add_pr(gmx_load1_pr(x+sciz+1), shZ_S);
        iz_S2          = gmx_add_pr(gmx_load1_pr(x+sciz+2), shZ_S);
        iz_S3          = gmx_add_pr(gmx_load1_pr(x+sciz+3), shZ_S);

        if (do_coul)
        {
            iq_S0      = gmx_set1_pr(facel*q[sci]);
            iq_S1      = gmx_set1_pr(facel*q[sci+1]);
            iq_S2      = gmx_set1_pr(facel*q[sci+2]);
            iq_S3      = gmx_set1_pr(facel*q[sci+3]);
        }

#ifdef LJ_COMB_LB
        hsig_i_S0      = gmx_load1_pr(ljc+sci2+0);
        hsig_i_S1      = gmx_load1_pr(ljc+sci2+1);
        hsig_i_S2      = gmx_load1_pr(ljc+sci2+2);
        hsig_i_S3      = gmx_load1_pr(ljc+sci2+3);
        seps_i_S0      = gmx_load1_pr(ljc+sci2+STRIDE+0);
        seps_i_S1      = gmx_load1_pr(ljc+sci2+STRIDE+1);
        seps_i_S2      = gmx_load1_pr(ljc+sci2+STRIDE+2);
        seps_i_S3      = gmx_load1_pr(ljc+sci2+STRIDE+3);
#else
#ifdef LJ_COMB_GEOM
        c6s_S0         = gmx_load1_pr(ljc+sci2+0);
        c6s_S1         = gmx_load1_pr(ljc+sci2+1);
        if (!half_LJ)
        {
            c6s_S2     = gmx_load1_pr(ljc+sci2+2);
            c6s_S3     = gmx_load1_pr(ljc+sci2+3);
        }
        c12s_S0        = gmx_load1_pr(ljc+sci2+STRIDE+0);
        c12s_S1        = gmx_load1_pr(ljc+sci2+STRIDE+1);
        if (!half_LJ)
        {
            c12s_S2    = gmx_load1_pr(ljc+sci2+STRIDE+2);
            c12s_S3    = gmx_load1_pr(ljc+sci2+STRIDE+3);
        }
#else
        nbfp0     = nbfp_ptr + type[sci  ]*nbat->ntype*nbfp_stride;
        nbfp1     = nbfp_ptr + type[sci+1]*nbat->ntype*nbfp_stride;
        if (!half_LJ)
        {
            nbfp2 = nbfp_ptr + type[sci+2]*nbat->ntype*nbfp_stride;
            nbfp3 = nbfp_ptr + type[sci+3]*nbat->ntype*nbfp_stride;
        }
#endif
#endif

        /* Zero the potential energy for this list */
        Vvdwtot_S        = gmx_setzero_pr();
        vctot_S          = gmx_setzero_pr();

        /* Clear i atom forces */
        fix_S0           = gmx_setzero_pr();
        fix_S1           = gmx_setzero_pr();
        fix_S2           = gmx_setzero_pr();
        fix_S3           = gmx_setzero_pr();
        fiy_S0           = gmx_setzero_pr();
        fiy_S1           = gmx_setzero_pr();
        fiy_S2           = gmx_setzero_pr();
        fiy_S3           = gmx_setzero_pr();
        fiz_S0           = gmx_setzero_pr();
        fiz_S1           = gmx_setzero_pr();
        fiz_S2           = gmx_setzero_pr();
        fiz_S3           = gmx_setzero_pr();

        cjind = cjind0;

        /* Currently all kernels use (at least half) LJ */
#define CALC_LJ
        if (half_LJ)
        {
#define CALC_COULOMB
#define HALF_LJ
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
            }
#undef HALF_LJ
#undef CALC_COULOMB
        }
        else if (do_coul)
        {
#define CALC_COULOMB
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
            }
#undef CALC_COULOMB
        }
        else
        {
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_4xn_inner.h"
            }
        }
#undef CALC_LJ
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
#define gmx_load_pr4  _mm_load_ps
#define gmx_store_pr4 _mm_store_ps
#define gmx_add_pr4   _mm_add_ps
#else
#define gmx_load_pr4  _mm256_load_pd
#define gmx_store_pr4 _mm256_store_pd
#define gmx_add_pr4   _mm256_add_pd
#endif
        GMX_MM_TRANSPOSE_SUM4_PR(fix_S0, fix_S1, fix_S2, fix_S3, fix_S);
        gmx_store_pr4(f+scix, gmx_add_pr4(fix_S, gmx_load_pr4(f+scix)));

        GMX_MM_TRANSPOSE_SUM4_PR(fiy_S0, fiy_S1, fiy_S2, fiy_S3, fiy_S);
        gmx_store_pr4(f+sciy, gmx_add_pr4(fiy_S, gmx_load_pr4(f+sciy)));

        GMX_MM_TRANSPOSE_SUM4_PR(fiz_S0, fiz_S1, fiz_S2, fiz_S3, fiz_S);
        gmx_store_pr4(f+sciz, gmx_add_pr4(fiz_S, gmx_load_pr4(f+sciz)));

#ifdef CALC_SHIFTFORCES
        gmx_store_pr4(shf, fix_S);
        fshift[ish3+0] += SUM_SIMD4(shf);
        gmx_store_pr4(shf, fiy_S);
        fshift[ish3+1] += SUM_SIMD4(shf);
        gmx_store_pr4(shf, fiz_S);
        fshift[ish3+2] += SUM_SIMD4(shf);
#endif
#else
        GMX_MM_TRANSPOSE_SUM2_PD(fix_S0, fix_S1, fix0_S);
        _mm_store_pd(f+scix, _mm_add_pd(fix0_S, _mm_load_pd(f+scix)));
        GMX_MM_TRANSPOSE_SUM2_PD(fix_S2, fix_S3, fix2_S);
        _mm_store_pd(f+scix+2, _mm_add_pd(fix2_S, _mm_load_pd(f+scix+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiy_S0, fiy_S1, fiy0_S);
        _mm_store_pd(f+sciy, _mm_add_pd(fiy0_S, _mm_load_pd(f+sciy)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiy_S2, fiy_S3, fiy2_S);
        _mm_store_pd(f+sciy+2, _mm_add_pd(fiy2_S, _mm_load_pd(f+sciy+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiz_S0, fiz_S1, fiz0_S);
        _mm_store_pd(f+sciz, _mm_add_pd(fiz0_S, _mm_load_pd(f+sciz)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiz_S2, fiz_S3, fiz2_S);
        _mm_store_pd(f+sciz+2, _mm_add_pd(fiz2_S, _mm_load_pd(f+sciz+2)));

#ifdef CALC_SHIFTFORCES
        _mm_store_pd(shf, _mm_add_pd(fix0_S, fix2_S));
        fshift[ish3+0] += shf[0] + shf[1];
        _mm_store_pd(shf, _mm_add_pd(fiy0_S, fiy2_S));
        fshift[ish3+1] += shf[0] + shf[1];
        _mm_store_pd(shf, _mm_add_pd(fiz0_S, fiz2_S));
        fshift[ish3+2] += shf[0] + shf[1];
#endif
#endif

#ifdef CALC_ENERGIES
        if (do_coul)
        {
            gmx_store_pr(tmpsum, vctot_S);
            *Vc += SUM_SIMD(tmpsum);
        }

        gmx_store_pr(tmpsum, Vvdwtot_S);
        *Vvdw += SUM_SIMD(tmpsum);
#endif

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}


#undef gmx_load_pr4
#undef gmx_store_pr4
#undef gmx_store_pr4

#undef CALC_SHIFTFORCES

#undef UNROLLI
#undef UNROLLJ
#undef STRIDE
#undef TAB_FDV0
#undef NBFP_STRIDE

#undef GMX_USE_HALF_WIDTH_SIMD_HERE
