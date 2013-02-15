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

/* GMX_MM256_HERE should be set before including this file */
#include "gmx_simd_macros.h"

#define SUM_SIMD4(x) (x[0]+x[1]+x[2]+x[3])

#define UNROLLI    NBNXN_CPU_CLUSTER_I_SIZE
#define UNROLLJ    (GMX_SIMD_WIDTH_HERE/2)

#if defined GMX_MM256_HERE
#define STRIDE     4
#endif

#ifdef GMX_MM256_HERE
#ifndef GMX_DOUBLE
/* single precision 2x(4+4) kernel */
#define SUM_SIMD(x) (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])
#define TAB_FDV0
#else
#error "unsupported kernel configuration"
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
NBK_FUNC_NAME(nbnxn_kernel_simd_2xnn, noener)
#else
#ifndef ENERGY_GROUPS
NBK_FUNC_NAME(nbnxn_kernel_simd_2xnn, ener)
#else
NBK_FUNC_NAME(nbnxn_kernel_simd_2xnn, energrp)
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

    gmx_mm_pr  shX_SSE;
    gmx_mm_pr  shY_SSE;
    gmx_mm_pr  shZ_SSE;
    gmx_mm_pr  ix_SSE0, iy_SSE0, iz_SSE0;
    gmx_mm_pr  ix_SSE2, iy_SSE2, iz_SSE2;
    gmx_mm_pr  fix_SSE0, fiy_SSE0, fiz_SSE0;
    gmx_mm_pr  fix_SSE2, fiy_SSE2, fiz_SSE2;
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
    __m128     fix_SSE, fiy_SSE, fiz_SSE;
#else
    __m256d    fix_SSE, fiy_SSE, fiz_SSE;
#endif
#else
    __m128d    fix0_SSE, fiy0_SSE, fiz0_SSE;
    __m128d    fix2_SSE, fiy2_SSE, fiz2_SSE;
#endif

    /* AVX: use floating point masks, as there are no integer instructions */
    gmx_mm_pr  mask0 = _mm256_castsi256_ps(_mm256_set_epi32( 0x0080, 0x0040, 0x0020, 0x0010, 0x0008, 0x0004, 0x0002, 0x0001 ));
    gmx_mm_pr  mask2 = _mm256_castsi256_ps(_mm256_set_epi32( 0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100 ));

    gmx_mm_pr  diag_jmi_SSE;
#if UNROLLI == UNROLLJ
    gmx_mm_pr  diag_SSE0, diag_SSE2;
#else
    gmx_mm_pr  diag0_SSE0, diag0_SSE2;
    gmx_mm_pr  diag1_SSE0, diag1_SSE2;
#endif

    gmx_mm_pr  zero_SSE = gmx_set1_pr(0);

    gmx_mm_pr  one_SSE = gmx_set1_pr(1.0);
    gmx_mm_pr  iq_SSE0 = gmx_setzero_pr();
    gmx_mm_pr  iq_SSE2 = gmx_setzero_pr();
    gmx_mm_pr  mrc_3_SSE;
#ifdef CALC_ENERGIES
    gmx_mm_pr  hrc_3_SSE, moh_rc_SSE;
#endif

#ifdef CALC_COUL_TAB
    /* Coulomb table variables */
    gmx_mm_pr   invtsp_SSE;
    const real *tab_coul_F;
#ifndef TAB_FDV0
    const real *tab_coul_V;
#endif
#ifdef GMX_MM256_HERE
    int        ti0_array[2*GMX_SIMD_WIDTH_HERE-1], *ti0;
    int        ti2_array[2*GMX_SIMD_WIDTH_HERE-1], *ti2;
#endif
#ifdef CALC_ENERGIES
    gmx_mm_pr  mhalfsp_SSE;
#endif
#endif

#ifdef CALC_COUL_EWALD
    gmx_mm_pr beta2_SSE, beta_SSE;
#endif

#if defined CALC_ENERGIES && (defined CALC_COUL_EWALD || defined CALC_COUL_TAB)
    gmx_mm_pr  sh_ewald_SSE;
#endif

#ifdef LJ_COMB_LB
    const real *ljc;

    gmx_mm_pr   hsig_i_SSE0, seps_i_SSE0;
    gmx_mm_pr   hsig_i_SSE2, seps_i_SSE2;
#else
#ifdef FIX_LJ_C
    real        pvdw_array[2*UNROLLI*UNROLLJ+3];
    real       *pvdw_c6, *pvdw_c12;
    gmx_mm_pr   c6_SSE0, c12_SSE0;
    gmx_mm_pr   c6_SSE2, c12_SSE2;
#endif

#ifdef LJ_COMB_GEOM
    const real *ljc;

    gmx_mm_pr   c6s_SSE0, c12s_SSE0;
    gmx_mm_pr   c6s_SSE1, c12s_SSE1;
    gmx_mm_pr   c6s_SSE2 = gmx_setzero_pr(), c12s_SSE2 = gmx_setzero_pr();
    gmx_mm_pr   c6s_SSE3 = gmx_setzero_pr(), c12s_SSE3 = gmx_setzero_pr();
#endif
#endif /* LJ_COMB_LB */

    gmx_mm_pr  vctotSSE, VvdwtotSSE;
    gmx_mm_pr  sixthSSE, twelvethSSE;

    gmx_mm_pr  avoid_sing_SSE;
    gmx_mm_pr  rc2_SSE;
#ifdef VDW_CUTOFF_CHECK
    gmx_mm_pr  rcvdw2_SSE;
#endif

#ifdef CALC_ENERGIES
    gmx_mm_pr  sh_invrc6_SSE, sh_invrc12_SSE;

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
    diag_jmi_SSE = gmx_load_pr(nbat->simd_2xnn_diag);
    /* Generate all the diagonal masks as comparison results */
#if UNROLLI == UNROLLJ
    diag_SSE0    = gmx_cmplt_pr(zero_SSE, diag_jmi_SSE);
    diag_jmi_SSE = gmx_sub_pr(diag_jmi_SSE, one_SSE);
    diag_jmi_SSE = gmx_sub_pr(diag_jmi_SSE, one_SSE);
    diag_SSE2    = gmx_cmplt_pr(zero_SSE, diag_jmi_SSE);
#else
#if 2*UNROLLI == UNROLLJ
    diag0_SSE0 = gmx_cmplt_pr(diag_i_SSE, diag_j_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag0_SSE2 = gmx_cmplt_pr(diag_i_SSE, diag_j_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag1_SSE0 = gmx_cmplt_pr(diag_i_SSE, diag_j_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag_i_SSE = gmx_add_pr(diag_i_SSE, one_SSE);
    diag1_SSE2 = gmx_cmplt_pr(diag_i_SSE, diag_j_SSE);
#endif
#endif

#ifdef CALC_COUL_TAB
#ifdef GMX_MM256_HERE
    /* Generate aligned table index pointers */
    ti0 = (int *)(((size_t)(ti0_array+GMX_SIMD_WIDTH_HERE-1)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(int)-1))));
    ti2 = (int *)(((size_t)(ti2_array+GMX_SIMD_WIDTH_HERE-1)) & (~((size_t)(GMX_SIMD_WIDTH_HERE*sizeof(int)-1))));
#endif

    invtsp_SSE  = gmx_set1_pr(ic->tabq_scale);
#ifdef CALC_ENERGIES
    mhalfsp_SSE = gmx_set1_pr(-0.5/ic->tabq_scale);
#endif

#ifdef TAB_FDV0
    tab_coul_F = ic->tabq_coul_FDV0;
#else
    tab_coul_F = ic->tabq_coul_F;
    tab_coul_V = ic->tabq_coul_V;
#endif
#endif /* CALC_COUL_TAB */

#ifdef CALC_COUL_EWALD
    beta2_SSE = gmx_set1_pr(ic->ewaldcoeff*ic->ewaldcoeff);
    beta_SSE  = gmx_set1_pr(ic->ewaldcoeff);
#endif

#if (defined CALC_COUL_TAB || defined CALC_COUL_EWALD) && defined CALC_ENERGIES
    sh_ewald_SSE = gmx_set1_pr(ic->sh_ewald);
#endif

    q                   = nbat->q;
    type                = nbat->type;
    facel               = ic->epsfac;
    shiftvec            = shift_vec[0];
    x                   = nbat->x;

    avoid_sing_SSE = gmx_set1_pr(NBNXN_AVOID_SING_R2_INC);

    /* The kernel either supports rcoulomb = rvdw or rcoulomb >= rvdw */
    rc2_SSE    = gmx_set1_pr(ic->rcoulomb*ic->rcoulomb);
#ifdef VDW_CUTOFF_CHECK
    rcvdw2_SSE = gmx_set1_pr(ic->rvdw*ic->rvdw);
#endif

#ifdef CALC_ENERGIES
    sixthSSE    = gmx_set1_pr(1.0/6.0);
    twelvethSSE = gmx_set1_pr(1.0/12.0);

    sh_invrc6_SSE  = gmx_set1_pr(ic->sh_invrc6);
    sh_invrc12_SSE = gmx_set1_pr(ic->sh_invrc6*ic->sh_invrc6);
#endif

    mrc_3_SSE = gmx_set1_pr(-2*ic->k_rf);

#ifdef CALC_ENERGIES
    hrc_3_SSE = gmx_set1_pr(ic->k_rf);

    moh_rc_SSE = gmx_set1_pr(-ic->c_rf);
#endif

#ifdef CALC_ENERGIES
    tmpsum = (real *)(((size_t)(tmpsum_array+7)) & (~((size_t)31)));
#endif
#ifdef CALC_SHIFTFORCES
    shf = (real *)(((size_t)(shf_array+7)) & (~((size_t)31)));
#endif

#ifdef FIX_LJ_C
    pvdw_c6  = (real *)(((size_t)(pvdw_array+3)) & (~((size_t)15)));
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
    c6_SSE0            = gmx_load_pr(pvdw_c6 +0*UNROLLJ);
    c6_SSE1            = gmx_load_pr(pvdw_c6 +1*UNROLLJ);
    c6_SSE2            = gmx_load_pr(pvdw_c6 +2*UNROLLJ);
    c6_SSE3            = gmx_load_pr(pvdw_c6 +3*UNROLLJ);

    c12_SSE0           = gmx_load_pr(pvdw_c12+0*UNROLLJ);
    c12_SSE1           = gmx_load_pr(pvdw_c12+1*UNROLLJ);
    c12_SSE2           = gmx_load_pr(pvdw_c12+2*UNROLLJ);
    c12_SSE3           = gmx_load_pr(pvdw_c12+3*UNROLLJ);
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

        shX_SSE = gmx_load1_pr(shiftvec+ish3);
        shY_SSE = gmx_load1_pr(shiftvec+ish3+1);
        shZ_SSE = gmx_load1_pr(shiftvec+ish3+2);

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

#define gmx_load2_hpr(x)  _mm256_insertf128_ps(gmx_load1_pr(x), gmx_load1_hpr(x+1), 1)

        /* Load i atom data */
        sciy             = scix + STRIDE;
        sciz             = sciy + STRIDE;
        ix_SSE0          = gmx_add_pr(gmx_load2_hpr(x+scix), shX_SSE);
        ix_SSE2          = gmx_add_pr(gmx_load2_hpr(x+scix+2), shX_SSE);
        iy_SSE0          = gmx_add_pr(gmx_load2_hpr(x+sciy), shY_SSE);
        iy_SSE2          = gmx_add_pr(gmx_load2_hpr(x+sciy+2), shY_SSE);
        iz_SSE0          = gmx_add_pr(gmx_load2_hpr(x+sciz), shZ_SSE);
        iz_SSE2          = gmx_add_pr(gmx_load2_hpr(x+sciz+2), shZ_SSE);

        if (do_coul)
        {
            gmx_mm_pr facel_SSE;

            facel_SSE    = gmx_set1_pr(facel);

            iq_SSE0      = gmx_mul_pr(facel_SSE, gmx_load2_hpr(q+sci));
            iq_SSE2      = gmx_mul_pr(facel_SSE, gmx_load2_hpr(q+sci+2));
        }

#ifdef LJ_COMB_LB
        hsig_i_SSE0      = gmx_load2_hpr(ljc+sci2+0);
        hsig_i_SSE2      = gmx_load2_hpr(ljc+sci2+2);
        seps_i_SSE0      = gmx_load2_hpr(ljc+sci2+STRIDE+0);
        seps_i_SSE2      = gmx_load2_hpr(ljc+sci2+STRIDE+2);
#else
#ifdef LJ_COMB_GEOM
        c6s_SSE0         = gmx_load2_hpr(ljc+sci2+0);
        if (!half_LJ)
        {
            c6s_SSE2     = gmx_load2_hpr(ljc+sci2+2);
        }
        c12s_SSE0        = gmx_load2_hpr(ljc+sci2+STRIDE+0);
        if (!half_LJ)
        {
            c12s_SSE2    = gmx_load2_hpr(ljc+sci2+STRIDE+2);
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
        VvdwtotSSE       = gmx_setzero_pr();
        vctotSSE         = gmx_setzero_pr();

        /* Clear i atom forces */
        fix_SSE0           = gmx_setzero_pr();
        fix_SSE2           = gmx_setzero_pr();
        fiy_SSE0           = gmx_setzero_pr();
        fiy_SSE2           = gmx_setzero_pr();
        fiz_SSE0           = gmx_setzero_pr();
        fiz_SSE2           = gmx_setzero_pr();

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
#include "nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_2xnn_inner.h"
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
#include "nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_2xnn_inner.h"
            }
#undef CALC_COULOMB
        }
        else
        {
#define CHECK_EXCLS
            while (cjind < cjind1 && nbl->cj[cjind].excl != SIMD_MASK_ALL)
            {
#include "nbnxn_kernel_simd_2xnn_inner.h"
                cjind++;
            }
#undef CHECK_EXCLS
            for (; (cjind < cjind1); cjind++)
            {
#include "nbnxn_kernel_simd_2xnn_inner.h"
            }
        }
#undef CALC_LJ
        ninner += cjind1 - cjind0;

        /* Add accumulated i-forces to the force array */
#if UNROLLJ >= 4
#ifndef GMX_DOUBLE
#define gmx_load_ps4  _mm_load_ps
#define gmx_store_ps4 _mm_store_ps
#define gmx_add_ps4   _mm_add_ps
#else
#define gmx_load_ps4  _mm256_load_pd
#define gmx_store_ps4 _mm256_store_pd
#define gmx_add_ps4   _mm256_add_pd
#endif
        GMX_MM_TRANSPOSE_SUM4H_PR(fix_SSE0, fix_SSE2, fix_SSE);
        gmx_store_ps4(f+scix, gmx_add_ps4(fix_SSE, gmx_load_ps4(f+scix)));

        GMX_MM_TRANSPOSE_SUM4H_PR(fiy_SSE0, fiy_SSE2, fiy_SSE);
        gmx_store_ps4(f+sciy, gmx_add_ps4(fiy_SSE, gmx_load_ps4(f+sciy)));

        GMX_MM_TRANSPOSE_SUM4H_PR(fiz_SSE0, fiz_SSE2, fiz_SSE);
        gmx_store_ps4(f+sciz, gmx_add_ps4(fiz_SSE, gmx_load_ps4(f+sciz)));

#ifdef CALC_SHIFTFORCES
        gmx_store_ps4(shf, fix_SSE);
        fshift[ish3+0] += SUM_SIMD4(shf);
        gmx_store_ps4(shf, fiy_SSE);
        fshift[ish3+1] += SUM_SIMD4(shf);
        gmx_store_ps4(shf, fiz_SSE);
        fshift[ish3+2] += SUM_SIMD4(shf);
#endif
#else
        GMX_MM_TRANSPOSE_SUM2_PD(fix_SSE0, fix_SSE1, fix0_SSE);
        _mm_store_pd(f+scix, _mm_add_pd(fix0_SSE, _mm_load_pd(f+scix)));
        GMX_MM_TRANSPOSE_SUM2_PD(fix_SSE2, fix_SSE3, fix2_SSE);
        _mm_store_pd(f+scix+2, _mm_add_pd(fix2_SSE, _mm_load_pd(f+scix+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiy_SSE0, fiy_SSE1, fiy0_SSE);
        _mm_store_pd(f+sciy, _mm_add_pd(fiy0_SSE, _mm_load_pd(f+sciy)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiy_SSE2, fiy_SSE3, fiy2_SSE);
        _mm_store_pd(f+sciy+2, _mm_add_pd(fiy2_SSE, _mm_load_pd(f+sciy+2)));

        GMX_MM_TRANSPOSE_SUM2_PD(fiz_SSE0, fiz_SSE1, fiz0_SSE);
        _mm_store_pd(f+sciz, _mm_add_pd(fiz0_SSE, _mm_load_pd(f+sciz)));
        GMX_MM_TRANSPOSE_SUM2_PD(fiz_SSE2, fiz_SSE3, fiz2_SSE);
        _mm_store_pd(f+sciz+2, _mm_add_pd(fiz2_SSE, _mm_load_pd(f+sciz+2)));

#ifdef CALC_SHIFTFORCES
        _mm_store_pd(shf, _mm_add_pd(fix0_SSE, fix2_SSE));
        fshift[ish3+0] += shf[0] + shf[1];
        _mm_store_pd(shf, _mm_add_pd(fiy0_SSE, fiy2_SSE));
        fshift[ish3+1] += shf[0] + shf[1];
        _mm_store_pd(shf, _mm_add_pd(fiz0_SSE, fiz2_SSE));
        fshift[ish3+2] += shf[0] + shf[1];
#endif
#endif

#ifdef CALC_ENERGIES
        if (do_coul)
        {
            gmx_store_pr(tmpsum, vctotSSE);
            *Vc += SUM_SIMD(tmpsum);
        }

        gmx_store_pr(tmpsum, VvdwtotSSE);
        *Vvdw += SUM_SIMD(tmpsum);
#endif

        /* Outer loop uses 6 flops/iteration */
    }

#ifdef COUNT_PAIRS
    printf("atom pairs %d\n", npair);
#endif
}

#undef gmx_load2_hpr

#undef gmx_load_ps4
#undef gmx_store_ps4
#undef gmx_store_ps4

#undef CALC_SHIFTFORCES

#undef UNROLLI
#undef UNROLLJ
#undef STRIDE
#undef TAB_FDV0
#undef NBFP_STRIDE
