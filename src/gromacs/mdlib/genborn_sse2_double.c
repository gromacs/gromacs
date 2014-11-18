/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <math.h>
#include <string.h>

#include "gromacs/domdec/domdec.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/legacyheaders/genborn.h"
#include "gromacs/legacyheaders/names.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/smalloc.h"

/* Only compile this file if SSE2 intrinsics are available */
#if 0 && defined (GMX_SIMD_X86_SSE2_OR_HIGHER)
#include "genborn_sse2_double.h"

#include <emmintrin.h>
#include <gmx_sse2_double.h>

int
calc_gb_rad_still_sse2_double(t_commrec *cr, t_forcerec *fr,
                              int natoms, gmx_localtop_t *top,
                              double *x, t_nblist *nl,
                              gmx_genborn_t *born)
{
    int           i, k, n, ii, is3, ii3, nj0, nj1, offset;
    int           jnrA, jnrB, j3A, j3B;
    int          *mdtype;
    double        shX, shY, shZ;
    int          *jjnr;
    double       *shiftvec;

    double        gpi_ai, gpi2;
    double        factor;
    double       *gb_radius;
    double       *vsolv;
    double       *work;
    double       *dadx;

    __m128d       ix, iy, iz;
    __m128d       jx, jy, jz;
    __m128d       dx, dy, dz;
    __m128d       tx, ty, tz;
    __m128d       rsq, rinv, rinv2, rinv4, rinv6;
    __m128d       ratio, gpi, rai, raj, vai, vaj, rvdw;
    __m128d       ccf, dccf, theta, cosq, term, sinq, res, prod, prod_ai, tmp;
    __m128d       mask, icf4, icf6, mask_cmp;

    const __m128d half   = _mm_set1_pd(0.5);
    const __m128d three  = _mm_set1_pd(3.0);
    const __m128d one    = _mm_set1_pd(1.0);
    const __m128d two    = _mm_set1_pd(2.0);
    const __m128d zero   = _mm_set1_pd(0.0);
    const __m128d four   = _mm_set1_pd(4.0);

    const __m128d still_p5inv  = _mm_set1_pd(STILL_P5INV);
    const __m128d still_pip5   = _mm_set1_pd(STILL_PIP5);
    const __m128d still_p4     = _mm_set1_pd(STILL_P4);

    factor  = 0.5 * ONE_4PI_EPS0;

    gb_radius = born->gb_radius;
    vsolv     = born->vsolv;
    work      = born->gpol_still_work;
    jjnr      = nl->jjnr;
    shiftvec  = fr->shift_vec[0];
    dadx      = fr->dadx;

    jnrA = jnrB = 0;
    jx   = _mm_setzero_pd();
    jy   = _mm_setzero_pd();
    jz   = _mm_setzero_pd();

    n = 0;

    for (i = 0; i < natoms; i++)
    {
        work[i] = 0;
    }

    for (i = 0; i < nl->nri; i++)
    {
        ii     = nl->iinr[i];
        ii3    = ii*3;
        is3    = 3*nl->shift[i];
        shX    = shiftvec[is3];
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];
        nj1    = nl->jindex[i+1];

        ix     = _mm_set1_pd(shX+x[ii3+0]);
        iy     = _mm_set1_pd(shY+x[ii3+1]);
        iz     = _mm_set1_pd(shZ+x[ii3+2]);


        /* Polarization energy for atom ai */
        gpi    = _mm_setzero_pd();

        rai     = _mm_load1_pd(gb_radius+ii);
        prod_ai = _mm_set1_pd(STILL_P4*vsolv[ii]);

        for (k = nj0; k < nj1-1; k += 2)
        {
            jnrA        = jjnr[k];
            jnrB        = jjnr[k+1];

            j3A         = 3*jnrA;
            j3B         = 3*jnrB;

            GMX_MM_LOAD_1RVEC_2POINTERS_PD(x+j3A, x+j3B, jx, jy, jz);

            GMX_MM_LOAD_2VALUES_PD(gb_radius+jnrA, gb_radius+jnrB, raj);
            GMX_MM_LOAD_2VALUES_PD(vsolv+jnrA, vsolv+jnrB, vaj);

            dx          = _mm_sub_pd(ix, jx);
            dy          = _mm_sub_pd(iy, jy);
            dz          = _mm_sub_pd(iz, jz);

            rsq         = gmx_mm_calc_rsq_pd(dx, dy, dz);
            rinv        = gmx_mm_invsqrt_pd(rsq);
            rinv2       = _mm_mul_pd(rinv, rinv);
            rinv4       = _mm_mul_pd(rinv2, rinv2);
            rinv6       = _mm_mul_pd(rinv4, rinv2);

            rvdw        = _mm_add_pd(rai, raj);
            ratio       = _mm_mul_pd(rsq, gmx_mm_inv_pd( _mm_mul_pd(rvdw, rvdw)));

            mask_cmp    = _mm_cmple_pd(ratio, still_p5inv);

            /* gmx_mm_sincos_pd() is quite expensive, so avoid calculating it if we can! */
            if (0 == _mm_movemask_pd(mask_cmp) )
            {
                /* if ratio>still_p5inv for ALL elements */
                ccf         = one;
                dccf        = _mm_setzero_pd();
            }
            else
            {
                ratio       = _mm_min_pd(ratio, still_p5inv);
                theta       = _mm_mul_pd(ratio, still_pip5);
                gmx_mm_sincos_pd(theta, &sinq, &cosq);
                term        = _mm_mul_pd(half, _mm_sub_pd(one, cosq));
                ccf         = _mm_mul_pd(term, term);
                dccf        = _mm_mul_pd(_mm_mul_pd(two, term),
                                         _mm_mul_pd(sinq, theta));
            }

            prod        = _mm_mul_pd(still_p4, vaj);
            icf4        = _mm_mul_pd(ccf, rinv4);
            icf6        = _mm_mul_pd( _mm_sub_pd( _mm_mul_pd(four, ccf), dccf), rinv6);

            GMX_MM_INCREMENT_2VALUES_PD(work+jnrA, work+jnrB, _mm_mul_pd(prod_ai, icf4));

            gpi           = _mm_add_pd(gpi, _mm_mul_pd(prod, icf4) );

            _mm_store_pd(dadx, _mm_mul_pd(prod, icf6));
            dadx += 2;
            _mm_store_pd(dadx, _mm_mul_pd(prod_ai, icf6));
            dadx += 2;
        }

        if (k < nj1)
        {
            jnrA        = jjnr[k];

            j3A         = 3*jnrA;

            GMX_MM_LOAD_1RVEC_1POINTER_PD(x+j3A, jx, jy, jz);

            GMX_MM_LOAD_1VALUE_PD(gb_radius+jnrA, raj);
            GMX_MM_LOAD_1VALUE_PD(vsolv+jnrA, vaj);

            dx          = _mm_sub_sd(ix, jx);
            dy          = _mm_sub_sd(iy, jy);
            dz          = _mm_sub_sd(iz, jz);

            rsq         = gmx_mm_calc_rsq_pd(dx, dy, dz);
            rinv        = gmx_mm_invsqrt_pd(rsq);
            rinv2       = _mm_mul_sd(rinv, rinv);
            rinv4       = _mm_mul_sd(rinv2, rinv2);
            rinv6       = _mm_mul_sd(rinv4, rinv2);

            rvdw        = _mm_add_sd(rai, raj);
            ratio       = _mm_mul_sd(rsq, gmx_mm_inv_pd( _mm_mul_pd(rvdw, rvdw)));

            mask_cmp    = _mm_cmple_sd(ratio, still_p5inv);

            /* gmx_mm_sincos_pd() is quite expensive, so avoid calculating it if we can! */
            if (0 == _mm_movemask_pd(mask_cmp) )
            {
                /* if ratio>still_p5inv for ALL elements */
                ccf         = one;
                dccf        = _mm_setzero_pd();
            }
            else
            {
                ratio       = _mm_min_sd(ratio, still_p5inv);
                theta       = _mm_mul_sd(ratio, still_pip5);
                gmx_mm_sincos_pd(theta, &sinq, &cosq);
                term        = _mm_mul_sd(half, _mm_sub_sd(one, cosq));
                ccf         = _mm_mul_sd(term, term);
                dccf        = _mm_mul_sd(_mm_mul_sd(two, term),
                                         _mm_mul_sd(sinq, theta));
            }

            prod        = _mm_mul_sd(still_p4, vaj);
            icf4        = _mm_mul_sd(ccf, rinv4);
            icf6        = _mm_mul_sd( _mm_sub_sd( _mm_mul_sd(four, ccf), dccf), rinv6);

            GMX_MM_INCREMENT_1VALUE_PD(work+jnrA, _mm_mul_sd(prod_ai, icf4));

            gpi           = _mm_add_sd(gpi, _mm_mul_sd(prod, icf4) );

            _mm_store_pd(dadx, _mm_mul_pd(prod, icf6));
            dadx += 2;
            _mm_store_pd(dadx, _mm_mul_pd(prod_ai, icf6));
            dadx += 2;
        }
        gmx_mm_update_1pot_pd(gpi, work+ii);
    }

    /* Sum up the polarization energy from other nodes */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, work);
    }

    /* Compute the radii */
    for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
    {
        if (born->use[i] != 0)
        {
            gpi_ai           = born->gpol[i] + work[i]; /* add gpi to the initial pol energy gpi_ai*/
            gpi2             = gpi_ai * gpi_ai;
            born->bRad[i]    = factor*gmx_invsqrt(gpi2);
            fr->invsqrta[i]  = gmx_invsqrt(born->bRad[i]);
        }
    }

    /* Extra (local) communication required for DD */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_spread_real(cr->dd, born->bRad);
        dd_atom_spread_real(cr->dd, fr->invsqrta);
    }

    return 0;
}


int
calc_gb_rad_hct_obc_sse2_double(t_commrec *cr, t_forcerec * fr, int natoms, gmx_localtop_t *top,
                                double *x, t_nblist *nl, gmx_genborn_t *born, t_mdatoms *md, int gb_algorithm)
{
    int           i, ai, k, n, ii, ii3, is3, nj0, nj1, at0, at1, offset;
    int           jnrA, jnrB;
    int           j3A, j3B;
    double        shX, shY, shZ;
    double        rr, rr_inv, rr_inv2, sum_tmp, sum, sum2, sum3, gbr;
    double        sum_ai2, sum_ai3, tsum, tchain, doffset;
    double       *obc_param;
    double       *gb_radius;
    double       *work;
    int        *  jjnr;
    double       *dadx;
    double       *shiftvec;
    double        min_rad, rad;

    __m128d       ix, iy, iz, jx, jy, jz;
    __m128d       dx, dy, dz, t1, t2, t3, t4;
    __m128d       rsq, rinv, r;
    __m128d       rai, rai_inv, raj, raj_inv, rai_inv2, sk, sk2, lij, dlij, duij;
    __m128d       uij, lij2, uij2, lij3, uij3, diff2;
    __m128d       lij_inv, sk2_inv, prod, log_term, tmp, tmp_sum;
    __m128d       sum_ai, tmp_ai, sk_ai, sk_aj, sk2_ai, sk2_aj, sk2_rinv;
    __m128d       dadx1, dadx2;
    __m128d       logterm;
    __m128d       mask;
    __m128d       obc_mask1, obc_mask2, obc_mask3;

    __m128d       oneeighth   = _mm_set1_pd(0.125);
    __m128d       onefourth   = _mm_set1_pd(0.25);

    const __m128d half  = _mm_set1_pd(0.5);
    const __m128d three = _mm_set1_pd(3.0);
    const __m128d one   = _mm_set1_pd(1.0);
    const __m128d two   = _mm_set1_pd(2.0);
    const __m128d zero  = _mm_set1_pd(0.0);
    const __m128d neg   = _mm_set1_pd(-1.0);

    /* Set the dielectric offset */
    doffset   = born->gb_doffset;
    gb_radius = born->gb_radius;
    obc_param = born->param;
    work      = born->gpol_hct_work;
    jjnr      = nl->jjnr;
    dadx      = fr->dadx;
    shiftvec  = fr->shift_vec[0];

    jx        = _mm_setzero_pd();
    jy        = _mm_setzero_pd();
    jz        = _mm_setzero_pd();

    jnrA = jnrB = 0;

    for (i = 0; i < born->nr; i++)
    {
        work[i] = 0;
    }

    for (i = 0; i < nl->nri; i++)
    {
        ii     = nl->iinr[i];
        ii3    = ii*3;
        is3    = 3*nl->shift[i];
        shX    = shiftvec[is3];
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];
        nj1    = nl->jindex[i+1];

        ix     = _mm_set1_pd(shX+x[ii3+0]);
        iy     = _mm_set1_pd(shY+x[ii3+1]);
        iz     = _mm_set1_pd(shZ+x[ii3+2]);

        rai     = _mm_load1_pd(gb_radius+ii);
        rai_inv = gmx_mm_inv_pd(rai);

        sum_ai = _mm_setzero_pd();

        sk_ai  = _mm_load1_pd(born->param+ii);
        sk2_ai = _mm_mul_pd(sk_ai, sk_ai);

        for (k = nj0; k < nj1-1; k += 2)
        {
            jnrA        = jjnr[k];
            jnrB        = jjnr[k+1];

            j3A         = 3*jnrA;
            j3B         = 3*jnrB;

            GMX_MM_LOAD_1RVEC_2POINTERS_PD(x+j3A, x+j3B, jx, jy, jz);
            GMX_MM_LOAD_2VALUES_PD(gb_radius+jnrA, gb_radius+jnrB, raj);
            GMX_MM_LOAD_2VALUES_PD(obc_param+jnrA, obc_param+jnrB, sk_aj);

            dx    = _mm_sub_pd(ix, jx);
            dy    = _mm_sub_pd(iy, jy);
            dz    = _mm_sub_pd(iz, jz);

            rsq         = gmx_mm_calc_rsq_pd(dx, dy, dz);

            rinv        = gmx_mm_invsqrt_pd(rsq);
            r           = _mm_mul_pd(rsq, rinv);

            /* Compute raj_inv aj1-4 */
            raj_inv     = gmx_mm_inv_pd(raj);

            /* Evaluate influence of atom aj -> ai */
            t1            = _mm_add_pd(r, sk_aj);
            t2            = _mm_sub_pd(r, sk_aj);
            t3            = _mm_sub_pd(sk_aj, r);
            obc_mask1     = _mm_cmplt_pd(rai, t1);
            obc_mask2     = _mm_cmplt_pd(rai, t2);
            obc_mask3     = _mm_cmplt_pd(rai, t3);

            uij           = gmx_mm_inv_pd(t1);
            lij           = _mm_or_pd(   _mm_and_pd(obc_mask2, gmx_mm_inv_pd(t2)),
                                         _mm_andnot_pd(obc_mask2, rai_inv));
            dlij          = _mm_and_pd(one, obc_mask2);
            uij2          = _mm_mul_pd(uij, uij);
            uij3          = _mm_mul_pd(uij2, uij);
            lij2          = _mm_mul_pd(lij, lij);
            lij3          = _mm_mul_pd(lij2, lij);

            diff2         = _mm_sub_pd(uij2, lij2);
            lij_inv       = gmx_mm_invsqrt_pd(lij2);
            sk2_aj        = _mm_mul_pd(sk_aj, sk_aj);
            sk2_rinv      = _mm_mul_pd(sk2_aj, rinv);
            prod          = _mm_mul_pd(onefourth, sk2_rinv);

            logterm       = gmx_mm_log_pd(_mm_mul_pd(uij, lij_inv));

            t1            = _mm_sub_pd(lij, uij);
            t2            = _mm_mul_pd(diff2,
                                       _mm_sub_pd(_mm_mul_pd(onefourth, r),
                                                  prod));
            t3            = _mm_mul_pd(half, _mm_mul_pd(rinv, logterm));
            t1            = _mm_add_pd(t1, _mm_add_pd(t2, t3));
            t4            = _mm_mul_pd(two, _mm_sub_pd(rai_inv, lij));
            t4            = _mm_and_pd(t4, obc_mask3);
            t1            = _mm_mul_pd(half, _mm_add_pd(t1, t4));

            sum_ai        = _mm_add_pd(sum_ai, _mm_and_pd(t1, obc_mask1) );

            t1            = _mm_add_pd(_mm_mul_pd(half, lij2),
                                       _mm_mul_pd(prod, lij3));
            t1            = _mm_sub_pd(t1,
                                       _mm_mul_pd(onefourth,
                                                  _mm_add_pd(_mm_mul_pd(lij, rinv),
                                                             _mm_mul_pd(lij3, r))));
            t2            = _mm_mul_pd(onefourth,
                                       _mm_add_pd(_mm_mul_pd(uij, rinv),
                                                  _mm_mul_pd(uij3, r)));
            t2            = _mm_sub_pd(t2,
                                       _mm_add_pd(_mm_mul_pd(half, uij2),
                                                  _mm_mul_pd(prod, uij3)));
            t3            = _mm_mul_pd(_mm_mul_pd(onefourth, logterm),
                                       _mm_mul_pd(rinv, rinv));
            t3            = _mm_sub_pd(t3,
                                       _mm_mul_pd(_mm_mul_pd(diff2, oneeighth),
                                                  _mm_add_pd(one,
                                                             _mm_mul_pd(sk2_rinv, rinv))));
            t1            = _mm_mul_pd(rinv,
                                       _mm_add_pd(_mm_mul_pd(dlij, t1),
                                                  _mm_add_pd(t2, t3)));

            dadx1         = _mm_and_pd(t1, obc_mask1);

            /* Evaluate influence of atom ai -> aj */
            t1            = _mm_add_pd(r, sk_ai);
            t2            = _mm_sub_pd(r, sk_ai);
            t3            = _mm_sub_pd(sk_ai, r);
            obc_mask1     = _mm_cmplt_pd(raj, t1);
            obc_mask2     = _mm_cmplt_pd(raj, t2);
            obc_mask3     = _mm_cmplt_pd(raj, t3);

            uij           = gmx_mm_inv_pd(t1);
            lij           = _mm_or_pd(   _mm_and_pd(obc_mask2, gmx_mm_inv_pd(t2)),
                                         _mm_andnot_pd(obc_mask2, raj_inv));
            dlij          = _mm_and_pd(one, obc_mask2);
            uij2          = _mm_mul_pd(uij, uij);
            uij3          = _mm_mul_pd(uij2, uij);
            lij2          = _mm_mul_pd(lij, lij);
            lij3          = _mm_mul_pd(lij2, lij);

            diff2         = _mm_sub_pd(uij2, lij2);
            lij_inv       = gmx_mm_invsqrt_pd(lij2);
            sk2_rinv      = _mm_mul_pd(sk2_ai, rinv);
            prod          = _mm_mul_pd(onefourth, sk2_rinv);

            logterm       = gmx_mm_log_pd(_mm_mul_pd(uij, lij_inv));

            t1            = _mm_sub_pd(lij, uij);
            t2            = _mm_mul_pd(diff2,
                                       _mm_sub_pd(_mm_mul_pd(onefourth, r),
                                                  prod));
            t3            = _mm_mul_pd(half, _mm_mul_pd(rinv, logterm));
            t1            = _mm_add_pd(t1, _mm_add_pd(t2, t3));
            t4            = _mm_mul_pd(two, _mm_sub_pd(raj_inv, lij));
            t4            = _mm_and_pd(t4, obc_mask3);
            t1            = _mm_mul_pd(half, _mm_add_pd(t1, t4));

            GMX_MM_INCREMENT_2VALUES_PD(work+jnrA, work+jnrB, _mm_and_pd(t1, obc_mask1));

            t1            = _mm_add_pd(_mm_mul_pd(half, lij2),
                                       _mm_mul_pd(prod, lij3));
            t1            = _mm_sub_pd(t1,
                                       _mm_mul_pd(onefourth,
                                                  _mm_add_pd(_mm_mul_pd(lij, rinv),
                                                             _mm_mul_pd(lij3, r))));
            t2            = _mm_mul_pd(onefourth,
                                       _mm_add_pd(_mm_mul_pd(uij, rinv),
                                                  _mm_mul_pd(uij3, r)));
            t2            = _mm_sub_pd(t2,
                                       _mm_add_pd(_mm_mul_pd(half, uij2),
                                                  _mm_mul_pd(prod, uij3)));
            t3            = _mm_mul_pd(_mm_mul_pd(onefourth, logterm),
                                       _mm_mul_pd(rinv, rinv));
            t3            = _mm_sub_pd(t3,
                                       _mm_mul_pd(_mm_mul_pd(diff2, oneeighth),
                                                  _mm_add_pd(one,
                                                             _mm_mul_pd(sk2_rinv, rinv))));
            t1            = _mm_mul_pd(rinv,
                                       _mm_add_pd(_mm_mul_pd(dlij, t1),
                                                  _mm_add_pd(t2, t3)));

            dadx2         = _mm_and_pd(t1, obc_mask1);

            _mm_store_pd(dadx, dadx1);
            dadx += 2;
            _mm_store_pd(dadx, dadx2);
            dadx += 2;
        } /* end normal inner loop */

        if (k < nj1)
        {
            jnrA        = jjnr[k];

            j3A         = 3*jnrA;

            GMX_MM_LOAD_1RVEC_1POINTER_PD(x+j3A, jx, jy, jz);
            GMX_MM_LOAD_1VALUE_PD(gb_radius+jnrA, raj);
            GMX_MM_LOAD_1VALUE_PD(obc_param+jnrA, sk_aj);

            dx    = _mm_sub_sd(ix, jx);
            dy    = _mm_sub_sd(iy, jy);
            dz    = _mm_sub_sd(iz, jz);

            rsq         = gmx_mm_calc_rsq_pd(dx, dy, dz);

            rinv        = gmx_mm_invsqrt_pd(rsq);
            r           = _mm_mul_sd(rsq, rinv);

            /* Compute raj_inv aj1-4 */
            raj_inv     = gmx_mm_inv_pd(raj);

            /* Evaluate influence of atom aj -> ai */
            t1            = _mm_add_sd(r, sk_aj);
            t2            = _mm_sub_sd(r, sk_aj);
            t3            = _mm_sub_sd(sk_aj, r);
            obc_mask1     = _mm_cmplt_sd(rai, t1);
            obc_mask2     = _mm_cmplt_sd(rai, t2);
            obc_mask3     = _mm_cmplt_sd(rai, t3);

            uij           = gmx_mm_inv_pd(t1);
            lij           = _mm_or_pd(_mm_and_pd(obc_mask2, gmx_mm_inv_pd(t2)),
                                      _mm_andnot_pd(obc_mask2, rai_inv));
            dlij          = _mm_and_pd(one, obc_mask2);
            uij2          = _mm_mul_sd(uij, uij);
            uij3          = _mm_mul_sd(uij2, uij);
            lij2          = _mm_mul_sd(lij, lij);
            lij3          = _mm_mul_sd(lij2, lij);

            diff2         = _mm_sub_sd(uij2, lij2);
            lij_inv       = gmx_mm_invsqrt_pd(lij2);
            sk2_aj        = _mm_mul_sd(sk_aj, sk_aj);
            sk2_rinv      = _mm_mul_sd(sk2_aj, rinv);
            prod          = _mm_mul_sd(onefourth, sk2_rinv);

            logterm       = gmx_mm_log_pd(_mm_mul_sd(uij, lij_inv));

            t1            = _mm_sub_sd(lij, uij);
            t2            = _mm_mul_sd(diff2,
                                       _mm_sub_sd(_mm_mul_pd(onefourth, r),
                                                  prod));
            t3            = _mm_mul_sd(half, _mm_mul_sd(rinv, logterm));
            t1            = _mm_add_sd(t1, _mm_add_sd(t2, t3));
            t4            = _mm_mul_sd(two, _mm_sub_sd(rai_inv, lij));
            t4            = _mm_and_pd(t4, obc_mask3);
            t1            = _mm_mul_sd(half, _mm_add_sd(t1, t4));

            sum_ai        = _mm_add_sd(sum_ai, _mm_and_pd(t1, obc_mask1) );

            t1            = _mm_add_sd(_mm_mul_sd(half, lij2),
                                       _mm_mul_sd(prod, lij3));
            t1            = _mm_sub_sd(t1,
                                       _mm_mul_sd(onefourth,
                                                  _mm_add_sd(_mm_mul_sd(lij, rinv),
                                                             _mm_mul_sd(lij3, r))));
            t2            = _mm_mul_sd(onefourth,
                                       _mm_add_sd(_mm_mul_sd(uij, rinv),
                                                  _mm_mul_sd(uij3, r)));
            t2            = _mm_sub_sd(t2,
                                       _mm_add_sd(_mm_mul_sd(half, uij2),
                                                  _mm_mul_sd(prod, uij3)));
            t3            = _mm_mul_sd(_mm_mul_sd(onefourth, logterm),
                                       _mm_mul_sd(rinv, rinv));
            t3            = _mm_sub_sd(t3,
                                       _mm_mul_sd(_mm_mul_sd(diff2, oneeighth),
                                                  _mm_add_sd(one,
                                                             _mm_mul_sd(sk2_rinv, rinv))));
            t1            = _mm_mul_sd(rinv,
                                       _mm_add_sd(_mm_mul_sd(dlij, t1),
                                                  _mm_add_pd(t2, t3)));

            dadx1         = _mm_and_pd(t1, obc_mask1);

            /* Evaluate influence of atom ai -> aj */
            t1            = _mm_add_sd(r, sk_ai);
            t2            = _mm_sub_sd(r, sk_ai);
            t3            = _mm_sub_sd(sk_ai, r);
            obc_mask1     = _mm_cmplt_sd(raj, t1);
            obc_mask2     = _mm_cmplt_sd(raj, t2);
            obc_mask3     = _mm_cmplt_sd(raj, t3);

            uij           = gmx_mm_inv_pd(t1);
            lij           = _mm_or_pd(   _mm_and_pd(obc_mask2, gmx_mm_inv_pd(t2)),
                                         _mm_andnot_pd(obc_mask2, raj_inv));
            dlij          = _mm_and_pd(one, obc_mask2);
            uij2          = _mm_mul_sd(uij, uij);
            uij3          = _mm_mul_sd(uij2, uij);
            lij2          = _mm_mul_sd(lij, lij);
            lij3          = _mm_mul_sd(lij2, lij);

            diff2         = _mm_sub_sd(uij2, lij2);
            lij_inv       = gmx_mm_invsqrt_pd(lij2);
            sk2_rinv      = _mm_mul_sd(sk2_ai, rinv);
            prod          = _mm_mul_sd(onefourth, sk2_rinv);

            logterm       = gmx_mm_log_pd(_mm_mul_sd(uij, lij_inv));

            t1            = _mm_sub_sd(lij, uij);
            t2            = _mm_mul_sd(diff2,
                                       _mm_sub_sd(_mm_mul_sd(onefourth, r),
                                                  prod));
            t3            = _mm_mul_sd(half, _mm_mul_sd(rinv, logterm));
            t1            = _mm_add_sd(t1, _mm_add_sd(t2, t3));
            t4            = _mm_mul_sd(two, _mm_sub_sd(raj_inv, lij));
            t4            = _mm_and_pd(t4, obc_mask3);
            t1            = _mm_mul_sd(half, _mm_add_sd(t1, t4));

            GMX_MM_INCREMENT_1VALUE_PD(work+jnrA, _mm_and_pd(t1, obc_mask1));

            t1            = _mm_add_sd(_mm_mul_sd(half, lij2),
                                       _mm_mul_sd(prod, lij3));
            t1            = _mm_sub_sd(t1,
                                       _mm_mul_sd(onefourth,
                                                  _mm_add_sd(_mm_mul_sd(lij, rinv),
                                                             _mm_mul_sd(lij3, r))));
            t2            = _mm_mul_sd(onefourth,
                                       _mm_add_sd(_mm_mul_sd(uij, rinv),
                                                  _mm_mul_sd(uij3, r)));
            t2            = _mm_sub_sd(t2,
                                       _mm_add_sd(_mm_mul_sd(half, uij2),
                                                  _mm_mul_sd(prod, uij3)));
            t3            = _mm_mul_sd(_mm_mul_sd(onefourth, logterm),
                                       _mm_mul_sd(rinv, rinv));
            t3            = _mm_sub_sd(t3,
                                       _mm_mul_sd(_mm_mul_sd(diff2, oneeighth),
                                                  _mm_add_sd(one,
                                                             _mm_mul_sd(sk2_rinv, rinv))));
            t1            = _mm_mul_sd(rinv,
                                       _mm_add_sd(_mm_mul_sd(dlij, t1),
                                                  _mm_add_sd(t2, t3)));

            dadx2         = _mm_and_pd(t1, obc_mask1);

            _mm_store_pd(dadx, dadx1);
            dadx += 2;
            _mm_store_pd(dadx, dadx2);
            dadx += 2;
        }
        gmx_mm_update_1pot_pd(sum_ai, work+ii);

    }

    /* Parallel summations */
    if (DOMAINDECOMP(cr))
    {
        dd_atom_sum_real(cr->dd, work);
    }

    if (gb_algorithm == egbHCT)
    {
        /* HCT */
        for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
        {
            if (born->use[i] != 0)
            {
                rr      = top->atomtypes.gb_radius[md->typeA[i]]-doffset;
                sum     = 1.0/rr - work[i];
                min_rad = rr + doffset;
                rad     = 1.0/sum;

                born->bRad[i]   = rad > min_rad ? rad : min_rad;
                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);
            }
        }

        /* Extra communication required for DD */
        if (DOMAINDECOMP(cr))
        {
            dd_atom_spread_real(cr->dd, born->bRad);
            dd_atom_spread_real(cr->dd, fr->invsqrta);
        }
    }
    else
    {
        /* OBC */
        for (i = 0; i < fr->natoms_force; i++) /* PELA born->nr */
        {
            if (born->use[i] != 0)
            {
                rr      = top->atomtypes.gb_radius[md->typeA[i]];
                rr_inv2 = 1.0/rr;
                rr      = rr-doffset;
                rr_inv  = 1.0/rr;
                sum     = rr * work[i];
                sum2    = sum  * sum;
                sum3    = sum2 * sum;

                tsum          = tanh(born->obc_alpha*sum-born->obc_beta*sum2+born->obc_gamma*sum3);
                born->bRad[i] = rr_inv - tsum*rr_inv2;
                born->bRad[i] = 1.0 / born->bRad[i];

                fr->invsqrta[i] = gmx_invsqrt(born->bRad[i]);

                tchain         = rr * (born->obc_alpha-2*born->obc_beta*sum+3*born->obc_gamma*sum2);
                born->drobc[i] = (1.0-tsum*tsum)*tchain*rr_inv2;
            }
        }
        /* Extra (local) communication required for DD */
        if (DOMAINDECOMP(cr))
        {
            dd_atom_spread_real(cr->dd, born->bRad);
            dd_atom_spread_real(cr->dd, fr->invsqrta);
            dd_atom_spread_real(cr->dd, born->drobc);
        }
    }



    return 0;
}


int
calc_gb_chainrule_sse2_double(int natoms, t_nblist *nl, double *dadx, double *dvda,
                              double *x, double *f, double *fshift, double *shiftvec,
                              int gb_algorithm, gmx_genborn_t *born, t_mdatoms *md)
{
    int           i, k, n, ii, jnr, ii3, is3, nj0, nj1, n0, n1;
    int           jnrA, jnrB;
    int           j3A, j3B;
    int        *  jjnr;

    double        rbi, shX, shY, shZ;
    double       *rb;

    __m128d       ix, iy, iz;
    __m128d       jx, jy, jz;
    __m128d       fix, fiy, fiz;
    __m128d       dx, dy, dz;
    __m128d       tx, ty, tz;

    __m128d       rbai, rbaj, f_gb, f_gb_ai;
    __m128d       xmm1, xmm2, xmm3;

    const __m128d two = _mm_set1_pd(2.0);

    rb     = born->work;

    jjnr   = nl->jjnr;

    /* Loop to get the proper form for the Born radius term, sse style */
    n0 = 0;
    n1 = natoms;

    if (gb_algorithm == egbSTILL)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = (2 * rbi * rbi * dvda[i])/ONE_4PI_EPS0;
        }
    }
    else if (gb_algorithm == egbHCT)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * dvda[i];
        }
    }
    else if (gb_algorithm == egbOBC)
    {
        for (i = n0; i < n1; i++)
        {
            rbi   = born->bRad[i];
            rb[i] = rbi * rbi * born->drobc[i] * dvda[i];
        }
    }

    jz = _mm_setzero_pd();

    n = j3A = j3B = 0;

    for (i = 0; i < nl->nri; i++)
    {
        ii     = nl->iinr[i];
        ii3    = ii*3;
        is3    = 3*nl->shift[i];
        shX    = shiftvec[is3];
        shY    = shiftvec[is3+1];
        shZ    = shiftvec[is3+2];
        nj0    = nl->jindex[i];
        nj1    = nl->jindex[i+1];

        ix     = _mm_set1_pd(shX+x[ii3+0]);
        iy     = _mm_set1_pd(shY+x[ii3+1]);
        iz     = _mm_set1_pd(shZ+x[ii3+2]);

        rbai   = _mm_load1_pd(rb+ii);
        fix    = _mm_setzero_pd();
        fiy    = _mm_setzero_pd();
        fiz    = _mm_setzero_pd();


        for (k = nj0; k < nj1-1; k += 2)
        {
            jnrA        = jjnr[k];
            jnrB        = jjnr[k+1];

            j3A         = 3*jnrA;
            j3B         = 3*jnrB;

            GMX_MM_LOAD_1RVEC_2POINTERS_PD(x+j3A, x+j3B, jx, jy, jz);

            dx          = _mm_sub_pd(ix, jx);
            dy          = _mm_sub_pd(iy, jy);
            dz          = _mm_sub_pd(iz, jz);

            GMX_MM_LOAD_2VALUES_PD(rb+jnrA, rb+jnrB, rbaj);

            /* load chain rule terms for j1-4 */
            f_gb        = _mm_load_pd(dadx);
            dadx       += 2;
            f_gb_ai     = _mm_load_pd(dadx);
            dadx       += 2;

            /* calculate scalar force */
            f_gb    = _mm_mul_pd(f_gb, rbai);
            f_gb_ai = _mm_mul_pd(f_gb_ai, rbaj);
            f_gb    = _mm_add_pd(f_gb, f_gb_ai);

            tx     = _mm_mul_pd(f_gb, dx);
            ty     = _mm_mul_pd(f_gb, dy);
            tz     = _mm_mul_pd(f_gb, dz);

            fix    = _mm_add_pd(fix, tx);
            fiy    = _mm_add_pd(fiy, ty);
            fiz    = _mm_add_pd(fiz, tz);

            GMX_MM_DECREMENT_1RVEC_2POINTERS_PD(f+j3A, f+j3B, tx, ty, tz);
        }

        /*deal with odd elements */
        if (k < nj1)
        {
            jnrA        = jjnr[k];
            j3A         = 3*jnrA;

            GMX_MM_LOAD_1RVEC_1POINTER_PD(x+j3A, jx, jy, jz);

            dx          = _mm_sub_sd(ix, jx);
            dy          = _mm_sub_sd(iy, jy);
            dz          = _mm_sub_sd(iz, jz);

            GMX_MM_LOAD_1VALUE_PD(rb+jnrA, rbaj);

            /* load chain rule terms */
            f_gb        = _mm_load_pd(dadx);
            dadx       += 2;
            f_gb_ai     = _mm_load_pd(dadx);
            dadx       += 2;

            /* calculate scalar force */
            f_gb    = _mm_mul_sd(f_gb, rbai);
            f_gb_ai = _mm_mul_sd(f_gb_ai, rbaj);
            f_gb    = _mm_add_sd(f_gb, f_gb_ai);

            tx     = _mm_mul_sd(f_gb, dx);
            ty     = _mm_mul_sd(f_gb, dy);
            tz     = _mm_mul_sd(f_gb, dz);

            fix    = _mm_add_sd(fix, tx);
            fiy    = _mm_add_sd(fiy, ty);
            fiz    = _mm_add_sd(fiz, tz);

            GMX_MM_DECREMENT_1RVEC_1POINTER_PD(f+j3A, tx, ty, tz);
        }

        /* fix/fiy/fiz now contain four partial force terms, that all should be
         * added to the i particle forces and shift forces.
         */
        gmx_mm_update_iforce_1atom_pd(&fix, &fiy, &fiz, f+ii3, fshift+is3);
    }

    return 0;
}

#else
/* keep compiler happy */
int genborn_sse2_dummy;

#endif /* SSE2 intrinsics available */
