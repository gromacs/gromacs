/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/math/invertmatrix.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/pbc-simd.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

using namespace gmx; // TODO: Remove when this file is moved into gmx namespace

typedef struct
{
    real   mO;
    real   mH;
    real   wh;
    real   dOH;
    real   dHH;
    real   ra;
    real   rb;
    real   rc;
    real   irc2;
    /* For projection */
    real   imO;
    real   imH;
    real   invdOH;
    real   invdHH;
    matrix invmat;
} settleparam_t;

typedef struct gmx_settledata
{
    settleparam_t massw;    /* Parameters for SETTLE for coordinates */
    settleparam_t mass1;    /* Parameters with all masses 1, for forces */

    int           nsettle;  /* The number of settles on our rank */
    int          *ow1;      /* Index to OW1 atoms, size nsettle + SIMD padding */
    int          *hw2;      /* Index to HW2 atoms, size nsettle + SIMD padding */
    int          *hw3;      /* Index to HW3 atoms, size nsettle + SIMD padding */
    real         *virfac;   /* Virial factor 0 or 1, size nsettle + SIMD pad. */
    int           nalloc;   /* Allocation size of ow1, hw2, hw3, virfac */

    bool          bUseSimd; /* Use SIMD intrinsics code, if possible */
} t_gmx_settledata;


static void init_proj_matrix(real invmO, real invmH, real dOH, real dHH,
                             matrix inverseCouplingMatrix)
{
    /* We normalize the inverse masses with invmO for the matrix inversion.
     * so we can keep using masses of almost zero for frozen particles,
     * without running out of the float range in invertMatrix.
     */
    double invmORelative = 1.0;
    double invmHRelative = invmH/static_cast<double>(invmO);
    double distanceRatio = dHH/static_cast<double>(dOH);

    /* Construct the constraint coupling matrix */
    matrix mat;
    mat[0][0] = invmORelative + invmHRelative;
    mat[0][1] = invmORelative*(1.0 - 0.5*gmx::square(distanceRatio));
    mat[0][2] = invmHRelative*0.5*distanceRatio;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = invmHRelative + invmHRelative;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    invertMatrix(mat, inverseCouplingMatrix);

    msmul(inverseCouplingMatrix, 1/invmO, inverseCouplingMatrix);
}

static void settleparam_init(settleparam_t *p,
                             real mO, real mH, real invmO, real invmH,
                             real dOH, real dHH)
{
    /* We calculate parameters in double precision to minimize errors.
     * The velocity correction applied during SETTLE coordinate constraining
     * introduces a systematic error of approximately 1 bit per atom,
     * depending on what the compiler does with the code.
     */
    double wohh;

    p->mO     = mO;
    p->mH     = mH;
    wohh      = mO + 2.0*mH;
    p->wh     = mH/wohh;
    p->dOH    = dOH;
    p->dHH    = dHH;
    double rc = dHH/2.0;
    double ra = 2.0*mH*std::sqrt(dOH*dOH - rc*rc)/wohh;
    p->rb     = std::sqrt(dOH*dOH - rc*rc) - ra;
    p->rc     = rc;
    p->ra     = ra;
    p->irc2   = 1.0/dHH;

    /* For projection: inverse masses and coupling matrix inversion */
    p->imO    = invmO;
    p->imH    = invmH;

    p->invdOH = 1.0/dOH;
    p->invdHH = 1.0/dHH;

    init_proj_matrix(invmO, invmH, dOH, dHH, p->invmat);

    if (debug)
    {
        fprintf(debug, "wh =%g, rc = %g, ra = %g\n",
                p->wh, p->rc, p->ra);
        fprintf(debug, "rb = %g, irc2 = %g, dHH = %g, dOH = %g\n",
                p->rb, p->irc2, p->dHH, p->dOH);
    }
}

gmx_settledata_t settle_init(const gmx_mtop_t *mtop)
{
    /* Check that we have only one settle type */
    int                   settle_type = -1;
    gmx_mtop_ilistloop_t  iloop       = gmx_mtop_ilistloop_init(mtop);
    t_ilist              *ilist;
    int                   nmol;
    const int             nral1       = 1 + NRAL(F_SETTLE);
    while (gmx_mtop_ilistloop_next(iloop, &ilist, &nmol))
    {
        for (int i = 0; i < ilist[F_SETTLE].nr; i += nral1)
        {
            if (settle_type == -1)
            {
                settle_type = ilist[F_SETTLE].iatoms[i];
            }
            else if (ilist[F_SETTLE].iatoms[i] != settle_type)
            {
                gmx_fatal(FARGS,
                          "The [molecules] section of your topology specifies more than one block of\n"
                          "a [moleculetype] with a [settles] block. Only one such is allowed.\n"
                          "If you are trying to partition your solvent into different *groups*\n"
                          "(e.g. for freezing, T-coupling, etc.), you are using the wrong approach. Index\n"
                          "files specify groups. Otherwise, you may wish to change the least-used\n"
                          "block of molecules with SETTLE constraints into 3 normal constraints.");
            }
        }
    }
    GMX_RELEASE_ASSERT(settle_type >= 0, "settle_init called without settles");

    gmx_settledata_t settled;

    snew(settled, 1);

    /* We will not initialize the normal SETTLE parameters here yet,
     * since the atom (inv)masses can depend on the integrator and
     * free-energy perturbation. We set mO=-1 to trigger later initialization.
     */
    settled->massw.mO = -1;

    real dOH = mtop->ffparams.iparams[settle_type].settle.doh;
    real dHH = mtop->ffparams.iparams[settle_type].settle.dhh;
    settleparam_init(&settled->mass1, 1.0, 1.0, 1.0, 1.0, dOH, dHH);

    settled->ow1    = nullptr;
    settled->hw2    = nullptr;
    settled->hw3    = nullptr;
    settled->virfac = nullptr;
    settled->nalloc = 0;

    /* Without SIMD configured, this bool is not used */
    settled->bUseSimd = (getenv("GMX_DISABLE_SIMD_KERNELS") == nullptr);

    return settled;
}

void settle_free(gmx_settledata_t settled)
{
    sfree_aligned(settled->ow1);
    sfree_aligned(settled->hw2);
    sfree_aligned(settled->hw3);
    sfree_aligned(settled->virfac);
    sfree(settled);
}

void settle_set_constraints(gmx_settledata_t  settled,
                            const t_ilist    *il_settle,
                            const t_mdatoms  *mdatoms)
{
#if GMX_SIMD_HAVE_REAL
    const int pack_size = GMX_SIMD_REAL_WIDTH;
#else
    const int pack_size = 1;
#endif

    const int nral1     = 1 + NRAL(F_SETTLE);
    int       nsettle   = il_settle->nr/nral1;
    settled->nsettle = nsettle;

    if (nsettle > 0)
    {
        const t_iatom *iatoms = il_settle->iatoms;

        /* Here we initialize the normal SETTLE parameters */
        if (settled->massw.mO < 0)
        {
            int firstO = iatoms[1];
            int firstH = iatoms[2];
            settleparam_init(&settled->massw,
                             mdatoms->massT[firstO],
                             mdatoms->massT[firstH],
                             mdatoms->invmass[firstO],
                             mdatoms->invmass[firstH],
                             settled->mass1.dOH,
                             settled->mass1.dHH);
        }

        if (nsettle + pack_size > settled->nalloc)
        {
            settled->nalloc = over_alloc_dd(nsettle + pack_size);
            sfree_aligned(settled->ow1);
            sfree_aligned(settled->hw2);
            sfree_aligned(settled->hw3);
            sfree_aligned(settled->virfac);
            snew_aligned(settled->ow1, settled->nalloc, 64);
            snew_aligned(settled->hw2, settled->nalloc, 64);
            snew_aligned(settled->hw3, settled->nalloc, 64);
            snew_aligned(settled->virfac, settled->nalloc, 64);
        }

        for (int i = 0; i < nsettle; i++)
        {
            settled->ow1[i]    = iatoms[i*nral1 + 1];
            settled->hw2[i]    = iatoms[i*nral1 + 2];
            settled->hw3[i]    = iatoms[i*nral1 + 3];
            /* We should avoid double counting of virial contributions for
             * SETTLEs that appear in multiple DD domains, so we only count
             * the contribution on the home range of the oxygen atom.
             */
            settled->virfac[i] = (iatoms[i*nral1 + 1] < mdatoms->homenr ? 1 : 0);
        }

        /* Pack the index array to the full SIMD width with copies from
         * the last normal entry, but with no virial contribution.
         */
        int end_packed = ((nsettle + pack_size - 1)/pack_size)*pack_size;
        for (int i = nsettle; i < end_packed; i++)
        {
            settled->ow1[i]    = settled->ow1[nsettle - 1];
            settled->hw2[i]    = settled->hw2[nsettle - 1];
            settled->hw3[i]    = settled->hw3[nsettle - 1];
            settled->virfac[i] = 0;
        }
    }
}

void settle_proj(gmx_settledata_t settled, int econq,
                 int nsettle, t_iatom iatoms[],
                 const t_pbc *pbc,
                 rvec x[],
                 rvec *der, rvec *derp,
                 int calcvir_atom_end, tensor vir_r_m_dder)
{
    /* Settle for projection out constraint components
     * of derivatives of the coordinates.
     * Berk Hess 2008-1-10
     */

    settleparam_t *p;
    real           imO, imH, dOH, dHH, invdOH, invdHH;
    matrix         invmat;
    int            i, m, m2, ow1, hw2, hw3;
    rvec           roh2, roh3, rhh, dc, fc;

    calcvir_atom_end *= DIM;

    if (econq == econqForce)
    {
        p = &settled->mass1;
    }
    else
    {
        p = &settled->massw;
    }
    imO    = p->imO;
    imH    = p->imH;
    copy_mat(p->invmat, invmat);
    dOH    = p->dOH;
    dHH    = p->dHH;
    invdOH = p->invdOH;
    invdHH = p->invdHH;

#ifdef PRAGMAS
#pragma ivdep
#endif

    const int nral1 = 1 + NRAL(F_SETTLE);

    for (i = 0; i < nsettle; i++)
    {
        ow1 = iatoms[i*nral1 + 1];
        hw2 = iatoms[i*nral1 + 2];
        hw3 = iatoms[i*nral1 + 3];

        if (pbc == nullptr)
        {
            rvec_sub(x[ow1], x[hw2], roh2);
            rvec_sub(x[ow1], x[hw3], roh3);
            rvec_sub(x[hw2], x[hw3], rhh);
        }
        else
        {
            pbc_dx_aiuc(pbc, x[ow1], x[hw2], roh2);
            pbc_dx_aiuc(pbc, x[ow1], x[hw3], roh3);
            pbc_dx_aiuc(pbc, x[hw2], x[hw3], rhh);
        }
        svmul(invdOH, roh2, roh2);
        svmul(invdOH, roh3, roh3);
        svmul(invdHH, rhh, rhh);
        /* 18 flops */

        /* Determine the projections of der on the bonds */
        clear_rvec(dc);
        for (m = 0; m < DIM; m++)
        {
            dc[0] += (der[ow1][m] - der[hw2][m])*roh2[m];
            dc[1] += (der[ow1][m] - der[hw3][m])*roh3[m];
            dc[2] += (der[hw2][m] - der[hw3][m])*rhh [m];
        }
        /* 27 flops */

        /* Determine the correction for the three bonds */
        mvmul(invmat, dc, fc);
        /* 15 flops */

        /* Subtract the corrections from derp */
        for (m = 0; m < DIM; m++)
        {
            derp[ow1][m] -= imO*( fc[0]*roh2[m] + fc[1]*roh3[m]);
            derp[hw2][m] -= imH*(-fc[0]*roh2[m] + fc[2]*rhh [m]);
            derp[hw3][m] -= imH*(-fc[1]*roh3[m] - fc[2]*rhh [m]);
        }

        /* 45 flops */

        if (ow1 < calcvir_atom_end)
        {
            /* Determining r \dot m der is easy,
             * since fc contains the mass weighted corrections for der.
             */

            for (m = 0; m < DIM; m++)
            {
                for (m2 = 0; m2 < DIM; m2++)
                {
                    vir_r_m_dder[m][m2] +=
                        dOH*roh2[m]*roh2[m2]*fc[0] +
                        dOH*roh3[m]*roh3[m2]*fc[1] +
                        dHH*rhh [m]*rhh [m2]*fc[2];
                }
            }
        }
    }
}


/* The actual settle code, templated for real/SimdReal and for optimization */
template<typename T, typename TypeBool, int packSize,
         typename TypePbc,
         bool bCorrectVelocity,
         bool bCalcVirial>
static void settleTemplate(const gmx_settledata_t settled,
                           int settleStart, int settleEnd,
                           const TypePbc pbc,
                           const real *x, real *xprime,
                           real invdt, real * gmx_restrict v,
                           tensor vir_r_m_dr,
                           bool *bErrorHasOccurred)
{
    /* ******************************************************************* */
    /*                                                                  ** */
    /*    Original code by Shuichi Miyamoto, last update Oct. 1, 1992   ** */
    /*                                                                  ** */
    /*    Algorithm changes by Berk Hess:                               ** */
    /*    2004-07-15 Convert COM to double precision to avoid drift     ** */
    /*    2006-10-16 Changed velocity update to use differences         ** */
    /*    2012-09-24 Use oxygen as reference instead of COM             ** */
    /*    2016-02    Complete rewrite of the code for SIMD              ** */
    /*                                                                  ** */
    /*    Reference for the SETTLE algorithm                            ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992).    ** */
    /*                                                                  ** */
    /* ******************************************************************* */

    assert(settleStart % packSize == 0);
    assert(settleEnd   % packSize == 0);

    TypeBool       bError = TypeBool(false);

    settleparam_t *p    = &settled->massw;
    T              wh   = T(p->wh);
    T              rc   = T(p->rc);
    T              ra   = T(p->ra);
    T              rb   = T(p->rb);
    T              irc2 = T(p->irc2);
    T              mO   = T(p->mO);
    T              mH   = T(p->mH);

    T              almost_zero = T(1e-12);

    T              sum_r_m_dr[DIM][DIM];

    if (bCalcVirial)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            for (int d = 0; d < DIM; d++)
            {
                sum_r_m_dr[d2][d] = T(0);
            }
        }
    }

    for (int i = settleStart; i < settleEnd; i += packSize)
    {
        /* Here we pad up to packSize with copies from the last valid entry.
         * This gives correct results, since we store (not increment) all
         * output, so we store the same output multiple times.
         */
        const int *ow1 = settled->ow1 + i;
        const int *hw2 = settled->hw2 + i;
        const int *hw3 = settled->hw3 + i;

        T          x_ow1[DIM], x_hw2[DIM], x_hw3[DIM];

        gatherLoadUTranspose<3>(x, ow1, &x_ow1[XX], &x_ow1[YY], &x_ow1[ZZ]);
        gatherLoadUTranspose<3>(x, hw2, &x_hw2[XX], &x_hw2[YY], &x_hw2[ZZ]);
        gatherLoadUTranspose<3>(x, hw3, &x_hw3[XX], &x_hw3[YY], &x_hw3[ZZ]);

        T xprime_ow1[DIM], xprime_hw2[DIM], xprime_hw3[DIM];

        gatherLoadUTranspose<3>(xprime, ow1, &xprime_ow1[XX], &xprime_ow1[YY], &xprime_ow1[ZZ]);
        gatherLoadUTranspose<3>(xprime, hw2, &xprime_hw2[XX], &xprime_hw2[YY], &xprime_hw2[ZZ]);
        gatherLoadUTranspose<3>(xprime, hw3, &xprime_hw3[XX], &xprime_hw3[YY], &xprime_hw3[ZZ]);

        T dist21[DIM], dist31[DIM];
        T doh2[DIM], doh3[DIM];
        T sh_hw2[DIM], sh_hw3[DIM];

        pbc_dx_aiuc(pbc, x_hw2, x_ow1, dist21);

        pbc_dx_aiuc(pbc, x_hw3, x_ow1, dist31);

        /* Tedious way of doing pbc */
        pbc_dx_aiuc(pbc, xprime_hw2, xprime_ow1, doh2);
        for (int d = 0; d < DIM; d++)
        {
            sh_hw2[d]     = xprime_hw2[d] - (xprime_ow1[d] + doh2[d]);
            xprime_hw2[d] = xprime_hw2[d] - sh_hw2[d];
        }
        pbc_dx_aiuc(pbc, xprime_hw3, xprime_ow1, doh3);
        for (int d = 0; d < DIM; d++)
        {
            sh_hw3[d]     = xprime_hw3[d] - (xprime_ow1[d] + doh3[d]);
            xprime_hw3[d] = xprime_hw3[d] - sh_hw3[d];
        }

        /* Not calculating the center of mass using the oxygen position
         * and the O-H distances, as done below, will make SETTLE
         * the largest source of energy drift for simulations of water,
         * as then the oxygen coordinate is multiplied by 0.89 at every step,
         * which can then transfer a systematic rounding to the oxygen velocity.
         */
        T a1[DIM], com[DIM];
        for (int d = 0; d < DIM; d++)
        {
            a1[d]  = -(doh2[d] + doh3[d]) * wh;
            com[d] = xprime_ow1[d] - a1[d];
        }
        T b1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            b1[d]  = xprime_hw2[d] - com[d];
        }
        T c1[DIM];
        for (int d = 0; d < DIM; d++)
        {
            c1[d]  = xprime_hw3[d] - com[d];
        }
        /* 15 flops */

        T xakszd = dist21[YY] * dist31[ZZ] - dist21[ZZ] * dist31[YY];
        T yakszd = dist21[ZZ] * dist31[XX] - dist21[XX] * dist31[ZZ];
        T zakszd = dist21[XX] * dist31[YY] - dist21[YY] * dist31[XX];
        T xaksxd = a1[YY] * zakszd - a1[ZZ] * yakszd;
        T yaksxd = a1[ZZ] * xakszd - a1[XX] * zakszd;
        T zaksxd = a1[XX] * yakszd - a1[YY] * xakszd;
        T xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        T yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        T zaksyd = xakszd * yaksxd - yakszd * xaksxd;
        /* 27 flops */

        T axlng = gmx::invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        T aylng = gmx::invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        T azlng = gmx::invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        T trns1[DIM], trns2[DIM], trns3[DIM];

        trns1[XX] = xaksxd * axlng;
        trns2[XX] = yaksxd * axlng;
        trns3[XX] = zaksxd * axlng;
        trns1[YY] = xaksyd * aylng;
        trns2[YY] = yaksyd * aylng;
        trns3[YY] = zaksyd * aylng;
        trns1[ZZ] = xakszd * azlng;
        trns2[ZZ] = yakszd * azlng;
        trns3[ZZ] = zakszd * azlng;
        /* 24 flops */

        T b0d[2], c0d[2];

        for (int d = 0; d < 2; d++)
        {
            b0d[d] = trns1[d] * dist21[XX] + trns2[d] * dist21[YY] + trns3[d] * dist21[ZZ];
            c0d[d] = trns1[d] * dist31[XX] + trns2[d] * dist31[YY] + trns3[d] * dist31[ZZ];
        }

        T a1d_z, b1d[DIM], c1d[DIM];

        a1d_z = trns1[ZZ] * a1[XX] + trns2[ZZ] * a1[YY] + trns3[ZZ] * a1[ZZ];
        for (int d = 0; d < DIM; d++)
        {
            b1d[d] = trns1[d] * b1[XX] + trns2[d] * b1[YY] + trns3[d] * b1[ZZ];
            c1d[d] = trns1[d] * c1[XX] + trns2[d] * c1[YY] + trns3[d] * c1[ZZ];
        }
        /* 65 flops */

        T tmp, tmp2;

        T sinphi = a1d_z * gmx::invsqrt(ra*ra);
        tmp2     = 1.0 - sinphi * sinphi;

        /* If tmp2 gets close to or beyond zero we have severly distorted
         * water molecules and we should terminate the simulation.
         * Below we take the max with almost_zero to continue the loop.
         */
        bError   = bError || (tmp2 <= almost_zero);

        tmp2     = max(tmp2, almost_zero);
        tmp      = gmx::invsqrt(tmp2);
        T cosphi = tmp2*tmp;
        T sinpsi = (b1d[ZZ] - c1d[ZZ]) * irc2 * tmp;
        tmp2     = 1.0 - sinpsi * sinpsi;

        T cospsi = tmp2*gmx::invsqrt(tmp2);
        /* 46 flops */

        T a2d_y  =  ra * cosphi;
        T b2d_x  = -rc * cospsi;
        T t1     = -rb * cosphi;
        T t2     =  rc * sinpsi * sinphi;
        T b2d_y  =  t1 - t2;
        T c2d_y  =  t1 + t2;
        /* 7 flops */

        /*     --- Step3  al,be,ga            --- */
        T alpha  = b2d_x * (b0d[XX] - c0d[XX]) + b0d[YY] * b2d_y + c0d[YY] * c2d_y;
        T beta   = b2d_x * (c0d[YY] - b0d[YY]) + b0d[XX] * b2d_y + c0d[XX] * c2d_y;
        T gamma  = b0d[XX] * b1d[YY] - b1d[XX] * b0d[YY] + c0d[XX] * c1d[YY] - c1d[XX] * c0d[YY];
        T al2be2 = alpha * alpha + beta * beta;
        tmp2     = (al2be2 - gamma * gamma);
        T sinthe = (alpha * gamma - beta * tmp2*gmx::invsqrt(tmp2)) * gmx::invsqrt(al2be2*al2be2);
        /* 47 flops */

        /*  --- Step4  A3' --- */
        tmp2     = 1.0 - sinthe * sinthe;
        T costhe = tmp2*gmx::invsqrt(tmp2);

        T a3d[DIM], b3d[DIM], c3d[DIM];

        a3d[XX]  = -a2d_y * sinthe;
        a3d[YY]  = a2d_y * costhe;
        a3d[ZZ]  = a1d_z;
        b3d[XX]  = b2d_x * costhe - b2d_y * sinthe;
        b3d[YY]  = b2d_x * sinthe + b2d_y * costhe;
        b3d[ZZ]  = b1d[ZZ];
        c3d[XX]  = -b2d_x * costhe - c2d_y * sinthe;
        c3d[YY]  = -b2d_x * sinthe + c2d_y * costhe;
        c3d[ZZ]  = c1d[ZZ];
        /* 26 flops */

        /*    --- Step5  A3 --- */
        T a3[DIM], b3[DIM], c3[DIM];

        a3[XX] = trns1[XX]*a3d[XX] + trns1[YY]*a3d[YY] + trns1[ZZ]*a3d[ZZ];
        a3[YY] = trns2[XX]*a3d[XX] + trns2[YY]*a3d[YY] + trns2[ZZ]*a3d[ZZ];
        a3[ZZ] = trns3[XX]*a3d[XX] + trns3[YY]*a3d[YY] + trns3[ZZ]*a3d[ZZ];
        b3[XX] = trns1[XX]*b3d[XX] + trns1[YY]*b3d[YY] + trns1[ZZ]*b3d[ZZ];
        b3[YY] = trns2[XX]*b3d[XX] + trns2[YY]*b3d[YY] + trns2[ZZ]*b3d[ZZ];
        b3[ZZ] = trns3[XX]*b3d[XX] + trns3[YY]*b3d[YY] + trns3[ZZ]*b3d[ZZ];
        c3[XX] = trns1[XX]*c3d[XX] + trns1[YY]*c3d[YY] + trns1[ZZ]*c3d[ZZ];
        c3[YY] = trns2[XX]*c3d[XX] + trns2[YY]*c3d[YY] + trns2[ZZ]*c3d[ZZ];
        c3[ZZ] = trns3[XX]*c3d[XX] + trns3[YY]*c3d[YY] + trns3[ZZ]*c3d[ZZ];
        /* 45 flops */

        /* Compute and store the corrected new coordinate */
        for (int d = 0; d < DIM; d++)
        {
            xprime_ow1[d] = com[d] + a3[d];
        }
        for (int d = 0; d < DIM; d++)
        {
            xprime_hw2[d] = com[d] + b3[d] + sh_hw2[d];;
        }
        for (int d = 0; d < DIM; d++)
        {
            xprime_hw3[d] = com[d] + c3[d] + sh_hw3[d];
        }
        /* 9 flops + 6 pbc flops */

        transposeScatterStoreU<3>(xprime, ow1, xprime_ow1[XX], xprime_ow1[YY], xprime_ow1[ZZ]);
        transposeScatterStoreU<3>(xprime, hw2, xprime_hw2[XX], xprime_hw2[YY], xprime_hw2[ZZ]);
        transposeScatterStoreU<3>(xprime, hw3, xprime_hw3[XX], xprime_hw3[YY], xprime_hw3[ZZ]);

        // cppcheck-suppress duplicateExpression
        if (bCorrectVelocity || bCalcVirial)
        {
            T da[DIM], db[DIM], dc[DIM];
            for (int d = 0; d < DIM; d++)
            {
                da[d] = a3[d] - a1[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                db[d] = b3[d] - b1[d];
            }
            for (int d = 0; d < DIM; d++)
            {
                dc[d] = c3[d] - c1[d];
            }
            /* 9 flops */

            if (bCorrectVelocity)
            {
                T v_ow1[DIM], v_hw2[DIM], v_hw3[DIM];

                gatherLoadUTranspose<3>(v, ow1, &v_ow1[XX], &v_ow1[YY], &v_ow1[ZZ]);
                gatherLoadUTranspose<3>(v, hw2, &v_hw2[XX], &v_hw2[YY], &v_hw2[ZZ]);
                gatherLoadUTranspose<3>(v, hw3, &v_hw3[XX], &v_hw3[YY], &v_hw3[ZZ]);

                /* Add the position correction divided by dt to the velocity */
                for (int d = 0; d < DIM; d++)
                {
                    v_ow1[d] = gmx::fma(da[d], invdt, v_ow1[d]);
                }
                for (int d = 0; d < DIM; d++)
                {
                    v_hw2[d] = gmx::fma(db[d], invdt, v_hw2[d]);
                }
                for (int d = 0; d < DIM; d++)
                {
                    v_hw3[d] = gmx::fma(dc[d], invdt, v_hw3[d]);
                }
                /* 3*6 flops */

                transposeScatterStoreU<3>(v, ow1, v_ow1[XX], v_ow1[YY], v_ow1[ZZ]);
                transposeScatterStoreU<3>(v, hw2, v_hw2[XX], v_hw2[YY], v_hw2[ZZ]);
                transposeScatterStoreU<3>(v, hw3, v_hw3[XX], v_hw3[YY], v_hw3[ZZ]);
            }

            if (bCalcVirial)
            {
                /* Filter out the non-local settles */
                T filter = load<T>(settled->virfac + i);
                T mOf    = filter*mO;
                T mHf    = filter*mH;

                T mdo[DIM], mdb[DIM], mdc[DIM];

                for (int d = 0; d < DIM; d++)
                {
                    mdb[d] = mHf*db[d];
                    mdc[d] = mHf*dc[d];
                    mdo[d] = mOf*da[d] + mdb[d] + mdc[d];
                }

                for (int d2 = 0; d2 < DIM; d2++)
                {
                    for (int d = 0; d < DIM; d++)
                    {
                        sum_r_m_dr[d2][d] = sum_r_m_dr[d2][d] -
                            (x_ow1[d2]*mdo[d] +
                             dist21[d2]*mdb[d] +
                             dist31[d2]*mdc[d]);
                    }
                }
                /* 71 flops */
            }
        }
    }

    if (bCalcVirial)
    {
        for (int d2 = 0; d2 < DIM; d2++)
        {
            for (int d = 0; d < DIM; d++)
            {
                vir_r_m_dr[d2][d] += reduce(sum_r_m_dr[d2][d]);
            }
        }
    }

    *bErrorHasOccurred = anyTrue(bError);
}

/* Wrapper template function that divides the settles over threads
 * and instantiates the core template with instantiated booleans.
 */
template<typename T, typename TypeBool, int packSize, typename TypePbc>
static void settleTemplateWrapper(gmx_settledata_t settled,
                                  int nthread, int thread,
                                  TypePbc pbc,
                                  const real x[], real xprime[],
                                  real invdt, real *v,
                                  bool bCalcVirial, tensor vir_r_m_dr,
                                  bool *bErrorHasOccurred)
{
    /* We need to assign settles to threads in groups of pack_size */
    int numSettlePacks = (settled->nsettle + packSize - 1)/packSize;
    /* Round the end value up to give thread 0 more work */
    int settleStart    = ((numSettlePacks* thread      + nthread - 1)/nthread)*packSize;
    int settleEnd      = ((numSettlePacks*(thread + 1) + nthread - 1)/nthread)*packSize;

    if (v != nullptr)
    {
        if (!bCalcVirial)
        {
            settleTemplate<T, TypeBool, packSize,
                           TypePbc,
                           true,
                           false>
                (settled, settleStart, settleEnd,
                pbc,
                x, xprime,
                invdt, v,
                nullptr,
                bErrorHasOccurred);
        }
        else
        {
            settleTemplate<T, TypeBool, packSize,
                           TypePbc,
                           true,
                           true>
                (settled, settleStart, settleEnd,
                pbc,
                x, xprime,
                invdt, v,
                vir_r_m_dr,
                bErrorHasOccurred);
        }
    }
    else
    {
        if (!bCalcVirial)
        {
            settleTemplate<T, TypeBool, packSize,
                           TypePbc,
                           false,
                           false>
                (settled, settleStart, settleEnd,
                pbc,
                x, xprime,
                invdt, v,
                nullptr,
                bErrorHasOccurred);
        }
        else
        {
            settleTemplate<T, TypeBool, packSize,
                           TypePbc,
                           false,
                           true>
                (settled, settleStart, settleEnd,
                pbc,
                x, xprime,
                invdt, v,
                vir_r_m_dr,
                bErrorHasOccurred);
        }
    }
}

void csettle(gmx_settledata_t settled,
             int nthread, int thread,
             const t_pbc *pbc,
             const real x[], real xprime[],
             real invdt, real *v,
             bool bCalcVirial, tensor vir_r_m_dr,
             bool *bErrorHasOccurred)
{
#if GMX_SIMD_HAVE_REAL
    if (settled->bUseSimd)
    {
        /* Convert the pbc struct for SIMD */
        GMX_ALIGNED(real, GMX_SIMD_REAL_WIDTH) pbcSimd[9*GMX_SIMD_REAL_WIDTH];
        set_pbc_simd(pbc, pbcSimd);

        settleTemplateWrapper<SimdReal, SimdBool, GMX_SIMD_REAL_WIDTH,
                              const real *>(settled,
                                            nthread, thread,
                                            pbcSimd,
                                            x, xprime,
                                            invdt,
                                            v,
                                            bCalcVirial, vir_r_m_dr,
                                            bErrorHasOccurred);
    }
    else
#endif
    {
        /* This construct is needed because pbc_dx_aiuc doesn't accept pbc=NULL */
        t_pbc        pbcNo;
        const t_pbc *pbcNonNull;

        if (pbc != nullptr)
        {
            pbcNonNull = pbc;
        }
        else
        {
            set_pbc(&pbcNo, epbcNONE, nullptr);
            pbcNonNull = &pbcNo;
        }

        settleTemplateWrapper<real, bool, 1,
                              const t_pbc *>(settled,
                                             nthread, thread,
                                             pbcNonNull,
                                             x, xprime,
                                             invdt,
                                             v,
                                             bCalcVirial, vir_r_m_dr,
                                             bErrorHasOccurred);
    }
}
