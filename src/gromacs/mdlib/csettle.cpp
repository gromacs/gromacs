/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#include <stdio.h>

#include "gromacs/legacyheaders/constr.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

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
    settleparam_t massw;
    settleparam_t mass1;
} t_gmx_settledata;


static void init_proj_matrix(settleparam_t *p,
                             real invmO, real invmH, real dOH, real dHH)
{
    real   imOn, imHn;
    matrix mat;

    p->imO = invmO;
    p->imH = invmH;
    /* We normalize the inverse masses with imO for the matrix inversion.
     * so we can keep using masses of almost zero for frozen particles,
     * without running out of the float range in m_inv.
     */
    imOn = 1;
    imHn = p->imH/p->imO;

    /* Construct the constraint coupling matrix */
    mat[0][0] = imOn + imHn;
    mat[0][1] = imOn*(1 - 0.5*dHH*dHH/(dOH*dOH));
    mat[0][2] = imHn*0.5*dHH/dOH;
    mat[1][1] = mat[0][0];
    mat[1][2] = mat[0][2];
    mat[2][2] = imHn + imHn;
    mat[1][0] = mat[0][1];
    mat[2][0] = mat[0][2];
    mat[2][1] = mat[1][2];

    m_inv(mat, p->invmat);

    msmul(p->invmat, 1/p->imO, p->invmat);

    p->invdOH = 1/dOH;
    p->invdHH = 1/dHH;
}

static void settleparam_init(settleparam_t *p,
                             real mO, real mH, real invmO, real invmH,
                             real dOH, real dHH)
{
    double wohh;

    p->mO     = mO;
    p->mH     = mH;
    wohh      = mO + 2.0*mH;
    p->wh     = mH/wohh;
    p->dOH    = dOH;
    p->dHH    = dHH;
    p->rc     = dHH/2.0;
    p->ra     = 2.0*mH*sqrt(dOH*dOH - p->rc*p->rc)/wohh;
    p->rb     = sqrt(dOH*dOH - p->rc*p->rc) - p->ra;
    p->irc2   = 1.0/dHH;

    /* For projection: connection matrix inversion */
    init_proj_matrix(p, invmO, invmH, dOH, dHH);

    if (debug)
    {
        fprintf(debug, "wh =%g, rc = %g, ra = %g\n",
                p->wh, p->rc, p->ra);
        fprintf(debug, "rb = %g, irc2 = %g, dHH = %g, dOH = %g\n",
                p->rb, p->irc2, p->dHH, p->dOH);
    }
}

gmx_settledata_t settle_init(real mO, real mH, real invmO, real invmH,
                             real dOH, real dHH)
{
    gmx_settledata_t settled;

    snew(settled, 1);

    settleparam_init(&settled->massw, mO, mH, invmO, invmH, dOH, dHH);

    settleparam_init(&settled->mass1, 1.0, 1.0, 1.0, 1.0, dOH, dHH);

    return settled;
}

#ifdef DEBUG
static void check_cons(FILE *fp, char *title, real x[], int OW1, int HW2, int HW3)
{
    rvec dOH1, dOH2, dHH;
    int  m;

    for (m = 0; (m < DIM); m++)
    {
        dOH1[m] = x[OW1+m]-x[HW2+m];
        dOH2[m] = x[OW1+m]-x[HW3+m];
        dHH[m]  = x[HW2+m]-x[HW3+m];
    }
    fprintf(fp, "%10s, OW1=%3d, HW2=%3d, HW3=%3d,  dOH1: %8.3f, dOH2: %8.3f, dHH: %8.3f\n",
            title, OW1/DIM, HW2/DIM, HW3/DIM, norm(dOH1), norm(dOH2), norm(dHH));
}
#endif


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

    for (i = 0; i < nsettle; i++)
    {
        ow1 = iatoms[i*4+1];
        hw2 = iatoms[i*4+2];
        hw3 = iatoms[i*4+3];

        if (pbc == NULL)
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


void csettle(gmx_settledata_t settled,
             int nsettle, t_iatom iatoms[],
             const t_pbc *pbc,
             real b4[], real after[],
             real invdt, real *v, int CalcVirAtomEnd,
             tensor vir_r_m_dr,
             int *error)
{
    /* ***************************************************************** */
    /*                                                               ** */
    /*    Subroutine : setlep - reset positions of TIP3P waters      ** */
    /*    Author : Shuichi Miyamoto                                  ** */
    /*    Date of last update : Oct. 1, 1992                         ** */
    /*                                                               ** */
    /*    Reference for the SETTLE algorithm                         ** */
    /*           S. Miyamoto et al., J. Comp. Chem., 13, 952 (1992). ** */
    /*                                                               ** */
    /* ***************************************************************** */

    /* Initialized data */
    settleparam_t *p;
    real           wh, ra, rb, rc, irc2;
    real           mO, mH;

    /* Local variables */
    real     gama, beta, alpa, xcom, ycom, zcom, al2be2, tmp, tmp2;
    real     axlng, aylng, azlng, trns11, trns21, trns31, trns12, trns22,
             trns32, trns13, trns23, trns33, cosphi, costhe, sinphi, sinthe,
             cospsi, xaksxd, yaksxd, xakszd, yakszd, zakszd, zaksxd, xaksyd,
             xb0, yb0, zb0, xc0, yc0, zc0, xa1;
    real     ya1, za1, xb1, yb1;
    real     zb1, xc1, yc1, zc1, yaksyd, zaksyd, sinpsi, xa3, ya3, za3,
             xb3, yb3, zb3, xc3, yc3, zc3, xb0d, yb0d, xc0d, yc0d,
             za1d, xb1d, yb1d, zb1d, xc1d, yc1d, zc1d, ya2d, xb2d, yb2d, yc2d,
             xa3d, ya3d, za3d, xb3d, yb3d, zb3d, xc3d, yc3d, zc3d;
    real     t1, t2;
    real     dax, day, daz, dbx, dby, dbz, dcx, dcy, dcz;
    real     mdax, mday, mdaz, mdbx, mdby, mdbz, mdcx, mdcy, mdcz;

    gmx_bool bOK;

    int      i, ow1, hw2, hw3;

    rvec     dx, sh_hw2 = {0, 0, 0}, sh_hw3 = {0, 0, 0};
    rvec     doh2, doh3;
    int      is;

    *error = -1;

    CalcVirAtomEnd *= 3;

    p     = &settled->massw;
    wh    = p->wh;
    rc    = p->rc;
    ra    = p->ra;
    rb    = p->rb;
    irc2  = p->irc2;
    mO    = p->mO;
    mH    = p->mH;

#ifdef PRAGMAS
#pragma ivdep
#endif
    for (i = 0; i < nsettle; ++i)
    {
        bOK = TRUE;
        /*    --- Step1  A1' ---      */
        ow1 = iatoms[i*4+1] * 3;
        hw2 = iatoms[i*4+2] * 3;
        hw3 = iatoms[i*4+3] * 3;
        if (pbc == NULL)
        {
            xb0 = b4[hw2 + XX] - b4[ow1 + XX];
            yb0 = b4[hw2 + YY] - b4[ow1 + YY];
            zb0 = b4[hw2 + ZZ] - b4[ow1 + ZZ];
            xc0 = b4[hw3 + XX] - b4[ow1 + XX];
            yc0 = b4[hw3 + YY] - b4[ow1 + YY];
            zc0 = b4[hw3 + ZZ] - b4[ow1 + ZZ];
            /* 6 flops */

            rvec_sub(after+hw2, after+ow1, doh2);
            rvec_sub(after+hw3, after+ow1, doh3);
            /* 6 flops */
        }
        else
        {
            pbc_dx_aiuc(pbc, b4+hw2, b4+ow1, dx);
            xb0 = dx[XX];
            yb0 = dx[YY];
            zb0 = dx[ZZ];
            pbc_dx_aiuc(pbc, b4+hw3, b4+ow1, dx);
            xc0 = dx[XX];
            yc0 = dx[YY];
            zc0 = dx[ZZ];

            /* Tedious way of doing pbc */
            is = pbc_dx_aiuc(pbc, after+hw2, after+ow1, doh2);
            if (is == CENTRAL)
            {
                clear_rvec(sh_hw2);
            }
            else
            {
                sh_hw2[XX] = after[hw2 + XX] - (after[ow1 + XX] + doh2[XX]);
                sh_hw2[YY] = after[hw2 + YY] - (after[ow1 + YY] + doh2[YY]);
                sh_hw2[ZZ] = after[hw2 + ZZ] - (after[ow1 + ZZ] + doh2[ZZ]);
                rvec_dec(after+hw2, sh_hw2);
            }
            is = pbc_dx_aiuc(pbc, after+hw3, after+ow1, doh3);
            if (is == CENTRAL)
            {
                clear_rvec(sh_hw3);
            }
            else
            {
                sh_hw3[XX] = after[hw3 + XX] - (after[ow1 + XX] + doh3[XX]);
                sh_hw3[YY] = after[hw3 + YY] - (after[ow1 + YY] + doh3[YY]);
                sh_hw3[ZZ] = after[hw3 + ZZ] - (after[ow1 + ZZ] + doh3[ZZ]);
                rvec_dec(after+hw3, sh_hw3);
            }
        }

        /* Not calculating the center of mass using the oxygen position
         * and the O-H distances, as done below, will make SETTLE
         * the largest source of energy drift for simulations of water,
         * as then the oxygen coordinate is multiplied by 0.89 at every step,
         * which can then transfer a systematic rounding to the oxygen velocity.
         */
        xa1 = -(doh2[XX] + doh3[XX]) * wh;
        ya1 = -(doh2[YY] + doh3[YY]) * wh;
        za1 = -(doh2[ZZ] + doh3[ZZ]) * wh;

        xcom = after[ow1 + XX] - xa1;
        ycom = after[ow1 + YY] - ya1;
        zcom = after[ow1 + ZZ] - za1;

        xb1 = after[hw2 + XX] - xcom;
        yb1 = after[hw2 + YY] - ycom;
        zb1 = after[hw2 + ZZ] - zcom;
        xc1 = after[hw3 + XX] - xcom;
        yc1 = after[hw3 + YY] - ycom;
        zc1 = after[hw3 + ZZ] - zcom;
        /* 15 flops */

        xakszd = yb0 * zc0 - zb0 * yc0;
        yakszd = zb0 * xc0 - xb0 * zc0;
        zakszd = xb0 * yc0 - yb0 * xc0;
        xaksxd = ya1 * zakszd - za1 * yakszd;
        yaksxd = za1 * xakszd - xa1 * zakszd;
        zaksxd = xa1 * yakszd - ya1 * xakszd;
        xaksyd = yakszd * zaksxd - zakszd * yaksxd;
        yaksyd = zakszd * xaksxd - xakszd * zaksxd;
        zaksyd = xakszd * yaksxd - yakszd * xaksxd;
        /* 27 flops */

        axlng = gmx_invsqrt(xaksxd * xaksxd + yaksxd * yaksxd + zaksxd * zaksxd);
        aylng = gmx_invsqrt(xaksyd * xaksyd + yaksyd * yaksyd + zaksyd * zaksyd);
        azlng = gmx_invsqrt(xakszd * xakszd + yakszd * yakszd + zakszd * zakszd);

        trns11 = xaksxd * axlng;
        trns21 = yaksxd * axlng;
        trns31 = zaksxd * axlng;
        trns12 = xaksyd * aylng;
        trns22 = yaksyd * aylng;
        trns32 = zaksyd * aylng;
        trns13 = xakszd * azlng;
        trns23 = yakszd * azlng;
        trns33 = zakszd * azlng;
        /* 24 flops */

        xb0d = trns11 * xb0 + trns21 * yb0 + trns31 * zb0;
        yb0d = trns12 * xb0 + trns22 * yb0 + trns32 * zb0;
        xc0d = trns11 * xc0 + trns21 * yc0 + trns31 * zc0;
        yc0d = trns12 * xc0 + trns22 * yc0 + trns32 * zc0;
        /*
           xa1d = trns11 * xa1 + trns21 * ya1 + trns31 * za1;
           ya1d = trns12 * xa1 + trns22 * ya1 + trns32 * za1;
         */
        za1d = trns13 * xa1 + trns23 * ya1 + trns33 * za1;
        xb1d = trns11 * xb1 + trns21 * yb1 + trns31 * zb1;
        yb1d = trns12 * xb1 + trns22 * yb1 + trns32 * zb1;
        zb1d = trns13 * xb1 + trns23 * yb1 + trns33 * zb1;
        xc1d = trns11 * xc1 + trns21 * yc1 + trns31 * zc1;
        yc1d = trns12 * xc1 + trns22 * yc1 + trns32 * zc1;
        zc1d = trns13 * xc1 + trns23 * yc1 + trns33 * zc1;
        /* 65 flops */

        sinphi = za1d * gmx_invsqrt(ra*ra);
        tmp    = 1.0 - sinphi * sinphi;
        if (tmp <= 0)
        {
            bOK = FALSE;
        }
        else
        {
            tmp2   = gmx_invsqrt(tmp);
            cosphi = tmp*tmp2;
            sinpsi = (zb1d - zc1d) * irc2 * tmp2;
            tmp2   = 1.0 - sinpsi * sinpsi;
            if (tmp2 <= 0)
            {
                bOK = FALSE;
            }
            else
            {
                cospsi = tmp2*gmx_invsqrt(tmp2);
            }
        }
        /* 46 flops */

        if (bOK)
        {
            ya2d =  ra * cosphi;
            xb2d = -rc * cospsi;
            t1   = -rb * cosphi;
            t2   =  rc * sinpsi * sinphi;
            yb2d =  t1 - t2;
            yc2d =  t1 + t2;
            /* 7 flops */

            /*     --- Step3  al,be,ga            --- */
            alpa   = xb2d * (xb0d - xc0d) + yb0d * yb2d + yc0d * yc2d;
            beta   = xb2d * (yc0d - yb0d) + xb0d * yb2d + xc0d * yc2d;
            gama   = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d;
            al2be2 = alpa * alpa + beta * beta;
            tmp2   = (al2be2 - gama * gama);
            sinthe = (alpa * gama - beta * tmp2*gmx_invsqrt(tmp2)) * gmx_invsqrt(al2be2*al2be2);
            /* 47 flops */

            /*  --- Step4  A3' --- */
            tmp2   = 1.0 - sinthe * sinthe;
            costhe = tmp2*gmx_invsqrt(tmp2);
            xa3d   = -ya2d * sinthe;
            ya3d   = ya2d * costhe;
            za3d   = za1d;
            xb3d   = xb2d * costhe - yb2d * sinthe;
            yb3d   = xb2d * sinthe + yb2d * costhe;
            zb3d   = zb1d;
            xc3d   = -xb2d * costhe - yc2d * sinthe;
            yc3d   = -xb2d * sinthe + yc2d * costhe;
            zc3d   = zc1d;
            /* 26 flops */

            /*    --- Step5  A3 --- */
            xa3 = trns11 * xa3d + trns12 * ya3d + trns13 * za3d;
            ya3 = trns21 * xa3d + trns22 * ya3d + trns23 * za3d;
            za3 = trns31 * xa3d + trns32 * ya3d + trns33 * za3d;
            xb3 = trns11 * xb3d + trns12 * yb3d + trns13 * zb3d;
            yb3 = trns21 * xb3d + trns22 * yb3d + trns23 * zb3d;
            zb3 = trns31 * xb3d + trns32 * yb3d + trns33 * zb3d;
            xc3 = trns11 * xc3d + trns12 * yc3d + trns13 * zc3d;
            yc3 = trns21 * xc3d + trns22 * yc3d + trns23 * zc3d;
            zc3 = trns31 * xc3d + trns32 * yc3d + trns33 * zc3d;
            /* 45 flops */
            after[ow1]     = xcom + xa3;
            after[ow1 + 1] = ycom + ya3;
            after[ow1 + 2] = zcom + za3;
            after[hw2]     = xcom + xb3;
            after[hw2 + 1] = ycom + yb3;
            after[hw2 + 2] = zcom + zb3;
            after[hw3]     = xcom + xc3;
            after[hw3 + 1] = ycom + yc3;
            after[hw3 + 2] = zcom + zc3;
            /* 9 flops */

            if (pbc != NULL)
            {
                rvec_inc(after+hw2, sh_hw2);
                rvec_inc(after+hw3, sh_hw3);
            }

            dax = xa3 - xa1;
            day = ya3 - ya1;
            daz = za3 - za1;
            dbx = xb3 - xb1;
            dby = yb3 - yb1;
            dbz = zb3 - zb1;
            dcx = xc3 - xc1;
            dcy = yc3 - yc1;
            dcz = zc3 - zc1;
            /* 9 flops, counted with the virial */

            if (v != NULL)
            {
                v[ow1]     += dax*invdt;
                v[ow1 + 1] += day*invdt;
                v[ow1 + 2] += daz*invdt;
                v[hw2]     += dbx*invdt;
                v[hw2 + 1] += dby*invdt;
                v[hw2 + 2] += dbz*invdt;
                v[hw3]     += dcx*invdt;
                v[hw3 + 1] += dcy*invdt;
                v[hw3 + 2] += dcz*invdt;
                /* 3*6 flops */
            }

            if (ow1 < CalcVirAtomEnd)
            {
                mdax                = mO*dax;
                mday                = mO*day;
                mdaz                = mO*daz;
                mdbx                = mH*dbx;
                mdby                = mH*dby;
                mdbz                = mH*dbz;
                mdcx                = mH*dcx;
                mdcy                = mH*dcy;
                mdcz                = mH*dcz;
                vir_r_m_dr[XX][XX] -= b4[ow1  ]*mdax + (b4[ow1  ]+xb0)*mdbx + (b4[ow1  ]+xc0)*mdcx;
                vir_r_m_dr[XX][YY] -= b4[ow1  ]*mday + (b4[ow1  ]+xb0)*mdby + (b4[ow1  ]+xc0)*mdcy;
                vir_r_m_dr[XX][ZZ] -= b4[ow1  ]*mdaz + (b4[ow1  ]+xb0)*mdbz + (b4[ow1  ]+xc0)*mdcz;
                vir_r_m_dr[YY][XX] -= b4[ow1+1]*mdax + (b4[ow1+1]+yb0)*mdbx + (b4[ow1+1]+yc0)*mdcx;
                vir_r_m_dr[YY][YY] -= b4[ow1+1]*mday + (b4[ow1+1]+yb0)*mdby + (b4[ow1+1]+yc0)*mdcy;
                vir_r_m_dr[YY][ZZ] -= b4[ow1+1]*mdaz + (b4[ow1+1]+yb0)*mdbz + (b4[ow1+1]+yc0)*mdcz;
                vir_r_m_dr[ZZ][XX] -= b4[ow1+2]*mdax + (b4[ow1+2]+zb0)*mdbx + (b4[ow1+2]+zc0)*mdcx;
                vir_r_m_dr[ZZ][YY] -= b4[ow1+2]*mday + (b4[ow1+2]+zb0)*mdby + (b4[ow1+2]+zc0)*mdcy;
                vir_r_m_dr[ZZ][ZZ] -= b4[ow1+2]*mdaz + (b4[ow1+2]+zb0)*mdbz + (b4[ow1+2]+zc0)*mdcz;
                /* 3*24 - 9 flops */
            }
        }
        else
        {
            *error = i;
        }
#ifdef DEBUG
        if (debug)
        {
            check_cons(debug, "settle", after, ow1, hw2, hw3);
        }
#endif
    }
}
