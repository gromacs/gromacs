/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gmx_arpack.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "gmx_blas.h"
#include "gmx_lapack.h"

#ifndef DOXYGEN

static void F77_FUNC(dstqrb, DSTQRB)(int* n, double* d__, double* e, double* z__, double* work, int* info)
{
    int    i__1, i__2;
    double d__1, d__2;
    int    c__0  = 0;
    int    c__1  = 1;
    double c_b31 = 1.;

    double b, c__, f, g;
    int    i__, j, k, l, m;
    double p, r__, s;
    int    l1, ii, mm, lm1, mm1, nm1;
    double rt1, rt2, eps;
    int    lsv;
    double tst, eps2;
    int    lend, jtot, lendm1, lendp1, iscale;

    int    lendsv, nmaxit, icompz;
    double ssfmax, ssfmin, safmin, minval, safmax, anorm;


    --work;
    --z__;
    --e;
    --d__;

    *info = 0;

    icompz = 2;

    if (*n == 0)
    {
        return;
    }

    if (*n == 1)
    {
        z__[1] = 1.;
        return;
    }

    eps = GMX_DOUBLE_EPS;

    d__1   = eps;
    eps2   = d__1 * d__1;
    minval = GMX_DOUBLE_MIN;
    safmin = minval / GMX_DOUBLE_EPS;
    safmax = 1. / safmin;
    ssfmax = std::sqrt(safmax) / 3.;
    ssfmin = std::sqrt(safmin) / eps2;

    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        z__[j] = 0.;
    }
    z__[*n] = 1.;

    nmaxit = *n * 30;
    jtot   = 0;

    l1  = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n)
    {
        goto L160;
    }
    if (l1 > 1)
    {
        e[l1 - 1] = 0.;
    }
    if (l1 <= nm1)
    {
        i__1 = nm1;
        for (m = l1; m <= i__1; ++m)
        {
            tst = std::abs(e[m]);
            if (tst == 0.)
            {
                goto L30;
            }
            if (tst <= std::sqrt(std::abs(d__[m])) * std::sqrt(std::abs(d__[m + 1])) * eps)
            {
                e[m] = 0.;
                goto L30;
            }
        }
    }
    m = *n;

L30:
    l      = l1;
    lsv    = l;
    lend   = m;
    lendsv = lend;
    l1     = m + 1;
    if (lend == l)
    {
        goto L10;
    }

    i__1   = lend - l + 1;
    anorm  = F77_FUNC(dlanst, DLANST)("i", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm == 0.)
    {
        goto L10;
    }
    if (anorm > ssfmax)
    {
        iscale = 1;
        i__1   = lend - l + 1;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, info);
        i__1 = lend - l;
        F77_FUNC(dlascl, DLASCL)("g", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, info);
    }
    else if (anorm < ssfmin)
    {
        iscale = 2;
        i__1   = lend - l + 1;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, info);
        i__1 = lend - l;
        F77_FUNC(dlascl, DLASCL)("g", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, info);
    }

    if (std::abs(d__[lend]) < std::abs(d__[l]))
    {
        lend = lsv;
        l    = lendsv;
    }

    if (lend > l)
    {

    L40:
        if (l != lend)
        {
            lendm1 = lend - 1;
            i__1   = lendm1;
            for (m = l; m <= i__1; ++m)
            {
                d__2 = std::abs(e[m]);
                tst  = d__2 * d__2;
                if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m + 1]) + safmin)
                {
                    goto L60;
                }
            }
        }

        m = lend;

    L60:
        if (m < lend)
        {
            e[m] = 0.;
        }
        p = d__[l];
        if (m == l)
        {
            goto L80;
        }

        if (m == l + 1)
        {
            if (icompz > 0)
            {
                F77_FUNC(dlaev2, DLAEV2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
                work[l]          = c__;
                work[*n - 1 + l] = s;

                tst        = z__[l + 1];
                z__[l + 1] = c__ * tst - s * z__[l];
                z__[l]     = s * tst + c__ * z__[l];
            }
            else
            {
                F77_FUNC(dlae2, DLAE2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
            }
            d__[l]     = rt1;
            d__[l + 1] = rt2;
            e[l]       = 0.;
            l += 2;
            if (l <= lend)
            {
                goto L40;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;

        g   = (d__[l + 1] - p) / (e[l] * 2.);
        r__ = F77_FUNC(dlapy2, DLAPY2)(&g, &c_b31);
        g   = d__[m] - p + e[l] / (g + ((g > 0) ? r__ : -r__));

        s   = 1.;
        c__ = 1.;
        p   = 0.;

        mm1  = m - 1;
        i__1 = l;
        for (i__ = mm1; i__ >= i__1; --i__)
        {
            f = s * e[i__];
            b = c__ * e[i__];
            F77_FUNC(dlartg, DLARTG)(&g, &f, &c__, &s, &r__);
            if (i__ != m - 1)
            {
                e[i__ + 1] = r__;
            }
            g            = d__[i__ + 1] - p;
            r__          = (d__[i__] - g) * s + c__ * 2. * b;
            p            = s * r__;
            d__[i__ + 1] = g + p;
            g            = c__ * r__ - b;

            if (icompz > 0)
            {
                work[i__]          = c__;
                work[*n - 1 + i__] = -s;
            }
        }

        if (icompz > 0)
        {
            mm = m - l + 1;

            F77_FUNC(dlasr, DLASR)
            ("r", "v", "b", &c__1, &mm, &work[l], &work[*n - 1 + l], &z__[l], &c__1);
        }

        d__[l] -= p;
        e[l] = g;
        goto L40;

    L80:
        d__[l] = p;

        ++l;
        if (l <= lend)
        {
            goto L40;
        }
        goto L140;
    }
    else
    {

    L90:
        if (l != lend)
        {
            lendp1 = lend + 1;
            i__1   = lendp1;
            for (m = l; m >= i__1; --m)
            {
                d__2 = std::abs(e[m - 1]);
                tst  = d__2 * d__2;
                if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m - 1]) + safmin)
                {
                    goto L110;
                }
            }
        }

        m = lend;

    L110:
        if (m > lend)
        {
            e[m - 1] = 0.;
        }
        p = d__[l];
        if (m == l)
        {
            goto L130;
        }

        if (m == l - 1)
        {
            if (icompz > 0)
            {
                F77_FUNC(dlaev2, DLAEV2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s);

                tst        = z__[l];
                z__[l]     = c__ * tst - s * z__[l - 1];
                z__[l - 1] = s * tst + c__ * z__[l - 1];
            }
            else
            {
                F77_FUNC(dlae2, DLAE2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
            }
            d__[l - 1] = rt1;
            d__[l]     = rt2;
            e[l - 1]   = 0.;
            l += -2;
            if (l >= lend)
            {
                goto L90;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;


        g   = (d__[l - 1] - p) / (e[l - 1] * 2.);
        r__ = F77_FUNC(dlapy2, DLAPY2)(&g, &c_b31);
        g   = d__[m] - p + e[l - 1] / (g + ((g > 0) ? r__ : -r__));

        s   = 1.;
        c__ = 1.;
        p   = 0.;

        lm1  = l - 1;
        i__1 = lm1;
        for (i__ = m; i__ <= i__1; ++i__)
        {
            f = s * e[i__];
            b = c__ * e[i__];
            F77_FUNC(dlartg, DLARTG)(&g, &f, &c__, &s, &r__);
            if (i__ != m)
            {
                e[i__ - 1] = r__;
            }
            g        = d__[i__] - p;
            r__      = (d__[i__ + 1] - g) * s + c__ * 2. * b;
            p        = s * r__;
            d__[i__] = g + p;
            g        = c__ * r__ - b;

            if (icompz > 0)
            {
                work[i__]          = c__;
                work[*n - 1 + i__] = s;
            }
        }

        if (icompz > 0)
        {
            mm = l - m + 1;

            F77_FUNC(dlasr, DLASR)
            ("r", "v", "f", &c__1, &mm, &work[m], &work[*n - 1 + m], &z__[m], &c__1);
        }

        d__[l] -= p;
        e[lm1] = g;
        goto L90;

    L130:
        d__[l] = p;

        --l;
        if (l >= lend)
        {
            goto L90;
        }
        goto L140;
    }

L140:
    if (iscale == 1)
    {
        i__1 = lendsv - lsv + 1;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], n, info);
        i__1 = lendsv - lsv;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }
    else if (iscale == 2)
    {
        i__1 = lendsv - lsv + 1;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], n, info);
        i__1 = lendsv - lsv;
        F77_FUNC(dlascl, DLASCL)
        ("g", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }

    if (jtot < nmaxit)
    {
        goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        if (e[i__] != 0.)
        {
            ++(*info);
        }
    }
    goto L190;

L160:
    if (icompz == 0)
    {

        F77_FUNC(dlasrt, DLASRT)("i", n, &d__[1], info);
    }
    else
    {

        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii)
        {
            i__  = ii - 1;
            k    = i__;
            p    = d__[i__];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j)
            {
                if (d__[j] < p)
                {
                    k = j;
                    p = d__[j];
                }
            }
            if (k != i__)
            {
                d__[k]   = d__[i__];
                d__[i__] = p;

                p        = z__[k];
                z__[k]   = z__[i__];
                z__[i__] = p;
            }
        }
    }

L190:
    return;
}

static void F77_FUNC(dgetv0, DGETV0)(int*        ido,
                                     const char* bmat,
                                     int gmx_unused* itry,
                                     int*            initv,
                                     int*            n,
                                     int*            j,
                                     double*         v,
                                     int*            ldv,
                                     double*         resid,
                                     double*         rnorm,
                                     int*            ipntr,
                                     double*         workd,
                                     int*            iwork,
                                     int*            ierr)
{
    int    c__1  = 1;
    double c_b22 = 1.;
    double c_b24 = 0.;
    double c_b27 = -1.;
    int    v_dim1, v_offset, i__1;

    int jj;
    int idist;

    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;
    --iwork;

    if (*ido == 0)
    {

        *ierr    = 0;
        iwork[7] = 0;
        iwork[5] = 0;
        iwork[6] = 0;

        if (!(*initv))
        {
            idist = 2;
            F77_FUNC(dlarnv, DLARNV)(&idist, &iwork[1], n, &resid[1]);
        }

        if (*bmat == 'G')
        {
            ipntr[1] = 1;
            ipntr[2] = *n + 1;
            F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
            *ido = -1;
            goto L9000;
        }
    }

    if (iwork[5] == 1)
    {
        goto L20;
    }

    if (iwork[6] == 1)
    {
        goto L40;
    }

    iwork[5] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(dcopy, DCOPY)(n, &workd[*n + 1], &c__1, &resid[1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L20:


    iwork[5] = 0;
    if (*bmat == 'G')
    {
        workd[*n * 3 + 4] = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        workd[*n * 3 + 4] = std::sqrt(std::abs(workd[*n * 3 + 4]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 4] = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }
    *rnorm = workd[*n * 3 + 4];

    if (*j == 1)
    {
        goto L50;
    }
    iwork[6] = 1;
L30:

    i__1 = *j - 1;
    F77_FUNC(dgemv, DGEMV)
    ("T", n, &i__1, &c_b22, &v[v_offset], ldv, &workd[1], &c__1, &c_b24, &workd[*n + 1], &c__1);
    i__1 = *j - 1;
    F77_FUNC(dgemv, DGEMV)
    ("N", n, &i__1, &c_b27, &v[v_offset], ldv, &workd[*n + 1], &c__1, &c_b22, &resid[1], &c__1);

    if (*bmat == 'G')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L40:

    if (*bmat == 'G')
    {
        *rnorm = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        *rnorm = std::sqrt(std::abs(*rnorm));
    }
    else if (*bmat == 'I')
    {
        *rnorm = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }

    if (*rnorm > workd[*n * 3 + 4] * .717F)
    {
        goto L50;
    }

    ++iwork[7];
    if (iwork[7] <= 1)
    {

        workd[*n * 3 + 4] = *rnorm;
        goto L30;
    }
    else
    {

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.;
        }
        *rnorm = 0.;
        *ierr  = -1;
    }

L50:

    *ido = 99;

L9000:
    return;
}


static void F77_FUNC(dsapps, DSAPPS)(int*    n,
                                     int*    kev,
                                     int*    np,
                                     double* shift,
                                     double* v,
                                     int*    ldv,
                                     double* h__,
                                     int*    ldh,
                                     double* resid,
                                     double* q,
                                     int*    ldq,
                                     double* workd)
{
    double c_b4  = 0.;
    double c_b5  = 1.;
    double c_b14 = -1.;
    int    c__1  = 1;
    int    h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    double c__, f, g;
    int    i__, j;
    double r__, s, a1, a2, a3, a4;
    int    jj;
    double big;
    int    iend, itop;
    double epsmch;
    int    istart, kplusp;

    --workd;
    --resid;
    --shift;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1   = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    epsmch = GMX_DOUBLE_EPS;
    itop   = 1;


    kplusp = *kev + *np;

    F77_FUNC(dlaset, DLASET)("All", &kplusp, &kplusp, &c_b4, &c_b5, &q[q_offset], ldq);

    if (*np == 0)
    {
        goto L9000;
    }

    i__1 = *np;
    for (jj = 1; jj <= i__1; ++jj)
    {

        istart = itop;

    L20:

        i__2 = kplusp - 1;
        for (i__ = istart; i__ <= i__2; ++i__)
        {
            big = std::abs(h__[i__ + (h_dim1 * 2)]) + std::abs(h__[i__ + 1 + (h_dim1 * 2)]);
            if (h__[i__ + 1 + h_dim1] <= epsmch * big)
            {
                h__[i__ + 1 + h_dim1] = 0.;
                iend                  = i__;
                goto L40;
            }
        }
        iend = kplusp;
    L40:

        if (istart < iend)
        {

            f = h__[istart + (h_dim1 << 1)] - shift[jj];
            g = h__[istart + 1 + h_dim1];
            F77_FUNC(dlartg, DLARTG)(&f, &g, &c__, &s, &r__);

            a1 = c__ * h__[istart + (h_dim1 << 1)] + s * h__[istart + 1 + h_dim1];
            a2 = c__ * h__[istart + 1 + h_dim1] + s * h__[istart + 1 + (h_dim1 << 1)];
            a4 = c__ * h__[istart + 1 + (h_dim1 << 1)] - s * h__[istart + 1 + h_dim1];
            a3 = c__ * h__[istart + 1 + h_dim1] - s * h__[istart + (h_dim1 << 1)];
            h__[istart + (h_dim1 << 1)]     = c__ * a1 + s * a2;
            h__[istart + 1 + (h_dim1 << 1)] = c__ * a4 - s * a3;
            h__[istart + 1 + h_dim1]        = c__ * a3 + s * a4;

            i__3 = istart + jj;
            i__2 = (i__3 < kplusp) ? i__3 : kplusp;
            for (j = 1; j <= i__2; ++j)
            {
                a1 = c__ * q[j + istart * q_dim1] + s * q[j + (istart + 1) * q_dim1];
                q[j + (istart + 1) * q_dim1] =
                        -s * q[j + istart * q_dim1] + c__ * q[j + (istart + 1) * q_dim1];
                q[j + istart * q_dim1] = a1;
            }

            i__2 = iend - 1;
            for (i__ = istart + 1; i__ <= i__2; ++i__)
            {

                f = h__[i__ + h_dim1];
                g = s * h__[i__ + 1 + h_dim1];

                h__[i__ + 1 + h_dim1] = c__ * h__[i__ + 1 + h_dim1];
                F77_FUNC(dlartg, DLARTG)(&f, &g, &c__, &s, &r__);

                if (r__ < 0.)
                {
                    r__ = -r__;
                    c__ = -c__;
                    s   = -s;
                }

                h__[i__ + h_dim1] = r__;

                a1 = c__ * h__[i__ + (h_dim1 << 1)] + s * h__[i__ + 1 + h_dim1];
                a2 = c__ * h__[i__ + 1 + h_dim1] + s * h__[i__ + 1 + (h_dim1 << 1)];
                a3 = c__ * h__[i__ + 1 + h_dim1] - s * h__[i__ + (h_dim1 << 1)];
                a4 = c__ * h__[i__ + 1 + (h_dim1 << 1)] - s * h__[i__ + 1 + h_dim1];

                h__[i__ + (h_dim1 << 1)]     = c__ * a1 + s * a2;
                h__[i__ + 1 + (h_dim1 << 1)] = c__ * a4 - s * a3;
                h__[i__ + 1 + h_dim1]        = c__ * a3 + s * a4;

                i__4 = j + jj;
                i__3 = (i__4 < kplusp) ? i__4 : kplusp;
                for (j = 1; j <= i__3; ++j)
                {
                    a1 = c__ * q[j + i__ * q_dim1] + s * q[j + (i__ + 1) * q_dim1];
                    q[j + (i__ + 1) * q_dim1] = -s * q[j + i__ * q_dim1] + c__ * q[j + (i__ + 1) * q_dim1];
                    q[j + i__ * q_dim1] = a1;
                }
            }
        }

        istart = iend + 1;

        if (h__[iend + h_dim1] < 0.)
        {
            h__[iend + h_dim1] = -h__[iend + h_dim1];
            F77_FUNC(dscal, DSCAL)(&kplusp, &c_b14, &q[iend * q_dim1 + 1], &c__1);
        }

        if (iend < kplusp)
        {
            goto L20;
        }

        i__2 = kplusp - 1;
        for (i__ = itop; i__ <= i__2; ++i__)
        {
            if (h__[i__ + 1 + h_dim1] > 0.)
            {
                goto L90;
            }
            ++itop;
        }

    L90:;
    }

    i__1 = kplusp - 1;
    for (i__ = itop; i__ <= i__1; ++i__)
    {
        big = std::abs(h__[i__ + (h_dim1 * 2)]) + std::abs(h__[i__ + 1 + (h_dim1 * 2)]);
        if (h__[i__ + 1 + h_dim1] <= epsmch * big)
        {
            h__[i__ + 1 + h_dim1] = 0.;
        }
    }

    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(dgemv, DGEMV)
        ("N", n, &kplusp, &c_b5, &v[v_offset], ldv, &q[(*kev + 1) * q_dim1 + 1], &c__1, &c_b4, &workd[*n + 1], &c__1);
    }

    i__1 = *kev;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = kplusp - i__ + 1;
        F77_FUNC(dgemv, DGEMV)
        ("N", n, &i__2, &c_b5, &v[v_offset], ldv, &q[(*kev - i__ + 1) * q_dim1 + 1], &c__1, &c_b4, &workd[1], &c__1);
        F77_FUNC(dcopy, DCOPY)(n, &workd[1], &c__1, &v[(kplusp - i__ + 1) * v_dim1 + 1], &c__1);
    }

    F77_FUNC(dlacpy, DLACPY)("All", n, kev, &v[(*np + 1) * v_dim1 + 1], ldv, &v[v_offset], ldv);

    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(dcopy, DCOPY)(n, &workd[*n + 1], &c__1, &v[(*kev + 1) * v_dim1 + 1], &c__1);
    }

    F77_FUNC(dscal, DSCAL)(n, &q[kplusp + *kev * q_dim1], &resid[1], &c__1);
    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(daxpy, DAXPY)
        (n, &h__[*kev + 1 + h_dim1], &v[(*kev + 1) * v_dim1 + 1], &c__1, &resid[1], &c__1);
    }


L9000:
    return;
}


static void F77_FUNC(dsortr, DSORTR)(const char* which, int* apply, int* n, double* x1, double* x2)
{
    int i__1;

    int    i__, j, igap;
    double temp;


    igap = *n / 2;

    if (!std::strncmp(which, "SA", 2))
    {

    L10:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L20:

            if (j < 0)
            {
                goto L30;
            }

            if (x1[j] < x1[j + igap])
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L30;
            }
            j -= igap;
            goto L20;
        L30:;
        }
        igap /= 2;
        goto L10;
    }
    else if (!std::strncmp(which, "SM", 2))
    {

    L40:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L50:

            if (j < 0)
            {
                goto L60;
            }

            if (std::abs(x1[j]) < std::abs(x1[j + igap]))
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L60;
            }
            j -= igap;
            goto L50;
        L60:;
        }
        igap /= 2;
        goto L40;
    }
    else if (!std::strncmp(which, "LA", 2))
    {

    L70:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L80:

            if (j < 0)
            {
                goto L90;
            }

            if (x1[j] > x1[j + igap])
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L90;
            }
            j -= igap;
            goto L80;
        L90:;
        }
        igap /= 2;
        goto L70;
    }
    else if (!std::strncmp(which, "LM", 2))
    {


    L100:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L110:

            if (j < 0)
            {
                goto L120;
            }

            if (std::abs(x1[j]) > std::abs(x1[j + igap]))
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L120;
            }
            j -= igap;
            goto L110;
        L120:;
        }
        igap /= 2;
        goto L100;
    }

L9000:
    return;
}


static void F77_FUNC(dsesrt,
                     DSESRT)(const char* which, int* apply, int* n, double* x, int* na, double* a, int* lda)
{
    int a_dim1, a_offset, i__1;
    int c__1 = 1;

    int    i__, j, igap;
    double temp;

    a_dim1   = *lda;
    a_offset = 1 + a_dim1 * 0;
    a -= a_offset;

    igap = *n / 2;

    if (!std::strncmp(which, "SA", 2))
    {

    L10:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L20:

            if (j < 0)
            {
                goto L30;
            }

            if (x[j] < x[j + igap])
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(dswap, DSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L30;
            }
            j -= igap;
            goto L20;
        L30:;
        }
        igap /= 2;
        goto L10;
    }
    else if (!std::strncmp(which, "SM", 2))
    {

    L40:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L50:

            if (j < 0)
            {
                goto L60;
            }

            if (std::abs(x[j]) < std::abs(x[j + igap]))
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(dswap, DSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L60;
            }
            j -= igap;
            goto L50;
        L60:;
        }
        igap /= 2;
        goto L40;
    }
    else if (!std::strncmp(which, "LA", 2))
    {

    L70:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L80:

            if (j < 0)
            {
                goto L90;
            }

            if (x[j] > x[j + igap])
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(dswap, DSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L90;
            }
            j -= igap;
            goto L80;
        L90:;
        }
        igap /= 2;
        goto L70;
    }
    else if (!std::strncmp(which, "LM", 2))
    {

    L100:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L110:

            if (j < 0)
            {
                goto L120;
            }

            if (std::abs(x[j]) > std::abs(x[j + igap]))
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(dswap, DSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L120;
            }
            j -= igap;
            goto L110;
        L120:;
        }
        igap /= 2;
        goto L100;
    }

L9000:
    return;
}


static void F77_FUNC(dsgets,
                     DSGETS)(int* ishift, const char* which, int* kev, int* np, double* ritz, double* bounds, double* shifts)
{
    int c__1 = 1;
    int i__1, i__2;
    int kevd2;

    --shifts;
    --bounds;
    --ritz;

    if (!std::strncmp(which, "BE", 2))
    {
        i__1 = *kev + *np;
        F77_FUNC(dsortr, DSORTR)("LA", &c__1, &i__1, &ritz[1], &bounds[1]);
        kevd2 = *kev / 2;
        if (*kev > 1)
        {
            i__1 = (kevd2 < *np) ? kevd2 : *np;
            i__2 = (kevd2 > *np) ? kevd2 : *np;
            F77_FUNC(dswap, DSWAP)(&i__1, &ritz[1], &c__1, &ritz[i__2 + 1], &c__1);
            i__1 = (kevd2 < *np) ? kevd2 : *np;
            i__2 = (kevd2 > *np) ? kevd2 : *np;
            F77_FUNC(dswap, DSWAP)(&i__1, &bounds[1], &c__1, &bounds[i__2 + 1], &c__1);
        }
    }
    else
    {
        i__1 = *kev + *np;
        F77_FUNC(dsortr, DSORTR)(which, &c__1, &i__1, &ritz[1], &bounds[1]);
    }

    if (*ishift == 1 && *np > 0)
    {

        F77_FUNC(dsortr, DSORTR)("SM", &c__1, np, &bounds[1], &ritz[1]);
        F77_FUNC(dcopy, DCOPY)(np, &ritz[1], &c__1, &shifts[1], &c__1);
    }


    return;
}


static void F77_FUNC(dsconv, DSCONV)(int* n, double* ritz, double* bounds, double* tol, int* nconv)
{
    double c_b3 = 2 / 3.;
    int    i__1;
    double d__2, d__3;

    int    i__;
    double eps23, temp;

    --bounds;
    --ritz;

    eps23 = GMX_DOUBLE_EPS;
    eps23 = std::pow(eps23, c_b3);

    *nconv = 0;
    i__1   = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {

        d__2 = eps23;
        d__3 = std::abs(ritz[i__]);
        temp = (d__2 > d__3) ? d__2 : d__3;
        if (bounds[i__] <= *tol * temp)
        {
            ++(*nconv);
        }
    }

    return;
}


static void F77_FUNC(dseigt, DSEIGT)(double* rnorm,
                                     int*    n,
                                     double* h__,
                                     int*    ldh,
                                     double* eig,
                                     double* bounds,
                                     double* workl,
                                     int*    ierr)
{
    int c__1 = 1;
    int h_dim1, h_offset, i__1;

    int k;


    --workl;
    --bounds;
    --eig;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    F77_FUNC(dcopy, DCOPY)(n, &h__[(h_dim1 << 1) + 1], &c__1, &eig[1], &c__1);
    i__1 = *n - 1;
    F77_FUNC(dcopy, DCOPY)(&i__1, &h__[h_dim1 + 2], &c__1, &workl[1], &c__1);
    F77_FUNC(dstqrb, DSTQRB)(n, &eig[1], &workl[1], &bounds[1], &workl[*n + 1], ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }

    i__1 = *n;
    for (k = 1; k <= i__1; ++k)
    {
        bounds[k] = *rnorm * std::abs(bounds[k]);
    }


L9000:
    return;
}


static void F77_FUNC(dsaitr, DSAITR)(int*        ido,
                                     const char* bmat,
                                     int*        n,
                                     int*        k,
                                     int*        np,
                                     int*        mode,
                                     double*     resid,
                                     double*     rnorm,
                                     double*     v,
                                     int*        ldv,
                                     double*     h__,
                                     int*        ldh,
                                     int*        ipntr,
                                     double*     workd,
                                     int*        iwork,
                                     int*        info)
{

    int    c__0  = 0;
    int    c__1  = 1;
    double c_b18 = 1.;
    double c_b42 = 0.;
    double c_b50 = -1.;

    int    h_dim1, h_offset, v_dim1, v_offset, i__1;
    int    i__, jj;
    double temp1;
    int    infol;
    double safmin, minval;


    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipntr;
    --iwork;
    minval = GMX_DOUBLE_MIN;
    safmin = minval / GMX_DOUBLE_EPS;

    if (*ido == 0)
    {
        *info    = 0;
        iwork[5] = 0;
        iwork[6] = 0;
        iwork[4] = 0;
        iwork[2] = 0;
        iwork[3] = 0;

        iwork[12] = *k + 1;

        iwork[8]  = 1;
        iwork[9]  = iwork[8] + *n;
        iwork[10] = iwork[9] + *n;
    }

    if (iwork[5] == 1)
    {
        goto L50;
    }
    if (iwork[6] == 1)
    {
        goto L60;
    }
    if (iwork[2] == 1)
    {
        goto L70;
    }
    if (iwork[3] == 1)
    {
        goto L90;
    }
    if (iwork[4] == 1)
    {
        goto L30;
    }

L1000:


    if (*rnorm > 0.)
    {
        goto L40;
    }

    iwork[11] = 1;
L20:
    iwork[4] = 1;
    *ido     = 0;
L30:

    F77_FUNC(dgetv0, DGETV0)
    (ido,
     bmat,
     &iwork[11],
     &c__0,
     n,
     &iwork[12],
     &v[v_offset],
     ldv,
     &resid[1],
     rnorm,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     &iwork[7]);
    if (*ido != 99)
    {
        goto L9000;
    }
    if (iwork[7] < 0)
    {
        ++iwork[11];
        if (iwork[11] <= 3)
        {
            goto L20;
        }

        *info = iwork[12] - 1;
        *ido  = 99;
        goto L9000;
    }

L40:

    F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &v[iwork[12] * v_dim1 + 1], &c__1);
    if (*rnorm >= safmin)
    {
        temp1 = 1. / *rnorm;
        F77_FUNC(dscal, DSCAL)(n, &temp1, &v[iwork[12] * v_dim1 + 1], &c__1);
        F77_FUNC(dscal, DSCAL)(n, &temp1, &workd[iwork[8]], &c__1);
    }
    else
    {

        F77_FUNC(dlascl, DLASCL)
        ("General", &i__, &i__, rnorm, &c_b18, n, &c__1, &v[iwork[12] * v_dim1 + 1], n, &infol);
        F77_FUNC(dlascl, DLASCL)
        ("General", &i__, &i__, rnorm, &c_b18, n, &c__1, &workd[iwork[8]], n, &infol);
    }

    iwork[5] = 1;
    F77_FUNC(dcopy, DCOPY)(n, &v[iwork[12] * v_dim1 + 1], &c__1, &workd[iwork[10]], &c__1);
    ipntr[1] = iwork[10];
    ipntr[2] = iwork[9];
    ipntr[3] = iwork[8];
    *ido     = 1;

    goto L9000;
L50:


    iwork[5] = 0;

    F77_FUNC(dcopy, DCOPY)(n, &workd[iwork[9]], &c__1, &resid[1], &c__1);

    if (*mode == 2)
    {
        goto L65;
    }
    if (*bmat == 'G')
    {
        iwork[6] = 1;
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L60:

    iwork[6] = 0;

L65:
    if (*mode == 2)
    {

        workd[*n * 3 + 3] = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[iwork[10]], &c__1);
        workd[*n * 3 + 3] = std::sqrt(std::abs(workd[*n * 3 + 3]));
    }
    else if (*bmat == 'G')
    {
        workd[*n * 3 + 3] = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        workd[*n * 3 + 3] = std::sqrt(std::abs(workd[*n * 3 + 3]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 3] = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }

    if (*mode != 2)
    {
        F77_FUNC(dgemv, DGEMV)
        ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[8]], &c__1, &c_b42, &workd[iwork[9]], &c__1);
    }
    else
    {
        F77_FUNC(dgemv, DGEMV)
        ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[10]], &c__1, &c_b42, &workd[iwork[9]], &c__1);
    }

    F77_FUNC(dgemv, DGEMV)
    ("N", n, &iwork[12], &c_b50, &v[v_offset], ldv, &workd[iwork[9]], &c__1, &c_b18, &resid[1], &c__1);

    h__[iwork[12] + (h_dim1 << 1)] = workd[iwork[9] + iwork[12] - 1];
    if (iwork[12] == 1 || iwork[4] == 1)
    {
        h__[iwork[12] + h_dim1] = 0.;
    }
    else
    {
        h__[iwork[12] + h_dim1] = *rnorm;
    }

    iwork[2] = 1;
    iwork[1] = 0;

    if (*bmat == 'G')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[iwork[9]], &c__1);
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L70:

    iwork[2] = 0;

    if (*bmat == 'G')
    {
        *rnorm = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        *rnorm = std::sqrt(std::abs(*rnorm));
    }
    else if (*bmat == 'I')
    {
        *rnorm = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }

    if (*rnorm > workd[*n * 3 + 3] * .717F)
    {
        goto L100;
    }

L80:

    F77_FUNC(dgemv, DGEMV)
    ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[8]], &c__1, &c_b42, &workd[iwork[9]], &c__1);

    F77_FUNC(dgemv, DGEMV)
    ("N", n, &iwork[12], &c_b50, &v[v_offset], ldv, &workd[iwork[9]], &c__1, &c_b18, &resid[1], &c__1);

    if (iwork[12] == 1 || iwork[4] == 1)
    {
        h__[iwork[12] + h_dim1] = 0.;
    }
    h__[iwork[12] + (h_dim1 << 1)] += workd[iwork[9] + iwork[12] - 1];

    iwork[3] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[iwork[9]], &c__1);
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L90:


    if (*bmat == 'G')
    {
        workd[*n * 3 + 2] = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        workd[*n * 3 + 2] = std::sqrt(std::abs(workd[*n * 3 + 2]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 2] = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }


    if (workd[*n * 3 + 2] > *rnorm * .717F)
    {

        *rnorm = workd[*n * 3 + 2];
    }
    else
    {

        *rnorm = workd[*n * 3 + 2];
        ++iwork[1];
        if (iwork[1] <= 1)
        {
            goto L80;
        }

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.;
        }
        *rnorm = 0.;
    }

L100:

    iwork[4] = 0;
    iwork[3] = 0;

    if (h__[iwork[12] + h_dim1] < 0.)
    {
        h__[iwork[12] + h_dim1] = -h__[iwork[12] + h_dim1];
        if (iwork[12] < *k + *np)
        {
            F77_FUNC(dscal, DSCAL)(n, &c_b50, &v[(iwork[12] + 1) * v_dim1 + 1], &c__1);
        }
        else
        {
            F77_FUNC(dscal, DSCAL)(n, &c_b50, &resid[1], &c__1);
        }
    }

    ++iwork[12];
    if (iwork[12] > *k + *np)
    {
        *ido = 99;


        goto L9000;
    }

    goto L1000;

L9000:
    return;
}


static void F77_FUNC(dsaup2, DSAUP2)(int*        ido,
                                     const char* bmat,
                                     int*        n,
                                     const char* which,
                                     int*        nev,
                                     int*        np,
                                     double*     tol,
                                     double*     resid,
                                     int*        mode,
                                     int gmx_unused* iupd,
                                     int*            ishift,
                                     int*            mxiter,
                                     double*         v,
                                     int*            ldv,
                                     double*         h__,
                                     int*            ldh,
                                     double*         ritz,
                                     double*         bounds,
                                     double*         q,
                                     int*            ldq,
                                     double*         workl,
                                     int*            ipntr,
                                     double*         workd,
                                     int*            iwork,
                                     int*            info)
{
    double c_b3 = 2 / 3.;
    int    c__1 = 1;
    int    c__0 = 0;

    int    h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, i__3;
    double d__2, d__3;
    int    j;
    double eps23;
    int    ierr;
    double temp;
    int    nevd2;
    int    nevm2;
    int    nevbef;
    char   wprime[2];
    int    nptemp;

    --workd;
    --resid;
    --workl;
    --bounds;
    --ritz;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1   = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --ipntr;
    --iwork;
    eps23 = GMX_DOUBLE_EPS;
    eps23 = std::pow(eps23, c_b3);

    if (*ido == 0)
    {

        iwork[41] = 1;
        iwork[42] = 3;
        iwork[43] = 5;
        iwork[44] = 7;

        iwork[9]  = *nev;
        iwork[10] = *np;

        iwork[7] = iwork[9] + iwork[10];
        iwork[8] = 0;
        iwork[6] = 0;

        iwork[2] = 1;
        iwork[4] = 0;
        iwork[5] = 0;
        iwork[1] = 0;

        if (*info != 0)
        {

            iwork[3] = 1;
            *info    = 0;
        }
        else
        {
            iwork[3] = 0;
        }
    }

    if (iwork[2] == 1)
    {
        F77_FUNC(dgetv0, DGETV0)
        (ido,
         bmat,
         &c__1,
         &iwork[3],
         n,
         &c__1,
         &v[v_offset],
         ldv,
         &resid[1],
         &workd[*n * 3 + 1],
         &ipntr[1],
         &workd[1],
         &iwork[41],
         info);

        if (*ido != 99)
        {
            goto L9000;
        }

        if (workd[*n * 3 + 1] == 0.)
        {

            *info = -9;
            goto L1200;
        }
        iwork[2] = 0;
        *ido     = 0;
    }

    if (iwork[4] == 1)
    {
        goto L20;
    }

    if (iwork[5] == 1)
    {
        goto L50;
    }

    if (iwork[1] == 1)
    {
        goto L100;
    }

    F77_FUNC(dsaitr, DSAITR)
    (ido,
     bmat,
     n,
     &c__0,
     &iwork[9],
     mode,
     &resid[1],
     &workd[*n * 3 + 1],
     &v[v_offset],
     ldv,
     &h__[h_offset],
     ldh,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {

        *np     = *info;
        *mxiter = iwork[6];
        *info   = -9999;
        goto L1200;
    }

L1000:

    ++iwork[6];


    *ido = 0;
L20:
    iwork[4] = 1;

    F77_FUNC(dsaitr, DSAITR)
    (ido,
     bmat,
     n,
     nev,
     np,
     mode,
     &resid[1],
     &workd[*n * 3 + 1],
     &v[v_offset],
     ldv,
     &h__[h_offset],
     ldh,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {

        *np     = *info;
        *mxiter = iwork[6];
        *info   = -9999;
        goto L1200;
    }
    iwork[4] = 0;

    F77_FUNC(dseigt, DSEIGT)
    (&workd[*n * 3 + 1], &iwork[7], &h__[h_offset], ldh, &ritz[1], &bounds[1], &workl[1], &ierr);

    if (ierr != 0)
    {
        *info = -8;
        goto L1200;
    }

    F77_FUNC(dcopy, DCOPY)(&iwork[7], &ritz[1], &c__1, &workl[iwork[7] + 1], &c__1);
    F77_FUNC(dcopy, DCOPY)(&iwork[7], &bounds[1], &c__1, &workl[(iwork[7] << 1) + 1], &c__1);

    *nev = iwork[9];
    *np  = iwork[10];
    F77_FUNC(dsgets, DSGETS)(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);

    F77_FUNC(dcopy, DCOPY)(nev, &bounds[*np + 1], &c__1, &workl[*np + 1], &c__1);
    F77_FUNC(dsconv, DSCONV)(nev, &ritz[*np + 1], &workl[*np + 1], tol, &iwork[8]);

    nptemp = *np;
    i__1   = nptemp;
    for (j = 1; j <= i__1; ++j)
    {
        if (bounds[j] == 0.)
        {
            --(*np);
            ++(*nev);
        }
    }

    if (iwork[8] >= iwork[9] || iwork[6] > *mxiter || *np == 0)
    {

        if (!std::strncmp(which, "BE", 2))
        {

            std::strncpy(wprime, "SA", 2);
            F77_FUNC(dsortr, DSORTR)(wprime, &c__1, &iwork[7], &ritz[1], &bounds[1]);
            nevd2 = *nev / 2;
            nevm2 = *nev - nevd2;
            if (*nev > 1)
            {
                i__1 = (nevd2 < *np) ? nevd2 : *np;
                i__2 = iwork[7] - nevd2 + 1, i__3 = iwork[7] - *np + 1;
                F77_FUNC(dswap, DSWAP)
                (&i__1, &ritz[nevm2 + 1], &c__1, &ritz[((i__2 > i__3) ? i__2 : i__3)], &c__1);
                i__1 = (nevd2 < *np) ? nevd2 : *np;
                i__2 = iwork[7] - nevd2 + 1, i__3 = iwork[7] - *np;
                F77_FUNC(dswap, DSWAP)
                (&i__1, &bounds[nevm2 + 1], &c__1, &bounds[((i__2 > i__3) ? i__2 : i__3) + 1], &c__1);
            }
        }
        else
        {

            if (!std::strncmp(which, "LM", 2))
            {
                std::strncpy(wprime, "SM", 2);
            }
            if (!std::strncmp(which, "SM", 2))
            {
                std::strncpy(wprime, "LM", 2);
            }
            if (!std::strncmp(which, "LA", 2))
            {
                std::strncpy(wprime, "SA", 2);
            }
            if (!std::strncmp(which, "SA", 2))
            {
                std::strncpy(wprime, "LA", 2);
            }

            F77_FUNC(dsortr, DSORTR)(wprime, &c__1, &iwork[7], &ritz[1], &bounds[1]);
        }

        i__1 = iwork[9];
        for (j = 1; j <= i__1; ++j)
        {
            d__2 = eps23;
            d__3 = std::abs(ritz[j]);
            temp = (d__2 > d__3) ? d__2 : d__3;
            bounds[j] /= temp;
        }

        std::strncpy(wprime, "LA", 2);
        F77_FUNC(dsortr, DSORTR)(wprime, &c__1, &iwork[9], &bounds[1], &ritz[1]);

        i__1 = iwork[9];
        for (j = 1; j <= i__1; ++j)
        {
            d__2 = eps23;
            d__3 = std::abs(ritz[j]);
            temp = (d__2 > d__3) ? d__2 : d__3;
            bounds[j] *= temp;
        }

        if (!std::strncmp(which, "BE", 2))
        {

            std::strncpy(wprime, "LA", 2);
            F77_FUNC(dsortr, DSORTR)(wprime, &c__1, &iwork[8], &ritz[1], &bounds[1]);
        }
        else
        {
            F77_FUNC(dsortr, DSORTR)(which, &c__1, &iwork[8], &ritz[1], &bounds[1]);
        }

        h__[h_dim1 + 1] = workd[*n * 3 + 1];


        if (iwork[6] > *mxiter && iwork[8] < *nev)
        {
            *info = 1;
        }
        if (*np == 0 && iwork[8] < iwork[9])
        {
            *info = 2;
        }

        *np = iwork[8];
        goto L1100;
    }
    else if (iwork[8] < *nev && *ishift == 1)
    {
        nevbef = *nev;
        i__1 = iwork[8], i__2 = *np / 2;
        *nev += (i__1 < i__2) ? i__1 : i__2;
        if (*nev == 1 && iwork[7] >= 6)
        {
            *nev = iwork[7] / 2;
        }
        else if (*nev == 1 && iwork[7] > 2)
        {
            *nev = 2;
        }
        *np = iwork[7] - *nev;


        if (nevbef < *nev)
        {
            F77_FUNC(dsgets, DSGETS)(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);
        }
    }


    if (*ishift == 0)
    {

        iwork[5] = 1;
        *ido     = 3;
        goto L9000;
    }

L50:

    iwork[5] = 0;

    if (*ishift == 0)
    {
        F77_FUNC(dcopy, DCOPY)(np, &workl[1], &c__1, &ritz[1], &c__1);
    }

    F77_FUNC(dsapps, DSAPPS)
    (n, nev, np, &ritz[1], &v[v_offset], ldv, &h__[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workd[1]);

    iwork[1] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(dcopy, DCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L100:

    if (*bmat == 'G')
    {
        workd[*n * 3 + 1] = F77_FUNC(ddot, DDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        workd[*n * 3 + 1] = std::sqrt(std::abs(workd[*n * 3 + 1]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 1] = F77_FUNC(dnrm2, DNRM2)(n, &resid[1], &c__1);
    }
    iwork[1] = 0;

    goto L1000;

L1100:

    *mxiter = iwork[6];
    *nev    = iwork[8];

L1200:
    *ido = 99;

L9000:
    return;
}


void F77_FUNC(dsaupd, DSAUPD)(int*        ido,
                              const char* bmat,
                              int*        n,
                              const char* which,
                              int*        nev,
                              double*     tol,
                              double*     resid,
                              int*        ncv,
                              double*     v,
                              int*        ldv,
                              int*        iparam,
                              int*        ipntr,
                              double*     workd,
                              int*        iwork,
                              double*     workl,
                              int*        lworkl,
                              int*        info)
{
    int v_dim1, v_offset, i__1, i__2;
    int j;

    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --iwork;
    --workl;

    if (*ido == 0)
    {


        iwork[2]  = 0;
        iwork[5]  = iparam[1];
        iwork[10] = iparam[3];
        iwork[12] = iparam[4];

        iwork[6]  = 1;
        iwork[11] = iparam[7];


        if (*n <= 0)
        {
            iwork[2] = -1;
        }
        else if (*nev <= 0)
        {
            iwork[2] = -2;
        }
        else if (*ncv <= *nev || *ncv > *n)
        {
            iwork[2] = -3;
        }


        iwork[15] = *ncv - *nev;

        if (iwork[10] <= 0)
        {
            iwork[2] = -4;
        }
        if (std::strncmp(which, "LM", 2) && std::strncmp(which, "SM", 2) && std::strncmp(which, "LA", 2)
            && std::strncmp(which, "SA", 2) && std::strncmp(which, "BE", 2))
        {
            iwork[2] = -5;
        }
        if (*bmat != 'I' && *bmat != 'G')
        {
            iwork[2] = -6;
        }

        i__1 = *ncv;
        if (*lworkl < i__1 * i__1 + (*ncv << 3))
        {
            iwork[2] = -7;
        }
        if (iwork[11] < 1 || iwork[11] > 5)
        {
            iwork[2] = -10;
        }
        else if (iwork[11] == 1 && *bmat == 'G')
        {
            iwork[2] = -11;
        }
        else if (iwork[5] < 0 || iwork[5] > 1)
        {
            iwork[2] = -12;
        }
        else if (*nev == 1 && !std::strncmp(which, "BE", 2))
        {
            iwork[2] = -13;
        }

        if (iwork[2] != 0)
        {
            *info = iwork[2];
            *ido  = 99;
            goto L9000;
        }

        if (iwork[12] <= 0)
        {
            iwork[12] = 1;
        }
        if (*tol <= 0.)
        {
            *tol = GMX_DOUBLE_EPS;
        }

        iwork[15] = *ncv - *nev;
        iwork[13] = *nev;
        i__2      = *ncv;
        i__1      = i__2 * i__2 + (*ncv << 3);
        for (j = 1; j <= i__1; ++j)
        {
            workl[j] = 0.;
        }

        iwork[8]  = *ncv;
        iwork[9]  = *ncv;
        iwork[3]  = 1;
        iwork[16] = iwork[3] + (iwork[8] << 1);
        iwork[1]  = iwork[16] + *ncv;
        iwork[4]  = iwork[1] + *ncv;
        i__1      = *ncv;
        iwork[7]  = iwork[4] + i__1 * i__1;
        iwork[14] = iwork[7] + *ncv * 3;

        ipntr[4]  = iwork[14];
        ipntr[5]  = iwork[3];
        ipntr[6]  = iwork[16];
        ipntr[7]  = iwork[1];
        ipntr[11] = iwork[7];
    }

    F77_FUNC(dsaup2, DSAUP2)
    (ido,
     bmat,
     n,
     which,
     &iwork[13],
     &iwork[15],
     tol,
     &resid[1],
     &iwork[11],
     &iwork[6],
     &iwork[5],
     &iwork[10],
     &v[v_offset],
     ldv,
     &workl[iwork[3]],
     &iwork[8],
     &workl[iwork[16]],
     &workl[iwork[1]],
     &workl[iwork[4]],
     &iwork[9],
     &workl[iwork[7]],
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido == 3)
    {
        iparam[8] = iwork[15];
    }
    if (*ido != 99)
    {
        goto L9000;
    }

    iparam[3] = iwork[10];
    iparam[5] = iwork[15];

    if (*info < 0)
    {
        goto L9000;
    }
    if (*info == 2)
    {
        *info = 3;
    }

L9000:

    return;
}


void F77_FUNC(dseupd, DSEUPD)(int*        rvec,
                              const char* howmny,
                              int*        select,
                              double*     d__,
                              double*     z__,
                              int*        ldz,
                              double*     sigma,
                              const char* bmat,
                              int*        n,
                              const char* which,
                              int*        nev,
                              double*     tol,
                              double*     resid,
                              int*        ncv,
                              double*     v,
                              int*        ldv,
                              int*        iparam,
                              int*        ipntr,
                              double*     workd,
                              double*     workl,
                              int*        lworkl,
                              int*        info)
{
    double c_b21  = 2 / 3.;
    int    c__1   = 1;
    double c_b102 = 1.;
    int    v_dim1, v_offset, z_dim1, z_offset, i__1;
    double d__1, d__2, d__3;

    int    j, k, ih, iq, iw, ibd, ihb, ihd, ldh, ilg, ldq, ism, irz;
    int    mode;
    double eps23;
    int    ierr;
    double temp;
    int    next;
    char   type__[6];
    int    ritz;
    int    reord;
    int    nconv;
    double rnorm;
    double bnorm2;
    double thres1 = 0, thres2 = 0;
    int    bounds;
    double tempbnd;
    int    leftptr, rghtptr;


    --workd;
    --resid;
    z_dim1   = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --d__;
    --select;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --workl;

    mode  = iparam[7];
    nconv = iparam[5];
    *info = 0;

    if (nconv == 0)
    {
        goto L9000;
    }
    ierr = 0;

    if (nconv <= 0)
    {
        ierr = -14;
    }
    if (*n <= 0)
    {
        ierr = -1;
    }
    if (*nev <= 0)
    {
        ierr = -2;
    }
    if (*ncv <= *nev || *ncv > *n)
    {
        ierr = -3;
    }
    if (std::strncmp(which, "LM", 2) && std::strncmp(which, "SM", 2) && std::strncmp(which, "LA", 2)
        && std::strncmp(which, "SA", 2) && std::strncmp(which, "BE", 2))
    {
        ierr = -5;
    }
    if (*bmat != 'I' && *bmat != 'G')
    {
        ierr = -6;
    }
    if (*howmny != 'A' && *howmny != 'P' && *howmny != 'S' && *rvec)
    {
        ierr = -15;
    }
    if (*rvec && *howmny == 'S')
    {
        ierr = -16;
    }
    i__1 = *ncv;
    if (*rvec && *lworkl < i__1 * i__1 + (*ncv << 3))
    {
        ierr = -7;
    }

    if (mode == 1 || mode == 2)
    {
        std::strncpy(type__, "REGULR", 6);
    }
    else if (mode == 3)
    {
        std::strncpy(type__, "SHIFTI", 6);
    }
    else if (mode == 4)
    {
        std::strncpy(type__, "BUCKLE", 6);
    }
    else if (mode == 5)
    {
        std::strncpy(type__, "CAYLEY", 6);
    }
    else
    {
        ierr = -10;
    }
    if (mode == 1 && *bmat == 'G')
    {
        ierr = -11;
    }
    if (*nev == 1 && !std::strncmp(which, "BE", 2))
    {
        ierr = -12;
    }

    if (ierr != 0)
    {
        *info = ierr;
        goto L9000;
    }

    ih        = ipntr[5];
    ritz      = ipntr[6];
    bounds    = ipntr[7];
    ldh       = *ncv;
    ldq       = *ncv;
    ihd       = bounds + ldh;
    ihb       = ihd + ldh;
    iq        = ihb + ldh;
    iw        = iq + ldh * *ncv;
    next      = iw + (*ncv << 1);
    ipntr[4]  = next;
    ipntr[8]  = ihd;
    ipntr[9]  = ihb;
    ipntr[10] = iq;

    irz = ipntr[11] + *ncv;
    ibd = irz + *ncv;


    eps23 = GMX_DOUBLE_EPS;
    eps23 = std::pow(eps23, c_b21);

    rnorm = workl[ih];
    if (*bmat == 'I')
    {
        bnorm2 = rnorm;
    }
    else if (*bmat == 'G')
    {
        bnorm2 = F77_FUNC(dnrm2, DNRM2)(n, &workd[1], &c__1);
    }

    if (*rvec)
    {

        if (!std::strncmp(which, "LM", 2) || !std::strncmp(which, "SM", 2)
            || !std::strncmp(which, "LA", 2) || !std::strncmp(which, "SA", 2))
        {
        }
        else if (!std::strncmp(which, "BE", 2))
        {


            ism = (*nev > nconv) ? *nev : nconv;
            ism /= 2;
            ilg    = ism + 1;
            thres1 = workl[ism];
            thres2 = workl[ilg];
        }

        reord = 0;
        i__1  = *ncv - 1;
        for (j = 0; j <= i__1; ++j)
        {
            select[j + 1] = 0;
            if (!std::strncmp(which, "LM", 2))
            {
                if (std::abs(workl[irz + j]) >= std::abs(thres1))
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "SM", 2))
            {
                if (std::abs(workl[irz + j]) <= std::abs(thres1))
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "LA", 2))
            {
                if (workl[irz + j] >= thres1)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "SA", 2))
            {
                if (workl[irz + j] <= thres1)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "BE", 2))
            {
                if (workl[irz + j] <= thres1 || workl[irz + j] >= thres2)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            if (j + 1 > nconv)
            {
                reord = select[j + 1] || reord;
            }
        }

        i__1 = *ncv - 1;
        F77_FUNC(dcopy, DCOPY)(&i__1, &workl[ih + 1], &c__1, &workl[ihb], &c__1);
        F77_FUNC(dcopy, DCOPY)(ncv, &workl[ih + ldh], &c__1, &workl[ihd], &c__1);

        F77_FUNC(dsteqr, DSTEQR)
        ("Identity", ncv, &workl[ihd], &workl[ihb], &workl[iq], &ldq, &workl[iw], &ierr);

        if (ierr != 0)
        {
            *info = -8;
            goto L9000;
        }


        if (reord)
        {

            leftptr = 1;
            rghtptr = *ncv;

            if (*ncv == 1)
            {
                goto L30;
            }

        L20:
            if (select[leftptr])
            {

                ++leftptr;
            }
            else if (!select[rghtptr])
            {

                --rghtptr;
            }
            else
            {

                temp                     = workl[ihd + leftptr - 1];
                workl[ihd + leftptr - 1] = workl[ihd + rghtptr - 1];
                workl[ihd + rghtptr - 1] = temp;
                F77_FUNC(dcopy, DCOPY)
                (ncv, &workl[iq + *ncv * (leftptr - 1)], &c__1, &workl[iw], &c__1);
                F77_FUNC(dcopy, DCOPY)
                (ncv, &workl[iq + *ncv * (rghtptr - 1)], &c__1, &workl[iq + *ncv * (leftptr - 1)], &c__1);
                F77_FUNC(dcopy, DCOPY)
                (ncv, &workl[iw], &c__1, &workl[iq + *ncv * (rghtptr - 1)], &c__1);
                ++leftptr;
                --rghtptr;
            }

            if (leftptr < rghtptr)
            {
                goto L20;
            }

        L30:;
        }

        F77_FUNC(dcopy, DCOPY)(&nconv, &workl[ihd], &c__1, &d__[1], &c__1);
    }
    else
    {

        F77_FUNC(dcopy, DCOPY)(&nconv, &workl[ritz], &c__1, &d__[1], &c__1);
        F77_FUNC(dcopy, DCOPY)(ncv, &workl[ritz], &c__1, &workl[ihd], &c__1);
    }
    if (!std::strncmp(type__, "REGULR", 6))
    {

        if (*rvec)
        {
            F77_FUNC(dsesrt, DSESRT)("LA", rvec, &nconv, &d__[1], ncv, &workl[iq], &ldq);
        }
        else
        {
            F77_FUNC(dcopy, DCOPY)(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
        }
    }
    else
    {

        F77_FUNC(dcopy, DCOPY)(ncv, &workl[ihd], &c__1, &workl[iw], &c__1);
        if (!std::strncmp(type__, "SHIFTI", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = 1. / workl[ihd + k - 1] + *sigma;
            }
        }
        else if (!std::strncmp(type__, "BUCKLE", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = *sigma * workl[ihd + k - 1] / (workl[ihd + k - 1] - 1.);
            }
        }
        else if (!std::strncmp(type__, "CAYLEY", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = *sigma * (workl[ihd + k - 1] + 1.) / (workl[ihd + k - 1] - 1.);
            }
        }

        F77_FUNC(dcopy, DCOPY)(&nconv, &workl[ihd], &c__1, &d__[1], &c__1);
        F77_FUNC(dsortr, DSORTR)("LA", &c__1, &nconv, &workl[ihd], &workl[iw]);
        if (*rvec)
        {
            F77_FUNC(dsesrt, DSESRT)("LA", rvec, &nconv, &d__[1], ncv, &workl[iq], &ldq);
        }
        else
        {
            F77_FUNC(dcopy, DCOPY)(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
            d__1 = bnorm2 / rnorm;
            F77_FUNC(dscal, DSCAL)(ncv, &d__1, &workl[ihb], &c__1);
            F77_FUNC(dsortr, DSORTR)("LA", &c__1, &nconv, &d__[1], &workl[ihb]);
        }
    }

    if (*rvec && *howmny == 'A')
    {

        F77_FUNC(dgeqr2, DGEQR2)
        (ncv, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb], &ierr);

        F77_FUNC(dorm2r, DORM2R)
        ("Right",
         "Notranspose",
         n,
         ncv,
         &nconv,
         &workl[iq],
         &ldq,
         &workl[iw + *ncv],
         &v[v_offset],
         ldv,
         &workd[*n + 1],
         &ierr);
        F77_FUNC(dlacpy, DLACPY)("All", n, &nconv, &v[v_offset], ldv, &z__[z_offset], ldz);

        i__1 = *ncv - 1;
        for (j = 1; j <= i__1; ++j)
        {
            workl[ihb + j - 1] = 0.;
        }
        workl[ihb + *ncv - 1] = 1.;
        F77_FUNC(dorm2r, DORM2R)
        ("Left", "Transpose", ncv, &c__1, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb], ncv, &temp, &ierr);
    }
    else if (*rvec && *howmny == 'S')
    {
    }

    if (!std::strncmp(type__, "REGULR", 6) && *rvec)
    {

        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            workl[ihb + j - 1] = rnorm * std::abs(workl[ihb + j - 1]);
        }
    }
    else if (std::strncmp(type__, "REGULR", 6) && *rvec)
    {

        F77_FUNC(dscal, DSCAL)(ncv, &bnorm2, &workl[ihb], &c__1);
        if (!std::strncmp(type__, "SHIFTI", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                d__2               = workl[iw + k - 1];
                workl[ihb + k - 1] = std::abs(workl[ihb + k - 1]) / (d__2 * d__2);
            }
        }
        else if (!std::strncmp(type__, "BUCKLE", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                d__2               = workl[iw + k - 1] - 1.;
                workl[ihb + k - 1] = *sigma * std::abs(workl[ihb + k - 1]) / (d__2 * d__2);
            }
        }
        else if (!std::strncmp(type__, "CAYLEY", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihb + k - 1] =
                        std::abs(workl[ihb + k - 1] / workl[iw + k - 1] * (workl[iw + k - 1] - 1.));
            }
        }
    }

    if (*rvec && (!std::strncmp(type__, "SHIFTI", 6) || !std::strncmp(type__, "CAYLEY", 6)))
    {

        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iq + k * ldq + *ncv - 1] / workl[iw + k];
        }
    }
    else if (*rvec && !std::strncmp(type__, "BUCKLE", 6))
    {

        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iq + k * ldq + *ncv - 1] / (workl[iw + k] - 1.);
        }
    }

    if (std::strncmp(type__, "REGULR", 6))
    {
        F77_FUNC(dger, DGER)
        (n, &nconv, &c_b102, &resid[1], &c__1, &workl[iw], &c__1, &z__[z_offset], ldz);
    }

L9000:

    return;
}


/* Selected single precision arpack routines */


static void F77_FUNC(sstqrb, SSTQRB)(int* n, float* d__, float* e, float* z__, float* work, int* info)
{
    int   i__1, i__2;
    float d__1, d__2;
    int   c__0  = 0;
    int   c__1  = 1;
    float c_b31 = 1.;

    float b, c__, f, g;
    int   i__, j, k, l, m;
    float p, r__, s;
    int   l1, ii, mm, lm1, mm1, nm1;
    float rt1, rt2, eps;
    int   lsv;
    float tst, eps2;
    int   lend, jtot, lendm1, lendp1, iscale;

    int   lendsv, nmaxit, icompz;
    float ssfmax, ssfmin, safmin, minval, safmax, anorm;


    --work;
    --z__;
    --e;
    --d__;

    *info = 0;

    icompz = 2;

    if (*n == 0)
    {
        return;
    }

    if (*n == 1)
    {
        z__[1] = 1.;
        return;
    }

    eps = GMX_FLOAT_EPS;

    d__1   = eps;
    eps2   = d__1 * d__1;
    minval = GMX_FLOAT_MIN;
    safmin = minval / GMX_FLOAT_EPS;
    safmax = 1. / safmin;
    ssfmax = std::sqrt(safmax) / 3.;
    ssfmin = std::sqrt(safmin) / eps2;

    i__1 = *n - 1;
    for (j = 1; j <= i__1; ++j)
    {
        z__[j] = 0.;
    }
    z__[*n] = 1.;

    nmaxit = *n * 30;
    jtot   = 0;

    l1  = 1;
    nm1 = *n - 1;

L10:
    if (l1 > *n)
    {
        goto L160;
    }
    if (l1 > 1)
    {
        e[l1 - 1] = 0.;
    }
    if (l1 <= nm1)
    {
        i__1 = nm1;
        for (m = l1; m <= i__1; ++m)
        {
            tst = std::abs(e[m]);
            if (tst == 0.)
            {
                goto L30;
            }
            if (tst <= std::sqrt(std::abs(d__[m])) * std::sqrt(std::abs(d__[m + 1])) * eps)
            {
                e[m] = 0.;
                goto L30;
            }
        }
    }
    m = *n;

L30:
    l      = l1;
    lsv    = l;
    lend   = m;
    lendsv = lend;
    l1     = m + 1;
    if (lend == l)
    {
        goto L10;
    }

    i__1   = lend - l + 1;
    anorm  = F77_FUNC(slanst, SLANST)("i", &i__1, &d__[l], &e[l]);
    iscale = 0;
    if (anorm == 0.)
    {
        goto L10;
    }
    if (anorm > ssfmax)
    {
        iscale = 1;
        i__1   = lend - l + 1;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &d__[l], n, info);
        i__1 = lend - l;
        F77_FUNC(slascl, SLASCL)("g", &c__0, &c__0, &anorm, &ssfmax, &i__1, &c__1, &e[l], n, info);
    }
    else if (anorm < ssfmin)
    {
        iscale = 2;
        i__1   = lend - l + 1;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &d__[l], n, info);
        i__1 = lend - l;
        F77_FUNC(slascl, SLASCL)("g", &c__0, &c__0, &anorm, &ssfmin, &i__1, &c__1, &e[l], n, info);
    }

    if (std::abs(d__[lend]) < std::abs(d__[l]))
    {
        lend = lsv;
        l    = lendsv;
    }

    if (lend > l)
    {

    L40:
        if (l != lend)
        {
            lendm1 = lend - 1;
            i__1   = lendm1;
            for (m = l; m <= i__1; ++m)
            {
                d__2 = std::abs(e[m]);
                tst  = d__2 * d__2;
                if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m + 1]) + safmin)
                {
                    goto L60;
                }
            }
        }

        m = lend;

    L60:
        if (m < lend)
        {
            e[m] = 0.;
        }
        p = d__[l];
        if (m == l)
        {
            goto L80;
        }

        if (m == l + 1)
        {
            if (icompz > 0)
            {
                F77_FUNC(slaev2, SLAEV2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2, &c__, &s);
                work[l]          = c__;
                work[*n - 1 + l] = s;

                tst        = z__[l + 1];
                z__[l + 1] = c__ * tst - s * z__[l];
                z__[l]     = s * tst + c__ * z__[l];
            }
            else
            {
                F77_FUNC(slae2, SLAE2)(&d__[l], &e[l], &d__[l + 1], &rt1, &rt2);
            }
            d__[l]     = rt1;
            d__[l + 1] = rt2;
            e[l]       = 0.;
            l += 2;
            if (l <= lend)
            {
                goto L40;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;

        g   = (d__[l + 1] - p) / (e[l] * 2.);
        r__ = F77_FUNC(slapy2, SLAPY2)(&g, &c_b31);
        g   = d__[m] - p + e[l] / (g + ((g > 0) ? r__ : -r__));

        s   = 1.;
        c__ = 1.;
        p   = 0.;

        mm1  = m - 1;
        i__1 = l;
        for (i__ = mm1; i__ >= i__1; --i__)
        {
            f = s * e[i__];
            b = c__ * e[i__];
            F77_FUNC(slartg, SLARTG)(&g, &f, &c__, &s, &r__);
            if (i__ != m - 1)
            {
                e[i__ + 1] = r__;
            }
            g            = d__[i__ + 1] - p;
            r__          = (d__[i__] - g) * s + c__ * 2. * b;
            p            = s * r__;
            d__[i__ + 1] = g + p;
            g            = c__ * r__ - b;

            if (icompz > 0)
            {
                work[i__]          = c__;
                work[*n - 1 + i__] = -s;
            }
        }

        if (icompz > 0)
        {
            mm = m - l + 1;

            F77_FUNC(slasr, SLASR)
            ("r", "v", "b", &c__1, &mm, &work[l], &work[*n - 1 + l], &z__[l], &c__1);
        }

        d__[l] -= p;
        e[l] = g;
        goto L40;

    L80:
        d__[l] = p;

        ++l;
        if (l <= lend)
        {
            goto L40;
        }
        goto L140;
    }
    else
    {

    L90:
        if (l != lend)
        {
            lendp1 = lend + 1;
            i__1   = lendp1;
            for (m = l; m >= i__1; --m)
            {
                d__2 = std::abs(e[m - 1]);
                tst  = d__2 * d__2;
                if (tst <= eps2 * std::abs(d__[m]) * std::abs(d__[m - 1]) + safmin)
                {
                    goto L110;
                }
            }
        }

        m = lend;

    L110:
        if (m > lend)
        {
            e[m - 1] = 0.;
        }
        p = d__[l];
        if (m == l)
        {
            goto L130;
        }

        if (m == l - 1)
        {
            if (icompz > 0)
            {
                F77_FUNC(slaev2, SLAEV2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2, &c__, &s);

                tst        = z__[l];
                z__[l]     = c__ * tst - s * z__[l - 1];
                z__[l - 1] = s * tst + c__ * z__[l - 1];
            }
            else
            {
                F77_FUNC(slae2, SLAE2)(&d__[l - 1], &e[l - 1], &d__[l], &rt1, &rt2);
            }
            d__[l - 1] = rt1;
            d__[l]     = rt2;
            e[l - 1]   = 0.;
            l += -2;
            if (l >= lend)
            {
                goto L90;
            }
            goto L140;
        }

        if (jtot == nmaxit)
        {
            goto L140;
        }
        ++jtot;


        g   = (d__[l - 1] - p) / (e[l - 1] * 2.);
        r__ = F77_FUNC(slapy2, SLAPY2)(&g, &c_b31);
        g   = d__[m] - p + e[l - 1] / (g + ((g > 0) ? r__ : -r__));

        s   = 1.;
        c__ = 1.;
        p   = 0.;

        lm1  = l - 1;
        i__1 = lm1;
        for (i__ = m; i__ <= i__1; ++i__)
        {
            f = s * e[i__];
            b = c__ * e[i__];
            F77_FUNC(slartg, SLARTG)(&g, &f, &c__, &s, &r__);
            if (i__ != m)
            {
                e[i__ - 1] = r__;
            }
            g        = d__[i__] - p;
            r__      = (d__[i__ + 1] - g) * s + c__ * 2. * b;
            p        = s * r__;
            d__[i__] = g + p;
            g        = c__ * r__ - b;

            if (icompz > 0)
            {
                work[i__]          = c__;
                work[*n - 1 + i__] = s;
            }
        }

        if (icompz > 0)
        {
            mm = l - m + 1;

            F77_FUNC(slasr, SLASR)
            ("r", "v", "f", &c__1, &mm, &work[m], &work[*n - 1 + m], &z__[m], &c__1);
        }

        d__[l] -= p;
        e[lm1] = g;
        goto L90;

    L130:
        d__[l] = p;

        --l;
        if (l >= lend)
        {
            goto L90;
        }
        goto L140;
    }

L140:
    if (iscale == 1)
    {
        i__1 = lendsv - lsv + 1;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &d__[lsv], n, info);
        i__1 = lendsv - lsv;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &ssfmax, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }
    else if (iscale == 2)
    {
        i__1 = lendsv - lsv + 1;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &d__[lsv], n, info);
        i__1 = lendsv - lsv;
        F77_FUNC(slascl, SLASCL)
        ("g", &c__0, &c__0, &ssfmin, &anorm, &i__1, &c__1, &e[lsv], n, info);
    }

    if (jtot < nmaxit)
    {
        goto L10;
    }
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        if (e[i__] != 0.)
        {
            ++(*info);
        }
    }
    goto L190;

L160:
    if (icompz == 0)
    {

        F77_FUNC(slasrt, SLASRT)("i", n, &d__[1], info);
    }
    else
    {

        i__1 = *n;
        for (ii = 2; ii <= i__1; ++ii)
        {
            i__  = ii - 1;
            k    = i__;
            p    = d__[i__];
            i__2 = *n;
            for (j = ii; j <= i__2; ++j)
            {
                if (d__[j] < p)
                {
                    k = j;
                    p = d__[j];
                }
            }
            if (k != i__)
            {
                d__[k]   = d__[i__];
                d__[i__] = p;

                p        = z__[k];
                z__[k]   = z__[i__];
                z__[i__] = p;
            }
        }
    }

L190:
    return;
}

static void F77_FUNC(sgetv0, SGETV0)(int*        ido,
                                     const char* bmat,
                                     int gmx_unused* itry,
                                     int*            initv,
                                     int*            n,
                                     int*            j,
                                     float*          v,
                                     int*            ldv,
                                     float*          resid,
                                     float*          rnorm,
                                     int*            ipntr,
                                     float*          workd,
                                     int*            iwork,
                                     int*            ierr)
{
    int   c__1  = 1;
    float c_b22 = 1.;
    float c_b24 = 0.;
    float c_b27 = -1.;
    int   v_dim1, v_offset, i__1;

    int jj;
    int idist;

    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --ipntr;
    --iwork;

    if (*ido == 0)
    {

        *ierr    = 0;
        iwork[7] = 0;
        iwork[5] = 0;
        iwork[6] = 0;

        if (!(*initv))
        {
            idist = 2;
            F77_FUNC(slarnv, SLARNV)(&idist, &iwork[1], n, &resid[1]);
        }

        if (*bmat == 'G')
        {
            ipntr[1] = 1;
            ipntr[2] = *n + 1;
            F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
            *ido = -1;
            goto L9000;
        }
    }

    if (iwork[5] == 1)
    {
        goto L20;
    }

    if (iwork[6] == 1)
    {
        goto L40;
    }

    iwork[5] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(scopy, SCOPY)(n, &workd[*n + 1], &c__1, &resid[1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L20:


    iwork[5] = 0;
    if (*bmat == 'G')
    {
        workd[*n * 3 + 4] = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        workd[*n * 3 + 4] = std::sqrt(std::abs(workd[*n * 3 + 4]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 4] = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }
    *rnorm = workd[*n * 3 + 4];

    if (*j == 1)
    {
        goto L50;
    }
    iwork[6] = 1;
L30:

    i__1 = *j - 1;
    F77_FUNC(sgemv, SGEMV)
    ("T", n, &i__1, &c_b22, &v[v_offset], ldv, &workd[1], &c__1, &c_b24, &workd[*n + 1], &c__1);
    i__1 = *j - 1;
    F77_FUNC(sgemv, SGEMV)
    ("N", n, &i__1, &c_b27, &v[v_offset], ldv, &workd[*n + 1], &c__1, &c_b22, &resid[1], &c__1);

    if (*bmat == 'G')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;
        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L40:

    if (*bmat == 'G')
    {
        *rnorm = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        *rnorm = std::sqrt(std::abs(*rnorm));
    }
    else if (*bmat == 'I')
    {
        *rnorm = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }

    if (*rnorm > workd[*n * 3 + 4] * .717F)
    {
        goto L50;
    }

    ++iwork[7];
    if (iwork[7] <= 1)
    {

        workd[*n * 3 + 4] = *rnorm;
        goto L30;
    }
    else
    {

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.;
        }
        *rnorm = 0.;
        *ierr  = -1;
    }

L50:

    *ido = 99;

L9000:
    return;
}


static void F77_FUNC(ssapps, SSAPPS)(int*   n,
                                     int*   kev,
                                     int*   np,
                                     float* shift,
                                     float* v,
                                     int*   ldv,
                                     float* h__,
                                     int*   ldh,
                                     float* resid,
                                     float* q,
                                     int*   ldq,
                                     float* workd)
{
    float c_b4  = 0.;
    float c_b5  = 1.;
    float c_b14 = -1.;
    int   c__1  = 1;
    int   h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, i__3, i__4;
    float c__, f, g;
    int   i__, j;
    float r__, s, a1, a2, a3, a4;
    int   jj;
    float big;
    int   iend, itop;
    float epsmch;
    int   istart, kplusp;

    --workd;
    --resid;
    --shift;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1   = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;

    epsmch = GMX_FLOAT_EPS;
    itop   = 1;


    kplusp = *kev + *np;

    F77_FUNC(slaset, SLASET)("All", &kplusp, &kplusp, &c_b4, &c_b5, &q[q_offset], ldq);

    if (*np == 0)
    {
        goto L9000;
    }

    i__1 = *np;
    for (jj = 1; jj <= i__1; ++jj)
    {

        istart = itop;

    L20:

        i__2 = kplusp - 1;
        for (i__ = istart; i__ <= i__2; ++i__)
        {
            big = std::abs(h__[i__ + (h_dim1 * 2)]) + std::abs(h__[i__ + 1 + (h_dim1 * 2)]);
            if (h__[i__ + 1 + h_dim1] <= epsmch * big)
            {
                h__[i__ + 1 + h_dim1] = 0.;
                iend                  = i__;
                goto L40;
            }
        }
        iend = kplusp;
    L40:

        if (istart < iend)
        {

            f = h__[istart + (h_dim1 << 1)] - shift[jj];
            g = h__[istart + 1 + h_dim1];
            F77_FUNC(slartg, SLARTG)(&f, &g, &c__, &s, &r__);

            a1 = c__ * h__[istart + (h_dim1 << 1)] + s * h__[istart + 1 + h_dim1];
            a2 = c__ * h__[istart + 1 + h_dim1] + s * h__[istart + 1 + (h_dim1 << 1)];
            a4 = c__ * h__[istart + 1 + (h_dim1 << 1)] - s * h__[istart + 1 + h_dim1];
            a3 = c__ * h__[istart + 1 + h_dim1] - s * h__[istart + (h_dim1 << 1)];
            h__[istart + (h_dim1 << 1)]     = c__ * a1 + s * a2;
            h__[istart + 1 + (h_dim1 << 1)] = c__ * a4 - s * a3;
            h__[istart + 1 + h_dim1]        = c__ * a3 + s * a4;

            i__3 = istart + jj;
            i__2 = (i__3 < kplusp) ? i__3 : kplusp;
            for (j = 1; j <= i__2; ++j)
            {
                a1 = c__ * q[j + istart * q_dim1] + s * q[j + (istart + 1) * q_dim1];
                q[j + (istart + 1) * q_dim1] =
                        -s * q[j + istart * q_dim1] + c__ * q[j + (istart + 1) * q_dim1];
                q[j + istart * q_dim1] = a1;
            }

            i__2 = iend - 1;
            for (i__ = istart + 1; i__ <= i__2; ++i__)
            {

                f = h__[i__ + h_dim1];
                g = s * h__[i__ + 1 + h_dim1];

                h__[i__ + 1 + h_dim1] = c__ * h__[i__ + 1 + h_dim1];
                F77_FUNC(slartg, SLARTG)(&f, &g, &c__, &s, &r__);

                if (r__ < 0.)
                {
                    r__ = -r__;
                    c__ = -c__;
                    s   = -s;
                }

                h__[i__ + h_dim1] = r__;

                a1 = c__ * h__[i__ + (h_dim1 << 1)] + s * h__[i__ + 1 + h_dim1];
                a2 = c__ * h__[i__ + 1 + h_dim1] + s * h__[i__ + 1 + (h_dim1 << 1)];
                a3 = c__ * h__[i__ + 1 + h_dim1] - s * h__[i__ + (h_dim1 << 1)];
                a4 = c__ * h__[i__ + 1 + (h_dim1 << 1)] - s * h__[i__ + 1 + h_dim1];

                h__[i__ + (h_dim1 << 1)]     = c__ * a1 + s * a2;
                h__[i__ + 1 + (h_dim1 << 1)] = c__ * a4 - s * a3;
                h__[i__ + 1 + h_dim1]        = c__ * a3 + s * a4;

                i__4 = j + jj;
                i__3 = (i__4 < kplusp) ? i__4 : kplusp;
                for (j = 1; j <= i__3; ++j)
                {
                    a1 = c__ * q[j + i__ * q_dim1] + s * q[j + (i__ + 1) * q_dim1];
                    q[j + (i__ + 1) * q_dim1] = -s * q[j + i__ * q_dim1] + c__ * q[j + (i__ + 1) * q_dim1];
                    q[j + i__ * q_dim1] = a1;
                }
            }
        }

        istart = iend + 1;

        if (h__[iend + h_dim1] < 0.)
        {
            h__[iend + h_dim1] = -h__[iend + h_dim1];
            F77_FUNC(sscal, SSCAL)(&kplusp, &c_b14, &q[iend * q_dim1 + 1], &c__1);
        }

        if (iend < kplusp)
        {
            goto L20;
        }

        i__2 = kplusp - 1;
        for (i__ = itop; i__ <= i__2; ++i__)
        {
            if (h__[i__ + 1 + h_dim1] > 0.)
            {
                goto L90;
            }
            ++itop;
        }

    L90:;
    }

    i__1 = kplusp - 1;
    for (i__ = itop; i__ <= i__1; ++i__)
    {
        big = std::abs(h__[i__ + (h_dim1 * 2)]) + std::abs(h__[i__ + 1 + (h_dim1 * 2)]);
        if (h__[i__ + 1 + h_dim1] <= epsmch * big)
        {
            h__[i__ + 1 + h_dim1] = 0.;
        }
    }

    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(sgemv, SGEMV)
        ("N", n, &kplusp, &c_b5, &v[v_offset], ldv, &q[(*kev + 1) * q_dim1 + 1], &c__1, &c_b4, &workd[*n + 1], &c__1);
    }

    i__1 = *kev;
    for (i__ = 1; i__ <= i__1; ++i__)
    {
        i__2 = kplusp - i__ + 1;
        F77_FUNC(sgemv, SGEMV)
        ("N", n, &i__2, &c_b5, &v[v_offset], ldv, &q[(*kev - i__ + 1) * q_dim1 + 1], &c__1, &c_b4, &workd[1], &c__1);
        F77_FUNC(scopy, SCOPY)(n, &workd[1], &c__1, &v[(kplusp - i__ + 1) * v_dim1 + 1], &c__1);
    }

    F77_FUNC(slacpy, SLACPY)("All", n, kev, &v[(*np + 1) * v_dim1 + 1], ldv, &v[v_offset], ldv);

    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(scopy, SCOPY)(n, &workd[*n + 1], &c__1, &v[(*kev + 1) * v_dim1 + 1], &c__1);
    }

    F77_FUNC(sscal, SSCAL)(n, &q[kplusp + *kev * q_dim1], &resid[1], &c__1);
    if (h__[*kev + 1 + h_dim1] > 0.)
    {
        F77_FUNC(saxpy, SAXPY)
        (n, &h__[*kev + 1 + h_dim1], &v[(*kev + 1) * v_dim1 + 1], &c__1, &resid[1], &c__1);
    }


L9000:
    return;
}


static void F77_FUNC(ssortr, SSORTR)(const char* which, int* apply, int* n, float* x1, float* x2)
{
    int i__1;

    int   i__, j, igap;
    float temp;


    igap = *n / 2;

    if (!std::strncmp(which, "SA", 2))
    {

    L10:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L20:

            if (j < 0)
            {
                goto L30;
            }

            if (x1[j] < x1[j + igap])
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L30;
            }
            j -= igap;
            goto L20;
        L30:;
        }
        igap /= 2;
        goto L10;
    }
    else if (!std::strncmp(which, "SM", 2))
    {

    L40:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L50:

            if (j < 0)
            {
                goto L60;
            }

            if (std::abs(x1[j]) < std::abs(x1[j + igap]))
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L60;
            }
            j -= igap;
            goto L50;
        L60:;
        }
        igap /= 2;
        goto L40;
    }
    else if (!std::strncmp(which, "LA", 2))
    {

    L70:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L80:

            if (j < 0)
            {
                goto L90;
            }

            if (x1[j] > x1[j + igap])
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L90;
            }
            j -= igap;
            goto L80;
        L90:;
        }
        igap /= 2;
        goto L70;
    }
    else if (!std::strncmp(which, "LM", 2))
    {


    L100:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L110:

            if (j < 0)
            {
                goto L120;
            }

            if (std::abs(x1[j]) > std::abs(x1[j + igap]))
            {
                temp         = x1[j];
                x1[j]        = x1[j + igap];
                x1[j + igap] = temp;
                if (*apply)
                {
                    temp         = x2[j];
                    x2[j]        = x2[j + igap];
                    x2[j + igap] = temp;
                }
            }
            else
            {
                goto L120;
            }
            j -= igap;
            goto L110;
        L120:;
        }
        igap /= 2;
        goto L100;
    }

L9000:
    return;
}


static void F77_FUNC(ssesrt,
                     SSESRT)(const char* which, int* apply, int* n, float* x, int* na, float* a, int* lda)
{
    int a_dim1, a_offset, i__1;
    int c__1 = 1;

    int   i__, j, igap;
    float temp;

    a_dim1   = *lda;
    a_offset = 1 + a_dim1 * 0;
    a -= a_offset;

    igap = *n / 2;

    if (!std::strncmp(which, "SA", 2))
    {

    L10:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L20:

            if (j < 0)
            {
                goto L30;
            }

            if (x[j] < x[j + igap])
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(sswap, SSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L30;
            }
            j -= igap;
            goto L20;
        L30:;
        }
        igap /= 2;
        goto L10;
    }
    else if (!std::strncmp(which, "SM", 2))
    {

    L40:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L50:

            if (j < 0)
            {
                goto L60;
            }

            if (std::abs(x[j]) < std::abs(x[j + igap]))
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(sswap, SSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L60;
            }
            j -= igap;
            goto L50;
        L60:;
        }
        igap /= 2;
        goto L40;
    }
    else if (!std::strncmp(which, "LA", 2))
    {

    L70:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L80:

            if (j < 0)
            {
                goto L90;
            }

            if (x[j] > x[j + igap])
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(sswap, SSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L90;
            }
            j -= igap;
            goto L80;
        L90:;
        }
        igap /= 2;
        goto L70;
    }
    else if (!std::strncmp(which, "LM", 2))
    {

    L100:
        if (igap == 0)
        {
            goto L9000;
        }
        i__1 = *n - 1;
        for (i__ = igap; i__ <= i__1; ++i__)
        {
            j = i__ - igap;
        L110:

            if (j < 0)
            {
                goto L120;
            }

            if (std::abs(x[j]) > std::abs(x[j + igap]))
            {
                temp        = x[j];
                x[j]        = x[j + igap];
                x[j + igap] = temp;
                if (*apply)
                {
                    F77_FUNC(sswap, SSWAP)
                    (na, &a[j * a_dim1 + 1], &c__1, &a[(j + igap) * a_dim1 + 1], &c__1);
                }
            }
            else
            {
                goto L120;
            }
            j -= igap;
            goto L110;
        L120:;
        }
        igap /= 2;
        goto L100;
    }

L9000:
    return;
}


static void F77_FUNC(ssgets,
                     SSGETS)(int* ishift, const char* which, int* kev, int* np, float* ritz, float* bounds, float* shifts)
{
    int c__1 = 1;
    int i__1, i__2;
    int kevd2;

    --shifts;
    --bounds;
    --ritz;

    if (!std::strncmp(which, "BE", 2))
    {
        i__1 = *kev + *np;
        F77_FUNC(ssortr, SSORTR)("LA", &c__1, &i__1, &ritz[1], &bounds[1]);
        kevd2 = *kev / 2;
        if (*kev > 1)
        {
            i__1 = (kevd2 < *np) ? kevd2 : *np;
            i__2 = (kevd2 > *np) ? kevd2 : *np;
            F77_FUNC(sswap, SSWAP)(&i__1, &ritz[1], &c__1, &ritz[i__2 + 1], &c__1);
            i__1 = (kevd2 < *np) ? kevd2 : *np;
            i__2 = (kevd2 > *np) ? kevd2 : *np;
            F77_FUNC(sswap, SSWAP)(&i__1, &bounds[1], &c__1, &bounds[i__2 + 1], &c__1);
        }
    }
    else
    {
        i__1 = *kev + *np;
        F77_FUNC(ssortr, SSORTR)(which, &c__1, &i__1, &ritz[1], &bounds[1]);
    }

    if (*ishift == 1 && *np > 0)
    {

        F77_FUNC(ssortr, SSORTR)("SM", &c__1, np, &bounds[1], &ritz[1]);
        F77_FUNC(scopy, SCOPY)(np, &ritz[1], &c__1, &shifts[1], &c__1);
    }


    return;
}


static void F77_FUNC(ssconv, SSCONV)(int* n, float* ritz, float* bounds, float* tol, int* nconv)
{
    float c_b3 = 2 / 3.;
    int   i__1;
    float d__2, d__3;

    int   i__;
    float eps23, temp;

    --bounds;
    --ritz;

    eps23 = GMX_FLOAT_EPS;
    eps23 = std::pow(eps23, c_b3);

    *nconv = 0;
    i__1   = *n;
    for (i__ = 1; i__ <= i__1; ++i__)
    {

        d__2 = eps23;
        d__3 = std::abs(ritz[i__]);
        temp = (d__2 > d__3) ? d__2 : d__3;
        if (bounds[i__] <= *tol * temp)
        {
            ++(*nconv);
        }
    }

    return;
}


static void F77_FUNC(
        sseigt,
        SSEIGT)(float* rnorm, int* n, float* h__, int* ldh, float* eig, float* bounds, float* workl, int* ierr)
{
    int c__1 = 1;
    int h_dim1, h_offset, i__1;

    int k;


    --workl;
    --bounds;
    --eig;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;

    F77_FUNC(scopy, SCOPY)(n, &h__[(h_dim1 << 1) + 1], &c__1, &eig[1], &c__1);
    i__1 = *n - 1;
    F77_FUNC(scopy, SCOPY)(&i__1, &h__[h_dim1 + 2], &c__1, &workl[1], &c__1);
    F77_FUNC(sstqrb, SSTQRB)(n, &eig[1], &workl[1], &bounds[1], &workl[*n + 1], ierr);
    if (*ierr != 0)
    {
        goto L9000;
    }

    i__1 = *n;
    for (k = 1; k <= i__1; ++k)
    {
        bounds[k] = *rnorm * std::abs(bounds[k]);
    }


L9000:
    return;
}


static void F77_FUNC(ssaitr, SSAITR)(int*        ido,
                                     const char* bmat,
                                     int*        n,
                                     int*        k,
                                     int*        np,
                                     int*        mode,
                                     float*      resid,
                                     float*      rnorm,
                                     float*      v,
                                     int*        ldv,
                                     float*      h__,
                                     int*        ldh,
                                     int*        ipntr,
                                     float*      workd,
                                     int*        iwork,
                                     int*        info)
{

    int   c__0  = 0;
    int   c__1  = 1;
    float c_b18 = 1.;
    float c_b42 = 0.;
    float c_b50 = -1.;

    int   h_dim1, h_offset, v_dim1, v_offset, i__1;
    int   i__, jj;
    float temp1;
    int   infol;
    float safmin, minval;


    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --ipntr;
    --iwork;
    minval = GMX_FLOAT_MIN;
    safmin = minval / GMX_FLOAT_EPS;

    if (*ido == 0)
    {
        *info    = 0;
        iwork[5] = 0;
        iwork[6] = 0;
        iwork[4] = 0;
        iwork[2] = 0;
        iwork[3] = 0;

        iwork[12] = *k + 1;

        iwork[8]  = 1;
        iwork[9]  = iwork[8] + *n;
        iwork[10] = iwork[9] + *n;
    }

    if (iwork[5] == 1)
    {
        goto L50;
    }
    if (iwork[6] == 1)
    {
        goto L60;
    }
    if (iwork[2] == 1)
    {
        goto L70;
    }
    if (iwork[3] == 1)
    {
        goto L90;
    }
    if (iwork[4] == 1)
    {
        goto L30;
    }

L1000:


    if (*rnorm > 0.)
    {
        goto L40;
    }

    iwork[11] = 1;
L20:
    iwork[4] = 1;
    *ido     = 0;
L30:

    F77_FUNC(sgetv0, sgetv0)
    (ido,
     bmat,
     &iwork[11],
     &c__0,
     n,
     &iwork[12],
     &v[v_offset],
     ldv,
     &resid[1],
     rnorm,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     &iwork[7]);
    if (*ido != 99)
    {
        goto L9000;
    }
    if (iwork[7] < 0)
    {
        ++iwork[11];
        if (iwork[11] <= 3)
        {
            goto L20;
        }

        *info = iwork[12] - 1;
        *ido  = 99;
        goto L9000;
    }

L40:

    F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &v[iwork[12] * v_dim1 + 1], &c__1);
    if (*rnorm >= safmin)
    {
        temp1 = 1. / *rnorm;
        F77_FUNC(sscal, SSCAL)(n, &temp1, &v[iwork[12] * v_dim1 + 1], &c__1);
        F77_FUNC(sscal, SSCAL)(n, &temp1, &workd[iwork[8]], &c__1);
    }
    else
    {

        F77_FUNC(slascl, SLASCL)
        ("General", &i__, &i__, rnorm, &c_b18, n, &c__1, &v[iwork[12] * v_dim1 + 1], n, &infol);
        F77_FUNC(slascl, SLASCL)
        ("General", &i__, &i__, rnorm, &c_b18, n, &c__1, &workd[iwork[8]], n, &infol);
    }

    iwork[5] = 1;
    F77_FUNC(scopy, SCOPY)(n, &v[iwork[12] * v_dim1 + 1], &c__1, &workd[iwork[10]], &c__1);
    ipntr[1] = iwork[10];
    ipntr[2] = iwork[9];
    ipntr[3] = iwork[8];
    *ido     = 1;

    goto L9000;
L50:


    iwork[5] = 0;

    F77_FUNC(scopy, SCOPY)(n, &workd[iwork[9]], &c__1, &resid[1], &c__1);

    if (*mode == 2)
    {
        goto L65;
    }
    if (*bmat == 'G')
    {
        iwork[6] = 1;
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L60:

    iwork[6] = 0;

L65:
    if (*mode == 2)
    {

        workd[*n * 3 + 3] = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[iwork[10]], &c__1);
        workd[*n * 3 + 3] = std::sqrt(std::abs(workd[*n * 3 + 3]));
    }
    else if (*bmat == 'G')
    {
        workd[*n * 3 + 3] = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        workd[*n * 3 + 3] = std::sqrt(std::abs(workd[*n * 3 + 3]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 3] = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }

    if (*mode != 2)
    {
        F77_FUNC(sgemv, SGEMV)
        ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[8]], &c__1, &c_b42, &workd[iwork[9]], &c__1);
    }
    else
    {
        F77_FUNC(sgemv, SGEMV)
        ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[10]], &c__1, &c_b42, &workd[iwork[9]], &c__1);
    }

    F77_FUNC(sgemv, SGEMV)
    ("N", n, &iwork[12], &c_b50, &v[v_offset], ldv, &workd[iwork[9]], &c__1, &c_b18, &resid[1], &c__1);

    h__[iwork[12] + (h_dim1 << 1)] = workd[iwork[9] + iwork[12] - 1];
    if (iwork[12] == 1 || iwork[4] == 1)
    {
        h__[iwork[12] + h_dim1] = 0.;
    }
    else
    {
        h__[iwork[12] + h_dim1] = *rnorm;
    }

    iwork[2] = 1;
    iwork[1] = 0;

    if (*bmat == 'G')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[iwork[9]], &c__1);
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L70:

    iwork[2] = 0;

    if (*bmat == 'G')
    {
        *rnorm = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        *rnorm = std::sqrt(std::abs(*rnorm));
    }
    else if (*bmat == 'I')
    {
        *rnorm = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }

    if (*rnorm > workd[*n * 3 + 3] * .717F)
    {
        goto L100;
    }

L80:

    F77_FUNC(sgemv, SGEMV)
    ("T", n, &iwork[12], &c_b18, &v[v_offset], ldv, &workd[iwork[8]], &c__1, &c_b42, &workd[iwork[9]], &c__1);

    F77_FUNC(sgemv, SGEMV)
    ("N", n, &iwork[12], &c_b50, &v[v_offset], ldv, &workd[iwork[9]], &c__1, &c_b18, &resid[1], &c__1);

    if (iwork[12] == 1 || iwork[4] == 1)
    {
        h__[iwork[12] + h_dim1] = 0.;
    }
    h__[iwork[12] + (h_dim1 << 1)] += workd[iwork[9] + iwork[12] - 1];

    iwork[3] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[iwork[9]], &c__1);
        ipntr[1] = iwork[9];
        ipntr[2] = iwork[8];
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
    }
L90:


    if (*bmat == 'G')
    {
        workd[*n * 3 + 2] = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[iwork[8]], &c__1);
        workd[*n * 3 + 2] = std::sqrt(std::abs(workd[*n * 3 + 2]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 2] = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }


    if (workd[*n * 3 + 2] > *rnorm * .717F)
    {

        *rnorm = workd[*n * 3 + 2];
    }
    else
    {

        *rnorm = workd[*n * 3 + 2];
        ++iwork[1];
        if (iwork[1] <= 1)
        {
            goto L80;
        }

        i__1 = *n;
        for (jj = 1; jj <= i__1; ++jj)
        {
            resid[jj] = 0.;
        }
        *rnorm = 0.;
    }

L100:

    iwork[4] = 0;
    iwork[3] = 0;

    if (h__[iwork[12] + h_dim1] < 0.)
    {
        h__[iwork[12] + h_dim1] = -h__[iwork[12] + h_dim1];
        if (iwork[12] < *k + *np)
        {
            F77_FUNC(sscal, SSCAL)(n, &c_b50, &v[(iwork[12] + 1) * v_dim1 + 1], &c__1);
        }
        else
        {
            F77_FUNC(sscal, SSCAL)(n, &c_b50, &resid[1], &c__1);
        }
    }

    ++iwork[12];
    if (iwork[12] > *k + *np)
    {
        *ido = 99;


        goto L9000;
    }

    goto L1000;

L9000:
    return;
}


static void F77_FUNC(ssaup2, SSAUP2)(int*        ido,
                                     const char* bmat,
                                     int*        n,
                                     const char* which,
                                     int*        nev,
                                     int*        np,
                                     float*      tol,
                                     float*      resid,
                                     int*        mode,
                                     int gmx_unused* iupd,
                                     int*            ishift,
                                     int*            mxiter,
                                     float*          v,
                                     int*            ldv,
                                     float*          h__,
                                     int*            ldh,
                                     float*          ritz,
                                     float*          bounds,
                                     float*          q,
                                     int*            ldq,
                                     float*          workl,
                                     int*            ipntr,
                                     float*          workd,
                                     int*            iwork,
                                     int*            info)
{
    float c_b3 = 2 / 3.;
    int   c__1 = 1;
    int   c__0 = 0;

    int   h_dim1, h_offset, q_dim1, q_offset, v_dim1, v_offset, i__1, i__2, i__3;
    float d__2, d__3;
    int   j;
    float eps23;
    int   ierr;
    float temp;
    int   nevd2;
    int   nevm2;
    int   nevbef;
    char  wprime[2];
    int   nptemp;

    --workd;
    --resid;
    --workl;
    --bounds;
    --ritz;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    h_dim1   = *ldh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    q_dim1   = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --ipntr;
    --iwork;
    eps23 = GMX_FLOAT_EPS;
    eps23 = std::pow(eps23, c_b3);

    if (*ido == 0)
    {

        iwork[41] = 1;
        iwork[42] = 3;
        iwork[43] = 5;
        iwork[44] = 7;

        iwork[9]  = *nev;
        iwork[10] = *np;

        iwork[7] = iwork[9] + iwork[10];
        iwork[8] = 0;
        iwork[6] = 0;

        iwork[2] = 1;
        iwork[4] = 0;
        iwork[5] = 0;
        iwork[1] = 0;

        if (*info != 0)
        {

            iwork[3] = 1;
            *info    = 0;
        }
        else
        {
            iwork[3] = 0;
        }
    }

    if (iwork[2] == 1)
    {
        F77_FUNC(sgetv0, SGETV0)
        (ido,
         bmat,
         &c__1,
         &iwork[3],
         n,
         &c__1,
         &v[v_offset],
         ldv,
         &resid[1],
         &workd[*n * 3 + 1],
         &ipntr[1],
         &workd[1],
         &iwork[41],
         info);

        if (*ido != 99)
        {
            goto L9000;
        }

        if (workd[*n * 3 + 1] == 0.)
        {

            *info = -9;
            goto L1200;
        }
        iwork[2] = 0;
        *ido     = 0;
    }

    if (iwork[4] == 1)
    {
        goto L20;
    }

    if (iwork[5] == 1)
    {
        goto L50;
    }

    if (iwork[1] == 1)
    {
        goto L100;
    }

    F77_FUNC(ssaitr, SSAITR)
    (ido,
     bmat,
     n,
     &c__0,
     &iwork[9],
     mode,
     &resid[1],
     &workd[*n * 3 + 1],
     &v[v_offset],
     ldv,
     &h__[h_offset],
     ldh,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {

        *np     = *info;
        *mxiter = iwork[6];
        *info   = -9999;
        goto L1200;
    }

L1000:

    ++iwork[6];


    *ido = 0;
L20:
    iwork[4] = 1;

    F77_FUNC(ssaitr, SSAITR)
    (ido,
     bmat,
     n,
     nev,
     np,
     mode,
     &resid[1],
     &workd[*n * 3 + 1],
     &v[v_offset],
     ldv,
     &h__[h_offset],
     ldh,
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido != 99)
    {
        goto L9000;
    }

    if (*info > 0)
    {

        *np     = *info;
        *mxiter = iwork[6];
        *info   = -9999;
        goto L1200;
    }
    iwork[4] = 0;

    F77_FUNC(sseigt, SSEIGT)
    (&workd[*n * 3 + 1], &iwork[7], &h__[h_offset], ldh, &ritz[1], &bounds[1], &workl[1], &ierr);

    if (ierr != 0)
    {
        *info = -8;
        goto L1200;
    }

    F77_FUNC(scopy, SCOPY)(&iwork[7], &ritz[1], &c__1, &workl[iwork[7] + 1], &c__1);
    F77_FUNC(scopy, SCOPY)(&iwork[7], &bounds[1], &c__1, &workl[(iwork[7] << 1) + 1], &c__1);

    *nev = iwork[9];
    *np  = iwork[10];
    F77_FUNC(ssgets, SSGETS)(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);

    F77_FUNC(scopy, SCOPY)(nev, &bounds[*np + 1], &c__1, &workl[*np + 1], &c__1);
    F77_FUNC(ssconv, SSCONV)(nev, &ritz[*np + 1], &workl[*np + 1], tol, &iwork[8]);


    nptemp = *np;
    i__1   = nptemp;
    for (j = 1; j <= i__1; ++j)
    {
        if (bounds[j] == 0.)
        {
            --(*np);
            ++(*nev);
        }
    }

    if (iwork[8] >= iwork[9] || iwork[6] > *mxiter || *np == 0)
    {

        if (!std::strncmp(which, "BE", 2))
        {

            std::strncpy(wprime, "SA", 2);
            F77_FUNC(ssortr, SSORTR)(wprime, &c__1, &iwork[7], &ritz[1], &bounds[1]);
            nevd2 = *nev / 2;
            nevm2 = *nev - nevd2;
            if (*nev > 1)
            {
                i__1 = (nevd2 < *np) ? nevd2 : *np;
                i__2 = iwork[7] - nevd2 + 1, i__3 = iwork[7] - *np + 1;
                F77_FUNC(sswap, SSWAP)
                (&i__1, &ritz[nevm2 + 1], &c__1, &ritz[((i__2 > i__3) ? i__2 : i__3)], &c__1);
                i__1 = (nevd2 < *np) ? nevd2 : *np;
                i__2 = iwork[7] - nevd2 + 1, i__3 = iwork[7] - *np;
                F77_FUNC(sswap, SSWAP)
                (&i__1, &bounds[nevm2 + 1], &c__1, &bounds[((i__2 > i__3) ? i__2 : i__3) + 1], &c__1);
            }
        }
        else
        {

            if (!std::strncmp(which, "LM", 2))
            {
                std::strncpy(wprime, "SM", 2);
            }
            if (!std::strncmp(which, "SM", 2))
            {
                std::strncpy(wprime, "LM", 2);
            }
            if (!std::strncmp(which, "LA", 2))
            {
                std::strncpy(wprime, "SA", 2);
            }
            if (!std::strncmp(which, "SA", 2))
            {
                std::strncpy(wprime, "LA", 2);
            }

            F77_FUNC(ssortr, SSORTR)(wprime, &c__1, &iwork[7], &ritz[1], &bounds[1]);
        }

        i__1 = iwork[9];
        for (j = 1; j <= i__1; ++j)
        {
            d__2 = eps23;
            d__3 = std::abs(ritz[j]);
            temp = (d__2 > d__3) ? d__2 : d__3;
            bounds[j] /= temp;
        }

        std::strncpy(wprime, "LA", 2);
        F77_FUNC(ssortr, SSORTR)(wprime, &c__1, &iwork[9], &bounds[1], &ritz[1]);

        i__1 = iwork[9];
        for (j = 1; j <= i__1; ++j)
        {
            d__2 = eps23;
            d__3 = std::abs(ritz[j]);
            temp = (d__2 > d__3) ? d__2 : d__3;
            bounds[j] *= temp;
        }

        if (!std::strncmp(which, "BE", 2))
        {

            std::strncpy(wprime, "LA", 2);
            F77_FUNC(ssortr, SSORTR)(wprime, &c__1, &iwork[8], &ritz[1], &bounds[1]);
        }
        else
        {
            F77_FUNC(ssortr, SSORTR)(which, &c__1, &iwork[8], &ritz[1], &bounds[1]);
        }

        h__[h_dim1 + 1] = workd[*n * 3 + 1];


        if (iwork[6] > *mxiter && iwork[8] < *nev)
        {
            *info = 1;
        }
        if (*np == 0 && iwork[8] < iwork[9])
        {
            *info = 2;
        }

        *np = iwork[8];
        goto L1100;
    }
    else if (iwork[8] < *nev && *ishift == 1)
    {
        nevbef = *nev;
        i__1 = iwork[8], i__2 = *np / 2;
        *nev += (i__1 < i__2) ? i__1 : i__2;
        if (*nev == 1 && iwork[7] >= 6)
        {
            *nev = iwork[7] / 2;
        }
        else if (*nev == 1 && iwork[7] > 2)
        {
            *nev = 2;
        }
        *np = iwork[7] - *nev;


        if (nevbef < *nev)
        {
            F77_FUNC(ssgets, SSGETS)(ishift, which, nev, np, &ritz[1], &bounds[1], &workl[1]);
        }
    }


    if (*ishift == 0)
    {

        iwork[5] = 1;
        *ido     = 3;
        goto L9000;
    }

L50:

    iwork[5] = 0;

    if (*ishift == 0)
    {
        F77_FUNC(scopy, SCOPY)(np, &workl[1], &c__1, &ritz[1], &c__1);
    }

    F77_FUNC(ssapps, SSAPPS)
    (n, nev, np, &ritz[1], &v[v_offset], ldv, &h__[h_offset], ldh, &resid[1], &q[q_offset], ldq, &workd[1]);

    iwork[1] = 1;
    if (*bmat == 'G')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[*n + 1], &c__1);
        ipntr[1] = *n + 1;
        ipntr[2] = 1;
        *ido     = 2;

        goto L9000;
    }
    else if (*bmat == 'I')
    {
        F77_FUNC(scopy, SCOPY)(n, &resid[1], &c__1, &workd[1], &c__1);
    }

L100:

    if (*bmat == 'G')
    {
        workd[*n * 3 + 1] = F77_FUNC(sdot, SDOT)(n, &resid[1], &c__1, &workd[1], &c__1);
        workd[*n * 3 + 1] = std::sqrt(std::abs(workd[*n * 3 + 1]));
    }
    else if (*bmat == 'I')
    {
        workd[*n * 3 + 1] = F77_FUNC(snrm2, SNRM2)(n, &resid[1], &c__1);
    }
    iwork[1] = 0;

    goto L1000;

L1100:

    *mxiter = iwork[6];
    *nev    = iwork[8];

L1200:
    *ido = 99;

L9000:
    return;
}


void F77_FUNC(ssaupd, SSAUPD)(int*        ido,
                              const char* bmat,
                              int*        n,
                              const char* which,
                              int*        nev,
                              float*      tol,
                              float*      resid,
                              int*        ncv,
                              float*      v,
                              int*        ldv,
                              int*        iparam,
                              int*        ipntr,
                              float*      workd,
                              int*        iwork,
                              float*      workl,
                              int*        lworkl,
                              int*        info)
{
    int v_dim1, v_offset, i__1, i__2;
    int j;

    --workd;
    --resid;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --iwork;
    --workl;

    if (*ido == 0)
    {


        iwork[2]  = 0;
        iwork[5]  = iparam[1];
        iwork[10] = iparam[3];
        iwork[12] = iparam[4];

        iwork[6]  = 1;
        iwork[11] = iparam[7];


        if (*n <= 0)
        {
            iwork[2] = -1;
        }
        else if (*nev <= 0)
        {
            iwork[2] = -2;
        }
        else if (*ncv <= *nev || *ncv > *n)
        {
            iwork[2] = -3;
        }


        iwork[15] = *ncv - *nev;

        if (iwork[10] <= 0)
        {
            iwork[2] = -4;
        }
        if (std::strncmp(which, "LM", 2) && std::strncmp(which, "SM", 2) && std::strncmp(which, "LA", 2)
            && std::strncmp(which, "SA", 2) && std::strncmp(which, "BE", 2))
        {
            iwork[2] = -5;
        }
        if (*bmat != 'I' && *bmat != 'G')
        {
            iwork[2] = -6;
        }

        i__1 = *ncv;
        if (*lworkl < i__1 * i__1 + (*ncv << 3))
        {
            iwork[2] = -7;
        }
        if (iwork[11] < 1 || iwork[11] > 5)
        {
            iwork[2] = -10;
        }
        else if (iwork[11] == 1 && *bmat == 'G')
        {
            iwork[2] = -11;
        }
        else if (iwork[5] < 0 || iwork[5] > 1)
        {
            iwork[2] = -12;
        }
        else if (*nev == 1 && !std::strncmp(which, "BE", 2))
        {
            iwork[2] = -13;
        }

        if (iwork[2] != 0)
        {
            *info = iwork[2];
            *ido  = 99;
            goto L9000;
        }

        if (iwork[12] <= 0)
        {
            iwork[12] = 1;
        }
        if (*tol <= 0.)
        {
            *tol = GMX_FLOAT_EPS;
        }

        iwork[15] = *ncv - *nev;
        iwork[13] = *nev;
        i__2      = *ncv;
        i__1      = i__2 * i__2 + (*ncv << 3);
        for (j = 1; j <= i__1; ++j)
        {
            workl[j] = 0.;
        }

        iwork[8]  = *ncv;
        iwork[9]  = *ncv;
        iwork[3]  = 1;
        iwork[16] = iwork[3] + (iwork[8] << 1);
        iwork[1]  = iwork[16] + *ncv;
        iwork[4]  = iwork[1] + *ncv;
        i__1      = *ncv;
        iwork[7]  = iwork[4] + i__1 * i__1;
        iwork[14] = iwork[7] + *ncv * 3;

        ipntr[4]  = iwork[14];
        ipntr[5]  = iwork[3];
        ipntr[6]  = iwork[16];
        ipntr[7]  = iwork[1];
        ipntr[11] = iwork[7];
    }

    F77_FUNC(ssaup2, SSAUP2)
    (ido,
     bmat,
     n,
     which,
     &iwork[13],
     &iwork[15],
     tol,
     &resid[1],
     &iwork[11],
     &iwork[6],
     &iwork[5],
     &iwork[10],
     &v[v_offset],
     ldv,
     &workl[iwork[3]],
     &iwork[8],
     &workl[iwork[16]],
     &workl[iwork[1]],
     &workl[iwork[4]],
     &iwork[9],
     &workl[iwork[7]],
     &ipntr[1],
     &workd[1],
     &iwork[21],
     info);

    if (*ido == 3)
    {
        iparam[8] = iwork[15];
    }
    if (*ido != 99)
    {
        goto L9000;
    }

    iparam[3] = iwork[10];
    iparam[5] = iwork[15];

    if (*info < 0)
    {
        goto L9000;
    }
    if (*info == 2)
    {
        *info = 3;
    }

L9000:

    return;
}


void F77_FUNC(sseupd, SSEUPD)(int*        rvec,
                              const char* howmny,
                              int*        select,
                              float*      d__,
                              float*      z__,
                              int*        ldz,
                              float*      sigma,
                              const char* bmat,
                              int*        n,
                              const char* which,
                              int*        nev,
                              float*      tol,
                              float*      resid,
                              int*        ncv,
                              float*      v,
                              int*        ldv,
                              int*        iparam,
                              int*        ipntr,
                              float*      workd,
                              float*      workl,
                              int*        lworkl,
                              int*        info)
{
    float c_b21  = 2 / 3.;
    int   c__1   = 1;
    float c_b102 = 1.;
    int   v_dim1, v_offset, z_dim1, z_offset, i__1;
    float d__1, d__2, d__3;

    int   j, k, ih, iq, iw, ibd, ihb, ihd, ldh, ilg, ldq, ism, irz;
    int   mode;
    float eps23;
    int   ierr;
    float temp;
    int   next;
    char  type__[6];
    int   ritz;
    int   reord;
    int   nconv;
    float rnorm;
    float bnorm2;
    float thres1 = 0, thres2 = 0;
    int   bounds;
    float tempbnd;
    int   leftptr, rghtptr;


    --workd;
    --resid;
    z_dim1   = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --d__;
    --select;
    v_dim1   = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --iparam;
    --ipntr;
    --workl;

    mode  = iparam[7];
    nconv = iparam[5];
    *info = 0;

    if (nconv == 0)
    {
        goto L9000;
    }
    ierr = 0;

    if (nconv <= 0)
    {
        ierr = -14;
    }
    if (*n <= 0)
    {
        ierr = -1;
    }
    if (*nev <= 0)
    {
        ierr = -2;
    }
    if (*ncv <= *nev || *ncv > *n)
    {
        ierr = -3;
    }
    if (std::strncmp(which, "LM", 2) && std::strncmp(which, "SM", 2) && std::strncmp(which, "LA", 2)
        && std::strncmp(which, "SA", 2) && std::strncmp(which, "BE", 2))
    {
        ierr = -5;
    }
    if (*bmat != 'I' && *bmat != 'G')
    {
        ierr = -6;
    }
    if (*howmny != 'A' && *howmny != 'P' && *howmny != 'S' && *rvec)
    {
        ierr = -15;
    }
    if (*rvec && *howmny == 'S')
    {
        ierr = -16;
    }
    i__1 = *ncv;
    if (*rvec && *lworkl < i__1 * i__1 + (*ncv << 3))
    {
        ierr = -7;
    }

    if (mode == 1 || mode == 2)
    {
        std::strncpy(type__, "REGULR", 6);
    }
    else if (mode == 3)
    {
        std::strncpy(type__, "SHIFTI", 6);
    }
    else if (mode == 4)
    {
        std::strncpy(type__, "BUCKLE", 6);
    }
    else if (mode == 5)
    {
        std::strncpy(type__, "CAYLEY", 6);
    }
    else
    {
        ierr = -10;
    }
    if (mode == 1 && *bmat == 'G')
    {
        ierr = -11;
    }
    if (*nev == 1 && !std::strncmp(which, "BE", 2))
    {
        ierr = -12;
    }

    if (ierr != 0)
    {
        *info = ierr;
        goto L9000;
    }

    ih        = ipntr[5];
    ritz      = ipntr[6];
    bounds    = ipntr[7];
    ldh       = *ncv;
    ldq       = *ncv;
    ihd       = bounds + ldh;
    ihb       = ihd + ldh;
    iq        = ihb + ldh;
    iw        = iq + ldh * *ncv;
    next      = iw + (*ncv << 1);
    ipntr[4]  = next;
    ipntr[8]  = ihd;
    ipntr[9]  = ihb;
    ipntr[10] = iq;

    irz = ipntr[11] + *ncv;
    ibd = irz + *ncv;


    eps23 = GMX_FLOAT_EPS;
    eps23 = std::pow(eps23, c_b21);

    rnorm = workl[ih];
    if (*bmat == 'I')
    {
        bnorm2 = rnorm;
    }
    else if (*bmat == 'G')
    {
        bnorm2 = F77_FUNC(snrm2, SNRM2)(n, &workd[1], &c__1);
    }

    if (*rvec)
    {

        if (!std::strncmp(which, "LM", 2) || !std::strncmp(which, "SM", 2)
            || !std::strncmp(which, "LA", 2) || !std::strncmp(which, "SA", 2))
        {
        }
        else if (!std::strncmp(which, "BE", 2))
        {


            ism = (*nev > nconv) ? *nev : nconv;
            ism /= 2;
            ilg    = ism + 1;
            thres1 = workl[ism];
            thres2 = workl[ilg];
        }

        reord = 0;
        i__1  = *ncv - 1;
        for (j = 0; j <= i__1; ++j)
        {
            select[j + 1] = 0;
            if (!std::strncmp(which, "LM", 2))
            {
                if (std::abs(workl[irz + j]) >= std::abs(thres1))
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "SM", 2))
            {
                if (std::abs(workl[irz + j]) <= std::abs(thres1))
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "LA", 2))
            {
                if (workl[irz + j] >= thres1)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "SA", 2))
            {
                if (workl[irz + j] <= thres1)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            else if (!std::strncmp(which, "BE", 2))
            {
                if (workl[irz + j] <= thres1 || workl[irz + j] >= thres2)
                {
                    d__2    = eps23;
                    d__3    = std::abs(workl[irz + j]);
                    tempbnd = (d__2 > d__3) ? d__2 : d__3;
                    if (workl[ibd + j] <= *tol * tempbnd)
                    {
                        select[j + 1] = 1;
                    }
                }
            }
            if (j + 1 > nconv)
            {
                reord = select[j + 1] || reord;
            }
        }

        i__1 = *ncv - 1;
        F77_FUNC(scopy, SCOPY)(&i__1, &workl[ih + 1], &c__1, &workl[ihb], &c__1);
        F77_FUNC(scopy, SCOPY)(ncv, &workl[ih + ldh], &c__1, &workl[ihd], &c__1);

        F77_FUNC(ssteqr, SSTEQR)
        ("Identity", ncv, &workl[ihd], &workl[ihb], &workl[iq], &ldq, &workl[iw], &ierr);

        if (ierr != 0)
        {
            *info = -8;
            goto L9000;
        }


        if (reord)
        {

            leftptr = 1;
            rghtptr = *ncv;

            if (*ncv == 1)
            {
                goto L30;
            }

        L20:
            if (select[leftptr])
            {

                ++leftptr;
            }
            else if (!select[rghtptr])
            {

                --rghtptr;
            }
            else
            {

                temp                     = workl[ihd + leftptr - 1];
                workl[ihd + leftptr - 1] = workl[ihd + rghtptr - 1];
                workl[ihd + rghtptr - 1] = temp;
                F77_FUNC(scopy, SCOPY)
                (ncv, &workl[iq + *ncv * (leftptr - 1)], &c__1, &workl[iw], &c__1);
                F77_FUNC(scopy, SCOPY)
                (ncv, &workl[iq + *ncv * (rghtptr - 1)], &c__1, &workl[iq + *ncv * (leftptr - 1)], &c__1);
                F77_FUNC(scopy, SCOPY)
                (ncv, &workl[iw], &c__1, &workl[iq + *ncv * (rghtptr - 1)], &c__1);
                ++leftptr;
                --rghtptr;
            }

            if (leftptr < rghtptr)
            {
                goto L20;
            }

        L30:;
        }

        F77_FUNC(scopy, SCOPY)(&nconv, &workl[ihd], &c__1, &d__[1], &c__1);
    }
    else
    {

        F77_FUNC(scopy, SCOPY)(&nconv, &workl[ritz], &c__1, &d__[1], &c__1);
        F77_FUNC(scopy, SCOPY)(ncv, &workl[ritz], &c__1, &workl[ihd], &c__1);
    }
    if (!std::strncmp(type__, "REGULR", 6))
    {

        if (*rvec)
        {
            F77_FUNC(ssesrt, SSESRT)("LA", rvec, &nconv, &d__[1], ncv, &workl[iq], &ldq);
        }
        else
        {
            F77_FUNC(scopy, SCOPY)(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
        }
    }
    else
    {

        F77_FUNC(scopy, SCOPY)(ncv, &workl[ihd], &c__1, &workl[iw], &c__1);
        if (!std::strncmp(type__, "SHIFTI", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = 1. / workl[ihd + k - 1] + *sigma;
            }
        }
        else if (!std::strncmp(type__, "BUCKLE", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = *sigma * workl[ihd + k - 1] / (workl[ihd + k - 1] - 1.);
            }
        }
        else if (!std::strncmp(type__, "CAYLEY", 6))
        {
            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihd + k - 1] = *sigma * (workl[ihd + k - 1] + 1.) / (workl[ihd + k - 1] - 1.);
            }
        }

        F77_FUNC(scopy, SCOPY)(&nconv, &workl[ihd], &c__1, &d__[1], &c__1);
        F77_FUNC(ssortr, SSORTR)("LA", &c__1, &nconv, &workl[ihd], &workl[iw]);
        if (*rvec)
        {
            F77_FUNC(ssesrt, SSESRT)("LA", rvec, &nconv, &d__[1], ncv, &workl[iq], &ldq);
        }
        else
        {
            F77_FUNC(scopy, SCOPY)(ncv, &workl[bounds], &c__1, &workl[ihb], &c__1);
            d__1 = bnorm2 / rnorm;
            F77_FUNC(sscal, SSCAL)(ncv, &d__1, &workl[ihb], &c__1);
            F77_FUNC(ssortr, SSORTR)("LA", &c__1, &nconv, &d__[1], &workl[ihb]);
        }
    }

    if (*rvec && *howmny == 'A')
    {

        F77_FUNC(sgeqr2, SGEQR2)
        (ncv, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb], &ierr);

        F77_FUNC(sorm2r, SORM2R)
        ("Right",
         "Notranspose",
         n,
         ncv,
         &nconv,
         &workl[iq],
         &ldq,
         &workl[iw + *ncv],
         &v[v_offset],
         ldv,
         &workd[*n + 1],
         &ierr);
        F77_FUNC(slacpy, SLACPY)("All", n, &nconv, &v[v_offset], ldv, &z__[z_offset], ldz);

        i__1 = *ncv - 1;
        for (j = 1; j <= i__1; ++j)
        {
            workl[ihb + j - 1] = 0.;
        }
        workl[ihb + *ncv - 1] = 1.;
        F77_FUNC(sorm2r, SORM2R)
        ("Left", "Transpose", ncv, &c__1, &nconv, &workl[iq], &ldq, &workl[iw + *ncv], &workl[ihb], ncv, &temp, &ierr);
    }
    else if (*rvec && *howmny == 'S')
    {
    }

    if (!std::strncmp(type__, "REGULR", 6) && *rvec)
    {

        i__1 = *ncv;
        for (j = 1; j <= i__1; ++j)
        {
            workl[ihb + j - 1] = rnorm * std::abs(workl[ihb + j - 1]);
        }
    }
    else if (std::strncmp(type__, "REGULR", 6) && *rvec)
    {

        F77_FUNC(sscal, SSCAL)(ncv, &bnorm2, &workl[ihb], &c__1);
        if (!std::strncmp(type__, "SHIFTI", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                d__2               = workl[iw + k - 1];
                workl[ihb + k - 1] = std::abs(workl[ihb + k - 1]) / (d__2 * d__2);
            }
        }
        else if (!std::strncmp(type__, "BUCKLE", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                d__2               = workl[iw + k - 1] - 1.;
                workl[ihb + k - 1] = *sigma * std::abs(workl[ihb + k - 1]) / (d__2 * d__2);
            }
        }
        else if (!std::strncmp(type__, "CAYLEY", 6))
        {

            i__1 = *ncv;
            for (k = 1; k <= i__1; ++k)
            {
                workl[ihb + k - 1] =
                        std::abs(workl[ihb + k - 1] / workl[iw + k - 1] * (workl[iw + k - 1] - 1.));
            }
        }
    }

    if (*rvec && (!std::strncmp(type__, "SHIFTI", 6) || !std::strncmp(type__, "CAYLEY", 6)))
    {

        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iq + k * ldq + *ncv - 1] / workl[iw + k];
        }
    }
    else if (*rvec && !std::strncmp(type__, "BUCKLE", 6))
    {

        i__1 = nconv - 1;
        for (k = 0; k <= i__1; ++k)
        {
            workl[iw + k] = workl[iq + k * ldq + *ncv - 1] / (workl[iw + k] - 1.);
        }
    }

    if (std::strncmp(type__, "REGULR", 6))
    {
        F77_FUNC(sger, SGER)
        (n, &nconv, &c_b102, &resid[1], &c__1, &workl[iw], &c__1, &z__[z_offset], ldz);
    }

L9000:

    return;
}

#endif
