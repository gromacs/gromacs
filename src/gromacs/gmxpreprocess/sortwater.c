/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#include "sortwater.h"

#include <stdlib.h>

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/math/vec.h"
#include "gromacs/random/random.h"
#include "gromacs/utility/smalloc.h"

static rvec   *xptr, box_1;
static int     nwat;
static matrix  BOX;
static ivec    NBOX;

void randwater(int astart, int nwater, int nwatom, rvec x[], rvec v[],
               gmx_rng_t rng)
{
    int  i, j, wi, wj, *tab;
    rvec buf;

    snew(tab, nwater);
    for (i = 0; (i < nwater); i++)
    {
        tab[i] = i;
    }
    for (j = 0; (j < 23*nwater); j++)
    {
        wi = (int) (nwater*gmx_rng_uniform_real(rng)) % nwater;
        do
        {
            wj = (int) (nwater*gmx_rng_uniform_real(rng)) % nwater;
        }
        while (wi == wj);
        wi = astart+wi*nwatom;
        wj = astart+wj*nwatom;
        /* Swap coords for wi and wj */
        for (i = 0; (i < nwatom); i++)
        {
            copy_rvec(x[wi+i], buf);
            copy_rvec(x[wj+i], x[wi+i]);
            copy_rvec(buf, x[wj+i]);
            if (v)
            {
                copy_rvec(v[wi+i], buf);
                copy_rvec(v[wj+i], v[wi+i]);
                copy_rvec(buf, v[wj+i]);
            }
        }
    }
    sfree(tab);
}

static int rvcomp(const void *a, const void *b)
{
    int aa, bb;

    aa = nwat*(*((int *)a));
    bb = nwat*(*((int *)b));
    if (xptr[aa][XX] < xptr[bb][XX])
    {
        return -1;
    }
    else if (xptr[aa][XX] > xptr[bb][XX])
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

static int block_index(rvec x, ivec nbox)
{
    ivec ixyz;
    int  m;

    for (m = 0; (m < DIM); m++)
    {
        ixyz[m] = ((int)((1+x[m]*box_1[m])*nbox[m])) % nbox[m];
    }
    return ixyz[XX]*nbox[YY]*nbox[ZZ]+ixyz[YY]*nbox[ZZ]+ixyz[ZZ];
}

static int blockcomp(const void *a, const void *b)
{
    int aa, bb, aind, bind;

    aa   = nwat*(*((int *)a));
    bb   = nwat*(*((int *)b));

    aind = block_index(xptr[aa], NBOX);
    bind = block_index(xptr[bb], NBOX);

    if (aind == bind)
    {
        if (xptr[aa][XX] < xptr[bb][XX])
        {
            return -1;
        }
        else if (xptr[aa][XX] > xptr[bb][XX])
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        return aind-bind;
    }
}

static void lo_sortwater(int astart, int nwater, int nwatom, rvec x[], rvec v[],
                         gmx_bool bBlock)
{
    int   i, j, i0, rvi;
    int  *rvindex;
    rvec *tmp;

    /* Sort indices to rvecs */
    snew(rvindex, nwater);
    for (i = 0; (i < nwater); i++)
    {
        rvindex[i] = i;
    }
    xptr = x+astart;
    nwat = nwatom;

    qsort(rvindex, nwater, sizeof(rvindex[0]), bBlock ? blockcomp : rvcomp);
    if (debug)
    {
        for (i = 0; (i < nwater); i++)
        {
            rvi = rvindex[i]*nwatom;
            fprintf(debug, "rvindex[%5d] = %5d (x = %8.3f  %8.3f  %8.3f)\n",
                    i, rvi, x[astart+rvi][XX], x[astart+rvi][YY], x[astart+rvi][ZZ]);
        }
    }
    snew(tmp, nwater*nwatom);

    for (i = 0; (i < nwater); i++)
    {
        i0 = astart+nwatom*rvindex[i];
        for (j = 0; (j < nwatom); j++)
        {
            copy_rvec(x[i0+j], tmp[nwatom*i+j]);
        }
    }
    for (i = 0; (i < nwater*nwatom); i++)
    {
        copy_rvec(tmp[i], x[astart+i]);
    }

    for (i = 0; (i < nwater); i++)
    {
        i0 = astart+nwatom*rvindex[i];
        for (j = 0; (j < nwatom); j++)
        {
            copy_rvec(v[i0+j], tmp[nwatom*i+j]);
        }
    }
    for (i = 0; (i < nwater*nwatom); i++)
    {
        copy_rvec(tmp[i], v[astart+i]);
    }

    sfree(tmp);
    sfree(rvindex);
}

void sortwater(int astart, int nwater, int nwatom, rvec x[], rvec v[])
{
    lo_sortwater(astart, nwater, nwatom, x, v, FALSE);
}

static void factorize(int nn, int fac[])
{
    int i, n = nn;

    for (i = 0; (i <= n); i++)
    {
        fac[i] = 0;
    }
    fac[1] = 1;
    for (i = 2; (i <= n); )
    {
        if ((n % i) == 0)
        {
            fac[i]++;
            n = n/i;
        }
        else
        {
            i++;
        }
    }
    if (debug)
    {
        fprintf(debug, "Factorizing %d into primes:\n", nn);
        for (i = 2; (i <= nn); i++)
        {
            if (fac[i])
            {
                fprintf(debug, "%d ^ %d\n", i, fac[i]);
            }
        }
    }
}

static int ipow(int base, int exp)
{
    int i, ip;

    for (ip = 1, i = 0; (i < exp); i++)
    {
        ip *= base;
    }
    return ip;
}

static int iv_comp(const void *a, const void *b)
{
    int *ia, *ib;

    ia = (int *)a;
    ib = (int *)b;
    if (ia[XX] != ib[XX])
    {
        return (ia[XX] - ib[XX]);
    }
    else if (ia[YY] != ib[YY])
    {
        return (ia[YY] - ib[YY]);
    }
    else
    {
        return (ia[ZZ] - ib[ZZ]);
    }
}

static int add_bb(ivec BB[], int n, ivec b)
{
#define SWPX(vv, xx, yy) { int tmp; tmp = vv[xx]; vv[xx] = vv[yy]; vv[yy] = tmp; }
    copy_ivec(b, BB[n++]); /* x y z */
    SWPX(b, XX, YY);
    copy_ivec(b, BB[n++]); /* y x z */
    SWPX(b, XX, ZZ);
    copy_ivec(b, BB[n++]); /* z x y */
    SWPX(b, XX, YY);
    copy_ivec(b, BB[n++]); /* x z y */
    SWPX(b, XX, ZZ);
    copy_ivec(b, BB[n++]); /* y z x */
    SWPX(b, XX, YY);
    copy_ivec(b, BB[n++]); /* z y x */
    SWPX(b, XX, ZZ);       /* Back to normal */
#undef SWPX
    return n;
}

static real box_weight(ivec nbox, matrix box)
{
    rvec lx;
    int  m;

    /* Calculate area of subbox */
    for (m = 0; (m < DIM); m++)
    {
        lx[m] = box[m][m]/nbox[m];
    }
    return 2*(lx[XX]*lx[YY]+lx[XX]*lx[ZZ]+lx[YY]*lx[ZZ]);
}

static int w_comp(const void *a, const void *b)
{
    int *ia, *ib;
    real wa, wb;

    ia = (int *)a;
    ib = (int *)b;

    wa = box_weight(ia, BOX);
    wb = box_weight(ib, BOX);
    if (fabs(wa - wb) < 1e-4)
    {
        return (iiprod(ia, ia) - iiprod(ib, ib));
    }
    else if (wa < wb)
    {
        return -1;
    }
    else
    {
        return 1;
    }
}

static void buildbox(int nnode, ivec nbox, matrix box)
{
    ivec *BB, bxyz;
    int   i, j, m, n, n3, ny, *fx, *fy, nbb;

    n3 = ipow(nnode, 3)*6;
    snew(BB, n3);
    nbb = 0;
    snew(fx, nnode+1);
    snew(fy, nnode+1);
    factorize(nnode, fx);
    for (i = 0; (i <= nnode); i++)
    {
        for (m = 1; (m <= fx[i]); m++)
        {
            bxyz[XX] = ipow(i, m);
            ny       = nnode/bxyz[XX];
            factorize(ny, fy);
            for (j = 0; (j <= ny); j++)
            {
                for (n = 1; (n <= fy[j]); n++)
                {
                    bxyz[YY] = ipow(j, n);
                    bxyz[ZZ] = ny/bxyz[YY];
                    if (bxyz[ZZ] > 0)
                    {
                        nbb = add_bb(BB, nbb, bxyz);
                    }
                }
            }
        }
    }
    /* Sort boxes and remove doubles */
    qsort(BB, nbb, sizeof(BB[0]), iv_comp);
    j = 0;
    for (i = 1; (i < nbb); i++)
    {
        if ((BB[i][XX] != BB[j][XX]) ||
            (BB[i][YY] != BB[j][YY]) ||
            (BB[i][ZZ] != BB[j][ZZ]))
        {
            j++;
            copy_ivec(BB[i], BB[j]);
        }
    }
    nbb = ++j;
    /* Sort boxes according to weight */
    copy_mat(box, BOX);
    qsort(BB, nbb, sizeof(BB[0]), w_comp);
    for (i = 0; (i < nbb); i++)
    {
        fprintf(stderr, "nbox = %2d %2d %2d [ prod %3d ] area = %12.5f (nm^2)\n",
                BB[i][XX], BB[i][YY], BB[i][ZZ],
                BB[i][XX]*BB[i][YY]*BB[i][ZZ],
                box_weight(BB[i], box));
    }
    copy_ivec(BB[0], nbox);
    sfree(BB);
    sfree(fy);
    sfree(fx);
}

void mkcompact(int astart, int nwater, int nwatom, rvec x[], rvec v[],
               int nnode, matrix box)
{
    /* Make a compact configuration for each processor.
     * Divide the computational box in near cubic boxes and spread them
     * evenly over processors.
     */
/*   ivec nbox; */
    int  m;

    if (nnode <= 1)
    {
        return;
    }

    buildbox(nnode, NBOX, box);
    /* copy_ivec(nbox,NBOX); */
    for (m = 0; (m < DIM); m++)
    {
        box_1[m] = 1.0/box[m][m];
    }

    lo_sortwater(astart, nwater, nwatom, x, v, TRUE);
}
