/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2012,2014, by the GROMACS development team, led by
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

#include "qsort_threadsafe.h"

#include <stdlib.h>

static void
qsort_swapfunc(void *        a,
               void *        b,
               size_t        n,
               int           swaptype)
{
    int  *  ia;
    int  *  ib;
    int     itmp;

    char *  ca;
    char *  cb;
    char    ctmp;

    if (swaptype <= 1)
    {
        ia = (int *)a;
        ib = (int *)b;
        for (; n > 0; ia += 1, ib += 1, n -= sizeof(int))
        {
            itmp      = *ia;
            *ia       = *ib;
            *ib       = itmp;
        }
    }
    else
    {
        ca = (char *)a;
        cb = (char *)b;
        for (; n > 0; ca += 1, cb += 1, n -= 1)
        {
            ctmp       = *ca;
            *ca        = *cb;
            *cb        = ctmp;
        }
    }
}


static void *
qsort_med3(void *        a,
           void *        b,
           void *        c,
           int          (*compar) (const void *a, const void *b))
{
    if (compar(a, b) < 0)
    {
        if (compar(b, c) < 0)
        {
            return b;
        }
        else if (compar(a, c) < 0)
        {
            return c;
        }
        else
        {
            return a;
        }
    }
    else
    {
        if (compar(b, c) > 0)
        {
            return b;
        }
        else if (compar(a, c) > 0)
        {
            return c;
        }
        else
        {
            return a;
        }
    }
}


void
gmx_qsort(void *           base,
          size_t           nmemb,
          size_t           size,
          int            (*compar)(const void *, const void *))
{
#define QSORT_EXCH(a, b, t) (t = a, a = b, b = t)
#define QSORT_SWAP(a, b) swaptype != 0 ? qsort_swapfunc(a, b, size, swaptype) : \
    (void)QSORT_EXCH(*(int *)(a), *(int *)(b), t)

    char  *pa, *pb, *pc, *pd, *pl, *pm, *pn, *pv, *cbase;
    int    r, swaptype;
    int    t, v;
    size_t s, st;

    cbase = (char *)base;

    swaptype = (size_t)(cbase - (char *)0) % sizeof(int) || size % sizeof(int) ? 2 : size == sizeof(int) ? 0 : 1;

    if (nmemb < 7)
    {
        /* Insertion sort on smallest arrays */
        for (pm = cbase + size; pm < cbase + nmemb*size; pm += size)
        {
            for (pl = pm; (pl > cbase) && compar((void *)(pl-size), (void *) pl) > 0; pl -= size)
            {
                QSORT_SWAP(pl, pl-size);
            }
        }
        return;
    }

    /* Small arrays, middle element */
    pm = cbase + (nmemb/2)*size;

    if (nmemb > 7)
    {
        pl = cbase;
        pn = cbase + (nmemb-1)*size;
        if (nmemb > 40)
        {
            /* Big arrays, pseudomedian of 9 */
            s  = (nmemb/8)*size;
            pl = (char *)qsort_med3((void *)pl, (void *)((size_t)pl+s), (void *)((size_t)pl+2*s), compar);
            pm = (char *)qsort_med3((void *)((size_t)pm-s), (void *)pm, (void *)((size_t)pm+s), compar);
            pn = (char *)qsort_med3((void *)((size_t)pn-2*s), (void *)((size_t)pn-s), (void *)pn, compar);
        }
        /* Mid-size, med of 3 */
        pm = (char *)qsort_med3((void *)pl, (void *)pm, (void *)pn, compar);
    }

    /* pv points to partition value */
    if (swaptype != 0)
    {
        pv = cbase;
        QSORT_SWAP(pv, pm);
    }
    else
    {
        pv = (char*)(void*)&v;
        v  = *(int *)pm;
    }

    pa = pb = cbase;
    pc = pd = cbase + (nmemb-1)*size;

    for (;; )
    {
        while (pb <= pc && (r = compar((void *)pb, (void *) pv)) <= 0)
        {
            if (r == 0)
            {
                QSORT_SWAP(pa, pb);
                pa += size;
            }
            pb += size;
        }
        while (pc >= pb && (r = compar((void *)pc, (void *) pv)) >= 0)
        {
            if (r == 0)
            {
                QSORT_SWAP(pc, pd);
                pd -= size;
            }
            pc -= size;
        }
        if (pb > pc)
        {
            break;
        }
        QSORT_SWAP(pb, pc);
        pb += size;
        pc -= size;
    }
    pn = cbase + nmemb*size;

    s  = pa-cbase;
    st = pb-pa;
    if (st < s)
    {
        s = st;
    }

    if (s > 0)
    {
        qsort_swapfunc(cbase, pb-s, s, swaptype);
    }

    s  = pd-pc;
    st = pn-pd-size;
    if (st < s)
    {
        s = st;
    }

    if (s > 0)
    {
        qsort_swapfunc(pb, pn-s, s, swaptype);
    }

    if ((s = pb-pa) > size)
    {
        gmx_qsort(cbase, s/size, size, compar);
    }

    if ((s = pd-pc) > size)
    {
        gmx_qsort(pn-s, s/size, size, compar);
    }

#undef QSORT_EXCH
#undef QSORT_SWAP

    return;
}
