/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
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

#include "config.h"

#include <cmath>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/crosscorr.h"
#include "gromacs/correlationfunctions/expfit.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"

#define max_hx 7
typedef int t_hx[max_hx];
#define NRHXTYPES max_hx
const char *hxtypenames[NRHXTYPES] =
{"n-n", "n-n+1", "n-n+2", "n-n+3", "n-n+4", "n-n+5", "n-n>6"};
#define MAXHH 4

static const int NOTSET = -49297;

/* -----------------------------------------*/

enum {
    gr0,  gr1,    grI,  grNR
};
enum {
    hbNo, hbDist, hbHB, hbNR, hbR2
};
static const unsigned char c_acceptorMask = (1 << 0);
static const unsigned char c_donorMask    = (1 << 1);
static const unsigned char c_inGroupMask  = (1 << 2);


static const char *grpnames[grNR] = {"0", "1", "I" };

static gmx_bool    bDebug = FALSE;

#define HB_NO 0
#define HB_YES 1<<0
#define HB_INS 1<<1
#define HB_YESINS HB_YES|HB_INS
#define HB_NR (1<<2)
#define MAXHYDRO 4

#define ISHB(h)    ((h) & 2)
#define ISDIST(h)  ((h) & 1)
#define ISDIST2(h) ((h) & 4)
#define ISACC(h)   ((h) & c_acceptorMask)
#define ISDON(h)   ((h) & c_donorMask)
#define ISINGRP(h) ((h) & c_inGroupMask)

typedef struct {
    int      nr;
    int      maxnr;
    int     *atoms;
} t_ncell;

typedef struct {
    t_ncell d[grNR];
    t_ncell a[grNR];
} t_gridcell;

typedef int     t_icell[grNR];
typedef int h_id[MAXHYDRO];

typedef struct {
    int      history[MAXHYDRO];
    /* Has this hbond existed ever? If so as hbDist or hbHB or both.
     * Result is stored as a bitmap (1 = hbDist) || (2 = hbHB)
     */
    /* Bitmask array which tells whether a hbond is present
     * at a given time. Either of these may be NULL
     */
    int            n0;                 /* First frame a HB was found     */
    int            nframes, maxframes; /* Amount of frames in this hbond */
    unsigned int **h;
    unsigned int **g;
    /* See Xu and Berne, JPCB 105 (2001), p. 11929. We define the
     * function g(t) = [1-h(t)] H(t) where H(t) is one when the donor-
     * acceptor distance is less than the user-specified distance (typically
     * 0.35 nm).
     */
} t_hbond;

typedef struct {
    int      nra, max_nra;
    int     *acc;         /* Atom numbers of the acceptors     */
    int     *grp;         /* Group index                       */
    int     *aptr;        /* Map atom number to acceptor index */
} t_acceptors;

typedef struct {
    int       nrd, max_nrd;
    int      *don;               /* Atom numbers of the donors         */
    int      *grp;               /* Group index                        */
    int      *dptr;              /* Map atom number to donor index     */
    int      *nhydro;            /* Number of hydrogens for each donor */
    h_id     *hydro;             /* The atom numbers of the hydrogens  */
    h_id     *nhbonds;           /* The number of HBs per H at current */
} t_donors;

typedef struct {
    gmx_bool        bHBmap, bDAnr;
    int             wordlen;
    /* The following arrays are nframes long */
    int             nframes, max_frames, maxhydro;
    int            *nhb, *ndist;
    h_id           *n_bound;
    real           *time;
    t_icell        *danr;
    t_hx           *nhx;
    /* These structures are initialized from the topology at start up */
    t_donors        d;
    t_acceptors     a;
    /* This holds a matrix with all possible hydrogen bonds */
    int             nrhb, nrdist;
    t_hbond      ***hbmap;
} t_hbdata;

/* Changed argument 'bMerge' into 'oneHB' below,
 * since -contact should cause maxhydro to be 1,
 * not just -merge.
 * - Erik Marklund May 29, 2006
 */

static t_hbdata *mk_hbdata(gmx_bool bHBmap, gmx_bool bDAnr, gmx_bool oneHB)
{
    t_hbdata *hb;

    snew(hb, 1);
    hb->wordlen = 8*sizeof(unsigned int);
    hb->bHBmap  = bHBmap;
    hb->bDAnr   = bDAnr;
    if (oneHB)
    {
        hb->maxhydro = 1;
    }
    else
    {
        hb->maxhydro = MAXHYDRO;
    }
    return hb;
}

static void mk_hbmap(t_hbdata *hb)
{
    int  i, j;

    snew(hb->hbmap, hb->d.nrd);
    for (i = 0; (i < hb->d.nrd); i++)
    {
        snew(hb->hbmap[i], hb->a.nra);
        if (hb->hbmap[i] == nullptr)
        {
            gmx_fatal(FARGS, "Could not allocate enough memory for hbmap");
        }
        for (j = 0; (j > hb->a.nra); j++)
        {
            hb->hbmap[i][j] = nullptr;
        }
    }
}

static void add_frames(t_hbdata *hb, int nframes)
{
    if (nframes >= hb->max_frames)
    {
        hb->max_frames += 4096;
        srenew(hb->time, hb->max_frames);
        srenew(hb->nhb, hb->max_frames);
        srenew(hb->ndist, hb->max_frames);
        srenew(hb->n_bound, hb->max_frames);
        srenew(hb->nhx, hb->max_frames);
        if (hb->bDAnr)
        {
            srenew(hb->danr, hb->max_frames);
        }
    }
    hb->nframes = nframes;
}

#define OFFSET(frame) (frame / 32)
#define MASK(frame)   (1 << (frame % 32))

static void _set_hb(unsigned int hbexist[], unsigned int frame, gmx_bool bValue)
{
    if (bValue)
    {
        hbexist[OFFSET(frame)] |= MASK(frame);
    }
    else
    {
        hbexist[OFFSET(frame)] &= ~MASK(frame);
    }
}

static gmx_bool is_hb(unsigned int hbexist[], int frame)
{
    return ((hbexist[OFFSET(frame)] & MASK(frame)) != 0) ? 1 : 0;
}

static void set_hb(t_hbdata *hb, int id, int ih, int ia, int frame, int ihb)
{
    unsigned int *ghptr = nullptr;

    if (ihb == hbHB)
    {
        ghptr = hb->hbmap[id][ia]->h[ih];
    }
    else if (ihb == hbDist)
    {
        ghptr = hb->hbmap[id][ia]->g[ih];
    }
    else
    {
        gmx_fatal(FARGS, "Incomprehensible iValue %d in set_hb", ihb);
    }

    _set_hb(ghptr, frame-hb->hbmap[id][ia]->n0, TRUE);
}

static void add_ff(t_hbdata *hbd, int id, int h, int ia, int frame, int ihb)
{
    int         i, j, n;
    t_hbond    *hb       = hbd->hbmap[id][ia];
    int         maxhydro = std::min(hbd->maxhydro, hbd->d.nhydro[id]);
    int         wlen     = hbd->wordlen;
    int         delta    = 32*wlen;

    if (!hb->h[0])
    {
        hb->n0        = frame;
        hb->maxframes = delta;
        for (i = 0; (i < maxhydro); i++)
        {
            snew(hb->h[i], hb->maxframes/wlen);
            snew(hb->g[i], hb->maxframes/wlen);
        }
    }
    else
    {
        hb->nframes = frame-hb->n0;
        /* We need a while loop here because hbonds may be returning
         * after a long time.
         */
        while (hb->nframes >= hb->maxframes)
        {
            n = hb->maxframes + delta;
            for (i = 0; (i < maxhydro); i++)
            {
                srenew(hb->h[i], n/wlen);
                srenew(hb->g[i], n/wlen);
                for (j = hb->maxframes/wlen; (j < n/wlen); j++)
                {
                    hb->h[i][j] = 0;
                    hb->g[i][j] = 0;
                }
            }

            hb->maxframes = n;
        }
    }
    if (frame >= 0)
    {
        set_hb(hbd, id, h, ia, frame, ihb);
    }

}

static void inc_nhbonds(t_donors *ddd, int d, int h)
{
    int j;
    int dptr = ddd->dptr[d];

    for (j = 0; (j < ddd->nhydro[dptr]); j++)
    {
        if (ddd->hydro[dptr][j] == h)
        {
            ddd->nhbonds[dptr][j]++;
            break;
        }
    }
    if (j == ddd->nhydro[dptr])
    {
        gmx_fatal(FARGS, "No such hydrogen %d on donor %d\n", h+1, d+1);
    }
}

static int _acceptor_index(t_acceptors *a, int grp, int i,
                           const char *file, int line)
{
    int ai = a->aptr[i];

    if (a->grp[ai] != grp)
    {
        if (debug && bDebug)
        {
            fprintf(debug, "Acc. group inconsist.. grp[%d] = %d, grp = %d (%s, %d)\n",
                    ai, a->grp[ai], grp, file, line);
        }
        return NOTSET;
    }
    else
    {
        return ai;
    }
}
#define acceptor_index(a, grp, i) _acceptor_index(a, grp, i, __FILE__, __LINE__)

static int _donor_index(t_donors *d, int grp, int i, const char *file, int line)
{
    int di = d->dptr[i];

    if (di == NOTSET)
    {
        return NOTSET;
    }

    if (d->grp[di] != grp)
    {
        if (debug && bDebug)
        {
            fprintf(debug, "Don. group inconsist.. grp[%d] = %d, grp = %d (%s, %d)\n",
                    di, d->grp[di], grp, file, line);
        }
        return NOTSET;
    }
    else
    {
        return di;
    }
}
#define donor_index(d, grp, i) _donor_index(d, grp, i, __FILE__, __LINE__)

static gmx_bool isInterchangable(t_hbdata *hb, int d, int a, int grpa, int grpd)
{
    /* g_hbond doesn't allow overlapping groups */
    if (grpa != grpd)
    {
        return FALSE;
    }
    return
        donor_index(&hb->d, grpd, a) != NOTSET
        && acceptor_index(&hb->a, grpa, d) != NOTSET;
}


static void add_hbond(t_hbdata *hb, int d, int a, int h, int grpd, int grpa,
                      int frame, gmx_bool bMerge, int ihb, gmx_bool bContact)
{
    int      k, id, ia, hh;
    gmx_bool daSwap = FALSE;

    if ((id = hb->d.dptr[d]) == NOTSET)
    {
        gmx_fatal(FARGS, "No donor atom %d", d+1);
    }
    else if (grpd != hb->d.grp[id])
    {
        gmx_fatal(FARGS, "Inconsistent donor groups, %d instead of %d, atom %d",
                  grpd, hb->d.grp[id], d+1);
    }
    if ((ia = hb->a.aptr[a]) == NOTSET)
    {
        gmx_fatal(FARGS, "No acceptor atom %d", a+1);
    }
    else if (grpa != hb->a.grp[ia])
    {
        gmx_fatal(FARGS, "Inconsistent acceptor groups, %d instead of %d, atom %d",
                  grpa, hb->a.grp[ia], a+1);
    }

    if (bMerge)
    {

        if (isInterchangable(hb, d, a, grpd, grpa) && d > a)
        /* Then swap identity so that the id of d is lower then that of a.
         *
         * This should really be redundant by now, as is_hbond() now ought to return
         * hbNo in the cases where this conditional is TRUE. */
        {
            daSwap = TRUE;
            k      = d;
            d      = a;
            a      = k;

            /* Now repeat donor/acc check. */
            if ((id = hb->d.dptr[d]) == NOTSET)
            {
                gmx_fatal(FARGS, "No donor atom %d", d+1);
            }
            else if (grpd != hb->d.grp[id])
            {
                gmx_fatal(FARGS, "Inconsistent donor groups, %d instead of %d, atom %d",
                          grpd, hb->d.grp[id], d+1);
            }
            if ((ia = hb->a.aptr[a]) == NOTSET)
            {
                gmx_fatal(FARGS, "No acceptor atom %d", a+1);
            }
            else if (grpa != hb->a.grp[ia])
            {
                gmx_fatal(FARGS, "Inconsistent acceptor groups, %d instead of %d, atom %d",
                          grpa, hb->a.grp[ia], a+1);
            }
        }
    }

    if (hb->hbmap)
    {
        /* Loop over hydrogens to find which hydrogen is in this particular HB */
        if ((ihb == hbHB) && !bMerge && !bContact)
        {
            for (k = 0; (k < hb->d.nhydro[id]); k++)
            {
                if (hb->d.hydro[id][k] == h)
                {
                    break;
                }
            }
            if (k == hb->d.nhydro[id])
            {
                gmx_fatal(FARGS, "Donor %d does not have hydrogen %d (a = %d)",
                          d+1, h+1, a+1);
            }
        }
        else
        {
            k = 0;
        }

        if (hb->bHBmap)
        {

#pragma omp critical
            {
                try
                {
                    if (hb->hbmap[id][ia] == nullptr)
                    {
                        snew(hb->hbmap[id][ia], 1);
                        snew(hb->hbmap[id][ia]->h, hb->maxhydro);
                        snew(hb->hbmap[id][ia]->g, hb->maxhydro);
                    }
                    add_ff(hb, id, k, ia, frame, ihb);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
        }

        /* Strange construction with frame >=0 is a relic from old code
         * for selected hbond analysis. It may be necessary again if that
         * is made to work again.
         */
        if (frame >= 0)
        {
            hh = hb->hbmap[id][ia]->history[k];
            if (ihb == hbHB)
            {
                hb->nhb[frame]++;
                if (!(ISHB(hh)))
                {
                    hb->hbmap[id][ia]->history[k] = hh | 2;
                    hb->nrhb++;
                }
            }
            else
            {
                if (ihb == hbDist)
                {
                    hb->ndist[frame]++;
                    if (!(ISDIST(hh)))
                    {
                        hb->hbmap[id][ia]->history[k] = hh | 1;
                        hb->nrdist++;
                    }
                }
            }
        }
    }
    else
    {
        if (frame >= 0)
        {
            if (ihb == hbHB)
            {
                hb->nhb[frame]++;
            }
            else
            {
                if (ihb == hbDist)
                {
                    hb->ndist[frame]++;
                }
            }
        }
    }
    if (bMerge && daSwap)
    {
        h = hb->d.hydro[id][0];
    }
    /* Increment number if HBonds per H */
    if (ihb == hbHB && !bContact)
    {
        inc_nhbonds(&(hb->d), d, h);
    }
}

static char *mkatomname(const t_atoms *atoms, int i)
{
    static char buf[32];
    int         rnr;

    rnr = atoms->atom[i].resind;
    sprintf(buf, "%4s%d%-4s",
            *atoms->resinfo[rnr].name, atoms->resinfo[rnr].nr, *atoms->atomname[i]);

    return buf;
}

static void gen_datable(int *index, int isize, unsigned char *datable, int natoms)
{
    /* Generates table of all atoms and sets the ingroup bit for atoms in index[] */
    int i;

    for (i = 0; i < isize; i++)
    {
        if (index[i] >= natoms)
        {
            gmx_fatal(FARGS, "Atom has index %d larger than number of atoms %d.", index[i], natoms);
        }
        datable[index[i]] |= c_inGroupMask;
    }
}

static void clear_datable_grp(unsigned char *datable, int size)
{
    /* Clears group information from the table */
    int        i;
    if (size > 0)
    {
        for (i = 0; i < size; i++)
        {
            datable[i] &= ~c_inGroupMask;
        }
    }
}

static void add_acc(t_acceptors *a, int ia, int grp)
{
    if (a->nra >= a->max_nra)
    {
        a->max_nra += 16;
        srenew(a->acc, a->max_nra);
        srenew(a->grp, a->max_nra);
    }
    a->grp[a->nra]   = grp;
    a->acc[a->nra++] = ia;
}

static void search_acceptors(const t_topology *top, int isize,
                             const int *index, t_acceptors *a, int grp,
                             gmx_bool bNitAcc,
                             gmx_bool bContact, gmx_bool bDoIt, unsigned char *datable)
{
    int i, n;

    if (bDoIt)
    {
        for (i = 0; (i < isize); i++)
        {
            n = index[i];
            if ((bContact ||
                 (((*top->atoms.atomname[n])[0] == 'O') ||
                  (bNitAcc && ((*top->atoms.atomname[n])[0] == 'N')))) &&
                ISINGRP(datable[n]))
            {
                datable[n] |= c_acceptorMask;
                add_acc(a, n, grp);
            }
        }
    }
    snew(a->aptr, top->atoms.nr);
    for (i = 0; (i < top->atoms.nr); i++)
    {
        a->aptr[i] = NOTSET;
    }
    for (i = 0; (i < a->nra); i++)
    {
        a->aptr[a->acc[i]] = i;
    }
}

static void add_h2d(int id, int ih, t_donors *ddd)
{
    int i;

    for (i = 0; (i < ddd->nhydro[id]); i++)
    {
        if (ddd->hydro[id][i] == ih)
        {
            printf("Hm. This isn't the first time I found this donor (%d,%d)\n",
                   ddd->don[id], ih);
            break;
        }
    }
    if (i == ddd->nhydro[id])
    {
        if (ddd->nhydro[id] >= MAXHYDRO)
        {
            gmx_fatal(FARGS, "Donor %d has more than %d hydrogens!",
                      ddd->don[id], MAXHYDRO);
        }
        ddd->hydro[id][i] = ih;
        ddd->nhydro[id]++;
    }
}

static void add_dh(t_donors *ddd, int id, int ih, int grp, unsigned char *datable)
{
    int i;

    if (!datable || ISDON(datable[id]))
    {
        if (ddd->dptr[id] == NOTSET)   /* New donor */
        {
            i             = ddd->nrd;
            ddd->dptr[id] = i;
        }
        else
        {
            i = ddd->dptr[id];
        }

        if (i == ddd->nrd)
        {
            if (ddd->nrd >= ddd->max_nrd)
            {
                ddd->max_nrd += 128;
                srenew(ddd->don, ddd->max_nrd);
                srenew(ddd->nhydro, ddd->max_nrd);
                srenew(ddd->hydro, ddd->max_nrd);
                srenew(ddd->nhbonds, ddd->max_nrd);
                srenew(ddd->grp, ddd->max_nrd);
            }
            ddd->don[ddd->nrd]    = id;
            ddd->nhydro[ddd->nrd] = 0;
            ddd->grp[ddd->nrd]    = grp;
            ddd->nrd++;
        }
        else
        {
            ddd->don[i] = id;
        }
        add_h2d(i, ih, ddd);
    }
    else
    if (datable)
    {
        printf("Warning: Atom %d is not in the d/a-table!\n", id);
    }
}

static void search_donors(const t_topology *top, int isize, const int *index,
                          t_donors *ddd, int grp, gmx_bool bContact, gmx_bool bDoIt,
                          unsigned char *datable)
{
    int            i, j;
    t_functype     func_type;
    int            nr1, nr2, nr3;

    if (!ddd->dptr)
    {
        snew(ddd->dptr, top->atoms.nr);
        for (i = 0; (i < top->atoms.nr); i++)
        {
            ddd->dptr[i] = NOTSET;
        }
    }

    if (bContact)
    {
        if (bDoIt)
        {
            for (i = 0; (i < isize); i++)
            {
                datable[index[i]] |= c_donorMask;
                add_dh(ddd, index[i], -1, grp, datable);
            }
        }
    }
    else
    {
        for (func_type = 0; (func_type < F_NRE); func_type++)
        {
            const t_ilist *interaction = &(top->idef.il[func_type]);
            if (func_type == F_POSRES || func_type == F_FBPOSRES)
            {
                /* The ilist looks strange for posre. Bug in grompp?
                 * We don't need posre interactions for hbonds anyway.*/
                continue;
            }
            for (i = 0; i < interaction->nr;
                 i += interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1)
            {
                /* next function */
                if (func_type != top->idef.functype[interaction->iatoms[i]])
                {
                    fprintf(stderr, "Error in func_type %s",
                            interaction_function[func_type].longname);
                    continue;
                }

                /* check out this functype */
                if (func_type == F_SETTLE)
                {
                    nr1 = interaction->iatoms[i+1];
                    nr2 = interaction->iatoms[i+2];
                    nr3 = interaction->iatoms[i+3];

                    if (ISINGRP(datable[nr1]))
                    {
                        if (ISINGRP(datable[nr2]))
                        {
                            datable[nr1] |= c_donorMask;
                            add_dh(ddd, nr1, nr1+1, grp, datable);
                        }
                        if (ISINGRP(datable[nr3]))
                        {
                            datable[nr1] |= c_donorMask;
                            add_dh(ddd, nr1, nr1+2, grp, datable);
                        }
                    }
                }
                else if (IS_CHEMBOND(func_type))
                {
                    for (j = 0; j < 2; j++)
                    {
                        nr1 = interaction->iatoms[i+1+j];
                        nr2 = interaction->iatoms[i+2-j];
                        if ((*top->atoms.atomname[nr1][0] == 'H') &&
                            ((*top->atoms.atomname[nr2][0] == 'O') ||
                             (*top->atoms.atomname[nr2][0] == 'N')) &&
                            ISINGRP(datable[nr1]) && ISINGRP(datable[nr2]))
                        {
                            datable[nr2] |= c_donorMask;
                            add_dh(ddd, nr2, nr1, grp, datable);
                        }
                    }
                }
            }
        }
#ifdef SAFEVSITES
        for (func_type = 0; func_type < F_NRE; func_type++)
        {
            interaction = &top->idef.il[func_type];
            for (i = 0; i < interaction->nr;
                 i += interaction_function[top->idef.functype[interaction->iatoms[i]]].nratoms+1)
            {
                /* next function */
                if (func_type != top->idef.functype[interaction->iatoms[i]])
                {
                    gmx_incons("function type in search_donors");
                }

                if (interaction_function[func_type].flags & IF_VSITE)
                {
                    nr1 = interaction->iatoms[i+1];
                    if (*top->atoms.atomname[nr1][0]  == 'H')
                    {
                        nr2  = nr1-1;
                        stop = FALSE;
                        while (!stop && ( *top->atoms.atomname[nr2][0] == 'H'))
                        {
                            if (nr2)
                            {
                                nr2--;
                            }
                            else
                            {
                                stop = TRUE;
                            }
                        }
                        if (!stop && ( ( *top->atoms.atomname[nr2][0] == 'O') ||
                                       ( *top->atoms.atomname[nr2][0] == 'N') ) &&
                            ISINGRP(datable[nr1]) && ISINGRP(datable[nr2]))
                        {
                            datable[nr2] |= c_donorMask;
                            add_dh(ddd, nr2, nr1, grp, datable);
                        }
                    }
                }
            }
        }
#endif
    }
}

static t_gridcell ***init_grid(gmx_bool bBox, rvec box[], real rcut, ivec ngrid)
{
    t_gridcell ***grid;
    int           i, y, z;

    if (bBox)
    {
        for (i = 0; i < DIM; i++)
        {
            ngrid[i] = static_cast<int>(box[i][i]/(1.2*rcut));
        }
    }

    if (!bBox || (ngrid[XX] < 3) || (ngrid[YY] < 3) || (ngrid[ZZ] < 3) )
    {
        for (i = 0; i < DIM; i++)
        {
            ngrid[i] = 1;
        }
    }
    else
    {
        printf("\nWill do grid-seach on %dx%dx%d grid, rcut=%g\n",
               ngrid[XX], ngrid[YY], ngrid[ZZ], rcut);
    }
    snew(grid, ngrid[ZZ]);
    for (z = 0; z < ngrid[ZZ]; z++)
    {
        snew((grid)[z], ngrid[YY]);
        for (y = 0; y < ngrid[YY]; y++)
        {
            snew((grid)[z][y], ngrid[XX]);
        }
    }
    return grid;
}

static void reset_nhbonds(t_donors *ddd)
{
    int i, j;

    for (i = 0; (i < ddd->nrd); i++)
    {
        for (j = 0; (j < MAXHH); j++)
        {
            ddd->nhbonds[i][j] = 0;
        }
    }
}

static void pbc_correct_gem(rvec dx, matrix box, rvec hbox);
static void pbc_in_gridbox(rvec dx, matrix box);

static void build_grid(t_hbdata *hb, rvec x[], rvec xshell,
                       gmx_bool bBox, matrix box, rvec hbox,
                       real rcut, real rshell,
                       ivec ngrid, t_gridcell ***grid)
{
    int         i, m, gr, xi, yi, zi, nr;
    int        *ad;
    ivec        grididx;
    rvec        invdelta, dshell;
    t_ncell    *newgrid;
    gmx_bool    bDoRshell, bInShell, bAcc;
    real        rshell2 = 0;
    int         gx, gy, gz;
    int         dum = -1;

    bDoRshell = (rshell > 0);
    rshell2   = gmx::square(rshell);
    bInShell  = TRUE;

#define DBB(x) if (debug && bDebug) fprintf(debug, "build_grid, line %d, %s = %d\n", __LINE__,#x, x)
    DBB(dum);
    for (m = 0; m < DIM; m++)
    {
        hbox[m] = box[m][m]*0.5;
        if (bBox)
        {
            invdelta[m] = ngrid[m]/box[m][m];
            if (1/invdelta[m] < rcut)
            {
                gmx_fatal(FARGS, "Your computational box has shrunk too much.\n"
                          "%s can not handle this situation, sorry.\n",
                          gmx::getProgramContext().displayName());
            }
        }
        else
        {
            invdelta[m] = 0;
        }
    }
    grididx[XX] = 0;
    grididx[YY] = 0;
    grididx[ZZ] = 0;
    DBB(dum);
    /* resetting atom counts */
    for (gr = 0; (gr < grNR); gr++)
    {
        for (zi = 0; zi < ngrid[ZZ]; zi++)
        {
            for (yi = 0; yi < ngrid[YY]; yi++)
            {
                for (xi = 0; xi < ngrid[XX]; xi++)
                {
                    grid[zi][yi][xi].d[gr].nr = 0;
                    grid[zi][yi][xi].a[gr].nr = 0;
                }
            }
        }
        DBB(dum);

        /* put atoms in grid cells */
        for (bAcc = FALSE; (bAcc <= TRUE); bAcc++)
        {
            if (bAcc)
            {
                nr = hb->a.nra;
                ad = hb->a.acc;
            }
            else
            {
                nr = hb->d.nrd;
                ad = hb->d.don;
            }
            DBB(bAcc);
            for (i = 0; (i < nr); i++)
            {
                /* check if we are inside the shell */
                /* if bDoRshell=FALSE then bInShell=TRUE always */
                DBB(i);
                if (bDoRshell)
                {
                    bInShell = TRUE;
                    rvec_sub(x[ad[i]], xshell, dshell);
                    if (bBox)
                    {
                        gmx_bool bDone = FALSE;
                        while (!bDone)
                        {
                            bDone = TRUE;
                            for (m = DIM-1; m >= 0 && bInShell; m--)
                            {
                                if (dshell[m] < -hbox[m])
                                {
                                    bDone = FALSE;
                                    rvec_inc(dshell, box[m]);
                                }
                                if (dshell[m] >= hbox[m])
                                {
                                    bDone      = FALSE;
                                    dshell[m] -= 2*hbox[m];
                                }
                            }
                        }
                        for (m = DIM-1; m >= 0 && bInShell; m--)
                        {
                            /* if we're outside the cube, we're outside the sphere also! */
                            if ( (dshell[m] > rshell) || (-dshell[m] > rshell) )
                            {
                                bInShell = FALSE;
                            }
                        }
                    }
                    /* if we're inside the cube, check if we're inside the sphere */
                    if (bInShell)
                    {
                        bInShell = norm2(dshell) < rshell2;
                    }
                }
                DBB(i);
                if (bInShell)
                {
                    if (bBox)
                    {
                        pbc_in_gridbox(x[ad[i]], box);

                        for (m = DIM-1; m >= 0; m--)
                        {   /* determine grid index of atom */
                            grididx[m] = static_cast<int>(x[ad[i]][m]*invdelta[m]);
                            grididx[m] = (grididx[m]+ngrid[m]) % ngrid[m];
                        }
                    }

                    gx = grididx[XX];
                    gy = grididx[YY];
                    gz = grididx[ZZ];
                    range_check(gx, 0, ngrid[XX]);
                    range_check(gy, 0, ngrid[YY]);
                    range_check(gz, 0, ngrid[ZZ]);
                    DBB(gx);
                    DBB(gy);
                    DBB(gz);
                    /* add atom to grid cell */
                    if (bAcc)
                    {
                        newgrid = &(grid[gz][gy][gx].a[gr]);
                    }
                    else
                    {
                        newgrid = &(grid[gz][gy][gx].d[gr]);
                    }
                    if (newgrid->nr >= newgrid->maxnr)
                    {
                        newgrid->maxnr += 10;
                        DBB(newgrid->maxnr);
                        srenew(newgrid->atoms, newgrid->maxnr);
                    }
                    DBB(newgrid->nr);
                    newgrid->atoms[newgrid->nr] = ad[i];
                    newgrid->nr++;
                }
            }
        }
    }
}

static void count_da_grid(ivec ngrid, t_gridcell ***grid, t_icell danr)
{
    int gr, xi, yi, zi;

    for (gr = 0; (gr < grNR); gr++)
    {
        danr[gr] = 0;
        for (zi = 0; zi < ngrid[ZZ]; zi++)
        {
            for (yi = 0; yi < ngrid[YY]; yi++)
            {
                for (xi = 0; xi < ngrid[XX]; xi++)
                {
                    danr[gr] += grid[zi][yi][xi].d[gr].nr;
                }
            }
        }
    }
}

/* The grid loop.
 * Without a box, the grid is 1x1x1, so all loops are 1 long.
 * With a rectangular box (bTric==FALSE) all loops are 3 long.
 * With a triclinic box all loops are 3 long, except when a cell is
 * located next to one of the box edges which is not parallel to the
 * x/y-plane, in that case all cells in a line or layer are searched.
 * This could be implemented slightly more efficient, but the code
 * would get much more complicated.
 */
static gmx_inline gmx_bool grid_loop_begin(int n, int x, gmx_bool bTric, gmx_bool bEdge)
{
    return ((n == 1) ? x : bTric && bEdge ? 0     : (x-1));
}
static gmx_inline gmx_bool grid_loop_end(int n, int x, gmx_bool bTric, gmx_bool bEdge)
{
    return ((n == 1) ? x : bTric && bEdge ? (n-1) : (x+1));
}
static gmx_inline int grid_mod(int j, int n)
{
    return (j+n) % (n);
}

static void dump_grid(FILE *fp, ivec ngrid, t_gridcell ***grid)
{
    int gr, x, y, z, sum[grNR];

    fprintf(fp, "grid %dx%dx%d\n", ngrid[XX], ngrid[YY], ngrid[ZZ]);
    for (gr = 0; gr < grNR; gr++)
    {
        sum[gr] = 0;
        fprintf(fp, "GROUP %d (%s)\n", gr, grpnames[gr]);
        for (z = 0; z < ngrid[ZZ]; z += 2)
        {
            fprintf(fp, "Z=%d,%d\n", z, z+1);
            for (y = 0; y < ngrid[YY]; y++)
            {
                for (x = 0; x < ngrid[XX]; x++)
                {
                    fprintf(fp, "%3d", grid[x][y][z].d[gr].nr);
                    sum[gr] += grid[z][y][x].d[gr].nr;
                    fprintf(fp, "%3d", grid[x][y][z].a[gr].nr);
                    sum[gr] += grid[z][y][x].a[gr].nr;

                }
                fprintf(fp, " | ");
                if ( (z+1) < ngrid[ZZ])
                {
                    for (x = 0; x < ngrid[XX]; x++)
                    {
                        fprintf(fp, "%3d", grid[z+1][y][x].d[gr].nr);
                        sum[gr] += grid[z+1][y][x].d[gr].nr;
                        fprintf(fp, "%3d", grid[z+1][y][x].a[gr].nr);
                        sum[gr] += grid[z+1][y][x].a[gr].nr;
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }
    fprintf(fp, "TOTALS:");
    for (gr = 0; gr < grNR; gr++)
    {
        fprintf(fp, " %d=%d", gr, sum[gr]);
    }
    fprintf(fp, "\n");
}

/* New GMX record! 5 * in a row. Congratulations!
 * Sorry, only four left.
 */
static void free_grid(ivec ngrid, t_gridcell ****grid)
{
    int           y, z;
    t_gridcell ***g = *grid;

    for (z = 0; z < ngrid[ZZ]; z++)
    {
        for (y = 0; y < ngrid[YY]; y++)
        {
            sfree(g[z][y]);
        }
        sfree(g[z]);
    }
    sfree(g);
    g = nullptr;
}

static void pbc_correct_gem(rvec dx, matrix box, rvec hbox)
{
    int      m;
    gmx_bool bDone = FALSE;
    while (!bDone)
    {
        bDone = TRUE;
        for (m = DIM-1; m >= 0; m--)
        {
            if (dx[m] < -hbox[m])
            {
                bDone = FALSE;
                rvec_inc(dx, box[m]);
            }
            if (dx[m] >= hbox[m])
            {
                bDone = FALSE;
                rvec_dec(dx, box[m]);
            }
        }
    }
}

static void pbc_in_gridbox(rvec dx, matrix box)
{
    int      m;
    gmx_bool bDone = FALSE;
    while (!bDone)
    {
        bDone = TRUE;
        for (m = DIM-1; m >= 0; m--)
        {
            if (dx[m] < 0)
            {
                bDone = FALSE;
                rvec_inc(dx, box[m]);
            }
            if (dx[m] >= box[m][m])
            {
                bDone = FALSE;
                rvec_dec(dx, box[m]);
            }
        }
    }
}

/* Added argument r2cut, changed contact and implemented
 * use of second cut-off.
 * - Erik Marklund, June 29, 2006
 */
static int is_hbond(t_hbdata *hb, int grpd, int grpa, int d, int a,
                    real rcut, real r2cut, real ccut,
                    rvec x[], gmx_bool bBox, matrix box, rvec hbox,
                    real *d_ha, real *ang, gmx_bool bDA, int *hhh,
                    gmx_bool bContact, gmx_bool bMerge)
{
    int      h, hh, id;
    rvec     r_da, r_ha, r_dh;
    real     rc2, r2c2, rda2, rha2, ca;
    gmx_bool HAinrange = FALSE; /* If !bDA. Needed for returning hbDist in a correct way. */
    gmx_bool daSwap    = FALSE;

    if (d == a)
    {
        return hbNo;
    }

    if (((id = donor_index(&hb->d, grpd, d)) == NOTSET) ||
        (acceptor_index(&hb->a, grpa, a) == NOTSET))
    {
        return hbNo;
    }

    rc2  = rcut*rcut;
    r2c2 = r2cut*r2cut;

    rvec_sub(x[d], x[a], r_da);
    /* Insert projection code here */

    if (bMerge && d > a && isInterchangable(hb, d, a, grpd, grpa))
    {
        /* Then this hbond/contact will be found again, or it has already been found. */
        /*return hbNo;*/
    }
    if (bBox)
    {
        if (d > a && bMerge && isInterchangable(hb, d, a, grpd, grpa)) /* acceptor is also a donor and vice versa? */
        {                                                              /* return hbNo; */
            daSwap = TRUE;                                             /* If so, then their history should be filed with donor and acceptor swapped. */
        }
        pbc_correct_gem(r_da, box, hbox);
    }
    rda2 = iprod(r_da, r_da);

    if (bContact)
    {
        if (daSwap && grpa == grpd)
        {
            return hbNo;
        }
        if (rda2 <= rc2)
        {
            return hbHB;
        }
        else if (rda2 < r2c2)
        {
            return hbDist;
        }
        else
        {
            return hbNo;
        }
    }
    *hhh = NOTSET;

    if (bDA && (rda2 > rc2))
    {
        return hbNo;
    }

    for (h = 0; (h < hb->d.nhydro[id]); h++)
    {
        hh   = hb->d.hydro[id][h];
        rha2 = rc2+1;
        if (!bDA)
        {
            rvec_sub(x[hh], x[a], r_ha);
            if (bBox)
            {
                pbc_correct_gem(r_ha, box, hbox);
            }
            rha2 = iprod(r_ha, r_ha);
        }

        if (bDA || (rha2 <= rc2))
        {
            rvec_sub(x[d], x[hh], r_dh);
            if (bBox)
            {
                pbc_correct_gem(r_dh, box, hbox);
            }

            if (!bDA)
            {
                HAinrange = TRUE;
            }
            ca = cos_angle(r_dh, r_da);
            /* if angle is smaller, cos is larger */
            if (ca >= ccut)
            {
                *hhh  = hh;
                *d_ha = std::sqrt(bDA ? rda2 : rha2);
                *ang  = std::acos(ca);
                return hbHB;
            }
        }
    }
    if (bDA || HAinrange)
    {
        return hbDist;
    }
    else
    {
        return hbNo;
    }
}

/* Merging is now done on the fly, so do_merge is most likely obsolete now.
 * Will do some more testing before removing the function entirely.
 * - Erik Marklund, MAY 10 2010 */
static void do_merge(t_hbdata *hb, int ntmp,
                     unsigned int htmp[], unsigned int gtmp[],
                     t_hbond *hb0, t_hbond *hb1)
{
    /* Here we need to make sure we're treating periodicity in
     * the right way for the geminate recombination kinetics. */

    int       m, mm, n00, n01, nn0, nnframes;

    /* Decide where to start from when merging */
    n00      = hb0->n0;
    n01      = hb1->n0;
    nn0      = std::min(n00, n01);
    nnframes = std::max(n00 + hb0->nframes, n01 + hb1->nframes) - nn0;
    /* Initiate tmp arrays */
    for (m = 0; (m < ntmp); m++)
    {
        htmp[m] = 0;
        gtmp[m] = 0;
    }
    /* Fill tmp arrays with values due to first HB */
    /* Once again '<' had to be replaced with '<='
       to catch the last frame in which the hbond
       appears.
       - Erik Marklund, June 1, 2006 */
    for (m = 0; (m <= hb0->nframes); m++)
    {
        mm       = m+n00-nn0;
        htmp[mm] = is_hb(hb0->h[0], m);
    }
    for (m = 0; (m <= hb0->nframes); m++)
    {
        mm       = m+n00-nn0;
        gtmp[mm] = is_hb(hb0->g[0], m);
    }
    /* Next HB */
    for (m = 0; (m <= hb1->nframes); m++)
    {
        mm       = m+n01-nn0;
        htmp[mm] = htmp[mm] || is_hb(hb1->h[0], m);
        gtmp[mm] = gtmp[mm] || is_hb(hb1->g[0], m);
    }
    /* Reallocate target array */
    if (nnframes > hb0->maxframes)
    {
        srenew(hb0->h[0], 4+nnframes/hb->wordlen);
        srenew(hb0->g[0], 4+nnframes/hb->wordlen);
    }

    /* Copy temp array to target array */
    for (m = 0; (m <= nnframes); m++)
    {
        _set_hb(hb0->h[0], m, htmp[m]);
        _set_hb(hb0->g[0], m, gtmp[m]);
    }

    /* Set scalar variables */
    hb0->n0        = nn0;
    hb0->maxframes = nnframes;
}

static void merge_hb(t_hbdata *hb, gmx_bool bTwo, gmx_bool bContact)
{
    int           i, inrnew, indnew, j, ii, jj, id, ia, ntmp;
    unsigned int *htmp, *gtmp;
    t_hbond      *hb0, *hb1;

    inrnew = hb->nrhb;
    indnew = hb->nrdist;

    /* Check whether donors are also acceptors */
    printf("Merging hbonds with Acceptor and Donor swapped\n");

    ntmp = 2*hb->max_frames;
    snew(gtmp, ntmp);
    snew(htmp, ntmp);
    for (i = 0; (i < hb->d.nrd); i++)
    {
        fprintf(stderr, "\r%d/%d", i+1, hb->d.nrd);
        fflush(stderr);
        id = hb->d.don[i];
        ii = hb->a.aptr[id];
        for (j = 0; (j < hb->a.nra); j++)
        {
            ia = hb->a.acc[j];
            jj = hb->d.dptr[ia];
            if ((id != ia) && (ii != NOTSET) && (jj != NOTSET) &&
                (!bTwo || (hb->d.grp[i] != hb->a.grp[j])))
            {
                hb0 = hb->hbmap[i][j];
                hb1 = hb->hbmap[jj][ii];
                if (hb0 && hb1 && ISHB(hb0->history[0]) && ISHB(hb1->history[0]))
                {
                    do_merge(hb, ntmp, htmp, gtmp, hb0, hb1);
                    if (ISHB(hb1->history[0]))
                    {
                        inrnew--;
                    }
                    else if (ISDIST(hb1->history[0]))
                    {
                        indnew--;
                    }
                    else
                    if (bContact)
                    {
                        gmx_incons("No contact history");
                    }
                    else
                    {
                        gmx_incons("Neither hydrogen bond nor distance");
                    }
                    sfree(hb1->h[0]);
                    sfree(hb1->g[0]);
                    hb1->h[0]       = nullptr;
                    hb1->g[0]       = nullptr;
                    hb1->history[0] = hbNo;
                }
            }
        }
    }
    fprintf(stderr, "\n");
    printf("- Reduced number of hbonds from %d to %d\n", hb->nrhb, inrnew);
    printf("- Reduced number of distances from %d to %d\n", hb->nrdist, indnew);
    hb->nrhb   = inrnew;
    hb->nrdist = indnew;
    sfree(gtmp);
    sfree(htmp);
}

static void do_nhb_dist(FILE *fp, t_hbdata *hb, real t)
{
    int  i, j, k, n_bound[MAXHH], nbtot;

    /* Set array to 0 */
    for (k = 0; (k < MAXHH); k++)
    {
        n_bound[k] = 0;
    }
    /* Loop over possible donors */
    for (i = 0; (i < hb->d.nrd); i++)
    {
        for (j = 0; (j < hb->d.nhydro[i]); j++)
        {
            n_bound[hb->d.nhbonds[i][j]]++;
        }
    }
    fprintf(fp, "%12.5e", t);
    nbtot = 0;
    for (k = 0; (k < MAXHH); k++)
    {
        fprintf(fp, "  %8d", n_bound[k]);
        nbtot += n_bound[k]*k;
    }
    fprintf(fp, "  %8d\n", nbtot);
}

static void do_hblife(const char *fn, t_hbdata *hb, gmx_bool bMerge, gmx_bool bContact,
                      const gmx_output_env_t *oenv)
{
    FILE          *fp;
    const char    *leg[] = { "p(t)", "t p(t)" };
    int           *histo;
    int            i, j, j0, k, m, nh, ihb, ohb, nhydro, ndump = 0;
    int            nframes = hb->nframes;
    unsigned int **h;
    real           t, x1, dt;
    double         sum, integral;
    t_hbond       *hbh;

    snew(h, hb->maxhydro);
    snew(histo, nframes+1);
    /* Total number of hbonds analyzed here */
    for (i = 0; (i < hb->d.nrd); i++)
    {
        for (k = 0; (k < hb->a.nra); k++)
        {
            hbh = hb->hbmap[i][k];
            if (hbh)
            {
                if (bMerge)
                {
                    if (hbh->h[0])
                    {
                        h[0]   = hbh->h[0];
                        nhydro = 1;
                    }
                    else
                    {
                        nhydro = 0;
                    }
                }
                else
                {
                    nhydro = 0;
                    for (m = 0; (m < hb->maxhydro); m++)
                    {
                        if (hbh->h[m])
                        {
                            h[nhydro++] = bContact ? hbh->g[m] : hbh->h[m];
                        }
                    }
                }
                for (nh = 0; (nh < nhydro); nh++)
                {
                    ohb = 0;
                    j0  = 0;

                    for (j = 0; (j <= hbh->nframes); j++)
                    {
                        ihb      = is_hb(h[nh], j);
                        if (debug && (ndump < 10))
                        {
                            fprintf(debug, "%5d  %5d\n", j, ihb);
                        }
                        if (ihb != ohb)
                        {
                            if (ihb)
                            {
                                j0 = j;
                            }
                            else
                            {
                                histo[j-j0]++;
                            }
                            ohb = ihb;
                        }
                    }
                    ndump++;
                }
            }
        }
    }
    fprintf(stderr, "\n");
    if (bContact)
    {
        fp = xvgropen(fn, "Uninterrupted contact lifetime", output_env_get_xvgr_tlabel(oenv), "()", oenv);
    }
    else
    {
        fp = xvgropen(fn, "Uninterrupted hydrogen bond lifetime", output_env_get_xvgr_tlabel(oenv), "()",
                      oenv);
    }

    xvgr_legend(fp, asize(leg), leg, oenv);
    j0 = nframes-1;
    while ((j0 > 0) && (histo[j0] == 0))
    {
        j0--;
    }
    sum = 0;
    for (i = 0; (i <= j0); i++)
    {
        sum += histo[i];
    }
    dt       = hb->time[1]-hb->time[0];
    sum      = dt*sum;
    integral = 0;
    for (i = 1; (i <= j0); i++)
    {
        t  = hb->time[i] - hb->time[0] - 0.5*dt;
        x1 = t*histo[i]/sum;
        fprintf(fp, "%8.3f  %10.3e  %10.3e\n", t, histo[i]/sum, x1);
        integral += x1;
    }
    integral *= dt;
    xvgrclose(fp);
    printf("%s lifetime = %.2f ps\n", bContact ? "Contact" : "HB", integral);
    printf("Note that the lifetime obtained in this manner is close to useless\n");
    printf("Use the -ac option instead and check the Forward lifetime\n");
    please_cite(stdout, "Spoel2006b");
    sfree(h);
    sfree(histo);
}

static void dump_ac(t_hbdata *hb, gmx_bool oneHB, int nDump)
{
    FILE     *fp;
    int       i, j, k, m, nd, ihb, idist;
    int       nframes = hb->nframes;
    gmx_bool  bPrint;
    t_hbond  *hbh;

    if (nDump <= 0)
    {
        return;
    }
    fp = gmx_ffopen("debug-ac.xvg", "w");
    for (j = 0; (j < nframes); j++)
    {
        fprintf(fp, "%10.3f", hb->time[j]);
        for (i = nd = 0; (i < hb->d.nrd) && (nd < nDump); i++)
        {
            for (k = 0; (k < hb->a.nra) && (nd < nDump); k++)
            {
                bPrint = FALSE;
                ihb    = idist = 0;
                hbh    = hb->hbmap[i][k];
                if (oneHB)
                {
                    if (hbh->h[0])
                    {
                        ihb    = is_hb(hbh->h[0], j);
                        idist  = is_hb(hbh->g[0], j);
                        bPrint = TRUE;
                    }
                }
                else
                {
                    for (m = 0; (m < hb->maxhydro) && !ihb; m++)
                    {
                        ihb   = ihb   || ((hbh->h[m]) && is_hb(hbh->h[m], j));
                        idist = idist || ((hbh->g[m]) && is_hb(hbh->g[m], j));
                    }
                    /* This is not correct! */
                    /* What isn't correct? -Erik M */
                    bPrint = TRUE;
                }
                if (bPrint)
                {
                    fprintf(fp, "  %1d-%1d", ihb, idist);
                    nd++;
                }
            }
        }
        fprintf(fp, "\n");
    }
    gmx_ffclose(fp);
}

static real calc_dg(real tau, real temp)
{
    real kbt;

    kbt = BOLTZ*temp;
    if (tau <= 0)
    {
        return -666;
    }
    else
    {
        return kbt*std::log(kbt*tau/PLANCK);
    }
}

typedef struct {
    int   n0, n1, nparams, ndelta;
    real  kkk[2];
    real *t, *ct, *nt, *kt, *sigma_ct, *sigma_nt, *sigma_kt;
} t_luzar;

static real compute_weighted_rates(int n, real t[], real ct[], real nt[],
                                   real kt[], real sigma_ct[], real sigma_nt[],
                                   real sigma_kt[], real *k, real *kp,
                                   real *sigma_k, real *sigma_kp,
                                   real fit_start)
{
#define NK 10
    int      i, j;
    t_luzar  tl;
    real     kkk = 0, kkp = 0, kk2 = 0, kp2 = 0, chi2;

    *sigma_k  = 0;
    *sigma_kp = 0;

    for (i = 0; (i < n); i++)
    {
        if (t[i] >= fit_start)
        {
            break;
        }
    }
    tl.n0       = i;
    tl.n1       = n;
    tl.nparams  = 2;
    tl.ndelta   = 1;
    tl.t        = t;
    tl.ct       = ct;
    tl.nt       = nt;
    tl.kt       = kt;
    tl.sigma_ct = sigma_ct;
    tl.sigma_nt = sigma_nt;
    tl.sigma_kt = sigma_kt;
    tl.kkk[0]   = *k;
    tl.kkk[1]   = *kp;

    chi2      = 0; /*optimize_luzar_parameters(debug, &tl, 1000, 1e-3); */
    *k        = tl.kkk[0];
    *kp       = tl.kkk[1] = *kp;
    tl.ndelta = NK;
    for (j = 0; (j < NK); j++)
    {
        kkk += tl.kkk[0];
        kkp += tl.kkk[1];
        kk2 += gmx::square(tl.kkk[0]);
        kp2 += gmx::square(tl.kkk[1]);
        tl.n0++;
    }
    *sigma_k  = std::sqrt(kk2/NK - gmx::square(kkk/NK));
    *sigma_kp = std::sqrt(kp2/NK - gmx::square(kkp/NK));

    return chi2;
}

void analyse_corr(int n, real t[], real ct[], real nt[], real kt[],
                  real sigma_ct[], real sigma_nt[], real sigma_kt[],
                  real fit_start, real temp)
{
    int        i0, i;
    real       k = 1, kp = 1, kow = 1;
    real       Q = 0, chi2, tau_hb, dtau, tau_rlx, e_1, sigma_k, sigma_kp, ddg;
    double     tmp, sn2 = 0, sc2 = 0, sk2 = 0, scn = 0, sck = 0, snk = 0;
    gmx_bool   bError = (sigma_ct != nullptr) && (sigma_nt != nullptr) && (sigma_kt != nullptr);

    for (i0 = 0; (i0 < n-2) && ((t[i0]-t[0]) < fit_start); i0++)
    {
        ;
    }
    if (i0 < n-2)
    {
        for (i = i0; (i < n); i++)
        {
            sc2 += gmx::square(ct[i]);
            sn2 += gmx::square(nt[i]);
            sk2 += gmx::square(kt[i]);
            sck += ct[i]*kt[i];
            snk += nt[i]*kt[i];
            scn += ct[i]*nt[i];
        }
        printf("Hydrogen bond thermodynamics at T = %g K\n", temp);
        tmp = (sn2*sc2-gmx::square(scn));
        if ((tmp > 0) && (sn2 > 0))
        {
            k    = (sn2*sck-scn*snk)/tmp;
            kp   = (k*scn-snk)/sn2;
            if (bError)
            {
                chi2 = 0;
                for (i = i0; (i < n); i++)
                {
                    chi2 += gmx::square(k*ct[i]-kp*nt[i]-kt[i]);
                }
                compute_weighted_rates(n, t, ct, nt, kt, sigma_ct, sigma_nt,
                                       sigma_kt, &k, &kp,
                                       &sigma_k, &sigma_kp, fit_start);
                Q   = 0; /* quality_of_fit(chi2, 2);*/
                ddg = BOLTZ*temp*sigma_k/k;
                printf("Fitting paramaters chi^2 = %10g, Quality of fit = %10g\n",
                       chi2, Q);
                printf("The Rate and Delta G are followed by an error estimate\n");
                printf("----------------------------------------------------------\n"
                       "Type      Rate (1/ps)  Sigma Time (ps)  DG (kJ/mol)  Sigma\n");
                printf("Forward    %10.3f %6.2f   %8.3f  %10.3f %6.2f\n",
                       k, sigma_k, 1/k, calc_dg(1/k, temp), ddg);
                ddg = BOLTZ*temp*sigma_kp/kp;
                printf("Backward   %10.3f %6.2f   %8.3f  %10.3f %6.2f\n",
                       kp, sigma_kp, 1/kp, calc_dg(1/kp, temp), ddg);
            }
            else
            {
                chi2 = 0;
                for (i = i0; (i < n); i++)
                {
                    chi2 += gmx::square(k*ct[i]-kp*nt[i]-kt[i]);
                }
                printf("Fitting parameters chi^2 = %10g\nQ = %10g\n",
                       chi2, Q);
                printf("--------------------------------------------------\n"
                       "Type      Rate (1/ps) Time (ps)  DG (kJ/mol)  Chi^2\n");
                printf("Forward    %10.3f   %8.3f  %10.3f  %10g\n",
                       k, 1/k, calc_dg(1/k, temp), chi2);
                printf("Backward   %10.3f   %8.3f  %10.3f\n",
                       kp, 1/kp, calc_dg(1/kp, temp));
            }
        }
        if (sc2 > 0)
        {
            kow  = 2*sck/sc2;
            printf("One-way    %10.3f   %s%8.3f  %10.3f\n",
                   kow, bError ? "       " : "", 1/kow, calc_dg(1/kow, temp));
        }
        else
        {
            printf(" - Numerical problems computing HB thermodynamics:\n"
                   "sc2 = %g  sn2 = %g  sk2 = %g sck = %g snk = %g scn = %g\n",
                   sc2, sn2, sk2, sck, snk, scn);
        }
        /* Determine integral of the correlation function */
        tau_hb = evaluate_integral(n, t, ct, nullptr, (t[n-1]-t[0])/2, &dtau);
        printf("Integral   %10.3f   %s%8.3f  %10.3f\n", 1/tau_hb,
               bError ? "       " : "", tau_hb, calc_dg(tau_hb, temp));
        e_1 = std::exp(-1.0);
        for (i = 0; (i < n-2); i++)
        {
            if ((ct[i] > e_1) && (ct[i+1] <= e_1))
            {
                break;
            }
        }
        if (i < n-2)
        {
            /* Determine tau_relax from linear interpolation */
            tau_rlx = t[i]-t[0] + (e_1-ct[i])*(t[i+1]-t[i])/(ct[i+1]-ct[i]);
            printf("Relaxation %10.3f   %8.3f  %s%10.3f\n", 1/tau_rlx,
                   tau_rlx, bError ? "       " : "",
                   calc_dg(tau_rlx, temp));
        }
    }
    else
    {
        printf("Correlation functions too short to compute thermodynamics\n");
    }
}

void compute_derivative(int nn, real x[], real y[], real dydx[])
{
    int j;

    /* Compute k(t) = dc(t)/dt */
    for (j = 1; (j < nn-1); j++)
    {
        dydx[j] = (y[j+1]-y[j-1])/(x[j+1]-x[j-1]);
    }
    /* Extrapolate endpoints */
    dydx[0]    = 2*dydx[1]   -  dydx[2];
    dydx[nn-1] = 2*dydx[nn-2] - dydx[nn-3];
}

static void normalizeACF(real *ct, real *gt, int nhb, int len)
{
    real ct_fac, gt_fac = 0;
    int  i;

    /* Xu and Berne use the same normalization constant */

    ct_fac = 1.0/ct[0];
    if (nhb != 0)
    {
        gt_fac = 1.0/nhb;
    }

    printf("Normalization for c(t) = %g for gh(t) = %g\n", ct_fac, gt_fac);
    for (i = 0; i < len; i++)
    {
        ct[i] *= ct_fac;
        if (gt != nullptr)
        {
            gt[i] *= gt_fac;
        }
    }
}

static void do_hbac(const char *fn, t_hbdata *hb,
                    int nDump, gmx_bool bMerge, gmx_bool bContact, real fit_start,
                    real temp, gmx_bool R2, const gmx_output_env_t *oenv,
                    int nThreads)
{
    FILE          *fp;
    int            i, j, k, m, ihb, idist, n2, nn;

    const char    *legLuzar[] = {
        "Ac\\sfin sys\\v{}\\z{}(t)",
        "Ac(t)",
        "Cc\\scontact,hb\\v{}\\z{}(t)",
        "-dAc\\sfs\\v{}\\z{}/dt"
    };
    gmx_bool       bNorm = FALSE;
    double         nhb   = 0;
    real          *rhbex = nullptr, *ht, *gt, *ght, *dght, *kt;
    real          *ct, tail, tail2, dtail, *cct;
    const real     tol     = 1e-3;
    int            nframes = hb->nframes;
    unsigned int **h       = nullptr, **g = nullptr;
    int            nh, nhbonds, nhydro;
    t_hbond       *hbh;
    int            acType;
    int           *dondata      = nullptr;

    enum {
        AC_NONE, AC_NN, AC_GEM, AC_LUZAR
    };

    const bool bOMP = GMX_OPENMP;

    printf("Doing autocorrelation ");

    acType = AC_LUZAR;
    printf("according to the theory of Luzar and Chandler.\n");
    fflush(stdout);

    /* build hbexist matrix in reals for autocorr */
    /* Allocate memory for computing ACF (rhbex) and aggregating the ACF (ct) */
    n2 = 1;
    while (n2 < nframes)
    {
        n2 *= 2;
    }

    nn = nframes/2;

    if (acType != AC_NN || bOMP)
    {
        snew(h, hb->maxhydro);
        snew(g, hb->maxhydro);
    }

    /* Dump hbonds for debugging */
    dump_ac(hb, bMerge || bContact, nDump);

    /* Total number of hbonds analyzed here */
    nhbonds = 0;

    if (acType != AC_LUZAR && bOMP)
    {
        nThreads = std::min((nThreads <= 0) ? INT_MAX : nThreads, gmx_omp_get_max_threads());

        gmx_omp_set_num_threads(nThreads);
        snew(dondata, nThreads);
        for (i = 0; i < nThreads; i++)
        {
            dondata[i] = -1;
        }
        printf("ACF calculations parallelized with OpenMP using %i threads.\n"
               "Expect close to linear scaling over this donor-loop.\n", nThreads);
        fflush(stdout);
        fprintf(stderr, "Donors: [thread no]\n");
        {
            char tmpstr[7];
            for (i = 0; i < nThreads; i++)
            {
                snprintf(tmpstr, 7, "[%i]", i);
                fprintf(stderr, "%-7s", tmpstr);
            }
        }
        fprintf(stderr, "\n");
    }


    /* Build the ACF */
    snew(rhbex, 2*n2);
    snew(ct, 2*n2);
    snew(gt, 2*n2);
    snew(ht, 2*n2);
    snew(ght, 2*n2);
    snew(dght, 2*n2);

    snew(kt, nn);
    snew(cct, nn);

    for (i = 0; (i < hb->d.nrd); i++)
    {
        for (k = 0; (k < hb->a.nra); k++)
        {
            nhydro = 0;
            hbh    = hb->hbmap[i][k];

            if (hbh)
            {
                if (bMerge || bContact)
                {
                    if (ISHB(hbh->history[0]))
                    {
                        h[0]   = hbh->h[0];
                        g[0]   = hbh->g[0];
                        nhydro = 1;
                    }
                }
                else
                {
                    for (m = 0; (m < hb->maxhydro); m++)
                    {
                        if (bContact ? ISDIST(hbh->history[m]) : ISHB(hbh->history[m]))
                        {
                            g[nhydro] = hbh->g[m];
                            h[nhydro] = hbh->h[m];
                            nhydro++;
                        }
                    }
                }

                int nf = hbh->nframes;
                for (nh = 0; (nh < nhydro); nh++)
                {
                    int nrint = bContact ? hb->nrdist : hb->nrhb;
                    if ((((nhbonds+1) % 10) == 0) || (nhbonds+1 == nrint))
                    {
                        fprintf(stderr, "\rACF %d/%d", nhbonds+1, nrint);
                        fflush(stderr);
                    }
                    nhbonds++;
                    for (j = 0; (j < nframes); j++)
                    {
                        if (j <= nf)
                        {
                            ihb   = is_hb(h[nh], j);
                            idist = is_hb(g[nh], j);
                        }
                        else
                        {
                            ihb = idist = 0;
                        }
                        rhbex[j] = ihb;
                        /* For contacts: if a second cut-off is provided, use it,
                         * otherwise use g(t) = 1-h(t) */
                        if (!R2 && bContact)
                        {
                            gt[j]  = 1-ihb;
                        }
                        else
                        {
                            gt[j]  = idist*(1-ihb);
                        }
                        ht[j]    = rhbex[j];
                        nhb     += ihb;
                    }

                    /* The autocorrelation function is normalized after summation only */
                    low_do_autocorr(nullptr, oenv, nullptr, nframes, 1, -1, &rhbex, hb->time[1]-hb->time[0],
                                    eacNormal, 1, FALSE, bNorm, FALSE, 0, -1, 0);

                    /* Cross correlation analysis for thermodynamics */
                    for (j = nframes; (j < n2); j++)
                    {
                        ht[j] = 0;
                        gt[j] = 0;
                    }

                    cross_corr(n2, ht, gt, dght);

                    for (j = 0; (j < nn); j++)
                    {
                        ct[j]  += rhbex[j];
                        ght[j] += dght[j];
                    }
                }
            }
        }
    }
    fprintf(stderr, "\n");
    sfree(h);
    sfree(g);
    normalizeACF(ct, ght, static_cast<int>(nhb), nn);

    /* Determine tail value for statistics */
    tail  = 0;
    tail2 = 0;
    for (j = nn/2; (j < nn); j++)
    {
        tail  += ct[j];
        tail2 += ct[j]*ct[j];
    }
    tail  /= (nn - nn/2);
    tail2 /= (nn - nn/2);
    dtail  = std::sqrt(tail2-tail*tail);

    /* Check whether the ACF is long enough */
    if (dtail > tol)
    {
        printf("\nWARNING: Correlation function is probably not long enough\n"
               "because the standard deviation in the tail of C(t) > %g\n"
               "Tail value (average C(t) over second half of acf): %g +/- %g\n",
               tol, tail, dtail);
    }
    for (j = 0; (j < nn); j++)
    {
        cct[j] = ct[j];
        ct[j]  = (cct[j]-tail)/(1-tail);
    }
    /* Compute negative derivative k(t) = -dc(t)/dt */
    compute_derivative(nn, hb->time, ct, kt);
    for (j = 0; (j < nn); j++)
    {
        kt[j] = -kt[j];
    }


    if (bContact)
    {
        fp = xvgropen(fn, "Contact Autocorrelation", output_env_get_xvgr_tlabel(oenv), "C(t)", oenv);
    }
    else
    {
        fp = xvgropen(fn, "Hydrogen Bond Autocorrelation", output_env_get_xvgr_tlabel(oenv), "C(t)", oenv);
    }
    xvgr_legend(fp, asize(legLuzar), legLuzar, oenv);


    for (j = 0; (j < nn); j++)
    {
        fprintf(fp, "%10g  %10g  %10g  %10g  %10g\n",
                hb->time[j]-hb->time[0], ct[j], cct[j], ght[j], kt[j]);
    }
    xvgrclose(fp);

    analyse_corr(nn, hb->time, ct, ght, kt, nullptr, nullptr, nullptr,
                 fit_start, temp);

    do_view(oenv, fn, nullptr);
    sfree(rhbex);
    sfree(ct);
    sfree(gt);
    sfree(ht);
    sfree(ght);
    sfree(dght);
    sfree(cct);
    sfree(kt);
}

static void init_hbframe(t_hbdata *hb, int nframes, real t)
{
    int i;

    hb->time[nframes]   = t;
    hb->nhb[nframes]    = 0;
    hb->ndist[nframes]  = 0;
    for (i = 0; (i < max_hx); i++)
    {
        hb->nhx[nframes][i] = 0;
    }
}

static FILE *open_donor_properties_file(const char             *fn,
                                        t_hbdata               *hb,
                                        const gmx_output_env_t *oenv)
{
    FILE       *fp    = nullptr;
    const char *leg[] = { "Nbound", "Nfree" };

    if (!fn || !hb)
    {
        return nullptr;
    }

    fp = xvgropen(fn, "Donor properties", output_env_get_xvgr_tlabel(oenv), "Number", oenv);
    xvgr_legend(fp, asize(leg), leg, oenv);

    return fp;
}

static void analyse_donor_properties(FILE *fp, t_hbdata *hb, int nframes, real t)
{
    int i, j, k, nbound, nb, nhtot;

    if (!fp || !hb)
    {
        return;
    }
    nbound = 0;
    nhtot  = 0;
    for (i = 0; (i < hb->d.nrd); i++)
    {
        for (k = 0; (k < hb->d.nhydro[i]); k++)
        {
            nb = 0;
            nhtot++;
            for (j = 0; (j < hb->a.nra) && (nb == 0); j++)
            {
                if (hb->hbmap[i][j] && hb->hbmap[i][j]->h[k] &&
                    is_hb(hb->hbmap[i][j]->h[k], nframes))
                {
                    nb = 1;
                }
            }
            nbound += nb;
        }
    }
    fprintf(fp, "%10.3e  %6d  %6d\n", t, nbound, nhtot-nbound);
}

static void dump_hbmap(t_hbdata *hb,
                       int nfile, t_filenm fnm[], gmx_bool bTwo,
                       gmx_bool bContact, int isize[], int *index[], char *grpnames[],
                       const t_atoms *atoms)
{
    FILE    *fp, *fplog;
    int      ddd, hhh, aaa, i, j, k, m, grp;
    char     ds[32], hs[32], as[32];
    gmx_bool first;

    fp = opt2FILE("-hbn", nfile, fnm, "w");
    if (opt2bSet("-g", nfile, fnm))
    {
        fplog = gmx_ffopen(opt2fn("-g", nfile, fnm), "w");
        fprintf(fplog, "# %10s  %12s  %12s\n", "Donor", "Hydrogen", "Acceptor");
    }
    else
    {
        fplog = nullptr;
    }
    for (grp = gr0; grp <= (bTwo ? gr1 : gr0); grp++)
    {
        fprintf(fp, "[ %s ]", grpnames[grp]);
        for (i = 0; i < isize[grp]; i++)
        {
            fprintf(fp, (i%15) ? " " : "\n");
            fprintf(fp, " %4d", index[grp][i]+1);
        }
        fprintf(fp, "\n");

        if (!bContact)
        {
            fprintf(fp, "[ donors_hydrogens_%s ]\n", grpnames[grp]);
            for (i = 0; (i < hb->d.nrd); i++)
            {
                if (hb->d.grp[i] == grp)
                {
                    for (j = 0; (j < hb->d.nhydro[i]); j++)
                    {
                        fprintf(fp, " %4d %4d", hb->d.don[i]+1,
                                hb->d.hydro[i][j]+1);
                    }
                    fprintf(fp, "\n");
                }
            }
            first = TRUE;
            fprintf(fp, "[ acceptors_%s ]", grpnames[grp]);
            for (i = 0; (i < hb->a.nra); i++)
            {
                if (hb->a.grp[i] == grp)
                {
                    fprintf(fp, (i%15 && !first) ? " " : "\n");
                    fprintf(fp, " %4d", hb->a.acc[i]+1);
                    first = FALSE;
                }
            }
            fprintf(fp, "\n");
        }
    }
    if (bTwo)
    {
        fprintf(fp, bContact ? "[ contacts_%s-%s ]\n" :
                "[ hbonds_%s-%s ]\n", grpnames[0], grpnames[1]);
    }
    else
    {
        fprintf(fp, bContact ? "[ contacts_%s ]" : "[ hbonds_%s ]\n", grpnames[0]);
    }

    for (i = 0; (i < hb->d.nrd); i++)
    {
        ddd = hb->d.don[i];
        for (k = 0; (k < hb->a.nra); k++)
        {
            aaa = hb->a.acc[k];
            for (m = 0; (m < hb->d.nhydro[i]); m++)
            {
                if (hb->hbmap[i][k] && ISHB(hb->hbmap[i][k]->history[m]))
                {
                    sprintf(ds, "%s", mkatomname(atoms, ddd));
                    sprintf(as, "%s", mkatomname(atoms, aaa));
                    if (bContact)
                    {
                        fprintf(fp, " %6d %6d\n", ddd+1, aaa+1);
                        if (fplog)
                        {
                            fprintf(fplog, "%12s  %12s\n", ds, as);
                        }
                    }
                    else
                    {
                        hhh = hb->d.hydro[i][m];
                        sprintf(hs, "%s", mkatomname(atoms, hhh));
                        fprintf(fp, " %6d %6d %6d\n", ddd+1, hhh+1, aaa+1);
                        if (fplog)
                        {
                            fprintf(fplog, "%12s  %12s  %12s\n", ds, hs, as);
                        }
                    }
                }
            }
        }
    }
    gmx_ffclose(fp);
    if (fplog)
    {
        gmx_ffclose(fplog);
    }
}

/* sync_hbdata() updates the parallel t_hbdata p_hb using hb as template.
 * It mimics add_frames() and init_frame() to some extent. */
static void sync_hbdata(t_hbdata *p_hb, int nframes)
{
    if (nframes >= p_hb->max_frames)
    {
        p_hb->max_frames += 4096;
        srenew(p_hb->nhb,   p_hb->max_frames);
        srenew(p_hb->ndist, p_hb->max_frames);
        srenew(p_hb->n_bound, p_hb->max_frames);
        srenew(p_hb->nhx, p_hb->max_frames);
        if (p_hb->bDAnr)
        {
            srenew(p_hb->danr, p_hb->max_frames);
        }
        std::memset(&(p_hb->nhb[nframes]),   0, sizeof(int) * (p_hb->max_frames-nframes));
        std::memset(&(p_hb->ndist[nframes]), 0, sizeof(int) * (p_hb->max_frames-nframes));
        p_hb->nhb[nframes]   = 0;
        p_hb->ndist[nframes] = 0;

    }
    p_hb->nframes = nframes;

    std::memset(&(p_hb->nhx[nframes]), 0, sizeof(int)*max_hx); /* zero the helix count for this frame */
}

int gmx_hbond(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] computes and analyzes hydrogen bonds. Hydrogen bonds are",
        "determined based on cutoffs for the angle Hydrogen - Donor - Acceptor",
        "(zero is extended) and the distance Donor - Acceptor",
        "(or Hydrogen - Acceptor using [TT]-noda[tt]).",
        "OH and NH groups are regarded as donors, O is an acceptor always,",
        "N is an acceptor by default, but this can be switched using",
        "[TT]-nitacc[tt]. Dummy hydrogen atoms are assumed to be connected",
        "to the first preceding non-hydrogen atom.[PAR]",

        "You need to specify two groups for analysis, which must be either",
        "identical or non-overlapping. All hydrogen bonds between the two",
        "groups are analyzed.[PAR]",

        "If you set [TT]-shell[tt], you will be asked for an additional index group",
        "which should contain exactly one atom. In this case, only hydrogen",
        "bonds between atoms within the shell distance from the one atom are",
        "considered.[PAR]",

        "With option -ac, rate constants for hydrogen bonding can be derived with the",
        "model of Luzar and Chandler (Nature 379:55, 1996; J. Chem. Phys. 113:23, 2000).",
        "If contact kinetics are analyzed by using the -contact option, then",
        "n(t) can be defined as either all pairs that are not within contact distance r at time t",
        "(corresponding to leaving the -r2 option at the default value 0) or all pairs that",
        "are within distance r2 (corresponding to setting a second cut-off value with option -r2).",
        "See mentioned literature for more details and definitions."
        "[PAR]",

        /*    "It is also possible to analyse specific hydrogen bonds with",
              "[TT]-sel[tt]. This index file must contain a group of atom triplets",
              "Donor Hydrogen Acceptor, in the following way::",
           "",
           "[ selected ]",
           "     20    21    24",
           "     25    26    29",
           "      1     3     6",
           "",
           "Note that the triplets need not be on separate lines.",
           "Each atom triplet specifies a hydrogen bond to be analyzed,",
           "note also that no check is made for the types of atoms.[PAR]",
         */

        "[BB]Output:[bb]",
        "",
        " * [TT]-num[tt]:  number of hydrogen bonds as a function of time.",
        " * [TT]-ac[tt]:   average over all autocorrelations of the existence",
        "   functions (either 0 or 1) of all hydrogen bonds.",
        " * [TT]-dist[tt]: distance distribution of all hydrogen bonds.",
        " * [TT]-ang[tt]:  angle distribution of all hydrogen bonds.",
        " * [TT]-hx[tt]:   the number of n-n+i hydrogen bonds as a function of time",
        "   where n and n+i stand for residue numbers and i ranges from 0 to 6.",
        "   This includes the n-n+3, n-n+4 and n-n+5 hydrogen bonds associated",
        "   with helices in proteins.",
        " * [TT]-hbn[tt]:  all selected groups, donors, hydrogens and acceptors",
        "   for selected groups, all hydrogen bonded atoms from all groups and",
        "   all solvent atoms involved in insertion.",
        " * [TT]-hbm[tt]:  existence matrix for all hydrogen bonds over all",
        "   frames, this also contains information on solvent insertion",
        "   into hydrogen bonds. Ordering is identical to that in [TT]-hbn[tt]",
        "   index file.",
        " * [TT]-dan[tt]: write out the number of donors and acceptors analyzed for",
        "   each timeframe. This is especially useful when using [TT]-shell[tt].",
        " * [TT]-nhbdist[tt]: compute the number of HBonds per hydrogen in order to",
        "   compare results to Raman Spectroscopy.",
        "",
        "Note: options [TT]-ac[tt], [TT]-life[tt], [TT]-hbn[tt] and [TT]-hbm[tt]",
        "require an amount of memory proportional to the total numbers of donors",
        "times the total number of acceptors in the selected group(s)."
    };

    static real        acut     = 30, abin = 1, rcut = 0.35, r2cut = 0, rbin = 0.005, rshell = -1;
    static real        maxnhb   = 0, fit_start = 1, fit_end = 60, temp = 298.15;
    static gmx_bool    bNitAcc  = TRUE, bDA = TRUE, bMerge = TRUE;
    static int         nDump    = 0;
    static int         nThreads = 0;

    static gmx_bool    bContact     = FALSE;

    /* options */
    t_pargs     pa [] = {
        { "-a",    FALSE,  etREAL, {&acut},
          "Cutoff angle (degrees, Hydrogen - Donor - Acceptor)" },
        { "-r",    FALSE,  etREAL, {&rcut},
          "Cutoff radius (nm, X - Acceptor, see next option)" },
        { "-da",   FALSE,  etBOOL, {&bDA},
          "Use distance Donor-Acceptor (if TRUE) or Hydrogen-Acceptor (FALSE)" },
        { "-r2",   FALSE,  etREAL, {&r2cut},
          "Second cutoff radius. Mainly useful with [TT]-contact[tt] and [TT]-ac[tt]"},
        { "-abin", FALSE,  etREAL, {&abin},
          "Binwidth angle distribution (degrees)" },
        { "-rbin", FALSE,  etREAL, {&rbin},
          "Binwidth distance distribution (nm)" },
        { "-nitacc", FALSE, etBOOL, {&bNitAcc},
          "Regard nitrogen atoms as acceptors" },
        { "-contact", FALSE, etBOOL, {&bContact},
          "Do not look for hydrogen bonds, but merely for contacts within the cut-off distance" },
        { "-shell", FALSE, etREAL, {&rshell},
          "when > 0, only calculate hydrogen bonds within # nm shell around "
          "one particle" },
        { "-fitstart", FALSE, etREAL, {&fit_start},
          "Time (ps) from which to start fitting the correlation functions in order to obtain the forward and backward rate constants for HB breaking and formation. With [TT]-gemfit[tt] we suggest [TT]-fitstart 0[tt]" },
        { "-fitend", FALSE, etREAL, {&fit_end},
          "Time (ps) to which to stop fitting the correlation functions in order to obtain the forward and backward rate constants for HB breaking and formation (only with [TT]-gemfit[tt])" },
        { "-temp",  FALSE, etREAL, {&temp},
          "Temperature (K) for computing the Gibbs energy corresponding to HB breaking and reforming" },
        { "-dump",  FALSE, etINT, {&nDump},
          "Dump the first N hydrogen bond ACFs in a single [REF].xvg[ref] file for debugging" },
        { "-max_hb", FALSE, etREAL, {&maxnhb},
          "Theoretical maximum number of hydrogen bonds used for normalizing HB autocorrelation function. Can be useful in case the program estimates it wrongly" },
        { "-merge", FALSE, etBOOL, {&bMerge},
          "H-bonds between the same donor and acceptor, but with different hydrogen are treated as a single H-bond. Mainly important for the ACF." },
#if GMX_OPENMP
        { "-nthreads", FALSE, etINT, {&nThreads},
          "Number of threads used for the parallel loop over autocorrelations. nThreads <= 0 means maximum number of threads. Requires linking with OpenMP. The number of threads is limited by the number of cores (before OpenMP v.3 ) or environment variable OMP_THREAD_LIMIT (OpenMP v.3)"},
#endif
    };
    const char *bugs[] = {
        "The option [TT]-sel[tt] that used to work on selected hbonds is out of order, and therefore not available for the time being."
    };
    t_filenm    fnm[] = {
        { efTRX, "-f",   nullptr,     ffREAD  },
        { efTPR, nullptr,   nullptr,     ffREAD  },
        { efNDX, nullptr,   nullptr,     ffOPTRD },
        /*    { efNDX, "-sel", "select", ffOPTRD },*/
        { efXVG, "-num", "hbnum",  ffWRITE },
        { efLOG, "-g",   "hbond",  ffOPTWR },
        { efXVG, "-ac",  "hbac",   ffOPTWR },
        { efXVG, "-dist", "hbdist", ffOPTWR },
        { efXVG, "-ang", "hbang",  ffOPTWR },
        { efXVG, "-hx",  "hbhelix", ffOPTWR },
        { efNDX, "-hbn", "hbond",  ffOPTWR },
        { efXPM, "-hbm", "hbmap",  ffOPTWR },
        { efXVG, "-don", "donor",  ffOPTWR },
        { efXVG, "-dan", "danum",  ffOPTWR },
        { efXVG, "-life", "hblife", ffOPTWR },
        { efXVG, "-nhbdist", "nhbdist", ffOPTWR }

    };
#define NFILE asize(fnm)

    char                  hbmap [HB_NR] = { ' ',    'o',      '-',       '*' };
    const char           *hbdesc[HB_NR] = { "None", "Present", "Inserted", "Present & Inserted" };
    t_rgb                 hbrgb [HB_NR] = { {1, 1, 1}, {1, 0, 0},   {0, 0, 1},    {1, 0, 1} };

    t_trxstatus          *status;
    int                   trrStatus = 1;
    t_topology            top;
    t_pargs              *ppa;
    int                   npargs, natoms, nframes = 0, shatom;
    int                  *isize;
    char                **grpnames;
    int                 **index;
    rvec                 *x, hbox;
    matrix                box;
    real                  t, ccut, dist = 0.0, ang = 0.0;
    double                max_nhb, aver_nhb, aver_dist;
    int                   h = 0, i = 0, j, k = 0, ogrp, nsel;
    int                   xi, yi, zi, ai;
    int                   xj, yj, zj, aj, xjj, yjj, zjj;
    gmx_bool              bSelected, bHBmap, bStop, bTwo, bBox, bTric;
    int                  *adist, *rdist;
    int                   grp, nabin, nrbin, resdist, ihb;
    char                **leg;
    t_hbdata             *hb;
    FILE                 *fp, *fpnhb = nullptr, *donor_properties = nullptr;
    t_gridcell         ***grid;
    t_ncell              *icell, *jcell;
    ivec                  ngrid;
    unsigned char        *datable;
    gmx_output_env_t     *oenv;
    int                   ii, hh, actual_nThreads;
    int                   threadNr = 0;
    gmx_bool              bParallel;
    gmx_bool              bEdge_yjj, bEdge_xjj;

    t_hbdata            **p_hb    = nullptr;                      /* one per thread, then merge after the frame loop */
    int                 **p_adist = nullptr, **p_rdist = nullptr; /* a histogram for each thread. */

    const bool            bOMP = GMX_OPENMP;

    npargs = asize(pa);
    ppa    = add_acf_pargs(&npargs, pa);

    if (!parse_common_args(&argc, argv, PCA_CAN_TIME | PCA_TIME_UNIT, NFILE, fnm, npargs,
                           ppa, asize(desc), desc, asize(bugs), bugs, &oenv))
    {
        sfree(ppa);
        return 0;
    }

    /* process input */
    bSelected = FALSE;
    ccut      = std::cos(acut*DEG2RAD);

    if (bContact)
    {
        if (bSelected)
        {
            gmx_fatal(FARGS, "Can not analyze selected contacts.");
        }
        if (!bDA)
        {
            gmx_fatal(FARGS, "Can not analyze contact between H and A: turn off -noda");
        }
    }

    /* Initiate main data structure! */
    bHBmap = (opt2bSet("-ac", NFILE, fnm) ||
              opt2bSet("-life", NFILE, fnm) ||
              opt2bSet("-hbn", NFILE, fnm) ||
              opt2bSet("-hbm", NFILE, fnm));

    if (opt2bSet("-nhbdist", NFILE, fnm))
    {
        const char *leg[MAXHH+1] = { "0 HBs", "1 HB", "2 HBs", "3 HBs", "Total" };
        fpnhb = xvgropen(opt2fn("-nhbdist", NFILE, fnm),
                         "Number of donor-H with N HBs", output_env_get_xvgr_tlabel(oenv), "N", oenv);
        xvgr_legend(fpnhb, asize(leg), leg, oenv);
    }

    hb = mk_hbdata(bHBmap, opt2bSet("-dan", NFILE, fnm), bMerge || bContact);

    /* get topology */
    t_inputrec      irInstance;
    t_inputrec     *ir = &irInstance;
    read_tpx_top(ftp2fn(efTPR, NFILE, fnm), ir, box, &natoms, nullptr, nullptr, &top);

    snew(grpnames, grNR);
    snew(index, grNR);
    snew(isize, grNR);
    /* Make Donor-Acceptor table */
    snew(datable, top.atoms.nr);

    if (bSelected)
    {
        /* analyze selected hydrogen bonds */
        printf("Select group with selected atoms:\n");
        get_index(&(top.atoms), opt2fn("-sel", NFILE, fnm),
                  1, &nsel, index, grpnames);
        if (nsel % 3)
        {
            gmx_fatal(FARGS, "Number of atoms in group '%s' not a multiple of 3\n"
                      "and therefore cannot contain triplets of "
                      "Donor-Hydrogen-Acceptor", grpnames[0]);
        }
        bTwo = FALSE;

        for (i = 0; (i < nsel); i += 3)
        {
            int dd = index[0][i];
            int aa = index[0][i+2];
            /* int */ hh = index[0][i+1];
            add_dh (&hb->d, dd, hh, i, datable);
            add_acc(&hb->a, aa, i);
            /* Should this be here ? */
            snew(hb->d.dptr, top.atoms.nr);
            snew(hb->a.aptr, top.atoms.nr);
            add_hbond(hb, dd, aa, hh, gr0, gr0, 0, bMerge, 0, bContact);
        }
        printf("Analyzing %d selected hydrogen bonds from '%s'\n",
               isize[0], grpnames[0]);
    }
    else
    {
        /* analyze all hydrogen bonds: get group(s) */
        printf("Specify 2 groups to analyze:\n");
        get_index(&(top.atoms), ftp2fn_null(efNDX, NFILE, fnm),
                  2, isize, index, grpnames);

        /* check if we have two identical or two non-overlapping groups */
        bTwo = isize[0] != isize[1];
        for (i = 0; (i < isize[0]) && !bTwo; i++)
        {
            bTwo = index[0][i] != index[1][i];
        }
        if (bTwo)
        {
            printf("Checking for overlap in atoms between %s and %s\n",
                   grpnames[0], grpnames[1]);

            gen_datable(index[0], isize[0], datable, top.atoms.nr);

            for (i = 0; i < isize[1]; i++)
            {
                if (ISINGRP(datable[index[1][i]]))
                {
                    gmx_fatal(FARGS, "Partial overlap between groups '%s' and '%s'",
                              grpnames[0], grpnames[1]);
                }
            }
        }
        if (bTwo)
        {
            printf("Calculating %s "
                   "between %s (%d atoms) and %s (%d atoms)\n",
                   bContact ? "contacts" : "hydrogen bonds",
                   grpnames[0], isize[0], grpnames[1], isize[1]);
        }
        else
        {
            fprintf(stderr, "Calculating %s in %s (%d atoms)\n",
                    bContact ? "contacts" : "hydrogen bonds", grpnames[0], isize[0]);
        }
    }
    sfree(datable);

    /* search donors and acceptors in groups */
    snew(datable, top.atoms.nr);
    for (i = 0; (i < grNR); i++)
    {
        if ( ((i == gr0) && !bSelected ) ||
             ((i == gr1) && bTwo ))
        {
            gen_datable(index[i], isize[i], datable, top.atoms.nr);
            if (bContact)
            {
                search_acceptors(&top, isize[i], index[i], &hb->a, i,
                                 bNitAcc, TRUE, (bTwo && (i == gr0)) || !bTwo, datable);
                search_donors   (&top, isize[i], index[i], &hb->d, i,
                                 TRUE, (bTwo && (i == gr1)) || !bTwo, datable);
            }
            else
            {
                search_acceptors(&top, isize[i], index[i], &hb->a, i, bNitAcc, FALSE, TRUE, datable);
                search_donors   (&top, isize[i], index[i], &hb->d, i, FALSE, TRUE, datable);
            }
            if (bTwo)
            {
                clear_datable_grp(datable, top.atoms.nr);
            }
        }
    }
    sfree(datable);
    printf("Found %d donors and %d acceptors\n", hb->d.nrd, hb->a.nra);
    /*if (bSelected)
       snew(donors[gr0D], dons[gr0D].nrd);*/

    donor_properties = open_donor_properties_file(opt2fn_null("-don", NFILE, fnm), hb, oenv);

    if (bHBmap)
    {
        printf("Making hbmap structure...");
        /* Generate hbond data structure */
        mk_hbmap(hb);
        printf("done.\n");
    }

    /* check input */
    bStop = FALSE;
    if (hb->d.nrd + hb->a.nra == 0)
    {
        printf("No Donors or Acceptors found\n");
        bStop = TRUE;
    }
    if (!bStop)
    {
        if (hb->d.nrd == 0)
        {
            printf("No Donors found\n");
            bStop = TRUE;
        }
        if (hb->a.nra == 0)
        {
            printf("No Acceptors found\n");
            bStop = TRUE;
        }
    }
    if (bStop)
    {
        gmx_fatal(FARGS, "Nothing to be done");
    }

    shatom = 0;
    if (rshell > 0)
    {
        int      shisz;
        int     *shidx;
        char    *shgrpnm;
        /* get index group with atom for shell */
        do
        {
            printf("Select atom for shell (1 atom):\n");
            get_index(&(top.atoms), ftp2fn_null(efNDX, NFILE, fnm),
                      1, &shisz, &shidx, &shgrpnm);
            if (shisz != 1)
            {
                printf("group contains %d atoms, should be 1 (one)\n", shisz);
            }
        }
        while (shisz != 1);
        shatom = shidx[0];
        printf("Will calculate hydrogen bonds within a shell "
               "of %g nm around atom %i\n", rshell, shatom+1);
    }

    /* Analyze trajectory */
    natoms = read_first_x(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &t, &x, box);
    if (natoms > top.atoms.nr)
    {
        gmx_fatal(FARGS, "Topology (%d atoms) does not match trajectory (%d atoms)",
                  top.atoms.nr, natoms);
    }

    bBox  = (ir->ePBC != epbcNONE);
    grid  = init_grid(bBox, box, (rcut > r2cut) ? rcut : r2cut, ngrid);
    nabin = static_cast<int>(acut/abin);
    nrbin = static_cast<int>(rcut/rbin);
    snew(adist, nabin+1);
    snew(rdist, nrbin+1);

#if !GMX_OPENMP
#define __ADIST adist
#define __RDIST rdist
#define __HBDATA hb
#else /* GMX_OPENMP ================================================== \
       * Set up the OpenMP stuff,                                       |
       * like the number of threads and such                            |
       * Also start the parallel loop.                                  |
       */
#define __ADIST p_adist[threadNr]
#define __RDIST p_rdist[threadNr]
#define __HBDATA p_hb[threadNr]
#endif
    if (bOMP)
    {
        bParallel = !bSelected;

        if (bParallel)
        {
            actual_nThreads = std::min((nThreads <= 0) ? INT_MAX : nThreads, gmx_omp_get_max_threads());

            gmx_omp_set_num_threads(actual_nThreads);
            printf("Frame loop parallelized with OpenMP using %i threads.\n", actual_nThreads);
            fflush(stdout);
        }
        else
        {
            actual_nThreads = 1;
        }

        snew(p_hb,    actual_nThreads);
        snew(p_adist, actual_nThreads);
        snew(p_rdist, actual_nThreads);
        for (i = 0; i < actual_nThreads; i++)
        {
            snew(p_hb[i], 1);
            snew(p_adist[i], nabin+1);
            snew(p_rdist[i], nrbin+1);

            p_hb[i]->max_frames = 0;
            p_hb[i]->nhb        = nullptr;
            p_hb[i]->ndist      = nullptr;
            p_hb[i]->n_bound    = nullptr;
            p_hb[i]->time       = nullptr;
            p_hb[i]->nhx        = nullptr;

            p_hb[i]->bHBmap     = hb->bHBmap;
            p_hb[i]->bDAnr      = hb->bDAnr;
            p_hb[i]->wordlen    = hb->wordlen;
            p_hb[i]->nframes    = hb->nframes;
            p_hb[i]->maxhydro   = hb->maxhydro;
            p_hb[i]->danr       = hb->danr;
            p_hb[i]->d          = hb->d;
            p_hb[i]->a          = hb->a;
            p_hb[i]->hbmap      = hb->hbmap;
            p_hb[i]->time       = hb->time; /* This may need re-syncing at every frame. */

            p_hb[i]->nrhb   = 0;
            p_hb[i]->nrdist = 0;
        }
    }

    /* Make a thread pool here,
     * instead of forking anew at every frame. */

#pragma omp parallel \
    firstprivate(i) \
    private(j, h, ii, hh, \
    xi, yi, zi, xj, yj, zj, threadNr, \
    dist, ang, icell, jcell, \
    grp, ogrp, ai, aj, xjj, yjj, zjj, \
    ihb, resdist, \
    k, bTric, \
    bEdge_xjj, bEdge_yjj) \
    default(shared)
    {                           /* Start of parallel region */
#if !defined __clang_analyzer__ // clang complains about unused value.
        threadNr = gmx_omp_get_thread_num();
#endif

        do
        {

            bTric = bBox && TRICLINIC(box);

            if (bOMP)
            {
                try
                {
                    sync_hbdata(p_hb[threadNr], nframes);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
#pragma omp single
            {
                try
                {
                    build_grid(hb, x, x[shatom], bBox, box, hbox, (rcut > r2cut) ? rcut : r2cut,
                               rshell, ngrid, grid);
                    reset_nhbonds(&(hb->d));

                    if (debug && bDebug)
                    {
                        dump_grid(debug, ngrid, grid);
                    }

                    add_frames(hb, nframes);
                    init_hbframe(hb, nframes, output_env_conv_time(oenv, t));

                    if (hb->bDAnr)
                    {
                        count_da_grid(ngrid, grid, hb->danr[nframes]);
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            } /* omp single */

            if (bOMP)
            {
                p_hb[threadNr]->time = hb->time; /* This pointer may have changed. */
            }

            if (bSelected)
            {

#pragma omp single
                {
                    try
                    {
                        /* Do not parallelize this just yet. */
                        /* int ii; */
                        for (ii = 0; (ii < nsel); ii++)
                        {
                            int dd = index[0][i];
                            int aa = index[0][i+2];
                            /* int */ hh = index[0][i+1];
                            ihb          = is_hbond(hb, ii, ii, dd, aa, rcut, r2cut, ccut, x, bBox, box,
                                                    hbox, &dist, &ang, bDA, &h, bContact, bMerge);

                            if (ihb)
                            {
                                /* add to index if not already there */
                                /* Add a hbond */
                                add_hbond(hb, dd, aa, hh, ii, ii, nframes, bMerge, ihb, bContact);
                            }
                        }
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                } /* omp single */
            }     /* if (bSelected) */
            else
            {
                /* The outer grid loop will have to do for now. */
#pragma omp for schedule(dynamic)
                for (xi = 0; xi < ngrid[XX]; xi++)
                {
                    try
                    {
                        for (yi = 0; (yi < ngrid[YY]); yi++)
                        {
                            for (zi = 0; (zi < ngrid[ZZ]); zi++)
                            {

                                /* loop over donor groups gr0 (always) and gr1 (if necessary) */
                                for (grp = gr0; (grp <= (bTwo ? gr1 : gr0)); grp++)
                                {
                                    icell = &(grid[zi][yi][xi].d[grp]);

                                    if (bTwo)
                                    {
                                        ogrp = 1-grp;
                                    }
                                    else
                                    {
                                        ogrp = grp;
                                    }

                                    /* loop over all hydrogen atoms from group (grp)
                                     * in this gridcell (icell)
                                     */
                                    for (ai = 0; (ai < icell->nr); ai++)
                                    {
                                        i  = icell->atoms[ai];

                                        /* loop over all adjacent gridcells (xj,yj,zj) */
                                        for (zjj = grid_loop_begin(ngrid[ZZ], zi, bTric, FALSE);
                                             zjj <= grid_loop_end(ngrid[ZZ], zi, bTric, FALSE);
                                             zjj++)
                                        {
                                            zj        = grid_mod(zjj, ngrid[ZZ]);
                                            bEdge_yjj = (zj == 0) || (zj == ngrid[ZZ] - 1);
                                            for (yjj = grid_loop_begin(ngrid[YY], yi, bTric, bEdge_yjj);
                                                 yjj <= grid_loop_end(ngrid[YY], yi, bTric, bEdge_yjj);
                                                 yjj++)
                                            {
                                                yj        = grid_mod(yjj, ngrid[YY]);
                                                bEdge_xjj =
                                                    (yj == 0) || (yj == ngrid[YY] - 1) ||
                                                    (zj == 0) || (zj == ngrid[ZZ] - 1);
                                                for (xjj = grid_loop_begin(ngrid[XX], xi, bTric, bEdge_xjj);
                                                     xjj <= grid_loop_end(ngrid[XX], xi, bTric, bEdge_xjj);
                                                     xjj++)
                                                {
                                                    xj    = grid_mod(xjj, ngrid[XX]);
                                                    jcell = &(grid[zj][yj][xj].a[ogrp]);
                                                    /* loop over acceptor atoms from other group (ogrp)
                                                     * in this adjacent gridcell (jcell)
                                                     */
                                                    for (aj = 0; (aj < jcell->nr); aj++)
                                                    {
                                                        j = jcell->atoms[aj];

                                                        /* check if this once was a h-bond */
                                                        ihb  = is_hbond(__HBDATA, grp, ogrp, i, j, rcut, r2cut, ccut, x, bBox, box,
                                                                        hbox, &dist, &ang, bDA, &h, bContact, bMerge);

                                                        if (ihb)
                                                        {
                                                            /* add to index if not already there */
                                                            /* Add a hbond */
                                                            add_hbond(__HBDATA, i, j, h, grp, ogrp, nframes, bMerge, ihb, bContact);

                                                            /* make angle and distance distributions */
                                                            if (ihb == hbHB && !bContact)
                                                            {
                                                                if (dist > rcut)
                                                                {
                                                                    gmx_fatal(FARGS, "distance is higher than what is allowed for an hbond: %f", dist);
                                                                }
                                                                ang *= RAD2DEG;
                                                                __ADIST[static_cast<int>( ang/abin)]++;
                                                                __RDIST[static_cast<int>(dist/rbin)]++;
                                                                if (!bTwo)
                                                                {
                                                                    if (donor_index(&hb->d, grp, i) == NOTSET)
                                                                    {
                                                                        gmx_fatal(FARGS, "Invalid donor %d", i);
                                                                    }
                                                                    if (acceptor_index(&hb->a, ogrp, j) == NOTSET)
                                                                    {
                                                                        gmx_fatal(FARGS, "Invalid acceptor %d", j);
                                                                    }
                                                                    resdist = std::abs(top.atoms.atom[i].resind-top.atoms.atom[j].resind);
                                                                    if (resdist >= max_hx)
                                                                    {
                                                                        resdist = max_hx-1;
                                                                    }
                                                                    __HBDATA->nhx[nframes][resdist]++;
                                                                }
                                                            }

                                                        }
                                                    } /* for aj  */
                                                }     /* for xjj */
                                            }         /* for yjj */
                                        }             /* for zjj */
                                    }                 /* for ai  */
                                }                     /* for grp */
                            }                         /* for xi,yi,zi */
                        }
                    }
                    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
                }
            } /* if (bSelected) {...} else */


            /* Better wait for all threads to finnish using x[] before updating it. */
            k = nframes;
#pragma omp barrier
#pragma omp critical
            {
                try
                {
                    /* Sum up histograms and counts from p_hb[] into hb */
                    if (bOMP)
                    {
                        hb->nhb[k]   += p_hb[threadNr]->nhb[k];
                        hb->ndist[k] += p_hb[threadNr]->ndist[k];
                        for (j = 0; j < max_hx; j++)
                        {
                            hb->nhx[k][j]  += p_hb[threadNr]->nhx[k][j];
                        }
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

            /* Here are a handful of single constructs
             * to share the workload a bit. The most
             * important one is of course the last one,
             * where there's a potential bottleneck in form
             * of slow I/O.                    */
#pragma omp barrier
#pragma omp single
            {
                try
                {
                    analyse_donor_properties(donor_properties, hb, k, t);
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

#pragma omp single
            {
                try
                {
                    if (fpnhb)
                    {
                        do_nhb_dist(fpnhb, hb, t);
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

#pragma omp single
            {
                try
                {
                    trrStatus = (read_next_x(oenv, status, &t, x, box));
                    nframes++;
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

#pragma omp barrier
        }
        while (trrStatus);

        if (bOMP)
        {
#pragma omp critical
            {
                hb->nrhb   += p_hb[threadNr]->nrhb;
                hb->nrdist += p_hb[threadNr]->nrdist;
            }

            /* Free parallel datastructures */
            sfree(p_hb[threadNr]->nhb);
            sfree(p_hb[threadNr]->ndist);
            sfree(p_hb[threadNr]->nhx);

#pragma omp for
            for (i = 0; i < nabin; i++)
            {
                try
                {
                    for (j = 0; j < actual_nThreads; j++)
                    {

                        adist[i] += p_adist[j][i];
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }
#pragma omp for
            for (i = 0; i <= nrbin; i++)
            {
                try
                {
                    for (j = 0; j < actual_nThreads; j++)
                    {
                        rdist[i] += p_rdist[j][i];
                    }
                }
                GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
            }

            sfree(p_adist[threadNr]);
            sfree(p_rdist[threadNr]);
        }
    } /* End of parallel region */
    if (bOMP)
    {
        sfree(p_adist);
        sfree(p_rdist);
    }

    if (nframes < 2 && (opt2bSet("-ac", NFILE, fnm) || opt2bSet("-life", NFILE, fnm)))
    {
        gmx_fatal(FARGS, "Cannot calculate autocorrelation of life times with less than two frames");
    }

    free_grid(ngrid, &grid);

    close_trx(status);

    if (donor_properties)
    {
        xvgrclose(donor_properties);
    }

    if (fpnhb)
    {
        xvgrclose(fpnhb);
    }

    /* Compute maximum possible number of different hbonds */
    if (maxnhb > 0)
    {
        max_nhb = maxnhb;
    }
    else
    {
        max_nhb = 0.5*(hb->d.nrd*hb->a.nra);
    }

    if (bHBmap)
    {
        if (hb->nrhb == 0)
        {
            printf("No %s found!!\n", bContact ? "contacts" : "hydrogen bonds");
        }
        else
        {
            printf("Found %d different %s in trajectory\n"
                   "Found %d different atom-pairs within %s distance\n",
                   hb->nrhb, bContact ? "contacts" : "hydrogen bonds",
                   hb->nrdist, (r2cut > 0) ? "second cut-off" : "hydrogen bonding");

            if (bMerge)
            {
                merge_hb(hb, bTwo, bContact);
            }

            if (opt2bSet("-hbn", NFILE, fnm))
            {
                dump_hbmap(hb, NFILE, fnm, bTwo, bContact, isize, index, grpnames, &top.atoms);
            }

            /* Moved the call to merge_hb() to a line BEFORE dump_hbmap
             * to make the -hbn and -hmb output match eachother.
             * - Erik Marklund, May 30, 2006 */
        }
    }
    /* Print out number of hbonds and distances */
    aver_nhb  = 0;
    aver_dist = 0;
    fp        = xvgropen(opt2fn("-num", NFILE, fnm), bContact ? "Contacts" :
                         "Hydrogen Bonds", output_env_get_xvgr_tlabel(oenv), "Number", oenv);
    snew(leg, 2);
    snew(leg[0], STRLEN);
    snew(leg[1], STRLEN);
    sprintf(leg[0], "%s", bContact ? "Contacts" : "Hydrogen bonds");
    sprintf(leg[1], "Pairs within %g nm", (r2cut > 0) ? r2cut : rcut);
    xvgr_legend(fp, 2, (const char**)leg, oenv);
    sfree(leg[1]);
    sfree(leg[0]);
    sfree(leg);
    for (i = 0; (i < nframes); i++)
    {
        fprintf(fp, "%10g  %10d  %10d\n", hb->time[i], hb->nhb[i], hb->ndist[i]);
        aver_nhb  += hb->nhb[i];
        aver_dist += hb->ndist[i];
    }
    xvgrclose(fp);
    aver_nhb  /= nframes;
    aver_dist /= nframes;
    /* Print HB distance distribution */
    if (opt2bSet("-dist", NFILE, fnm))
    {
        int sum;

        sum = 0;
        for (i = 0; i < nrbin; i++)
        {
            sum += rdist[i];
        }

        fp = xvgropen(opt2fn("-dist", NFILE, fnm),
                      "Hydrogen Bond Distribution",
                      bDA ?
                      "Donor - Acceptor Distance (nm)" :
                      "Hydrogen - Acceptor Distance (nm)", "", oenv);
        for (i = 0; i < nrbin; i++)
        {
            fprintf(fp, "%10g %10g\n", (i+0.5)*rbin, rdist[i]/(rbin*sum));
        }
        xvgrclose(fp);
    }

    /* Print HB angle distribution */
    if (opt2bSet("-ang", NFILE, fnm))
    {
        long sum;

        sum = 0;
        for (i = 0; i < nabin; i++)
        {
            sum += adist[i];
        }

        fp = xvgropen(opt2fn("-ang", NFILE, fnm),
                      "Hydrogen Bond Distribution",
                      "Hydrogen - Donor - Acceptor Angle (\\SO\\N)", "", oenv);
        for (i = 0; i < nabin; i++)
        {
            fprintf(fp, "%10g %10g\n", (i+0.5)*abin, adist[i]/(abin*sum));
        }
        xvgrclose(fp);
    }

    /* Print HB in alpha-helix */
    if (opt2bSet("-hx", NFILE, fnm))
    {
        fp = xvgropen(opt2fn("-hx", NFILE, fnm),
                      "Hydrogen Bonds", output_env_get_xvgr_tlabel(oenv), "Count", oenv);
        xvgr_legend(fp, NRHXTYPES, hxtypenames, oenv);
        for (i = 0; i < nframes; i++)
        {
            fprintf(fp, "%10g", hb->time[i]);
            for (j = 0; j < max_hx; j++)
            {
                fprintf(fp, " %6d", hb->nhx[i][j]);
            }
            fprintf(fp, "\n");
        }
        xvgrclose(fp);
    }

    printf("Average number of %s per timeframe %.3f out of %g possible\n",
           bContact ? "contacts" : "hbonds",
           bContact ? aver_dist : aver_nhb, max_nhb);

    /* Do Autocorrelation etc. */
    if (hb->bHBmap)
    {
        /*
           Added support for -contact in ac and hbm calculations below.
           - Erik Marklund, May 29, 2006
         */
        if (opt2bSet("-ac", NFILE, fnm) || opt2bSet("-life", NFILE, fnm))
        {
            please_cite(stdout, "Spoel2006b");
        }
        if (opt2bSet("-ac", NFILE, fnm))
        {
            do_hbac(opt2fn("-ac", NFILE, fnm), hb, nDump,
                    bMerge, bContact, fit_start, temp, r2cut > 0, oenv,
                    nThreads);
        }
        if (opt2bSet("-life", NFILE, fnm))
        {
            do_hblife(opt2fn("-life", NFILE, fnm), hb, bMerge, bContact, oenv);
        }
        if (opt2bSet("-hbm", NFILE, fnm))
        {
            t_matrix mat;
            int      id, ia, hh, x, y;
            mat.flags = mat.y0 = 0;

            if ((nframes > 0) && (hb->nrhb > 0))
            {
                mat.nx = nframes;
                mat.ny = hb->nrhb;

                snew(mat.matrix, mat.nx);
                for (x = 0; (x < mat.nx); x++)
                {
                    snew(mat.matrix[x], mat.ny);
                }
                y = 0;
                for (id = 0; (id < hb->d.nrd); id++)
                {
                    for (ia = 0; (ia < hb->a.nra); ia++)
                    {
                        for (hh = 0; (hh < hb->maxhydro); hh++)
                        {
                            if (hb->hbmap[id][ia])
                            {
                                if (ISHB(hb->hbmap[id][ia]->history[hh]))
                                {
                                    for (x = 0; (x <= hb->hbmap[id][ia]->nframes); x++)
                                    {
                                        int nn0 = hb->hbmap[id][ia]->n0;
                                        range_check(y, 0, mat.ny);
                                        mat.matrix[x+nn0][y] = is_hb(hb->hbmap[id][ia]->h[hh], x);
                                    }
                                    y++;
                                }
                            }
                        }
                    }
                }
                mat.axis_x = hb->time;
                snew(mat.axis_y, mat.ny);
                for (j = 0; j < mat.ny; j++)
                {
                    mat.axis_y[j] = j;
                }
                sprintf(mat.title, bContact ? "Contact Existence Map" :
                        "Hydrogen Bond Existence Map");
                sprintf(mat.legend, bContact ? "Contacts" : "Hydrogen Bonds");
                sprintf(mat.label_x, "%s", output_env_get_xvgr_tlabel(oenv).c_str());
                sprintf(mat.label_y, bContact ? "Contact Index" : "Hydrogen Bond Index");
                mat.bDiscrete = TRUE;
                mat.nmap      = 2;
                snew(mat.map, mat.nmap);
                for (i = 0; i < mat.nmap; i++)
                {
                    mat.map[i].code.c1 = hbmap[i];
                    mat.map[i].desc    = hbdesc[i];
                    mat.map[i].rgb     = hbrgb[i];
                }
                fp = opt2FILE("-hbm", NFILE, fnm, "w");
                write_xpm_m(fp, mat);
                gmx_ffclose(fp);
                for (x = 0; x < mat.nx; x++)
                {
                    sfree(mat.matrix[x]);
                }
                sfree(mat.axis_y);
                sfree(mat.matrix);
                sfree(mat.map);
            }
            else
            {
                fprintf(stderr, "No hydrogen bonds/contacts found. No hydrogen bond map will be printed.\n");
            }
        }
    }

    if (hb->bDAnr)
    {
        int    i, j, nleg;
        char **legnames;
        char   buf[STRLEN];

#define USE_THIS_GROUP(j) ( (j == gr0) || (bTwo && (j == gr1)) )

        fp = xvgropen(opt2fn("-dan", NFILE, fnm),
                      "Donors and Acceptors", output_env_get_xvgr_tlabel(oenv), "Count", oenv);
        nleg = (bTwo ? 2 : 1)*2;
        snew(legnames, nleg);
        i = 0;
        for (j = 0; j < grNR; j++)
        {
            if (USE_THIS_GROUP(j) )
            {
                sprintf(buf, "Donors %s", grpnames[j]);
                legnames[i++] = gmx_strdup(buf);
                sprintf(buf, "Acceptors %s", grpnames[j]);
                legnames[i++] = gmx_strdup(buf);
            }
        }
        if (i != nleg)
        {
            gmx_incons("number of legend entries");
        }
        xvgr_legend(fp, nleg, (const char**)legnames, oenv);
        for (i = 0; i < nframes; i++)
        {
            fprintf(fp, "%10g", hb->time[i]);
            for (j = 0; (j < grNR); j++)
            {
                if (USE_THIS_GROUP(j) )
                {
                    fprintf(fp, " %6d", hb->danr[i][j]);
                }
            }
            fprintf(fp, "\n");
        }
        xvgrclose(fp);
    }

    return 0;
}
