/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
#include <stdlib.h>
#include <string.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/cmat.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/legacyheaders/viewit.h"
#include "gromacs/linearalgebra/eigensolver.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/random/random.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

/* print to two file pointers at once (i.e. stderr and log) */
static gmx_inline
void lo_ffprintf(FILE *fp1, FILE *fp2, const char *buf)
{
    fprintf(fp1, "%s", buf);
    fprintf(fp2, "%s", buf);
}

/* just print a prepared buffer to fp1 and fp2 */
static gmx_inline
void ffprintf(FILE *fp1, FILE *fp2, const char *buf)
{
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with one argument, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_d(FILE *fp1, FILE *fp2, char *buf, const char *fmt, int arg)
{
    sprintf(buf, fmt, arg);
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with one argument, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_g(FILE *fp1, FILE *fp2, char *buf, const char *fmt, real arg)
{
    sprintf(buf, fmt, arg);
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with one argument, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_s(FILE *fp1, FILE *fp2, char *buf, const char *fmt, const char *arg)
{
    sprintf(buf, fmt, arg);
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with two arguments, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_dd(FILE *fp1, FILE *fp2, char *buf, const char *fmt, int arg1, int arg2)
{
    sprintf(buf, fmt, arg1, arg2);
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with two arguments, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_gg(FILE *fp1, FILE *fp2, char *buf, const char *fmt, real arg1, real arg2)
{
    sprintf(buf, fmt, arg1, arg2);
    lo_ffprintf(fp1, fp2, buf);
}

/* prepare buffer with two arguments, then print to fp1 and fp2 */
static gmx_inline
void ffprintf_ss(FILE *fp1, FILE *fp2, char *buf, const char *fmt, const char *arg1, const char *arg2)
{
    sprintf(buf, fmt, arg1, arg2);
    lo_ffprintf(fp1, fp2, buf);
}

typedef struct {
    int  ncl;
    int *cl;
} t_clusters;

typedef struct {
    int  nr;
    int *nb;
} t_nnb;

void cp_index(int nn, int from[], int to[])
{
    int i;

    for (i = 0; (i < nn); i++)
    {
        to[i] = from[i];
    }
}

void mc_optimize(FILE *log, t_mat *m, real *time,
                 int maxiter, int nrandom,
                 int seed, real kT,
                 const char *conv, output_env_t oenv)
{
    FILE      *fp = NULL;
    real       ecur, enext, emin, prob, enorm;
    int        i, j, iswap, jswap, nn, nuphill = 0;
    gmx_rng_t  rng;
    t_mat     *minimum;

    if (m->n1 != m->nn)
    {
        fprintf(stderr, "Can not do Monte Carlo optimization with a non-square matrix.\n");
        return;
    }
    printf("\nDoing Monte Carlo optimization to find the smoothest trajectory\n");
    printf("by reordering the frames to minimize the path between the two structures\n");
    printf("that have the largest pairwise RMSD.\n");

    iswap = jswap = -1;
    enorm = m->mat[0][0];
    for (i = 0; (i < m->n1); i++)
    {
        for (j = 0; (j < m->nn); j++)
        {
            if (m->mat[i][j] > enorm)
            {
                enorm   = m->mat[i][j];
                iswap   = i;
                jswap   = j;
            }
        }
    }
    if ((iswap == -1) || (jswap == -1))
    {
        fprintf(stderr, "Matrix contains identical values in all fields\n");
        return;
    }
    swap_rows(m, 0, iswap);
    swap_rows(m, m->n1-1, jswap);
    emin = ecur = mat_energy(m);
    printf("Largest distance %g between %d and %d. Energy: %g.\n",
           enorm, iswap, jswap, emin);

    rng = gmx_rng_init(seed);
    nn  = m->nn;

    /* Initiate and store global minimum */
    minimum     = init_mat(nn, m->b1D);
    minimum->nn = nn;
    copy_t_mat(minimum, m);

    if (NULL != conv)
    {
        fp = xvgropen(conv, "Convergence of the MC optimization",
                      "Energy", "Step", oenv);
    }
    for (i = 0; (i < maxiter); i++)
    {
        /* Generate new swapping candidates */
        do
        {
            iswap = 1+(nn-2)*gmx_rng_uniform_real(rng);
            jswap = 1+(nn-2)*gmx_rng_uniform_real(rng);
        }
        while ((iswap == jswap) || (iswap >= nn-1) || (jswap >= nn-1));

        /* Apply swap and compute energy */
        swap_rows(m, iswap, jswap);
        enext = mat_energy(m);

        /* Compute probability */
        prob = 0;
        if ((enext < ecur) || (i < nrandom))
        {
            prob = 1;
            if (enext < emin)
            {
                /* Store global minimum */
                copy_t_mat(minimum, m);
                emin = enext;
            }
        }
        else if (kT > 0)
        {
            /* Try Monte Carlo step */
            prob = exp(-(enext-ecur)/(enorm*kT));
        }

        if ((prob == 1) || (gmx_rng_uniform_real(rng) < prob))
        {
            if (enext > ecur)
            {
                nuphill++;
            }

            fprintf(log, "Iter: %d Swapped %4d and %4d (energy: %g prob: %g)\n",
                    i, iswap, jswap, enext, prob);
            if (NULL != fp)
            {
                fprintf(fp, "%6d  %10g\n", i, enext);
            }
            ecur = enext;
        }
        else
        {
            swap_rows(m, jswap, iswap);
        }
    }
    fprintf(log, "%d uphill steps were taken during optimization\n",
            nuphill);

    /* Now swap the matrix to get it into global minimum mode */
    copy_t_mat(m, minimum);

    fprintf(log, "Global minimum energy %g\n", mat_energy(minimum));
    fprintf(log, "Global minimum energy %g\n", mat_energy(m));
    fprintf(log, "Swapped time and frame indices and RMSD to next neighbor:\n");
    for (i = 0; (i < m->nn); i++)
    {
        fprintf(log, "%10g  %5d  %10g\n",
                time[m->m_ind[i]],
                m->m_ind[i],
                (i < m->nn-1) ? m->mat[m->m_ind[i]][m->m_ind[i+1]] : 0);
    }

    if (NULL != fp)
    {
        xvgrclose(fp);
    }
}

static void calc_dist(int nind, rvec x[], real **d)
{
    int      i, j;
    real    *xi;
    rvec     dx;

    for (i = 0; (i < nind-1); i++)
    {
        xi = x[i];
        for (j = i+1; (j < nind); j++)
        {
            /* Should use pbc_dx when analysing multiple molecueles,
             * but the box is not stored for every frame.
             */
            rvec_sub(xi, x[j], dx);
            d[i][j] = norm(dx);
        }
    }
}

static real rms_dist(int isize, real **d, real **d_r)
{
    int  i, j;
    real r, r2;

    r2 = 0.0;
    for (i = 0; (i < isize-1); i++)
    {
        for (j = i+1; (j < isize); j++)
        {
            r   = d[i][j]-d_r[i][j];
            r2 += r*r;
        }
    }
    r2 /= (isize*(isize-1))/2;

    return sqrt(r2);
}

static int rms_dist_comp(const void *a, const void *b)
{
    t_dist *da, *db;

    da = (t_dist *)a;
    db = (t_dist *)b;

    if (da->dist - db->dist < 0)
    {
        return -1;
    }
    else if (da->dist - db->dist > 0)
    {
        return 1;
    }
    return 0;
}

static int clust_id_comp(const void *a, const void *b)
{
    t_clustid *da, *db;

    da = (t_clustid *)a;
    db = (t_clustid *)b;

    return da->clust - db->clust;
}

static int nrnb_comp(const void *a, const void *b)
{
    t_nnb *da, *db;

    da = (t_nnb *)a;
    db = (t_nnb *)b;

    /* return the b-a, we want highest first */
    return db->nr - da->nr;
}

void gather(t_mat *m, real cutoff, t_clusters *clust)
{
    t_clustid *c;
    t_dist    *d;
    int        i, j, k, nn, cid, n1, diff;
    gmx_bool   bChange;

    /* First we sort the entries in the RMSD matrix */
    n1 = m->nn;
    nn = ((n1-1)*n1)/2;
    snew(d, nn);
    for (i = k = 0; (i < n1); i++)
    {
        for (j = i+1; (j < n1); j++, k++)
        {
            d[k].i    = i;
            d[k].j    = j;
            d[k].dist = m->mat[i][j];
        }
    }
    if (k != nn)
    {
        gmx_incons("gather algortihm");
    }
    qsort(d, nn, sizeof(d[0]), rms_dist_comp);

    /* Now we make a cluster index for all of the conformations */
    c = new_clustid(n1);

    /* Now we check the closest structures, and equalize their cluster numbers */
    fprintf(stderr, "Linking structures ");
    do
    {
        fprintf(stderr, "*");
        bChange = FALSE;
        for (k = 0; (k < nn) && (d[k].dist < cutoff); k++)
        {
            diff = c[d[k].j].clust - c[d[k].i].clust;
            if (diff)
            {
                bChange = TRUE;
                if (diff > 0)
                {
                    c[d[k].j].clust = c[d[k].i].clust;
                }
                else
                {
                    c[d[k].i].clust = c[d[k].j].clust;
                }
            }
        }
    }
    while (bChange);
    fprintf(stderr, "\nSorting and renumbering clusters\n");
    /* Sort on cluster number */
    qsort(c, n1, sizeof(c[0]), clust_id_comp);

    /* Renumber clusters */
    cid = 1;
    for (k = 1; k < n1; k++)
    {
        if (c[k].clust != c[k-1].clust)
        {
            c[k-1].clust = cid;
            cid++;
        }
        else
        {
            c[k-1].clust = cid;
        }
    }
    c[k-1].clust = cid;
    if (debug)
    {
        for (k = 0; (k < n1); k++)
        {
            fprintf(debug, "Cluster index for conformation %d: %d\n",
                    c[k].conf, c[k].clust);
        }
    }
    clust->ncl = cid;
    for (k = 0; k < n1; k++)
    {
        clust->cl[c[k].conf] = c[k].clust;
    }

    sfree(c);
    sfree(d);
}

gmx_bool jp_same(int **nnb, int i, int j, int P)
{
    gmx_bool bIn;
    int      k, ii, jj, pp;

    bIn = FALSE;
    for (k = 0; nnb[i][k] >= 0; k++)
    {
        bIn = bIn || (nnb[i][k] == j);
    }
    if (!bIn)
    {
        return FALSE;
    }

    bIn = FALSE;
    for (k = 0; nnb[j][k] >= 0; k++)
    {
        bIn = bIn || (nnb[j][k] == i);
    }
    if (!bIn)
    {
        return FALSE;
    }

    pp = 0;
    for (ii = 0; nnb[i][ii] >= 0; ii++)
    {
        for (jj = 0; nnb[j][jj] >= 0; jj++)
        {
            if ((nnb[i][ii] == nnb[j][jj]) && (nnb[i][ii] != -1))
            {
                pp++;
            }
        }
    }

    return (pp >= P);
}

static void jarvis_patrick(int n1, real **mat, int M, int P,
                           real rmsdcut, t_clusters *clust)
{
    t_dist     *row;
    t_clustid  *c;
    int       **nnb;
    int         i, j, k, cid, diff, max;
    gmx_bool    bChange;
    real      **mcpy = NULL;

    if (rmsdcut < 0)
    {
        rmsdcut = 10000;
    }

    /* First we sort the entries in the RMSD matrix row by row.
     * This gives us the nearest neighbor list.
     */
    snew(nnb, n1);
    snew(row, n1);
    for (i = 0; (i < n1); i++)
    {
        for (j = 0; (j < n1); j++)
        {
            row[j].j    = j;
            row[j].dist = mat[i][j];
        }
        qsort(row, n1, sizeof(row[0]), rms_dist_comp);
        if (M > 0)
        {
            /* Put the M nearest neighbors in the list */
            snew(nnb[i], M+1);
            for (j = k = 0; (k < M) && (j < n1) && (mat[i][row[j].j] < rmsdcut); j++)
            {
                if (row[j].j  != i)
                {
                    nnb[i][k]  = row[j].j;
                    k++;
                }
            }
            nnb[i][k] = -1;
        }
        else
        {
            /* Put all neighbors nearer than rmsdcut in the list */
            max = 0;
            k   = 0;
            for (j = 0; (j < n1) && (mat[i][row[j].j] < rmsdcut); j++)
            {
                if (row[j].j != i)
                {
                    if (k >= max)
                    {
                        max += 10;
                        srenew(nnb[i], max);
                    }
                    nnb[i][k] = row[j].j;
                    k++;
                }
            }
            if (k == max)
            {
                srenew(nnb[i], max+1);
            }
            nnb[i][k] = -1;
        }
    }
    sfree(row);
    if (debug)
    {
        fprintf(debug, "Nearest neighborlist. M = %d, P = %d\n", M, P);
        for (i = 0; (i < n1); i++)
        {
            fprintf(debug, "i:%5d nbs:", i);
            for (j = 0; nnb[i][j] >= 0; j++)
            {
                fprintf(debug, "%5d[%5.3f]", nnb[i][j], mat[i][nnb[i][j]]);
            }
            fprintf(debug, "\n");
        }
    }

    c = new_clustid(n1);
    fprintf(stderr, "Linking structures ");
    /* Use mcpy for temporary storage of booleans */
    mcpy = mk_matrix(n1, n1, FALSE);
    for (i = 0; i < n1; i++)
    {
        for (j = i+1; j < n1; j++)
        {
            mcpy[i][j] = jp_same(nnb, i, j, P);
        }
    }
    do
    {
        fprintf(stderr, "*");
        bChange = FALSE;
        for (i = 0; i < n1; i++)
        {
            for (j = i+1; j < n1; j++)
            {
                if (mcpy[i][j])
                {
                    diff = c[j].clust - c[i].clust;
                    if (diff)
                    {
                        bChange = TRUE;
                        if (diff > 0)
                        {
                            c[j].clust = c[i].clust;
                        }
                        else
                        {
                            c[i].clust = c[j].clust;
                        }
                    }
                }
            }
        }
    }
    while (bChange);

    fprintf(stderr, "\nSorting and renumbering clusters\n");
    /* Sort on cluster number */
    qsort(c, n1, sizeof(c[0]), clust_id_comp);

    /* Renumber clusters */
    cid = 1;
    for (k = 1; k < n1; k++)
    {
        if (c[k].clust != c[k-1].clust)
        {
            c[k-1].clust = cid;
            cid++;
        }
        else
        {
            c[k-1].clust = cid;
        }
    }
    c[k-1].clust = cid;
    clust->ncl   = cid;
    for (k = 0; k < n1; k++)
    {
        clust->cl[c[k].conf] = c[k].clust;
    }
    if (debug)
    {
        for (k = 0; (k < n1); k++)
        {
            fprintf(debug, "Cluster index for conformation %d: %d\n",
                    c[k].conf, c[k].clust);
        }
    }

/* Again, I don't see the point in this... (AF) */
/*   for(i=0; (i<n1); i++) { */
/*     for(j=0; (j<n1); j++) */
/*       mcpy[c[i].conf][c[j].conf] = mat[i][j]; */
/*   } */
/*   for(i=0; (i<n1); i++) { */
/*     for(j=0; (j<n1); j++) */
/*       mat[i][j] = mcpy[i][j]; */
/*   } */
    done_matrix(n1, &mcpy);

    sfree(c);
    for (i = 0; (i < n1); i++)
    {
        sfree(nnb[i]);
    }
    sfree(nnb);
}

static void dump_nnb (FILE *fp, const char *title, int n1, t_nnb *nnb)
{
    int i, j;

    /* dump neighbor list */
    fprintf(fp, "%s", title);
    for (i = 0; (i < n1); i++)
    {
        fprintf(fp, "i:%5d #:%5d nbs:", i, nnb[i].nr);
        for (j = 0; j < nnb[i].nr; j++)
        {
            fprintf(fp, "%5d", nnb[i].nb[j]);
        }
        fprintf(fp, "\n");
    }
}

static void gromos(int n1, real **mat, real rmsdcut, t_clusters *clust)
{
    t_dist *row;
    t_nnb  *nnb;
    int     i, j, k, j1, max;

    /* Put all neighbors nearer than rmsdcut in the list */
    fprintf(stderr, "Making list of neighbors within cutoff ");
    snew(nnb, n1);
    snew(row, n1);
    for (i = 0; (i < n1); i++)
    {
        max = 0;
        k   = 0;
        /* put all neighbors within cut-off in list */
        for (j = 0; j < n1; j++)
        {
            if (mat[i][j] < rmsdcut)
            {
                if (k >= max)
                {
                    max += 10;
                    srenew(nnb[i].nb, max);
                }
                nnb[i].nb[k] = j;
                k++;
            }
        }
        /* store nr of neighbors, we'll need that */
        nnb[i].nr = k;
        if (i%(1+n1/100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (i*100+1)/n1);
        }
    }
    fprintf(stderr, "%3d%%\n", 100);
    sfree(row);

    /* sort neighbor list on number of neighbors, largest first */
    qsort(nnb, n1, sizeof(nnb[0]), nrnb_comp);

    if (debug)
    {
        dump_nnb(debug, "Nearest neighborlist after sort.\n", n1, nnb);
    }

    /* turn first structure with all its neighbors (largest) into cluster
       remove them from pool of structures and repeat for all remaining */
    fprintf(stderr, "Finding clusters %4d", 0);
    /* cluster id's start at 1: */
    k = 1;
    while (nnb[0].nr)
    {
        /* set cluster id (k) for first item in neighborlist */
        for (j = 0; j < nnb[0].nr; j++)
        {
            clust->cl[nnb[0].nb[j]] = k;
        }
        /* mark as done */
        nnb[0].nr = 0;
        sfree(nnb[0].nb);

        /* adjust number of neighbors for others, taking removals into account: */
        for (i = 1; i < n1 && nnb[i].nr; i++)
        {
            j1 = 0;
            for (j = 0; j < nnb[i].nr; j++)
            {
                /* if this neighbor wasn't removed */
                if (clust->cl[nnb[i].nb[j]] == 0)
                {
                    /* shift the rest (j1<=j) */
                    nnb[i].nb[j1] = nnb[i].nb[j];
                    /* next */
                    j1++;
                }
            }
            /* now j1 is the new number of neighbors */
            nnb[i].nr = j1;
        }
        /* sort again on nnb[].nr, because we have new # neighbors: */
        /* but we only need to sort upto i, i.e. when nnb[].nr>0 */
        qsort(nnb, i, sizeof(nnb[0]), nrnb_comp);

        fprintf(stderr, "\b\b\b\b%4d", k);
        /* new cluster id */
        k++;
    }
    fprintf(stderr, "\n");
    sfree(nnb);
    if (debug)
    {
        fprintf(debug, "Clusters (%d):\n", k);
        for (i = 0; i < n1; i++)
        {
            fprintf(debug, " %3d", clust->cl[i]);
        }
        fprintf(debug, "\n");
    }

    clust->ncl = k-1;
}

rvec **read_whole_trj(const char *fn, int isize, atom_id index[], int skip,
                      int *nframe, real **time, const output_env_t oenv, gmx_bool bPBC, gmx_rmpbc_t gpbc)
{
    rvec       **xx, *x;
    matrix       box;
    real         t;
    int          i, i0, j, max_nf;
    int          natom;
    t_trxstatus *status;


    max_nf = 0;
    xx     = NULL;
    *time  = NULL;
    natom  = read_first_x(oenv, &status, fn, &t, &x, box);
    i      = 0;
    i0     = 0;
    do
    {
        if (bPBC)
        {
            gmx_rmpbc(gpbc, natom, box, x);
        }
        if (i0 >= max_nf)
        {
            max_nf += 10;
            srenew(xx, max_nf);
            srenew(*time, max_nf);
        }
        if ((i % skip) == 0)
        {
            snew(xx[i0], isize);
            /* Store only the interesting atoms */
            for (j = 0; (j < isize); j++)
            {
                copy_rvec(x[index[j]], xx[i0][j]);
            }
            (*time)[i0] = t;
            i0++;
        }
        i++;
    }
    while (read_next_x(oenv, status, &t, x, box));
    fprintf(stderr, "Allocated %lu bytes for frames\n",
            (unsigned long) (max_nf*isize*sizeof(**xx)));
    fprintf(stderr, "Read %d frames from trajectory %s\n", i0, fn);
    *nframe = i0;
    sfree(x);

    return xx;
}

static int plot_clusters(int nf, real **mat, t_clusters *clust,
                         int minstruct)
{
    int  i, j, ncluster, ci;
    int *cl_id, *nstruct, *strind;

    snew(cl_id, nf);
    snew(nstruct, nf);
    snew(strind, nf);
    for (i = 0; i < nf; i++)
    {
        strind[i] = 0;
        cl_id[i]  = clust->cl[i];
        nstruct[cl_id[i]]++;
    }
    ncluster = 0;
    for (i = 0; i < nf; i++)
    {
        if (nstruct[i] >= minstruct)
        {
            ncluster++;
            for (j = 0; (j < nf); j++)
            {
                if (cl_id[j] == i)
                {
                    strind[j] = ncluster;
                }
            }
        }
    }
    ncluster++;
    fprintf(stderr, "There are %d clusters with at least %d conformations\n",
            ncluster, minstruct);

    for (i = 0; (i < nf); i++)
    {
        ci = cl_id[i];
        for (j = 0; j < i; j++)
        {
            if ((ci == cl_id[j]) && (nstruct[ci] >= minstruct))
            {
                /* color different clusters with different colors, as long as
                   we don't run out of colors */
                mat[i][j] = strind[i];
            }
            else
            {
                mat[i][j] = 0;
            }
        }
    }
    sfree(strind);
    sfree(nstruct);
    sfree(cl_id);

    return ncluster;
}

static void mark_clusters(int nf, real **mat, real val, t_clusters *clust)
{
    int i, j, v;

    for (i = 0; i < nf; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (clust->cl[i] == clust->cl[j])
            {
                mat[i][j] = val;
            }
            else
            {
                mat[i][j] = 0;
            }
        }
    }
}

static char *parse_filename(const char *fn, int maxnr)
{
    int         i;
    char       *fnout;
    const char *ext;
    char        buf[STRLEN];

    if (strchr(fn, '%'))
    {
        gmx_fatal(FARGS, "will not number filename %s containing '%c'", fn, '%');
    }
    /* number of digits needed in numbering */
    i = (int)(log(maxnr)/log(10)) + 1;
    /* split fn and ext */
    ext = strrchr(fn, '.');
    if (!ext)
    {
        gmx_fatal(FARGS, "cannot separate extension in filename %s", fn);
    }
    ext++;
    /* insert e.g. '%03d' between fn and ext */
    sprintf(buf, "%s%%0%dd.%s", fn, i, ext);
    snew(fnout, strlen(buf)+1);
    strcpy(fnout, buf);

    return fnout;
}

static void ana_trans(t_clusters *clust, int nf,
                      const char *transfn, const char *ntransfn, FILE *log,
                      t_rgb rlo, t_rgb rhi, const output_env_t oenv)
{
    FILE  *fp;
    real **trans, *axis;
    int   *ntrans;
    int    i, ntranst, maxtrans;
    char   buf[STRLEN];

    snew(ntrans, clust->ncl);
    snew(trans, clust->ncl);
    snew(axis, clust->ncl);
    for (i = 0; i < clust->ncl; i++)
    {
        axis[i] = i+1;
        snew(trans[i], clust->ncl);
    }
    ntranst  = 0;
    maxtrans = 0;
    for (i = 1; i < nf; i++)
    {
        if (clust->cl[i] != clust->cl[i-1])
        {
            ntranst++;
            ntrans[clust->cl[i-1]-1]++;
            ntrans[clust->cl[i]-1]++;
            trans[clust->cl[i-1]-1][clust->cl[i]-1]++;
            maxtrans = max(maxtrans, trans[clust->cl[i]-1][clust->cl[i-1]-1]);
        }
    }
    ffprintf_dd(stderr, log, buf, "Counted %d transitions in total, "
                "max %d between two specific clusters\n", ntranst, maxtrans);
    if (transfn)
    {
        fp = gmx_ffopen(transfn, "w");
        i  = min(maxtrans+1, 80);
        write_xpm(fp, 0, "Cluster Transitions", "# transitions",
                  "from cluster", "to cluster",
                  clust->ncl, clust->ncl, axis, axis, trans,
                  0, maxtrans, rlo, rhi, &i);
        gmx_ffclose(fp);
    }
    if (ntransfn)
    {
        fp = xvgropen(ntransfn, "Cluster Transitions", "Cluster #", "# transitions",
                      oenv);
        for (i = 0; i < clust->ncl; i++)
        {
            fprintf(fp, "%5d %5d\n", i+1, ntrans[i]);
        }
        xvgrclose(fp);
    }
    sfree(ntrans);
    for (i = 0; i < clust->ncl; i++)
    {
        sfree(trans[i]);
    }
    sfree(trans);
    sfree(axis);
}

static void analyze_clusters(int nf, t_clusters *clust, real **rmsd,
                             int natom, t_atoms *atoms, rvec *xtps,
                             real *mass, rvec **xx, real *time,
                             int ifsize, atom_id *fitidx,
                             int iosize, atom_id *outidx,
                             const char *trxfn, const char *sizefn,
                             const char *transfn, const char *ntransfn,
                             const char *clustidfn, gmx_bool bAverage,
                             int write_ncl, int write_nst, real rmsmin,
                             gmx_bool bFit, FILE *log, t_rgb rlo, t_rgb rhi,
                             const output_env_t oenv)
{
    FILE        *size_fp = NULL;
    char         buf[STRLEN], buf1[40], buf2[40], buf3[40], *trxsfn;
    t_trxstatus *trxout  = NULL;
    t_trxstatus *trxsout = NULL;
    int          i, i1, cl, nstr, *structure, first = 0, midstr;
    gmx_bool    *bWrite = NULL;
    real         r, clrmsd, midrmsd;
    rvec        *xav = NULL;
    matrix       zerobox;

    clear_mat(zerobox);

    ffprintf_d(stderr, log, buf, "\nFound %d clusters\n\n", clust->ncl);
    trxsfn = NULL;
    if (trxfn)
    {
        /* do we write all structures? */
        if (write_ncl)
        {
            trxsfn = parse_filename(trxfn, max(write_ncl, clust->ncl));
            snew(bWrite, nf);
        }
        ffprintf_ss(stderr, log, buf, "Writing %s structure for each cluster to %s\n",
                    bAverage ? "average" : "middle", trxfn);
        if (write_ncl)
        {
            /* find out what we want to tell the user:
               Writing [all structures|structures with rmsd > %g] for
               {all|first %d} clusters {with more than %d structures} to %s     */
            if (rmsmin > 0.0)
            {
                sprintf(buf1, "structures with rmsd > %g", rmsmin);
            }
            else
            {
                sprintf(buf1, "all structures");
            }
            buf2[0] = buf3[0] = '\0';
            if (write_ncl >= clust->ncl)
            {
                if (write_nst == 0)
                {
                    sprintf(buf2, "all ");
                }
            }
            else
            {
                sprintf(buf2, "the first %d ", write_ncl);
            }
            if (write_nst)
            {
                sprintf(buf3, " with more than %d structures", write_nst);
            }
            sprintf(buf, "Writing %s for %sclusters%s to %s\n", buf1, buf2, buf3, trxsfn);
            ffprintf(stderr, log, buf);
        }

        /* Prepare a reference structure for the orientation of the clusters  */
        if (bFit)
        {
            reset_x(ifsize, fitidx, natom, NULL, xtps, mass);
        }
        trxout = open_trx(trxfn, "w");
        /* Calculate the average structure in each cluster,               *
         * all structures are fitted to the first struture of the cluster */
        snew(xav, natom);
    }

    if (transfn || ntransfn)
    {
        ana_trans(clust, nf, transfn, ntransfn, log, rlo, rhi, oenv);
    }

    if (clustidfn)
    {
        FILE *fp = xvgropen(clustidfn, "Clusters", output_env_get_xvgr_tlabel(oenv), "Cluster #", oenv);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(fp, "@    s0 symbol 2\n");
            fprintf(fp, "@    s0 symbol size 0.2\n");
            fprintf(fp, "@    s0 linestyle 0\n");
        }
        for (i = 0; i < nf; i++)
        {
            fprintf(fp, "%8g %8d\n", time[i], clust->cl[i]);
        }
        xvgrclose(fp);
    }
    if (sizefn)
    {
        size_fp = xvgropen(sizefn, "Cluster Sizes", "Cluster #", "# Structures", oenv);
        if (output_env_get_print_xvgr_codes(oenv))
        {
            fprintf(size_fp, "@g%d type %s\n", 0, "bar");
        }
    }
    snew(structure, nf);
    fprintf(log, "\n%3s | %3s  %4s | %6s %4s | cluster members\n",
            "cl.", "#st", "rmsd", "middle", "rmsd");
    for (cl = 1; cl <= clust->ncl; cl++)
    {
        /* prepare structures (fit, middle, average) */
        if (xav)
        {
            for (i = 0; i < natom; i++)
            {
                clear_rvec(xav[i]);
            }
        }
        nstr = 0;
        for (i1 = 0; i1 < nf; i1++)
        {
            if (clust->cl[i1] == cl)
            {
                structure[nstr] = i1;
                nstr++;
                if (trxfn && (bAverage || write_ncl) )
                {
                    if (bFit)
                    {
                        reset_x(ifsize, fitidx, natom, NULL, xx[i1], mass);
                    }
                    if (nstr == 1)
                    {
                        first = i1;
                    }
                    else if (bFit)
                    {
                        do_fit(natom, mass, xx[first], xx[i1]);
                    }
                    if (xav)
                    {
                        for (i = 0; i < natom; i++)
                        {
                            rvec_inc(xav[i], xx[i1][i]);
                        }
                    }
                }
            }
        }
        if (sizefn)
        {
            fprintf(size_fp, "%8d %8d\n", cl, nstr);
        }
        clrmsd  = 0;
        midstr  = 0;
        midrmsd = 10000;
        for (i1 = 0; i1 < nstr; i1++)
        {
            r = 0;
            if (nstr > 1)
            {
                for (i = 0; i < nstr; i++)
                {
                    if (i < i1)
                    {
                        r += rmsd[structure[i]][structure[i1]];
                    }
                    else
                    {
                        r += rmsd[structure[i1]][structure[i]];
                    }
                }
                r /= (nstr - 1);
            }
            if (r < midrmsd)
            {
                midstr  = structure[i1];
                midrmsd = r;
            }
            clrmsd += r;
        }
        clrmsd /= nstr;

        /* dump cluster info to logfile */
        if (nstr > 1)
        {
            sprintf(buf1, "%6.3f", clrmsd);
            if (buf1[0] == '0')
            {
                buf1[0] = ' ';
            }
            sprintf(buf2, "%5.3f", midrmsd);
            if (buf2[0] == '0')
            {
                buf2[0] = ' ';
            }
        }
        else
        {
            sprintf(buf1, "%5s", "");
            sprintf(buf2, "%5s", "");
        }
        fprintf(log, "%3d | %3d %s | %6g%s |", cl, nstr, buf1, time[midstr], buf2);
        for (i = 0; i < nstr; i++)
        {
            if ((i % 7 == 0) && i)
            {
                sprintf(buf, "\n%3s | %3s  %4s | %6s %4s |", "", "", "", "", "");
            }
            else
            {
                buf[0] = '\0';
            }
            i1 = structure[i];
            fprintf(log, "%s %6g", buf, time[i1]);
        }
        fprintf(log, "\n");

        /* write structures to trajectory file(s) */
        if (trxfn)
        {
            if (write_ncl)
            {
                for (i = 0; i < nstr; i++)
                {
                    bWrite[i] = FALSE;
                }
            }
            if (cl < write_ncl+1 && nstr > write_nst)
            {
                /* Dump all structures for this cluster */
                /* generate numbered filename (there is a %d in trxfn!) */
                sprintf(buf, trxsfn, cl);
                trxsout = open_trx(buf, "w");
                for (i = 0; i < nstr; i++)
                {
                    bWrite[i] = TRUE;
                    if (rmsmin > 0.0)
                    {
                        for (i1 = 0; i1 < i && bWrite[i]; i1++)
                        {
                            if (bWrite[i1])
                            {
                                bWrite[i] = rmsd[structure[i1]][structure[i]] > rmsmin;
                            }
                        }
                    }
                    if (bWrite[i])
                    {
                        write_trx(trxsout, iosize, outidx, atoms, i, time[structure[i]], zerobox,
                                  xx[structure[i]], NULL, NULL);
                    }
                }
                close_trx(trxsout);
            }
            /* Dump the average structure for this cluster */
            if (bAverage)
            {
                for (i = 0; i < natom; i++)
                {
                    svmul(1.0/nstr, xav[i], xav[i]);
                }
            }
            else
            {
                for (i = 0; i < natom; i++)
                {
                    copy_rvec(xx[midstr][i], xav[i]);
                }
                if (bFit)
                {
                    reset_x(ifsize, fitidx, natom, NULL, xav, mass);
                }
            }
            if (bFit)
            {
                do_fit(natom, mass, xtps, xav);
            }
            r = cl;
            write_trx(trxout, iosize, outidx, atoms, cl, time[midstr], zerobox, xav, NULL, NULL);
        }
    }
    /* clean up */
    if (trxfn)
    {
        close_trx(trxout);
        sfree(xav);
        if (write_ncl)
        {
            sfree(bWrite);
        }
    }
    sfree(structure);
    if (trxsfn)
    {
        sfree(trxsfn);
    }

    if (size_fp)
    {
        xvgrclose(size_fp);
    }
}

static void convert_mat(t_matrix *mat, t_mat *rms)
{
    int i, j;

    rms->n1 = mat->nx;
    matrix2real(mat, rms->mat);
    /* free input xpm matrix data */
    for (i = 0; i < mat->nx; i++)
    {
        sfree(mat->matrix[i]);
    }
    sfree(mat->matrix);

    for (i = 0; i < mat->nx; i++)
    {
        for (j = i; j < mat->nx; j++)
        {
            rms->sumrms += rms->mat[i][j];
            rms->maxrms  = max(rms->maxrms, rms->mat[i][j]);
            if (j != i)
            {
                rms->minrms = min(rms->minrms, rms->mat[i][j]);
            }
        }
    }
    rms->nn = mat->nx;
}

int gmx_cluster(int argc, char *argv[])
{
    const char        *desc[] = {
        "[THISMODULE] can cluster structures using several different methods.",
        "Distances between structures can be determined from a trajectory",
        "or read from an [REF].xpm[ref] matrix file with the [TT]-dm[tt] option.",
        "RMS deviation after fitting or RMS deviation of atom-pair distances",
        "can be used to define the distance between structures.[PAR]",

        "single linkage: add a structure to a cluster when its distance to any",
        "element of the cluster is less than [TT]cutoff[tt].[PAR]",

        "Jarvis Patrick: add a structure to a cluster when this structure",
        "and a structure in the cluster have each other as neighbors and",
        "they have a least [TT]P[tt] neighbors in common. The neighbors",
        "of a structure are the M closest structures or all structures within",
        "[TT]cutoff[tt].[PAR]",

        "Monte Carlo: reorder the RMSD matrix using Monte Carlo such that",
        "the order of the frames is using the smallest possible increments.",
        "With this it is possible to make a smooth animation going from one",
        "structure to another with the largest possible (e.g.) RMSD between",
        "them, however the intermediate steps should be as small as possible.",
        "Applications could be to visualize a potential of mean force",
        "ensemble of simulations or a pulling simulation. Obviously the user",
        "has to prepare the trajectory well (e.g. by not superimposing frames).",
        "The final result can be inspect visually by looking at the matrix",
        "[REF].xpm[ref] file, which should vary smoothly from bottom to top.[PAR]",

        "diagonalization: diagonalize the RMSD matrix.[PAR]",

        "gromos: use algorithm as described in Daura [IT]et al.[it]",
        "([IT]Angew. Chem. Int. Ed.[it] [BB]1999[bb], [IT]38[it], pp 236-240).",
        "Count number of neighbors using cut-off, take structure with",
        "largest number of neighbors with all its neighbors as cluster",
        "and eliminate it from the pool of clusters. Repeat for remaining",
        "structures in pool.[PAR]",

        "When the clustering algorithm assigns each structure to exactly one",
        "cluster (single linkage, Jarvis Patrick and gromos) and a trajectory",
        "file is supplied, the structure with",
        "the smallest average distance to the others or the average structure",
        "or all structures for each cluster will be written to a trajectory",
        "file. When writing all structures, separate numbered files are made",
        "for each cluster.[PAR]",

        "Two output files are always written:",
        "",
        " * [TT]-o[tt] writes the RMSD values in the upper left half of the matrix",
        "   and a graphical depiction of the clusters in the lower right half",
        "   When [TT]-minstruct[tt] = 1 the graphical depiction is black",
        "   when two structures are in the same cluster.",
        "   When [TT]-minstruct[tt] > 1 different colors will be used for each",
        "   cluster.",
        " * [TT]-g[tt] writes information on the options used and a detailed list",
        "   of all clusters and their members.",
        "",

        "Additionally, a number of optional output files can be written:",
        "",
        " * [TT]-dist[tt] writes the RMSD distribution.",
        " * [TT]-ev[tt] writes the eigenvectors of the RMSD matrix",
        "   diagonalization.",
        " * [TT]-sz[tt] writes the cluster sizes.",
        " * [TT]-tr[tt] writes a matrix of the number transitions between",
        "   cluster pairs.",
        " * [TT]-ntr[tt] writes the total number of transitions to or from",
        "   each cluster.",
        " * [TT]-clid[tt] writes the cluster number as a function of time.",
        " * [TT]-cl[tt] writes average (with option [TT]-av[tt]) or central",
        "   structure of each cluster or writes numbered files with cluster members",
        "   for a selected set of clusters (with option [TT]-wcl[tt], depends on",
        "   [TT]-nst[tt] and [TT]-rmsmin[tt]). The center of a cluster is the",
        "   structure with the smallest average RMSD from all other structures",
        "   of the cluster.",
    };

    FILE              *fp, *log;
    int                nf, i, i1, i2, j;
    gmx_int64_t        nrms = 0;

    matrix             box;
    rvec              *xtps, *usextps, *x1, **xx = NULL;
    const char        *fn, *trx_out_fn;
    t_clusters         clust;
    t_mat             *rms, *orig = NULL;
    real              *eigenvalues;
    t_topology         top;
    int                ePBC;
    t_atoms            useatoms;
    t_matrix          *readmat = NULL;
    real              *eigenvectors;

    int                isize = 0, ifsize = 0, iosize = 0;
    atom_id           *index = NULL, *fitidx, *outidx;
    char              *grpname;
    real               rmsd, **d1, **d2, *time = NULL, time_invfac, *mass = NULL;
    char               buf[STRLEN], buf1[80], title[STRLEN];
    gmx_bool           bAnalyze, bUseRmsdCut, bJP_RMSD = FALSE, bReadMat, bReadTraj, bPBC = TRUE;

    int                method, ncluster = 0;
    static const char *methodname[] = {
        NULL, "linkage", "jarvis-patrick", "monte-carlo",
        "diagonalization", "gromos", NULL
    };
    enum {
        m_null, m_linkage, m_jarvis_patrick,
        m_monte_carlo, m_diagonalize, m_gromos, m_nr
    };
    /* Set colors for plotting: white = zero RMS, black = maximum */
    static t_rgb rlo_top  = { 1.0, 1.0, 1.0 };
    static t_rgb rhi_top  = { 0.0, 0.0, 0.0 };
    static t_rgb rlo_bot  = { 1.0, 1.0, 1.0 };
    static t_rgb rhi_bot  = { 0.0, 0.0, 1.0 };
    static int   nlevels  = 40, skip = 1;
    static real  scalemax = -1.0, rmsdcut = 0.1, rmsmin = 0.0;
    gmx_bool     bRMSdist = FALSE, bBinary = FALSE, bAverage = FALSE, bFit = TRUE;
    static int   niter    = 10000, nrandom = 0, seed = 1993, write_ncl = 0, write_nst = 1, minstruct = 1;
    static real  kT       = 1e-3;
    static int   M        = 10, P = 3;
    output_env_t oenv;
    gmx_rmpbc_t  gpbc = NULL;

    t_pargs      pa[] = {
        { "-dista", FALSE, etBOOL, {&bRMSdist},
          "Use RMSD of distances instead of RMS deviation" },
        { "-nlevels", FALSE, etINT,  {&nlevels},
          "Discretize RMSD matrix in this number of levels" },
        { "-cutoff", FALSE, etREAL, {&rmsdcut},
          "RMSD cut-off (nm) for two structures to be neighbor" },
        { "-fit",   FALSE, etBOOL, {&bFit},
          "Use least squares fitting before RMSD calculation" },
        { "-max",   FALSE, etREAL, {&scalemax},
          "Maximum level in RMSD matrix" },
        { "-skip",  FALSE, etINT,  {&skip},
          "Only analyze every nr-th frame" },
        { "-av",    FALSE, etBOOL, {&bAverage},
          "Write average iso middle structure for each cluster" },
        { "-wcl",   FALSE, etINT,  {&write_ncl},
          "Write the structures for this number of clusters to numbered files" },
        { "-nst",   FALSE, etINT,  {&write_nst},
          "Only write all structures if more than this number of structures per cluster" },
        { "-rmsmin", FALSE, etREAL, {&rmsmin},
          "minimum rms difference with rest of cluster for writing structures" },
        { "-method", FALSE, etENUM, {methodname},
          "Method for cluster determination" },
        { "-minstruct", FALSE, etINT, {&minstruct},
          "Minimum number of structures in cluster for coloring in the [REF].xpm[ref] file" },
        { "-binary", FALSE, etBOOL, {&bBinary},
          "Treat the RMSD matrix as consisting of 0 and 1, where the cut-off "
          "is given by [TT]-cutoff[tt]" },
        { "-M",     FALSE, etINT,  {&M},
          "Number of nearest neighbors considered for Jarvis-Patrick algorithm, "
          "0 is use cutoff" },
        { "-P",     FALSE, etINT,  {&P},
          "Number of identical nearest neighbors required to form a cluster" },
        { "-seed",  FALSE, etINT,  {&seed},
          "Random number seed for Monte Carlo clustering algorithm: <= 0 means generate" },
        { "-niter", FALSE, etINT,  {&niter},
          "Number of iterations for MC" },
        { "-nrandom", FALSE, etINT,  {&nrandom},
          "The first iterations for MC may be done complete random, to shuffle the frames" },
        { "-kT",    FALSE, etREAL, {&kT},
          "Boltzmann weighting factor for Monte Carlo optimization "
          "(zero turns off uphill steps)" },
        { "-pbc", FALSE, etBOOL,
          { &bPBC }, "PBC check" }
    };
    t_filenm     fnm[] = {
        { efTRX, "-f",     NULL,        ffOPTRD },
        { efTPS, "-s",     NULL,        ffOPTRD },
        { efNDX, NULL,     NULL,        ffOPTRD },
        { efXPM, "-dm",   "rmsd",       ffOPTRD },
        { efXPM, "-om",   "rmsd-raw",   ffWRITE },
        { efXPM, "-o",    "rmsd-clust", ffWRITE },
        { efLOG, "-g",    "cluster",    ffWRITE },
        { efXVG, "-dist", "rmsd-dist",  ffOPTWR },
        { efXVG, "-ev",   "rmsd-eig",   ffOPTWR },
        { efXVG, "-conv", "mc-conv",    ffOPTWR },
        { efXVG, "-sz",   "clust-size", ffOPTWR},
        { efXPM, "-tr",   "clust-trans", ffOPTWR},
        { efXVG, "-ntr",  "clust-trans", ffOPTWR},
        { efXVG, "-clid", "clust-id",   ffOPTWR},
        { efTRX, "-cl",   "clusters.pdb", ffOPTWR }
    };
#define NFILE asize(fnm)

    if (!parse_common_args(&argc, argv,
                           PCA_CAN_VIEW | PCA_CAN_TIME | PCA_TIME_UNIT,
                           NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL,
                           &oenv))
    {
        return 0;
    }

    /* parse options */
    bReadMat   = opt2bSet("-dm", NFILE, fnm);
    bReadTraj  = opt2bSet("-f", NFILE, fnm) || !bReadMat;
    if (opt2parg_bSet("-av", asize(pa), pa) ||
        opt2parg_bSet("-wcl", asize(pa), pa) ||
        opt2parg_bSet("-nst", asize(pa), pa) ||
        opt2parg_bSet("-rmsmin", asize(pa), pa) ||
        opt2bSet("-cl", NFILE, fnm) )
    {
        trx_out_fn = opt2fn("-cl", NFILE, fnm);
    }
    else
    {
        trx_out_fn = NULL;
    }
    if (bReadMat && output_env_get_time_factor(oenv) != 1)
    {
        fprintf(stderr,
                "\nWarning: assuming the time unit in %s is %s\n",
                opt2fn("-dm", NFILE, fnm), output_env_get_time_unit(oenv));
    }
    if (trx_out_fn && !bReadTraj)
    {
        fprintf(stderr, "\nWarning: "
                "cannot write cluster structures without reading trajectory\n"
                "         ignoring option -cl %s\n", trx_out_fn);
    }

    method = 1;
    while (method < m_nr && gmx_strcasecmp(methodname[0], methodname[method]) != 0)
    {
        method++;
    }
    if (method == m_nr)
    {
        gmx_fatal(FARGS, "Invalid method");
    }

    bAnalyze = (method == m_linkage || method == m_jarvis_patrick ||
                method == m_gromos );

    /* Open log file */
    log = ftp2FILE(efLOG, NFILE, fnm, "w");

    fprintf(stderr, "Using %s method for clustering\n", methodname[0]);
    fprintf(log, "Using %s method for clustering\n", methodname[0]);

    /* check input and write parameters to log file */
    bUseRmsdCut = FALSE;
    if (method == m_jarvis_patrick)
    {
        bJP_RMSD = (M == 0) || opt2parg_bSet("-cutoff", asize(pa), pa);
        if ((M < 0) || (M == 1))
        {
            gmx_fatal(FARGS, "M (%d) must be 0 or larger than 1", M);
        }
        if (M < 2)
        {
            sprintf(buf1, "Will use P=%d and RMSD cutoff (%g)", P, rmsdcut);
            bUseRmsdCut = TRUE;
        }
        else
        {
            if (P >= M)
            {
                gmx_fatal(FARGS, "Number of neighbors required (P) must be less than M");
            }
            if (bJP_RMSD)
            {
                sprintf(buf1, "Will use P=%d, M=%d and RMSD cutoff (%g)", P, M, rmsdcut);
                bUseRmsdCut = TRUE;
            }
            else
            {
                sprintf(buf1, "Will use P=%d, M=%d", P, M);
            }
        }
        ffprintf_s(stderr, log, buf, "%s for determining the neighbors\n\n", buf1);
    }
    else /* method != m_jarvis */
    {
        bUseRmsdCut = ( bBinary || method == m_linkage || method == m_gromos );
    }
    if (bUseRmsdCut && method != m_jarvis_patrick)
    {
        fprintf(log, "Using RMSD cutoff %g nm\n", rmsdcut);
    }
    if (method == m_monte_carlo)
    {
        fprintf(log, "Using %d iterations\n", niter);
    }

    if (skip < 1)
    {
        gmx_fatal(FARGS, "skip (%d) should be >= 1", skip);
    }

    /* get input */
    if (bReadTraj)
    {
        /* don't read mass-database as masses (and top) are not used */
        read_tps_conf(ftp2fn(efTPS, NFILE, fnm), buf, &top, &ePBC, &xtps, NULL, box,
                      TRUE);
        if (bPBC)
        {
            gpbc = gmx_rmpbc_init(&top.idef, ePBC, top.atoms.nr);
        }

        fprintf(stderr, "\nSelect group for least squares fit%s:\n",
                bReadMat ? "" : " and RMSD calculation");
        get_index(&(top.atoms), ftp2fn_null(efNDX, NFILE, fnm),
                  1, &ifsize, &fitidx, &grpname);
        if (trx_out_fn)
        {
            fprintf(stderr, "\nSelect group for output:\n");
            get_index(&(top.atoms), ftp2fn_null(efNDX, NFILE, fnm),
                      1, &iosize, &outidx, &grpname);
            /* merge and convert both index groups: */
            /* first copy outidx to index. let outidx refer to elements in index */
            snew(index, iosize);
            isize = iosize;
            for (i = 0; i < iosize; i++)
            {
                index[i]  = outidx[i];
                outidx[i] = i;
            }
            /* now lookup elements from fitidx in index, add them if necessary
               and also let fitidx refer to elements in index */
            for (i = 0; i < ifsize; i++)
            {
                j = 0;
                while (j < isize && index[j] != fitidx[i])
                {
                    j++;
                }
                if (j >= isize)
                {
                    /* slow this way, but doesn't matter much */
                    isize++;
                    srenew(index, isize);
                }
                index[j]  = fitidx[i];
                fitidx[i] = j;
            }
        }
        else /* !trx_out_fn */
        {
            isize = ifsize;
            snew(index, isize);
            for (i = 0; i < ifsize; i++)
            {
                index[i]  = fitidx[i];
                fitidx[i] = i;
            }
        }
    }

    if (bReadTraj)
    {
        /* Loop over first coordinate file */
        fn = opt2fn("-f", NFILE, fnm);

        xx = read_whole_trj(fn, isize, index, skip, &nf, &time, oenv, bPBC, gpbc);
        output_env_conv_times(oenv, nf, time);
        if (!bRMSdist || bAnalyze)
        {
            /* Center all frames on zero */
            snew(mass, isize);
            for (i = 0; i < ifsize; i++)
            {
                mass[fitidx[i]] = top.atoms.atom[index[fitidx[i]]].m;
            }
            if (bFit)
            {
                for (i = 0; i < nf; i++)
                {
                    reset_x(ifsize, fitidx, isize, NULL, xx[i], mass);
                }
            }
        }
        if (bPBC)
        {
            gmx_rmpbc_done(gpbc);
        }
    }

    if (bReadMat)
    {
        fprintf(stderr, "Reading rms distance matrix ");
        read_xpm_matrix(opt2fn("-dm", NFILE, fnm), &readmat);
        fprintf(stderr, "\n");
        if (readmat[0].nx != readmat[0].ny)
        {
            gmx_fatal(FARGS, "Matrix (%dx%d) is not square",
                      readmat[0].nx, readmat[0].ny);
        }
        if (bReadTraj && bAnalyze && (readmat[0].nx != nf))
        {
            gmx_fatal(FARGS, "Matrix size (%dx%d) does not match the number of "
                      "frames (%d)", readmat[0].nx, readmat[0].ny, nf);
        }

        nf = readmat[0].nx;
        sfree(time);
        time        = readmat[0].axis_x;
        time_invfac = output_env_get_time_invfactor(oenv);
        for (i = 0; i < nf; i++)
        {
            time[i] *= time_invfac;
        }

        rms = init_mat(readmat[0].nx, method == m_diagonalize);
        convert_mat(&(readmat[0]), rms);

        nlevels = readmat[0].nmap;
    }
    else   /* !bReadMat */
    {
        rms  = init_mat(nf, method == m_diagonalize);
        nrms = ((gmx_int64_t)nf*((gmx_int64_t)nf-1))/2;
        if (!bRMSdist)
        {
            fprintf(stderr, "Computing %dx%d RMS deviation matrix\n", nf, nf);
            /* Initialize work array */
            snew(x1, isize);
            for (i1 = 0; i1 < nf; i1++)
            {
                for (i2 = i1+1; i2 < nf; i2++)
                {
                    for (i = 0; i < isize; i++)
                    {
                        copy_rvec(xx[i1][i], x1[i]);
                    }
                    if (bFit)
                    {
                        do_fit(isize, mass, xx[i2], x1);
                    }
                    rmsd = rmsdev(isize, mass, xx[i2], x1);
                    set_mat_entry(rms, i1, i2, rmsd);
                }
                nrms -= (gmx_int64_t) (nf-i1-1);
                fprintf(stderr, "\r# RMSD calculations left: " "%"GMX_PRId64 "   ", nrms);
            }
            sfree(x1);
        }
        else /* bRMSdist */
        {
            fprintf(stderr, "Computing %dx%d RMS distance deviation matrix\n", nf, nf);

            /* Initiate work arrays */
            snew(d1, isize);
            snew(d2, isize);
            for (i = 0; (i < isize); i++)
            {
                snew(d1[i], isize);
                snew(d2[i], isize);
            }
            for (i1 = 0; i1 < nf; i1++)
            {
                calc_dist(isize, xx[i1], d1);
                for (i2 = i1+1; (i2 < nf); i2++)
                {
                    calc_dist(isize, xx[i2], d2);
                    set_mat_entry(rms, i1, i2, rms_dist(isize, d1, d2));
                }
                nrms -= (nf-i1-1);
                fprintf(stderr, "\r# RMSD calculations left: " "%"GMX_PRId64 "   ", nrms);
            }
            /* Clean up work arrays */
            for (i = 0; (i < isize); i++)
            {
                sfree(d1[i]);
                sfree(d2[i]);
            }
            sfree(d1);
            sfree(d2);
        }
        fprintf(stderr, "\n\n");
    }
    ffprintf_gg(stderr, log, buf, "The RMSD ranges from %g to %g nm\n",
                rms->minrms, rms->maxrms);
    ffprintf_g(stderr, log, buf, "Average RMSD is %g\n", 2*rms->sumrms/(nf*(nf-1)));
    ffprintf_d(stderr, log, buf, "Number of structures for matrix %d\n", nf);
    ffprintf_g(stderr, log, buf, "Energy of the matrix is %g.\n", mat_energy(rms));
    if (bUseRmsdCut && (rmsdcut < rms->minrms || rmsdcut > rms->maxrms) )
    {
        fprintf(stderr, "WARNING: rmsd cutoff %g is outside range of rmsd values "
                "%g to %g\n", rmsdcut, rms->minrms, rms->maxrms);
    }
    if (bAnalyze && (rmsmin < rms->minrms) )
    {
        fprintf(stderr, "WARNING: rmsd minimum %g is below lowest rmsd value %g\n",
                rmsmin, rms->minrms);
    }
    if (bAnalyze && (rmsmin > rmsdcut) )
    {
        fprintf(stderr, "WARNING: rmsd minimum %g is above rmsd cutoff %g\n",
                rmsmin, rmsdcut);
    }

    /* Plot the rmsd distribution */
    rmsd_distribution(opt2fn("-dist", NFILE, fnm), rms, oenv);

    if (bBinary)
    {
        for (i1 = 0; (i1 < nf); i1++)
        {
            for (i2 = 0; (i2 < nf); i2++)
            {
                if (rms->mat[i1][i2] < rmsdcut)
                {
                    rms->mat[i1][i2] = 0;
                }
                else
                {
                    rms->mat[i1][i2] = 1;
                }
            }
        }
    }

    snew(clust.cl, nf);
    switch (method)
    {
        case m_linkage:
            /* Now sort the matrix and write it out again */
            gather(rms, rmsdcut, &clust);
            break;
        case m_diagonalize:
            /* Do a diagonalization */
            snew(eigenvalues, nf);
            snew(eigenvectors, nf*nf);
            memcpy(eigenvectors, rms->mat[0], nf*nf*sizeof(real));
            eigensolver(eigenvectors, nf, 0, nf, eigenvalues, rms->mat[0]);
            sfree(eigenvectors);

            fp = xvgropen(opt2fn("-ev", NFILE, fnm), "RMSD matrix Eigenvalues",
                          "Eigenvector index", "Eigenvalues (nm\\S2\\N)", oenv);
            for (i = 0; (i < nf); i++)
            {
                fprintf(fp, "%10d  %10g\n", i, eigenvalues[i]);
            }
            xvgrclose(fp);
            break;
        case m_monte_carlo:
            orig     = init_mat(rms->nn, FALSE);
            orig->nn = rms->nn;
            copy_t_mat(orig, rms);
            mc_optimize(log, rms, time, niter, nrandom, seed, kT,
                        opt2fn_null("-conv", NFILE, fnm), oenv);
            break;
        case m_jarvis_patrick:
            jarvis_patrick(rms->nn, rms->mat, M, P, bJP_RMSD ? rmsdcut : -1, &clust);
            break;
        case m_gromos:
            gromos(rms->nn, rms->mat, rmsdcut, &clust);
            break;
        default:
            gmx_fatal(FARGS, "DEATH HORROR unknown method \"%s\"", methodname[0]);
    }

    if (method == m_monte_carlo || method == m_diagonalize)
    {
        fprintf(stderr, "Energy of the matrix after clustering is %g.\n",
                mat_energy(rms));
    }

    if (bAnalyze)
    {
        if (minstruct > 1)
        {
            ncluster = plot_clusters(nf, rms->mat, &clust, minstruct);
        }
        else
        {
            mark_clusters(nf, rms->mat, rms->maxrms, &clust);
        }
        init_t_atoms(&useatoms, isize, FALSE);
        snew(usextps, isize);
        useatoms.resinfo = top.atoms.resinfo;
        for (i = 0; i < isize; i++)
        {
            useatoms.atomname[i]    = top.atoms.atomname[index[i]];
            useatoms.atom[i].resind = top.atoms.atom[index[i]].resind;
            useatoms.nres           = max(useatoms.nres, useatoms.atom[i].resind+1);
            copy_rvec(xtps[index[i]], usextps[i]);
        }
        useatoms.nr = isize;
        analyze_clusters(nf, &clust, rms->mat, isize, &useatoms, usextps, mass, xx, time,
                         ifsize, fitidx, iosize, outidx,
                         bReadTraj ? trx_out_fn : NULL,
                         opt2fn_null("-sz", NFILE, fnm),
                         opt2fn_null("-tr", NFILE, fnm),
                         opt2fn_null("-ntr", NFILE, fnm),
                         opt2fn_null("-clid", NFILE, fnm),
                         bAverage, write_ncl, write_nst, rmsmin, bFit, log,
                         rlo_bot, rhi_bot, oenv);
    }
    gmx_ffclose(log);

    if (bBinary && !bAnalyze)
    {
        /* Make the clustering visible */
        for (i2 = 0; (i2 < nf); i2++)
        {
            for (i1 = i2+1; (i1 < nf); i1++)
            {
                if (rms->mat[i1][i2])
                {
                    rms->mat[i1][i2] = rms->maxrms;
                }
            }
        }
    }

    fp = opt2FILE("-o", NFILE, fnm, "w");
    fprintf(stderr, "Writing rms distance/clustering matrix ");
    if (bReadMat)
    {
        write_xpm(fp, 0, readmat[0].title, readmat[0].legend, readmat[0].label_x,
                  readmat[0].label_y, nf, nf, readmat[0].axis_x, readmat[0].axis_y,
                  rms->mat, 0.0, rms->maxrms, rlo_top, rhi_top, &nlevels);
    }
    else
    {
        sprintf(buf, "Time (%s)", output_env_get_time_unit(oenv));
        sprintf(title, "RMS%sDeviation / Cluster Index",
                bRMSdist ? " Distance " : " ");
        if (minstruct > 1)
        {
            write_xpm_split(fp, 0, title, "RMSD (nm)", buf, buf,
                            nf, nf, time, time, rms->mat, 0.0, rms->maxrms, &nlevels,
                            rlo_top, rhi_top, 0.0, (real) ncluster,
                            &ncluster, TRUE, rlo_bot, rhi_bot);
        }
        else
        {
            write_xpm(fp, 0, title, "RMSD (nm)", buf, buf,
                      nf, nf, time, time, rms->mat, 0.0, rms->maxrms,
                      rlo_top, rhi_top, &nlevels);
        }
    }
    fprintf(stderr, "\n");
    gmx_ffclose(fp);
    if (NULL != orig)
    {
        fp = opt2FILE("-om", NFILE, fnm, "w");
        sprintf(buf, "Time (%s)", output_env_get_time_unit(oenv));
        sprintf(title, "RMS%sDeviation", bRMSdist ? " Distance " : " ");
        write_xpm(fp, 0, title, "RMSD (nm)", buf, buf,
                  nf, nf, time, time, orig->mat, 0.0, orig->maxrms,
                  rlo_top, rhi_top, &nlevels);
        gmx_ffclose(fp);
        done_mat(&orig);
        sfree(orig);
    }
    /* now show what we've done */
    do_view(oenv, opt2fn("-o", NFILE, fnm), "-nxy");
    do_view(oenv, opt2fn_null("-sz", NFILE, fnm), "-nxy");
    if (method == m_diagonalize)
    {
        do_view(oenv, opt2fn_null("-ev", NFILE, fnm), "-nxy");
    }
    do_view(oenv, opt2fn("-dist", NFILE, fnm), "-nxy");
    if (bAnalyze)
    {
        do_view(oenv, opt2fn_null("-tr", NFILE, fnm), "-nxy");
        do_view(oenv, opt2fn_null("-ntr", NFILE, fnm), "-nxy");
        do_view(oenv, opt2fn_null("-clid", NFILE, fnm), "-nxy");
    }
    do_view(oenv, opt2fn_null("-conv", NFILE, fnm), NULL);

    return 0;
}
