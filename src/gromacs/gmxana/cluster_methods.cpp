/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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

#include "cluster_methods.h"

#include <cmath>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/cmat.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformintdistribution.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

void mc_optimize(FILE* log, t_mat* m, real* time, int maxiter, int nrandom, int seed, real kT, const char* conv, gmx_output_env_t* oenv)
{
    FILE*  fp = nullptr;
    real   ecur, enext, emin, prob, enorm;
    int    i, j, iswap, jswap, nn, nuphill = 0;
    t_mat* minimum;

    if (seed == 0)
    {
        seed = static_cast<int>(gmx::makeRandomSeed());
    }
    gmx::DefaultRandomEngine rng(seed);

    if (m->n1 != m->nn)
    {
        fprintf(stderr, "Can not do Monte Carlo optimization with a non-square matrix.\n");
        return;
    }
    printf("\nDoing Monte Carlo optimization to find the smoothest trajectory\n");
    printf("by reordering the frames to minimize the path between the two structures\n");
    printf("that have the largest pairwise RMSD.\n");
    printf("Using random seed %d.\n", seed);

    iswap = jswap = -1;
    enorm         = m->mat[0][0];
    for (i = 0; (i < m->n1); i++)
    {
        for (j = 0; (j < m->nn); j++)
        {
            if (m->mat[i][j] > enorm)
            {
                enorm = m->mat[i][j];
                iswap = i;
                jswap = j;
            }
        }
    }
    if ((iswap == -1) || (jswap == -1))
    {
        fprintf(stderr, "Matrix contains identical values in all fields\n");
        return;
    }
    swap_rows(m, 0, iswap);
    swap_rows(m, m->n1 - 1, jswap);
    emin = ecur = mat_energy(m);
    printf("Largest distance %g between %d and %d. Energy: %g.\n", enorm, iswap, jswap, emin);

    nn = m->nn;

    /* Initiate and store global minimum */
    minimum     = init_mat(nn, m->b1D);
    minimum->nn = nn;
    copy_t_mat(minimum, m);

    if (nullptr != conv)
    {
        fp = xvgropen(conv, "Convergence of the MC optimization", "Energy", "Step", oenv);
    }

    gmx::UniformIntDistribution<int>   intDistNN(1, nn - 2); // [1,nn-2]
    gmx::UniformRealDistribution<real> realDistOne;          // [0,1)

    for (i = 0; (i < maxiter); i++)
    {
        /* Generate new swapping candidates */
        do
        {
            iswap = intDistNN(rng);
            jswap = intDistNN(rng);
        } while ((iswap == jswap) || (iswap >= nn - 1) || (jswap >= nn - 1));

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
            prob = std::exp(-(enext - ecur) / (enorm * kT));
        }

        if (prob == 1 || realDistOne(rng) < prob)
        {
            if (enext > ecur)
            {
                nuphill++;
            }

            fprintf(log, "Iter: %d Swapped %4d and %4d (energy: %g prob: %g)\n", i, iswap, jswap, enext, prob);
            if (nullptr != fp)
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
    fprintf(log, "%d uphill steps were taken during optimization\n", nuphill);

    /* Now swap the matrix to get it into global minimum mode */
    copy_t_mat(m, minimum);

    fprintf(log, "Global minimum energy %g\n", mat_energy(minimum));
    fprintf(log, "Global minimum energy %g\n", mat_energy(m));
    fprintf(log, "Swapped time and frame indices and RMSD to next neighbor:\n");
    for (i = 0; (i < m->nn); i++)
    {
        fprintf(log,
                "%10g  %5d  %10g\n",
                time[m->m_ind[i]],
                m->m_ind[i],
                (i < m->nn - 1) ? m->mat[m->m_ind[i]][m->m_ind[i + 1]] : 0);
    }

    if (nullptr != fp)
    {
        xvgrclose(fp);
    }
}

static bool rms_dist_comp(const t_dist& a, const t_dist& b)
{
    return a.dist < b.dist;
}

static bool clust_id_comp(const t_clustid& a, const t_clustid& b)
{
    return a.clust < b.clust;
}

static bool nrnb_comp(const t_nnb& a, const t_nnb& b)
{
    /* return b<a, we want highest first */
    return b.nr < a.nr;
}

void gather(t_mat* m, real cutoff, t_clusters* clust)
{
    t_clustid* c;
    t_dist*    d;
    int        i, j, k, nn, cid, n1, diff;
    gmx_bool   bChange;

    /* First we sort the entries in the RMSD matrix */
    n1 = m->nn;
    nn = ((n1 - 1) * n1) / 2;
    snew(d, nn);
    for (i = k = 0; (i < n1); i++)
    {
        for (j = i + 1; (j < n1); j++, k++)
        {
            d[k].i    = i;
            d[k].j    = j;
            d[k].dist = m->mat[i][j];
        }
    }
    if (k != nn)
    {
        gmx_incons("gather algorithm");
    }
    std::sort(d, d + nn, rms_dist_comp);

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
    } while (bChange);
    fprintf(stderr, "\nSorting and renumbering clusters\n");
    /* Sort on cluster number */
    std::sort(c, c + n1, clust_id_comp);

    /* Renumber clusters */
    cid = 1;
    for (k = 1; k < n1; k++)
    {
        if (c[k].clust != c[k - 1].clust)
        {
            c[k - 1].clust = cid;
            cid++;
        }
        else
        {
            c[k - 1].clust = cid;
        }
    }
    c[k - 1].clust = cid;
    if (debug)
    {
        for (k = 0; (k < n1); k++)
        {
            fprintf(debug, "Cluster index for conformation %d: %d\n", c[k].conf, c[k].clust);
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

static gmx_bool jp_same(int** nnb, int i, int j, int P)
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

void jarvis_patrick(int n1, real** mat, int M, int P, real rmsdcut, t_clusters* clust)
{
    t_dist*    row;
    t_clustid* c;
    int**      nnb;
    int        i, j, k, cid, diff, maxval;
    gmx_bool   bChange;
    real**     mcpy = nullptr;

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
        std::sort(row, row + n1, rms_dist_comp);
        if (M > 0)
        {
            /* Put the M nearest neighbors in the list */
            snew(nnb[i], M + 1);
            for (j = k = 0; (k < M) && (j < n1) && (mat[i][row[j].j] < rmsdcut); j++)
            {
                if (row[j].j != i)
                {
                    nnb[i][k] = row[j].j;
                    k++;
                }
            }
            nnb[i][k] = -1;
        }
        else
        {
            /* Put all neighbors nearer than rmsdcut in the list */
            maxval = 0;
            k      = 0;
            for (j = 0; (j < n1) && (mat[i][row[j].j] < rmsdcut); j++)
            {
                if (row[j].j != i)
                {
                    if (k >= maxval)
                    {
                        maxval += 10;
                        srenew(nnb[i], maxval);
                    }
                    nnb[i][k] = row[j].j;
                    k++;
                }
            }
            if (k == maxval)
            {
                srenew(nnb[i], maxval + 1);
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
        for (j = i + 1; j < n1; j++)
        {
            mcpy[i][j] = static_cast<real>(jp_same(nnb, i, j, P));
        }
    }
    do
    {
        fprintf(stderr, "*");
        bChange = FALSE;
        for (i = 0; i < n1; i++)
        {
            for (j = i + 1; j < n1; j++)
            {
                if (mcpy[i][j] != 0.0F)
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
    } while (bChange);

    fprintf(stderr, "\nSorting and renumbering clusters\n");
    /* Sort on cluster number */
    std::sort(c, c + n1, clust_id_comp);

    /* Renumber clusters */
    cid = 1;
    for (k = 1; k < n1; k++)
    {
        if (c[k].clust != c[k - 1].clust)
        {
            c[k - 1].clust = cid;
            cid++;
        }
        else
        {
            c[k - 1].clust = cid;
        }
    }
    c[k - 1].clust = cid;
    clust->ncl     = cid;
    for (k = 0; k < n1; k++)
    {
        clust->cl[c[k].conf] = c[k].clust;
    }
    if (debug)
    {
        for (k = 0; (k < n1); k++)
        {
            fprintf(debug, "Cluster index for conformation %d: %d\n", c[k].conf, c[k].clust);
        }
    }

    done_matrix(n1, &mcpy);

    sfree(c);
    for (i = 0; (i < n1); i++)
    {
        sfree(nnb[i]);
    }
    sfree(nnb);
}

static void dump_nnb(FILE* fp, const char* title, int n1, t_nnb* nnb)
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

void gromos(int n1, real** mat, real rmsdcut, t_clusters* clust)
{
    t_nnb* nnb;
    int    i, j, k, j1, maxval;

    /* Put all neighbors nearer than rmsdcut in the list */
    fprintf(stderr, "Making list of neighbors within cutoff ");
    snew(nnb, n1);
    for (i = 0; (i < n1); i++)
    {
        maxval = 0;
        k      = 0;
        /* put all neighbors within cut-off in list */
        for (j = 0; j < n1; j++)
        {
            if (mat[i][j] < rmsdcut)
            {
                if (k >= maxval)
                {
                    maxval += 10;
                    srenew(nnb[i].nb, maxval);
                }
                nnb[i].nb[k] = j;
                k++;
            }
        }
        /* store nr of neighbors, we'll need that */
        nnb[i].nr = k;
        if (i % (1 + n1 / 100) == 0)
        {
            fprintf(stderr, "%3d%%\b\b\b\b", (i * 100 + 1) / n1);
        }
    }
    fprintf(stderr, "%3d%%\n", 100);

    /* sort neighbor list on number of neighbors, largest first */
    std::sort(nnb, nnb + n1, nrnb_comp);

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
        std::sort(nnb, nnb + i, nrnb_comp);

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

    clust->ncl = k - 1;
}
