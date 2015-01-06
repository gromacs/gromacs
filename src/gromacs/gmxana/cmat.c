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

#include "cmat.h"

#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

t_mat *init_mat(int n1, gmx_bool b1D)
{
    t_mat *m;

    snew(m, 1);
    m->n1     = n1;
    m->nn     = 0;
    m->b1D    = b1D;
    m->maxrms = 0;
    m->minrms = 1e20;
    m->sumrms = 0;
    m->mat    = mk_matrix(n1, n1, b1D);

    snew(m->erow, n1);
    snew(m->m_ind, n1);
    reset_index(m);

    return m;
}

void copy_t_mat(t_mat *dst, t_mat *src)
{
    int i, j;

    if (dst->nn != src->nn)
    {
        fprintf(stderr, "t_mat structures not identical in size dst %d src %d\n", dst->nn, src->nn);
        return;
    }
    dst->maxrms = src->maxrms;
    dst->minrms = src->minrms;
    dst->sumrms = src->sumrms;
    for (i = 0; (i < src->nn); i++)
    {
        for (j = 0; (j < src->nn); j++)
        {
            dst->mat[i][j] = src->mat[i][j];
        }
        dst->erow[i]  = src->erow[i];
        dst->m_ind[i] = src->m_ind[i];
    }
}

void enlarge_mat(t_mat *m, int deltan)
{
    int i, j;

    srenew(m->erow, m->nn+deltan);
    srenew(m->m_ind, m->nn+deltan);
    srenew(m->mat, m->nn+deltan);

    /* Reallocate existing rows in the matrix, and set them to zero */
    for (i = 0; (i < m->nn); i++)
    {
        srenew(m->mat[i], m->nn+deltan);
        for (j = m->nn; (j < m->nn+deltan); j++)
        {
            m->mat[i][j] = 0;
        }
    }
    /* Allocate new rows of the matrix, set energies to zero */
    for (i = m->nn; (i < m->nn+deltan); i++)
    {
        m->erow[i]  = 0;
        snew(m->mat[i], m->nn+deltan);
    }
    m->nn += deltan;
}

void reset_index(t_mat *m)
{
    int i;

    for (i = 0; (i < m->n1); i++)
    {
        m->m_ind[i] = i;
    }
}

void set_mat_entry(t_mat *m, int i, int j, real val)
{
    m->mat[i][j] = m->mat[j][i] = val;
    m->maxrms    = max(m->maxrms, val);
    if (j != i)
    {
        m->minrms  = min(m->minrms, val);
    }
    m->sumrms   += val;
    m->nn        = max(m->nn, max(j+1, i+1));
}

void done_mat(t_mat **m)
{
    done_matrix((*m)->n1, &((*m)->mat));
    sfree((*m)->m_ind);
    sfree((*m)->erow);
    sfree(*m);
    *m = NULL;
}

real mat_energy(t_mat *m)
{
    int  j;
    real emat = 0;

    for (j = 0; (j < m->nn-1); j++)
    {
        emat += sqr(m->mat[j][j+1]);
    }
    return emat;
}

void swap_rows(t_mat *m, int iswap, int jswap)
{
    real *tmp, ttt;
    int   i, itmp;

    /* Swap indices */
    itmp            = m->m_ind[iswap];
    m->m_ind[iswap] = m->m_ind[jswap];
    m->m_ind[jswap] = itmp;

    /* Swap rows (since the matrix is an array of pointers) */
    tmp           = m->mat[iswap];
    m->mat[iswap] = m->mat[jswap];
    m->mat[jswap] = tmp;

    /* Swap columns */
    for (i = 0; (i < m->nn); i++)
    {
        ttt              = m->mat[i][iswap];
        m->mat[i][iswap] = m->mat[i][jswap];
        m->mat[i][jswap] = ttt;
    }
}

void swap_mat(t_mat *m)
{
    t_mat *tmp;
    int    i, j;

    tmp = init_mat(m->nn, FALSE);
    for (i = 0; (i < m->nn); i++)
    {
        for (j = 0; (j < m->nn); j++)
        {
            tmp->mat[m->m_ind[i]][m->m_ind[j]] = m->mat[i][j];
        }
    }
    /*tmp->mat[i][j] =  m->mat[m->m_ind[i]][m->m_ind[j]]; */
    for (i = 0; (i < m->nn); i++)
    {
        for (j = 0; (j < m->nn); j++)
        {
            m->mat[i][j] = tmp->mat[i][j];
        }
    }
    done_mat(&tmp);
}

void low_rmsd_dist(const char *fn, real maxrms, int nn, real **mat,
                   const output_env_t oenv)
{
    FILE   *fp;
    int     i, j, *histo, x;
    real    fac;

    fac = 100/maxrms;
    snew(histo, 101);
    for (i = 0; i < nn; i++)
    {
        for (j = i+1; j < nn; j++)
        {
            x = (int)(fac*mat[i][j]+0.5);
            if (x <= 100)
            {
                histo[x]++;
            }
        }
    }

    fp = xvgropen(fn, "RMS Distribution", "RMS (nm)", "a.u.", oenv);
    for (i = 0; (i < 101); i++)
    {
        fprintf(fp, "%10g  %10d\n", i/fac, histo[i]);
    }
    xvgrclose(fp);
    sfree(histo);
}

void rmsd_distribution(const char *fn, t_mat *rms, const output_env_t oenv)
{
    low_rmsd_dist(fn, rms->maxrms, rms->nn, rms->mat, oenv);
}

t_clustid *new_clustid(int n1)
{
    t_clustid *c;
    int        i;

    snew(c, n1);
    for (i = 0; (i < n1); i++)
    {
        c[i].conf  = i;
        c[i].clust = i;
    }
    return c;
}
