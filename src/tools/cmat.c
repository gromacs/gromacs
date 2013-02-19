/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "cmat.h"
#include "smalloc.h"
#include "macros.h"
#include "xvgr.h"
#include "matio.h"
#include "futil.h"

t_mat *init_mat(int n1, gmx_bool b1D)
{
    t_mat *m;

    snew(m, 1);
    m->n1     = n1;
    m->nn     = 0;
    m->b1D    = b1D;
    m->emat   = 0;
    m->maxrms = 0;
    m->minrms = 1e20;
    m->sumrms = 0;
    m->nn     = 0;
    m->mat    = mk_matrix(n1, n1, b1D);

    snew(m->erow, n1);
    snew(m->m_ind, n1);
    reset_index(m);

    return m;
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

real row_energy(int nn, int row, real *mat)
{
    real re = 0;
    int  i;

    for (i = 0; (i < nn); i++)
    {
        re += abs(i-row)*mat[i];
    }
    return re/nn;
}

real mat_energy(t_mat *m)
{
    real re, retot;
    int  j, jj;

    retot = 0;
    for (j = 0; (j < m->nn); j++)
    {
        jj         = m->m_ind[j];
        re         = row_energy(m->nn, jj, m->mat[j]);
        m->erow[j] = re;
        retot     += re;
    }
    m->emat = retot/m->nn;
    return m->emat;
}

void swap_rows(t_mat *m, int isw, int jsw)
{
    real *tmp, ttt;
    int   i;

    /* Swap rows */
    tmp         = m->mat[isw];
    m->mat[isw] = m->mat[jsw];
    m->mat[jsw] = tmp;
    /* Swap columns */
    for (i = 0; (i < m->nn); i++)
    {
        ttt            = m->mat[isw][i];
        m->mat[isw][i] = m->mat[jsw][i];
        m->mat[jsw][i] = ttt;
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
    ffclose(fp);
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
