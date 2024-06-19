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

#include "cmat.h"

#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/fileio/matio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/functions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"

t_mat* init_mat(int n1, gmx_bool b1D)
{
    t_mat* m;

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

void copy_t_mat(t_mat* dst, t_mat* src)
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

void reset_index(t_mat* m)
{
    int i;

    for (i = 0; (i < m->n1); i++)
    {
        m->m_ind[i] = i;
    }
}

void set_mat_entry(t_mat* m, int i, int j, real val)
{
    m->mat[i][j] = m->mat[j][i] = val;
    m->maxrms                   = std::max(m->maxrms, val);
    if (j != i)
    {
        m->minrms = std::min(m->minrms, val);
    }
    m->sumrms += val;
    m->nn = std::max(m->nn, std::max(j + 1, i + 1));
}

void done_mat(t_mat** m)
{
    done_matrix((*m)->n1, &((*m)->mat));
    sfree((*m)->m_ind);
    sfree((*m)->erow);
    sfree(*m);
    *m = nullptr;
}

real mat_energy(t_mat* m)
{
    int  j;
    real emat = 0;

    for (j = 0; (j < m->nn - 1); j++)
    {
        emat += gmx::square(m->mat[j][j + 1]);
    }
    return emat;
}

void swap_rows(t_mat* m, int iswap, int jswap)
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

void low_rmsd_dist(const char* fn, real maxrms, int nn, real** mat, const gmx_output_env_t* oenv)
{
    FILE* fp;
    int   i, j, *histo, x;
    real  fac;

    fac = 100 / maxrms;
    snew(histo, 101);
    for (i = 0; i < nn; i++)
    {
        for (j = i + 1; j < nn; j++)
        {
            x = gmx::roundToInt(fac * mat[i][j]);
            if (x <= 100)
            {
                histo[x]++;
            }
        }
    }

    fp = xvgropen(fn, "RMS Distribution", "RMS (nm)", "counts", oenv);
    for (i = 0; (i < 101); i++)
    {
        fprintf(fp, "%10g  %10d\n", i / fac, histo[i]);
    }
    xvgrclose(fp);
    sfree(histo);
}

void rmsd_distribution(const char* fn, t_mat* rms, const gmx_output_env_t* oenv)
{
    low_rmsd_dist(fn, rms->maxrms, rms->nn, rms->mat, oenv);
}

t_clustid* new_clustid(int n1)
{
    t_clustid* c;
    int        i;

    snew(c, n1);
    for (i = 0; (i < n1); i++)
    {
        c[i].conf  = i;
        c[i].clust = i;
    }
    return c;
}
