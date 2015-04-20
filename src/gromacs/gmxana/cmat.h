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

#ifndef _cmat_h
#define _cmat_h

#include "gromacs/legacyheaders/typedefs.h"

typedef struct {
    int  i, j;
    real dist;
} t_dist;

typedef struct {
    int  conf, clust;
} t_clustid;

typedef struct {
    int      n1, nn;
    int     *m_ind;
    gmx_bool b1D;
    real     minrms, maxrms, sumrms;
    real    *erow;
    real   **mat;
} t_mat;

/* The matrix is indexed using the matrix index */
#define EROW(m, i)  m->erow[i]

extern t_mat *init_mat(int n1, gmx_bool b1D);

extern void copy_t_mat(t_mat *dst, t_mat *src);

extern void enlarge_mat(t_mat *m, int deltan);

extern void reset_index(t_mat *m);

extern void swap_rows(t_mat *m, int iswap, int jswap);

extern void set_mat_entry(t_mat *m, int i, int j, real val);

extern void done_mat(t_mat **m);

extern real mat_energy(t_mat *mat);

extern void swap_mat(t_mat *m);

extern void low_rmsd_dist(const char *fn, real maxrms, int nn, real **mat,
                          const output_env_t oenv);

extern void rmsd_distribution(const char *fn, t_mat *m, const output_env_t oenv);

extern t_clustid *new_clustid(int n1);

#endif
