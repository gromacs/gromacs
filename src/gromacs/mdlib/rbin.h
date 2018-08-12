/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014,2015,2018, by the GROMACS development team, led by
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
#ifndef GMX_MDLIB_RBIN_H
#define GMX_MDLIB_RBIN_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct t_commrec;

typedef struct {
    int     nreal;
    int     maxreal;
    double *rbuf;
} t_bin;

t_bin *mk_bin();
/* Create a real bin */

void destroy_bin(t_bin *b);
/* Destroy the bin structure */

void reset_bin(t_bin *b);
/* Reset number of entries to zero */

int add_binr(t_bin *b, int nr, const real r[]);
int add_binr(t_bin *b, gmx::ArrayRef<const real> r);
int add_bind(t_bin *b, int nr, const double r[]);
/* Add reals to the bin. Returns index */

void sum_bin(t_bin *b, const t_commrec *cr);
/* Globally sum the reals in the bin */

void extract_binr(t_bin *b, int index, int nr, real r[]);
void extract_binr(t_bin *b, int index, gmx::ArrayRef<real> r);
void extract_bind(t_bin *b, int index, int nr, double r[]);
/* Extract values from the bin, starting from index (see add_bin) */

#endif
