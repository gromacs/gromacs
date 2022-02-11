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
#ifndef GMX_MDLIB_RBIN_H
#define GMX_MDLIB_RBIN_H

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/real.h"

struct t_commrec;

typedef struct
{
    int     nreal;
    int     maxreal;
    double* rbuf;
} t_bin;

t_bin* mk_bin();
/* Create a real bin */

void destroy_bin(t_bin* b);
/* Destroy the bin structure */

void reset_bin(t_bin* b);
/* Reset number of entries to zero */

int add_binr(t_bin* b, int nr, const real r[]);
int add_binr(t_bin* b, gmx::ArrayRef<const real> r);
int add_bind(t_bin* b, int nr, const double r[]);
int add_bind(t_bin* b, gmx::ArrayRef<const double> r);
/* Add reals to the bin. Returns index */

void sum_bin(t_bin* b, const t_commrec* cr);
/* Globally sum the reals in the bin */

void extract_binr(t_bin* b, int index, int nr, real r[]);
void extract_binr(t_bin* b, int index, gmx::ArrayRef<real> r);
void extract_bind(t_bin* b, int index, int nr, double r[]);
void extract_bind(t_bin* b, int index, gmx::ArrayRef<double> r);
/* Extract values from the bin, starting from index (see add_bin) */

#endif
