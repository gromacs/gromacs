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

#ifndef GMX_TOOLS_COMPARE_H
#define GMX_TOOLS_COMPARE_H

struct gmx_groups_t;
struct gmx_output_env_t;
struct t_inputrec;
struct t_state;
struct t_topology;

/* Routines for comparing data structures from non-trajectory binary
   file formats (e.g. as used by gmx check). */

void
cmp_real(FILE *fp, const char *s, int index, real i1, real i2, real ftol, real abstol);

void cmp_top(FILE *fp, t_topology *t1, t_topology *t2, real ftol, real abstol);

void cmp_groups(FILE *fp, gmx_groups_t *g0, gmx_groups_t *g1,
                int natoms0, int natoms1);

void cmp_inputrec(FILE *fp, t_inputrec *ir1, t_inputrec *ir2, real ftol, real abstol);

void comp_pull_AB(FILE *fp, pull_params_t *pull, real ftol, real abstol);

void comp_state(t_state *st1, t_state *st2, gmx_bool bRMSD, real ftol, real abstol);

void
comp_trx(const gmx_output_env_t *oenv, const char *fn1, const char *fn2,
         gmx_bool bRMSD, real ftol, real abstol);
/* Compare two binary trajectory files */

void
comp_enx(const char *fn1, const char *fn2, real ftol, real abstol,
         const char *lastener);
/* Compare two binary energy files */

#endif
