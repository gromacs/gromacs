/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_FILEIO_GROIO_H
#define GMX_FILEIO_GROIO_H

#include <stdio.h>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_mtop_t;
struct t_atoms;
struct t_topology;
struct t_trxframe;

void get_coordnum(const char *infile, int *natoms);
void gmx_gro_read_conf(const char *infile,
                       struct t_topology *top, rvec x[], rvec *v, matrix box);

gmx_bool gro_next_x_or_v(FILE *status, struct t_trxframe *fr);
int gro_first_x_or_v(FILE *status, struct t_trxframe *fr);
/* read first/next x and/or v frame from gro file */

void write_hconf_indexed_p(FILE *out, const char *title, const t_atoms *atoms,
                           int nx, const int index[], int ndec,
                           const rvec *x, const rvec *v, const matrix box);

void write_hconf_mtop(FILE *out, const char *title, struct gmx_mtop_t *mtop, int pr,
                      const rvec *x, const rvec *v, const matrix box);

void write_hconf_p(FILE *out, const char *title, const t_atoms *atoms, int ndec,
                   const rvec *x, const rvec *v, const matrix box);
/* Write a Gromos file with precision ndec: number of decimal places in x,
 * v has one place more. */

void write_conf_p(const char *outfile, const char *title,
                  const t_atoms *atoms, int pr,
                  const rvec *x, const rvec *v, const matrix box);

#ifdef __cplusplus
}
#endif

#endif
