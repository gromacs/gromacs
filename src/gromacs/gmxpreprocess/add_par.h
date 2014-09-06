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

#ifndef GMX_GMXPREPROCESS_ADD_PAR_H
#define GMX_GMXPREPROCESS_ADD_PAR_H

#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

void add_param(t_params *ps, int ai, int aj, real *c, char *s);

void add_imp_param(t_params *ps, int ai, int aj, int ak, int al,
                   real c0, real c1, char *s);

void add_dih_param(t_params *ps, int ai, int aj, int ak, int al,
                   real c0, real c1, real c2, char *s);

void add_cmap_param(t_params *ps, int ai, int aj, int ak, int al, int am,
                    char *s);

void add_vsite2_atoms(t_params *ps, int ai, int aj, int ak);

void add_vsite3_atoms(t_params *ps, int ai, int aj, int ak, int al,
                      gmx_bool bSwapParity);

void add_vsite2_param(t_params *ps, int ai, int aj, int ak, real c0);

void add_vsite3_param(t_params *ps, int ai, int aj, int ak, int al,
                      real c0, real c1);

void add_vsite4_atoms(t_params *ps, int ai, int aj, int ak, int al,
                      int am);

int search_jtype(t_restp *rp, char *name, gmx_bool bFirstRes);

#ifdef __cplusplus
}
#endif

#endif
