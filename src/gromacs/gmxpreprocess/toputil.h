/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_TOPUTIL_H
#define GMX_GMXPREPROCESS_TOPUTIL_H

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"

#ifdef __cplusplus
extern "C" {
#endif

/* UTILITIES */

int name2index(char *str, char ***typenames, int ntypes);

void pr_alloc (int extra, t_params *pr);

void set_p_string(t_param *p, const char *s);

void cp_param(t_param *dest, t_param *src);

void add_param_to_list(t_params *list, t_param *b);

/* INITIATE */

void init_plist(t_params plist[]);

void init_molinfo(t_molinfo *mol);

/* FREE */
void done_mi(t_molinfo *mi);

/* PRINTING */

void print_blocka(FILE *out, const char *szName, const char *szIndex,
                  const char *szA, t_blocka *block);

void print_atoms(FILE *out, gpp_atomtype_t atype, t_atoms *at, int *cgnr,
                 gmx_bool bRTPresname);

void print_bondeds(FILE *out, int natoms, directive d,
                   int ftype, int fsubtype, t_params plist[]);

void print_excl(FILE *out, int natoms, t_excls excls[]);

#ifdef __cplusplus
}
#endif

#endif
