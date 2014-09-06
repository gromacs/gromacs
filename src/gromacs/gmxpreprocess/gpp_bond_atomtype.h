/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_GPP_BONDATOMTYPE_H
#define GMX_GMXPREPROCESS_GPP_BONDATOMTYPE_H

#include <stdio.h>

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gpp_bondatomtype *t_bond_atomtype;

int get_bond_atomtype_type(char *str, t_bond_atomtype at);
/* Return atomtype corresponding to case-insensitive str
   or NOTSET if not found */

char *get_bond_atomtype_name(int nt, t_bond_atomtype at);
/* Return name corresponding to atomtype nt, or NULL if not found */

t_bond_atomtype init_bond_atomtype(void);
/* Return a new atomtype structure */

void done_bond_atomtype(t_bond_atomtype *at);
/* Free the memory in the structure */

void add_bond_atomtype(t_bond_atomtype at, struct t_symtab *tab,
                       char *name);
/* Add a complete new atom type to an existing atomtype structure */

#ifdef __cplusplus
}
#endif

#endif
