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

#ifndef GMX_GMXPREPROCESS_PGUTIL_H
#define GMX_GMXPREPROCESS_PGUTIL_H

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

/* Search an atom in array of pointers to strings, starting from start
 * if type starts with '-' then searches backwards from start.
 * bondtype is only used for printing the error/warning string,
 * when bondtype="check" no error/warning is issued.
 * When bAllowMissing=FALSE an fatal error is issued, otherwise a warning.
 */
atom_id search_atom(const char *type, int start,
                    t_atoms *atoms,
                    const char *bondtype, gmx_bool bAllowMissing);

/* Similar to search_atom, but this routine searches for the specified
 * atom in residue resind.
 */
atom_id
search_res_atom(const char *type, int resind,
                t_atoms *atoms,
                const char *bondtype, gmx_bool bAllowMissing);


void set_at(t_atom *at, real m, real q, int type, int resind);


#ifdef __cplusplus
}
#endif

#endif
