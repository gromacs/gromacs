/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_RESIDUETYPES_H
#define GMX_TOPOLOGY_RESIDUETYPES_H

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct gmx_residuetype_t gmx_residuetype_t;

int
gmx_residuetype_init(gmx_residuetype_t **rt);

int
gmx_residuetype_destroy(gmx_residuetype_t *rt);

int
gmx_residuetype_get_type(gmx_residuetype_t *rt, const char *resname, const char **p_restype);

int
gmx_residuetype_add(gmx_residuetype_t *rt, const char *newresname, const char *newrestype);

int
gmx_residuetype_get_alltypes(gmx_residuetype_t   *rt,
                             const char        ***p_typenames,
                             int                 *ntypes);

gmx_bool
gmx_residuetype_is_protein(gmx_residuetype_t *rt, const char *resnm);

gmx_bool
gmx_residuetype_is_dna(gmx_residuetype_t *rt, const char *resnm);

gmx_bool
gmx_residuetype_is_rna(gmx_residuetype_t *rt, const char *resnm);

int
gmx_residuetype_get_size(gmx_residuetype_t *rt);

int
gmx_residuetype_get_index(gmx_residuetype_t *rt, const char *resnm);

const char *
gmx_residuetype_get_name(gmx_residuetype_t *rt, int index);

#ifdef __cplusplus
}
#endif

#endif
