/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
#ifndef GMX_TOPOLOGY_ATOMPROP_H
#define GMX_TOPOLOGY_ATOMPROP_H

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Abstract type for the atom property database */
typedef struct gmx_atomprop *gmx_atomprop_t;

enum {
    epropMass, epropVDW, epropDGsol, epropElectroneg, epropElement,
    epropNR
};

gmx_atomprop_t gmx_atomprop_init(void);
/* Initializes and returns the atom properties struct */

void gmx_atomprop_destroy(gmx_atomprop_t aps);
/* Get rid of memory after use */

char *gmx_atomprop_element(gmx_atomprop_t aps, int atomnumber);

int gmx_atomprop_atomnumber(gmx_atomprop_t aps, const char *element);

gmx_bool gmx_atomprop_query(gmx_atomprop_t aps,
                            int eprop, const char *resnm, const char *atomnm,
                            real *value);
/* Extract a value from the database. Returns TRUE on succes,
 * FALSE otherwise. In the latter case, value is a deafult value.
 * The first time this function is called for this property
 * the database will be read.
 */

#ifdef __cplusplus
}
#endif

#endif
