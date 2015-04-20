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

#ifndef GMX_GMXPREPROCESS_RESALL_H
#define GMX_GMXPREPROCESS_RESALL_H

#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/gmxpreprocess/hackblock.h"
#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

char *search_rtp(const char *key, int nrtp, t_restp rtp[]);
/* Search for an entry in the rtp database, returns the rtp residue name.
 * A mismatch of one character is allowed, if there is only one nearly
 * matching entry in the database, a warning will be generated.
 */

t_restp *get_restp(const char *rtpname, int nrtp, t_restp rtp[]);
/* Return the entry in the rtp database with rtp name rtpname.
 * Generates a fatal error when rtpname is not found.
 */

gpp_atomtype_t read_atype(const char *ffdir, struct t_symtab *tab);
/* read atom type database(s) */

void read_resall(char *resdb, int *nrtp, t_restp **rtp,
                 gpp_atomtype_t atype, struct t_symtab *tab,
                 gmx_bool bAllowOverrideRTP);
/* read rtp database, append to the existing database */

void print_resall(FILE *out, int nrtp, t_restp rtp[],
                  gpp_atomtype_t atype);
/* write rtp database */
#ifdef __cplusplus
}
#endif

#endif
