/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _h_db_h
#define _h_db_h
#include "visibility.h"
#include "sysstuff.h"
#include "hackblock.h"

/* number of control atoms for each 'add' type */
extern const int ncontrol[];

/* functions for the h-database */

extern void read_ab(char *line, const char *fn, t_hack *ab);
/* Read one add block */

GMX_LIBGMXPREPROCESS_EXPORT
extern int read_h_db(const char *ffdir, t_hackblock **ah);
/* Read the database from hdb file(s) in ffdir or current dir */

extern void print_ab(FILE *out, t_hack *ab, char *nname);
/* print one add block */

extern void print_h_db(FILE *out, int nh, t_hackblock ah[]);
/* Print the database to file */

extern int compaddh(const void *a, const void *b);

extern t_hackblock *search_h_db(int nh, t_hackblock ah[], char *key);
/* Search for an entry in the database */



#endif  /* _h_db_h */
