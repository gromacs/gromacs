/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012, by the GROMACS development team, led by
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

#ifndef _topio_h
#define _topio_h
#include "visibility.h"
#include "typedefs.h"
#include "readir.h"
#include "grompp.h"
#include "gpp_atomtype.h"

GMX_LIBGMXPREPROCESS_EXPORT
extern double check_mol(gmx_mtop_t *mtop,warninp_t wi);
/* Check mass and charge */

GMX_LIBGMXPREPROCESS_EXPORT
extern char **do_top(gmx_bool         bVerbose,
		     const char   *topfile,
		     const char   *topppfile,
		     t_gromppopts *opts,
		     gmx_bool         bZero,
		     t_symtab     *symtab,
		     t_params     plist[],
		     int          *combination_rule,
		     double       *repulsion_power,
		     real         *fudgeQQ,
		     gpp_atomtype_t atype,
		     int          *nrmols,
		     t_molinfo    **molinfo,
		     t_inputrec   *ir,
		     int          *nmolblock,
		     gmx_molblock_t **molblock,
		     gmx_bool         bGB,
		     warninp_t    wi);


/* This routine expects sys->molt[m].ilist to be of size F_NRE and ordered. */
GMX_LIBGMXPREPROCESS_EXPORT
void generate_qmexcl(gmx_mtop_t *sys,t_inputrec *ir,warninp_t    wi);

#endif	/* _topio_h */
