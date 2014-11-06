/*
 * This source file is part of the Alexandria project.
 *
 * Copyright (C) 2014 David van der Spoel and Paul J. van Maaren
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GENTOP_NM2TYPE_H
#define GENTOP_NM2TYPE_H

#include <stdio.h>
#include "gromacs/topology/atomprop.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "poldata.h"
#include "gentop_vsite.h"

extern int nm2type(FILE *fp, const char *molname,
                   gmx_poldata_t pd, gmx_atomprop_t aps,
                   t_symtab *tab, t_atoms *atoms, gmx_bool bRing[],
                   double bondorder[],
                   gpp_atomtype_t atype, int *nbonds, t_params *bond,
                   char **smname,
                   rvec x[], t_pbc *pbc, real th_toler, real phi_toler,
                   alexandria::GentopVsites gvt);
/* Try to determine the atomtype (force field dependent) for the atoms
 * with help of the bond list and the coordinates!
 */

gpp_atomtype_t set_atom_type(FILE *fp, const char *molname,
                             t_symtab *tab, t_atoms *atoms, t_params *bonds,
                             int nbonds[], gmx_bool bRing[], double bondorder[],
                             gmx_poldata_t pd, gmx_atomprop_t aps,
                             rvec x[], t_pbc *pbc, real th_toler,
                             real ph_toler, alexandria::GentopVsites gvt);


#endif
