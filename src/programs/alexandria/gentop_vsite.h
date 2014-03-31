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
#ifndef GENTOP_VSITE_H
#define GENTOP_VSITE_H

#include "gromacs/gmxpreprocess/grompp.h"
#include "poldata.h"

enum {
    egvtNO, egvtLINEAR, egvtPLANAR, egvtRING_PLANAR, egvtALL, egvtNR
};

typedef struct gentop_vsite *gentop_vsite_t;

extern gentop_vsite_t gentop_vsite_init(int egvt);

extern void gentop_vsite_done(gentop_vsite_t *gvt);

extern void gentop_vsite_add_linear(gentop_vsite_t gvt, int a1, int a2, int a3);

extern void gentop_vsite_add_planar(gentop_vsite_t gvt, int a1, int a2, int a3, int a4, int nbonds[]);

extern void gentop_vsite_add_ring_planar(gentop_vsite_t gvt, int natom,
                                         int a[], int nbonds[]);

extern void gentop_vsite_generate_special(gentop_vsite_t gvt, gmx_bool bGenVsites,
                                          t_atoms *atoms, rvec **x,
                                          t_params plist[],
                                          t_symtab *symtab, gpp_atomtype_t atype,
                                          t_excls **excls, gmx_poldata_t pd);

#endif
