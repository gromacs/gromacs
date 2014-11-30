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
#ifndef GENTOP_CORE_H
#define GENTOP_CORE_H

#include <stdio.h>
#include <vector>
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/gmxpreprocess/gpp_nextnb.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "plistwrapper.h"
#include "poldata.h"

void calc_angles_dihs(t_params *ang, t_params *dih, rvec x[], gmx_bool bPBC, matrix box);

real calc_dip(t_atoms *atoms, rvec x[]);

void dump_hybridization(FILE *fp, t_atoms *atoms, int nbonds[]);

void reset_q(t_atoms *atoms);

void print_rtp(const char *filenm, const char *title, t_atoms *atoms,
               t_params plist[], int cgnr[], int nbts, int bts[]);

void symmetrize_charges(gmx_bool bQsym,
                        t_atoms *atoms, 
                        alexandria::PlistWrapperIterator bonds,
                        gmx_poldata_t pd,
                        gmx_atomprop_t aps, const char *symm_string,
                        std::vector<int> &sym_charges);

enum eChargeGroup {
    ecgAtom, ecgGroup, ecgNeutral, ecgNR
};

int *generate_charge_groups(eChargeGroup cgtp, t_atoms *atoms,
                            std::vector<alexandria::PlistWrapper> &pw,
                            bool bUsePDBcharge,
                            real *qtot, real *mtot);

void sort_on_charge_groups(int *cgnr, t_atoms *atoms, 
                           std::vector<alexandria::PlistWrapper> &pw,
                           rvec x[], t_excls excls[],
                           const char *ndxout,
                           int nmol);

#endif
