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
#ifndef MOLPROP_UTIL_H
#define MOLPROP_UTIL_H

#include "gromacs/gmxpreprocess/grompp-impl.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/topology/atomprop.h"
#include "molprop.h"
#include "poldata.h"
#include "molselect.h"

extern void generate_composition(std::vector<alexandria::MolProp> &mp, gmx_poldata_t pd);
extern void generate_formula(std::vector<alexandria::MolProp> &mp, gmx_atomprop_t ap);

alexandria::MolProp atoms_2_molprop(char *molname, int natoms, char **smnames,
                                    gmx_atomprop_t ap);

/* Return number of atoms, 0 means failure */
extern bool molprop_2_topology2(alexandria::MolProp mp, gmx_atomprop_t ap,
                                t_symtab *tab, const char *lot,
                                t_topology **top, const char *q_algorithm,
                                rvec **x, t_params plist[F_NRE],
                                int nexcl, t_excls **excls);

extern void merge_doubles(std::vector<alexandria::MolProp> &mp,
                          char *doubles, bool bForceMerge);

extern void merge_xml(int nfile, char **infiles,
                      std::vector<alexandria::MolProp> &mp,
                      char *outf, char *sorted, char *doubles,
                      gmx_atomprop_t ap, gmx_poldata_t pd,
                      bool bForceMerge);

typedef struct {
    int    n;
    int   *count;
    char **method, **basis, **type, **lot;
    int    nconf;
    char **conf;
} t_qmcount;

/* Check the available molprops to see what kind of calculations are stored in there */
extern t_qmcount *find_calculations(std::vector<alexandria::MolProp> mp,
                                    MolPropObservable mpo, const char *fc_str);

/*! \brief
 * Enumerated type for MolPropSort function
 *
 * \inpublicapi
 * \ingroup module_alexandria
 */
enum MolPropSortAlgorithm {
    MPSA_MOLNAME,
    MPSA_FORMULA,
    MPSA_COMPOSITION,
    MPSA_SELECTION,
    MPSA_NR
};

/*! \brief
 * Sorts a vector of molprops
 *
 * Function that uses the std::sort routine and can apply different sorting
 * keys.
 *
 * \param[inout]  mp        The vector of MolProp
 * \param[in]     mpsa      The algorithm used for sorting
 * \param[in]     apt       Database of atom properties
 * \param[in]     mgs       Optional structure containing selection criteria
 * \ingroup module_alexandria
 */
extern void MolPropSort(std::vector<alexandria::MolProp> &mp,
                        MolPropSortAlgorithm mpsa, gmx_atomprop_t apt,
                        gmx_molselect_t gms);

#endif
