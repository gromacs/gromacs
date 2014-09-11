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
#ifndef MOLPROP_TABLES_H
#define MOLPROP_TABLES_H

#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/topology/atomprop.h"
#include "poldata.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molselect.h"

//! Utility function converting float to char *
extern char *ftoa(double f);

//! Utility function converting int to char *
extern char *itoa(int f);

/*! \brief
 * Generates a LaTeX table containing the statistics (RMSD from experiment) of a calculated property per molecule category
 *
 * \param[out] fp   File pointer to write to
 * \param[in] eprop The property of choice
 * \param[in] mp    Array of molecules
 * \param[in] qmc
 * \param[in] iQM
 * \param[in] lot   Level of theory in 'A/B' form
 * \param[in] exp_type The type of experimental property
 * \param[in] outlier  Outlier indicates a level (in units of sigma, one standard deviation). Calculations that deviate more than this level from the experiment are not taken into account when computing statistics. Moreover, the outliers are printed to the standard error. If outlier is 0, no action is taken.
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform iQM to enum
 * \ingroup module_alexandria
 */
extern void gmx_molprop_stats_table(FILE *fp, MolPropObservable eprop,
                                    std::vector<alexandria::MolProp> mp,
                                    t_qmcount *qmc, char *exp_type,
                                    double outlier, gmx_molselect_t gms, iMolSelect ims);

/*! \brief
 * Generates a LaTeX table containing the composition in atoms of each molecule
 *
 * \param[out] fp   File pointer to write to
 * \param[in] mp    Array of molecules
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \ingroup module_alexandria
 */
extern void gmx_molprop_composition_table(FILE *fp,
                                          std::vector<alexandria::MolProp> mp,
                                          gmx_molselect_t gms, iMolSelect ims);

/*! \brief
 * Generates a LaTeX table containing the molecules in each molecule category
 *
 * \param[out] fp   File pointer to write to
 * \param[in] mp    Array of molecules
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \ingroup module_alexandria
 */
extern void gmx_molprop_category_table(FILE *fp,
                                       std::vector<alexandria::MolProp> mp,
                                       gmx_molselect_t gms, iMolSelect ims);

/*! \brief
 * Generates a LaTeX table containing properties for molecules from different sources
 *
 * \param[out] fp   File pointer to write to
 * \param[in] eprop The property of choice
 * \param[in] rel_toler Relative tolerance in the property to print it bold (i.e. indicating an outlier)
 * \param[in] abs_toler If non-zero, takes prevalence over rel_toler, and indicates the absolute tolerance for this property for designating as an outlier
 * \param[in] mp    Array of molecules
 * \param[in] qmc   Statistics of quantum calculations
 * \param[in] bPrintAll If set print also properties of molecules for which no experimental data is available
 * \param[in] bPrintBasis Print the basis set in the table header
 * \param[in] bPrintMultQ Print the multiplicity and total charge of the molecule
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \todo Introduce enum to switch between absolute and relative tolerance
 * \ingroup module_alexandria
 */
extern void gmx_molprop_prop_table(FILE *fp, MolPropObservable eprop,
                                   real rel_toler, real abs_toler,
                                   std::vector<alexandria::MolProp> mp,
                                   t_qmcount *qmc,
                                   bool bPrintAll, bool bPrintBasis,
                                   bool bPrintMultQ, gmx_molselect_t gms, iMolSelect ims);

/*! \brief
 * Generates a LaTeX table containing the atomtypes in Alexandria
 *
 * \param[out] fp   File pointer to write to
 * \param[in] bPolar Calling this with bPolar TRUE will print an atomtype table with polarizability information. With bPolar FALSE it will print the same table with EEM parameters.
 * \param[in] pd    Force field file
 * \param[in] mp    Array of molecules
 * \param[in] iQM   Selecting the reference between QM and Experimental
 * \param[in] lot   Level of theory in 'A/B' format
 * \param[in] exp_type The type of experimental property
 * \param[in] exp_type The type of experimental property
 * \param[in] oenv  Information for generating xvg files
 * \param[in] histo File name for histogram data
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform iQM to enum
 * \todo More explanation text
 * \ingroup module_alexandria
 */
extern void gmx_molprop_atomtype_table(FILE *fp, bool bPolar,
                                       gmx_poldata_t pd,
                                       std::vector<alexandria::MolProp> mp,
                                       char *lot, char *exp_type);

#endif
