/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_TABLES_H
#define MOLPROP_TABLES_H

#include "categories.h"
#include "molprop.h"
#include "molprop_util.h"

//! Utility function converting float to char *
extern char *ftoa(double f);

//! Utility function converting int to char *
extern char *itoa(int f);

namespace alexandria
{
/*! \brief
 * Generates a LaTeX table containing the statistics (RMSD from experiment) of a calculated property per molecule category
 *
 * \param[out] fp   File pointer to write to
 * \param[in] eprop The property of choice
 * \param[in] mp    Array of molecules
 * \param[in] qmc   QM calc statistics
 * \param[in] iQM
 * \param[in] lot   Level of theory in 'A/B' form
 * \param[in] exp_type The type of experimental property
 * \param[in] outlier  Outlier indicates a level (in units of sigma, one standard deviation). Calculations that deviate more than this level from the experiment are not taken into account when computing statistics. Moreover, the outliers are printed to the standard error. If outlier is 0, no action is taken.
 * \param[in] catmin Mininum number of molecules in a category for it to be printed
 * \param[in] cList Structure containing all categories and the number of molecules per category
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform iQM to enum
 * \ingroup module_alexandria
 */
void gmx_molprop_stats_table(FILE                 *fp,
                             MolPropObservable     eprop,
                             std::vector<MolProp> &mp,
                             const QmCount        &qmc,
                             char                 *exp_type,
                             double                outlier,
                             CategoryList          cList,
                             const MolSelect      &gms,
                             iMolSelect            ims);

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
void gmx_molprop_composition_table(FILE                 *fp,
                                   std::vector<MolProp>  mp,
                                   const MolSelect      &gms,
                                   iMolSelect            ims);

/*! \brief
 * Generates a LaTeX table containing the molecules in each molecule category
 *
 * \param[out] fp    File pointer to write to
 * \param[in] catmin Mininum number of molecules in a category for it to be printed
 * \param[in] cList  Structure containing all categories and the number of molecules per category
 * \ingroup module_alexandria
 */
void gmx_molprop_category_table(FILE                    *fp,
                                int                      catmin,
                                alexandria::CategoryList cList);

/*! \brief
 * Generates a LaTeX table containing properties for molecules from different sources
 *
 * \param[out] fp   File pointer to write to
 * \param[in] eprop The property of choice
 * \param[in] rel_toler Relative tolerance in the property to print it bold (i.e. indicating an outlier)
 * \param[in] abs_toler If non-zero, takes prevalence over rel_toler, and indicates the absolute tolerance for this property for designating as an outlier
 * \param[in] mp    Array of molecules
 * \param[in] qmc   Statistics of quantum calculations
 * \param[in] exp_type Which type of this property should we print?
 * \param[in] bPrintAll If set print also properties of molecules for which no experimental data is available
 * \param[in] bPrintBasis Print the basis set in the table header
 * \param[in] bPrintMultQ Print the multiplicity and total charge of the molecule
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \todo Introduce enum to switch between absolute and relative tolerance
 * \ingroup module_alexandria
 */
void gmx_molprop_prop_table(FILE                             *fp,
                            MolPropObservable                 eprop,
                            real                              rel_toler,
                            real                              abs_toler,
                            std::vector<alexandria::MolProp> &mp,
                            const QmCount                    &qmc,
                            const char                       *exp_type,
                            bool                              bPrintAll,
                            bool                              bPrintBasis,
                            bool                              bPrintMultQ,
                            const MolSelect                  &gms,
                            iMolSelect                        ims);

/*! \brief
 * Generates a LaTeX table containing the atomtypes in Alexandria
 *
 * \param[out] fp   File pointer to write to
 * \param[in] bPolar Calling this with bPolar TRUE will print an atomtype table with polarizability information. With bPolar FALSE it will print the same table with EEM parameters.
 * \param[in] pd    Force field data
 * \param[in] mp    Array of molecules
 * \param[in] lot   Level of theory in 'A/B' format
 * \param[in] exp_type The type of experimental property
 * \todo More explanation text
 * \ingroup module_alexandria
 */
void gmx_molprop_atomtype_table(FILE                                   *fp,
                                bool                                    bPolar,
                                const Poldata                          &pd,
                                const std::vector<alexandria::MolProp> &mp,
                                const char                             *lot,
                                const char                             *exp_type);

} // namespace alexandria

#endif
