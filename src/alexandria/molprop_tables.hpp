/*
 * $Id: molprop_tables.h,v 1.10 2009/05/29 15:01:18 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */

#ifndef _molprop_tables_hpp
#define _molprop_tables_hpp

#include "grompp.h"
#include "atomprop.h"
#include "poldata.h"
#include "molprop.hpp"
#include "molprop_util.hpp"
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
 * \param[in] ntot  
 * \param[in] qmc
 * \param[in] iQM
 * \param[in] lot   Level of theory in 'A/B' form
 * \param[in] outlier  Outlier indicates a level (in units of sigma, one standard deviation). Calculations that deviate more than this level from the experiment are not taken into account when computing statistics. Moreover, the outliers are printed to the standard error. If outlier is 0, no action is taken.
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \todo Transform iQM to enum
 * \ingroup module_alexandria
 */
extern void gmx_molprop_stats_table(FILE *fp,MolPropObservable eprop,
                                    std::vector<alexandria::MolProp> mp,
                                    int ntot,
                                    t_qmcount *qmc,int iQM,char *lot,
                                    double outlier,gmx_molselect_t gms,int ims);

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
                                          gmx_molselect_t gms,int ims);

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
                                       gmx_molselect_t gms,int ims);

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
 * \todo Replace gmx_bool by C++ bool
 * \ingroup module_alexandria
 */
extern void gmx_molprop_prop_table(FILE *fp,MolPropObservable eprop,
                                   real rel_toler,real abs_toler,
                                   std::vector<alexandria::MolProp> mp,
                                   t_qmcount *qmc,
                                   gmx_bool bPrintAll,gmx_bool bPrintBasis,
                                   gmx_bool bPrintMultQ,gmx_molselect_t gms,int ims);
				   
/*! \brief
 * Generates a LaTeX table containing the atomtypes in Alexandria
 *
 * \param[out] fp   File pointer to write to
 * \param[in] bPolar Calling this with bPolar TRUE will print an atomtype table with polarizability information. With bPolar FALSE it will print the same table with EEM parameters.
 * \param[in] npd   Number of force field structurs
 * \param[in] pd    Array of force field structures
 * \param[out] pd_aver Force field structure with average and standard deviation of parameters
 * \param[in] mp    Array of molecules
 * \param[in] iQM   Selecting the reference between QM and Experimental
 * \param[in] lot   Level of theory in 'A/B' format
 * \param[in] oenv  Information for generating xvg files
 * \param[in] histo File name for histogram data
 * \param[in] gms   Structure containing selections of which molecules to output
 * \param[in] ims   The actual selection of the right set
 * \todo Transform ims to enum
 * \todo Transform iQM to enum
 * \todo More explanation text
 * \ingroup module_alexandria
 */
extern void gmx_molprop_atomtype_table(FILE *fp,gmx_bool bPolar,
                                       int npd,gmx_poldata_t pd[],
                                       gmx_poldata_t pd_aver, /* Output! */
                                       std::vector<alexandria::MolProp> mp,
                                       int iQM,char *lot,output_env_t oenv,
                                       const char *histo);

#endif
