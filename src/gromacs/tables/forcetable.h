/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_TABLES_FORCETABLE_H
#define GMX_TABLES_FORCETABLE_H

/*! \libinternal \file
 * \brief
 * Old routines for table generation (will eventually be replaced)
 *
 * \inlibraryapi
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */

#include <cstdio>

#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/utility/real.h"

/*! \brief Flag to select user tables for make_tables */
#define GMX_MAKETABLES_FORCEUSER  (1<<0)
/*! \brief Flag to only make 1,4 pair tables for make_tables */
#define GMX_MAKETABLES_14ONLY     (1<<1)

/*! \brief Enumerated type to describe the interaction types in a table */
enum {
    etiCOUL,  //!< Coulomb
    etiLJ6,   //!< Dispersion
    etiLJ12,  //!< Repulsion
    etiNR     //!< Total number of interaction types
};

/*! \brief Function pointer to calculate the grid contribution for coulomb/LJ
 *
 * Used to tell table_spline3_fill_ewald_lr whether it
 * should calculate the grid contribution for electrostatics or LJ.
 */
typedef double (*real_space_grid_contribution_computer)(double, double);


/*! \brief Fill tables with the Ewald long-range force interaction
 *
 * Fill tables of ntab points with spacing dr with the ewald long-range
 * (mesh) force.
 * There are three separate tables with format FDV0, F, and V.
 * This function interpolates the Ewald mesh potential contribution
 * with coefficient beta using a quadratic spline.
 * The force can then be interpolated linearly.
 *
 * \param table_F    Force table
 * \param table_V    Potential table
 * \param table_FDV0 Combined table optimized for SIMD loads
 * \param ntab       Number of points in tables
 * \param dx         Spacing
 * \param beta       Ewald splitting paramter
 * \param v_lr       Pointer to function calculating real-space grid contribution
 */
void table_spline3_fill_ewald_lr(real                                 *table_F,
                                 real                                 *table_V,
                                 real                                 *table_FDV0,
                                 int                                   ntab,
                                 double                                dx,
                                 real                                  beta,
                                 real_space_grid_contribution_computer v_lr);

/*! \brief Compute scaling for the Ewald quadratic spline tables.
 *
 * \param ic  Pointer to interaction constant structure
 * \return The scaling factor
 */
real ewald_spline3_table_scale(const interaction_const_t *ic);

/*! \brief Return the real space grid contribution for Ewald
 *
 *  \param beta  Ewald splitting parameter
 *  \param r     Distance for which to calculate the real-space contrib
 *  \return      Real space grid contribution for Ewald electrostatics
 */
double v_q_ewald_lr(double beta, double r);

/*! \brief Return the real space grid contribution for LJ-Ewald
 *
 *  \param beta  Ewald splitting parameter
 *  \param r     Distance for which to calculate the real-space contrib
 *  \return      Real space grid contribution for Ewald Lennard-Jones interaction
 */
double v_lj_ewald_lr(double beta, double r);

/*! \brief Return tables for inner loops.
 *
 * \param fp     Log file pointer
 * \param ic     Non-bonded interaction constants
 * \param fn     File name from which to read user tables
 * \param rtab   Largest interaction distance to tabulate
 * \param flags  Flags to select table settings
 *
 * \return Pointer to inner loop table structure
 */
t_forcetable *make_tables(FILE *fp,
                          const interaction_const_t *ic,
                          const char *fn, real rtab, int flags);

/*! \brief Return a table for bonded interactions,
 *
 * \param  fplog   Pointer to log file
 * \param  fn      File name
 * \param  angle   Type of angle: bonds 0, angles 1, dihedrals 2
 * \return New bonded table datatype
 */
bondedtable_t make_bonded_table(FILE *fplog, const char *fn, int angle);

/*! \brief Return a table for GB calculations
 *
 * \param fr   Force record
 * \return     Pointer to new gb table structure
 */
t_forcetable *make_gb_table(const t_forcerec              *fr);

/*! \brief Construct and return tabulated dispersion and repulsion interactions
 *
 * This table can be used to compute long-range dispersion corrections */
t_forcetable *makeDispersionCorrectionTable(FILE *fp, const interaction_const_t *ic,
                                            real rtab, const char *tabfn);

#endif  /* GMX_TABLES_FORCETABLE_H */
