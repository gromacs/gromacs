/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015, by the GROMACS development team, led by
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
#ifndef GMX_TABLES_TABLES_H
#define GMX_TABLES_TABLES_H

#include "gromacs/legacyheaders/types/interaction_const.h"
#include "gromacs/legacyheaders/types/simple.h"

/*! \libinternal \brief
 * Function pointer to fill table.
 *
 * Function pointer used to tell table_spline3_fill_ewald_lr whether it
 * should calculate the grid contribution for electrostatics or LJ.
 * The first parameters is an array of constant used to evaluate
 * the function, the second parameter is the distance.
 */
typedef double (interaction_potential_function)(const double *, double);

/*! \libinternal\brief
 * Fill a table using a spline interpolation.
 *
 * This function interpolates the Ewald mesh potential contribution
 * with coefficient beta using a quadratic spline.
 * The force can then be interpolated linearly.
 * \param[out] table_F The array containing the force table
 * \param[out] table_V
 * \param[out] table_FDV0
 * \param[in]  ntab       Length of the table
 * \param[in]  dx         Spacing between the points
 * \param[in]  beta       The temperature factor
 * \param[in]  v_lr       A function pointer that can be used as an
 *                        alternative to the default Ewald long range.
 */
void table_spline3_fill_ewald_lr(real                           *table_F,
                                 real                           *table_V,
                                 real                           *table_FDV0,
                                 int                             ntab,
                                 double                          dx,
                                 real                            beta,
                                 interaction_potential_function  v_lr);

/*! \libinternal \brief
 * Return the scaling for the Ewald quadratic spline tables.
 *
 * \param[in] ic
 * \return the scaling factor.
 */
real ewald_spline3_table_scale(const interaction_const_t *ic);

/*! \libinternal\brief
 * Return the long range energy for Ewald Coulomb
 *
 * \param[in] beta Parameters to the long range function
 * \param[in] r    Distance
 * \return Long range energy
 */
double v_q_ewald_lr(const double *beta, double r);

/*! \libinternal\brief
 * Return the long range energy for Ewald Lennard Jones (r^{-6} only)
 *
 * \param[in] beta Parameters to the long range function
 * \param[in] r    Distance
 * \return Long range energy
 */
double v_lj_ewald_lr(const double *beta, double r);

#endif  /* GMX_TABLES_TABLES_H */
