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
/*! \libinternal\file
 * \brief Declares utilities for initializing long range correction tables
 *
 * \author Berk Hess <hess@kth.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \inlibraryapi
 */
#ifndef GMX_TABLES_LONGRANGETABLES_H
#define GMX_TABLES_LONGRANGETABLES_H

#include <vector>

/*! \brief
 * Compute the long range energy and force for Ewald Coulomb
 *
 * \param[in]  params Parameters to the long range function
 * \param[in]  r      Distance
 * \param[out] e      Energy
 * \param[out] f      Force
 */
void v_q_ewald_lr(const std::vector<double> &params,
                  double r, double *e, double *f);

/*! \brief
 * Compute the long range energy and force for Ewald Lennard Jones (r^{-6} only)
 *
 * \param[in]  params Parameters to the long range function
 * \param[in]  r      Distance
 * \param[out] e      Energy
 * \param[out] f      Force
 */
void v_lj_ewald_lr(const std::vector<double> &params,
                   double r, double *e, double *f);

/*! \brief Function to compute table scaling for long range corrections
 *
 * Distances are multiplied with the table scaling to obtain the index
 * in the tables when interpolating.
 * \param[in] type_coul       The electrostatics type
 * \param[in] ewaldcoeff_coul The coefficient for Coulomb
 * \param[in] r_coul          The Coulomb cut-off (switching distance)
 * \param[in] type_vdw        The Van der Waals type
 * \param[in] ewaldcoeff_vdw  The coefficient for Van der Waals
 * \param[in] r_vdw           The Van der Waals cut-off (switching distance)
 * \return The table scaling
 */
double longRangeCorrectionTableScale(int    type_coul,
                                     double ewaldcoeff_coul,
                                     double r_coul,
                                     int    type_vdw,
                                     double ewaldcoeff_vdw,
                                     double r_vdw);
/*! \brief Compute number of data points in the table
 *
 * \param[in] rtab  The length of the table
 * \param[in] scale The table scaling
 * \return The number of points needed.
 */
unsigned int longRangeCorrectionTableSize(double rtab,
                                          double scale);

#endif /* GMX_TABLES_TABLES_H */
