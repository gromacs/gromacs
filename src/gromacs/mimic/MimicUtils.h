/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief Provides utility functions for MiMiC QM/MM
 *
 * \author Viacheslav Bolnykh <v.bolnykh@hpc-leap.eu>
 * \ingroup module_utility
 */
#ifndef GMX_MIMIC_MIMICUTILS_H
#define GMX_MIMIC_MIMICUTILS_H

#include <vector>

/*! \brief Generates the list of QM atoms
 *
 * This function generates vector containing global IDs of QM atoms
 * based on the information stored in the QM/MM group (egcQMMM)
 *
 * \param[in]    mtop   Global topology object
 * \return              The list of global IDs of QM atoms
 */
std::vector<int> genQmmmIndices(const gmx_mtop_t &mtop);

/*! \brief Updates inter-molecular exclusion lists
 *
 * This function updates inter-molecular exclusions to exclude all
 * non-bonded interactions between QM atoms
 *
 * \param[inout]    excls   existing exclusions in local topology
 * \param[in]       ids     list of global QM atoms IDs
 */
void addExclusions(t_blocka               *excls,
                   const std::vector<int> &ids);

#endif //GMX_MIMIC_MIMICUTILS_H
