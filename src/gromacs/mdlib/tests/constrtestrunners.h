/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
 * \brief SHAKE and LINCS tests header.
 *
 * Contains description and constructor for the test data accumulating object,
 * declares CPU- and GPU-based functions used to apply SHAKE or LINCS on the
 * test data.
 *
 * \author Artem Zhmurov <zhmurov@gmail.com>
 * \ingroup module_mdlib
 */

#ifndef GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H
#define GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H

#include "constrtestdata.h"

struct t_pbc;

namespace gmx
{
namespace test
{

/*! \brief Apply SHAKE constraints to the test data.
 */
void applyShake(ConstraintsTestData* testData, t_pbc pbc);
/*! \brief Apply LINCS constraints to the test data.
 */
void applyLincs(ConstraintsTestData* testData, t_pbc pbc);
/*! \brief Apply CUDA version of LINCS constraints to the test data.
 *
 * All the data is copied to the GPU device, then LINCS is applied and
 * the resulting coordinates are copied back.
 */
void applyLincsCuda(ConstraintsTestData* testData, t_pbc pbc);

} // namespace test
} // namespace gmx

#endif // GMX_MDLIB_TESTS_CONSTRTESTRUNNERS_H
