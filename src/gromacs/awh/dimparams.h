/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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
 *
 * \brief
 * Declares the DimParams struct and AWH vector types.
 *
 * This class holds the physical information for a dimension
 * of the bias reaction-coordinate grid.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#ifndef GMX_AWH_DIMPARAMS_H
#define GMX_AWH_DIMPARAMS_H

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"

namespace gmx
{

//! The maximum dimensionality of the AWH coordinate.
static const int c_biasMaxNumDim = 4;

//! A real vector in AWH coordinate space.
typedef double awh_dvec[c_biasMaxNumDim];

//! An integer vector in AWH coordinate space.
typedef int awh_ivec[c_biasMaxNumDim];

/*! \internal \brief Constant parameters for each dimension of the coordinate.
 */
struct DimParams
{
    /*! \brief
     * Constructor.
     *
     * \param[in] conversionFactor  Conversion factor from user coordinate units to bias internal units (=DEG2RAD for angles).
     * \param[in] forceConstant     The harmonic force constant.
     * \param[in] beta              1/(k_B T).
     */
    DimParams(double conversionFactor,
              double forceConstant,
              double beta) :
        k(forceConstant),
        betak(beta*forceConstant),
        userCoordUnitsToInternal(conversionFactor)
    {
    };

    /*! \brief Convert internal coordinate units to external, user coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleInternalToUserInput(double value) const
    {
        return value/userCoordUnitsToInternal;
    }

    /*! \brief Convert external, user coordinate units to internal coordinate units.
     *
     * \param[in] value               Value to convert.
     * \returns the converted value.
     */
    double scaleUserInputToInternal(double value) const
    {
        return value*userCoordUnitsToInternal;
    }

    const double k;                        /**< Force constant (kJ/mol/nm^2) for each coordinate dimension. */
    const double betak;                    /**< Inverse variance (1/nm^2) for each coordinate dimension. */
    const double userCoordUnitsToInternal; /**< Conversion factor coordinate units. */
};

}      // namespace gmx

#endif /* GMX_AWH_DIMPARAMS_H */
