/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \libinternal \file
 * \brief
 * Declares structures for analytical or numerical input data to construct tables
 *
 * \inlibraryapi
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \ingroup module_tables
 */

#ifndef GMX_TABLES_TABLEINPUT_H
#define GMX_TABLES_TABLEINPUT_H

#include <functional>
#include <vector>

#include "gromacs/utility/arrayref.h"

namespace gmx
{

/*! \libinternal \brief Specification for analytical table function (name, function, derivative)
 */
struct AnalyticalSplineTableInput
{
    //NOLINTNEXTLINE(google-runtime-member-string-references)
    const std::string&            desc;       //!< \libinternal Brief description of function
    std::function<double(double)> function;   //!< \libinternal Analytical form of function
    std::function<double(double)> derivative; //!< \libinternal Analytical derivative
};

/*! \libinternal \brief Specification for vector table function (name, function, derivative, spacing)
 */
struct NumericalSplineTableInput
{
    //NOLINTNEXTLINE(google-runtime-member-string-references)
    const std::string&     desc;       //!< \libinternal Brief description of function
    ArrayRef<const double> function;   //!< \libinternal Vector with function values
    ArrayRef<const double> derivative; //!< \libinternal Vector with derivative values
    double                 spacing;    //!< \libinternal Distance between data points
};


} // namespace gmx


#endif // GMX_TABLES_TABLEINPUT_H
