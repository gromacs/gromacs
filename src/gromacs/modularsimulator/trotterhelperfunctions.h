/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 * \brief Defines helper functions used by the Trotter decomposition
 * algorithms (Nose-Hoover chains, MTTK)
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */
#ifndef GMX_MODULARSIMULATOR_TROTTERHELPERFUNCTIONS_H
#define GMX_MODULARSIMULATOR_TROTTERHELPERFUNCTIONS_H

#include "modularsimulatorinterfaces.h"

namespace gmx
{

/*! \brief Check whether two times are nearly equal
 *
 * Times are considered close if their absolute difference is smaller
 * than c_timePrecision.
 *
 * \param time           The test time
 * \param referenceTime  The reference time
 * \return bool          Whether the absolute difference is < c_timePrecision
 */
inline bool timesClose(Time time, Time referenceTime)
{
    /* Expected time precision
     * Times are typically incremented in the order of 1e-3 ps (1 fs), so
     * 1e-6 should be sufficiently tight.
     */
    constexpr real c_timePrecision = 1e-6;

    return (time - referenceTime) * (time - referenceTime) < c_timePrecision * c_timePrecision;
}

} // namespace gmx
#endif // GMX_MODULARSIMULATOR_TROTTERHELPERFUNCTIONS_H
