/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements functionality from any.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/any.h"

#include <cstdint>

#include <string>

#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{

//! \cond libapi
std::string simpleValueToString(const Any& value)
{
    if (value.isType<bool>())
    {
        // TODO: Consider if this would be better as yes/no instead of
        // true/false.
        return toString(value.cast<bool>());
    }
    else if (value.isType<float>())
    {
        return toString(value.cast<float>());
    }
    else if (value.isType<double>())
    {
        return toString(value.cast<double>());
    }
    else if (value.isType<int>())
    {
        return toString(value.cast<int>());
    }
    else if (value.isType<int64_t>())
    {
        return toString(value.cast<int64_t>());
    }
    else if (value.isType<std::string>())
    {
        return value.cast<std::string>();
    }
    GMX_RELEASE_ASSERT(false, "Unknown value type");
    return std::string();
}
//! \endcond

} // namespace gmx
