/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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
 * \brief
 * Implements functions in strconvert.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "strconvert.h"

#include <cerrno>
#include <cstdlib>

#include <limits>
#include <string>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

//! \cond libapi

bool boolFromString(const char* value)
{
    if (gmx_strcasecmp(value, "1") == 0 || gmx_strcasecmp(value, "yes") == 0
        || gmx_strcasecmp(value, "true") == 0)
    {
        return true;
    }
    if (gmx_strcasecmp(value, "0") == 0 || gmx_strcasecmp(value, "no") == 0
        || gmx_strcasecmp(value, "false") == 0)
    {
        return false;
    }
    GMX_THROW(InvalidInputError("Invalid value: '" + std::string(value)
                                + "'; supported values are: 1, 0, yes, no, true, false"));
}

int intFromString(const char* str)
{
    errno = 0;
    char*          endptr;
    const long int value = std::strtol(str, &endptr, 10);
    if (errno == ERANGE || value < std::numeric_limits<int>::min()
        || value > std::numeric_limits<int>::max())
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str)
                                    + "'; it causes an integer overflow"));
    }
    if (str[0] == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str) + "'; expected an integer"));
    }
    return value;
}

int64_t int64FromString(const char* str)
{
    errno = 0;
    char*         endptr;
    const int64_t value = str_to_int64_t(str, &endptr);
    if (errno == ERANGE)
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str)
                                    + "'; it causes an integer overflow"));
    }
    if (str[0] == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str) + "'; expected an integer"));
    }
    return value;
}

float floatFromString(const char* str)
{
    errno = 0;
    char*        endptr;
    const double value = std::strtod(str, &endptr);
    if (errno == ERANGE || value < -std::numeric_limits<float>::max()
        || value > std::numeric_limits<float>::max())
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str)
                                    + "'; it causes an overflow/underflow"));
    }
    if (str[0] == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str) + "'; expected a number"));
    }
    return value;
}

double doubleFromString(const char* str)
{
    errno = 0;
    char*        endptr;
    const double value = std::strtod(str, &endptr);
    if (errno == ERANGE)
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str)
                                    + "'; it causes an overflow/underflow"));
    }
    if (str[0] == '\0' || *endptr != '\0')
    {
        GMX_THROW(InvalidInputError("Invalid value: '" + std::string(str) + "'; expected a number"));
    }
    return value;
}

//! \endcond

} // namespace gmx
