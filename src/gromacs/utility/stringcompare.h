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
 * Declares utility functionst for string comparison.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRINGCOMPARE_H
#define GMX_UTILITY_STRINGCOMPARE_H

#include <string>

#include "gromacs/utility/cstringutil.h"

namespace gmx
{

//! \cond libapi
/*! \brief
 * Specifies how strings should be compared in various contexts.
 *
 * \ingroup module_utility
 */
enum class StringCompareType
{
    //! Only exact matches are accepted.
    Exact,
    //! Case-insensitive comparison.
    CaseInsensitive,
    //! Case-insensitive comparison that also ignores '-' and '_'.
    CaseAndDashInsensitive
};
//! \endcond

/*! \libinternal \brief
 * Compare object for std::string STL containers and algorithms that supports
 * run-time decision on how to compare.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class StringCompare
{
public:
    /*! \brief
     * Creates a comparer with the given type
     *
     * This is not explicit, which allows passing \ref StringCompareType
     * directly to, e.g., `std::map` constructors.
     */
    StringCompare(StringCompareType type = StringCompareType::Exact) : type_(type) {}

    //! The comparison operation.
    bool operator()(const std::string& a, const std::string& b) const
    {
        switch (type_)
        {
            case StringCompareType::Exact: return a < b;
            case StringCompareType::CaseInsensitive:
                return gmx_strcasecmp(a.c_str(), b.c_str()) < 0;
            case StringCompareType::CaseAndDashInsensitive:
                return gmx_strcasecmp_min(a.c_str(), b.c_str()) < 0;
        }
        return a < b;
    }

private:
    StringCompareType type_;
};

} // namespace gmx

#endif
