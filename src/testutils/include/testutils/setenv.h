/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Helper functions to have identical behavior of setenv and unsetenv
 * on Unix and Windows systems.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \inlibraryapi
 * \ingroup module_testutils
 */

#include "config.h"

#include <cstdlib>

#ifndef GMX_TESTUTILS_SETENV_H
#    define GMX_TESTUTILS_SETENV_H

namespace gmx
{
namespace test
{
//! Workaround to make setenv work on Windows
inline int gmxSetenv(const char* name, const char* value, int overwrite)
{
#    if GMX_NATIVE_WINDOWS
    if (!overwrite)
    {
        size_t size  = 0;
        int    error = getenv_s(&size, nullptr, 0, name);
        if (error != 0 || size != 0)
        {
            return error;
        }
    }
    return _putenv_s(name, value);
#    else
    return setenv(name, value, overwrite);
#    endif
}

//! Workaround to make unsetenv work on Windows
inline int gmxUnsetenv(const char* name)
{
#    if GMX_NATIVE_WINDOWS
    return _putenv_s(name, "");
#    else
    return unsetenv(name);
#    endif
}
} // namespace test
} // namespace gmx

#endif // GMX_TESTUTILS_SETENV_H
