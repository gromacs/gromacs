/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Implements helper routines for setting environment variables in tests
 *
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/setenv.h"

#include "config.h"

#include <cstdlib>

#include "gromacs/utility/gmxassert.h"

namespace gmx
{
namespace test
{

namespace
{

//! Polyfiller to make setenv work on Windows
void gmxSetenv(const char* name, const char* value)
{
#if GMX_NATIVE_WINDOWS
    int result = _putenv_s(name, value);
#else
    static constexpr int sc_overwrite = 1;
    int                  result       = setenv(name, value, sc_overwrite);
#endif
    GMX_RELEASE_ASSERT(result == 0, "Failed to set environment variable");
}

//! Polyfiller to make unsetenv work on Windows
void gmxUnsetenv(const char* name)
{
#if GMX_NATIVE_WINDOWS
    int result = _putenv_s(name, "");
#else
    int result = unsetenv(name);
#endif
    GMX_RELEASE_ASSERT(result == 0, "Failed to unset environment variable");
}
} // namespace

GmxEnvGuard::GmxEnvGuard(const char* const envVar, const char* const newValue)
{
    GMX_RELEASE_ASSERT(envVar != nullptr, "Can't have null here");
    const char* oldValue = std::getenv(envVar);
    envVar_              = envVar;
    oldValue_            = oldValue ? std::make_optional<std::string>(oldValue) : std::nullopt;
    if (newValue)
    {
        gmxSetenv(envVar, newValue);
    }
    else
    {
        gmxUnsetenv(envVar);
    }
}

GmxEnvGuard::~GmxEnvGuard()
{
    if (oldValue_.has_value())
    {
        gmxSetenv(envVar_.c_str(), oldValue_->c_str());
    }
    else
    {
        gmxUnsetenv(envVar_.c_str());
    }
}


} // namespace test
} // namespace gmx
