/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Implements gmx::IProgramContext and related methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/programcontext.h"

#include <cstddef>

#include <filesystem>

namespace gmx
{

namespace
{

//! \addtogroup module_utility
//! \{

/*! \brief
 * Default implementation of IProgramContext.
 *
 * This implementation is used if nothing has been set with
 * setProgramContext().
 *
 * Since it is constructed using a global initializer, it should not throw.
 */
class DefaultProgramContext : public IProgramContext
{
public:
    DefaultProgramContext() {}

    const char*            programName() const override { return "GROMACS"; }
    const char*            displayName() const override { return "GROMACS"; }
    std::filesystem::path  fullBinaryPath() const override { return ""; }
    InstallationPrefixInfo installationPrefix() const override
    {
        return InstallationPrefixInfo("", false);
    }
    const char* commandLine() const override { return ""; }
};

//! Global program info; stores the object set with setProgramContext().
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
const IProgramContext* g_programContext;
//! Default program context if nothing is set.
const DefaultProgramContext g_defaultContext;

//! \}

} // namespace

const IProgramContext& getProgramContext()
{
    if (g_programContext != nullptr)
    {
        return *g_programContext;
    }
    return g_defaultContext;
}

void setProgramContext(const IProgramContext* programContext)
{
    g_programContext = programContext;
}

} // namespace gmx
