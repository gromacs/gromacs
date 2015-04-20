/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::ProgramContextInterface and related methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "programcontext.h"

#include <cstddef>

namespace gmx
{

namespace
{

//! \addtogroup module_utility
//! \{

/*! \brief
 * Default implementation of ProgramContextInterface.
 *
 * This implementation is used if nothing has been set with
 * setProgramContext().
 *
 * Since it is constructed using a global initializer, it should not throw.
 */
class DefaultProgramContext : public ProgramContextInterface
{
    public:
        DefaultProgramContext() {}

        virtual const char *programName() const { return "GROMACS"; }
        virtual const char *displayName() const { return "GROMACS"; }
        virtual const char *fullBinaryPath() const { return ""; }
        virtual InstallationPrefixInfo installationPrefix() const
        {
            return InstallationPrefixInfo("", false);
        }
        virtual const char *commandLine() const { return ""; }
};

//! Global program info; stores the object set with setProgramContext().
const ProgramContextInterface *g_programContext;
//! Default program context if nothing is set.
const DefaultProgramContext    g_defaultContext;

//! \}

}   // namespace

const ProgramContextInterface &getProgramContext()
{
    if (g_programContext != NULL)
    {
        return *g_programContext;
    }
    return g_defaultContext;
}

void setProgramContext(const ProgramContextInterface *programContext)
{
    g_programContext = programContext;
}

} // namespace gmx
