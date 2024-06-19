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
 * Implements classes in stdiohelper.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/stdiohelper.h"

#include <cerrno>
#include <cstdio>

#include <filesystem>
#include <string>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/textwriter.h"

#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * StdioTestHelper
 */

StdioTestHelper::~StdioTestHelper()
{
    if (redirected)
    {
        fclose(stdin);
    }
}

void StdioTestHelper::redirectStringToStdin(const char* theString)
{
    const std::string fakeStdin = fileManager_.getTemporaryFilePath(".stdin").string();
    gmx::TextWriter::writeFileFromString(fakeStdin, theString);
    if (nullptr == std::freopen(fakeStdin.c_str(), "r", stdin))
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Failed to redirect a string to stdin"), "freopen", errno);
    }
    redirected = true;
}

} // namespace test
} // namespace gmx
