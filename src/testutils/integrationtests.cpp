/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Implements classes in integrationtests.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "integrationtests.h"

#include "testutils/testoptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/utility/file.h"
#include <stdlib.h>
#include <stdio.h>

namespace gmx
{
namespace test
{

/********************************************************************
 * IntegrationTestFixture
 */

std::string IntegrationTestFixture::s_maxBackup("-1");

//! \cond
GMX_TEST_OPTIONS(IntegrationTestOptions, options)
{
    options->addOption(StringOption("max-backup")
                           .store(&IntegrationTestFixture::s_maxBackup)
                           .description("Maximum number of backup files of old test output to write (-1 prevents backups being created)"));
}
//! \endcond

IntegrationTestFixture::IntegrationTestFixture()
{
    // TODO fix this when we have an encapsulation layer for handling
    // environment variables
#ifdef GMX_NATIVE_WINDOWS
    _putenv(("GMX_MAXBACKUP="+s_maxBackup).c_str());
#else
    setenv("GMX_MAXBACKUP", s_maxBackup.c_str(), true);
#endif
}

IntegrationTestFixture::~IntegrationTestFixture()
{
}

void
IntegrationTestFixture::redirectStringToStdin(const char* theString)
{
    std::string fakeStdin("fake-stdin");
    gmx::File::writeFileFromString(fakeStdin, theString);
    if (NULL == std::freopen(fakeStdin.c_str(), "r", stdin))
    {
        GMX_THROW_WITH_ERRNO(FileIOError("Failed to redirect a string to stdin"),
                             "freopen",
                             errno);
    }
}

} // namespace test
} // namespace gmx
