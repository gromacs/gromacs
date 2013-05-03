/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \ingroup module_integration_tests
 */
#include "integrationtests.h"
#include "integrationtests_impl.h"

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/utility/file.h"

namespace gmx
{
namespace test
{

/********************************************************************
 * IntegrationTestFixture
 */

IntegrationTestFixture::IntegrationTestFixture()
    : impl_(new Impl(this))
{
}

IntegrationTestFixture::~IntegrationTestFixture()
{
}

void
IntegrationTestFixture::createProgramCaller(std::string const& programName, ProgramCaller::cmain_func_ptr cmain)
{
    impl_->createProgramCaller(programName, cmain);
}

void
IntegrationTestFixture::setProgramOption(std::string const& programName, const char *optionName, std::string const& value)
{
    impl_->setProgramOption(programName, optionName, value);
}

void
IntegrationTestFixture::setProgramOption(std::string const& programName, const char *optionName, int const& value)
{
    impl_->setProgramOption(programName, optionName, gmx::formatString("%d", value));
}

void
IntegrationTestFixture::setProgramOption(std::string const& programName, const char *optionName, double const& value)
{
    impl_->setProgramOption(programName, optionName, gmx::formatString("%g", value));
}

template <class OptionType> void
IntegrationTestFixture::setProgramOption(std::string const& programName, const char *optionName, std::string const& value)
{
    impl_->setProgramOption<OptionType>(programName, optionName, value);
}

int
IntegrationTestFixture::runProgram(std::string const& programName)
{
    return impl_->getProgramCaller(programName)->run();
}

void
IntegrationTestFixture::redirectStringToStdin(const char* theString)
{
    std::string fakeStdin("fake-stdin");
    gmx::File::writeFileFromString(fakeStdin, theString);
    if (NULL == std::freopen(fakeStdin.c_str(), "r", stdin))
    {
        int cachedErrno = errno;
        GMX_THROW_WITH_ERRNO(FileIOError("Failed to redirect a string to stdin"),
                             "freopen",
                             cachedErrno
                             );
    }
}

// Template instantiation hints

template void IntegrationTestFixture::setProgramOption<IntegerOption>(std::string const& programName, const char *optionName, std::string const& value);
template void IntegrationTestFixture::setProgramOption<FileNameOption>(std::string const& programName, const char *optionName, std::string const& value);

} // namespace test
} // namespace gmx
