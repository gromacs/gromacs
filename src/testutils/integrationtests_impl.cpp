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
 * Implements pimpl class in integrationtests_impl.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_integration_tests
 */
#include "integrationtests_impl.h"

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "testutils/testoptions.h"

namespace gmx
{
namespace test
{
 
IntegrationTestFixture::Impl::Impl(IntegrationTestFixture *parent)
    : programCallers(),
      parent_(*parent)
{
}

template <class OptionType>
void
IntegrationTestFixture::Impl::defineProgramOption(std::string const& programName, const char *optionName, char const* description = "")
{
    try
    {
        getProgramCaller(programName)->addOption(OptionType(optionName).description(description));
    }
    catch(gmx::APIError &ex)
    {
        // TODO Is this useful, correct, etc.?
        ex.prependContext("While constructing call to " + std::string(programName)
                          + " to define option " + optionName);
    }
}

template <class OptionType> void
IntegrationTestFixture::Impl::setProgramOption(std::string const& programName, const char *optionName, std::string const& value)
{
    // Need the definition in a header file, so the templates
    // can be instantiated at link time.
    defineProgramOption<OptionType>(programName, optionName, "");
    setProgramOption(programName, optionName, value);
}

void
IntegrationTestFixture::Impl::setProgramOption(std::string const& programName, const char *optionName, std::string const& value)
{
    GMX_RELEASE_ASSERT(false != getProgramCaller(programName)->isDefined(optionName),
                       formatString("Option %s for program %s should have been defined before the attempt to set it",
                                    optionName, programName.c_str()).c_str());

    try
    {
        getProgramCaller(programName)->setOption(optionName, value);
    }
    catch(gmx::InvalidInputError &ex)
    {
        // TODO Is this useful, correct, etc.?
        ex.prependContext("While constructing call to " + std::string(programName)
                          + " to add option " + optionName
                          + " with value " + value);
    }
}

void
IntegrationTestFixture::Impl::createProgramCaller(std::string const& programName, ProgramCaller::cmain_func_ptr cmain)
{
    programCallers[programName] = boost::shared_ptr<ProgramCaller>(new ProgramCaller(programName.c_str(), cmain));
}

boost::shared_ptr<ProgramCaller>
IntegrationTestFixture::Impl::getProgramCaller(std::string const& programName)
{
    // Check that the caller for this program has been defined already
    if (programCallers.end() == programCallers.find(programName))
    {
        GMX_THROW(APIError("In IntegrationTestFixture::createProgramCaller(std::string const&), " + programName + " cannot be created automatically. Use createProgramCaller() before getProgramCaller()."));
    }
    return programCallers[programName];
}

// Template instantiation hints

template void IntegrationTestFixture::Impl::setProgramOption<FileNameOption>(std::string const& programName, const char *optionName, std::string const& value);
template void IntegrationTestFixture::Impl::setProgramOption<IntegerOption>(std::string const& programName, const char *optionName, std::string const& value);
template void IntegrationTestFixture::Impl::defineProgramOption<IntegerOption>(std::string const& programName, const char *optionName, char const* description);

} // namespace test
} // namespace gmx
