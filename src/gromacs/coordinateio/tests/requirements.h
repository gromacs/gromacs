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
/*!\file
 * \libinternal
 * \brief
 * Helpers and data for flag setting method.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_coordinateio
 */

#ifndef GMX_COORDINATEIO_TESTS_FLAG_H
#define GMX_COORDINATEIO_TESTS_FLAG_H

#include "gmxpre.h"

#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "gromacs/coordinateio/requirements.h"
#include "gromacs/coordinateio/tests/coordinate_test.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsassigner.h"
#include "gromacs/utility/any.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"

namespace gmx
{

template<typename>
class ArrayRef;

namespace test
{

//! Test enums to decide on how to use options framework
enum class TestEnums
{
    efTestString,
    efTestInt,
    efTestFloat
};

/*!\brief
 * Helper function to split string with multiple values into option settings.
 *
 * \param[in] assigner     Options assigner to use.
 * \param[in] optionValues String to split that contains option values.
 * \param[in] type         Determine which input time we are using.
 */
static void splitAndAddValues(OptionsAssigner* assigner, const std::string& optionValues, TestEnums type)
{
    std::vector<std::string> strings = splitString(optionValues);
    for (const auto& entry : strings)
    {
        switch (type)
        {
            case (TestEnums::efTestFloat):
            {
                real value = std::stod(entry);
                Any  var(value);
                assigner->appendValue(var);
                break;
            }
            case (TestEnums::efTestString):
            {
                assigner->appendValue(entry);
                break;
            }
            case (TestEnums::efTestInt):
            {
                int value = std::stoi(entry);
                Any var(value);
                assigner->appendValue(var);
                break;
            }
            default: GMX_THROW(InternalError("No such input type"));
        }
    }
}


/*! \libinternal \brief
 * Test fixture to test user inputs.
 */
class FlagTest : public gmx::test::CommandLineTestBase
{
public:
    FlagTest() { requirementsBuilder_.initOptions(&options_); }

    /*! \brief
     * Set value for options in flag setting object.
     *
     * \param[in] optionName   Name of the option to add value for.
     * \param[in] optionValues Values to set for option.
     * \param[in] options      Container for options.
     * \param[in] type         Need to know type of entries.
     */
    static void setModuleFlag(const std::string& optionName,
                              const std::string& optionValues,
                              Options*           options,
                              TestEnums          type)
    {
        OptionsAssigner assigner(options);
        assigner.start();
        assigner.startOption(optionName.c_str());
        splitAndAddValues(&assigner, optionValues, type);

        assigner.finishOption();
        assigner.finish();
    }

    //! Storage of requirements.
    OutputRequirementOptionDirector requirementsBuilder_;
    //! Options storage
    Options options_;
};

} // namespace test

} // namespace gmx

#endif
