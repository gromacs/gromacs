/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <gtest/gtest.h>

#include "gromacs/coordinateio/flags.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/variant.h"

#include "gromacs/coordinateio/tests/coordinate_test.h"

namespace gmx
{

namespace test
{

/*!\brief
 * Helper function to split string with multiple values into option settings.
 *
 * \param[in] assigner     Options assigner to use.
 * \param[in] optionValues String to split that contains option values.
 * \param[in] isFloat      Need to know if entries are floating point numbers.
 */
static void
splitAndAddValues(OptionsAssigner *assigner, const std::string &optionValues, bool isFloat)
{
    std::vector<std::string> strings = splitString(optionValues);
    for (const auto &entry : strings)
    {
        if (isFloat)
        {
            real    value = std::stod(entry);
            Variant var(value);
            assigner->appendValue(var);
        }
        else
        {
            assigner->appendValue(entry);
        }
    }
}


/*! \libinternal \brief
 * Test fixture to test flag selections.
 */
class FlagTest : public ModuleTest,
                 public ModuleSelection
{
    public:
        FlagTest()
        {
            userOption_.initFileOptions(getOption());
        }

        /*! \brief
         * Set value for options in flag setting object.
         *
         * \param[in] optionName   Name of the option to add value for.
         * \param[in] optionValues Values to set for option.
         * \param[in] options      Container for options.
         * \param[in] isFloat      Need to know if entries are floating point numbers.
         */
        void setModuleFlag(const std::string &optionName,
                           const std::string &optionValues,
                           Options           *options,
                           bool               isFloat)
        {
            OptionsAssigner assigner(options);
            assigner.start();
            assigner.startOption(optionName.c_str());
            splitAndAddValues(&assigner, optionValues, isFloat);

            assigner.finishOption();
            assigner.finish();
        }


        /*! \brief
         * Get the modules that flag should have registered.
         * Calls the checkOptions code before registering.
         */
        OutputAdapters getModules()
        {
            userOption_.checkOptions();
            return userOption_.registerModules(getTopAtoms());
        }
    private:
        //! Object for setting options
        CoordinateOutputUserOptions userOption_;
};

} // namespace test

} // namespace gmx

#endif
