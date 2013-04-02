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
 * Declares test fixture for mdrun tests
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_PROGRAM_CALLER_H
#define GMX_TESTUTILS_PROGRAM_CALLER_H

#include <string>
#include "gromacs/options/options.h"

namespace gmx
{

namespace test
{

/*! \internal \brief Class for encapsulating calls to GROMACS tools,
 * using simulated command-line arguments to their cmain functions.
 *
 * An object of this type can be instantiated directly. Any options
 * will then need to be defined, and then set, before the tool can
 * run. Example usage:
 *
 * // Any method in this class may throw std::bad_alloc if out of
 * memory.
 *
 * \ingroup module_integration_tests
 */
class ProgramCaller
{
    public:
        /*! \brief Helper typedef for C-style "main" functions
         */
        typedef int (*cmain_func_ptr)(int, char**);
        /*! \brief Construct a ProgramCaller object to call
         * _programName via function _cmain.
         */
        ProgramCaller(const char* _programName, cmain_func_ptr _cmain);
        /*! \brief Call the program with a command line constructed
         * from the current state of options.
         */
        int run();
        /*! \brief Defines an option for this program
         *
         * \param[in] settings Option description.
         * \returns   OptionInfo object for the created option (never NULL).
         * \throws    APIError if invalid option settings are provided.
         */
        OptionInfo *addOption(const AbstractOption &settings);
        //! Returns true if option \p name is defined.
        bool isDefined(const char *name) const;
        /*! \brief Sets the value for a previously defined option
         *
         * \param[in] settings Option description.
         * \throws    InvalidInputError if invalid option settings are provided.
         */
        void setOption(const char *optionName, std::string const& value);

    private:
        const char* programName;
        gmx::Options options;
        cmain_func_ptr cmain;
};

} // namespace test
} // namespace gmx

#endif
