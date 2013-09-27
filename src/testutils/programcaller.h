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

/*! \internal \brief Class for maintaining options for calling GROMACS
 *  tools.
 *
 * An object of this type can be instantiated directly. Any options
 * will then need to be defined, and then set, before the tool can
 * run.
 *
 * \todo Example usage:
 *
 * Any method in this class may throw std::bad_alloc if out of memory.
 *
 * \ingroup module_testutils
 */
class ProgramCaller
{
    public:
        //! Constructor
        ProgramCaller(std::string const &_programName);
        //! Copy constructor
        ProgramCaller(ProgramCaller &c);
        /*! \brief Defines an option for this program
         *
         * \param[in] settings  Option description.
         * \returns   OptionInfo object for the created option (never NULL).
         * \throws    APIError if invalid option settings are provided.
         */
        OptionInfo *addOption(const AbstractOption &settings);
        /*! \brief Gets a reference to the options
         *
         * \throws    none
         */
        gmx::Options &getOptions();
        /*! \brief Gets a const reference to the program name
         *
         * \throws    none
         */
        std::string const &getProgramName();

    private:
        std::string const programName_;
        gmx::Options      options_;
};

} // namespace test
} // namespace gmx

#endif
