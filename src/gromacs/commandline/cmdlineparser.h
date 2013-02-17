/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::CommandLineParser.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEPARSER_H
#define GMX_COMMANDLINE_CMDLINEPARSER_H

#include <string>
#include <vector>

#include "../utility/common.h"

namespace gmx
{

class Options;

/*! \brief
 * Implements command-line parsing for Options objects.
 *
 * Typical usage (without error checking):
 * \code
   gmx::Options options("name", "description");
   // Fill up options

   gmx::CommandLineParser(&options).parse(&argc, argv);
   options.finish();
 * \endcode
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class CommandLineParser
{
    public:
        /*! \brief
         * Creates a command-line parser that sets values for options.
         *
         * \param[in] options  Options object whose options should be set.
         * \throws  std::bad_alloc if out of memory.
         */
        CommandLineParser(Options *options);
        ~CommandLineParser();

        /*! \brief
         * Parses the command line.
         *
         * \throws  std::bad_alloc if out of memory.
         * \throws  InvalidInputError if any errors were detected in the input.
         *
         * All command-line arguments are parsed, and an aggregate exception
         * with all the detected errors is thrown in the end.
         *
         * Currently, the input parameters are not modified, but this may
         * change if/when support for parsing only part of the options is
         * implemented.
         */
        void parse(int *argc, char *argv[]);
        /*! \brief
         * Parses the command line from a std::vector.
         *
         * \param[in] commandLine  Array of command-line strings.
         * \throws  std::bad_alloc if out of memory.
         * \throws  InvalidInputError if any errors were detected in the input.
         *
         * \p commandLine should relate to the standard \c argv array
         * one-to-one.
         *
         * This method is provided for convenience for cases where the command
         * line needs to be stored before parsing.
         *
         * Currently, the input parameters are not modified, but this may
         * change if/when support for parsing only part of the options is
         * implemented.
         *
         * \see parse(int *, char *[])
         */
        void parse(std::vector<std::string> *commandLine);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
