/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::CommandLineParser.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEPARSER_H
#define GMX_COMMANDLINE_CMDLINEPARSER_H

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"

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
         * Makes the parser skip unknown options and keep them in \c argv.
         *
         * \param[in] bEnabled  Whether to skip and keep unknown options.
         * \returns   *this
         *
         * Setting this option to true has dual effect: unknown options are
         * silently skipped, and all recognized options are removed from
         * \c argc and \c argv in parse().  These effects should be easy to
         * separate into different flags if there is need for it.
         *
         * The default is false: unknown options result in exceptions and
         * \c argc and \c argv are not modified.
         *
         * Does not throw.
         */
        CommandLineParser &skipUnknown(bool bEnabled);

        /*! \brief
         * Parses the command line.
         *
         * \throws  std::bad_alloc if out of memory.
         * \throws  InvalidInputError if any errors were detected in the input.
         *
         * All command-line arguments are parsed, and an aggregate exception
         * with all the detected errors is thrown in the end.
         *
         * If skipUnknown() is false, the input parameters are not modified.
         * If skipUnknown() is true, recognized options and their values are
         * removed from the argument list.  \c argv[0] is never modified.
         */
        void parse(int *argc, char *argv[]);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
