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
 * Declares gmx::CommandLineHelpWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEHELPWRITER_H
#define GMX_COMMANDLINE_CMDLINEHELPWRITER_H

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class CommandLineHelpContext;
class Options;

template <typename T> class ConstArrayRef;

/*! \brief
 * Writes help information for Options in ascii format.
 *
 * \inpublicapi
 * \ingroup module_commandline
 */
class CommandLineHelpWriter
{
    public:
        /*! \brief
         * Creates an object that writer ascii-formatted help for Options.
         *
         * \param[in] options  Options for which help should be printed.
         */
        explicit CommandLineHelpWriter(const Options &options);
        ~CommandLineHelpWriter();

        /*! \brief
         * Sets whether long descriptions for sections are shown in the help.
         */
        CommandLineHelpWriter &setShowDescriptions(bool bShow);
        /*! \brief
         * Sets time unit to show in descriptions.
         *
         * \param[in] timeUnit  Time unit to show in descriptions.
         * \throws    std::bad_alloc if out of memory.
         *
         * For each time parameter, any "%t" in the description is replaced
         * with \p timeunit.
         * If not called, uses a default "ps".
         */
        CommandLineHelpWriter &setTimeUnitString(const char *timeUnit);

        /*! \brief
         * Sets the list of known bugs/limitations.
         *
         * \param[in] bugs  Array of bugs/limitations.
         *
         * Each entry in the input array identifies a separate issue.
         * The array passed should remain valid for the lifetime of the writer
         * object.
         */
        CommandLineHelpWriter &
        setKnownIssues(const ConstArrayRef<const char *> &bugs);

        /*! \brief
         * Writes the help.
         *
         * \param[in] context  Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        void writeHelp(const CommandLineHelpContext &context);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
