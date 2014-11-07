/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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
 * Declares gmx::ShellCompletionWriter.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_SHELLCOMPLETIONS_H
#define GMX_COMMANDLINE_SHELLCOMPLETIONS_H

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class CommandLineHelpContext;
class File;
class Options;

//! \cond internal
//! \addtogroup module_commandline
//! \{

//! Output format for ShellCompletionWriter.
enum ShellCompletionFormat
{
    eShellCompletionFormat_Bash         //!< Shell completions for bash.
};

//! \}
//! \endcond

class ShellCompletionWriter
{
    public:
        typedef std::vector<std::string> ModuleNameList;

        ShellCompletionWriter(const std::string     &binaryName,
                              ShellCompletionFormat  format);
        ~ShellCompletionWriter();

        File *outputFile();

        void startCompletions();
        void writeModuleCompletions(const char    *moduleName,
                                    const Options &options);
        void writeWrapperCompletions(const ModuleNameList &modules,
                                     const Options        &options);
        void finishCompletions();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
