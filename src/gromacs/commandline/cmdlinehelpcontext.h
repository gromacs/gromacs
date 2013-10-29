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
/*! \libinternal \file
 * \brief
 * Declares gmx::CommandLineHelpContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEHELPCONTEXT_H
#define GMX_COMMANDLINE_CMDLINEHELPCONTEXT_H

#include "../onlinehelp/helpwritercontext.h"
#include "../utility/common.h"

namespace gmx
{

/*! \libinternal \brief
 * Context information for writing out command-line help.
 *
 * This class wraps a HelpWriterContext, extending it with information specific
 * for command-line help export.  This way, code using only the routines in the
 * onlinehelp module is not exposed to extra features of the command-line help
 * export.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpContext
{
    public:
        /*! \brief
         * Creates the context.
         *
         * Wraps the constructor of HelpWriterContext.
         */
        CommandLineHelpContext(File *file, HelpOutputFormat format);
        ~CommandLineHelpContext();

        /*! \brief
         * Sets a display name for the module for which help is being written.
         *
         * \throws std::bad_alloc if out of memory.
         */
        void setModuleDisplayName(const std::string &name);
        /*! \brief
         * Sets the links to process in help output.
         *
         * A reference to \p links is stored until the CommandLineHelpContext
         * is destructed.  The caller is responsible for ensuring that the
         * links object remains valid long enough.
         */
        void setLinks(const HelpLinks &links);
        //! Sets whether hidden options should be shown in help output.
        void setShowHidden(bool bHidden);

        //! Returns the lower-level context for writing the help.
        const HelpWriterContext &writerContext() const;
        /*! \brief
         * Returns a display name for the module for which help is being written.
         *
         * Does not throw.
         */
        const char *moduleDisplayName() const;
        //! Returns whether hidden options should be shown in help output.
        bool showHidden() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

/*! \libinternal \brief
 * Helper for passing CommandLineHelpContext into parse_common_args().
 *
 * This class provides a mechanism to set and retrieve a global
 * CommandLineHelpContext object.  It is used to pass this object into
 * parse_common_args() from CommandLineModuleManager::runAsMainCMain() through
 * the main() function that is not aware of the wrapper binary mechanism.
 * It is not thread-safe because in this limited use case, it is always called
 * from a single-threaded context.
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class GlobalCommandLineHelpContext
{
    public:
        //! Returns the global context, or NULL if not set.
        static const CommandLineHelpContext *get();

        /*! \brief
         * Sets the global context for the scope.
         *
         * The global context is cleared when this object goes out of scope.
         *
         * It is an error to have more than one GlobalCommandLineHelpContext
         * object in existence at the same time.
         */
        explicit GlobalCommandLineHelpContext(const CommandLineHelpContext &context);
        //! Clears the global context.
        ~GlobalCommandLineHelpContext();
};

} // namespace gmx

#endif
