/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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

#include <memory>
#include <string>

#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class ShellCompletionWriter;
class TextWriter;

/*! \libinternal \brief
 * Context information for writing out command-line help.
 *
 * This class wraps a HelpWriterContext, extending it with information specific
 * for command-line help export.  This way, code using only the routines in the
 * onlinehelp module is not exposed to extra features of the command-line help
 * export.
 *
 * Copying a context works like with HelpWriterContext: the output file and
 * most state is shared.  However, setModuleDisplayName() and setShowHidden()
 * can be set independently for the child context.  Defaults for these options
 * are inherited from the parent.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpContext
{
public:
    /*! \brief
     * Creates a context for help export.
     *
     * Wraps the constructor of HelpWriterContext.
     */
    CommandLineHelpContext(TextWriter*        writer,
                           HelpOutputFormat   format,
                           const HelpLinks*   links,
                           const std::string& programName);
    //! Creates a context for a particular HelpWriterContext.
    explicit CommandLineHelpContext(const HelpWriterContext& writerContext);
    /*! \brief
     * Creates a context for shell completion.
     */
    explicit CommandLineHelpContext(ShellCompletionWriter* writer);
    //! Creates a copy of the context.
    explicit CommandLineHelpContext(const CommandLineHelpContext& other);
    //! Moves the context.
    CommandLineHelpContext(CommandLineHelpContext&& other) noexcept;
    //! Move-assigns the context.
    CommandLineHelpContext& operator=(CommandLineHelpContext&& other) noexcept;
    ~CommandLineHelpContext();

    /*! \brief
     * Sets a display name for the module for which help is being written.
     *
     * \throws std::bad_alloc if out of memory.
     */
    void setModuleDisplayName(const std::string& name);
    //! Sets whether hidden options should be shown in help output.
    void setShowHidden(bool bHidden);
    //! \copydoc HelpWriterContext::enterSubSection()
    void enterSubSection(const std::string& title);

    //! Returns the lower-level context for writing the help.
    const HelpWriterContext& writerContext() const;
    /*! \brief
     * Returns a display name for the module for which help is being written.
     *
     * Does not throw.
     */
    const char* moduleDisplayName() const;
    //! Returns whether hidden options should be shown in help output.
    bool showHidden() const;
    //! Returns whether this context is for exporting shell completions.
    bool isCompletionExport() const;
    /*! \brief
     * Returns the shell completion writer for this context.
     *
     * Can only be called if isCompletionExport() returns `true`.
     */
    ShellCompletionWriter& shellCompletionWriter() const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;

    GMX_DISALLOW_ASSIGN(CommandLineHelpContext);
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
    static const CommandLineHelpContext* get();

    /*! \brief
     * Sets the global context for the scope.
     *
     * The global context is cleared when this object goes out of scope.
     *
     * It is an error to have more than one GlobalCommandLineHelpContext
     * object in existence at the same time.
     */
    explicit GlobalCommandLineHelpContext(const CommandLineHelpContext& context);
    //! Clears the global context.
    ~GlobalCommandLineHelpContext();
};

} // namespace gmx

#endif
