/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::CommandLineHelpContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlinehelpcontext.h"

#include "gromacs/utility/gmxassert.h"

#include "shellcompletions.h"

namespace gmx
{

namespace
{

/*! \brief
 * Stores the global context set with GlobalCommandLineHelpContext.
 *
 * This is not protected by a mutex, since it is only used in command-line
 * start-up (i.e., single-threaded context), and is inherently not thread-safe.
 *
 * \ingroup module_commandline
 */
const CommandLineHelpContext *g_globalContext = NULL;

}   // namespace

/*! \internal \brief
 * Private implementation class for CommandLineHelpContext.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpContext::Impl
{
    public:
        //! Creates the implementation class and the low-level context.
        Impl(File *file, HelpOutputFormat format, const HelpLinks *links)
            : writerContext_(file, format, links), moduleDisplayName_("gmx"),
              completionWriter_(NULL), bHidden_(false)
        {
        }
        //! Creates an implementation class from a low-level context.
        explicit Impl(const HelpWriterContext &writerContext)
            : writerContext_(writerContext),
              completionWriter_(NULL), bHidden_(false)
        {
        }

        //! Wrapped lower-level context.
        HelpWriterContext      writerContext_;
        //! Display name for the module for which help is written.
        std::string            moduleDisplayName_;
        //! Shell completion writer (`NULL` if not doing completions).
        ShellCompletionWriter *completionWriter_;
        //! Whether hidden options should be shown in help output.
        bool                   bHidden_;
};

CommandLineHelpContext::CommandLineHelpContext(
        File *file, HelpOutputFormat format, const HelpLinks *links,
        const std::string &programName)
    : impl_(new Impl(file, format, links))
{
    impl_->writerContext_.setReplacement("[PROGRAM]", programName);
}

CommandLineHelpContext::CommandLineHelpContext(
        const HelpWriterContext &writerContext)
    : impl_(new Impl(writerContext))
{
}

CommandLineHelpContext::CommandLineHelpContext(
        ShellCompletionWriter *writer)
    : impl_(new Impl(writer->outputFile(), eHelpOutputFormat_Other, NULL))
{
    impl_->completionWriter_ = writer;
}

CommandLineHelpContext::CommandLineHelpContext(
        const CommandLineHelpContext &other)
    : impl_(new Impl(*other.impl_))
{
}

CommandLineHelpContext::~CommandLineHelpContext()
{
}

void CommandLineHelpContext::setModuleDisplayName(const std::string &name)
{
    impl_->writerContext_.setReplacement("[THISMODULE]", "[TT]" + name + "[tt]");
    impl_->moduleDisplayName_ = name;
}

void CommandLineHelpContext::setShowHidden(bool bHidden)
{
    impl_->bHidden_ = bHidden;
}

void CommandLineHelpContext::enterSubSection(const std::string &title)
{
    impl_->writerContext_.enterSubSection(title);
}

const HelpWriterContext &CommandLineHelpContext::writerContext() const
{
    return impl_->writerContext_;
}

const char *CommandLineHelpContext::moduleDisplayName() const
{
    return impl_->moduleDisplayName_.c_str();
}

bool CommandLineHelpContext::showHidden() const
{
    return impl_->bHidden_;
}

bool CommandLineHelpContext::isCompletionExport() const
{
    return impl_->completionWriter_ != NULL;
}

ShellCompletionWriter &CommandLineHelpContext::shellCompletionWriter() const
{
    GMX_RELEASE_ASSERT(isCompletionExport(),
                       "Invalid call when not writing shell completions");
    return *impl_->completionWriter_;
}

/********************************************************************
 * GlobalCommandLineHelpContext
 */

// static
const CommandLineHelpContext *GlobalCommandLineHelpContext::get()
{
    return g_globalContext;
}

GlobalCommandLineHelpContext::GlobalCommandLineHelpContext(
        const CommandLineHelpContext &context)
{
    GMX_RELEASE_ASSERT(g_globalContext == NULL,
                       "Global context set more than once");
    g_globalContext = &context;
}

GlobalCommandLineHelpContext::~GlobalCommandLineHelpContext()
{
    g_globalContext = NULL;
}

} // namespace gmx
