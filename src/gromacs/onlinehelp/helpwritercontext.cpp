/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
 * Implements gmx::HelpWriterContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "helpwritercontext.h"

#include <cctype>

#include <algorithm>
#include <string>

#include "gromacs/legacyheaders/smalloc.h"

#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/onlinehelp/wman.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

/*! \internal \brief
 * Stores the global context set with HelpWriterGlobalContext.
 *
 * This is not protected by a mutex, since it is only used in command-line
 * start-up (i.e., single-threaded context), and is inherently not thread-safe.
 */
const HelpWriterContext *g_globalContext = NULL;

/*! \internal \brief
 * Make the string uppercase.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all characters transformed to uppercase.
 * \throws    std::bad_alloc if out of memory.
 */
std::string toUpperCase(const std::string &text)
{
    std::string result(text);
    transform(result.begin(), result.end(), result.begin(), toupper);
    return result;
}

}   // namespace

/********************************************************************
 * HelpWriterContext::Impl
 */

/*! \internal \brief
 * Private implementation class for HelpWriterContext.
 *
 * \ingroup module_onlinehelp
 */
class HelpWriterContext::Impl
{
    public:
        //! Initializes the context with the given output file and format.
        explicit Impl(File *file, HelpOutputFormat format)
            : file_(*file), format_(format)
        {
        }

        //! Output file to which the help is written.
        File                   &file_;
        //! Output format for the help output.
        HelpOutputFormat        format_;
        //! Display name for the module for which help is written.
        std::string             moduleDisplayName_;
};

/********************************************************************
 * HelpWriterContext
 */

HelpWriterContext::HelpWriterContext(File *file, HelpOutputFormat format)
    : impl_(new Impl(file, format))
{
}

HelpWriterContext::~HelpWriterContext()
{
}

void HelpWriterContext::setModuleDisplayName(const std::string &name)
{
    impl_->moduleDisplayName_ = name;
}

HelpOutputFormat HelpWriterContext::outputFormat() const
{
    return impl_->format_;
}

File &HelpWriterContext::outputFile() const
{
    return impl_->file_;
}

const char *HelpWriterContext::moduleDisplayName() const
{
    return impl_->moduleDisplayName_.empty()
           ? "gmx"
           : impl_->moduleDisplayName_.c_str();
}

std::string HelpWriterContext::substituteMarkup(const std::string &text) const
{
    if (outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                          "This output format is not implemented"));
    }
    char            *resultStr = check_tty(text.c_str());
    scoped_ptr_sfree resultGuard(resultStr);
    return std::string(resultStr);
}

void HelpWriterContext::writeTitle(const std::string &title) const
{
    if (outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                          "This output format is not implemented"));
    }
    File &file = outputFile();
    file.writeLine(toUpperCase(title));
    file.writeLine();
}

void HelpWriterContext::writeTextBlock(const std::string &text) const
{
    if (outputFormat() != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                          "This output format is not implemented"));
    }
    TextLineWrapper wrapper;
    wrapper.settings().setLineLength(78);
    const char     *program = ProgramInfo::getInstance().programName().c_str();
    std::string     newText = replaceAll(text, "[PROGRAM]", program);
    outputFile().writeLine(wrapper.wrapToString(substituteMarkup(newText)));
}

/********************************************************************
 * HelpWriterGlobalContext
 */

// static
const HelpWriterContext *HelpWriterGlobalContext::get()
{
    return g_globalContext;
}

HelpWriterGlobalContext::HelpWriterGlobalContext(const HelpWriterContext &context)
{
    GMX_RELEASE_ASSERT(g_globalContext == NULL,
                       "Global context set more than once");
    g_globalContext = &context;
}

HelpWriterGlobalContext::~HelpWriterGlobalContext()
{
    g_globalContext = NULL;
}

} // namespace gmx
