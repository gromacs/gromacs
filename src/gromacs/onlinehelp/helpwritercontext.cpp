/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::HelpWriterContext.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_onlinehelp
 */
#include "helpwritercontext.h"

#include <cctype>

#include <algorithm>

#include <boost/shared_ptr.hpp>

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/wman.h"

#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programinfo.h"
#include "gromacs/utility/stringutil.h"

namespace
{

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

} // namespace

namespace gmx
{

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
        /*! \internal \brief
         * Shared, non-modifiable state for context objects.
         *
         * Contents of this structure are shared between all context objects
         * that are created from a common parent, e.g., by createSubsection().
         * This state should not be modified after construction.
         *
         * \ingroup module_onlinehelp
         */
        struct SharedState
        {
            //! Initializes the state with the given output file and format.
            SharedState(File *file, HelpOutputFormat format)
                : file_(*file), format_(format)
            {
            }

            //! Output file to which the help is written.
            File                   &file_;
            //! Output format for the help output.
            HelpOutputFormat        format_;
        };

        //! Smart pointer type for managing the shared state.
        typedef boost::shared_ptr<const SharedState> StatePointer;

        //! Initializes the context with the given state and section depth.
        Impl(const StatePointer &state, int sectionDepth)
            : state_(state), sectionDepth_(sectionDepth)
        {
        }

        //! Constant state shared by all child context objects.
        StatePointer    state_;
        //! Number of subsections above this context.
        int             sectionDepth_;

    private:
        GMX_DISALLOW_ASSIGN(Impl);
};

/********************************************************************
 * HelpWriterContext
 */

HelpWriterContext::HelpWriterContext(File *file, HelpOutputFormat format)
    : impl_(new Impl(Impl::StatePointer(new Impl::SharedState(file, format)), 0))
{
    if (format != eHelpOutputFormat_Console)
    {
        // TODO: Implement once the situation with Redmine issue #969 is more
        // clear.
        GMX_THROW(NotImplementedError(
                    "This output format is not implemented"));
    }
}

HelpWriterContext::HelpWriterContext(Impl *impl)
    : impl_(impl)
{
}

HelpWriterContext::HelpWriterContext(const HelpWriterContext &other)
    : impl_(new Impl(*other.impl_))
{
}

HelpWriterContext::~HelpWriterContext()
{
}

HelpWriterContext &HelpWriterContext::operator =(const HelpWriterContext &other)
{
    impl_.reset(new Impl(*other.impl_));
    return *this;
}

HelpOutputFormat HelpWriterContext::outputFormat() const
{
    return impl_->state_->format_;
}

File &HelpWriterContext::outputFile() const
{
    return impl_->state_->file_;
}

HelpWriterContext
HelpWriterContext::createSubsection(const std::string &title) const
{
    writeTitle(title);
    return HelpWriterContext(new Impl(impl_->state_, impl_->sectionDepth_ + 1));
}

std::string HelpWriterContext::substituteMarkup(const std::string &text) const
{
    char *resultStr = check_tty(text.c_str());
    scoped_ptr_sfree resultGuard(resultStr);
    return std::string(resultStr);
}

void HelpWriterContext::writeTitle(const std::string &title) const
{
    File &file = outputFile();
    file.writeLine(toUpperCase(title));
    file.writeLine();
}

void HelpWriterContext::writeTextBlock(const std::string &text) const
{
    TextLineWrapper wrapper;
    wrapper.settings().setLineLength(78);
    const char *program = ProgramInfo::getInstance().programName().c_str();
    std::string newText = replaceAll(text, "[PROGRAM]", program);
    outputFile().writeLine(wrapper.wrapToString(substituteMarkup(newText)));
}

} // namespace gmx
