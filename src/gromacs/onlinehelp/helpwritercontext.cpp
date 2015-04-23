/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements gmx::HelpWriterContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_onlinehelp
 */
#include "gmxpre.h"

#include "helpwritercontext.h"

#include <cctype>

#include <algorithm>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "gromacs/onlinehelp/helpformat.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

#include "rstparser.h"

namespace gmx
{

namespace
{

//! \internal \addtogroup module_onlinehelp
//! \{

//! Characters used for reStructuredText title underlining.
const char g_titleChars[] = "=-^*~+#'_.";

struct t_sandr
{
    const char *search;
    const char *replace;
};

/* The order of these arrays is significant. Text search and replace
 * for each element occurs in order, so earlier changes can induce
 * subsequent changes even though the original text might not appear
 * to invoke the latter changes.
 * TODO: Get rid of this behavior. It makes it very difficult to manage
 * replacements coming from multiple sources (e.g., hyperlinks).*/

//! List of replacements for console output.
const t_sandr sandrTty[] = {
    { "\\*", "*" },
    { "\\=", "=" },
    { "[REF]", "" },
    { "[ref]", "" },
    { "[TT]", "" },
    { "[tt]", "" },
    { "[BB]", "" },
    { "[bb]", "" },
    { "[IT]", "" },
    { "[it]", "" },
    { "[MATH]", "" },
    { "[math]", "" },
    { "[CHEVRON]", "<" },
    { "[chevron]", ">" },
    { "[MAG]", "|" },
    { "[mag]", "|" },
    { "[INT]", "integral" },
    { "[FROM]", " from " },
    { "[from]", "" },
    { "[TO]", " to " },
    { "[to]", " of" },
    { "[int]", "" },
    { "[SUM]", "sum" },
    { "[sum]", "" },
    { "[SUB]", "_" },
    { "[sub]", "" },
    { "[SQRT]", "sqrt(" },
    { "[sqrt]", ")" },
    { "[EXP]", "exp(" },
    { "[exp]", ")" },
    { "[LN]", "ln(" },
    { "[ln]", ")" },
    { "[LOG]", "log(" },
    { "[log]", ")" },
    { "[COS]", "cos(" },
    { "[cos]", ")" },
    { "[SIN]", "sin(" },
    { "[sin]", ")" },
    { "[TAN]", "tan(" },
    { "[tan]", ")" },
    { "[COSH]", "cosh(" },
    { "[cosh]", ")" },
    { "[SINH]", "sinh(" },
    { "[sinh]", ")" },
    { "[TANH]", "tanh(" },
    { "[tanh]", ")" },
    { "[PAR]", "\n\n" },
    { "[GRK]", "" },
    { "[grk]", "" }
};

//! List of replacements for reStructuredText output.
const t_sandr sandrRst[] = {
    { "[TT]", "``" },
    { "[tt]", "``" },
    { "[BB]", "**" },
    { "[bb]", "**" },
    { "[IT]", "*" },
    { "[it]", "*" },
    { "[MATH]", "" },
    { "[math]", "" },
    { "[CHEVRON]", "<" },
    { "[chevron]", ">" },
    { "[MAG]", "\\|" },
    { "[mag]", "\\|" },
    { "[INT]", "integral" },
    { "[FROM]", " from " },
    { "[from]", "" },
    { "[TO]", " to " },
    { "[to]", " of" },
    { "[int]", "" },
    { "[SUM]", "sum" },
    { "[sum]", "" },
    { "[SUB]", "_" },
    { "[sub]", "" },
    { "[SQRT]", "sqrt(" },
    { "[sqrt]", ")" },
    { "[EXP]", "exp(" },
    { "[exp]", ")" },
    { "[LN]", "ln(" },
    { "[ln]", ")" },
    { "[LOG]", "log(" },
    { "[log]", ")" },
    { "[COS]", "cos(" },
    { "[cos]", ")" },
    { "[SIN]", "sin(" },
    { "[sin]", ")" },
    { "[TAN]", "tan(" },
    { "[tan]", ")" },
    { "[COSH]", "cosh(" },
    { "[cosh]", ")" },
    { "[SINH]", "sinh(" },
    { "[sinh]", ")" },
    { "[TANH]", "tanh(" },
    { "[tanh]", ")" },
    { "[PAR]", "\n\n" },
    { "[GRK]", "" },
    { "[grk]", "" }
};

/*! \brief
 * Replaces all entries from a list of replacements.
 */
std::string repall(const std::string &s, int nsr, const t_sandr sa[])
{
    std::string result(s);
    for (int i = 0; i < nsr; ++i)
    {
        result = replaceAll(result, sa[i].search, sa[i].replace);
    }
    return result;
}

/*! \brief
 * Replaces all entries from a list of replacements.
 */
template <size_t nsr>
std::string repall(const std::string &s, const t_sandr (&sa)[nsr])
{
    return repall(s, nsr, sa);
}

/*! \brief
 * Custom output interface for HelpWriterContext::Impl::processMarkup().
 *
 * Provides an interface that is used to implement different types of output
 * from HelpWriterContext::Impl::processMarkup().
 */
class WrapperInterface
{
    public:
        virtual ~WrapperInterface() {}

        /*! \brief
         * Provides the wrapping settings.
         *
         * HelpWriterContext::Impl::processMarkup() may provide some default
         * values for the settings if they are not set; this is the reason the
         * return value is not const.
         */
        virtual TextLineWrapperSettings &settings() = 0;
        //! Appends the given string to output.
        virtual void wrap(const std::string &text)  = 0;
};

/*! \brief
 * Wraps markup output into a single string.
 */
class WrapperToString : public WrapperInterface
{
    public:
        //! Creates a wrapper with the given settings.
        explicit WrapperToString(const TextLineWrapperSettings &settings)
            : wrapper_(settings)
        {
        }

        virtual TextLineWrapperSettings &settings()
        {
            return wrapper_.settings();
        }
        virtual void wrap(const std::string &text)
        {
            result_.append(wrapper_.wrapToString(text));
        }
        //! Returns the result string.
        const std::string &result() const { return result_; }

    private:
        TextLineWrapper         wrapper_;
        std::string             result_;
};

/*! \brief
 * Wraps markup output into a vector of string (one line per element).
 */
class WrapperToVector : public WrapperInterface
{
    public:
        //! Creates a wrapper with the given settings.
        explicit WrapperToVector(const TextLineWrapperSettings &settings)
            : wrapper_(settings)
        {
        }

        virtual TextLineWrapperSettings &settings()
        {
            return wrapper_.settings();
        }
        virtual void wrap(const std::string &text)
        {
            const std::vector<std::string> &lines = wrapper_.wrapToVector(text);
            result_.insert(result_.end(), lines.begin(), lines.end());
        }
        //! Returns a vector with the output lines.
        const std::vector<std::string> &result() const { return result_; }

    private:
        TextLineWrapper          wrapper_;
        std::vector<std::string> result_;
};

/*! \brief
 * Makes the string uppercase.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all characters transformed to uppercase.
 * \throws    std::bad_alloc if out of memory.
 */
std::string toUpperCase(const std::string &text)
{
    std::string result(text);
    std::transform(result.begin(), result.end(), result.begin(), toupper);
    return result;
}

/*! \brief
 * Removes extra newlines from reStructuredText.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all sequences of more than two newlines replaced
 *     with just two newlines.
 * \throws    std::bad_alloc if out of memory.
 */
std::string removeExtraNewlinesRst(const std::string &text)
{
    // Start from 2, so that all newlines in the beginning get stripped off.
    int         newlineCount = 2;
    std::string result;
    result.reserve(text.length());
    for (size_t i = 0; i < text.length(); ++i)
    {
        if (text[i] == '\n')
        {
            ++newlineCount;
            if (newlineCount > 2)
            {
                continue;
            }
        }
        else
        {
            newlineCount = 0;
        }
        result.push_back(text[i]);
    }
    size_t last = result.find_last_not_of('\n');
    if (last != std::string::npos)
    {
        result.resize(last + 1);
    }
    return result;
}

//! \}

}   // namespace

/********************************************************************
 * HelpLinks::Impl
 */

/*! \internal \brief
 * Private implementation class for HelpLinks.
 *
 * \ingroup module_onlinehelp
 */
class HelpLinks::Impl
{
    public:
        struct LinkItem
        {
            LinkItem(const std::string &linkName,
                     const std::string &replacement)
                : linkName(linkName), replacement(replacement)
            {
            }
            std::string         linkName;
            std::string         replacement;
        };

        //! Shorthand for a list of links.
        typedef std::vector<LinkItem> LinkList;

        //! Initializes empty links with the given format.
        explicit Impl(HelpOutputFormat format) : format_(format)
        {
        }

        //! List of links.
        LinkList          links_;
        //! Output format for which the links are formatted.
        HelpOutputFormat  format_;
};

/********************************************************************
 * HelpLinks
 */

HelpLinks::HelpLinks(HelpOutputFormat format) : impl_(new Impl(format))
{
}

HelpLinks::~HelpLinks()
{
}

void HelpLinks::addLink(const std::string &linkName,
                        const std::string &targetName,
                        const std::string &displayName)
{
    std::string replacement;
    switch (impl_->format_)
    {
        case eHelpOutputFormat_Console:
            replacement = repall(displayName, sandrTty);
            break;
        case eHelpOutputFormat_Rst:
            replacement = targetName;
            break;
        default:
            GMX_RELEASE_ASSERT(false, "Output format not implemented for links");
    }
    impl_->links_.push_back(Impl::LinkItem(linkName, replacement));
}

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
        /*! \brief
         * Shared, non-modifiable state for context objects.
         *
         * Contents of this structure are shared between all context objects
         * that are created from a common parent.
         * This state should not be modified after construction.
         *
         * \ingroup module_onlinehelp
         */
        class SharedState
        {
            public:
                //! Initializes the state with the given parameters.
                SharedState(File *file, HelpOutputFormat format,
                            const HelpLinks *links)
                    : file_(*file), format_(format), links_(links)
                {
                }

                /*! \brief
                 * Returns a formatter for formatting options lists for console
                 * output.
                 *
                 * The formatter is lazily initialized on first access.
                 */
                TextTableFormatter &consoleOptionsFormatter() const
                {
                    GMX_RELEASE_ASSERT(format_ == eHelpOutputFormat_Console,
                                       "Accessing console formatter for non-console output");
                    if (!consoleOptionsFormatter_)
                    {
                        consoleOptionsFormatter_.reset(new TextTableFormatter());
                        consoleOptionsFormatter_->setFirstColumnIndent(1);
                        consoleOptionsFormatter_->addColumn(NULL, 7, false);
                        consoleOptionsFormatter_->addColumn(NULL, 18, false);
                        consoleOptionsFormatter_->addColumn(NULL, 16, false);
                        consoleOptionsFormatter_->addColumn(NULL, 28, false);
                    }
                    return *consoleOptionsFormatter_;
                }

                //! Output file to which the help is written.
                File                   &file_;
                //! Output format for the help output.
                HelpOutputFormat        format_;
                //! Links to use.
                const HelpLinks        *links_;

            private:
                //! Formatter for console output options.
                mutable boost::scoped_ptr<TextTableFormatter> consoleOptionsFormatter_;
        };

        struct ReplaceItem
        {
            ReplaceItem(const std::string &search,
                        const std::string &replace)
                : search(search), replace(replace)
            {
            }
            std::string         search;
            std::string         replace;
        };

        //! Smart pointer type for managing the shared state.
        typedef boost::shared_ptr<const SharedState> StatePointer;
        //! Shorthand for a list of markup/other replacements.
        typedef std::vector<ReplaceItem> ReplaceList;

        //! Initializes the context with the given state and section depth.
        Impl(const StatePointer &state, int sectionDepth)
            : state_(state), sectionDepth_(sectionDepth)
        {
        }

        //! Adds a new replacement.
        void addReplacement(const std::string &search,
                            const std::string &replace)
        {
            replacements_.push_back(ReplaceItem(search, replace));
        }

        //! Replaces links in a given string.
        std::string replaceLinks(const std::string &input) const;

        /*! \brief
         * Process markup and wrap lines within a block of text.
         *
         * \param[in] text     Text to process.
         * \param     wrapper  Object used to wrap the text.
         *
         * The \p wrapper should take care of either writing the text to output
         * or providing an interface for the caller to retrieve the output.
         */
        void processMarkup(const std::string &text,
                           WrapperInterface  *wrapper) const;

        //! Constant state shared by all child context objects.
        StatePointer            state_;
        //! List of markup/other replacements.
        ReplaceList             replacements_;
        //! Number of subsections above this context.
        int                     sectionDepth_;

    private:
        GMX_DISALLOW_ASSIGN(Impl);
};

std::string HelpWriterContext::Impl::replaceLinks(const std::string &input) const
{
    std::string result(input);
    if (state_->links_ != NULL)
    {
        HelpLinks::Impl::LinkList::const_iterator link;
        for (link  = state_->links_->impl_->links_.begin();
             link != state_->links_->impl_->links_.end(); ++link)
        {
            result = replaceAllWords(result, link->linkName, link->replacement);
        }
    }
    return result;
}

void HelpWriterContext::Impl::processMarkup(const std::string &text,
                                            WrapperInterface  *wrapper) const
{
    std::string result(text);
    for (ReplaceList::const_iterator i = replacements_.begin();
         i != replacements_.end(); ++i)
    {
        result = replaceAll(result, i->search, i->replace);
    }
    switch (state_->format_)
    {
        case eHelpOutputFormat_Console:
        {
            const int   baseFirstLineIndent = wrapper->settings().firstLineIndent();
            const int   baseIndent          = wrapper->settings().indent();
            result = repall(result, sandrTty);
            result = replaceLinks(result);
            std::string          paragraph;
            paragraph.reserve(result.length());
            RstParagraphIterator iter(result);
            while (iter.nextParagraph())
            {
                iter.getParagraphText(&paragraph);
                wrapper->settings().setFirstLineIndent(baseFirstLineIndent + iter.firstLineIndent());
                wrapper->settings().setIndent(baseIndent + iter.indent());
                wrapper->wrap(paragraph);
            }
            wrapper->settings().setFirstLineIndent(baseFirstLineIndent);
            wrapper->settings().setIndent(baseIndent);
            break;
        }
        case eHelpOutputFormat_Rst:
        {
            result = repall(result, sandrRst);
            result = replaceLinks(result);
            result = replaceAll(result, "[REF]", "");
            result = replaceAll(result, "[ref]", "");
            result = removeExtraNewlinesRst(result);
            wrapper->wrap(result);
            break;
        }
        default:
            GMX_THROW(InternalError("Invalid help output format"));
    }
}

/********************************************************************
 * HelpWriterContext
 */

HelpWriterContext::HelpWriterContext(File *file, HelpOutputFormat format)
    : impl_(new Impl(Impl::StatePointer(new Impl::SharedState(file, format, NULL)), 0))
{
}

HelpWriterContext::HelpWriterContext(File *file, HelpOutputFormat format,
                                     const HelpLinks *links)
    : impl_(new Impl(Impl::StatePointer(new Impl::SharedState(file, format, links)), 0))
{
    if (links != NULL)
    {
        GMX_RELEASE_ASSERT(links->impl_->format_ == format,
                           "Links must have the same output format as the context");
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

void HelpWriterContext::setReplacement(const std::string &search,
                                       const std::string &replace)
{
    impl_->addReplacement(search, replace);
}

HelpOutputFormat HelpWriterContext::outputFormat() const
{
    return impl_->state_->format_;
}

File &HelpWriterContext::outputFile() const
{
    return impl_->state_->file_;
}

void HelpWriterContext::enterSubSection(const std::string &title)
{
    GMX_RELEASE_ASSERT(impl_->sectionDepth_ - 1 < static_cast<int>(std::strlen(g_titleChars)),
                       "Too deeply nested subsections");
    writeTitle(title);
    ++impl_->sectionDepth_;
}

std::string
HelpWriterContext::substituteMarkupAndWrapToString(
        const TextLineWrapperSettings &settings, const std::string &text) const
{
    WrapperToString wrapper(settings);
    impl_->processMarkup(text, &wrapper);
    return wrapper.result();
}

std::vector<std::string>
HelpWriterContext::substituteMarkupAndWrapToVector(
        const TextLineWrapperSettings &settings, const std::string &text) const
{
    WrapperToVector wrapper(settings);
    impl_->processMarkup(text, &wrapper);
    return wrapper.result();
}

void HelpWriterContext::writeTitle(const std::string &title) const
{
    if (title.empty())
    {
        return;
    }
    File &file = outputFile();
    switch (outputFormat())
    {
        case eHelpOutputFormat_Console:
            file.writeLine(toUpperCase(title));
            file.writeLine();
            break;
        case eHelpOutputFormat_Rst:
            file.writeLine(title);
            file.writeLine(std::string(title.length(),
                                       g_titleChars[impl_->sectionDepth_]));
            break;
        default:
            GMX_THROW(NotImplementedError(
                              "This output format is not implemented"));
    }
}

void HelpWriterContext::writeTextBlock(const std::string &text) const
{
    TextLineWrapperSettings settings;
    if (outputFormat() == eHelpOutputFormat_Console)
    {
        settings.setLineLength(78);
    }
    outputFile().writeLine(substituteMarkupAndWrapToString(settings, text));
}

void HelpWriterContext::writeOptionListStart() const
{
}

void HelpWriterContext::writeOptionItem(const std::string &name,
                                        const std::string &value,
                                        const std::string &defaultValue,
                                        const std::string &info,
                                        const std::string &description) const
{
    File &file = outputFile();
    switch (outputFormat())
    {
        case eHelpOutputFormat_Console:
        {
            TextTableFormatter &formatter(impl_->state_->consoleOptionsFormatter());
            formatter.clear();
            formatter.addColumnLine(0, name);
            formatter.addColumnLine(1, value);
            if (!defaultValue.empty())
            {
                formatter.addColumnLine(2, "(" + defaultValue + ")");
            }
            if (!info.empty())
            {
                formatter.addColumnLine(3, "(" + info + ")");
            }
            TextLineWrapperSettings settings;
            settings.setIndent(11);
            settings.setLineLength(78);
            std::string formattedDescription
                = substituteMarkupAndWrapToString(settings, description);
            file.writeLine(formatter.formatRow());
            file.writeLine(formattedDescription);
            break;
        }
        case eHelpOutputFormat_Rst:
        {
            std::string args(value);
            if (!defaultValue.empty())
            {
                args.append(" (");
                args.append(defaultValue);
                args.append(")");
            }
            if (!info.empty())
            {
                args.append(" (");
                args.append(info);
                args.append(")");
            }
            file.writeLine(formatString("``%s`` %s", name.c_str(), args.c_str()));
            TextLineWrapperSettings settings;
            settings.setIndent(4);
            file.writeLine(substituteMarkupAndWrapToString(settings, description));
            break;
        }
        default:
            GMX_THROW(NotImplementedError(
                              "This output format is not implemented"));
    }
}

void HelpWriterContext::writeOptionListEnd() const
{
}

} // namespace gmx
