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
/*! \libinternal \file
 * \brief
 * Declares gmx::HelpWriterContext.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
#ifndef GMX_ONLINEHELP_HELPWRITERCONTEXT_H
#define GMX_ONLINEHELP_HELPWRITERCONTEXT_H

#include <string>
#include <vector>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class File;
class TextLineWrapperSettings;

/*! \cond libapi */
//! \libinternal Output format for help writing.
enum HelpOutputFormat
{
    eHelpOutputFormat_Console,  //!< Plain text directly on the console.
    eHelpOutputFormat_Rst,      //!< reStructuredText for online manual and man pages.
    eHelpOutputFormat_Other,    //!< Used for extensions in other modules.
    eHelpOutputFormat_NR        //!< Used for the number of output formats.
};
//! \endcond

/*! \libinternal \brief
 * Hyperlink data for writing out help.
 *
 * This class is separate from HelpWriterContext to allow constructing the list
 * of links once and reusing them across multiple help writer contexts.
 * This is used when exporting all the help from the wrapper binary to avoid
 * repeatedly constructing the same data structure for each help item.
 *
 * While the links are in principle independent of the output format, the
 * constructor takes the output format to be able to preformat the links,
 * avoiding repeated processing during markup substitution.  Could be hidden
 * behind the scenes in HelpWriterContext, but that would complicate the
 * implementation.
 *
 * \ingroup module_onlinehelp
 */
class HelpLinks
{
    public:
        /*! \brief
         * Initializes an empty links collection for the given output format.
         */
        explicit HelpLinks(HelpOutputFormat format);
        ~HelpLinks();

        /*! \brief
         * Adds a link.
         *
         * \param[in] linkName     Name of the link in input text.
         * \param[in] targetName   Hyperlink target.
         * \param[in] displayName  Text to show as the link.
         *
         * Any occurrence of \p linkName in the text passed to markup
         * substitution methods in HelpWriterContext is made into a hyperlink
         * to \p targetName if the markup format supports that.
         */
        void addLink(const std::string &linkName,
                     const std::string &targetName,
                     const std::string &displayName);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        //! Allows the context to use the links.
        friend class HelpWriterContext;
};

/*! \libinternal \brief
 * Context information for writing out help.
 *
 * The purpose of this class is to pass information about the output format to
 * methods that write help, and to abstract away most of the details of
 * different output formats.
 *
 * The state of a context object (excluding the fact that the output file is
 * written to) does not change after initial construction of the object.
 * Copying creates a context objects that share state with the source.
 *
 * \inlibraryapi
 * \ingroup module_onlinehelp
 */
class HelpWriterContext
{
    public:
        /*! \brief
         * Initializes a context with the given output file and format.
         *
         * \throws std::bad_alloc if out of memory.
         */
        HelpWriterContext(File *file, HelpOutputFormat format);
        /*! \brief
         * Initializes a context with the given output file, format and links.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * A reference to \p links is stored until the HelpWriterContext
         * is destructed.  The caller is responsible for ensuring that the
         * links object remains valid long enough.
         */
        HelpWriterContext(File *file, HelpOutputFormat format,
                          const HelpLinks *links);
        //! Creates a copy of the context.
        HelpWriterContext(const HelpWriterContext &other);
        ~HelpWriterContext();

        /*! \brief
         * Adds a string replacement for markup subsitution.
         *
         * \param[in] search   Text to replace in input.
         * \param[in] replace  Text that each occurrence of \p search is
         *     replaced with.
         * \throws std::bad_alloc if out of memory.
         *
         * \todo
         * Improve semantics if the same \p search item is set multiple
         * times.
         */
        void setReplacement(const std::string &search,
                            const std::string &replace);

        /*! \brief
         * Returns the active output format.
         *
         * Does not throw.
         */
        HelpOutputFormat outputFormat() const;
        /*! \brief
         * Returns the raw output file for writing the help.
         *
         * Using this file directly should be avoided, as it requires one to
         * have different code for each output format.
         * Using other methods in this class should be preferred.
         *
         * Does not throw.
         */
        File &outputFile() const;

        /*! \brief
         * Creates a subsection in the output help.
         *
         * \param[in] title  Title for the subsection.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Writes \p title using writeTitle() and makes any further
         * writeTitle() calls write headings one level deeper.
         *
         * Typical use for writing a subsection is to create a copy of the
         * context for the parent section, and then call enterSubSection() on
         * the copy.
         * The whole subsection should be written out using the returned
         * context before calling any further methods in the parent context.
         *
         * This method is only necessary if the subsection will contain further
         * subsections.  If there is only one level of subsections, it is
         * possible to use writeTitle() directly.
         */
        void enterSubSection(const std::string &title);

        /*! \brief
         * Substitutes markup used in help text and wraps lines.
         *
         * \param[in] settings Line wrapper settings.
         * \param[in] text     Text to substitute.
         * \returns   \p text with markup substituted and wrapped.
         * \throws    std::bad_alloc if out of memory.
         *
         * \see TextLineWrapper::wrapToString()
         */
        std::string
        substituteMarkupAndWrapToString(const TextLineWrapperSettings &settings,
                                        const std::string             &text) const;
        /*! \brief
         * Substitutes markup used in help text and wraps lines.
         *
         * \param[in] settings Line wrapper settings.
         * \param[in] text     Text to substitute.
         * \returns   \p text with markup substituted and wrapped as a list of
         *      lines.
         * \throws    std::bad_alloc if out of memory.
         *
         * \see TextLineWrapper::wrapToVector()
         */
        std::vector<std::string>
        substituteMarkupAndWrapToVector(const TextLineWrapperSettings &settings,
                                        const std::string             &text) const;
        /*! \brief
         * Writes a title for the current help topic.
         *
         * \param[in] title  Title to write.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        void writeTitle(const std::string &title) const;
        /*! \brief
         * Writes a formatted text block into the output.
         *
         * \param[in] text  Text to format.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Convenience function that calls substituteMarkupAndWrapToString()
         * and writes the result directly to the output file.
         */
        void writeTextBlock(const std::string &text) const;
        /*! \brief
         * Starts writing a list of options.
         *
         * Prints any necessary headers for a list of options formatted with
         * writeOptionItem().
         */
        void writeOptionListStart() const;
        /*! \brief
         * Writes an entry for a single option into the output.
         *
         * \param[in] name  Name of the option.
         * \param[in] value        Placeholder for option value.
         * \param[in] defaultValue Default value for the option.
         * \param[in] info         Additional (brief) info/attributes for the
         *      option.
         * \param[in] description  Full description of the option.
         */
        void writeOptionItem(
            const std::string &name, const std::string &value,
            const std::string &defaultValue, const std::string &info,
            const std::string &description) const;
        /*! \brief
         * Finishes writing a list of options.
         *
         * Prints any necessary footers for a list of options formatted with
         * writeOptionItem().
         */
        void writeOptionListEnd() const;

    private:
        class Impl;

        /*! \brief
         * Constructs a context object with the given implementation class.
         *
         * \param[in] impl  Implementation object.
         *
         * Does not throw.
         */
        explicit HelpWriterContext(Impl *impl);

        PrivateImplPointer<Impl> impl_;

        GMX_DISALLOW_ASSIGN(HelpWriterContext);
};

} // namespace gmx

#endif
