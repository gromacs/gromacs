/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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

#include "../utility/common.h"

namespace gmx
{

class File;
class TextLineWrapperSettings;

/*! \cond libapi */
//! \libinternal Output format for help writing.
enum HelpOutputFormat
{
    eHelpOutputFormat_Console,  //!< Plain text directly on the console.
    eHelpOutputFormat_Man,      //!< Man page.
    eHelpOutputFormat_Html,     //!< Html output for online manual.
    eHelpOutputFormat_NR        //!< Used for the number of output formats.
};
//! \endcond

/*! \libinternal \brief
 * Context information for writing out help.
 *
 * The purpose of this class is to pass information about the output format to
 * methods that write help, and to abstract away most of the details of
 * different output formats.
 * Additionally, it can keep other context information, although it currently
 * does not.  Such additional context information would be useful for
 * formatting links/references to other help topics.
 *
 * TODO: This class will need additional work as part of Redmine issue #969.
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
        ~HelpWriterContext();

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
         * Writes a formatted text block into the output.
         *
         * \param[in] settings Line wrapper settings.
         * \param[in] text     Text to format.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Convenience function that calls substituteMarkupAndWrapToString()
         * and writes the result directly to the output file.
         */
        void writeTextBlock(const TextLineWrapperSettings &settings,
                            const std::string             &text) const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
