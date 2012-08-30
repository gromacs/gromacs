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
/*! \libinternal \file
 * \brief
 * Declares gmx::HelpWriterContext.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
    eHelpOutputFormat_NR        //!< Used for the number of output formats.
};
//! \endcond

/*! \libinternal \brief
 * Context information for writing out help.
 *
 * The purpose of this class is to pass information about the output format to
 * methods that write help, and to abstract away most of the details of
 * different output formats.
 * Additionally, it can keep other context information.  Currently, it keeps
 * the section nesting depth, which is needed for properly formatting titles
 * for nested sections in many output formats.
 * Such additional context information could be useful also for
 * formatting links/references to other help topics.
 *
 * The state of a context object (excluding the fact that the output file is
 * written to) does not change after initial construction of the object.
 * Copy and assign create context objects that share state with the source.
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
        //! Creates a copy of the context.
        HelpWriterContext(const HelpWriterContext &other);
        ~HelpWriterContext();

        //! Assigns a context object.
        HelpWriterContext &operator =(const HelpWriterContext &other);

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
         * \returns   HelpWriterContext for writing the subsection.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Writes \p title using writeTitle() and returns a context object
         * that can be used to write the help for the subsection.
         * The whole subsection should be written out using the returned
         * context before calling any further methods in the parent context.
         *
         * This method is only necessary if the subsection will contain further
         * subsections.  If there is only one level of subsections, it is
         * possible to use writeTitle() directly.
         */
        HelpWriterContext createSubsection(const std::string &title) const;

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
                                        const std::string &text) const;
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
                                        const std::string &text) const;
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
         * \param[in] settins Line wrapper settings.
         * \param[in] text    Text to format.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         *
         * Convenience function that calls substituteMarkupAndWrapToString()
         * and writes the result directly to the output file.
         */
        void writeTextBlock(const TextLineWrapperSettings &settings,
                            const std::string &text) const;

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
};

} // namespace gmx

#endif
