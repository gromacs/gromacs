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
 * Declares gmx::OptionsAssigner.
 *
 * This header is only needed when implementing option parsers.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSASSIGNER_H
#define GMX_OPTIONS_OPTIONSASSIGNER_H

#include <string>

namespace gmx
{

class AbstractErrorReporter;

class Options;

/*! \libinternal \brief
 * Decorator class for assigning values to Options.
 *
 * This class extends the interface of an Options object by providing methods
 * to set values for options.  It also keeps track of necessary state variables
 * to assign values to options in subsections within the Options object.
 * Typical use (without error checking):
 * \code
gmx::options::Options options("name", "Title");
// Set up options

gmx::error::StandardReporter errors;
gmx::options::OptionsAssigner assigner(&options, &errors);
assigner.startOption("opt1");
assigner.appendValue("3");
assigner.startSubSection("section");
assigner.startOption("opt2"); // Now in the subsection
assigner.appendValue("yes");
assigner.finishSubSection()
assigner.startOption("opt3"); // Again in the main options
assigner.appendValue("2");
assigner.finish(); // At minimum, the return value of finish() should be checked.
 * \endcode
 *
 * As shown in the example, calling finishOption() or finishSubSection() is
 * optional; they are automatically called when appropriate by startOption(),
 * startSubSection(), and finish().
 * However, you need to call them explicitly if you want to act on the return
 * value: these calls do not influence the return value of
 * startOption() / startSubSection().
 * They do influence the return value of finish(), however.
 * The finish() method should always be called.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsAssigner
{
    public:
        /*! \brief
         * Creates an object that assigns to the given object.
         */
        OptionsAssigner(Options *options, AbstractErrorReporter *errors);
        ~OptionsAssigner();

        /*! \brief
         * Returns the error reporter object passed to the constructor.
         *
         * This method is provided for convenience such that users of this
         * class do not need to store a separate pointer to the error reporter
         * if they need to use it.
         */
        AbstractErrorReporter *errorReporter() const;

        /*! \brief
         * Sets the assigner to recognize boolean options with a "no" prefix.
         *
         * With this option set, \c startOption("noname") is interpreted as
         * \c startOption("name") followed by \c appendValue("no"), if there is
         * no option by the name "noname", but there is a boolean option with
         * name "name".
         *
         * By default, the prefix is not recognized.
         *
         * Can be set or cleared at any time, and will have effect on all
         * subsequent calls of startOption().
         */
        void setAcceptBooleanNoPrefix(bool enabled);
        /*! \brief
         * Sets the assigner to find options in non-active sections.
         *
         * By default, options are only looked for in the currently active
         * subsection.  With this option set, if no matching option is found in
         * the current section, a breadth-first search is performed, first on
         * all subsections of the current section, and then going up one level
         * at a time.  The first matching option is used, and the current
         * section is changed to the section that contains the matching option.
         *
         * Can be set or cleared at any time, and will have effect on all
         * subsequent calls of startOption().
         */
        void setNoStrictSectioning(bool enabled);

        /*! \brief
         * Start assigning values.
         *
         * \retval 0 on success.
         */
        int start();
        /*! \brief
         * Start assigning values to options in a subsection.
         *
         * \param[in] name  Name of the subsection to start assigning to.
         * \retval 0 if assignment can proceed.
         *
         * Does not call finishSubSection() automatically to enable nested
         * sections.
         */
        int startSubSection(const char *name);
        /*! \brief
         * Start assigning values for an option.
         *
         * \param[in] name  Name of the option to start assigning to.
         * \retval 0 if assignment can proceed.
         */
        int startOption(const char *name);
        /*! \brief
         * Appends a value to the value list of the current option.
         *
         * \param[in] value  String representation of the value to assign.
         * \retval 0 if assignment was successful.
         */
        int appendValue(const std::string &value);
        /*! \brief
         * Finish assigning values for the current option.
         *
         * \retval 0 if there were no errors in the assignment.
         *
         * This function returns non-zero only if the error could not have been
         * detected earlier, i.e., from the return value of appendValue().
         */
        int finishOption();
        /*! \brief
         * Finish assigning values to a subsection.
         *
         * \retval 0 for success.
         *
         * This function returns non-zero only if the error could not have been
         * detected earlier.
         */
        int finishSubSection();
        /*! \brief
         * Finish assigning options through the object.
         *
         * \retval 0 if there were no errors in the assignment.
         *
         * If an error was detected in any of the other calls to this class,
         * this function returns the error code of the first of such errors.
         */
        int finish();

    private:
        class Impl;

        Impl                   *_impl;

        // Disallow copy and assign.
        OptionsAssigner(const OptionsAssigner &);
        void operator =(const OptionsAssigner &);
};

} // namespace gmx

#endif
