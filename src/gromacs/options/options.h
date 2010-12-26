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
/*! \file
 * \brief
 * Declares gmx::Options.
 *
 * Together with basicoptions.h, this header forms the part of the public
 * API that most classes will use to provide options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONS_H
#define GMX_OPTIONS_OPTIONS_H

#include <string>

namespace gmx
{

class AbstractErrorReporter;

class AbstractOption;
class OptionsGlobalProperties;
class OptionsAssigner;
class OptionsIterator;

/*! \brief
 * Collection of options.
 *
 * This class provides a standard interface for implementing input options.
 * Standard usage is to write a method that creates an Options that is owned by
 * the object, populates it with supported options, and then returns it:
 * \code
// <as class attributes>
using gmx::Options;
Options      options("common", "Common Options");
std::string  arg1;
int          arg2;

// <populating>
using gmx::StringOption;
using gmx::IntegerOption;
options.addOption(StringOption("arg1").store(&arg1));
options.addOption(IntegerOption("arg2").store(&arg2));
return &options;
 * \endcode
 * The caller of that method can then use a parser implementation such as
 * CommandLineParser to provide values for the options.
 *
 * Header basicoptions.h provides declarations of several standard
 * option types for use with addOption().  Documentation of those classes
 * also give more examples of how to define options.
 *
 * In order to keep the public interface of this class simple and to reduce
 * build dependencies on objects that simply provide options, functionality
 * to assign values to options is provided by a separate OptionsAssigner class.
 * Similarly, functionality for looping over all options (e.g., for writing out
 * help) is provided by OptionsIterator.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class Options
{
    public:
        /*! \brief
         * Initializes the name and title of an option collection.
         *
         * \param[in] name  Single-word name.
         * \param[in] title Descriptive title.
         *
         * Copies the input strings.
         */
        Options(const char *name, const char *title);
        ~Options();

        //! Returns the short name of the option collection.
        const std::string &name() const;
        //! Returns the title of the option collection.
        const std::string &title() const;
        //! Returns the full description of the option collection.
        const std::string &description() const;

        //! Sets the full description of the option collection.
        void setDescription(const char *const *desc);
        //int addBugs(int nbugs, const char *const *bugs);

        /*! \brief
         * Adds an option collection as a subsection of this collection.
         *
         * The name() field of \p section is used as the name of the
         * subsection.
         *
         * For certain functionality to work properly, no options should
         * be added to the subsection after it has been added to another
         * collection.
         *
         * Only a pointer to the provided object is stored.  The caller is
         * responsible that the object exists for the lifetime of the
         * collection.
         * It is not possible to add the same Options object as a subsection to
         * several different Options.
         */
        void addSubSection(Options *section);
        /*! \brief
         * Adds a recognized option to the collection.
         *
         * The behavior is undefined if invalid settings are provided.
         * The current implementation asserts if it detects an error.
         *
         * See \link Options class documentation \endlink for example usage.
         */
        void addOption(const AbstractOption &settings);
        /*! \brief
         * Adds default options to this collection.
         *
         * Adds default options for altering global properties such as the time
         * unit into this collection.
         *
         * It is possible to call this method on a subsection of a collection.
         * Even in that case, the generated options set the properties for the
         * parent collection.
         */
        void addDefaultOptions();

        //! Returns true if option is set.
        bool isSet(const char *name) const;
        /*! \brief
         * Notifies the collection that all option values are assigned.
         *
         * \param[in] errors  Error reporter for reporting errors.
         * \retval 0 on success.
         *
         * This function should be called after no more option values are
         * to be assigned.  Values in storage variables are guaranteed to be
         * available only after this call.
         */
        int finish(AbstractErrorReporter *errors);

        /*! \brief
         * Returns the global property object for this collection.
         *
         * The caller should not store pointers or references to the object;
         * it can change if this object is added as a subsection into
         * another collection.
         *
         * \see OptionsGlobalProperties
         */
        OptionsGlobalProperties &globalProperties();
        //! \copydoc globalProperties()
        const OptionsGlobalProperties &globalProperties() const
        { return const_cast<Options *>(this)->globalProperties(); }

    private:
        class Impl;

        Impl                   *_impl;

        friend class Impl;

        //! Needed to be able to extend the interface of this object.
        friend class OptionsAssigner;
        //! Needed to be able to extend the interface of this object.
        friend class OptionsIterator;

        // Disallow copy and assign.
        Options(const Options &);
        void operator =(const Options &);
};

} // namespace gmx

#endif
