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
 * Declares gmx::OptionsVisitor interface and supporting classes.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSVISITOR_H
#define GMX_OPTIONS_OPTIONSVISITOR_H

#include <string>

namespace gmx
{

class AbstractOptionStorage;
class Options;

/*! \libinternal \brief
 * Wrapper class for accessing option information.
 *
 * This class isolates the details of the internal option implementation
 * from option visitors.
 *
 * \see OptionsVisitor
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionInfo
{
    public:
        /*! \brief
         * Wraps a given option object.
         */
        OptionInfo(const AbstractOptionStorage &option);

        //! Returns true if the option is a boolean option.
        bool isBoolean() const;
        //! Returns true if the option is a file name option.
        bool isFile() const;
        //! Returns true if the option is a hidden option.
        bool isHidden() const;
        //! Returns the name of the option.
        const std::string &name() const;
        //! Returns the description of the option.
        const std::string &description() const;
        //! Returns the type of the option as a string.
        const char *type() const;
        //! Returns the number of values given for the option.
        int valueCount() const;
        //! Returns the i'th value of the option as a string.
        std::string formatValue(int i) const;
        //! Returns all the values of the option as a single string.
        std::string formatValues() const;

    private:
        //! The wrapped option.
        const AbstractOptionStorage &_option;

        // Disallow copy and assign.
        OptionInfo(const OptionInfo &);
        void operator =(const OptionInfo &);
};

/*! \libinternal \brief
 * Pure interface for visiting options in a Options object.
 *
 * \see OptionsIterator
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsVisitor
{
    public:
        virtual ~OptionsVisitor() {}

        /*! \brief
         * Called for each subsection in Options.
         */
        virtual void visitSubSection(const Options &section) = 0;
        /*! \brief
         * Called for each option in Options.
         */
        virtual void visitOption(const OptionInfo &option) = 0;
};

/*! \libinternal \brief
 * Decorator class for visiting options in a Options object.
 *
 * This class provides an interface for looping through subsections and
 * options in a Options object.
 *
 * Typical use (loop over all options, iteratively descending into
 * subsections):
 * \code
class Visitor : public gmx::options::OptionsVisitor
{
    public:
        void visitSubSection(const Options &section)
        {
            OptionsIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOption(const OptionInfo &option)
        {
            // Do something.
        }
}

Visitor().visitSubSection(options);
 * \endcode
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsIterator
{
    public:
        /*! \brief
         * Creates an object for visiting options in a Options object.
         */
        OptionsIterator(const Options &options);

        /*! \brief
         * Visits each subsection in the wrapped Options object.
         */
        void acceptSubSections(OptionsVisitor *visitor) const;
        /*! \brief
         * Visits each option in the wrapped Options object.
         */
        void acceptOptions(OptionsVisitor *visitor) const;

    private:
        //! The wrapped Options object.
        const Options          &_options;

        // Disallow copy and assign.
        OptionsIterator(const OptionsIterator &);
        void operator =(const OptionsIterator &);
};

} // namespace gmx

#endif
