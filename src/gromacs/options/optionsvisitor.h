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

#include <cstddef>

#include <string>

#include "../utility/common.h"

#include "abstractoption.h"

namespace gmx
{

class Options;

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
 * Abstract base class for visiting options of a particular type.
 *
 * \see OptionsIterator
 * \see OptionsVisitor
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template <class InfoType>
class OptionsTypeVisitor : public OptionsVisitor
{
    public:
        virtual ~OptionsTypeVisitor() {}

        virtual void visitSubSection(const Options &section) = 0;
        /*! \brief
         * Called for each option of type \p InfoType.
         */
        virtual void visitOptionType(const InfoType &option) = 0;

    private:
        virtual void visitOption(const OptionInfo &option)
        {
            const InfoType *subtype = option.toType<InfoType>();
            if (subtype != NULL)
            {
                visitOptionType(*subtype);
            }
        }
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
   class Visitor : public gmx::OptionsVisitor
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
        explicit OptionsIterator(const Options &options);

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
        const Options          &options_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsIterator);
};

/*! \libinternal \brief
 * Pure interface for visiting options in a Options object, allowing
 * modifications.
 *
 * \see OptionsModifyingIterator
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsModifyingVisitor
{
    public:
        virtual ~OptionsModifyingVisitor() {}

        /*! \brief
         * Called for each subsection in Options.
         */
        virtual void visitSubSection(Options *section) = 0;
        /*! \brief
         * Called for each option in Options.
         */
        virtual void visitOption(OptionInfo *option) = 0;
};

/*! \libinternal \brief
 * Abstract base class for visiting options of a particular type, allowing
 * modifications.
 *
 * \see OptionsModifyingIterator
 * \see OptionsModifyingVisitor
 *
 * \inlibraryapi
 * \ingroup module_options
 */
template <class InfoType>
class OptionsModifyingTypeVisitor : public OptionsModifyingVisitor
{
    public:
        virtual ~OptionsModifyingTypeVisitor() {}

        virtual void visitSubSection(Options *section) = 0;
        /*! \brief
         * Called for each option of type \p InfoType.
         */
        virtual void visitOptionType(InfoType *option) = 0;

    private:
        virtual void visitOption(OptionInfo *option)
        {
            InfoType *subtype = option->toType<InfoType>();
            if (subtype != NULL)
            {
                visitOptionType(subtype);
            }
        }
};

/*! \libinternal \brief
 * Decorator class for visiting options in a Options object, allowing changes.
 *
 * This class works exactly like OptionsIterator, except that it uses
 * OptionsModifyingVisitor interface, which allows modifying the options.
 *
 * \see OptionsIterator
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsModifyingIterator
{
    public:
        /*! \brief
         * Creates an object for visiting options in a Options object.
         */
        explicit OptionsModifyingIterator(Options *options);

        /*! \brief
         * Visits each subsection in the wrapped Options object.
         */
        void acceptSubSections(OptionsModifyingVisitor *visitor) const;
        /*! \brief
         * Visits each option in the wrapped Options object.
         */
        void acceptOptions(OptionsModifyingVisitor *visitor) const;

    private:
        //! The wrapped Options object.
        Options                &options_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsModifyingIterator);
};

} // namespace gmx

#endif
