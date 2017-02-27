/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2014,2016,2017, by the GROMACS development team, led by
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
 * Declares gmx::OptionsVisitor interface and supporting classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSVISITOR_H
#define GMX_OPTIONS_OPTIONSVISITOR_H

#include <cstddef>

#include <string>

#include "gromacs/options/abstractoption.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class Options;
class OptionSectionInfo;

namespace internal
{
class OptionSectionImpl;
}

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

        //! Called for each section.
        virtual void visitSection(const OptionSectionInfo &section) = 0;
        //! Called for each option.
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

        virtual void visitSection(const OptionSectionInfo &section) = 0;
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
 * This class provides an interface for looping through sections and
 * options in a Options object.
 *
 * Typical use (loop over all options, iteratively descending into
 * subsections):
 * \code
   class Visitor : public gmx::OptionsVisitor
   {
       public:
           virtual void visitSection(const OptionSectionInfo &section)
           {
               OptionsIterator iterator(section);
               iterator.acceptSections(this);
               iterator.acceptOptions(this);
           }

           virtual void visitOption(const OptionInfo &option)
           {
               // Do something.
           }
   }

   Visitor().visitSection(options.rootSection());
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
         *
         * The created iterator iterates over the "root" section in the Options
         * object.
         */
        explicit OptionsIterator(const Options &options);
        //! Creates an object for visiting options in an options section.
        explicit OptionsIterator(const OptionSectionInfo &section);

        //! Visits each section in the wrapped section.
        void acceptSections(OptionsVisitor *visitor) const;
        //! Visits each option in the wrapped section.
        void acceptOptions(OptionsVisitor *visitor) const;

    private:
        //! The wrapped section object.
        const internal::OptionSectionImpl &section_;

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

        //! Called for each section.
        virtual void visitSection(OptionSectionInfo *section) = 0;
        //! Called for each option.
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

        virtual void visitSection(OptionSectionInfo *section) = 0;
        /*! \brief
         * Called for each option of type \p InfoType.
         */
        virtual void visitOptionType(InfoType *option) = 0;

    private:
        virtual void visitOption(OptionInfo *option)
        {
            InfoType *subtype = option->toType<InfoType>();
            if (subtype != nullptr)
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
         *
         * The created iterator iterates over the "root" section in the Options
         * object.
         */
        explicit OptionsModifyingIterator(Options *options);
        //! Creates an object for visiting options in an options section.
        explicit OptionsModifyingIterator(OptionSectionInfo *section);

        //! Visits each section in the wrapped section.
        void acceptSections(OptionsModifyingVisitor *visitor) const;
        //! Visits each option in the wrapped section.
        void acceptOptions(OptionsModifyingVisitor *visitor) const;

    private:
        //! The wrapped section object.
        internal::OptionSectionImpl &section_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsModifyingIterator);
};

} // namespace gmx

#endif
