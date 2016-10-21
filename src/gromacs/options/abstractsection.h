/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares base classes for declaring option sections.
 *
 * This header defines base classes for option section settings that are used
 * with IOptionsContainerWithSections::addSection().  These classes implement
 * the "named parameter" idiom for specifying section properties.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ABSTRACTSECTION_H
#define GMX_OPTIONS_ABSTRACTSECTION_H

#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/isectionstorage.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class IOptionSectionStorage;

namespace internal
{
class OptionSectionImpl;
}

/*! \brief
 * Base class for specifying option section properties.
 *
 * \ingroup module_options
 */
class AbstractOptionSection
{
    protected:
        //! \cond libapi
        //! Initializes option properties with the given name.
        explicit AbstractOptionSection(const char *name) : name_(name) {}

        /*! \brief
         * Creates a storage object corresponding to this section.
         *
         * Similar to AbstractOption::createStorage().
         */
        virtual IOptionSectionStorage *createStorage() const = 0;
        //! \endcond

    private:
        const char *name_;

        friend class internal::OptionSectionImpl;
};

/*! \brief
 * Base class for handles to option sections.
 *
 * This class implements the common functionality for adding options and
 * subsections to option sections.
 *
 * \ingroup module_options
 */
class AbstractOptionSectionHandle : public IOptionsContainerWithSections
{
    public:
        // From IOptionsContainer
        //! \copydoc IOptionsContainer::addGroup()
        virtual IOptionsContainer &addGroup();

    protected:
        //! \cond libapi
        /*! \brief
         * Returns the storage for a particular type of section.
         *
         * This is intended for use in derived class constructors, where the
         * handle needs access to the actual storage.  The handle should know
         * the type of storage created for the section type it deals with, so
         * the cast should always be successful.
         */
        template <typename StorageType>
        static StorageType *getStorage(internal::OptionSectionImpl *section)
        {
            IOptionSectionStorage *storage = getStorage(section);
            StorageType           *typedStorage
                = dynamic_cast<StorageType *>(storage);
            GMX_ASSERT(typedStorage != nullptr, "Mismatching section storage type");
            return typedStorage;
        }

        //! Wraps a given section storage object.
        explicit AbstractOptionSectionHandle(internal::OptionSectionImpl *section)
            : section_(section)
        {
        }
        //! \endcond

    private:
        // From IOptionsContainerWithSections
        virtual internal::OptionSectionImpl *
        addSectionImpl(const AbstractOptionSection &section);
        // From IOptionsContainer
        virtual OptionInfo *addOptionImpl(const AbstractOption &settings);

        /*! \brief
         * Implementation helper for the template method.
         *
         * This allows encapsulating the implementation within the source file.
         */
        static IOptionSectionStorage *getStorage(internal::OptionSectionImpl *section);

        internal::OptionSectionImpl *section_;
};

class AbstractOptionSectionInfo
{
    public:
        //! Wraps a given section storage object.
        explicit AbstractOptionSectionInfo(internal::OptionSectionImpl *section)
            : section_(*section)
        {
        }

        //! Returns the name of the section.
        const std::string &name() const;

        //! Returns the wrapped section storage object.
        internal::OptionSectionImpl       &section() { return section_; }
        //! Returns the wrapped section storage object.
        const internal::OptionSectionImpl &section() const { return section_; }

    private:
        internal::OptionSectionImpl &section_;

        GMX_DISALLOW_COPY_AND_ASSIGN(AbstractOptionSectionInfo);
};

} // namespace gmx

#endif
