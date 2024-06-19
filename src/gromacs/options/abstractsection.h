/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
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

#include <memory>
#include <string>

#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

#include "isectionstorage.h"

namespace gmx
{

class IOptionSectionStorage;
class AbstractOption;
class IOptionsContainer;
class OptionInfo;

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
    explicit AbstractOptionSection(const char* name) : name_(name) {}
    virtual ~AbstractOptionSection() {}
    /*! \brief
     * Creates a storage object corresponding to this section.
     *
     * Similar to AbstractOption::createStorage().
     */
    virtual std::unique_ptr<IOptionSectionStorage> createStorage() const = 0;
    //! \endcond

private:
    const char* name_;

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
    IOptionsContainer& addGroup() override;

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
    template<typename StorageType>
    static StorageType* getStorage(internal::OptionSectionImpl* section)
    {
        IOptionSectionStorage* storage      = getStorage(section);
        StorageType*           typedStorage = dynamic_cast<StorageType*>(storage);
        GMX_ASSERT(typedStorage != nullptr, "Mismatching section storage type");
        return typedStorage;
    }

    //! Wraps a given section storage object.
    explicit AbstractOptionSectionHandle(internal::OptionSectionImpl* section) : section_(section)
    {
    }
    //! \endcond

private:
    // From IOptionsContainerWithSections
    internal::OptionSectionImpl* addSectionImpl(const AbstractOptionSection& section) override;
    // From IOptionsContainer
    OptionInfo* addOptionImpl(const AbstractOption& settings) override;

    /*! \brief
     * Implementation helper for the template method.
     *
     * This allows encapsulating the implementation within the source file.
     */
    static IOptionSectionStorage* getStorage(internal::OptionSectionImpl* section);

    internal::OptionSectionImpl* section_;
};

class AbstractOptionSectionInfo
{
public:
    //! Wraps a given section storage object.
    explicit AbstractOptionSectionInfo(internal::OptionSectionImpl* section) : section_(*section) {}

    //! Returns the name of the section.
    const std::string& name() const;

    //! Returns the wrapped section storage object.
    internal::OptionSectionImpl& section() { return section_; }
    //! Returns the wrapped section storage object.
    const internal::OptionSectionImpl& section() const { return section_; }

private:
    internal::OptionSectionImpl& section_;

    GMX_DISALLOW_COPY_AND_ASSIGN(AbstractOptionSectionInfo);
};

} // namespace gmx

#endif
