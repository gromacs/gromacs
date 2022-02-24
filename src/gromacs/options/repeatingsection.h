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
 * Declares gmx::RepeatingOptionSection and related classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_REPEATINGSECTION_H
#define GMX_OPTIONS_REPEATINGSECTION_H

#include <memory>
#include <vector>

#include "gromacs/options/ioptionscontainerwithsections.h"

#include "abstractsection.h"
#include "isectionstorage.h"
#include "valuestore.h"

namespace gmx
{

template<class T>
class RepeatingOptionSectionHandle;
template<class T>
class RepeatingOptionSectionStorage;

/*! \brief
 * Declares an option section that creates a structure for each instance.
 *
 * \tparam  T  Structure that stores the values for an instance of the section.
 *
 * This class declares a section that can be specified multiple times, and
 * creates an instance of `T` for each instance.  Options within the section
 * can store their values in the created structure.
 *
 * \inpublicapi
 * \ingroup module_options
 */
template<class T>
class RepeatingOptionSection : public AbstractOptionSection
{
public:
    //! AbstractOptionSectionHandle corresponding to this option type.
    typedef RepeatingOptionSectionHandle<T> HandleType;

    //! Creates a section with the given name.
    explicit RepeatingOptionSection(const char* name) :
        AbstractOptionSection(name), values_(nullptr)
    {
    }

    //! Specifies a vector to receive the section structures.
    RepeatingOptionSection& storeVector(std::vector<T>* values)
    {
        values_ = values;
        return *this;
    }

private:
    std::unique_ptr<IOptionSectionStorage> createStorage() const override;

    std::vector<T>* values_;

    //! Allows accessing the properties when initializing the storage.
    friend class RepeatingOptionSectionStorage<T>;
};

/*! \internal
 * \brief
 * Implements handling of the structures that stores per-section values.
 *
 * \ingroup module_options
 */
template<class T>
class RepeatingOptionSectionStorage : public IOptionSectionStorage
{
public:
    //! Initializes the storage for given section properties.
    explicit RepeatingOptionSectionStorage(const RepeatingOptionSection<T>& section) :
        store_(new OptionValueStoreVector<T>(section.values_)), currentData_()
    {
    }

    void initStorage() override { defaultValues_ = currentData_; }
    void startSection() override { resetSection(); }
    void finishSection() override
    {
        store_->append(currentData_);
        resetSection();
    }

private:
    void resetSection() { currentData_ = defaultValues_; }

    //! Final storage location for the section structures.
    const std::unique_ptr<IOptionValueStore<T>> store_;
    T                                           defaultValues_;
    /*! \brief
     * Stores the values for the current in-process section.
     *
     * Options within the section store their values to this structure, and
     * they are then copied to the final storage when the section finishes.
     */
    T currentData_;

    //! Allows binding option storage to the current section data structure.
    friend class RepeatingOptionSectionHandle<T>;
};

template<class T>
std::unique_ptr<IOptionSectionStorage> RepeatingOptionSection<T>::createStorage() const
{

    return std::make_unique<RepeatingOptionSectionStorage<T>>(*this);
}

/*! \brief
 * Allows adding options to an RepeatingOptionSection.
 *
 * An instance of this class is returned from
 * IOptionsContainerWithSections::addSection(), and supports adding options and
 * subsections to a section created with OptionSection.
 *
 * Example:
 * \code
   struct SectionData { int value; }
   // as class attribute
   std::vector<SectionData> values;

   void MyClass::initOptions(gmx::IOptionsContainerWithSections *options)
   {
       auto sec = options->addSection(gmx::RepeatingOptionSection<SectionData>("sec").storeVector(&values));
       sec->addOption(gmx::IntegerOption("arg").store(&sec.bind().value));
   }
   \endcode
 *
 * \inpublicapi
 * \ingroup module_options
 */
template<class T>
class RepeatingOptionSectionHandle : public AbstractOptionSectionHandle
{
public:
    //! Wraps a given section storage object.
    explicit RepeatingOptionSectionHandle(internal::OptionSectionImpl* section) :
        AbstractOptionSectionHandle(section),
        storage_(getStorage<RepeatingOptionSectionStorage<T>>(section))
    {
    }

    /*! \brief
     * Supports storing option values within the per-section data structure.
     *
     * See class documentation for an example.
     */
    T& bind() { return storage_->currentData_; }

private:
    RepeatingOptionSectionStorage<T>* storage_;
};

} // namespace gmx

#endif
