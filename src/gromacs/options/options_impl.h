/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Declares private implementation class for gmx::Options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONS_IMPL_H
#define GMX_OPTIONS_OPTIONS_IMPL_H

#include <list>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/classhelpers.h"

#include "isectionstorage.h"

namespace gmx
{

namespace internal
{

/*! \internal
 * \brief
 * Internal implementation class for storing an option section.
 *
 * All options are stored within a section: the top-level contents of an
 * Options object are handled within an unnamed, "root" section.
 * This class handles the common functionality for all sections, related to
 * storing the options and subsections.  Functionality specific to a section
 * type is provided by IOptionSectionStorage.
 *
 * \ingroup module_options
 */
class OptionSectionImpl : public IOptionsContainerWithSections
{
public:
    /*! \internal \brief
     * Describes a group of options (see Options::addGroup()).
     *
     * \ingroup module_options
     */
    class Group : public IOptionsContainer
    {
    public:
        //! Convenience typedef for list of options.
        typedef std::vector<AbstractOptionStorage*> OptionList;
        //! Convenience typedef for list of subgroups.
        typedef std::list<Group> SubgroupList;

        //! Creates a group within the given Options.
        explicit Group(OptionSectionImpl* parent) : parent_(parent) {}

        // From IOptionsContainer
        IOptionsContainer& addGroup() override;
        OptionInfo*        addOptionImpl(const AbstractOption& settings) override;

        //! Containing options object.
        OptionSectionImpl* parent_;
        /*! \brief
         * List of options, in insertion order.
         *
         * Pointers in this container point to the objects managed by
         * Impl::optionsMap_.
         */
        OptionList options_;
        //! List of groups, in insertion order.
        SubgroupList subgroups_;
    };

    //! Smart pointer for managing an AbstractOptionStorage object.
    typedef std::unique_ptr<AbstractOptionStorage> AbstractOptionStoragePointer;
    //! Convenience typedef for a map that contains all the options.
    typedef std::map<std::string, AbstractOptionStoragePointer> OptionMap;
    //! Smart pointer for managing subsections.
    typedef std::unique_ptr<OptionSectionImpl> SectionPointer;
    //! Convenience typedef for a container for subsections.
    typedef std::vector<SectionPointer> SectionList;

    //! Creates storage for a new section.
    OptionSectionImpl(const OptionManagerContainer&          managers,
                      std::unique_ptr<IOptionSectionStorage> storage,
                      const char*                            name) :
        managers_(managers),
        storage_(std::move(storage)),
        info_(this),
        name_(name),
        rootGroup_(this),
        storageInitialized_(false)
    {
    }

    // From IOptionsContainerWithSections
    OptionSectionImpl* addSectionImpl(const AbstractOptionSection& section) override;

    // From IOptionsContainer
    IOptionsContainer& addGroup() override;
    OptionInfo*        addOptionImpl(const AbstractOption& settings) override;

    //! Returns section info object for this section.
    OptionSectionInfo& info() { return info_; }
    //! Returns section info object for this section.
    const OptionSectionInfo& info() const { return info_; }

    /*! \brief
     * Finds a subsection by name.
     *
     * \param[in] name  Name to search for.
     * \returns Pointer to the found subsection, or NULL if not found.
     *
     * Does not throw.
     */
    OptionSectionImpl* findSection(const char* name) const;
    /*! \brief
     * Finds an option by name.
     *
     * \param[in] name  Name to search for.
     * \returns Pointer to the found option, or NULL if not found.
     *
     * Does not throw.
     */
    AbstractOptionStorage* findOption(const char* name) const;

    /*! \brief
     * Called when entering the section.
     *
     * Calls AbstractOptionStorage::startSource() for all options.
     */
    void start();
    /*! \brief
     * Calls AbstractOptionStorage::finish() for all options.
     */
    void finish();

    //! Reference to the option managers in the parent Options object.
    const OptionManagerContainer& managers_;
    //! Type-specific storage object for this section.
    std::unique_ptr<IOptionSectionStorage> storage_;
    //! Info object for this section.
    OptionSectionInfo info_;
    //! Name of this section (empty and unused for the root section).
    std::string name_;
    /*! \brief
     * Group that contains all options (and subgroups).
     *
     * This is used to store the insertion order of options.
     */
    Group rootGroup_;
    //! Map from option names to options; owns the option storage objects.
    OptionMap optionMap_;
    //! List of subsections, in insertion order.
    SectionList subsections_;
    //! Whether initStorage() has been called for `storage_`.
    bool storageInitialized_;

    GMX_DISALLOW_COPY_AND_ASSIGN(OptionSectionImpl);
};

/*! \internal
 * \brief
 * Private implementation class for Options.
 *
 * Note that in addition to Options, the OptionsAssigner class also directly
 * accesses this class.
 *
 * \ingroup module_options
 */
class OptionsImpl
{
public:
    OptionsImpl();

    //! Option managers set for this collection.
    OptionManagerContainer managers_;
    //! Root section for this collection.
    OptionSectionImpl rootSection_;
};

} // namespace internal

} // namespace gmx

#endif
