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
/*! \file
 * \brief
 * Declares gmx::Options.
 *
 * Together with basicoptions.h, this header forms the part of the public
 * API that most classes will use to provide options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONS_H
#define GMX_OPTIONS_OPTIONS_H

#include <memory>
#include <string>

#include "gromacs/options/ioptionscontainerwithsections.h"

namespace gmx
{

class AbstractOption;
class OptionSection;
class OptionSectionInfo;
class OptionsAssigner;

namespace internal
{
class OptionsImpl;
}

/*! \brief
 * Base class for option managers.
 *
 * This class is used as a marker for all classes that are used with
 * Options::addManager().  It doesn't provide any methods, but only supports
 * transporting these classes through the Options collection into the
 * individual option implementation classes.
 *
 * The virtual destructor is present to make this class polymorphic, such that
 * `dynamic_cast` can be used when retrieving a manager of a certain type for
 * the individual options.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class IOptionManager
{
protected:
    virtual ~IOptionManager();
};

/*! \brief
 * Collection of options.
 *
 * See \ref module_options for an overview of how the options work.
 * The IOptionsContainerWithSections interface documents how to add options.
 *
 * In order to keep the public interface of this class simple, functionality
 * to assign values to options is provided by a separate OptionsAssigner class.
 * Similarly, functionality for looping over all options (e.g., for writing out
 * help) is provided by OptionsIterator.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class Options : public IOptionsContainerWithSections
{
public:
    //! Initializes an empty options root container.
    Options();
    ~Options() override;

    /*! \brief
     * Adds an option manager.
     *
     * \param    manager Manager to add.
     * \throws   std::bad_alloc if out of memory.
     *
     * Option managers are used by some types of options that require
     * interaction between different option instances (e.g., selection
     * options), or need to support globally set properties (e.g., a global
     * default file prefix).  Option objects can retrieve the pointer to
     * their manager when they are created, and the caller can alter the
     * behavior of the options through the manager.
     * See the individual managers for details.
     *
     * Caller is responsible for memory management of \p manager.
     * The Options object (and its contained options) only stores a
     * reference to the object.
     *
     * This method cannot be called after adding options or sections.
     */
    void addManager(IOptionManager* manager);

    // From IOptionsContainer
    IOptionsContainer& addGroup() override;

    //! Returns a handle to the root section.
    OptionSectionInfo& rootSection();
    //! Returns a handle to the root section.
    const OptionSectionInfo& rootSection() const;

    /*! \brief
     * Notifies the collection that all option values are assigned.
     *
     * \throws InvalidInputError if invalid user input is detected.
     *
     * This function should be called after no more option values are
     * to be assigned.  Values in storage variables are guaranteed to be
     * available only after this call, although in most cases, they are
     * available already during assignment.
     *
     * If invalid option values, e.g., missing required option, is detected
     * at this point, this function throws.  The thrown exception contains
     * information on all errors detected during the call.
     */
    void finish();

private:
    // From IOptionsContainerWithSections
    internal::OptionSectionImpl* addSectionImpl(const AbstractOptionSection& section) override;
    // From IOptionsContainer
    OptionInfo* addOptionImpl(const AbstractOption& settings) override;

    std::unique_ptr<internal::OptionsImpl> impl_;

    //! Needed to be able to extend the interface of this object.
    friend class OptionsAssigner;
};

} // namespace gmx

#endif
