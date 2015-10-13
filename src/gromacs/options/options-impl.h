/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/options/options.h"

namespace gmx
{

class AbstractOptionStorage;

namespace internal
{

/*! \internal
 * \brief
 * Private implementation class for Options.
 *
 * Note that in addition to Options, the OptionsAssigner and OptionsIterator
 * classes also directly access this class.
 *
 * \ingroup module_options
 */
class OptionsImpl
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
                typedef std::vector<AbstractOptionStorage *> OptionList;
                //! Convenience typedef for list of subgroups.
                typedef std::list<Group> SubgroupList;

                //! Creates a group within the given Options.
                explicit Group(OptionsImpl *parent) : parent_(parent) {}

                // From IOptionsContainer
                virtual IOptionsContainer &addGroup();
                virtual OptionInfo *addOption(const AbstractOption &settings);

                //! Containing options object.
                OptionsImpl  *parent_;
                /*! \brief
                 * List of options, in insertion order.
                 *
                 * Pointers in this container point to the objects managed by
                 * Impl::optionsMap_.
                 */
                OptionList    options_;
                //! List of groups, in insertion order.
                SubgroupList  subgroups_;
        };

        //! Smart pointer for managing an AbstractOptionStorage object.
        typedef std::unique_ptr<AbstractOptionStorage>
            AbstractOptionStoragePointer;
        //! Convenience type for list of sections.
        typedef std::vector<Options *> SubSectionList;
        //! Convenience typedef for a map that contains all the options.
        typedef std::map<std::string, AbstractOptionStoragePointer> OptionMap;

        //! Sets the name and title.
        OptionsImpl(const char *name, const char *title);

        /*! \brief
         * Finds a subsection by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found subsection, or NULL if not found.
         *
         * Does not throw.
         */
        Options *findSubSection(const char *name) const;
        /*! \brief
         * Finds an option by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found option, or NULL if not found.
         *
         * Does not throw.
         */
        AbstractOptionStorage *findOption(const char *name) const;

        /*! \brief
         * Calls AbstractOptionStorage::startSource() for all options,
         * including subsections.
         *
         * Does not throw.
         */
        void startSource();

        //! Name for the Options object.
        std::string             name_;
        /*! \brief
         * Option managers set for this collection.
         *
         * This is non-empty only for the top-level Options object.
         */
        OptionManagerContainer  managers_;
        /*! \brief
         * Group that contains all options (and subgroups).
         *
         * This is used to store the insertion order of options.
         */
        Group                   rootGroup_;
        //! Map from option names to options; owns the option storage objects.
        OptionMap               optionMap_;
        /*! \brief
         * List of subsections, in insertion order.
         *
         * This container contains only references to external objects; memory
         * management is performed elsewhere.
         */
        SubSectionList          subSections_;
        //! Options object that contains this object as a subsection, or NULL.
        Options                *parent_;
};

} // namespace internal

} // namespace gmx

#endif
