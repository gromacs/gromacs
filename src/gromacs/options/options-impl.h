/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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

#include <string>
#include <vector>

#include "gromacs/options/abstractoption.h"
#include "gromacs/options/optionmanagercontainer.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/uniqueptr.h"

namespace gmx
{

class AbstractOptionStorage;

/*! \internal
 * \brief
 * Private implementation class for Options.
 *
 * Note that in addition to Options, the OptionsAssigner and OptionsIterator
 * classes also directly access this class.
 *
 * \ingroup module_options
 */
class Options::Impl
{
    public:
        //! Smart pointer for managing an AbstractOptionStorage object.
        typedef gmx_unique_ptr<AbstractOptionStorage>::type
            AbstractOptionStoragePointer;
        //! Convenience type for list of sections.
        typedef std::vector<Options *> SubSectionList;
        //! Convenience type for list of options.
        typedef std::vector<AbstractOptionStoragePointer> OptionList;

        //! Sets the name and title.
        Impl(const char *name, const char *title);
        ~Impl();

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
        //! Description title for the Options object.
        std::string             title_;
        //! Full description for the Options object.
        std::string             description_;
        /*! \brief
         * Option managers set for this collection.
         *
         * This is non-empty only for the top-level Options object.
         */
        OptionManagerContainer  managers_;
        /*! \brief
         * List of subsections, in insertion order.
         *
         * This container contains only references to external objects; memory
         * management is performed elsewhere.
         */
        SubSectionList          subSections_;
        /*! \brief
         * List of options, in insertion order.
         *
         * All objects in this container are owned by this object, and are
         * freed in the destructor.
         */
        OptionList              options_;
        //! Options object that contains this object as a subsection, or NULL.
        Options                *parent_;
};

} // namespace gmx

#endif
