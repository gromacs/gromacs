/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Declares gmx::OptionManagerContainer.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONMANAGERCONTAINER_H
#define GMX_OPTIONS_OPTIONMANAGERCONTAINER_H

#include <vector>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class OptionManagerInterface;

/*! \libinternal
 * \brief
 * Container to keep managers added with Options::addManager() and pass them
 * to options.
 *
 * Consistency of the managers (e.g., that there is at most one manager of a
 * certain type) is only checked when the managers are accessed.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionManagerContainer
{
    public:
        OptionManagerContainer()
        {
        }

        //! Returns `true` if there are no managers.
        bool empty() const { return list_.empty(); }

        //! Adds a manager to the container.
        void add(OptionManagerInterface *manager)
        {
            list_.push_back(manager);
        }
        /*! \brief
         * Retrieves a manager of a certain type.
         *
         * \tparam  ManagerType  Type of manager to retrieve
         *     (should derive from OptionManagerInterface).
         * \returns The manager, or `NULL` if there is none.
         *
         * This method is used in AbstractOption::createStorage() to retrieve
         * a manager of a certain type for options that use a manager.
         *
         * The return value is `NULL` if there is no manager of the given type.
         * The caller needs to handle this (either by asserting, or by handling
         * the manager as optional).
         */
        template <class ManagerType>
        ManagerType *get() const
        {
            ManagerType *result = NULL;
            for (ListType::const_iterator i = list_.begin(); i != list_.end(); ++i)
            {
                ManagerType *curr = dynamic_cast<ManagerType *>(*i);
                if (curr != NULL)
                {
                    GMX_RELEASE_ASSERT(result == NULL,
                                       "More than one applicable option manager is set");
                    result = curr;
                }
            }
            return result;
        }

    private:
        //! Shorthand for the internal container type.
        typedef std::vector<OptionManagerInterface *> ListType;

        ListType  list_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionManagerContainer);
};

} // namespace gmx

#endif
