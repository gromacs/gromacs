/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 * Declares gmx::OptionsBehaviorCollection.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_BEHAVIORCOLLECTION_H
#define GMX_OPTIONS_BEHAVIORCOLLECTION_H

#include <memory>
#include <vector>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class IOptionsBehavior;
class Options;

//! Smart pointer for behaviors stored in OptionsBehaviorCollection.
typedef std::shared_ptr<IOptionsBehavior> OptionsBehaviorPointer;

/*! \libinternal \brief
 * Container for IOptionsBehavior objects.
 *
 * This class provides a container to keep IOptionsBehavior objects, and to
 * call the IOptionsBehavior methods for the contained objects.
 *
 * IOptionsBehavior methods are called for the contained objects in the same
 * order as in which the behaviors were inserted.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsBehaviorCollection
{
    public:
        /*! \brief
         * Constructs a container for storing behaviors associated with given
         * Options.
         *
         * Caller needs to ensure that provided Options remains in existence
         * while the container exists.
         */
        explicit OptionsBehaviorCollection(Options *options);
        ~OptionsBehaviorCollection();

        //! Adds a behavior to the collection.
        void addBehavior(const OptionsBehaviorPointer &behavior);
        //! Calls IOptionsBehavior::optionsFinishing() on all behaviors.
        void optionsFinishing();
        //! Calls IOptionsBehavior::optionsFinished() on all behaviors.
        void optionsFinished();

    private:
        Options                             *options_;
        std::vector<OptionsBehaviorPointer>  behaviors_;

        GMX_DISALLOW_COPY_AND_ASSIGN(OptionsBehaviorCollection);
};

} // namespace

#endif
