/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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
 * Declares gmx::scoped_ptr_sfree.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_SCOPED_PTR_SFREE_H
#define GMX_UTILITY_SCOPED_PTR_SFREE_H

#include "common.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

/*! \libinternal \brief
 * Stripped-down version of scoped_ptr that uses sfree().
 *
 * Currently only implements constructor from a pointer value and destructor;
 * other operations can be added if they become necessary.
 *
 * This class provides a very basic guard/smart pointer for C pointers.
 * Currently, boost::shared_ptr is used in a few locations that require more
 * flexibility, but is not suitable for all cases either.  A scoped_ptr with
 * deleter support would be a general enough implementation for all uses.
 * C++11 unique_ptr has this, but for non-C++11 support we need something else.
 *
 * Methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class scoped_ptr_sfree
{
    public:
        /*! \brief
         * Initializes a scoped_ptr that frees \p ptr on scope exit.
         *
         * \param[in] ptr  Pointer to use for initialization.
         */
        explicit scoped_ptr_sfree(void *ptr) : ptr_(ptr) {}
        //! Frees the pointer passed to the constructor.
        ~scoped_ptr_sfree() { sfree(ptr_); }

    private:
        void                   *ptr_;

        GMX_DISALLOW_COPY_AND_ASSIGN(scoped_ptr_sfree);
};

}      // namespace gmx

#endif
