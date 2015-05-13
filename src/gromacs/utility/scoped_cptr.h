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
 * Declares gmx::scoped_cptr and gmx::scoped_guard_sfree.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_SCOPED_PTR_SFREE_H
#define GMX_UTILITY_SCOPED_PTR_SFREE_H

#include "config.h"

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

//! sfree wrapper to be used as scoped_cptr deleter
template <class T>
inline void sfree_wrapper(T *p)
{
    sfree(p);
}

/*! \libinternal \brief
 * Stripped-down version of scoped_ptr that uses sfree() or custom deleter.
 *
 * Currently only implements some operations; other operations can be added
 * if they become necessary.
 * The presence of a release() method is not strictly according to `scoped_ptr`
 * design, but makes it easier to make existing C code exception-safe, and does
 * not really warrant a separate class for such a purpose.
 *
 * This class provides a basic guard/smart pointer for C pointers.
 *
 * Methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <class T, void D(T *) = sfree_wrapper>
class scoped_cptr
{
    public:
        /*! \brief
         * Initializes a scoped_ptr that frees \p ptr on scope exit.
         *
         * \param[in] ptr  Pointer to use for initialization.
         */
        explicit scoped_cptr(T *ptr = NULL) : ptr_(ptr) {}
        //! Frees the pointer passed to the constructor.
        ~scoped_cptr() { D(ptr_); }
        //! Returns the stored pointer.
        T *get() const { return ptr_; }
        //! Check for non-null pointer in boolean context.
#ifdef GMX_CXX11
        explicit
#endif
        operator bool () const { return ptr_ != 0; }
        //! Sets the pointer and frees previous pointer if necessary.
        void reset(T *ptr) { D(ptr_); ptr_ = ptr; }
        //! Clears the pointer without freeing the memory, and returns the old value.
        T *release() { T *ptr = ptr_; ptr_ = NULL; return ptr; }

    private:
        T                    *ptr_;

        GMX_DISALLOW_COPY_AND_ASSIGN(scoped_cptr);
};

//! Simple guard which calls sfree. See scoped_cptr for details.
typedef scoped_cptr<void> scoped_guard_sfree;

}      // namespace gmx

#endif
