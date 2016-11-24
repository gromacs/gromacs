/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014,2015,2016, by the GROMACS development team, led by
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
 * Declares gmx::unique_cptr and gmx::unique_guard_sfree.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_UNIQUE_PTR_SFREE_H
#define GMX_UTILITY_UNIQUE_PTR_SFREE_H

#include <memory>

#include "gromacs/utility/smalloc.h"

namespace gmx
{

//! sfree wrapper to be used as unique_cptr deleter
template <class T>
inline void sfree_wrapper(T *p)
{
    sfree(p);
}

/*! \libinternal \brief
 * unique-ptr which takes function pointer as template argument
 *
 * By default calls sfree as deleter.
 *
 * Methods in this class do not throw.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
template <class T, void D(T *) = sfree_wrapper>
class unique_cptr : public std::unique_ptr<T, void(*)(T*)>
{
    public:
        //! Initializes a unique_cptr to nulltpr.
        constexpr unique_cptr() : std::unique_ptr<T, void(*) (T*)>(nullptr, D) {};
        /*! \brief
         * Initializes a unique_cptr that frees \p ptr on destruction.
         *
         * \param[in] ptr  Pointer to use for initialization.
         */
        explicit unique_cptr(T* p) : std::unique_ptr<T, void(*) (T*)>(p, D) {};
};

//! Simple guard which calls sfree. See unique_cptr for details.
typedef unique_cptr<void> unique_guard_sfree;

}      // namespace gmx

#endif
