/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief Declares accessor policies for mdspan.
 *
 * Implement ways how to convert a linear offset from a pointer to memory access.
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 * \libinternal
 * \ingroup mdspan
 */
#ifndef MDSPAN_ACCESSOR_POLICY_H
#define MDSPAN_ACCESSOR_POLICY_H

#include <cstddef>

namespace gmx
{

/*! \libinternal \brief The most basic memory access model for mdspan.
 * \tparam ElemntType the type held in memory to be accessed
 */
template<class ElementType>
class accessor_basic
{
    public:
        //! Type of element to be accessed.
        using element_type  = ElementType;
        //! Pointer to element to be accessed.
        using pointer       = ElementType*;
        //! How to determine a memory offset, provided by self accessor.
        using offset_policy = accessor_basic;
        //! Type of references.
        using reference     = ElementType&;

        /*! \brief Shift a pointer by an offset.
         * \param[in] p Pointer to reference memory location.
         * \param[in] i offset from memory location.
         * \returns pointer to offset memory location.
         */
        constexpr typename offset_policy::pointer
        offset( pointer p, ptrdiff_t i ) const noexcept
        { return typename offset_policy::pointer(p+i); }

        /*! \brief Access element from an offset to given pointer.
         * \param[in] p Pointer to reference memory location.
         * \param[in] i offset from memory location.
         * \returns reference to element stored at offset from memory location.
         */
        constexpr reference access( pointer p, ptrdiff_t i ) const noexcept
        { return p[i]; }

        /*! \brief Decay pointer to pointer to ElementType.
         * NOTE This function does nothing, because it is the trivial implementaion of an accessor.
         * \returns input pointer as poointer to ElementType
         */
        constexpr ElementType* decay( pointer p ) const noexcept
        { return p; }
};

}      // namespace gmx
#endif /* end of include guard: MDSPAN_ACCESSOR_POLICY_H */
