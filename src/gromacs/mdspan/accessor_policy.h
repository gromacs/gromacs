/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
/*
 * This file is a modified version of original work of Sandia Corporation.
 * In the spirit of the original code, this particular file can be distributed
 * on the terms of Sandia Corporation.
 */
/*
 *                         Kokkos v. 2.0
 *               Copyright (2014) Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Kokkos is licensed under 3-clause BSD terms of use:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Christian R. Trott (crtrott@sandia.gov)
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
 * \tparam ElementType the type held in memory to be accessed
 */
template<class ElementType>
class accessor_basic
{
public:
    //! Type of element to be accessed.
    using element_type = ElementType;
    //! Pointer to element to be accessed.
    using pointer = ElementType*;
    //! How to determine a memory offset, provided by self accessor.
    using offset_policy = accessor_basic;
    //! Type of references.
    using reference = ElementType&;

    /*! \brief Shift a pointer by an offset.
     * \param[in] p Pointer to reference memory location.
     * \param[in] i offset from memory location.
     * \returns pointer to offset memory location.
     */
    constexpr typename offset_policy::pointer offset(pointer p, ptrdiff_t i) const noexcept
    {
        return typename offset_policy::pointer(p + i);
    }

    /*! \brief Access element from an offset to given pointer.
     * \param[in] p Pointer to reference memory location.
     * \param[in] i offset from memory location.
     * \returns reference to element stored at offset from memory location.
     */
    constexpr reference access(pointer p, ptrdiff_t i) const noexcept { return p[i]; }

    /*! \brief Decay pointer to pointer to ElementType.
     * NOTE This function does nothing, because it is the trivial implementation of an accessor.
     * \returns input pointer as pointer to ElementType
     */
    constexpr ElementType* decay(pointer p) const noexcept { return p; }
};

} // namespace gmx
#endif /* end of include guard: MDSPAN_ACCESSOR_POLICY_H */
