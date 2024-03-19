/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
/*
 * This file is a modified version of original work of Sandia Corporation.
 * In the spirit of the original code, this particular file can be distributed
 * on the terms of Sandia Corporation.
 */
/*
 *                          Kokkos v. 2.0
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
 * \brief Declares gmx::layout_right for mdspan.
 *
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup mdspan
 */
#ifndef MDSPAN_LAYOUTS_H
#define MDSPAN_LAYOUTS_H
#include <cstddef>

#include <type_traits>

namespace gmx
{

/*! \libinternal \brief Right-aligned array layout indexer.
 * Carries the mapping class performing the translation from multidimensional
 * index to one-dimensional number.
 */
class layout_right
{

public:
    /*! \libinternal \brief Mapping from multidimensional indices within extents to 1D index.
     * \tparam Extents the extents of the multidimensional integers for the mapping.
     */
    template<class Extents>
    class mapping
    {
    private:
        //! The extents.
        Extents m_extents;

    public:
        //! exposing the type of indices
        using index_type = ptrdiff_t;
        //! exposing the type of the extents
        using extents_type = Extents;
        //! Default constructor.
        constexpr mapping() noexcept = default;
        //! Default move constructor.
        constexpr mapping(mapping&&) noexcept = default;
        //! Default copy constructor.
        constexpr mapping(const mapping&) noexcept = default;
        //! Default move assignment
        mapping& operator=(mapping&&) noexcept = default;
        //! Default copy assignment
        mapping& operator=(const mapping&) noexcept = default;
        /*! \brief Construct mapping, setting extents
         * \param[in] ext the extents
         */
        constexpr mapping(const Extents& ext) noexcept : m_extents(ext) {}
        /*! \brief Return the extents.
         * \returns extents
         */
        constexpr const Extents& extents() const noexcept { return m_extents; }

    private:
        /* \brief End recursion helper function for static offset calculation.
         * \param[in] sum The accumulated offset over all dimensions
         * \returns The offset.
         */
        static constexpr index_type offset(const size_t /*r*/, const ptrdiff_t sum) { return sum; }

        /* \brief Statically calculate offset from index and extent.
         * For a multidimensional index (i0,i1,..,in), in a right memory
         * layout, 'i0' denotes the slowest moving dimension and
         * 'in' the fastest moving dimension.
         * The overall offset within extents N = (N0,..,Nn) is then
         * offest = i0 * N1 * .. * Nn + i1 * N2 * .. * Nn + in-1 * Nn + in
         *        = (((i0*N1+i1)*N2+i2)*N3+i3) ...
         * \param[in] r current rank
         * \param[in] sum current sum up to this rank
         * \param[in] i index
         * \oaram[in] indices The rest of the parameter pack.
         * \returns The offset.
         */
        template<class... Indices>
        constexpr index_type offset(const size_t r, ptrdiff_t sum, const index_type i, Indices... indices) const noexcept
        {
            return offset(r + 1, sum * m_extents.extent(r) + i, indices...);
        }

    public:
        /*! \brief Return the size of the underlying one-dimensional
         * data structure, so that the mapping is always valid.
         *
         * \returns number of span elements
         */
        constexpr index_type required_span_size() const noexcept
        {
            index_type size = 1;
            for (size_t r = 0; r < m_extents.rank(); r++)
            {
                size *= m_extents.extent(r);
            }
            return size;
        }

        /*! \brief Map the multidimensional indices to 1D.
         * Requires number of indicies have the same dimensionality as the mapping.
         * \tparam Indices type of the indices to be mapped
         * \param[in] indices the indices to be mapped
         * \returns One-dimensional integer index.
         */
        template<class... Indices>
        std::enable_if_t<sizeof...(Indices) == Extents::rank(), index_type> constexpr
        operator()(Indices... indices) const noexcept
        {
            return offset(0, 0, indices...);
        }

        //! Report that this mapping is always unique.
        static constexpr bool is_always_unique() noexcept { return true; }
        //! Report that this mapping is always contiguous.
        static constexpr bool is_always_contiguous() noexcept { return true; }
        //! Report that this mapping is always strided.
        static constexpr bool is_always_strided() noexcept { return true; }

        //! Report that this mapping is unique.
        constexpr bool is_unique() const noexcept { return true; }
        //! Report that this mapping is contiguous.
        constexpr bool is_contiguous() const noexcept { return true; }
        //! Report that this mapping is strided.
        constexpr bool is_strided() const noexcept { return true; }
        /*!\brief Return the stride of dimension r.
         * \param[in] R rank of the stride to be queried.
         * \returns the stride along dimension r.
         */
        constexpr index_type stride(const size_t R) const noexcept
        {
            ptrdiff_t stride = 1;
            for (size_t r = m_extents.rank() - 1; r > R; r--)
            {
                stride *= m_extents.extent(r);
            }
            return stride;
        }

    }; // class mapping

}; // class layout_right

} // namespace gmx
#endif /* end of include guard: MDSPAN_LAYOUTS_H */
