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
 * \brief Declares gmx::extents for mdspan.
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author Ronan Keryell <ronan.keryell@xilinx.com>
 * \author Carter Edwards <hedwards@nvidia.com>
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup mdspan
 */
#ifndef MDSPAN_EXTENTS_H
#define MDSPAN_EXTENTS_H

#include <cstddef>

#include <array>

namespace gmx
{

/*! \brief Define constant that signals dynamic extent.
 */
enum : std::ptrdiff_t
{
    dynamic_extent = -1
};

template<std::ptrdiff_t... StaticExtents>
class extents;

template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator==(const extents<LHS...>& lhs, const extents<RHS...>& rhs) noexcept;

template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator!=(const extents<LHS...>& lhs, const extents<RHS...>& rhs) noexcept;

namespace detail
{

template<int R, std::ptrdiff_t... StaticExtents>
struct extents_analyse;

/*! \libinternal \brief Enable querying extent of specific rank by splitting
 * a static extents off the variadic template arguments.
 *
 */
template<int R, std::ptrdiff_t E0, std::ptrdiff_t... StaticExtents>
struct extents_analyse<R, E0, StaticExtents...>
{

    //! The extent analysis of the next lower rank.
    using next_extents_analyse = extents_analyse<R - 1, StaticExtents...>;

    /*! \brief Accumulate the total rank from all extents.
     * \returns incremented rank of the next extent
     */
    static constexpr std::size_t rank() noexcept { return next_extents_analyse::rank() + 1; }
    /*! \brief Accumulate the dynamic rank from all extents.
     * This extent is static, so hand down query to the next extent analysis.
     * \returns the dynamic rank of the next extent analysis.
     */
    static constexpr std::size_t rank_dynamic() noexcept
    {
        return next_extents_analyse::rank_dynamic();
    }

    //! Store analysis of the next extent of next lower rank.
    next_extents_analyse next;

    //! Trivial constructor.
    constexpr extents_analyse() : next() {}

    /*! \brief Construction from dynamic extents hands the extents down
     * to the next extents analysis of lower rank.
     * \param[in] de dynamic extents
     */
    template<class... DynamicExtents>
    constexpr extents_analyse(DynamicExtents... de) : next(de...)
    {
    }

    /*! \brief Construct from an array of dynamic extentes and rank.
     * Hand down the dynamic rank parameters to the next extents analysis rank
     * \param[in] de dynamic extents
     * \param[in] r rank to read from the dynamic extent
     */
    template<std::size_t Rank>
    constexpr extents_analyse(const std::array<std::ptrdiff_t, Rank>& de, const std::size_t r) :
        next(de, r)
    {
    }

    //! Copy constructor.
    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse(extents_analyse<R, OtherStaticExtents...> rhs) : next(rhs.next)
    {
    }

    //! Assignment operator.
    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse& operator=(extents_analyse<R, OtherStaticExtents...> rhs)
    {
        next = rhs.next;
        return *this;
    }

    /*! \brief Report extent of dimension r.
     * \param[in] r the dimension to query
     * \returns the extent in dimension r.
     */
    constexpr std::ptrdiff_t extent(const std::size_t r) const noexcept
    {
        return (r == R) ? E0 : next.extent(r);
    }
    /*! \brief Report the static extent of dimension r.
     * \param[in] r the dimension to query
     * \returns the static extent in dimension r.
     */
    static constexpr std::ptrdiff_t static_extent(const std::size_t r) noexcept
    {
        return (r == R) ? E0 : next_extents_analyse::static_extent(r);
    }

    //! Returns the extent with the first dimension sliced off
    constexpr auto sliced_extents() const noexcept { return next; }
};

/*! \libinternal \brief Enable querying extent of specific rank by splitting
 * a dynamic extent off the variadic template arguments.
 */
template<int R, std::ptrdiff_t... StaticExtents>
struct extents_analyse<R, dynamic_extent, StaticExtents...>
{
    //! The extent analysis of the next lower rank.
    using next_extents_analyse = extents_analyse<R - 1, StaticExtents...>;
    /*! \brief Accumulate the total rank from all extents.
     * \returns incremented rank of the next extent
     */
    static constexpr std::size_t rank() noexcept { return next_extents_analyse::rank() + 1; }
    /*! \brief Accumulate the dynamic rank from all extents.
     * \returns the dynamic rank of the next extent analysis.
     */
    static constexpr std::size_t rank_dynamic() noexcept
    {
        return next_extents_analyse::rank_dynamic() + 1;
    }

    //! Store analysis of the next extent of next lower rank.
    next_extents_analyse next;
    //! The dynamic extent of this rank
    std::ptrdiff_t this_extent;

    //! Trivial constructor.
    extents_analyse() : next(), this_extent(0) {}

    /*! \brief Construction from dynamic extents hands the extents down
     * to the next extents analysis of lower rank.
     * \param[in] E the dynamic extent of this rank.
     * \param[in] de dynamic extents
     */
    template<class... DynamicExtents>
    extents_analyse(std::ptrdiff_t E, DynamicExtents... de) : next(de...), this_extent(E)
    {
    }

    /*! \brief Construct from an array of dynamic extentes and rank.
     * Hand down the dynamic rank parameters to the next extents analysis rank
     * \param[in] de dynamic extents
     * \param[in] r rank to read from the dynamic extent
     */
    template<std::size_t Rank>
    extents_analyse(const std::array<std::ptrdiff_t, Rank>& de, const std::size_t r) :
        next(de, r + 1), this_extent(de[r])
    {
    }

    //! Copy constructor.
    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse(extents_analyse<R, OtherStaticExtents...> rhs) :
        next(rhs.next), this_extent(rhs.extent(R))
    {
    }

    //! Assignment operator.
    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse& operator=(extents_analyse<R, OtherStaticExtents...> rhs)
    {
        next        = rhs.next;
        this_extent = rhs.extent(R);
        return *this;
    }

    /*! \brief Report extent of dimension r.
     * \param[in] r the dimension to query
     * \returns the extent in dimension r.
     */
    constexpr std::ptrdiff_t extent(const std::size_t r) const noexcept
    {
        return (r == R) ? this_extent : next.extent(r);
    }
    /*! \brief Report the static extent of dimension r.
     * \param[in] r the dimension to query
     * \returns the static extent in dimension r.
     */
    static constexpr std::ptrdiff_t static_extent(const std::size_t r) noexcept
    {
        return (r == R) ? dynamic_extent : next_extents_analyse::static_extent(r);
    }

    //! Returns the extent with the first dimension sliced off
    constexpr auto sliced_extents() const noexcept { return next; }
};

/*! \libinternal \brief Specialisation for rank 0 extents analysis.
 * Ends recursive rank analysis.
 */
template<>
struct extents_analyse<0>
{
    /*! \brief Rank of extent of rank 0.
     * \returns 0
     */
    static constexpr std::size_t rank() noexcept { return 0; }
    /*! \brief Dynamic rank of extent of rank 0.
     * \returns 0
     */
    static constexpr std::size_t rank_dynamic() noexcept { return 0; }

    //! Trivial constructor.
    constexpr extents_analyse() {}

    //! Construct from array and rank, doing nothing.
    template<std::size_t Rank>
    extents_analyse(const std::array<std::ptrdiff_t, Rank>& /*de*/, const std::size_t /*r*/)
    {
    }

    // extents_analyse & operator=(extents_analyse) = default;

    /*! \brief Extent of rank 0 is 1, ensuring that product of extents yields required size and not zero.
     * NOTE changed from ORNL reference implementation in making this static constexpr instead of constexpr .. const
     */
    static constexpr std::ptrdiff_t extent(const std::size_t /*r*/) noexcept { return 1; }

    //! Static extent of rank 0 is 1, ensuring that product of extents yields required size and not zero.
    static constexpr std::ptrdiff_t static_extent(const std::size_t /*r*/) noexcept { return 1; }
};

template<std::ptrdiff_t E0, std::ptrdiff_t... StaticExtents>
struct sliced_extents
{
    using type = extents<StaticExtents...>;
};
} // namespace detail

/*! \libinternal \brief Multidimensional extents with static and dynamic dimensions.
 *
 * Describes a multidimensional index space of rank R.
 * This is equivalent to the Cartesian product space of integer intervals
 * [0, N_0) x [0, N_1) x ... x [0,N_{R-1} )
 *
 * Confer to P0009r8 of the Library Evolution Working Group and mdspan.extents
 *
 * \tparam StaticExtents rank number of extents, where the dynamic_extent
 * constant for static extent is used to signal a dynamic extent.
 */
template<std::ptrdiff_t... StaticExtents>
class extents
{
private:
    using extents_analyse_t = detail::extents_analyse<sizeof...(StaticExtents), StaticExtents...>;

public:
    //! Type used to index elements.
    using index_type = std::ptrdiff_t;
    //! Trivial constructor
    constexpr extents() noexcept {}
    //! Move constructor
    constexpr extents(extents&&) noexcept = default;
    //! Copy constructor.
    constexpr extents(const extents&) noexcept = default;
    /*! \brief Construct with dynamic extents.
     *
     * Allows for extents(u,v,w..) syntax when setting dynamic extents
     *
     * \tparam IndexType type of index
     * \param[in] dn first dynamic index
     * \param[in] DynamicExtents parameter pack
     */
    template<class... IndexType>
    constexpr extents(std::ptrdiff_t dn, IndexType... DynamicExtents) noexcept :
        impl(dn, DynamicExtents...)
    {
        static_assert(1 + sizeof...(DynamicExtents) == rank_dynamic(), "");
    }

    /*! \brief Construct from array of dynamic extents.
     *
     * Allows for extents({{u,v,w..}}) syntax when setting dynamic extents
     *
     * \param[in] dynamic_extents array of dynamic rank size containing extents
     */
    constexpr extents(const std::array<std::ptrdiff_t, extents_analyse_t::rank_dynamic()> dynamic_extents) noexcept
        :
        impl(dynamic_extents, 0)
    {
    }

    //! Copy constructor
    template<std::ptrdiff_t... OtherStaticExtents>
    extents(const extents<OtherStaticExtents...>& other) : impl(other.impl)
    {
    }

    //! Default move assignment
    extents& operator=(extents&&) noexcept = default;
    //! Default copy assignment
    extents& operator=(const extents&) noexcept = default;
    //! Copy assignment
    template<std::ptrdiff_t... OtherStaticExtents>
    extents& operator=(const extents<OtherStaticExtents...>& other)
    {
        impl = other.impl;
        return *this;
    }
    //! Default destructor
    ~extents() = default;

    // [mdspan.extents.obs]
    /*! \brief The rank of the extent.
     * \returns the rank all extents together
     */
    static constexpr std::size_t rank() noexcept { return sizeof...(StaticExtents); }
    /*! \brief The rank of the dynamic extents.
     * \returns Only the dynamic extents.
     */
    static constexpr std::size_t rank_dynamic() noexcept
    {
        return extents_analyse_t::rank_dynamic();
    }
    /*! \brief The rank of the static extents.
     * \returns Only the static extents.
     */
    static constexpr index_type static_extent(std::size_t k) noexcept
    {
        return extents_analyse_t::static_extent(rank() - k);
    }
    /*! \brief The extent along a specific dimension.
     * \param[in] k the dimension
     * \returns the extent along that dimension
     */
    constexpr index_type extent(std::size_t k) const noexcept { return impl.extent(rank() - k); }
    //! Returns the extent with the first dimension sliced off
    constexpr auto sliced_extents() const noexcept
    {
        return typename detail::sliced_extents<StaticExtents...>::type(impl.sliced_extents());
    }

private:
    extents(extents_analyse_t o) : impl(o) {}
    //! For copy assignment, extents are friends of extents.
    template<std::ptrdiff_t...>
    friend class extents;
    //! The implementation class.
    extents_analyse_t impl;
};


/*! \brief Comparison operator.
 * \returns true if extents are equal
 */
template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator==(const extents<LHS...>& lhs, const extents<RHS...>& rhs) noexcept
{
    bool equal = lhs.rank() == rhs.rank();
    for (std::size_t r = 0; r < lhs.rank(); r++)
    {
        equal = equal && (lhs.extent(r) == rhs.extent(r));
    }
    return equal;
}

/*! \brief Check for non-equality.
 * \returns true if extents are unequal
 */
template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator!=(const extents<LHS...>& lhs, const extents<RHS...>& rhs) noexcept
{
    return !(lhs == rhs);
}

} // namespace gmx
#endif /* end of include guard: MDSPAN_EXTENTS_H */
