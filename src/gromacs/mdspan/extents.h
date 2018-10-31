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
#ifndef MDSPAN_EXTENTS_H
#define MDSPAN_EXTENTS_H

#include <cstddef> // std::ptrdiff_t

#include <array>   // std::array

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{

enum : std::ptrdiff_t {
    dynamic_extent = -1
};


// [mdspan.extents]
template< std::ptrdiff_t ... StaticExtents >
class extents;

// [mdspan.extents.compare]
template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
bool operator==(const extents<LHS...> &lhs,
                const extents<RHS...> &rhs) noexcept;

template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator!=(const extents<LHS...> &lhs,
                          const extents<RHS...> &rhs) noexcept;


} // namespace fundamentals_v3
} // namespace experimental
} // namespace std

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{
namespace detail
{
template< int R, std::ptrdiff_t ... StaticExtents >
struct extents_analyse;

template< int R, std::ptrdiff_t E0, std::ptrdiff_t ... StaticExtents >
struct extents_analyse<R, E0, StaticExtents...> {

    typedef extents_analyse<R-1, StaticExtents...> next_extents_analyse;

    static constexpr std::size_t rank() noexcept { return next_extents_analyse::rank()+1; }
    static constexpr std::size_t rank_dynamic() noexcept { return next_extents_analyse::rank_dynamic(); }

    next_extents_analyse next;

    extents_analyse() : next() {};

    template<class ... DynamicExtents>
    extents_analyse(DynamicExtents... de) : next(de ...) {}

    template<std::size_t Rank>
    extents_analyse(const array<std::ptrdiff_t, Rank> &de, const std::size_t r) : next(de, r) {}

    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse(extents_analyse<R, OtherStaticExtents...> rhs) : next(rhs.next) {}

    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse operator= (extents_analyse<R, OtherStaticExtents...> rhs)
    {
        next = rhs.next;
        return *this;
    }

    constexpr std::ptrdiff_t extent(const std::size_t r) const noexcept
    {
        return (r == R) ? E0 : next.extent(r); // TODO use early return and if instead of ternary operator after bumping requirements to C++14
    }
    static constexpr std::ptrdiff_t static_extent(const std::size_t r) noexcept
    {
        return (r == R) ? E0 : next_extents_analyse::static_extent(r); // TODO use early return and if instead of ternary operator after bumping requirements to C++14
    }
};

template< int R, std::ptrdiff_t ... StaticExtents >
struct extents_analyse<R, dynamic_extent, StaticExtents...> {
    typedef extents_analyse<R-1, StaticExtents...> next_extents_analyse;

    static constexpr std::size_t rank() noexcept { return next_extents_analyse::rank()+1; }
    static constexpr std::size_t rank_dynamic() noexcept { return next_extents_analyse::rank_dynamic()+1; }

    next_extents_analyse next;
    std::ptrdiff_t       this_extent;

    extents_analyse() : next(), this_extent(0) {}

    template<class ... DynamicExtents>
    extents_analyse(std::ptrdiff_t E, DynamicExtents... de) : next(de ...), this_extent(E) {}

    template<std::size_t Rank>
    extents_analyse(const array<std::ptrdiff_t, Rank> &de, const std::size_t r) : next(de, r+1), this_extent(de[r]) {}

    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse(extents_analyse<R, OtherStaticExtents...> rhs) : next(rhs.next), this_extent(rhs.extent(R)) {}

    template<std::ptrdiff_t... OtherStaticExtents>
    extents_analyse &operator= (extents_analyse<R, OtherStaticExtents...> rhs)
    {
        next        = rhs.next;
        this_extent = rhs.extent(R);
        return *this;
    }

    constexpr std::ptrdiff_t extent(const std::size_t r) const noexcept
    {
        return (r == R) ? this_extent : next.extent(r);
    }
    static constexpr std::ptrdiff_t static_extent(const std::size_t r) noexcept
    {
        return (r == R) ? dynamic_extent : next_extents_analyse::static_extent(r);
    }
};

template<>
struct extents_analyse<0> {
    static constexpr std::size_t rank() noexcept { return 0; }
    static constexpr std::size_t rank_dynamic() noexcept { return 0; }

    extents_analyse() {}

    template<std::size_t Rank>
    extents_analyse(const array<std::ptrdiff_t, Rank> &, const std::size_t) {}

    //extents_analyse & operator=(extents_analyse) = default;

    constexpr std::ptrdiff_t extent(const std::size_t) const noexcept
    {
        return 1;
    }
    static constexpr std::ptrdiff_t static_extent(const std::size_t) noexcept
    {
        return 1;
    }

};
}

template< std::ptrdiff_t ... StaticExtents >
class extents
{
    private:

        template< std::ptrdiff_t... > friend class extents;

        typedef detail::extents_analyse<sizeof ... (StaticExtents), StaticExtents...> extents_analyse_t;
        extents_analyse_t impl;
    public:

        using index_type = std::ptrdiff_t;

        constexpr extents() noexcept {}

        constexpr extents( extents && ) noexcept = default;

        constexpr extents( const extents & ) noexcept = default;

        template< class ... IndexType >
        constexpr extents( std::ptrdiff_t dn,
                           IndexType ...  DynamicExtents ) noexcept
            : impl( dn, DynamicExtents ... )
        { static_assert( 1+sizeof ... (DynamicExtents) == rank_dynamic(), "" ); }

        constexpr extents( const array<std::ptrdiff_t, extents_analyse_t::rank_dynamic()> dynamic_extents) noexcept
            : impl(dynamic_extents, 0) {}

        template<std::ptrdiff_t... OtherStaticExtents>
        extents( const extents<OtherStaticExtents...> &other )
            : impl( other.impl ) {}

        extents &operator= ( extents && ) noexcept = default;

        extents &operator= ( const extents & ) noexcept = default;

        template<std::ptrdiff_t... OtherStaticExtents>
        extents &operator= ( const extents<OtherStaticExtents...> &other )
        { impl = other.impl; return *this; }

        ~extents() = default;

        // [mdspan.extents.obs]

        static constexpr std::size_t rank() noexcept
        { return sizeof ... (StaticExtents); }

        static constexpr std::size_t rank_dynamic() noexcept
        { return extents_analyse_t::rank_dynamic(); }

        static constexpr index_type static_extent(std::size_t k) noexcept
        { return extents_analyse_t::static_extent(rank()-k); }

        constexpr index_type extent(std::size_t k) const noexcept
        { return impl.extent(rank()-k); }

};

// make constexpr again in C++14
template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
bool operator==(const extents<LHS...> &lhs,
                const extents<RHS...> &rhs) noexcept
{
    bool equal = lhs.rank() == rhs.rank();
    for (std::size_t r = 0; r < lhs.rank(); r++)
    {
        equal = equal && ( lhs.extent(r) == rhs.extent(r) );
    }
    return equal;
}

template<std::ptrdiff_t... LHS, std::ptrdiff_t... RHS>
constexpr bool operator!=(const extents<LHS...> &lhs,
                          const extents<RHS...> &rhs) noexcept
{
    return !(lhs == rhs);
}

}
}
}      // std::experimental::fundamentals_v3
#endif /* end of include guard: MDSPAN_EXTENTS_H */
