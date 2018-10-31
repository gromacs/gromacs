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
#ifndef MDSPAN_LAYOUTS_H
#define MDSPAN_LAYOUTS_H

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{

// [mdspan.layout]
class layout_right;
class layout_left;
class layout_stride;

}
}
}

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{

class layout_right
{

    public:
        template<class Extents>
        class mapping
        {
            private:

                Extents m_extents;

            public:

                using index_type   = ptrdiff_t;
                using extents_type = Extents;

                constexpr mapping() noexcept = default;

                constexpr mapping( mapping && ) noexcept = default;

                constexpr mapping( const mapping & ) noexcept = default;

                mapping &operator= ( mapping && ) noexcept = default;

                mapping &operator= ( const mapping & ) noexcept = default;

                constexpr mapping( const Extents &ext ) noexcept
                    : m_extents( ext ) {}

                constexpr const Extents &extents() const noexcept { return m_extents; }

            private:

                // ( ( ( ( i0 ) * N1 + i1 ) * N2 + i2 ) * N3 + i3 ) ...

                static constexpr index_type
                offset( const size_t, const ptrdiff_t sum)
                { return sum; }

                template<class ... Indices >
                inline constexpr index_type
                offset( const size_t r, ptrdiff_t sum, const index_type i, Indices... indices) const noexcept
                {
                    return offset( r+1, sum * m_extents.extent(r) + i, indices ...);
                }

            public:

                constexpr index_type required_span_size() const noexcept
                {
                    index_type size = 1;
                    for (size_t r = 0; r < m_extents.rank(); r++)
                    {
                        size *= m_extents.extent(r);
                    }
                    return size;
                }

                template<class ... Indices >
                constexpr
                typename enable_if<sizeof ... (Indices) == Extents::rank(), index_type>::type
                operator()( Indices ... indices ) const noexcept
                { return offset( 0, 0, indices ... ); }

                static constexpr bool is_always_unique()     noexcept { return true; }
                static constexpr bool is_always_contiguous() noexcept { return true; }
                static constexpr bool is_always_strided()    noexcept { return true; }

                constexpr bool is_unique()     const noexcept { return true; }
                constexpr bool is_contiguous() const noexcept { return true; }
                constexpr bool is_strided()    const noexcept { return true; }

                constexpr index_type stride(const size_t R) const noexcept
                {
                    ptrdiff_t stride_ = 1;
                    for (size_t r = m_extents.rank()-1; r > R; r--)
                    {
                        stride_ *= m_extents.extent(r);
                    }
                    return stride_;
                }

        }; // class mapping

};         // class layout_right

}
}
}   // experimental::fundamentals_v3

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{

class layout_left
{
    public:
        template<class Extents>
        class mapping
        {
            private:

                Extents m_extents;

            public:

                using index_type   = ptrdiff_t;
                using extents_type = Extents;

                constexpr mapping() noexcept = default;

                constexpr mapping( mapping && ) noexcept = default;

                constexpr mapping( const mapping & ) noexcept = default;

                mapping &operator= ( mapping && ) noexcept = default;

                mapping &operator= ( const mapping & ) noexcept = default;

                constexpr mapping( const Extents &ext ) noexcept
                    : m_extents( ext ) {}

                constexpr const Extents &extents() const noexcept { return m_extents; }

            private:

                // ( i0 + N0 * ( i1 + N1 * ( i2 + N2 * ( ... ) ) ) )

                static constexpr index_type
                offset( size_t ) noexcept
                { return 0; }

                template<class ... IndexType >
                constexpr index_type
                offset( const size_t r, index_type i, IndexType... indices ) const noexcept
                { return i + m_extents.extent(r) * offset( r+1, indices ... ); }

            public:

                constexpr index_type required_span_size() const noexcept
                {
                    ptrdiff_t size = 1;
                    for (size_t r = 0; r < m_extents.rank(); r++)
                    {
                        size *= m_extents.extent(r);
                    }
                    return size;
                }

                template<class ... Indices >
                constexpr
                typename enable_if<sizeof ... (Indices) == Extents::rank(), index_type>::type
                operator()( Indices ... indices ) const noexcept
                { return offset( 0, indices ... ); }

                static constexpr bool is_always_unique()     noexcept { return true; }
                static constexpr bool is_always_contiguous() noexcept { return true; }
                static constexpr bool is_always_strided()    noexcept { return true; }

                constexpr bool is_unique()     const noexcept { return true; }
                constexpr bool is_contiguous() const noexcept { return true; }
                constexpr bool is_strided()    const noexcept { return true; }

                constexpr index_type stride(const size_t R) const noexcept
                {
                    ptrdiff_t stride_ = 1;
                    for (size_t r = 0; r < R; r++)
                    {
                        stride_ *= m_extents.extent(r);
                    }
                    return stride_;
                }

        }; // class mapping

};         // class layout_left

}
}
}   // experimental::fundamentals_v3

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{

class layout_stride
{
    public:

        template<class Extents>
        class mapping
        {
            private:

                using stride_t = array<ptrdiff_t, Extents::rank()>;

                Extents   m_extents;
                stride_t  m_stride;
                int       m_contig;

            public:

                using index_type   = ptrdiff_t;
                using extents_type = Extents;

                constexpr mapping() noexcept = default;

                constexpr mapping( mapping && ) noexcept = default;

                constexpr mapping( const mapping & ) noexcept = default;

                mapping &operator= ( mapping && ) noexcept = default;

                mapping &operator= ( const mapping & ) noexcept = default;

                mapping( const Extents &ext, const stride_t &str ) noexcept
                    : m_extents(ext), m_stride(str), m_contig(1)
                {
                    int p[ Extents::rank() ? Extents::rank() : 1 ];

                    // Fill permutation such that
                    //   m_stride[ p[i] ] <= m_stride[ p[i+1] ]
                    //
                    for (size_t i = 0; i < Extents::rank(); ++i)
                    {

                        int j = i;

                        while (j && m_stride[i] < m_stride[ p[j-1] ])
                        { p[j] = p[j-1]; --j; }

                        p[j] = i;
                    }

                    for (size_t i = 1; i < Extents::rank(); ++i)
                    {
                        const int        j    = p[i-1];
                        const int        k    = p[i];
                        const index_type prev = m_stride[j] * m_extents.extent(j);
                        if (m_stride[k] != prev) { m_contig = 0; }
                    }
                }

                constexpr const Extents &extents() const noexcept { return m_extents; }

            private:

                // i0 * N0 + i1 * N1 + i2 * N2 + ...

                constexpr index_type
                offset(size_t) const noexcept
                { return 0; }

                template<class ... IndexType >
                constexpr index_type
                offset( const size_t K, const index_type i, IndexType... indices ) const noexcept
                { return i * m_stride[K] + offset(K+1, indices ...); }

            public:

                index_type required_span_size() const noexcept
                {
                    index_type max = 0;
                    for (size_t i = 0; i < Extents::rank(); ++i)
                    {
                        max += m_stride[i] * ( m_extents.extent(i) - 1 );
                    }
                    return max;
                }

                template<class ... Indices >
                constexpr
                typename enable_if<sizeof ... (Indices) == Extents::rank(), index_type>::type
                operator()( Indices ... indices ) const noexcept
                { return offset(0, indices ... ); }


                static constexpr bool is_always_unique()     noexcept { return true; }
                static constexpr bool is_always_contiguous() noexcept { return false; }
                static constexpr bool is_always_strided()    noexcept { return true; }

                constexpr bool is_unique()     const noexcept { return true; }
                constexpr bool is_contiguous() const noexcept { return m_contig; }
                constexpr bool is_strided()    const noexcept { return true; }

                constexpr index_type stride(size_t r) const noexcept
                { return m_stride[r]; }

        }; // class mapping

};         // class layout_stride

}
}
}      // experimental::fundamentals_v3

#endif /* end of include guard: MDSPAN_LAYOUTS_H */
