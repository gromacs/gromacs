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
#ifndef MDSPAN_SUBSPAN_H
#define MDSPAN_SUBSPAN_H

namespace std
{
namespace experimental
{
inline namespace fundamentals_v3
{
namespace detail
{

template<class ExtentsNew, class ExtentsOld, class ... SliceSpecifiers>
struct compose_new_extents;

template<ptrdiff_t ... ExtentsNew, ptrdiff_t E0, ptrdiff_t ... ExtentsOld, class ... SliceSpecifiers>
struct compose_new_extents<extents<ExtentsNew...>, extents<E0, ExtentsOld...>, all_type, SliceSpecifiers...> {
    typedef compose_new_extents<extents<ExtentsNew..., E0>, extents<ExtentsOld...>, SliceSpecifiers...> next_compose_new_extents;
    typedef typename next_compose_new_extents::extents_type extents_type;

    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents e,
                                                     array<ptrdiff_t, OrgExtents::rank()> &strides, ptrdiff_t &offset,
                                                     all_type, SliceSpecifiers... s, DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...);
    }
};
template<ptrdiff_t ... ExtentsNew, ptrdiff_t ... ExtentsOld, class ... SliceSpecifiers>
struct compose_new_extents<extents<ExtentsNew...>, extents<dynamic_extent, ExtentsOld...>, all_type, SliceSpecifiers...> {
    typedef compose_new_extents<extents<ExtentsNew..., dynamic_extent>, extents<ExtentsOld...>, SliceSpecifiers...> next_compose_new_extents;
    typedef typename next_compose_new_extents::extents_type extents_type;

    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents e,
                                                     array<ptrdiff_t, OrgExtents::rank()> &strides, ptrdiff_t &offset,
                                                     all_type, SliceSpecifiers... s, DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ..., e.extent(OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1));
    }
};

template<ptrdiff_t ... ExtentsNew, ptrdiff_t E0, ptrdiff_t ... ExtentsOld, class IT, class ... SliceSpecifiers>
struct compose_new_extents<extents<ExtentsNew...>, extents<E0, ExtentsOld...>, pair<IT, IT>, SliceSpecifiers...> {
    typedef compose_new_extents<extents<ExtentsNew..., dynamic_extent>, extents<ExtentsOld...>, SliceSpecifiers...> next_compose_new_extents;
    typedef typename next_compose_new_extents::extents_type extents_type;

    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents e, array<ptrdiff_t, OrgExtents::rank()> &strides, ptrdiff_t &offset,
                                                     pair<IT, IT> p, SliceSpecifiers ... s, DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        offset += p.first*strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ..., ptrdiff_t(p.second-p.first));
    }
};

template<ptrdiff_t ... ExtentsNew, ptrdiff_t E0, ptrdiff_t ... ExtentsOld, class IT, class ... SliceSpecifiers>
struct compose_new_extents<extents<ExtentsNew...>, extents<E0, ExtentsOld...>, IT, SliceSpecifiers...> {
    typedef compose_new_extents<extents<ExtentsNew...>, extents<ExtentsOld...>, SliceSpecifiers...> next_compose_new_extents;
    typedef typename next_compose_new_extents::extents_type extents_type;

    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents e, array<ptrdiff_t, OrgExtents::rank()> &strides, ptrdiff_t &offset,
                                                     const ptrdiff_t v, SliceSpecifiers... s, DynamicExtents... de)
    {
        offset += v*strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...);
    }
};

template<ptrdiff_t ... ExtentsNew>
struct compose_new_extents < extents<ExtentsNew...>, extents <>> {
    typedef extents<ExtentsNew...> extents_type;

    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents, array<ptrdiff_t, OrgExtents::rank()>, ptrdiff_t, DynamicExtents... de)
    {
        return extents_type(de ...);
    }
    template<class OrgExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents, array<ptrdiff_t, OrgExtents::rank()>, ptrdiff_t)
    {
        return extents_type();
    }
};


template<class Extents, class ... SliceSpecifiers>
struct subspan_deduce_extents {
    typedef typename compose_new_extents<extents<>, Extents, SliceSpecifiers...>::extents_type extents_type;
    typedef array<ptrdiff_t, Extents::rank()> stride_type;
    static constexpr extents_type create_sub_extents(const Extents e, stride_type &strides, ptrdiff_t &offset, SliceSpecifiers... s)
    {
        return compose_new_extents<extents<>, Extents, SliceSpecifiers...>::create_sub_extents(e, strides, offset, s ...);
    }
};

}

template<class ElementType, class Extents, class LayoutPolicy,
         class AccessorPolicy, class ... SliceSpecifiers>
basic_mdspan<ElementType, typename detail::subspan_deduce_extents<Extents, SliceSpecifiers...>::extents_type,
             layout_stride, typename AccessorPolicy::offset_policy >
subspan(const basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy> &src, SliceSpecifiers ... slices) noexcept
{
    typedef typename detail::subspan_deduce_extents<Extents, SliceSpecifiers...>::extents_type sub_extents_type;
    typedef typename AccessorPolicy::offset_policy sub_accessor_policy;
    typedef basic_mdspan<ElementType, sub_extents_type, layout_stride, sub_accessor_policy> sub_mdspan_type;

    array<ptrdiff_t, Extents::rank()> strides;
    for (size_t r = 0; r < Extents::rank(); r++)
    {
        strides[r] = src.stride(r);
    }

    ptrdiff_t        offset      = 0;
    sub_extents_type sub_extents = detail::subspan_deduce_extents<Extents, SliceSpecifiers...>::create_sub_extents(src.extents(), strides, offset, slices ...);

    array<ptrdiff_t, sub_extents_type::rank()> sub_strides;
    for (size_t r = 0; r < sub_extents_type::rank(); r++)
    {
        sub_strides[r] = strides[r];
    }

    typename AccessorPolicy::offset_policy::pointer ptr = src.accessor().offset(src.data(), offset);
    return sub_mdspan_type(ptr, typename sub_mdspan_type::mapping_type(sub_extents, sub_strides));
}

}
}
}

#endif /* end of include guard: MDSPAN_SUBSPAN_H */
