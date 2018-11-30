/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Declares gmx::mdspan
 *
 * \author Christian Trott <crtrott@sandia.gov>
 * \author Ronan Keryell <ronan.keryell@xilinx.com>
 * \author Carter Edwards <hedwards@nvidia.com>
 * \author David Hollman <dshollm@sandia.gov>
 * \author Christian Blau <cblau@gwdg.de>
 * \inlibraryapi
 * \ingroup module_mdspan
 */

#ifndef MDSPAN_SUBSPAN_H
#define MDSPAN_SUBSPAN_H

#include "gromacs/mdspan/mdspan.h"

namespace gmx
{

/*! \libinternal \brief
 * Tags access to all data along a rank.
 */
class all_type
{
    public:
        constexpr explicit all_type() = default;
};

namespace detail
{

/*! \internal \brief
 * Declare template class to compose new extents.
 */
template<class ExtentsNew, class ExtentsOld, class ... SliceSpecifiers>
struct compose_new_extents;

/*! \internal \brief
 * Slicing off last new extent from extents and updating strides and offset.
 *
 * \tparam ExtentsNew newly constructed extents of the subspan
 * \tparam E0 Currently spliced of extent
 */
template<
    ptrdiff_t ... ExtentsNew,
    ptrdiff_t E0,
    ptrdiff_t ... ExtentsOld,
    class ... SliceSpecifiers>
struct compose_new_extents<
    extents<ExtentsNew...>,
    extents<E0, ExtentsOld...>,
    all_type,
    SliceSpecifiers...>
{
    //! Type of the next new extents with the next new extent pulled out of the parameter pack
    using next_compose_new_extents = compose_new_extents<
                extents<ExtentsNew..., E0>,
                extents<ExtentsOld...>,
                SliceSpecifiers... >;
    //! Type of the new extents
    using extents_type = typename next_compose_new_extents::extents_type;

    /*! \brief Return the subspan extents while setting offset and strides.
     * \todo reenable constexpr with C++14
     */
    template<class OrgExtents, class ... DynamicExtents>
    static extents_type create_sub_extents(
            const OrgExtents e,
            std::array<ptrdiff_t, OrgExtents::rank()> &strides,
            ptrdiff_t *offset,
            all_type /*alltype*/,
            SliceSpecifiers... s,
            DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...);
    }
};

/*! \internal \brief
 * Template specialisation to compose new extents when the last slice in the
 * SliceSpecifier parameter pack is the whole extent.
 */
template<
    ptrdiff_t ... ExtentsNew,
    ptrdiff_t ... ExtentsOld,
    class ... SliceSpecifiers>
struct compose_new_extents<
    extents<ExtentsNew...>,
    extents<dynamic_extent, ExtentsOld...>,
    all_type,
    SliceSpecifiers...>
{
    //! \copydoc compose_new_extents::next_compose_new_extents
    using next_compose_new_extents = compose_new_extents<
                extents<ExtentsNew...,
                        dynamic_extent>,
                extents<ExtentsOld...>,
                SliceSpecifiers...>;

    //! \copydoc compose_new_extents::extents_type
    using extents_type = typename next_compose_new_extents::extents_type;

    /*! \brief
     * Return the subspan extents while setting offset and strides using all_type as last parameter.
     * \todo reenable constexpr with C++14
     */
    template <class OrgExtents, class ... DynamicExtents>
    static extents_type create_sub_extents(const OrgExtents e,
                                           std::array<ptrdiff_t, OrgExtents::rank()> &strides,
                                           ptrdiff_t *offset,
                                           all_type /*alltype*/,
                                           SliceSpecifiers... s,
                                           DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...,
                                                            e.extent(OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1));
    }
};


/*! \internal \brief
 * Template specialisation to compose new extents when the last slice in the
 * parameter pack is a pair of integral numbers.
 */
template<
    ptrdiff_t ... ExtentsNew,
    ptrdiff_t E0,
    ptrdiff_t ... ExtentsOld,
    class IT,
    class ... SliceSpecifiers>
struct compose_new_extents<
    extents<ExtentsNew...>,
    extents<E0, ExtentsOld...>,
    std::pair<IT, IT>,
    SliceSpecifiers...>
{
    //! \copydoc compose_new_extents::next_compose_new_extents
    using next_compose_new_extents = compose_new_extents<
                extents<ExtentsNew..., dynamic_extent>,
                extents<ExtentsOld...>,
                SliceSpecifiers...>;

    //! \copydoc compose_new_extents::extents_type
    using extents_type = typename next_compose_new_extents::extents_type;

    /*! \brief
     * Return the subspan extents while setting offset and strides using a pair of integral numbers as last parameter.
     * \todo reenable constexpr with C++14
     */
    template<class OrgExtents, class ... DynamicExtents>
    static extents_type create_sub_extents(const OrgExtents e,
                                           std::array<ptrdiff_t, OrgExtents::rank()> &strides,
                                           ptrdiff_t *offset,
                                           std::pair<IT, IT> p,
                                           SliceSpecifiers ... s,
                                           DynamicExtents... de)
    {
        strides[sizeof ... (ExtentsNew)] = strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        *offset += p.first*strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...,
                                                            ptrdiff_t(p.second-p.first));
    }
};

/*! \internal \brief
 * Template specialisation to compose new extents when the last slice in the
 * parameter pack is a single integral number.
 */
template<
    ptrdiff_t ... ExtentsNew,
    ptrdiff_t E0,
    ptrdiff_t ... ExtentsOld,
    class IT,
    class ... SliceSpecifiers>
struct compose_new_extents<
    extents<ExtentsNew...>,
    extents<E0, ExtentsOld...>,
    IT,
    SliceSpecifiers...>
{
    //! \copydoc compose_new_extents::next_compose_new_extents
    using next_compose_new_extents = compose_new_extents<
                extents<ExtentsNew...>,
                extents<ExtentsOld...>,
                SliceSpecifiers...>;

    //! \copydoc compose_new_extents::extents_type
    using extents_type = typename next_compose_new_extents::extents_type;

    /*! \brief
     * Return the subspan extents while setting offset and strides using an integral numbers as last parameter.
     * \todo reenable constexpr with C++14
     */
    template <class OrgExtents, class ... DynamicExtents>
    static extents_type create_sub_extents(const OrgExtents e,
                                           std::array<ptrdiff_t, OrgExtents::rank()> &strides,
                                           ptrdiff_t * offset,
                                           const ptrdiff_t v,
                                           SliceSpecifiers... s,
                                           DynamicExtents... de)
    {
        *offset += v*strides[OrgExtents::rank()-sizeof ... (SliceSpecifiers)-1];
        return next_compose_new_extents::create_sub_extents(e, strides, offset, s ..., de ...);
    }
};

/*! \internal \brief
 * Template specialisation to compose new extents when there are no more
 * slice specifiers in the paramter pack describing the subspan slices.
 */
template<ptrdiff_t ... ExtentsNew>
struct compose_new_extents< extents<ExtentsNew...>, extents<> >
{
    //! \copydoc compose_new_extents::extents_type
    using extents_type = extents<ExtentsNew...>;

    /*! \brief
     * Return the deduced subspan extents and setting the dynamic extents.
     */
    template<class OrgExtents, class ... DynamicExtents>
    static constexpr extents_type create_sub_extents(
            const OrgExtents /*e*/,
            std::array<ptrdiff_t, OrgExtents::rank()> /*strides*/,
            ptrdiff_t * /*offset*/,
            DynamicExtents... de)
    {
        return extents_type(de ...);
    }
    /*! \brief
     * Return the deduced subspan extents.
     */
    template<class OrgExtents>
    static constexpr extents_type create_sub_extents(const OrgExtents /*e*/,
                                                     std::array<ptrdiff_t, OrgExtents::rank()> /*strides*/,
                                                     ptrdiff_t * /*offset*/)
    {
        return extents_type();
    }
};

/*! \internal \brief
 * Deduce the extents of a subspan.
 */
template<class Extents, class ... SliceSpecifiers>
struct subspan_deduce_extents
{
    //! \copydoc compose_new_extents::extents_type
    using extents_type = typename compose_new_extents<extents<>, Extents, SliceSpecifiers...>::extents_type;
    //! Define stride stype for subspan, must have as many strides as there are ranks
    using stride_type  = std::array<ptrdiff_t, Extents::rank()>;
    //! create the subspan extents and set its stride and offset
    static constexpr extents_type create_sub_extents(const Extents e,
                                                     stride_type &strides, ptrdiff_t *offset, SliceSpecifiers... s)
    {
        return compose_new_extents<
            extents<>, Extents, SliceSpecifiers...>::create_sub_extents(e, strides, offset, s ...);
    }
};

}   // namespace detail

/*! \brief
 * Create a new basic_mdspan by slicing a given basic_mdspan.
 *
 * All template parameters but SliceSpecifiers are deduced from the basic_mdspan to be sliced.
 *
 * Slicing is specified with a single number, a pair of numbers or all_type.
 * A pair of numbers denotes the begin and end index of the slice.
 * Single number is equivalent to a pair of equal numbers,
 * all_type is equivalent in behaviour to a pair (0,extent-1)
 *
 * \note Layout of basic_mdspan slices is always strided.
 *
 * \tparam ElementType    type of element of the basic_mdspan to be sliced
 * \tparam Extents        extents of the basic_mdspan to be sliced
 * \tparam LayoutPolicy   layout policity of the basic_mdspan to be sliced
 * \tparam AccessorPolicy of the basic_mdspan to be sliced
 * \note the new subspan will retain only offset_policy from the basic_mdspan AccessorPolicy
 * \tparam SliceSpecifiers Parameter pack determining how the array shall be sliced
 *
 * \param[in] src    the span to be sliced
 * \param[in] slices define how the subspan shall be created, either pairs of
 *                   integral numbers, one integral number or all_type
 * \returns basic_mdspan that provides a sliced subspan view with a strided layout on the input basic_mdspan.
 */
template<class ElementType,
         class Extents,
         class LayoutPolicy,
         class AccessorPolicy,
         class ... SliceSpecifiers>
basic_mdspan< ElementType,
              typename detail::subspan_deduce_extents<Extents, SliceSpecifiers...>::extents_type,
              layout_stride,
              typename AccessorPolicy::offset_policy >
subspan(const basic_mdspan<ElementType, Extents, LayoutPolicy, AccessorPolicy> &src,
        SliceSpecifiers ... slices) noexcept
{
    // the extents of the subspan
    using sub_extents_type = typename detail::subspan_deduce_extents<Extents, SliceSpecifiers...>::extents_type;

    // the accessor policy of the subspan
    using sub_mdspan_type = basic_mdspan<
                ElementType, sub_extents_type, layout_stride, typename AccessorPolicy::offset_policy>;

    // copying the strides of the input basic_mdspan into an array
    std::array<ptrdiff_t, Extents::rank()> strides;
    for (size_t r = 0; r < Extents::rank(); r++)
    {
        strides[r] = src.stride(r);
    }

    // the memory offset of the first element in the subspan
    ptrdiff_t        offset = 0;
    // deduce the extents of the subspan from the input extents, strides and slice parameters
    sub_extents_type sub_extents = detail::subspan_deduce_extents<
            Extents, SliceSpecifiers...>::create_sub_extents(src.extents(), strides, &offset, slices ...);

    // copy the strides of the subspan
    std::array<ptrdiff_t, sub_extents_type::rank()> sub_strides;
    for (size_t r = 0; r < sub_extents_type::rank(); r++)
    {
        sub_strides[r] = strides[r];
    }

    // offset the pointer to the data by the subspan offset.
    typename AccessorPolicy::offset_policy::pointer ptr = src.accessor().offset(src.data(), offset);
    return sub_mdspan_type(ptr, typename sub_mdspan_type::mapping_type(sub_extents, sub_strides));
}

}      // namespace gmx

#endif /* end of include guard: MDSPAN_SUBSPAN_H */
