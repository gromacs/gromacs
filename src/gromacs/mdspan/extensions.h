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
/*! \libinternal
 * \file
 * \brief GROMACS extensions to mdspan.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_mdspan
 */

#ifndef GMX_MDSPAN_EXTENSIONS_H_
#define GMX_MDSPAN_EXTENSIONS_H_

#include "extents.h"

namespace gmx
{

/*! \brief
 * Free begin function addressing memory of a contiguously laid out basic_mdspan.
 *
 * \note Changing the elements that basic_mdspan views does not change
 *       the view itself, so a single begin that takes a const view suffices.
 */
template <class BasicMdspan>
constexpr typename std::enable_if<BasicMdspan::is_always_contiguous(),
                                  typename BasicMdspan::pointer>::type
begin(const BasicMdspan &basicMdspan)
{
    return basicMdspan.data();
}

/*! \brief
 * Free end function addressing memory of a contiguously laid out basic_mdspan.
 *
 * \note Changing the elements that basic_mdspan views does not change
 *       the view itself, so a single end that takes a const view suffices.
 */
template <class BasicMdspan>
constexpr typename std::enable_if<BasicMdspan::is_always_contiguous(),
                                  typename BasicMdspan::pointer>::type
end(const BasicMdspan &basicMdspan)
{
    return basicMdspan.data() + basicMdspan.mapping().required_span_size();
}

/*! \internal \brief
 * Drop front extent from multidimensional extents.
 * \tparam Extents mdspan
 * \tparam ExtentsAreDynamic type enabling compile time dropping of front extent
 *                           if extents after dropping the front one are purely static
 */
template<class Extents, typename ExtentsAreDynamic = void>
struct frontExtentDropper;

/*! \internal \brief Specialisation for dynamic extents after dropping the front extent.
 *
 * Keeps track of the dynamic extents that are only known during run time that
 * prevents this functionality to be constexpr as preferred for the static extents.
 *
 * Drops the first dimension of the template parameter extents, by matching
 * template parameters accordingly.
 *
 * \tparam E0      The first extent - to be dropped
 * \tparam Extents The leftover extents.
 */
template <ptrdiff_t E0, ptrdiff_t... Extents>
struct frontExtentDropper<extents<E0, Extents...>,
                          typename std::enable_if<extents<Extents...>::rank_dynamic() != 0>::type>
{
    //! The type of the extents with the first one dropped.
    using extents_type = extents<Extents...>;

    /*! \brief returns extents where the first extent is dropped,
     * maintaining the information about the dynamic extents.
     */
    static extents_type dropFirstExtent(const extents<E0, Extents...> fullExtents)
    {
        std::array<ptrdiff_t, extents_type::rank_dynamic()> dynamicExtents;
        size_t currentDynamicExtent = 0;
        for (size_t r = 0; r < extents_type::rank(); ++r)
        {
            // dropped the first extent of fullExtents in the result
            if (fullExtents.static_extent(r + 1) == dynamic_extent)
            {
                dynamicExtents[currentDynamicExtent] = fullExtents.extent(r + 1);
                ++currentDynamicExtent;
            }
        }
        return {dynamicExtents};
    }
};

/*! \internal \brief Specialisation for purely static extents after dropping the front extent.
 *
 * Allows compile time evaluation of extents.
 *
 * Drops the first dimension of the template parameter extents, by matching
 * template parameters accordingly.
 *
 * \tparam E0      The first extent - to be dropped
 * \tparam Extents The leftover extents.
 */
template <ptrdiff_t E0, ptrdiff_t... Extents>
struct frontExtentDropper<extents<E0, Extents...>, typename std::enable_if<extents<Extents...>::rank_dynamic() == 0>::type>
{
    using extents_type = extents<Extents...>;
    static constexpr extents_type dropFirstExtent(const extents<E0, Extents...> /*extents*/)
    {
        return {};
    }
};

}      // namespace gmx

#endif // GMX_MDSPAN_EXTENSIONS_H_
