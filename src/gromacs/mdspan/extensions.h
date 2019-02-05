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
 * \inlibraryapi
 * \ingroup module_mdspan
 */

#ifndef GMX_MDSPAN_EXTENSIONS_H_
#define GMX_MDSPAN_EXTENSIONS_H_

#include "gromacs/mdspan/mdspan.h"

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

}      // namespace gmx

#endif // GMX_MDSPAN_EXTENSIONS_H_
