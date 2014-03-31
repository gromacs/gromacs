/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2012,2013,2014, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements functions in selvalue.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

/*!
 * \param[out] val  Output structure
 *
 * The type of \p val is not touched.
 * Any contents of \p val are discarded without freeing.
 */
void
_gmx_selvalue_clear(gmx_ana_selvalue_t *val)
{
    val->nr     = 0;
    val->u.ptr  = NULL;
    val->nalloc = 0;
}

/*!
 * \param[in,out] val  Value structure to allocate.
 * \param[in]     n    Maximum number of values needed.
 * \returns       Zero on success.
 *
 * Reserves memory for the values within \p val to store at least \p n values,
 * of the type specified in the \p val structure.
 *
 * If the type is \ref POS_VALUE or \ref GROUP_VALUE, memory is reserved for
 * the data structures, but no memory is reserved inside these newly allocated
 * data structures.
 * Similarly, for \ref STR_VALUE values, the pointers are set to NULL.
 * For other values, the memory is uninitialized.
 */
int
_gmx_selvalue_reserve(gmx_ana_selvalue_t *val, int n)
{
    int  i;

    if (val->nalloc == -1)
    {
        return 0;
    }

    if (!val->u.ptr || val->nalloc < n)
    {
        switch (val->type)
        {
            case INT_VALUE:   srenew(val->u.i, n); break;
            case REAL_VALUE:  srenew(val->u.r, n); break;
            case STR_VALUE:
                srenew(val->u.s, n);
                for (i = val->nalloc; i < n; ++i)
                {
                    val->u.s[i] = NULL;
                }
                break;
            case POS_VALUE:
                GMX_RELEASE_ASSERT(val->u.ptr == NULL,
                                   "Reallocation of position values not supported");
                val->u.p = new gmx_ana_pos_t[n];
                break;
            case GROUP_VALUE:
                srenew(val->u.g, n);
                for (i = val->nalloc; i < n; ++i)
                {
                    gmx_ana_index_clear(&val->u.g[i]);
                }
                break;
            case NO_VALUE:    break;
        }
        val->nalloc = n;
    }
    return 0;
}

/*!
 * \param[in,out] val    Value structure to allocate.
 * \param[in]     ptr    Pointer where the values should be stored.
 * \returns       Zero on success.
 *
 * Automatic memory management is disabled for \p ptr, unless \p ptr is NULL.
 */
int
_gmx_selvalue_setstore(gmx_ana_selvalue_t *val, void *ptr)
{
    val->u.ptr  = ptr;
    val->nalloc = (ptr ? -1 : 0);
    return 0;
}

/*!
 * \param[in,out] val    Value structure to allocate.
 * \param[in]     ptr    Pointer where the values should be stored.
 * \param[in]     nalloc Number of values allocated for \p ptr.
 * \returns       Zero on success.
 */
int
_gmx_selvalue_setstore_alloc(gmx_ana_selvalue_t *val, void *ptr, int nalloc)
{
    val->u.ptr  = ptr;
    val->nalloc = nalloc;
    return 0;
}
