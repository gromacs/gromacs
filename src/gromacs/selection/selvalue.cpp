/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements functions in selvalue.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>

#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selvalue.h"

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
                srenew(val->u.p, n);
                for (i = val->nalloc; i < n; ++i)
                {
                    gmx_ana_pos_clear(&val->u.p[i]);
                }
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
