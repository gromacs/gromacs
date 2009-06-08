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
 * \brief Implementation of functions in selvalue.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>

#include <indexutil.h>
#include <position.h>
#include <selvalue.h>

/*!
 * \param[in,out] val  Value structure to allocate.
 * \param[in]     n    Maximum number of values needed.
 * \returns       Zero on success.
 *
 * Reserves memory for the values within \p val to store at least \p n values,
 * of the type specified in the \p val structure.
 *
 * If the type is \ref POS_VALUE or \ref GROUP_VALUE, only a single data
 * structure is reserved, independent of the value of \p n.
 * No memory is reserved inside these newly allocated data structures.
 */
int
_gmx_selvalue_reserve(gmx_ana_selvalue_t *val, int n)
{
    int  i;

    if (val->type == POS_VALUE || val->type == GROUP_VALUE)
    {
        n = 1;
    }
    if (!val->u.ptr || val->nalloc < n)
    {
        switch (val->type)
        {
            case INT_VALUE:   srenew(val->u.i, n); break;
            case REAL_VALUE:  srenew(val->u.r, n); break;
            case STR_VALUE:   srenew(val->u.s, n); break;
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
