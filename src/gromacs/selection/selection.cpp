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
 * Implementation of functions in selection.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <statutil.h>
#include <string2.h>
#include <xvgr.h>

#include "gromacs/selection/selection.h"
#include "position.h"
#include "selelem.h"
#include "selvalue.h"

namespace gmx
{

Selection::Selection(t_selelem *elem, const char *selstr)
{
    _sel.name = strdup(elem->name);
    _sel.selstr = strdup(selstr);
    gmx_ana_pos_clear(&_sel.p);

    if (elem->child->type == SEL_CONST)
    {
        gmx_ana_pos_copy(&_sel.p, elem->child->v.u.p, TRUE);
        _sel.bDynamic = FALSE;
    }
    else
    {
        t_selelem *child;

        child = elem->child;
        child->flags     &= ~SEL_ALLOCVAL;
        _gmx_selvalue_setstore(&child->v, &_sel.p);
        /* We should also skip any modifiers to determine the dynamic
         * status. */
        while (child->type == SEL_MODIFIER)
        {
            child = child->child;
            if (child->type == SEL_SUBEXPRREF)
            {
                child = child->child;
                /* Because most subexpression elements are created
                 * during compilation, we need to check for them
                 * explicitly here.
                 */
                if (child->type == SEL_SUBEXPR)
                {
                    child = child->child;
                }
            }
        }
        /* For variable references, we should skip the
         * SEL_SUBEXPRREF and SEL_SUBEXPR elements. */
        if (child->type == SEL_SUBEXPRREF)
        {
            child = child->child->child;
        }
        _sel.bDynamic = (child->child->flags & SEL_DYNAMIC);
    }
    /* The group will be set after compilation */
    _sel.m        = NULL;
    _sel.q        = NULL;
    _sel.g        = NULL;
    _sel.orgm     = NULL;
    _sel.orgq     = NULL;
    _sel.selelem  = elem;
    initCoveredFraction(CFRAC_NONE);
}

Selection::~Selection()
{
    sfree(_sel.name);
    sfree(_sel.selstr);
    gmx_ana_pos_deinit(&_sel.p);
    if (_sel.m != _sel.orgm)
    {
        sfree(_sel.m);
    }
    if (_sel.q != _sel.orgq)
    {
        sfree(_sel.q);
    }
    sfree(_sel.orgm);
    sfree(_sel.orgq);
}

void
Selection::printInfo() const
{
    fprintf(stderr, "\"%s\" (%d position%s, %d atom%s%s)", _sel.name,
            _sel.p.nr,     _sel.p.nr     == 1 ? "" : "s",
            _sel.g->isize, _sel.g->isize == 1 ? "" : "s",
            _sel.bDynamic ? ", dynamic" : "");
    fprintf(stderr, "\n");
}

/*!
 * \param[in] type Type of covered fraction required.
 * \returns   True if the covered fraction can be calculated for the selection.
 */
bool
Selection::initCoveredFraction(e_coverfrac_t type)
{
    gmx_ana_selection_t *sel = &_sel;

    sel->cfractype = type;
    if (type == CFRAC_NONE || !sel->selelem)
    {
        sel->bCFracDyn = FALSE;
    }
    else if (!_gmx_selelem_can_estimate_cover(sel->selelem))
    {
        sel->cfractype = CFRAC_NONE;
        sel->bCFracDyn = FALSE;
    }
    else
    {
        sel->bCFracDyn = TRUE;
    }
    sel->cfrac     = sel->bCFracDyn ? 0.0 : 1.0;
    sel->avecfrac  = sel->cfrac;
    return type == CFRAC_NONE || sel->cfractype != CFRAC_NONE;
}

} // namespace gmx
