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
 * Implements gmx::Selection.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <statutil.h>
#include <string2.h>
#include <xvgr.h>

#include "gromacs/selection/position.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selvalue.h"

#include "selelem.h"

namespace gmx
{

Selection::Selection(t_selelem *elem, const char *selstr)
{
    // TODO: This is not exception-safe if any called function throws.
    _sel.name = strdup(elem->name);
    _sel.selstr = strdup(selstr);
    gmx_ana_pos_clear(&_sel.p);

    if (elem->child->type == SEL_CONST)
    {
        gmx_ana_pos_copy(&_sel.p, elem->child->v.u.p, true);
        _sel.bDynamic = false;
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
Selection::printInfo(FILE *fp) const
{
    fprintf(fp, "\"%s\" (%d position%s, %d atom%s%s)", _sel.name,
            _sel.p.nr,     _sel.p.nr     == 1 ? "" : "s",
            _sel.g->isize, _sel.g->isize == 1 ? "" : "s",
            _sel.bDynamic ? ", dynamic" : "");
    fprintf(fp, "\n");
}


bool
Selection::initCoveredFraction(e_coverfrac_t type)
{
    gmx_ana_selection_t *sel = &_sel;

    sel->cfractype = type;
    if (type == CFRAC_NONE || !sel->selelem)
    {
        sel->bCFracDyn = false;
    }
    else if (!_gmx_selelem_can_estimate_cover(sel->selelem))
    {
        sel->cfractype = CFRAC_NONE;
        sel->bCFracDyn = false;
    }
    else
    {
        sel->bCFracDyn = true;
    }
    sel->cfrac     = sel->bCFracDyn ? 0.0 : 1.0;
    sel->avecfrac  = sel->cfrac;
    return type == CFRAC_NONE || sel->cfractype != CFRAC_NONE;
}


void
Selection::printDebugInfo(FILE *fp, int nmaxind) const
{
    fprintf(fp, "  ");
    printInfo(fp);
    fprintf(fp, "    ");
    gmx_ana_index_dump(fp, _sel.g, -1, nmaxind);

    fprintf(fp, "    Block (size=%d):", _sel.p.m.mapb.nr);
    if (!_sel.p.m.mapb.index)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        int n = _sel.p.m.mapb.nr;
        if (nmaxind >= 0 && n > nmaxind)
            n = nmaxind;
        for (int i = 0; i <= n; ++i)
            fprintf(fp, " %d", _sel.p.m.mapb.index[i]);
        if (n < _sel.p.m.mapb.nr)
            fprintf(fp, " ...");
    }
    fprintf(fp, "\n");

    int n = _sel.p.m.nr;
    if (nmaxind >= 0 && n > nmaxind)
        n = nmaxind;
    fprintf(fp, "    RefId:");
    if (!_sel.p.m.refid)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        for (int i = 0; i < n; ++i)
            fprintf(fp, " %d", _sel.p.m.refid[i]);
        if (n < _sel.p.m.nr)
            fprintf(fp, " ...");
    }
    fprintf(fp, "\n");

    fprintf(fp, "    MapId:");
    if (!_sel.p.m.mapid)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        for (int i = 0; i < n; ++i)
            fprintf(fp, " %d", _sel.p.m.mapid[i]);
        if (n < _sel.p.m.nr)
            fprintf(fp, " ...");
    }
    fprintf(fp, "\n");
}


void
Selection::initializeMassesAndCharges(const t_topology *top)
{
    snew(_sel.orgm, posCount());
    snew(_sel.orgq, posCount());
    for (int b = 0; b < posCount(); ++b)
    {
        _sel.orgq[b] = 0;
        if (top)
        {
            _sel.orgm[b] = 0;
            for (int i = _sel.p.m.mapb.index[b];
                     i < _sel.p.m.mapb.index[b+1];
                     ++i)
            {
                int index = _sel.p.g->index[i];
                _sel.orgm[b] += top->atoms.atom[index].m;
                _sel.orgq[b] += top->atoms.atom[index].q;
            }
        }
        else
        {
            _sel.orgm[b] = 1;
        }
    }
    if (isDynamic() && !hasFlag(efDynamicMask))
    {
        snew(_sel.m, posCount());
        snew(_sel.q, posCount());
        for (int b = 0; b < posCount(); ++b)
        {
            _sel.m[b] = _sel.orgm[b];
            _sel.q[b] = _sel.orgq[b];
        }
    }
    else
    {
        _sel.m = _sel.orgm;
        _sel.q = _sel.orgq;
    }
}

} // namespace gmx
