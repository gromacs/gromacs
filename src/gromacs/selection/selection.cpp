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
 * Implements classes in selection.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#include "selection.h"

#include "position.h"
#include "selelem.h"
#include "selvalue.h"

namespace gmx
{

namespace internal
{

SelectionData::SelectionData(SelectionTreeElement *elem,
                             const char *selstr)
    : name_(elem->name()), selectionText_(selstr),
      rootElement_(*elem), coveredFractionType_(CFRAC_NONE),
      coveredFraction_(1.0), averageCoveredFraction_(1.0),
      bDynamic_(false), bDynamicCoveredFraction_(false)
{
    gmx_ana_pos_clear(&rawPositions_);

    if (elem->child->type == SEL_CONST)
    {
        // TODO: This is not exception-safe if any called function throws.
        gmx_ana_pos_copy(&rawPositions_, elem->child->v.u.p, true);
    }
    else
    {
        SelectionTreeElementPointer child = elem->child;
        child->flags     &= ~SEL_ALLOCVAL;
        _gmx_selvalue_setstore(&child->v, &rawPositions_);
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
        bDynamic_ = (child->child->flags & SEL_DYNAMIC);
    }
    initCoveredFraction(CFRAC_NONE);
}


SelectionData::~SelectionData()
{
    gmx_ana_pos_deinit(&rawPositions_);
}


bool
SelectionData::initCoveredFraction(e_coverfrac_t type)
{
    coveredFractionType_ = type;
    if (type == CFRAC_NONE)
    {
        bDynamicCoveredFraction_ = false;
    }
    else if (!_gmx_selelem_can_estimate_cover(rootElement()))
    {
        coveredFractionType_ = CFRAC_NONE;
        bDynamicCoveredFraction_ = false;
    }
    else
    {
        bDynamicCoveredFraction_ = true;
    }
    coveredFraction_ = bDynamicCoveredFraction_ ? 0.0 : 1.0;
    averageCoveredFraction_ = coveredFraction_;
    return type == CFRAC_NONE || coveredFractionType_ != CFRAC_NONE;
}


void
SelectionData::initializeMassesAndCharges(const t_topology *top)
{
    posInfo_.reserve(posCount());
    for (int b = 0; b < posCount(); ++b)
    {
        real mass   = 1.0;
        real charge = 0.0;
        if (top != NULL)
        {
            mass = 0.0;
            for (int i = rawPositions_.m.mapb.index[b];
                 i < rawPositions_.m.mapb.index[b+1];
                 ++i)
            {
                int index = rawPositions_.g->index[i];
                mass   += top->atoms.atom[index].m;
                charge += top->atoms.atom[index].q;
            }
        }
        posInfo_.push_back(PositionInfo(mass, charge));
    }
    if (isDynamic() && !hasFlag(efSelection_DynamicMask))
    {
        originalPosInfo_ = posInfo_;
    }
}


void
SelectionData::refreshMassesAndCharges()
{
    if (!originalPosInfo_.empty())
    {
        posInfo_.clear();
        for (int i = 0; i < posCount(); ++i)
        {
            int refid  = rawPositions_.m.refid[i];
            posInfo_.push_back(originalPosInfo_[refid]);
        }
    }
}


void
SelectionData::updateCoveredFractionForFrame()
{
    if (isCoveredFractionDynamic())
    {
        real cfrac = _gmx_selelem_estimate_coverfrac(rootElement());
        coveredFraction_ = cfrac;
        averageCoveredFraction_ += cfrac;
    }
}


void
SelectionData::computeAverageCoveredFraction(int nframes)
{
    if (isCoveredFractionDynamic() && nframes > 0)
    {
        averageCoveredFraction_ /= nframes;
    }
}


void
SelectionData::restoreOriginalPositions()
{
    if (isDynamic())
    {
        gmx_ana_pos_t &p = rawPositions_;
        gmx_ana_index_copy(p.g, rootElement().v.u.g, false);
        p.g->name = NULL;
        gmx_ana_indexmap_update(&p.m, p.g, hasFlag(gmx::efSelection_DynamicMask));
        p.nr = p.m.nr;
        refreshMassesAndCharges();
    }
}

} // namespace internal


void
Selection::printInfo(FILE *fp) const
{
    fprintf(fp, "\"%s\" (%d position%s, %d atom%s%s)", name(),
            posCount(),  posCount()  == 1 ? "" : "s",
            atomCount(), atomCount() == 1 ? "" : "s",
            isDynamic() ? ", dynamic" : "");
    fprintf(fp, "\n");
}


void
Selection::printDebugInfo(FILE *fp, int nmaxind) const
{
    const gmx_ana_pos_t &p = data().rawPositions_;

    fprintf(fp, "  ");
    printInfo(fp);
    fprintf(fp, "    ");
    gmx_ana_index_dump(fp, p.g, -1, nmaxind);

    fprintf(fp, "    Block (size=%d):", p.m.mapb.nr);
    if (!p.m.mapb.index)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        int n = p.m.mapb.nr;
        if (nmaxind >= 0 && n > nmaxind)
        {
            n = nmaxind;
        }
        for (int i = 0; i <= n; ++i)
        {
            fprintf(fp, " %d", p.m.mapb.index[i]);
        }
        if (n < p.m.mapb.nr)
        {
            fprintf(fp, " ...");
        }
    }
    fprintf(fp, "\n");

    int n = posCount();
    if (nmaxind >= 0 && n > nmaxind)
    {
        n = nmaxind;
    }
    fprintf(fp, "    RefId:");
    if (!p.m.refid)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            fprintf(fp, " %d", p.m.refid[i]);
        }
        if (n < posCount())
        {
            fprintf(fp, " ...");
        }
    }
    fprintf(fp, "\n");

    fprintf(fp, "    MapId:");
    if (!p.m.mapid)
    {
        fprintf(fp, " (null)");
    }
    else
    {
        for (int i = 0; i < n; ++i)
        {
            fprintf(fp, " %d", p.m.mapid[i]);
        }
        if (n < posCount())
        {
            fprintf(fp, " ...");
        }
    }
    fprintf(fp, "\n");
}

} // namespace gmx
