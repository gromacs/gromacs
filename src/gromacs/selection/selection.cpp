/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements classes in selection.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "selection.h"

#include <string>

#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/position.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "selelem.h"
#include "selvalue.h"

namespace gmx
{

namespace internal
{

/********************************************************************
 * SelectionData
 */

SelectionData::SelectionData(SelectionTreeElement *elem,
                             const char           *selstr)
    : name_(elem->name()), selectionText_(selstr),
      rootElement_(*elem), coveredFractionType_(CFRAC_NONE),
      coveredFraction_(1.0), averageCoveredFraction_(1.0),
      bDynamic_(false), bDynamicCoveredFraction_(false)
{
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
            if (!child)
            {
                break;
            }
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
        if (child)
        {
            /* For variable references, we should skip the
             * SEL_SUBEXPRREF and SEL_SUBEXPR elements. */
            if (child->type == SEL_SUBEXPRREF)
            {
                child = child->child->child;
            }
            bDynamic_ = (child->child->flags & SEL_DYNAMIC);
        }
    }
    initCoveredFraction(CFRAC_NONE);
}


SelectionData::~SelectionData()
{
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
        coveredFractionType_     = CFRAC_NONE;
        bDynamicCoveredFraction_ = false;
    }
    else
    {
        bDynamicCoveredFraction_ = true;
    }
    coveredFraction_        = bDynamicCoveredFraction_ ? 0.0 : 1.0;
    averageCoveredFraction_ = coveredFraction_;
    return type == CFRAC_NONE || coveredFractionType_ != CFRAC_NONE;
}

namespace
{

/*! \brief
 * Helper function to compute total masses and charges for positions.
 *
 * \param[in]  top     Topology to take atom masses from.
 * \param[in]  pos     Positions to compute masses and charges for.
 * \param[out] masses  Output masses.
 * \param[out] charges Output charges.
 *
 * Does not throw if enough space has been reserved for the output vectors.
 */
void computeMassesAndCharges(const t_topology *top, const gmx_ana_pos_t &pos,
                             std::vector<real> *masses,
                             std::vector<real> *charges)
{
    GMX_ASSERT(top != NULL, "Should not have been called with NULL topology");
    masses->clear();
    charges->clear();
    for (int b = 0; b < pos.count(); ++b)
    {
        real mass   = 0.0;
        real charge = 0.0;
        for (int i = pos.m.mapb.index[b]; i < pos.m.mapb.index[b+1]; ++i)
        {
            const int index = pos.m.mapb.a[i];
            mass   += top->atoms.atom[index].m;
            charge += top->atoms.atom[index].q;
        }
        masses->push_back(mass);
        charges->push_back(charge);
    }
}

}       // namespace

void
SelectionData::refreshName()
{
    rootElement_.fillNameIfMissing(selectionText_.c_str());
    name_ = rootElement_.name();
}

void
SelectionData::initializeMassesAndCharges(const t_topology *top)
{
    GMX_ASSERT(posMass_.empty() && posCharge_.empty(),
               "Should not be called more than once");
    posMass_.reserve(posCount());
    posCharge_.reserve(posCount());
    if (top == NULL)
    {
        posMass_.resize(posCount(), 1.0);
        posCharge_.resize(posCount(), 0.0);
    }
    else
    {
        computeMassesAndCharges(top, rawPositions_, &posMass_, &posCharge_);
    }
}


void
SelectionData::refreshMassesAndCharges(const t_topology *top)
{
    if (top != NULL && isDynamic() && !hasFlag(efSelection_DynamicMask))
    {
        computeMassesAndCharges(top, rawPositions_, &posMass_, &posCharge_);
    }
}


void
SelectionData::updateCoveredFractionForFrame()
{
    if (isCoveredFractionDynamic())
    {
        real cfrac = _gmx_selelem_estimate_coverfrac(rootElement());
        coveredFraction_         = cfrac;
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
SelectionData::restoreOriginalPositions(const t_topology *top)
{
    if (isDynamic())
    {
        gmx_ana_pos_t &p = rawPositions_;
        gmx_ana_indexmap_update(&p.m, rootElement().v.u.g,
                                hasFlag(gmx::efSelection_DynamicMask));
        refreshMassesAndCharges(top);
    }
}

}   // namespace internal

/********************************************************************
 * Selection
 */

Selection::operator AnalysisNeighborhoodPositions() const
{
    AnalysisNeighborhoodPositions pos(data().rawPositions_.x,
                                      data().rawPositions_.count());
    if (hasOnlyAtoms())
    {
        pos.exclusionIds(atomIndices());
    }
    return pos;
}


void
Selection::setOriginalId(int i, int id)
{
    data().rawPositions_.m.mapid[i] = id;
    data().rawPositions_.m.orgid[i] = id;
}


int
Selection::initOriginalIdsToGroup(t_topology *top, e_index_t type)
{
    try
    {
        return gmx_ana_indexmap_init_orgid_group(&data().rawPositions_.m, top, type);
    }
    catch (const InconsistentInputError &)
    {
        GMX_ASSERT(type == INDEX_RES || type == INDEX_MOL,
                   "Expected that only grouping by residue/molecule would fail");
        std::string message =
            formatString("Cannot group selection '%s' into %s, because some "
                         "positions have atoms from more than one such group.",
                         name(), type == INDEX_MOL ? "molecules" : "residues");
        GMX_THROW(InconsistentInputError(message));
    }
}


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
    fprintf(fp, "    Group ");
    gmx_ana_index_t g;
    gmx_ana_index_set(&g, p.m.mapb.nra, p.m.mapb.a, 0);
    gmx_ana_index_dump(fp, &g, nmaxind);

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


/********************************************************************
 * SelectionPosition
 */

SelectionPosition::operator AnalysisNeighborhoodPositions() const
{
    AnalysisNeighborhoodPositions pos(sel_->rawPositions_.x,
                                      sel_->rawPositions_.count());
    if (sel_->hasOnlyAtoms())
    {
        // TODO: Move atomIndices() such that it can be reused here as well.
        pos.exclusionIds(constArrayRefFromArray<int>(sel_->rawPositions_.m.mapb.a,
                                                     sel_->rawPositions_.m.mapb.nra));
    }
    return pos.selectSingleFromArray(i_);
}

} // namespace gmx
