/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements gmx::PositionCalculationCollection and functions in poscalc.h.
 *
 * \todo
 * There is probably some room for optimization in the calculation of
 * positions with bases.
 * In particular, the current implementation may do a lot of unnecessary
 * copying.
 * The interface would need to be changed to make it possible to use the
 * same output positions for several calculations.
 *
 * \todo
 * The current algorithm for setting up base calculations could be improved
 * in cases when there are calculations that cannot use a common base but
 * still overlap partially (e.g., with three calculations A, B, and C
 * such that A could use both B and C as a base, but B and C cannot use the
 * same base).
 * Setting up the bases in an optimal manner in every possible situation can
 * be quite difficult unless several bases are allowed for one calculation,
 * but better heuristics could probably be implemented.
 * For best results, the setup should probably be postponed (at least
 * partially) to gmx_ana_poscalc_init_eval().
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "poscalc.h"

#include <string.h>

#include <algorithm>
#include <vector>

#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "centerofmass.h"

namespace gmx
{

/*! \internal \brief
 * Private implementation class for PositionCalculationCollection.
 *
 * \ingroup module_selection
 */
class PositionCalculationCollection::Impl
{
    public:
        Impl();
        ~Impl();

        /*! \brief
         * Inserts a position calculation structure into its collection.
         *
         * \param pc     Data structure to insert.
         * \param before Data structure before which to insert
         *   (NULL = insert at end).
         *
         * Inserts \p pc to its collection before \p before.
         * If \p before is NULL, \p pc is appended to the list.
         */
        void insertCalculation(gmx_ana_poscalc_t *pc, gmx_ana_poscalc_t *before);
        /*! \brief
         * Removes a position calculation structure from its collection.
         *
         * \param pc    Data structure to remove.
         *
         * Removes \p pc from its collection.
         */
        void removeCalculation(gmx_ana_poscalc_t *pc);

        /*! \copydoc PositionCalculationCollection::createCalculation()
         *
         * This method contains the actual implementation of the similarly
         * named method in the parent class.
         */
        gmx_ana_poscalc_t *createCalculation(e_poscalc_t type, int flags);

        /*! \brief
         * Maps given topology indices into frame indices.
         *
         * Only one position calculation at a time needs to access this (and
         * there are also other thread-unsafe constructs here), so a temporary
         * array is used to avoid repeated memory allocation.
         */
        ArrayRef<const int> getFrameIndices(int size, int index[])
        {
            if (mapToFrameAtoms_.empty())
            {
                return constArrayRefFromArray(index, size);
            }
            tmpFrameAtoms_.resize(size);
            for (int i = 0; i < size; ++i)
            {
                const int ii = index[i];
                GMX_ASSERT(ii >= 0 && ii <= static_cast<int>(mapToFrameAtoms_.size())
                           && mapToFrameAtoms_[ii] != -1,
                           "Invalid input atom index");
                tmpFrameAtoms_[i] = mapToFrameAtoms_[ii];
            }
            return tmpFrameAtoms_;
        }

        /*! \brief
         * Topology data.
         *
         * Can be NULL if none of the calculations require topology data or if
         * setTopology() has not been called.
         */
        const gmx_mtop_t         *top_;
        //! Pointer to the first data structure.
        gmx_ana_poscalc_t        *first_;
        //! Pointer to the last data structure.
        gmx_ana_poscalc_t        *last_;
        //! Whether the collection has been initialized for evaluation.
        bool                      bInit_;
        //! Mapping from topology atoms to frame atoms (one index for each topology atom).
        std::vector<int>          mapToFrameAtoms_;
        //! Working array for updating positions.
        std::vector<int>          tmpFrameAtoms_;
};

} // namespace gmx

/*! \internal \brief
 * Data structure for position calculation.
 */
struct gmx_ana_poscalc_t
{
    /*! \brief
     * Type of calculation.
     *
     * This field may differ from the type requested by the user, because
     * it is changed internally to the most effective calculation.
     * For example, if the user requests a COM calculation for residues
     * consisting of single atoms, it is simply set to POS_ATOM.
     * To provide a consistent interface to the user, the field \p itype
     * should be used when information should be given out.
     */
    e_poscalc_t               type;
    /*! \brief
     * Flags for calculation options.
     *
     * See \ref poscalc_flags "documentation of the flags".
     */
    int                       flags;

    /*! \brief
     * Type for the created indices.
     *
     * This field always agrees with the type that the user requested, but
     * may differ from \p type.
     */
    e_index_t                 itype;
    /*! \brief
     * Block data for the calculation.
     */
    t_blocka                  b;
    /*! \brief
     * Mapping from the blocks to the blocks of \p sbase.
     *
     * If \p sbase is NULL, this field is also.
     */
    int                      *baseid;
    /*! \brief
     * Maximum evaluation group.
     */
    gmx_ana_index_t           gmax;

    /** Position storage for calculations that are used as a base. */
    gmx_ana_pos_t            *p;

    /** true if the positions have been evaluated for the current frame. */
    bool                      bEval;
    /*! \brief
     * Base position data for this calculation.
     *
     * If not NULL, the centers required by this calculation have
     * already been calculated in \p sbase.
     * The structure pointed by \p sbase is always a static calculation.
     */
    gmx_ana_poscalc_t                        *sbase;
    /** Next structure in the linked list of calculations. */
    gmx_ana_poscalc_t                        *next;
    /** Previous structure in the linked list of calculations. */
    gmx_ana_poscalc_t                        *prev;
    /** Number of references to this structure. */
    int                                       refcount;
    /** Collection this calculation belongs to. */
    gmx::PositionCalculationCollection::Impl *coll;
};

const char * const gmx::PositionCalculationCollection::typeEnumValues[] = {
    "atom",
    "res_com",       "res_cog",
    "mol_com",       "mol_cog",
    "whole_res_com", "whole_res_cog",
    "whole_mol_com", "whole_mol_cog",
    "part_res_com",  "part_res_cog",
    "part_mol_com",  "part_mol_cog",
    "dyn_res_com",   "dyn_res_cog",
    "dyn_mol_com",   "dyn_mol_cog",
    nullptr,
};

/*! \brief
 * Returns the partition type for a given position type.
 *
 * \param [in] type  \c e_poscalc_t value to convert.
 * \returns    Corresponding \c e_indet_t.
 */
static e_index_t
index_type_for_poscalc(e_poscalc_t type)
{
    switch (type)
    {
        case POS_ATOM:    return INDEX_ATOM;
        case POS_RES:     return INDEX_RES;
        case POS_MOL:     return INDEX_MOL;
        case POS_ALL:     return INDEX_ALL;
        case POS_ALL_PBC: return INDEX_ALL;
    }
    return INDEX_UNKNOWN;
}

namespace gmx
{

namespace
{

//! Helper function for determining required topology information.
PositionCalculationCollection::RequiredTopologyInfo
requiredTopologyInfo(e_poscalc_t type, int flags)
{
    if (type != POS_ATOM)
    {
        if ((flags & POS_MASS) || (flags & POS_FORCES))
        {
            return PositionCalculationCollection::RequiredTopologyInfo::TopologyAndMasses;
        }
        if (type == POS_RES || type == POS_MOL)
        {
            return PositionCalculationCollection::RequiredTopologyInfo::Topology;
        }
    }
    return PositionCalculationCollection::RequiredTopologyInfo::None;
}

}   // namespace

// static
void
PositionCalculationCollection::typeFromEnum(const char *post,
                                            e_poscalc_t *type, int *flags)
{
    if (post[0] == 'a')
    {
        *type   = POS_ATOM;
        *flags &= ~(POS_MASS | POS_COMPLMAX | POS_COMPLWHOLE);
        return;
    }

    /* Process the prefix */
    const char *ptr = post;
    if (post[0] == 'w')
    {
        *flags &= ~POS_COMPLMAX;
        *flags |= POS_COMPLWHOLE;
        ptr     = post + 6;
    }
    else if (post[0] == 'p')
    {
        *flags &= ~POS_COMPLWHOLE;
        *flags |= POS_COMPLMAX;
        ptr     = post + 5;
    }
    else if (post[0] == 'd')
    {
        *flags &= ~(POS_COMPLMAX | POS_COMPLWHOLE);
        ptr     = post + 4;
    }

    if (ptr[0] == 'r')
    {
        *type = POS_RES;
    }
    else if (ptr[0] == 'm')
    {
        *type = POS_MOL;
    }
    else
    {
        GMX_THROW(InternalError("Unknown position calculation type"));
    }
    if (strlen(ptr) < 7)
    {
        GMX_THROW(InternalError("Unknown position calculation type"));
    }
    if (ptr[6] == 'm')
    {
        *flags |= POS_MASS;
    }
    else if (ptr[6] == 'g')
    {
        *flags &= ~POS_MASS;
    }
    else
    {
        GMX_THROW(InternalError("Unknown position calculation type"));
    }
}

// static
PositionCalculationCollection::RequiredTopologyInfo
PositionCalculationCollection::requiredTopologyInfoForType(const char *post,
                                                           bool        forces)
{
    e_poscalc_t  type;
    int          flags = (forces ? POS_FORCES : 0);
    PositionCalculationCollection::typeFromEnum(post, &type, &flags);
    return requiredTopologyInfo(type, flags);
}

/********************************************************************
 * PositionCalculationCollection::Impl
 */

PositionCalculationCollection::Impl::Impl()
    : top_(nullptr), first_(nullptr), last_(nullptr), bInit_(false)
{
}

PositionCalculationCollection::Impl::~Impl()
{
    // Loop backwards, because there can be internal references in that are
    // correctly handled by this direction.
    while (last_ != nullptr)
    {
        GMX_ASSERT(last_->refcount == 1,
                   "Dangling references to position calculations");
        gmx_ana_poscalc_free(last_);
    }
}

void
PositionCalculationCollection::Impl::insertCalculation(gmx_ana_poscalc_t *pc,
                                                       gmx_ana_poscalc_t *before)
{
    GMX_RELEASE_ASSERT(pc->coll == this, "Inconsistent collections");
    if (before == nullptr)
    {
        pc->next = nullptr;
        pc->prev = last_;
        if (last_ != nullptr)
        {
            last_->next = pc;
        }
        last_ = pc;
    }
    else
    {
        pc->prev     = before->prev;
        pc->next     = before;
        if (before->prev)
        {
            before->prev->next = pc;
        }
        before->prev = pc;
    }
    if (pc->prev == nullptr)
    {
        first_ = pc;
    }
}

void
PositionCalculationCollection::Impl::removeCalculation(gmx_ana_poscalc_t *pc)
{
    GMX_RELEASE_ASSERT(pc->coll == this, "Inconsistent collections");
    if (pc->prev != nullptr)
    {
        pc->prev->next = pc->next;
    }
    else if (pc == first_)
    {
        first_ = pc->next;
    }
    if (pc->next != nullptr)
    {
        pc->next->prev = pc->prev;
    }
    else if (pc == last_)
    {
        last_ = pc->prev;
    }
    pc->prev = pc->next = nullptr;
}

gmx_ana_poscalc_t *
PositionCalculationCollection::Impl::createCalculation(e_poscalc_t type, int flags)
{
    gmx_ana_poscalc_t *pc;

    snew(pc, 1);
    pc->type     = type;
    pc->itype    = index_type_for_poscalc(type);
    gmx_ana_poscalc_set_flags(pc, flags);
    pc->refcount = 1;
    pc->coll     = this;
    insertCalculation(pc, nullptr);
    return pc;
}


/********************************************************************
 * PositionCalculationCollection
 */

PositionCalculationCollection::PositionCalculationCollection()
    : impl_(new Impl)
{
}

PositionCalculationCollection::~PositionCalculationCollection()
{
}

void
PositionCalculationCollection::setTopology(const gmx_mtop_t *top)
{
    impl_->top_ = top;
}

void
PositionCalculationCollection::printTree(FILE *fp) const
{
    gmx_ana_poscalc_t *pc;
    int                i, j;

    fprintf(fp, "Position calculations:\n");
    i  = 1;
    pc = impl_->first_;
    while (pc)
    {
        fprintf(fp, "%2d ", i);
        switch (pc->type)
        {
            case POS_ATOM:    fprintf(fp, "ATOM");    break;
            case POS_RES:     fprintf(fp, "RES");     break;
            case POS_MOL:     fprintf(fp, "MOL");     break;
            case POS_ALL:     fprintf(fp, "ALL");     break;
            case POS_ALL_PBC: fprintf(fp, "ALL_PBC"); break;
        }
        if (pc->itype != index_type_for_poscalc(pc->type))
        {
            fprintf(fp, " (");
            switch (pc->itype)
            {
                case INDEX_UNKNOWN: fprintf(fp, "???");  break;
                case INDEX_ATOM:    fprintf(fp, "ATOM"); break;
                case INDEX_RES:     fprintf(fp, "RES");  break;
                case INDEX_MOL:     fprintf(fp, "MOL");  break;
                case INDEX_ALL:     fprintf(fp, "ALL");  break;
            }
            fprintf(fp, ")");
        }
        fprintf(fp, " flg=");
        if (pc->flags & POS_MASS)
        {
            fprintf(fp, "M");
        }
        if (pc->flags & POS_DYNAMIC)
        {
            fprintf(fp, "D");
        }
        if (pc->flags & POS_MASKONLY)
        {
            fprintf(fp, "A");
        }
        if (pc->flags & POS_COMPLMAX)
        {
            fprintf(fp, "Cm");
        }
        if (pc->flags & POS_COMPLWHOLE)
        {
            fprintf(fp, "Cw");
        }
        if (!pc->flags)
        {
            fprintf(fp, "0");
        }
        fprintf(fp, " nr=%d nra=%d", pc->b.nr, pc->b.nra);
        fprintf(fp, " refc=%d", pc->refcount);
        fprintf(fp, "\n");
        if (pc->gmax.nalloc_index > 0)
        {
            fprintf(fp, "   Group: ");
            if (pc->gmax.isize > 20)
            {
                fprintf(fp, " %d atoms", pc->gmax.isize);
            }
            else
            {
                for (j = 0; j < pc->gmax.isize; ++j)
                {
                    fprintf(fp, " %d", pc->gmax.index[j] + 1);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->b.nalloc_a > 0)
        {
            fprintf(fp, "   Atoms: ");
            if (pc->b.nra > 20)
            {
                fprintf(fp, " %d atoms", pc->b.nra);
            }
            else
            {
                for (j = 0; j < pc->b.nra; ++j)
                {
                    fprintf(fp, " %d", pc->b.a[j] + 1);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->b.nalloc_index > 0)
        {
            fprintf(fp, "   Blocks:");
            if (pc->b.nr > 20)
            {
                fprintf(fp, " %d pcs", pc->b.nr);
            }
            else
            {
                for (j = 0; j <= pc->b.nr; ++j)
                {
                    fprintf(fp, " %d", pc->b.index[j]);
                }
            }
            fprintf(fp, "\n");
        }
        if (pc->sbase)
        {
            gmx_ana_poscalc_t *base;

            fprintf(fp, "   Base: ");
            j    = 1;
            base = impl_->first_;
            while (base && base != pc->sbase)
            {
                ++j;
                base = base->next;
            }
            fprintf(fp, "%d", j);
            if (pc->baseid && pc->b.nr <= 20)
            {
                fprintf(fp, " id:");
                for (j = 0; j < pc->b.nr; ++j)
                {
                    fprintf(fp, " %d", pc->baseid[j]+1);
                }
            }
            fprintf(fp, "\n");
        }
        ++i;
        pc = pc->next;
    }
}

gmx_ana_poscalc_t *
PositionCalculationCollection::createCalculation(e_poscalc_t type, int flags)
{
    return impl_->createCalculation(type, flags);
}

gmx_ana_poscalc_t *
PositionCalculationCollection::createCalculationFromEnum(const char *post, int flags)
{
    e_poscalc_t  type;
    int          cflags = flags;
    typeFromEnum(post, &type, &cflags);
    return impl_->createCalculation(type, cflags);
}

void PositionCalculationCollection::getRequiredAtoms(gmx_ana_index_t *out) const
{
    gmx_ana_poscalc_t *pc     = impl_->first_;
    while (pc)
    {
        // Calculations with a base just copy positions from the base, so
        // those do not need to be considered in the check.
        if (!pc->sbase)
        {
            gmx_ana_index_t g;
            gmx_ana_index_set(&g, pc->b.nra, pc->b.a, 0);
            gmx_ana_index_union_unsorted(out, out, &g);
        }
        pc = pc->next;
    }
}

void PositionCalculationCollection::initEvaluation()
{
    if (impl_->bInit_)
    {
        return;
    }
    gmx_ana_poscalc_t *pc = impl_->first_;
    while (pc)
    {
        /* Initialize position storage for base calculations */
        if (pc->p)
        {
            gmx_ana_poscalc_init_pos(pc, pc->p);
        }
        /* Construct the mapping of the base positions */
        if (pc->sbase)
        {
            int                     bi, bj;

            snew(pc->baseid, pc->b.nr);
            for (bi = bj = 0; bi < pc->b.nr; ++bi, ++bj)
            {
                while (pc->sbase->b.a[pc->sbase->b.index[bj]] != pc->b.a[pc->b.index[bi]])
                {
                    ++bj;
                }
                pc->baseid[bi] = bj;
            }
        }
        /* Free the block data for dynamic calculations */
        if (pc->flags & POS_DYNAMIC)
        {
            if (pc->b.nalloc_index > 0)
            {
                sfree(pc->b.index);
                pc->b.nalloc_index = 0;
            }
            if (pc->b.nalloc_a > 0)
            {
                sfree(pc->b.a);
                pc->b.nalloc_a = 0;
            }
        }
        pc = pc->next;
    }
    impl_->bInit_ = true;
}

void PositionCalculationCollection::initFrame(const t_trxframe *fr)
{
    if (!impl_->bInit_)
    {
        initEvaluation();
    }
    /* Clear the evaluation flags */
    gmx_ana_poscalc_t *pc = impl_->first_;
    while (pc)
    {
        pc->bEval = false;
        pc        = pc->next;
    }
    if (fr->bIndex && fr->natoms > 0)
    {
        const int highestAtom = *std::max_element(fr->index, fr->index + fr->natoms);
        impl_->mapToFrameAtoms_.resize(highestAtom + 1);
        std::fill(impl_->mapToFrameAtoms_.begin(), impl_->mapToFrameAtoms_.end(), -1);
        for (int i = 0; i < fr->natoms; ++i)
        {
            impl_->mapToFrameAtoms_[fr->index[i]] = i;
        }
    }
    else
    {
        impl_->mapToFrameAtoms_.clear();
    }
}

} // namespace gmx

/*! \brief
 * Initializes position calculation using the maximum possible input index.
 *
 * \param[in,out] pc  Position calculation data structure.
 * \param[in]     g   Maximum index group for the calculation.
 * \param[in]     bBase Whether \p pc will be used as a base or not.
 *
 * \p bBase affects on how the \p pc->gmax field is initialized.
 */
static void
set_poscalc_maxindex(gmx_ana_poscalc_t *pc, gmx_ana_index_t *g, bool bBase)
{
    const gmx_mtop_t *top = pc->coll->top_;
    gmx_ana_index_make_block(&pc->b, top, g, pc->itype, pc->flags & POS_COMPLWHOLE);
    /* Set the type to POS_ATOM if the calculation in fact is such. */
    if (pc->b.nr == pc->b.nra)
    {
        pc->type   = POS_ATOM;
        pc->flags &= ~(POS_MASS | POS_COMPLMAX | POS_COMPLWHOLE);
    }
    /* Set the POS_COMPLWHOLE flag if the calculation in fact always uses
     * complete residues and molecules. */
    if (!(pc->flags & POS_COMPLWHOLE)
        && (!(pc->flags & POS_DYNAMIC) || (pc->flags & POS_COMPLMAX))
        && (pc->type == POS_RES || pc->type == POS_MOL)
        && gmx_ana_index_has_complete_elems(g, pc->itype, top))
    {
        pc->flags &= ~POS_COMPLMAX;
        pc->flags |= POS_COMPLWHOLE;
    }
    /* Setup the gmax field */
    if ((pc->flags & POS_COMPLWHOLE) && !bBase && pc->b.nra > g->isize)
    {
        gmx_ana_index_copy(&pc->gmax, g, true);
    }
    else
    {
        gmx_ana_index_set(&pc->gmax, pc->b.nra, pc->b.a, 0);
    }
}

/*! \brief
 * Checks whether a position calculation should use a base at all.
 *
 * \param[in] pc   Position calculation data to check.
 * \returns   true if \p pc can use a base and gets some benefit out of it,
 *   false otherwise.
 */
static bool
can_use_base(gmx_ana_poscalc_t *pc)
{
    /* For atoms, it should be faster to do a simple copy, so don't use a
     * base. */
    if (pc->type == POS_ATOM)
    {
        return false;
    }
    /* For dynamic selections that do not use completion, it is not possible
     * to use a base. */
    if ((pc->type == POS_RES || pc->type == POS_MOL)
        && (pc->flags & POS_DYNAMIC) && !(pc->flags & (POS_COMPLMAX | POS_COMPLWHOLE)))
    {
        return false;
    }
    /* Dynamic calculations for a single position cannot use a base. */
    if ((pc->type == POS_ALL || pc->type == POS_ALL_PBC)
        && (pc->flags & POS_DYNAMIC))
    {
        return false;
    }
    return true;
}

/*! \brief
 * Checks whether two position calculations should use a common base.
 *
 * \param[in]     pc1 Calculation 1 to check for.
 * \param[in]     pc2 Calculation 2 to check for.
 * \param[in]     g1  Index group structure that contains the atoms from
 *   \p pc1.
 * \param[in,out] g   Working space, should have enough allocated memory to
 *   contain the intersection of the atoms in \p pc1 and \p pc2.
 * \returns   true if the two calculations should be merged to use a common
 *   base, false otherwise.
 */
static bool
should_merge(gmx_ana_poscalc_t *pc1, gmx_ana_poscalc_t *pc2,
             gmx_ana_index_t *g1, gmx_ana_index_t *g)
{
    gmx_ana_index_t  g2;

    /* Do not merge calculations with different mass weighting. */
    if ((pc1->flags & POS_MASS) != (pc2->flags & POS_MASS))
    {
        return false;
    }
    /* Avoid messing up complete calculations. */
    if ((pc1->flags & POS_COMPLWHOLE) != (pc2->flags & POS_COMPLWHOLE))
    {
        return false;
    }
    /* Find the overlap between the calculations. */
    gmx_ana_index_set(&g2, pc2->b.nra, pc2->b.a, 0);
    gmx_ana_index_intersection(g, g1, &g2);
    /* Do not merge if there is no overlap. */
    if (g->isize == 0)
    {
        return false;
    }
    /* Full completion calculations always match if the type is correct. */
    if ((pc1->flags & POS_COMPLWHOLE) && (pc2->flags & POS_COMPLWHOLE)
        && pc1->type == pc2->type)
    {
        return true;
    }
    /* The calculations also match if the intersection consists of full
     * blocks. */
    if (gmx_ana_index_has_full_ablocks(g, &pc1->b)
        && gmx_ana_index_has_full_ablocks(g, &pc2->b))
    {
        return true;
    }
    return false;
}

/*! \brief
 * Creates a static base for position calculation.
 *
 * \param     pc  Data structure to copy.
 * \returns   Pointer to a newly allocated base for \p pc.
 *
 * Creates and returns a deep copy of \p pc, but clears the
 * \ref POS_DYNAMIC and \ref POS_MASKONLY flags.
 * The newly created structure is set as the base (\c gmx_ana_poscalc_t::sbase)
 * of \p pc and inserted in the collection before \p pc.
 */
static gmx_ana_poscalc_t *
create_simple_base(gmx_ana_poscalc_t *pc)
{
    gmx_ana_poscalc_t *base;
    int                flags;

    flags = pc->flags & ~(POS_DYNAMIC | POS_MASKONLY);
    base  = pc->coll->createCalculation(pc->type, flags);
    set_poscalc_maxindex(base, &pc->gmax, true);

    base->p = new gmx_ana_pos_t();

    pc->sbase = base;
    pc->coll->removeCalculation(base);
    pc->coll->insertCalculation(base, pc);

    return base;
}

/*! \brief
 * Merges a calculation into another calculation such that the new calculation
 * can be used as a base.
 *
 * \param[in,out] base Base calculation to merge to.
 * \param[in,out] pc   Position calculation to merge to \p base.
 *
 * After the call, \p base can be used as a base for \p pc (or any calculation
 * that used it as a base).
 * It is assumed that any overlap between \p base and \p pc is in complete
 * blocks, i.e., that the merge is possible.
 */
static void
merge_to_base(gmx_ana_poscalc_t *base, gmx_ana_poscalc_t *pc)
{
    gmx_ana_index_t  gp, gb, g;
    int              isize, bnr;
    int              i, bi, bj, bo;

    base->flags |= pc->flags & (POS_VELOCITIES | POS_FORCES);
    gmx_ana_index_set(&gp, pc->b.nra, pc->b.a, 0);
    gmx_ana_index_set(&gb, base->b.nra, base->b.a, 0);
    isize = gmx_ana_index_difference_size(&gp, &gb);
    if (isize > 0)
    {
        gmx_ana_index_clear(&g);
        gmx_ana_index_reserve(&g, base->b.nra + isize);
        /* Find the new blocks */
        gmx_ana_index_difference(&g, &gp, &gb);
        /* Count the blocks in g */
        i = bi = bnr = 0;
        while (i < g.isize)
        {
            while (pc->b.a[pc->b.index[bi]] != g.index[i])
            {
                ++bi;
            }
            i += pc->b.index[bi+1] - pc->b.index[bi];
            ++bnr;
            ++bi;
        }
        /* Merge the atoms into a temporary structure */
        gmx_ana_index_merge(&g, &gb, &g);
        /* Merge the blocks */
        srenew(base->b.index, base->b.nr + bnr + 1);
        i                   = g.isize - 1;
        bi                  = base->b.nr - 1;
        bj                  = pc->b.nr - 1;
        bo                  = base->b.nr + bnr - 1;
        base->b.index[bo+1] = i + 1;
        while (bo >= 0)
        {
            if (bi < 0 || base->b.a[base->b.index[bi+1]-1] != g.index[i])
            {
                i -= pc->b.index[bj+1] - pc->b.index[bj];
                --bj;
            }
            else
            {
                if (bj >= 0 && pc->b.a[pc->b.index[bj+1]-1] == g.index[i])
                {
                    --bj;
                }
                i -= base->b.index[bi+1] - base->b.index[bi];
                --bi;
            }
            base->b.index[bo] = i + 1;
            --bo;
        }
        base->b.nr           += bnr;
        base->b.nalloc_index += bnr;
        sfree(base->b.a);
        base->b.nra      = g.isize;
        base->b.a        = g.index;
        base->b.nalloc_a = g.isize;
        /* Refresh the gmax field */
        gmx_ana_index_set(&base->gmax, base->b.nra, base->b.a, 0);
    }
}

/*! \brief
 * Merges two bases into one.
 *
 * \param[in,out] tbase Base calculation to merge to.
 * \param[in]     mbase Base calculation to merge to \p tbase.
 *
 * After the call, \p mbase has been freed and \p tbase is used as the base
 * for all calculations that previously had \p mbase as their base.
 * It is assumed that any overlap between \p tbase and \p mbase is in complete
 * blocks, i.e., that the merge is possible.
 */
static void
merge_bases(gmx_ana_poscalc_t *tbase, gmx_ana_poscalc_t *mbase)
{
    gmx_ana_poscalc_t *pc;

    merge_to_base(tbase, mbase);
    mbase->coll->removeCalculation(mbase);
    /* Set tbase as the base for all calculations that had mbase */
    pc = tbase->coll->first_;
    while (pc)
    {
        if (pc->sbase == mbase)
        {
            pc->sbase = tbase;
            tbase->refcount++;
        }
        pc = pc->next;
    }
    /* Free mbase */
    mbase->refcount = 0;
    gmx_ana_poscalc_free(mbase);
}

/*! \brief
 * Setups the static base calculation for a position calculation.
 *
 * \param[in,out] pc   Position calculation to setup the base for.
 */
static void
setup_base(gmx_ana_poscalc_t *pc)
{
    gmx_ana_poscalc_t *base, *pbase, *next;
    gmx_ana_index_t    gp, g;

    /* Exit immediately if pc should not have a base. */
    if (!can_use_base(pc))
    {
        return;
    }

    gmx_ana_index_set(&gp, pc->b.nra, pc->b.a, 0);
    gmx_ana_index_clear(&g);
    gmx_ana_index_reserve(&g, pc->b.nra);
    pbase = pc;
    base  = pc->coll->first_;
    while (base)
    {
        /* Save the next calculation so that we can safely delete base */
        next = base->next;
        /* Skip pc, calculations that already have a base (we should match the
         * base instead), as well as calculations that should not have a base.
         * If the above conditions are met, check whether we should do a
         * merge.
         */
        if (base != pc && !base->sbase && can_use_base(base)
            && should_merge(pbase, base, &gp, &g))
        {
            /* Check whether this is the first base found */
            if (pbase == pc)
            {
                /* Create a real base if one is not present */
                if (!base->p)
                {
                    pbase = create_simple_base(base);
                }
                else
                {
                    pbase = base;
                }
                /* Make it a base for pc as well */
                merge_to_base(pbase, pc);
                pc->sbase = pbase;
                pbase->refcount++;
            }
            else /* This was not the first base */
            {
                if (!base->p)
                {
                    /* If it is not a real base, just make the new base as
                     * the base for it as well. */
                    merge_to_base(pbase, base);
                    base->sbase = pbase;
                    pbase->refcount++;
                }
                else
                {
                    /* If base is a real base, merge it with the new base
                     * and delete it. */
                    merge_bases(pbase, base);
                }
            }
            gmx_ana_index_set(&gp, pbase->b.nra, pbase->b.a, 0);
            gmx_ana_index_reserve(&g, pc->b.nra);
        }
        /* Proceed to the next unchecked calculation */
        base = next;
    }

    gmx_ana_index_deinit(&g);

    /* If no base was found, create one if one is required */
    if (!pc->sbase && (pc->flags & POS_DYNAMIC)
        && (pc->flags & (POS_COMPLMAX | POS_COMPLWHOLE)))
    {
        create_simple_base(pc);
    }
}

/*!
 * \param[in,out] pc    Position calculation data structure.
 * \param[in]     flags New flags.
 *
 * \p flags are added to the old flags.
 * If calculation type is \ref POS_ATOM, \ref POS_MASS is automatically
 * cleared.
 * If both \ref POS_DYNAMIC and \ref POS_MASKONLY are provided,
 * \ref POS_DYNAMIC is cleared.
 * If calculation type is not \ref POS_RES or \ref POS_MOL,
 * \ref POS_COMPLMAX and \ref POS_COMPLWHOLE are automatically cleared.
 */
void
gmx_ana_poscalc_set_flags(gmx_ana_poscalc_t *pc, int flags)
{
    if (pc->type == POS_ATOM)
    {
        flags &= ~POS_MASS;
    }
    if (flags & POS_MASKONLY)
    {
        flags &= ~POS_DYNAMIC;
    }
    if (pc->type != POS_RES && pc->type != POS_MOL)
    {
        flags &= ~(POS_COMPLMAX | POS_COMPLWHOLE);
    }
    pc->flags |= flags;
}

/*!
 * \param[in,out] pc  Position calculation data structure.
 * \param[in]     g   Maximum index group for the calculation.
 *
 * Subsequent calls to gmx_ana_poscalc_update() should use only subsets of
 * \p g for evaluation.
 *
 * The topology should have been set for the collection of which \p pc is
 * a member.
 */
void
gmx_ana_poscalc_set_maxindex(gmx_ana_poscalc_t *pc, gmx_ana_index_t *g)
{
    set_poscalc_maxindex(pc, g, false);
    setup_base(pc);
}

/*!
 * \param[in]  pc  Position calculation data structure.
 * \param[out] p   Output positions.
 *
 * Calls to gmx_ana_poscalc_update() using \p pc should use only positions
 * initialized with this function.
 * The \c p->g pointer is initialized to point to an internal group that
 * contains the maximum index group set with gmx_ana_poscalc_set_maxindex().
 */
void
gmx_ana_poscalc_init_pos(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p)
{
    gmx_ana_indexmap_init(&p->m, &pc->gmax, pc->coll->top_, pc->itype);
    /* Only do the static optimization when there is no completion */
    if (!(pc->flags & POS_DYNAMIC) && pc->b.nra == pc->gmax.isize)
    {
        gmx_ana_indexmap_set_static(&p->m, &pc->b);
    }
    gmx_ana_pos_reserve(p, p->m.mapb.nr, -1);
    if (pc->flags & POS_VELOCITIES)
    {
        gmx_ana_pos_reserve_velocities(p);
    }
    if (pc->flags & POS_FORCES)
    {
        gmx_ana_pos_reserve_forces(p);
    }
}

/*!
 * \param  pc  Position calculation data to be freed.
 *
 * The \p pc pointer is invalid after the call.
 */
void
gmx_ana_poscalc_free(gmx_ana_poscalc_t *pc)
{
    if (!pc)
    {
        return;
    }

    pc->refcount--;
    if (pc->refcount > 0)
    {
        return;
    }

    pc->coll->removeCalculation(pc);
    if (pc->b.nalloc_index > 0)
    {
        sfree(pc->b.index);
    }
    if (pc->b.nalloc_a > 0)
    {
        sfree(pc->b.a);
    }
    if (pc->flags & POS_COMPLWHOLE)
    {
        gmx_ana_index_deinit(&pc->gmax);
    }
    delete pc->p;
    if (pc->sbase)
    {
        gmx_ana_poscalc_free(pc->sbase);
        sfree(pc->baseid);
    }
    sfree(pc);
}

gmx::PositionCalculationCollection::RequiredTopologyInfo
gmx_ana_poscalc_required_topology_info(gmx_ana_poscalc_t *pc)
{
    return gmx::requiredTopologyInfo(pc->type, pc->flags);
}

/*!
 * \param[in]     pc   Position calculation data.
 * \param[in,out] p    Output positions, initialized previously with
 *   gmx_ana_poscalc_init_pos() using \p pc.
 * \param[in]     g    Index group to use for the update.
 * \param[in]     fr   Current frame.
 * \param[in]     pbc  PBC data, or NULL if no PBC should be used.
 *
 * gmx_ana_poscalc_init_frame() should be called for each frame before calling
 * this function.
 */
void
gmx_ana_poscalc_update(gmx_ana_poscalc_t *pc, gmx_ana_pos_t *p,
                       gmx_ana_index_t *g, t_trxframe *fr, const t_pbc *pbc)
{
    int  i, bi, bj;

    if (pc->bEval == true && !(pc->flags & POS_MASKONLY))
    {
        return;
    }
    if (pc->sbase)
    {
        gmx_ana_poscalc_update(pc->sbase, nullptr, nullptr, fr, pbc);
    }
    if (!p)
    {
        p = pc->p;
    }
    if (!g)
    {
        g = &pc->gmax;
    }

    /* Update the index map */
    if (pc->flags & POS_DYNAMIC)
    {
        gmx_ana_indexmap_update(&p->m, g, false);
    }
    else if (pc->flags & POS_MASKONLY)
    {
        gmx_ana_indexmap_update(&p->m, g, true);
        if (pc->bEval)
        {
            return;
        }
    }
    if (!(pc->flags & POS_DYNAMIC))
    {
        pc->bEval = true;
    }

    /* Evaluate the positions */
    if (pc->sbase)
    {
        /* TODO: It might be faster to evaluate the positions within this
         * loop instead of in the beginning. */
        if (pc->flags & POS_DYNAMIC)
        {
            for (bi = 0; bi < p->count(); ++bi)
            {
                bj = pc->baseid[p->m.refid[bi]];
                copy_rvec(pc->sbase->p->x[bj], p->x[bi]);
            }
            if (p->v)
            {
                for (bi = 0; bi < p->count(); ++bi)
                {
                    bj = pc->baseid[p->m.refid[bi]];
                    copy_rvec(pc->sbase->p->v[bj], p->v[bi]);
                }
            }
            if (p->f)
            {
                for (bi = 0; bi < p->count(); ++bi)
                {
                    bj = pc->baseid[p->m.refid[bi]];
                    copy_rvec(pc->sbase->p->f[bj], p->f[bi]);
                }
            }
        }
        else
        {
            for (bi = 0; bi < p->count(); ++bi)
            {
                bj = pc->baseid[bi];
                copy_rvec(pc->sbase->p->x[bj], p->x[bi]);
            }
            if (p->v)
            {
                for (bi = 0; bi < p->count(); ++bi)
                {
                    bj = pc->baseid[bi];
                    copy_rvec(pc->sbase->p->v[bj], p->v[bi]);
                }
            }
            if (p->f)
            {
                for (bi = 0; bi < p->count(); ++bi)
                {
                    bj = pc->baseid[bi];
                    copy_rvec(pc->sbase->p->f[bj], p->f[bi]);
                }
            }
        }
    }
    else /* pc->sbase is NULL */
    {
        if (pc->flags & POS_DYNAMIC)
        {
            pc->b.nr    = p->m.mapb.nr;
            pc->b.index = p->m.mapb.index;
            pc->b.nra   = g->isize;
            pc->b.a     = g->index;
        }
        if (p->v && !fr->bV)
        {
            for (i = 0; i < pc->b.nra; ++i)
            {
                clear_rvec(p->v[i]);
            }
        }
        if (p->f && !fr->bF)
        {
            for (i = 0; i < pc->b.nra; ++i)
            {
                clear_rvec(p->f[i]);
            }
        }
        gmx::ArrayRef<const int> index = pc->coll->getFrameIndices(pc->b.nra, pc->b.a);
        const gmx_mtop_t        *top   = pc->coll->top_;
        const bool               bMass = pc->flags & POS_MASS;
        switch (pc->type)
        {
            case POS_ATOM:
                for (i = 0; i < pc->b.nra; ++i)
                {
                    copy_rvec(fr->x[index[i]], p->x[i]);
                }
                if (p->v && fr->bV)
                {
                    for (i = 0; i < pc->b.nra; ++i)
                    {
                        copy_rvec(fr->v[index[i]], p->v[i]);
                    }
                }
                if (p->f && fr->bF)
                {
                    for (i = 0; i < pc->b.nra; ++i)
                    {
                        copy_rvec(fr->f[index[i]], p->f[i]);
                    }
                }
                break;
            case POS_ALL:
                gmx_calc_comg(top, fr->x, index.size(), index.data(), bMass, p->x[0]);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg(top, fr->v, index.size(), index.data(), bMass, p->v[0]);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_f(top, fr->f, index.size(), index.data(), bMass, p->f[0]);
                }
                break;
            case POS_ALL_PBC:
                gmx_calc_comg_pbc(top, fr->x, pbc, index.size(), index.data(), bMass, p->x[0]);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg(top, fr->v, index.size(), index.data(), bMass, p->v[0]);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_f(top, fr->f, index.size(), index.data(), bMass, p->f[0]);
                }
                break;
            default:
                // TODO: It would probably be better to do this without the type casts.
                gmx_calc_comg_block(top, fr->x, reinterpret_cast<t_block *>(&pc->b),
                                    index.data(), bMass, p->x);
                if (p->v && fr->bV)
                {
                    gmx_calc_comg_block(top, fr->v, reinterpret_cast<t_block *>(&pc->b),
                                        index.data(), bMass, p->v);
                }
                if (p->f && fr->bF)
                {
                    gmx_calc_comg_f_block(top, fr->f, reinterpret_cast<t_block *>(&pc->b),
                                          index.data(), bMass, p->f);
                }
                break;
        }
    }
}
