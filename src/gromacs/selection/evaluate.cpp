/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements functions in evaluate.h.
 *
 * \todo
 * One of the major bottlenecks for selection performance is that all the
 * evaluation is carried out for atoms.
 * There are several cases when the evaluation could be done for residues
 * or molecules instead, including keywords that select by residue and
 * cases where residue centers are used as reference positions.
 * Implementing this would require a mechanism for recognizing whether
 * something can be evaluated by residue/molecule instead by atom, and
 * converting selections by residue/molecule into selections by atom
 * when necessary.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "evaluate.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <memory>

#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "mempool.h"
#include "poscalc.h"
#include "selectioncollection_impl.h"
#include "selelem.h"
#include "selmethod.h"

using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;

namespace
{

/*! \brief
 * Reserves memory for a selection element from the evaluation memory pool.
 *
 * This class implements RAII semantics for allocating memory for selection
 * element values from a selection evaluation memory pool.
 *
 * \ingroup module_selection
 */
class MempoolSelelemReserver
{
public:
    //! Constructs a reserver without initial reservation.
    MempoolSelelemReserver() {}
    /*! \brief
     * Constructs a reserver with initial reservation.
     *
     * \param[in,out] sel    Selection element for which to reserve.
     * \param[in]     count  Number of values to reserve.
     *
     * \see reserve()
     */
    MempoolSelelemReserver(const SelectionTreeElementPointer& sel, int count)
    {
        reserve(sel, count);
    }
    //! Frees any memory allocated using this reserver.
    ~MempoolSelelemReserver()
    {
        if (sel_)
        {
            sel_->mempoolRelease();
        }
    }

    /*! \brief
     * Reserves memory for selection element values using this reserver.
     *
     * \param[in,out] sel    Selection element for which to reserve.
     * \param[in]     count  Number of values to reserve.
     *
     * Allocates space to store \p count output values in \p sel from the
     * memory pool associated with \p sel, or from the heap if there is no
     * memory pool.  Type of values to allocate is automatically determined
     * from \p sel.
     */
    void reserve(const SelectionTreeElementPointer& sel, int count)
    {
        GMX_RELEASE_ASSERT(!sel_, "Can only reserve one element with one instance");
        sel->mempoolReserve(count);
        sel_ = sel;
    }

private:
    SelectionTreeElementPointer sel_;
};

/*! \brief
 * Reserves memory for an index group from the evaluation memory pool.
 *
 * This class implements RAII semantics for allocating memory for an index
 * group from a selection evaluation memory pool.
 *
 * \ingroup module_selection
 */
class MempoolGroupReserver
{
public:
    /*! \brief
     * Creates a reserver associated with a given memory pool.
     *
     * \param    mp  Memory pool from which to reserve memory.
     */
    explicit MempoolGroupReserver(gmx_sel_mempool_t* mp) : mp_(mp), g_(nullptr) {}
    //! Frees any memory allocated using this reserver.
    ~MempoolGroupReserver()
    {
        if (g_ != nullptr)
        {
            _gmx_sel_mempool_free_group(mp_, g_);
        }
    }

    /*! \brief
     * Reserves memory for an index group using this reserver.
     *
     * \param[in,out] g      Index group to reserve.
     * \param[in]     count  Number of atoms to reserve space for.
     *
     * Allocates memory from the memory pool to store \p count atoms in
     * \p g.
     */
    void reserve(gmx_ana_index_t* g, int count)
    {
        GMX_RELEASE_ASSERT(g_ == nullptr, "Can only reserve one element with one instance");
        _gmx_sel_mempool_alloc_group(mp_, g, count);
        g_ = g;
    }

private:
    gmx_sel_mempool_t* mp_;
    gmx_ana_index_t*   g_;
};

/*! \brief
 * Assigns a temporary value for a selection element.
 *
 * This class implements RAII semantics for temporarily assigning the value
 * pointer of a selection element to point to a different location.
 *
 * \ingroup module_selection
 */
class SelelemTemporaryValueAssigner
{
public:
    //! Constructs an assigner without an initial assignment.
    SelelemTemporaryValueAssigner() : old_ptr_(nullptr), old_nalloc_(0) {}
    /*! \brief
     * Constructs an assigner with an initial assignment.
     *
     * \param[in,out] sel     Selection element for which to assign.
     * \param[in]     vsource Element to which \p sel values will point to.
     *
     * \see assign()
     */
    SelelemTemporaryValueAssigner(const SelectionTreeElementPointer& sel, const SelectionTreeElement& vsource)
    {
        assign(sel, vsource);
    }
    //! Undoes any temporary assignment done using this assigner.
    ~SelelemTemporaryValueAssigner()
    {
        if (sel_)
        {
            _gmx_selvalue_setstore_alloc(&sel_->v, old_ptr_, old_nalloc_);
        }
    }

    /*! \brief
     * Assigns a temporary value pointer.
     *
     * \param[in,out] sel     Selection element for which to assign.
     * \param[in]     vsource Element to which \p sel values will point to.
     *
     * Assigns the value pointer in \p sel to point to the values in
     * \p vsource, i.e., any access/modification to values in \p sel
     * actually accesses values in \p vsource.
     */
    void assign(const SelectionTreeElementPointer& sel, const SelectionTreeElement& vsource)
    {
        GMX_RELEASE_ASSERT(!sel_, "Can only assign one element with one instance");
        GMX_RELEASE_ASSERT(sel->v.type == vsource.v.type, "Mismatching selection value types");
        _gmx_selvalue_getstore_and_release(&sel->v, &old_ptr_, &old_nalloc_);
        _gmx_selvalue_setstore(&sel->v, vsource.v.u.ptr);
        sel_ = sel;
    }

private:
    SelectionTreeElementPointer sel_;
    void*                       old_ptr_;
    int                         old_nalloc_;
};

/*! \brief
 * Expands a value array from one-per-position to one-per-atom.
 *
 * \param[in,out] value  Array to expand.
 * \param[in,out] nr     Number of values in \p value.
 * \param[in]     pos    Position data.
 * \tparam        T      Value type of the array to expand.
 *
 * On input, \p value contains one value for each position in \p pos (and `*nr`
 * must match).  On output, \p value will contain a value for each atom used to
 * evaluate `pos`: each input value is replicated to all atoms that make up the
 * corresponding position.
 * The operation is done in-place.
 *
 * Does not throw.
 */
template<typename T>
void expandValueForPositions(T value[], int* nr, gmx_ana_pos_t* pos)
{
    GMX_RELEASE_ASSERT(*nr == pos->count(),
                       "Position update method did not return the correct number of values");
    *nr = pos->m.mapb.nra;
    // Loop over the positions in reverse order so that the expansion can be
    // done in-place (assumes that each position has at least one atom, which
    // should always be the case).
    int outputIndex = pos->m.mapb.nra;
    for (int i = pos->count() - 1; i >= 0; --i)
    {
        const int atomCount = pos->m.mapb.index[i + 1] - pos->m.mapb.index[i];
        outputIndex -= atomCount;
        GMX_ASSERT(outputIndex >= i, "In-place algorithm would overwrite data not yet used");
        std::fill(&value[outputIndex], &value[outputIndex + atomCount], value[i]);
    }
}

} // namespace

/*!
 * \param[in] fp       File handle to receive the output.
 * \param[in] evalfunc Function pointer to print.
 */
void _gmx_sel_print_evalfunc_name(FILE* fp, gmx::sel_evalfunc evalfunc)
{
    if (!evalfunc)
    {
        fprintf(fp, "none");
    }
    else if (evalfunc == &_gmx_sel_evaluate_root)
    {
        fprintf(fp, "root");
    }
    else if (evalfunc == &_gmx_sel_evaluate_static)
    {
        fprintf(fp, "static");
    }
    else if (evalfunc == &_gmx_sel_evaluate_subexpr_simple)
    {
        fprintf(fp, "subexpr_simple");
    }
    else if (evalfunc == &_gmx_sel_evaluate_subexpr_staticeval)
    {
        fprintf(fp, "subexpr_staticeval");
    }
    else if (evalfunc == &_gmx_sel_evaluate_subexpr)
    {
        fprintf(fp, "subexpr");
    }
    else if (evalfunc == &_gmx_sel_evaluate_subexprref_simple)
    {
        fprintf(fp, "ref_simple");
    }
    else if (evalfunc == &_gmx_sel_evaluate_subexprref)
    {
        fprintf(fp, "ref");
    }
    else if (evalfunc == &_gmx_sel_evaluate_method)
    {
        fprintf(fp, "method");
    }
    else if (evalfunc == &_gmx_sel_evaluate_modifier)
    {
        fprintf(fp, "mod");
    }
    else if (evalfunc == &_gmx_sel_evaluate_not)
    {
        fprintf(fp, "not");
    }
    else if (evalfunc == &_gmx_sel_evaluate_and)
    {
        fprintf(fp, "and");
    }
    else if (evalfunc == &_gmx_sel_evaluate_or)
    {
        fprintf(fp, "or");
    }
    else if (evalfunc == &_gmx_sel_evaluate_arithmetic)
    {
        fprintf(fp, "arithmetic");
    }
    else
    {
        fprintf(fp, "%p", reinterpret_cast<void*>(evalfunc));
    }
}

/*!
 * \param[out] data Evaluation data structure to initialize.
 * \param[in]  mp   Memory pool for intermediate evaluation values.
 * \param[in]  gall Index group with all the atoms.
 * \param[in]  top  Topology structure for evaluation.
 * \param[in]  fr   New frame for evaluation.
 * \param[in]  pbc  New PBC information for evaluation.
 */
void _gmx_sel_evaluate_init(gmx_sel_evaluate_t* data,
                            gmx_sel_mempool_t*  mp,
                            gmx_ana_index_t*    gall,
                            const gmx_mtop_t*   top,
                            t_trxframe*         fr,
                            t_pbc*              pbc)
{
    data->mp   = mp;
    data->gall = gall;
    data->top  = top;
    data->fr   = fr;
    data->pbc  = pbc;
}

/*! \brief
 * Recursively initializes the flags for evaluation.
 *
 * \param[in,out] sel Selection element to clear.
 *
 * The \ref SEL_INITFRAME flag is set for \ref SEL_EXPRESSION elements whose
 * method defines the \p init_frame callback (see sel_framefunc()), and
 * cleared for other elements.
 *
 * The \ref SEL_EVALFRAME flag is cleared for all elements.
 */
static void init_frame_eval(SelectionTreeElementPointer sel)
{
    while (sel)
    {
        sel->flags &= ~(SEL_INITFRAME | SEL_EVALFRAME);
        if (sel->type == SEL_EXPRESSION)
        {
            if (sel->u.expr.method && sel->u.expr.method->init_frame)
            {
                sel->flags |= SEL_INITFRAME;
            }
        }
        if (sel->child && sel->type != SEL_SUBEXPRREF)
        {
            init_frame_eval(sel->child);
        }
        sel = sel->next;
    }
}

namespace gmx
{

SelectionEvaluator::SelectionEvaluator() {}

/*!
 * \param[in,out] coll  The selection collection to evaluate.
 * \param[in] fr  Frame for which the evaluation should be carried out.
 * \param[in] pbc PBC data, or NULL if no PBC should be used.
 * \returns   0 on successful evaluation, a non-zero error code on error.
 *
 * This functions sets the global variables for topology, frame and PBC,
 * clears some information in the selection to initialize the evaluation
 * for a new frame, and evaluates \p sel and all the selections pointed by
 * the \p next pointers of \p sel.
 *
 * This is the only function that user code should call if they want to
 * evaluate a selection for a new frame.
 */
// NOLINTNEXTLINE readability-convert-member-functions-to-static
void SelectionEvaluator::evaluate(SelectionCollection* coll, t_trxframe* fr, t_pbc* pbc)
{
    gmx_ana_selcollection_t* sc = &coll->impl_->sc_;
    gmx_sel_evaluate_t       data;

    _gmx_sel_evaluate_init(&data, sc->mempool, &sc->gall, sc->top, fr, pbc);
    init_frame_eval(sc->root);
    SelectionTreeElementPointer sel = sc->root;
    while (sel)
    {
        /* Clear the evaluation group of subexpressions */
        if (sel->child && sel->child->type == SEL_SUBEXPR && sel->child->evaluate != nullptr)
        {
            sel->child->u.cgrp.isize = 0;
            /* Not strictly necessary, because the value will be overwritten
             * during first evaluation of the subexpression anyways, but we
             * clear the group for clarity. Note that this is _not_ done during
             * compilation because of some additional complexities involved
             * (see compiler.cpp), so it should not be relied upon in
             * _gmx_sel_evaluate_subexpr(). */
            if (sel->child->v.type == GROUP_VALUE)
            {
                sel->child->v.u.g->isize = 0;
            }
        }
        if (sel->evaluate)
        {
            sel->evaluate(&data, sel, nullptr);
        }
        sel = sel->next;
    }
    /* Update selection information */
    SelectionDataList::const_iterator isel;
    for (isel = sc->sel.begin(); isel != sc->sel.end(); ++isel)
    {
        internal::SelectionData& sel = **isel;
        sel.refreshMassesAndCharges(sc->top);
        sel.updateCoveredFractionForFrame();
    }
}

/*!
 * \param[in,out] coll  The selection collection to evaluate.
 * \param[in]     nframes Total number of frames.
 */
// NOLINTNEXTLINE readability-convert-member-functions-to-static
void SelectionEvaluator::evaluateFinal(SelectionCollection* coll, int nframes)
{
    gmx_ana_selcollection_t* sc = &coll->impl_->sc_;

    SelectionDataList::const_iterator isel;
    for (isel = sc->sel.begin(); isel != sc->sel.end(); ++isel)
    {
        internal::SelectionData& sel = **isel;
        sel.restoreOriginalPositions(sc->top);
        sel.computeAverageCoveredFraction(nframes);
    }
}

} // namespace gmx

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel  Selection element being evaluated.
 * \param[in] g    Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates each child of \p sel in \p g.
 */
void _gmx_sel_evaluate_children(gmx_sel_evaluate_t*                     data,
                                const gmx::SelectionTreeElementPointer& sel,
                                gmx_ana_index_t*                        g)
{
    SelectionTreeElementPointer child = sel->child;
    while (child)
    {
        if (child->evaluate)
        {
            child->evaluate(data, child, g);
        }
        child = child->next;
    }
}

void _gmx_sel_evaluate_root(gmx_sel_evaluate_t*                     data,
                            const gmx::SelectionTreeElementPointer& sel,
                            gmx_ana_index_t* /* g */)
{
    if (sel->u.cgrp.isize == 0 || !sel->child->evaluate)
    {
        return;
    }

    sel->child->evaluate(data, sel->child, sel->u.cgrp.isize < 0 ? nullptr : &sel->u.cgrp);
}

void _gmx_sel_evaluate_static(gmx_sel_evaluate_t* /* data */,
                              const gmx::SelectionTreeElementPointer& sel,
                              gmx_ana_index_t*                        g)
{
    if (sel->flags & SEL_UNSORTED)
    {
        gmx_ana_index_reserve(sel->v.u.g, sel->u.cgrp.isize);
        // This only works if g contains all the atoms, but that is currently
        // the only supported case.
        gmx_ana_index_copy(sel->v.u.g, &sel->u.cgrp, false);
    }
    else
    {
        gmx_ana_index_intersection(sel->v.u.g, &sel->u.cgrp, g);
    }
}


/*********************************************************************
 * SUBEXPRESSION EVALUATION
 *********************************************************************/

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel  Selection element being evaluated.
 * \param[in] g    Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates the child element (there should be exactly one) in \p g.
 * The compiler has taken care that the child actually stores the evaluated
 * value in the value pointer of this element.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_SUBEXPR elements that are used only once, and hence do not need
 * full subexpression handling.
 */
void _gmx_sel_evaluate_subexpr_simple(gmx_sel_evaluate_t*                     data,
                                      const gmx::SelectionTreeElementPointer& sel,
                                      gmx_ana_index_t*                        g)
{
    if (sel->child->evaluate)
    {
        sel->child->evaluate(data, sel->child, g);
    }
    sel->v.nr = sel->child->v.nr;
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel  Selection element being evaluated.
 * \param[in] g    Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * If this is the first call for this frame, evaluates the child element
 * there should be exactly one in \p g.
 * The compiler has taken care that the child actually stores the evaluated
 * value in the value pointer of this element.
 * Assumes that \p g is persistent for the duration of the whole evaluation.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_SUBEXPR elements that have a static evaluation group, and hence do
 * not need full subexpression handling.
 */
void _gmx_sel_evaluate_subexpr_staticeval(gmx_sel_evaluate_t*                     data,
                                          const gmx::SelectionTreeElementPointer& sel,
                                          gmx_ana_index_t*                        g)
{
    if (sel->u.cgrp.isize == 0)
    {
        sel->child->evaluate(data, sel->child, g);
        sel->v.nr = sel->child->v.nr;
        if (!g)
        {
            sel->u.cgrp.isize = -1;
        }
        else
        {
            gmx_ana_index_set(&sel->u.cgrp, g->isize, g->index, 0);
        }
    }
}

/*!
 * \param[in]  data  Data for the current frame.
 * \param[in]  sel   Selection element being evaluated.
 * \param[in]  g     Group for which \p sel should be evaluated.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Finds the part of \p g for which the subexpression
 * has not yet been evaluated by comparing \p g to \p sel->u.cgrp.
 * If the part is not empty, the child expression is evaluated for this
 * part, and the results merged to the old values of the child.
 * The value of \p sel itself is undefined after the call.
 *
 * \todo
 * The call to gmx_ana_index_difference() can take quite a lot of unnecessary
 * time if the subexpression is evaluated either several times for the same
 * group or for completely distinct groups.
 * However, in the majority of cases, these situations occur when
 * _gmx_sel_evaluate_subexpr_staticeval() can be used, so this should not be a
 * major problem.
 */
void _gmx_sel_evaluate_subexpr(gmx_sel_evaluate_t*                     data,
                               const gmx::SelectionTreeElementPointer& sel,
                               gmx_ana_index_t*                        g)
{
    gmx_ana_index_t gmiss;

    MempoolGroupReserver gmissreserver(data->mp);
    if (sel->u.cgrp.isize == 0)
    {
        {
            SelelemTemporaryValueAssigner assigner(sel->child, *sel);
            sel->child->evaluate(data, sel->child, g);
        }
        gmx_ana_index_copy(&sel->u.cgrp, g, false);
        gmiss.isize = 0;
    }
    else
    {
        gmissreserver.reserve(&gmiss, g->isize);
        gmx_ana_index_difference(&gmiss, g, &sel->u.cgrp);
    }
    if (gmiss.isize > 0)
    {
        MempoolSelelemReserver reserver(sel->child, gmiss.isize);
        /* Evaluate the missing values for the child */
        sel->child->evaluate(data, sel->child, &gmiss);
        /* Merge the missing values to the existing ones. */
        if (sel->v.type == GROUP_VALUE)
        {
            gmx_ana_index_merge(sel->v.u.g, sel->child->v.u.g, sel->v.u.g);
        }
        else
        {
            int i, j, k;

            i = sel->u.cgrp.isize - 1;
            j = gmiss.isize - 1;
            /* TODO: This switch is kind of ugly, but it may be difficult to
             * do this portably without C++ templates. */
            switch (sel->v.type)
            {
                case INT_VALUE:
                    for (k = sel->u.cgrp.isize + gmiss.isize - 1; k >= 0; k--)
                    {
                        if (i < 0 || (j >= 0 && sel->u.cgrp.index[i] < gmiss.index[j]))
                        {
                            sel->v.u.i[k] = sel->child->v.u.i[j--];
                        }
                        else
                        {
                            sel->v.u.i[k] = sel->v.u.i[i--];
                        }
                    }
                    break;

                case REAL_VALUE:
                    for (k = sel->u.cgrp.isize + gmiss.isize - 1; k >= 0; k--)
                    {
                        if (i < 0 || (j >= 0 && sel->u.cgrp.index[i] < gmiss.index[j]))
                        {
                            sel->v.u.r[k] = sel->child->v.u.r[j--];
                        }
                        else
                        {
                            sel->v.u.r[k] = sel->v.u.r[i--];
                        }
                    }
                    break;

                case STR_VALUE:
                    // Note: with the currently allowed syntax, this case is never
                    // reached.
                    for (k = sel->u.cgrp.isize + gmiss.isize - 1; k >= 0; k--)
                    {
                        if (i < 0 || (j >= 0 && sel->u.cgrp.index[i] < gmiss.index[j]))
                        {
                            sel->v.u.s[k] = sel->child->v.u.s[j--];
                        }
                        else
                        {
                            sel->v.u.s[k] = sel->v.u.s[i--];
                        }
                    }
                    break;

                case POS_VALUE:
                    /* TODO: Implement this */
                    GMX_THROW(gmx::NotImplementedError(
                            "position subexpressions not implemented properly"));

                case NO_VALUE:
                case GROUP_VALUE: GMX_THROW(gmx::InternalError("Invalid subexpression type"));
            }
        }
        gmx_ana_index_merge(&sel->u.cgrp, &sel->u.cgrp, &gmiss);
    }
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 for success.
 *
 * Sets the value pointers of the child and its child to point to the same
 * memory as the value pointer of this element to avoid copying, and then
 * evaluates evaluates the child.
 *
 * This function is used as gmx::SelectionTreeElement:evaluate for
 * \ref SEL_SUBEXPRREF elements for which the \ref SEL_SUBEXPR does not have
 * other references.
 */
void _gmx_sel_evaluate_subexprref_simple(gmx_sel_evaluate_t*                     data,
                                         const gmx::SelectionTreeElementPointer& sel,
                                         gmx_ana_index_t*                        g)
{
    if (g)
    {
        _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
        _gmx_selvalue_setstore_alloc(&sel->child->child->v, sel->v.u.ptr, sel->child->child->v.nalloc);
        sel->child->evaluate(data, sel->child, g);
    }
    sel->v.nr = sel->child->v.nr;
    if (sel->u.param)
    {
        sel->u.param->val.nr = sel->v.nr;
        if (sel->u.param->nvalptr)
        {
            *sel->u.param->nvalptr = sel->u.param->val.nr;
        }
    }
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * If the value type is \ref POS_VALUE, the value of the child is simply
 * copied to set the value of \p sel (the child subexpression should
 * already have been evaluated by its root).
 * If the value type is something else, the child is evaluated for the
 * group \p g, and the value of the child is then copied.
 * There should be only one child element.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_SUBEXPRREF elements.
 */
void _gmx_sel_evaluate_subexprref(gmx_sel_evaluate_t*                     data,
                                  const gmx::SelectionTreeElementPointer& sel,
                                  gmx_ana_index_t*                        g)
{
    int i, j;

    if (g != nullptr && sel->child->evaluate != nullptr)
    {
        sel->child->evaluate(data, sel->child, g);
    }
    const SelectionTreeElementPointer& expr = sel->child;
    switch (sel->v.type)
    {
        case INT_VALUE:
            if (!g)
            {
                sel->v.nr = expr->v.nr;
                memcpy(sel->v.u.i, expr->v.u.i, sel->v.nr * sizeof(*sel->v.u.i));
            }
            else
            {
                sel->v.nr = g->isize;
                /* Extract the values corresponding to g */
                for (i = j = 0; i < g->isize; ++i, ++j)
                {
                    while (sel->child->u.cgrp.index[j] < g->index[i])
                    {
                        ++j;
                    }
                    sel->v.u.i[i] = expr->v.u.i[j];
                }
            }
            break;

        case REAL_VALUE:
            if (!g)
            {
                sel->v.nr = expr->v.nr;
                memcpy(sel->v.u.r, expr->v.u.r, sel->v.nr * sizeof(*sel->v.u.r));
            }
            else
            {
                sel->v.nr = g->isize;
                /* Extract the values corresponding to g */
                for (i = j = 0; i < g->isize; ++i, ++j)
                {
                    while (sel->child->u.cgrp.index[j] < g->index[i])
                    {
                        ++j;
                    }
                    sel->v.u.r[i] = expr->v.u.r[j];
                }
            }
            break;

        case STR_VALUE:
            if (!g)
            {
                sel->v.nr = expr->v.nr;
                memcpy(sel->v.u.s, expr->v.u.s, sel->v.nr * sizeof(*sel->v.u.s));
            }
            else
            {
                sel->v.nr = g->isize;
                /* Extract the values corresponding to g */
                for (i = j = 0; i < g->isize; ++i, ++j)
                {
                    while (sel->child->u.cgrp.index[j] < g->index[i])
                    {
                        ++j;
                    }
                    sel->v.u.s[i] = expr->v.u.s[j];
                }
            }
            break;

        case POS_VALUE:
            /* Currently, there is no need to do anything fancy here,
             * but some future extensions may need a more flexible
             * implementation. */
            gmx_ana_pos_copy(sel->v.u.p, expr->v.u.p, false);
            break;

        case GROUP_VALUE:
            if (!g)
            {
                gmx_ana_index_copy(sel->v.u.g, expr->v.u.g, false);
            }
            else
            {
                gmx_ana_index_intersection(sel->v.u.g, expr->v.u.g, g);
            }
            break;

        default: /* should not be reached */
            GMX_THROW(gmx::InternalError("Invalid subexpression reference type"));
    }
    /* Store the number of values if needed */
    if (sel->u.param)
    {
        sel->u.param->val.nr = sel->v.nr;
        if (sel->u.param->nvalptr)
        {
            *sel->u.param->nvalptr = sel->u.param->val.nr;
        }
    }
}

/********************************************************************
 * METHOD EXPRESSION EVALUATION
 ********************************************************************/

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates each child of a \ref SEL_EXPRESSION element.
 * The value of \p sel is not touched.
 *
 * This function is not used as gmx::SelectionTreeElement::evaluate,
 * but is used internally.
 */
void _gmx_sel_evaluate_method_params(gmx_sel_evaluate_t*                     data,
                                     const gmx::SelectionTreeElementPointer& sel,
                                     gmx_ana_index_t*                        g)
{
    SelectionTreeElementPointer child = sel->child;
    while (child)
    {
        if (child->evaluate && !(child->flags & SEL_EVALFRAME))
        {
            if (child->flags & SEL_ATOMVAL)
            {
                child->evaluate(data, child, g);
            }
            else
            {
                child->flags |= SEL_EVALFRAME;
                child->evaluate(data, child, nullptr);
            }
        }
        child = child->next;
    }
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates all child selections (using _gmx_sel_evaluate_method_params())
 * to evaluate any parameter values.
 * If this is the first time this expression is evaluated for
 * the frame, sel_framefunc() callback is called if one is provided.
 * If a reference position calculation has been initialized for this element,
 * the positions are also updated, and sel_updatefunc_pos() is used to
 * evaluate the value. Otherwise, sel_updatefunc() is used.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_EXPRESSION elements.
 */
void _gmx_sel_evaluate_method(gmx_sel_evaluate_t*                     data,
                              const gmx::SelectionTreeElementPointer& sel,
                              gmx_ana_index_t*                        g)
{
    _gmx_sel_evaluate_method_params(data, sel, g);
    gmx::SelMethodEvalContext context(data->top, data->fr, data->pbc);
    if (sel->flags & SEL_INITFRAME)
    {
        sel->flags &= ~SEL_INITFRAME;
        sel->u.expr.method->init_frame(context, sel->u.expr.mdata);
    }
    if (sel->u.expr.pc)
    {
        gmx_ana_poscalc_update(sel->u.expr.pc, sel->u.expr.pos, g, data->fr, data->pbc);
        sel->u.expr.method->pupdate(context, sel->u.expr.pos, &sel->v, sel->u.expr.mdata);
        if ((sel->flags & SEL_ATOMVAL) && sel->v.nr < g->isize)
        {
            switch (sel->v.type)
            {
                case REAL_VALUE:
                    expandValueForPositions(sel->v.u.r, &sel->v.nr, sel->u.expr.pos);
                    break;
                default:
                    GMX_RELEASE_ASSERT(false, "Unimplemented value type for position update method");
            }
        }
    }
    else
    {
        sel->u.expr.method->update(context, g, &sel->v, sel->u.expr.mdata);
    }
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates all child selections (using _gmx_sel_evaluate_method_params())
 * to evaluate any parameter values.
 * If this is the first time this expression is evaluated for
 * the frame, sel_framefunc() callback is called if one is provided.
 * The modifier is then evaluated using sel_updatefunc_pos().
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_MODIFIER elements.
 */
void _gmx_sel_evaluate_modifier(gmx_sel_evaluate_t*                     data,
                                const gmx::SelectionTreeElementPointer& sel,
                                gmx_ana_index_t*                        g)
{
    _gmx_sel_evaluate_method_params(data, sel, g);
    gmx::SelMethodEvalContext context(data->top, data->fr, data->pbc);
    if (sel->flags & SEL_INITFRAME)
    {
        sel->flags &= ~SEL_INITFRAME;
        sel->u.expr.method->init_frame(context, sel->u.expr.mdata);
    }
    if (sel->child && sel->child->v.type != POS_VALUE)
    {
        GMX_THROW(gmx::NotImplementedError("Non-position valued modifiers not implemented"));
    }
    sel->u.expr.method->pupdate(context, nullptr, &sel->v, sel->u.expr.mdata);
}


/********************************************************************
 * BOOLEAN EXPRESSION EVALUATION
 ********************************************************************/

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates the child element (there should be only one) in the group
 * \p g, and then sets the value of \p sel to the complement of the
 * child value.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_BOOLEAN elements with \ref BOOL_NOT.
 */
void _gmx_sel_evaluate_not(gmx_sel_evaluate_t*                     data,
                           const gmx::SelectionTreeElementPointer& sel,
                           gmx_ana_index_t*                        g)
{
    MempoolSelelemReserver reserver(sel->child, g->isize);
    sel->child->evaluate(data, sel->child, g);
    gmx_ana_index_difference(sel->v.u.g, g, sel->child->v.u.g);
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Short-circuiting evaluation of logical AND expressions.
 *
 * Starts by evaluating the first child element in the group \p g.
 * The each following child element is evaluated in the intersection
 * of all the previous values until all children have been evaluated
 * or the intersection becomes empty.
 * The value of \p sel is set to the intersection of all the (evaluated)
 * child values.
 *
 * If the first child does not have an evaluation function, it is skipped
 * and the evaluation is started at the second child.
 * This happens if the first child is a constant expression and during
 * compilation it was detected that the evaluation group is always a subset
 * of the constant group
 * (currently, the compiler never detects this).
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_BOOLEAN elements with \ref BOOL_AND.
 */
void _gmx_sel_evaluate_and(gmx_sel_evaluate_t*                     data,
                           const gmx::SelectionTreeElementPointer& sel,
                           gmx_ana_index_t*                        g)
{
    SelectionTreeElementPointer child = sel->child;
    /* Skip the first child if it does not have an evaluation function. */
    if (!child->evaluate)
    {
        child = child->next;
    }
    /* Evaluate the first child */
    {
        MempoolSelelemReserver reserver(child, g->isize);
        child->evaluate(data, child, g);
        gmx_ana_index_copy(sel->v.u.g, child->v.u.g, false);
    }
    child = child->next;
    while (child && sel->v.u.g->isize > 0)
    {
        MempoolSelelemReserver reserver(child, sel->v.u.g->isize);
        child->evaluate(data, child, sel->v.u.g);
        gmx_ana_index_intersection(sel->v.u.g, sel->v.u.g, child->v.u.g);
        child = child->next;
    }
}

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Short-circuiting evaluation of logical OR expressions.
 *
 * Starts by evaluating the first child element in the group \p g.
 * For each subsequent child, finds the part of \p g that is not
 * included the value of any previous child, and evaluates the child
 * in that group until the last child is evaluated or all of \p g
 * is included in some child value.
 * The value of \p sel is set to the union of all the (evaluated)
 * child values.
 *
 * If the first child does not have an evaluation function, its value is
 * used without evaluation.
 * This happens if the first child is a constant expression, the selection
 * has been compiled, and the evaluation group is the same for each frame.
 * In this case, the compiler has taken care of that the child value is a
 * subset of \p g, making it unnecessary to evaluate it.
 *
 * This function is used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_BOOLEAN elements with \ref BOOL_OR.
 */
void _gmx_sel_evaluate_or(gmx_sel_evaluate_t*                     data,
                          const gmx::SelectionTreeElementPointer& sel,
                          gmx_ana_index_t*                        g)
{
    gmx_ana_index_t tmp, tmp2;

    SelectionTreeElementPointer child = sel->child;
    if (child->evaluate)
    {
        MempoolSelelemReserver reserver(child, g->isize);
        child->evaluate(data, child, g);
        gmx_ana_index_partition(sel->v.u.g, &tmp, g, child->v.u.g);
    }
    else
    {
        gmx_ana_index_partition(sel->v.u.g, &tmp, g, child->v.u.g);
    }
    child = child->next;
    while (child && tmp.isize > 0)
    {
        {
            MempoolSelelemReserver reserver(child, tmp.isize);
            child->evaluate(data, child, &tmp);
            gmx_ana_index_partition(&tmp, &tmp2, &tmp, child->v.u.g);
        }
        sel->v.u.g->isize += tmp.isize;
        tmp.isize = tmp2.isize;
        tmp.index = tmp2.index;
        child     = child->next;
    }
    gmx_ana_index_sort(sel->v.u.g);
}


/********************************************************************
 * ARITHMETIC EVALUATION
 ********************************************************************/

/*!
 * \param[in] data Data for the current frame.
 * \param[in] sel  Selection element being evaluated.
 * \param[in] g    Group for which \p sel should be evaluated.
 * \returns   0 on success, a non-zero error code on error.
 */
void _gmx_sel_evaluate_arithmetic(gmx_sel_evaluate_t*                     data,
                                  const gmx::SelectionTreeElementPointer& sel,
                                  gmx_ana_index_t*                        g)
{
    int  n, i, i1, i2;
    real lval, rval = 0., val = 0.;

    const SelectionTreeElementPointer& left  = sel->child;
    const SelectionTreeElementPointer& right = left->next;

    SelelemTemporaryValueAssigner assigner;
    MempoolSelelemReserver        reserver;
    if (left->mempool)
    {
        assigner.assign(left, *sel);
        if (right)
        {
            reserver.reserve(right, g->isize);
        }
    }
    else if (right && right->mempool)
    {
        assigner.assign(right, *sel);
    }
    _gmx_sel_evaluate_children(data, sel, g);

    n         = (sel->flags & SEL_SINGLEVAL) ? 1 : g->isize;
    sel->v.nr = n;

    bool bArithNeg = (sel->u.type == ARITH_NEG);
    GMX_ASSERT(right || bArithNeg, "Right operand cannot be null except for negations");
    for (i = i1 = i2 = 0; i < n; ++i)
    {
        lval = left->v.u.r[i1];
        if (!bArithNeg)
        {
            rval = right->v.u.r[i2];
        }
        switch (sel->u.type)
        {
            case ARITH_PLUS: val = lval + rval; break;
            case ARITH_MINUS: val = lval - rval; break;
            case ARITH_NEG: val = -lval; break;
            case ARITH_MULT: val = lval * rval; break;
            case ARITH_DIV: val = lval / rval; break;
            case ARITH_EXP: val = std::pow(lval, rval); break;
        }
        sel->v.u.r[i] = val;
        if (!(left->flags & SEL_SINGLEVAL))
        {
            ++i1;
        }
        if (!bArithNeg && !(right->flags & SEL_SINGLEVAL))
        {
            ++i2;
        }
    }
}
