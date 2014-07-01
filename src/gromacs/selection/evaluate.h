/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * \brief Evaluation functions for sel_evalfunc().
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 * Users should only use SelectionCollection::evaluate() to evaluate
 * selections.
 *
 * The functions defined in this header file are all the possible values
 * for the gmx::SelectionTreeElement::evaluate field (in addition to NULL).
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_EVALUATE_H
#define GMX_SELECTION_EVALUATE_H

#include "selelem.h"

struct gmx_ana_index_t;
struct gmx_sel_mempool_t;
struct t_pbc;
struct t_topology;
struct t_trxframe;

/*! \internal \brief
 * Data structure for passing information required during evaluation.
 */
struct gmx_sel_evaluate_t
{
    /** Memory pool for intermediate values. */
    gmx_sel_mempool_t        *mp;
    /** Index group that contains all the atoms. */
    gmx_ana_index_t          *gall;
    /** Topology information. */
    t_topology               *top;
    /** Current frame. */
    t_trxframe               *fr;
    /** PBC data. */
    t_pbc                    *pbc;
};

/*! \name Utility functions
 */
/*@{*/
/** Initializes an evaluation data structure. */
void
_gmx_sel_evaluate_init(gmx_sel_evaluate_t *data,
                       gmx_sel_mempool_t *mp, gmx_ana_index_t *gall,
                       t_topology *top, t_trxframe *fr, t_pbc *pbc);
/** Evaluates the children of a general selection element. */
void
_gmx_sel_evaluate_children(gmx_sel_evaluate_t                     *data,
                           const gmx::SelectionTreeElementPointer &sel,
                           gmx_ana_index_t                        *g);
/** Evaluates the children of a \ref SEL_EXPRESSION element. */
void
_gmx_sel_evaluate_method_params(gmx_sel_evaluate_t                     *data,
                                const gmx::SelectionTreeElementPointer &sel,
                                gmx_ana_index_t                        *g);
/*@}*/

/*! \name Misc. evaluation functions
 */
/*@{*/
/*! \brief
 * Evaluates a root selection element.
 *
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated
 *   (not used, can be NULL).
 * \returns   0 on success, a non-zero error code on error.
 *
 * Evaluates the first child element in the group defined by \p sel->u.cgrp.
 * If \p sel->u.cgrp is empty, nothing is done.
 * The value of \p sel is not touched (root elements do not evaluate to
 * values).
 *
 * This function can be used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_ROOT elements.
 */
void
_gmx_sel_evaluate_root(gmx_sel_evaluate_t                     *data,
                       const gmx::SelectionTreeElementPointer &sel,
                       gmx_ana_index_t                        *g);
/*! \brief
 * Evaluates a static group selection element.
 *
 * \param[in] data Data for the current frame.
 * \param[in] sel Selection element being evaluated.
 * \param[in] g   Group for which \p sel should be evaluated.
 * \returns   0 for success.
 *
 * Sets the value of \p sel to the intersection of \p g and \p sel->u.cgrp.
 *
 * This function can be used as gmx::SelectionTreeElement::evaluate for
 * \ref SEL_CONST elements with value type \ref GROUP_VALUE.
 */
void
_gmx_sel_evaluate_static(gmx_sel_evaluate_t                     *data,
                         const gmx::SelectionTreeElementPointer &sel,
                         gmx_ana_index_t                        *g);
/** Evaluates an arithmetic expression element. */
void
_gmx_sel_evaluate_arithmetic(gmx_sel_evaluate_t                     *data,
                             const gmx::SelectionTreeElementPointer &sel,
                             gmx_ana_index_t                        *g);
/*@}*/

/*! \name Subexpression evaluation functions
 */
/*@{*/
/** Evaluates a subexpression when there is only one reference. */
void
_gmx_sel_evaluate_subexpr_simple(gmx_sel_evaluate_t                     *data,
                                 const gmx::SelectionTreeElementPointer &sel,
                                 gmx_ana_index_t                        *g);
/** Evaluates a subexpression when the evaluation group is static. */
void
_gmx_sel_evaluate_subexpr_staticeval(gmx_sel_evaluate_t                     *data,
                                     const gmx::SelectionTreeElementPointer &sel,
                                     gmx_ana_index_t                        *g);
/** Evaluates a subexpression. */
void
_gmx_sel_evaluate_subexpr(gmx_sel_evaluate_t                     *data,
                          const gmx::SelectionTreeElementPointer &sel,
                          gmx_ana_index_t                        *g);
/** Evaluates a subexpression reference when there are no other references. */
void
_gmx_sel_evaluate_subexprref_simple(gmx_sel_evaluate_t                     *data,
                                    const gmx::SelectionTreeElementPointer &sel,
                                    gmx_ana_index_t                        *g);
/** Evaluates a subexpression reference. */
void
_gmx_sel_evaluate_subexprref(gmx_sel_evaluate_t                     *data,
                             const gmx::SelectionTreeElementPointer &sel,
                             gmx_ana_index_t                        *g);
/*@}*/

/*! \name Method evaluation functions
 */
/*@{*/

/** Evaluates a method expression. */
void
_gmx_sel_evaluate_method(gmx_sel_evaluate_t                     *data,
                         const gmx::SelectionTreeElementPointer &sel,
                         gmx_ana_index_t                        *g);
/** Evaluates a modifier expression. */
void
_gmx_sel_evaluate_modifier(gmx_sel_evaluate_t                     *data,
                           const gmx::SelectionTreeElementPointer &sel,
                           gmx_ana_index_t                        *g);
/*@}*/

/*! \name Boolean evaluation functions
 */
/*@{*/
/** Evaluates a boolean NOT element. */
void
_gmx_sel_evaluate_not(gmx_sel_evaluate_t                     *data,
                      const gmx::SelectionTreeElementPointer &sel,
                      gmx_ana_index_t                        *g);
/** Evaluates a boolean AND element with short-circuiting. */
void
_gmx_sel_evaluate_and(gmx_sel_evaluate_t                     *data,
                      const gmx::SelectionTreeElementPointer &sel,
                      gmx_ana_index_t                        *g);
/** Evaluates a boolean OR element with short-circuiting. */
void
_gmx_sel_evaluate_or(gmx_sel_evaluate_t                     *data,
                     const gmx::SelectionTreeElementPointer &sel,
                     gmx_ana_index_t                        *g);
/*@}*/

#endif
