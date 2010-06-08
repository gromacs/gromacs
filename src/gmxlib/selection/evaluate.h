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
 * \brief Evaluation functions for sel_evalfunc().
 *
 * This is an implementation header: there should be no need to use it outside
 * this directory.
 * Users should only use gmx_ana_selcollection_evaluate() declared in
 * \ref selection.h to evaluate selections.
 *
 * The functions defined in this header file are all the possible values
 * for the \c t_selelem::evaluate field (in addition to NULL).
 */
#ifndef SELECTION_EVALUATE_H
#define SELECTION_EVALUATE_H

#include <typedefs.h>

#include <indexutil.h>

#include "selelem.h"

struct gmx_sel_mempool_t;

/*! \internal \brief
 * Data structure for passing information required during evaluation.
 */
typedef struct gmx_sel_evaluate_t
{
    /** Memory pool for intermediate values. */
    struct gmx_sel_mempool_t *mp;
    /** Index group that contains all the atoms. */
    gmx_ana_index_t         *gall;
    /** Topology information. */
    t_topology              *top;
    /** Current frame. */
    t_trxframe              *fr;
    /** PBC data. */
    t_pbc                   *pbc;
} gmx_sel_evaluate_t;

/*! \name Utility functions
 */
/*@{*/
/** Initializes an evaluation data structure. */
void
_gmx_sel_evaluate_init(gmx_sel_evaluate_t *data,
                       struct gmx_sel_mempool_t *mp, gmx_ana_index_t *gall,
                       t_topology *top, t_trxframe *fr, t_pbc *pbc);
/** Evaluates the children of a general selection element. */
int
_gmx_sel_evaluate_children(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates the children of a \ref SEL_EXPRESSION element. */
int
_gmx_sel_evaluate_method_params(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/*@}*/

/*! \name Misc. evaluation functions
 */
/*@{*/
/** Evaluates a root selection element. */
int
_gmx_sel_evaluate_root(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a static group selection element. */
int
_gmx_sel_evaluate_static(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates an arithmetic expression element. */
int
_gmx_sel_evaluate_arithmetic(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/*@}*/

/*! \name Subexpression evaluation functions
 */
/*@{*/
/** Evaluates a subexpression when there is only one reference. */
int
_gmx_sel_evaluate_subexpr_simple(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a subexpression when the evaluation group is static. */
int
_gmx_sel_evaluate_subexpr_staticeval(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a subexpression. */
int
_gmx_sel_evaluate_subexpr(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a subexpression reference when there are no other references. */
int
_gmx_sel_evaluate_subexprref_simple(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a subexpression reference. */
int
_gmx_sel_evaluate_subexprref(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/*@}*/

/*! \name Method evaluation functions
 */
/*@{*/

/** Evaluates a method expression. */
int
_gmx_sel_evaluate_method(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a modifier expression. */
int
_gmx_sel_evaluate_modifier(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/*@}*/

/*! \name Boolean evaluation functions
 */
/*@{*/
/** Evaluates a boolean NOT element. */
int
_gmx_sel_evaluate_not(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a boolean AND element with short-circuiting. */
int
_gmx_sel_evaluate_and(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/** Evaluates a boolean OR element with short-circuiting. */
int
_gmx_sel_evaluate_or(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g);
/*@}*/

#endif
