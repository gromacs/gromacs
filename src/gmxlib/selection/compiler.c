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
 * \brief Selection compilation and optimization.
 *
 * \todo
 * Better error handling and memory management in error situations.
 * For example, memory handling in atom-valued method parameters within common
 * subexpressions may currently result in memory leaks.
 * Also, the main compilation function leaves the selection collection in
 * a bad state if an error occurs.
 *
 * \todo
 * The memory usage could be optimized.
 */
/*! \internal
 * \page selcompiler Selection compilation
 *
 * The compiler takes the selection element tree from the selection parser
 * (see \ref selparser) as input. The selection parser is quite independent of
 * selection evaluation details, and the compiler processes the tree to
 * conform to what the evaluation functions expect.
 * For better control and optimization possibilities, the compilation is
 * done on all selections simultaneously.
 * Hence, all the selections should be parsed before the compiler can be
 * called.
 *
 * The compiler initializes all fields in \c t_selelem not initialized by
 * the parser: \c t_selelem::v (some fields have already been initialized by
 * the parser), \c t_selelem::evaluate, and \c t_selelem::u (again, some
 * elements have been initialized in the parser).
 * The \c t_selelem::cdata field is used during the compilation to store
 * internal data, but the data is freed when the compiler returns.
 *
 * In addition to initializing the elements, the compiler reorganizes the tree
 * to simplify and optimize evaluation. The compiler also evaluates the static
 * parts of the selection: in the end of the compilation, static parts have
 * been replaced by the result of the evaluation.
 *
 * The compiler is called by calling gmx_ana_selcollection_compile().
 * This functions then does the compilation in several passes over the
 * \c t_selelem tree.
 *  -# Subexpressions are extracted: a separate root is created for each
 *     subexpression, and placed before the expression is first used.
 *     Currently, only variables and expressions used to evaluate parameter
 *     values are extracted, but common subexpression could also be detected
 *     here.
 *  -# A second pass with simple reordering and initialization is done:
 *    -# Boolean expressions are combined such that one element can evaluate,
 *       e.g., "A and B and C". The subexpressions in boolean expression are
 *       reordered such that static expressions come first without otherwise
 *       altering the relative order of the expressions.
 *    -# The \c t_selelem::evaluate field is set to the correct evaluation
 *       function from evaluate.h.
 *    -# The compiler data structure is allocated for each element, and
 *       the fields are initialized, with the exception of the contents of
 *       \c gmax and \c gmin fields.
 *    .
 *  -# The evaluation function of all elements is replaced with the
 *     analyze_static() function to be able to initialize the element before
 *     the actual evaluation function is called.
 *     The evaluation machinery is then called to initialize the whole tree,
 *     while simultaneously evaluating the static expressions.
 *     During the evaluation, track is kept of the smallest and largest
 *     possible selections, and these are stored in the internal compiler
 *     data structure for each element.
 *     To be able to do this for all possible values of dynamical expressions,
 *     special care needs to be taken with boolean expressions because they
 *     are short-circuiting. This is done through the
 *     \c t_compiler_data::bEvalMax flag, which makes dynamic child expressions
 *     of \c BOOL_OR expressions evaluate to empty groups, while subexpressions
 *     of \c BOOL_AND are evaluated to largest possible groups.
 *     Memory is also allocated to store the results of the evaluation.
 *     For each element, analyze_static() calls the actual evaluation function
 *     after the element has been properly initialized.
 *  -# Another evaluation pass is done over subexpressions with more than
 *     one reference to them. These cannot be completely processed during the
 *     first pass, because it is not known whether later references require
 *     additional evaluation of static expressions.
 *  -# Most of the processing is now done, and the next pass simply sets the
 *     evaluation group of root elements to the largest selection as determined
 *     in pass 3. Subexpressions that were evaluated to constants are no
 *     longer referenced at this time, and are removed.
 *  -# The next pass eliminates some unnecessary evaluation calls from
 *     subexpressions that are referenced only once, as well as initializing
 *     the position calculation data for selection method elements that require
 *     it. Compiler data is also freed as it is no longer needed.
 *  -# A final pass initializes the total masses and charges in the
 *     \c gmx_ana_selection_t data structures.
 *
 * The actual evaluation of the selection is described in the documentation
 * of the functions in evaluate.h.
 *
 * \todo
 * Some combinations of method parameter flags are not yet properly treated by
 * the compiler or the evaluation functions in evaluate.c. All the ones used by
 * currently implemented methods should work, but new combinations might not.
 *
 *
 * \section selcompiler_tree Element tree after compilation
 *
 * After the compilation, the selection element tree is suitable for
 * gmx_ana_selcollection_evaluate().
 * Enough memory has been allocated for \ref t_selelem::v
 * (and \ref t_selelem::cgrp for \ref SEL_SUBEXPR elements) to allow the
 * selection to be evaluated without allocating any memory.
 *
 *
 * \subsection selcompiler_tree_root Root elements
 *
 * The top level of the tree consists of a chain of \ref SEL_ROOT elements.
 * These are used for two purposes:
 *  -# A selection that should be evaluated.
 *     These elements appear in the same order as the selections in the input.
 *     For these elements, \ref t_selelem::v has been set to the maximum
 *     possible group that the selection can evaluate to, and
 *     \ref t_selelem::cgrp has been set to use a NULL group for evaluation.
 *  -# A subexpression that appears in one or more selections.
 *     Each selection that gives a value for a method parameter is a
 *     potential subexpression, as is any variable value.
 *     Only subexpressions that require evaluation for each frame are left
 *     after the selection is compiled.
 *     Each subexpression appears in the chain before any references to it.
 *     For these elements, \c t_selelem::cgrp has been set to the group
 *     that should be used to evaluate the subexpression.
 *     If \c t_selelem::cgrp is empty, the total evaluation group is not known
 *     in advance. If this is the case, \c t_selelem::evaluate is also NULL.
 *
 * The children of the \ref SEL_ROOT elements can be used to distinguish
 * the two types of root elements from each other; the rules are the same
 * as for the parsed tree (see \ref selparser_tree_root).
 * Subexpressions are treated as if they had been provided through variables.
 *
 * Selection names are stored as after parsing (see \ref selparser_tree_root).
 *
 *
 * \subsection selcompiler_tree_const Constant elements
 *
 * All (sub)selections that do not require particle positions have been
 * replaced with \ref SEL_CONST elements.
 * Constant elements from the parser are also retained if present in
 * dynamic parts of the selections.
 * Several constant elements with a NULL \c t_selelem::evaluate are left for
 * debugging purposes; of these, only the ones for \ref BOOL_OR expressions are
 * used during evaluation.
 *
 * The value is stored in \c t_selelem::v, and for group values with an
 * evaluation function set, also in \c t_selelem::cgrp.
 * For \ref GROUP_VALUE elements, unnecessary atoms (i.e., atoms that
 * could never be selected) have been removed from the value.
 *
 * \ref SEL_CONST elements have no children.
 *
 *
 * \subsection selcompiler_tree_method Method evaluation elements
 *
 * All selection methods that need to be evaluated dynamically are described
 * by a \ref SEL_EXPRESSION element. The \c t_selelem::method and
 * \c t_selelem::mdata fields have already been initialized by the parser,
 * and the compiler only calls the initialization functions in the method
 * data structure to do some additional initialization of these fields at
 * appropriate points. If the \c t_selelem::pc data field has been created by
 * the parser, the compiler initializes the data structure properly once the
 * required positions are known. If the \c t_selelem::pc field is NULL after
 * the parser, but the method provides only sel_updatefunc_pos(), an
 * appropriate position calculation data structure is created.
 * If \c t_selelem::pc is not NULL, \c t_selelem::pos is also initialized
 * to hold the positions calculated.
 *
 * Children of these elements are of type \ref SEL_SUBEXPRREF, and describe
 * parameter values that need to be evaluated for each frame. See the next
 * section for more details.
 * \ref SEL_CONST children can also appear, and stand for parameters that get
 * their value from a static expression. These elements are present only for
 * debugging purposes: they always have a NULL evaluation function.
 *
 *
 * \subsection selcompiler_tree_subexpr Subexpression elements
 *
 * As described in \ref selcompiler_tree_root, subexpressions are created
 * for each variable and each expression that gives a value to a selection
 * method parameter. As the only child of the \ref SEL_ROOT element,
 * these elements have a \ref SEL_SUBEXPR element. The \ref SEL_SUBEXPR
 * element has a single child, which evaluates the actual expression.
 * After compilation, only subexpressions that require particle positions
 * for evaluation are left.
 * For non-variable subexpression, automatic names have been generated to
 * help in debugging.
 *
 * For \ref SEL_SUBEXPR elements, memory has been allocated for
 * \c t_selelem::cgrp to store the group for which the expression has been
 * evaluated during the current frame.
 *
 * \ref SEL_SUBEXPRREF elements are used to describe references to
 * subexpressions. They have always a single child, which is the
 * \ref SEL_SUBEXPR element being referenced.
 *
 * If a subexpression is used only once and can be evaluated statically,
 * the evaluation has been optimized by setting the child of the
 * \ref SEL_SUBEXPR element to evaluate the value of \ref SEL_SUBEXPRREF
 * directly. In this case, the evaluation routines for the \ref SEL_SUBEXPRREF
 * and \ref SEL_SUBEXPR elements only propagate some status information,
 * but do not unnecessarily copy the values.
 *
 *
 * \subsection selcompiler_tree_bool Boolean elements
 *
 * \ref SEL_BOOLEAN elements have been merged such that one element
 * may carry out evaluation of more than one operation of the same type.
 * The static parts of the expressions have been evaluated, and are placed
 * in the first child. These are followed by the dynamic expressions, in the
 * order provided by the user.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <stdarg.h>

#include <smalloc.h>
#include <string2.h>
#include <vec.h>

#include <indexutil.h>
#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>

#include "evaluate.h"
#include "keywords.h"
#include "selcollection.h"
#include "selelem.h"

/*! \internal \brief
 * Internal data structure used by the compiler.
 */
typedef struct t_compiler_data
{
    /** The real evaluation method. */
    sel_evalfunc  evaluate; 
    /*! \brief
     * Whether the element is a method parameter.
     *
     * This flag is set for \ref SEL_SUBEXPR elements that are used to
     * evaluate non-atom-valued selection method parameters.
     */
    bool          bMethodParam;
    /*! \brief
     * TRUE if the whole subexpression should be treated static.
     *
     * This flag is always FALSE if \ref SEL_DYNAMIC is set for the element,
     * but it is also FALSE for static elements within common subexpressions.
     */
    bool          bStatic;
    /** TRUE if the subexpression will always be evaluated with the same group. */
    bool          bStaticEval;
    /** TRUE if the compiler evaluation routine should return the maximal selection. */
    bool          bEvalMax;
    /** Smallest selection that can be selected by the subexpression. */
    gmx_ana_index_t *gmin;
    /** Largest selection that can be selected by the subexpression. */
    gmx_ana_index_t *gmax;
    /** TRUE if memory has been allocated for \p gmin and \p gmax. */
    bool             bMinMaxAlloc;
} t_compiler_data;


/********************************************************************
 * COMPILER UTILITY FUNCTIONS
 ********************************************************************/

/*!
 * \param  sel Selection to free.
 *
 * This function only frees the data for the given selection, not its children.
 * It is safe to call the function when compiler data has not been allocated
 * or has already been freed; in such a case, nothing is done.
 */
void
_gmx_selelem_free_compiler_data(t_selelem *sel)
{
    if (sel->cdata)
    {
        sel->evaluate = sel->cdata->evaluate;
        if (sel->cdata->bMinMaxAlloc)
        {
            sel->cdata->gmin->name = NULL;
            sel->cdata->gmax->name = NULL;
            gmx_ana_index_deinit(sel->cdata->gmin);
            gmx_ana_index_deinit(sel->cdata->gmax);
            sfree(sel->cdata->gmin);
            sfree(sel->cdata->gmax);
        }
        sfree(sel->cdata);
    }
    sel->cdata = NULL;
}

/*! \brief
 * Allocates memory for storing the evaluated value of a selection element.
 *
 * \param     sel   Selection element to initialize
 * \param[in] isize Maximum evaluation group size.
 * \param[in] bChildEval TRUE if children have already been processed.
 * \returns   TRUE if the memory was allocated, FALSE if children need to
 *   be processed first.
 *
 * If called more than once, memory is (re)allocated to ensure that the
 * maximum of the \p isize values can be stored.
 */
static bool
alloc_selection_data(t_selelem *sel, int isize, bool bChildEval)
{
    int        nalloc;

    /* Find out the number of elements to allocate */
    if (sel->flags & SEL_SINGLEVAL)
    {
        nalloc = 1;
    }
    else if (sel->flags & SEL_ATOMVAL)
    {
        nalloc = isize;
    }
    else /* sel->flags should contain SEL_VARNUMVAL */
    {
        t_selelem *child;

        if (!bChildEval)
        {
            return FALSE;
        }
        child = (sel->type == SEL_SUBEXPRREF ? sel->child : sel);
        if (child->type == SEL_SUBEXPR)
        {
            child = child->child;
        }
        nalloc = (sel->v.type == POS_VALUE) ? child->v.u.p->nr : child->v.nr;
    }
    /* For positions, we actually want to allocate just a single structure
     * for nalloc positions. */
    if (sel->v.type == POS_VALUE)
    {
        isize  = nalloc;
        nalloc = 1;
    }
    /* Allocate memory for sel->v.u if needed */
    if ((sel->flags & SEL_ALLOCVAL)
        || (sel->type == SEL_SUBEXPRREF && sel->u.param
            && (sel->u.param->flags & (SPAR_VARNUM | SPAR_ATOMVAL))))
    {
        _gmx_selvalue_reserve(&sel->v, nalloc);
    }
    /* Reserve memory inside group and position structures if
     * SEL_ALLOCDATA is set. */
    if (sel->flags & SEL_ALLOCDATA)
    {
        if (sel->v.type == GROUP_VALUE)
        {
            gmx_ana_index_reserve(sel->v.u.g, isize);
        }
        else if (sel->v.type == POS_VALUE)
        {
            gmx_ana_pos_reserve(sel->v.u.p, isize, 0);
        }
    }
    return TRUE;
}

/*! \brief
 * Replace the evaluation function of each element in the subtree.
 *
 * \param     sel  Root of the selection subtree to process.
 * \param[in] eval The new evaluation function.
 */
static void
set_evaluation_function(t_selelem *sel, sel_evalfunc eval)
{
    sel->evaluate = eval;
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem *child = sel->child;
        while (child)
        {
            set_evaluation_function(child, eval);
            child = child->next;
        }
    }
}

/********************************************************************
 * SUBEXPRESSION EXTRACTION COMPILER PASS
 ********************************************************************/

/*! \brief
 * Creates a name with a running number for a subexpression.
 *
 * \param[in,out] sel The subexpression to be named.
 * \param[in]     i   Running number for the subexpression.
 *
 * The name of the selection becomes "SubExpr N", where N is \p i;
 * Memory is allocated for the name and the name is stored both in
 * \c t_selelem::name and \c t_selelem::u::cgrp::name; the latter
 * is freed by _gmx_selelem_free().
 */
static void
create_subexpression_name(t_selelem *sel, int i)
{
    int   len, ret;
    char *name;

    len = 8 + (int)log10(abs(i)) + 3;
    snew(name, len+1);
    /* FIXME: snprintf used to be used here for extra safety, but this
     * requires extra checking on Windows since it only provides a
     * non-C99-conforming implementation as _snprintf()... */
    ret = sprintf(name, "SubExpr %d", i);
    if (ret < 0 || ret > len)
    {
        sfree(name);
        name = NULL;
    }
    sel->name        = name;
    sel->u.cgrp.name = name;
}

/*! \brief
 * Processes and extracts subexpressions from a given selection subtree.
 *
 * \param   sel      Root of the subtree to process.
 * \param[in] gall     Index group that contains all the input atoms.
 * \param     subexprn Pointer to a subexpression counter.
 * \returns Pointer to a chain of subselections, or NULL if none were found.
 *
 * This function finds recursively all \ref SEL_SUBEXPRREF elements below
 * the given root element and ensures that their children are within
 * \ref SEL_SUBEXPR elements. It also creates a chain of \ref SEL_ROOT elements
 * that contain the subexpression as their children and returns the first
 * of these root elements.
 */
static t_selelem *
extract_item_subselections(t_selelem *sel, gmx_ana_index_t *gall, int *subexprn)
{
    t_selelem *root;
    t_selelem *subexpr;
    t_selelem *child;

    root = subexpr = NULL;
    child = sel->child;
    while (child)
    {
        if (!root)
        {
            root = subexpr = extract_item_subselections(child, gall, subexprn);
        }
        else
        {
            subexpr->next = extract_item_subselections(child, gall, subexprn);
        }
        while (subexpr && subexpr->next)
        {
            subexpr = subexpr->next;
        }
        /* The latter check excludes variable references */
        if (child->type == SEL_SUBEXPRREF && child->child->type != SEL_SUBEXPR)
        {
            /* Create the root element for the subexpression */
            if (!root)
            {
                root = subexpr = _gmx_selelem_create(SEL_ROOT);
            }
            else
            {
                subexpr->next = _gmx_selelem_create(SEL_ROOT);
                subexpr       = subexpr->next;
            }
            /* Set the evaluation group to all atoms */
            if (!(child->flags & SEL_ATOMVAL))
            {
                gmx_ana_index_set(&subexpr->u.cgrp, gall->isize, gall->index, NULL, 0);
            }
            /* Create the subexpression element */
            subexpr->child = _gmx_selelem_create(SEL_SUBEXPR);
            _gmx_selelem_set_vtype(subexpr->child, child->v.type);
            create_subexpression_name(subexpr->child, ++*subexprn);
            /* Move the actual subexpression under the created element */
            subexpr->child->child    = child->child;
            child->child             = subexpr->child;
            subexpr->child->refcount = 2;
            /* Set the flags for the created elements */
            subexpr->flags          |= (child->flags & SEL_VALFLAGMASK);
            subexpr->child->flags   |= (child->flags & SEL_VALFLAGMASK);
        }
        child = child->next;
    }

    return root;
}

/*! \brief
 * Extracts subexpressions of the selection chain.
 * 
 * \param   sel First selection in the whole selection chain.
 * \param[in] gall Index group that contains all the input atoms.
 * \returns The new first element for the chain.
 *
 * Finds all the subexpressions (and their subexpressions) in the
 * selection chain starting from \p sel and creates \ref SEL_SUBEXPR
 * elements for them.
 * \ref SEL_ROOT elements are also created for each subexpression
 * and inserted into the selection chain before the expressions that
 * refer to them.
 */
static t_selelem *
extract_subexpressions(t_selelem *sel, gmx_ana_index_t *gall)
{
    t_selelem   *root, *item, *next;
    int          subexprn;

    subexprn = 0;
    root = NULL;
    next = sel;
    while (next)
    {
        item = extract_item_subselections(next, gall, &subexprn);
        if (item)
        {
            if (!root)
            {
                root = item;
            }
            else
            {
                sel->next = item;
            }
            while (item->next)
            {
                item = item->next;
            }
            item->next = next;
        }
        else if (!root)
        {
            root = next;
        }
        sel = next;
        next = next->next;
    }
    return root;
}

/********************************************************************
 * BOOLEAN OPERATION REORDERING COMPILER PASS
 ********************************************************************/

/*! \brief
 * Removes redundant boolean selection elements.
 *
 * \param  sel Root of the selection subtree to optimize.
 *
 * This function merges similar boolean operations (e.g., (A or B) or C becomes
 * a single OR operation with three operands).
 */
static void
optimize_boolean_expressions(t_selelem *sel)
{
    t_selelem *child, *prev;

    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        prev  = NULL;
        child = sel->child;
        while (child)
        {
            optimize_boolean_expressions(child);
            /* Remove double negations */
            if (child->type == SEL_BOOLEAN && child->u.boolt == BOOL_NOT
                && child->child->type == SEL_BOOLEAN && child->child->u.boolt == BOOL_NOT)
            {
                /* Move the doubly negated expression up two levels */
                if (!prev)
                {
                    sel->child = child->child->child;
                    prev       = sel->child;
                }
                else
                {
                    prev->next = child->child->child;
                    prev       = prev->next;
                }
                child->child->child->next = child->next;
                /* Remove the two negations */
                child->child->child = NULL;
                child->next         = NULL;
                _gmx_selelem_free(child);
                child = prev;
            }
            prev  = child;
            child = child->next;
        }
    }
    if (sel->type != SEL_BOOLEAN || sel->u.boolt == BOOL_NOT)
    {
        return;
    }
    /* Merge subsequent binary operations */
    prev  = NULL;
    child = sel->child;
    while (child)
    {
        if (child->type == SEL_BOOLEAN && child->u.boolt == sel->u.boolt)
        {
            if (!prev)
            {
                sel->child = child->child;
                prev       = sel->child;
            }
            else
            {
                prev->next = child->child;
            }
            while (prev->next)
            {
                prev = prev->next;
            }
            prev->next = child->next;
            sfree(child->v.u.g);
            sfree(child);
            child = prev->next;
        }
        else
        {
            prev = child;
            child = child->next;
        }
    }
}

/*! \brief
 * Reorders children of boolean expressions such that static selections
 * come first.
 *
 * \param  sel Root of the selection subtree to reorder.
 *
 * The relative order of static expressions does not change.
 * The same is true for the dynamic expressions.
 */
static void
reorder_boolean_static_children(t_selelem *sel)
{
    t_selelem *child, *prev, *next;

    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            reorder_boolean_static_children(child);
            child = child->next;
        }
    }

    /* Reorder boolean expressions such that static selections come first */
    if (sel->type == SEL_BOOLEAN && (sel->flags & SEL_DYNAMIC))
    {
        t_selelem  start;

        start.next = sel->child;
        prev  = &start;
        child = &start;
        while (child->next)
        {
            /* child is the last handled static expression */
            /* prev is the last handled non-static expression */
            next = prev->next;
            while (next && (next->flags & SEL_DYNAMIC))
            {
                prev = next;
                next = next->next;
            }
            /* next is now the first static expression after child */
            if (!next)
            {
                break;
            }
            /* Reorder such that next comes after child */
            if (prev != child)
            {
                prev->next  = next->next;
                next->next  = child->next;
                child->next = next;
            }
            else
            {
                prev = prev->next;
            }
            /* Advance child by one */
            child = next;
        }

        sel->child = start.next;
    }
}

/********************************************************************
 * EVALUATION PREPARATION COMPILER PASS
 ********************************************************************/

/*! \brief
 * Initializes the evaluation groups for the selections.
 *
 * \param[in,out] sc   Selection collection data.
 *
 * The evaluation group of each \ref SEL_ROOT element corresponding to a
 * selection in \p sc is set to \p gall.
 */
static void
initialize_evalgrps(gmx_ana_selcollection_t *sc)
{
    t_selelem   *item;
    int          i;

    /* Initialize the output */
    for (i = 0; i < sc->nr; ++i)
    {
        item = sc->sel[i]->selelem;
        /* Set the evaluation group to all atoms */
        gmx_ana_index_set(&item->u.cgrp, sc->gall.isize, sc->gall.index,
                          item->u.cgrp.name, 0);
    }
}

/*! \brief
 * Prepares the selection (sub)tree for evaluation.
 *
 * \param[in,out] sel Root of the selection subtree to prepare.
 * \returns       TRUE on success, FALSE if any subexpression fails.
 *
 * This function sets the evaluation function (\c t_selelem::evaluate)
 * for the selection elements.
 * It also allocates memory for the \p sel->v.u.g or \p sel->v.u.p
 * structure if required.
 */
static bool
init_item_evaluation(t_selelem *sel)
{
    t_selelem         *child;

    /* Process children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            if (!init_item_evaluation(child))
            {
                return FALSE;
            }
            child = child->next;
        }
    }

    /* Make sure that the group/position structure is allocated */
    if (!sel->v.u.ptr && (sel->flags & SEL_ALLOCVAL))
    {
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_reserve(&sel->v, 1);
            sel->v.nr = 1;
        }
    }

    /* Set the evaluation function */
    switch (sel->type)
    {
        case SEL_CONST:
            if (sel->v.type == GROUP_VALUE)
            {
                sel->evaluate = &_gmx_sel_evaluate_static;
            }
            break;

        case SEL_EXPRESSION:
            sel->evaluate = &_gmx_sel_evaluate_method;
            break;

        case SEL_MODIFIER:
            if (sel->v.type != NO_VALUE)
            {
                sel->evaluate = &_gmx_sel_evaluate_modifier;
            }
            break;

        case SEL_BOOLEAN:
            switch (sel->u.boolt)
            {
                case BOOL_NOT: sel->evaluate = &_gmx_sel_evaluate_not; break;
                case BOOL_AND: sel->evaluate = &_gmx_sel_evaluate_and; break;
                case BOOL_OR:  sel->evaluate = &_gmx_sel_evaluate_or;  break;
                case BOOL_XOR:
                    gmx_impl("xor expressions not implemented");
                    return FALSE;
            }
            break;

        case SEL_ROOT:
            sel->evaluate = &_gmx_sel_evaluate_root;
            break;

        case SEL_SUBEXPR:
            sel->evaluate = &_gmx_sel_evaluate_subexpr;
            break;

        case SEL_SUBEXPRREF:
            sel->name     = sel->child->name;
            sel->evaluate = &_gmx_sel_evaluate_subexprref;
            break;
    }

    return TRUE;
}


/********************************************************************
 * COMPILER DATA INITIALIZATION PASS
 ********************************************************************/

/*! \brief
 * Allocates memory for the compiler data and initializes the structure.
 *
 * \param sel Root of the selection subtree to process.
 */
static void
init_item_compilerdata(t_selelem *sel)
{
    t_selelem   *child;

    /* Allocate the compiler data structure */
    snew(sel->cdata, 1);

    /* Store the real evaluation method because the compiler will replace it */
    sel->cdata->evaluate = sel->evaluate;

    /* Initialize the flags */
    sel->cdata->bMethodParam = FALSE;
    sel->cdata->bStatic      = !(sel->flags & SEL_DYNAMIC);
    sel->cdata->bStaticEval  = TRUE;
    sel->cdata->bEvalMax     = (sel->type == SEL_SUBEXPR ? TRUE : FALSE);
    /* Set the method parameter flag for non-atom-valued parameters */
    if (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER)
    {
        child = sel->child;
        while (child)
        {
            if (!(child->flags & SEL_ATOMVAL))
            {
                child->child->cdata->bMethodParam = TRUE;
            }
            child = child->next;
        }
    }

    /* Initialize children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            init_item_compilerdata(child);
            child = child->next;
        }
    }

    /* Determine whether we should evaluate the minimum or the maximum
     * for the children of this element. */
    if (sel->type == SEL_BOOLEAN)
    {
        bool  bEvalMax;

        bEvalMax = (sel->u.boolt == BOOL_AND);
        child = sel->child;
        while (child)
        {
            child->cdata->bEvalMax = bEvalMax;
            if (child->type == SEL_BOOLEAN && child->u.boolt == BOOL_NOT)
            {
                child->child->cdata->bEvalMax = !bEvalMax;
            }
            child = child->next;
        }
    }
    else if (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER
             || sel->type == SEL_SUBEXPR)
    {
        child = sel->child;
        while (child)
        {
            child->cdata->bEvalMax = TRUE;
            child = child->next;
        }
    }

    /* Initialize the minimum and maximum evaluation groups */
    sel->cdata->bMinMaxAlloc = FALSE;
    if (sel->type != SEL_ROOT && sel->v.type != NO_VALUE)
    {
        if (sel->type == SEL_SUBEXPR)
        {
            sel->cdata->gmin = sel->child->cdata->gmin;
            sel->cdata->gmax = sel->child->cdata->gmax;
        }
        else if (sel->v.type == GROUP_VALUE && sel->cdata->bStatic)
        {
            sel->cdata->gmin = sel->v.u.g;
            sel->cdata->gmax = sel->v.u.g;
        }
        else
        {
            sel->cdata->bMinMaxAlloc = TRUE;
            snew(sel->cdata->gmin, 1);
            snew(sel->cdata->gmax, 1);
        }
    }
}

/*! \brief
 * Initializes the static evaluation flag for a selection subtree.
 *
 * \param[in,out] sel  Root of the selection subtree to process.
 *
 * Sets the \c bStaticEval in the compiler data structure:
 * for any element for which the evaluation group may depend on the trajectory
 * frame, the flag is cleared.
 *
 * reorder_boolean_static_children() should have been called.
 */
static void
init_item_staticeval(t_selelem *sel)
{
    t_selelem   *child;

    /* Non-atom-valued method parameters should always have bStaticEval,
     * so don't do anything if a reference to them is encountered. */
    if (sel->type == SEL_SUBEXPRREF && sel->child->cdata->bMethodParam)
    {
        return;
    }

    /* Propagate the bStaticEval flag to children if it is not set */
    if (!sel->cdata->bStaticEval)
    {
        child = sel->child;
        while (child)
        {
            if ((sel->type != SEL_EXPRESSION && sel->type != SEL_MODIFIER)
                || (child->flags & SEL_ATOMVAL))
            {
                if (child->cdata->bStaticEval)
                {
                    child->cdata->bStaticEval = FALSE;
                    init_item_staticeval(child);
                }
            }
            child = child->next;
        }
    }
    else /* bStaticEval is set */
    {
        /* For boolean expressions, any expression after the first dynamic
         * expression should not have bStaticEval. */
        if (sel->type == SEL_BOOLEAN)
        {
            child = sel->child;
            while (child && !(child->flags & SEL_DYNAMIC))
            {
                child = child->next;
            }
            if (child)
            {
                child = child->next;
            }
            while (child)
            {
                child->cdata->bStaticEval = FALSE;
                child = child->next;
            }
        }

        /* Process the children */
        child = sel->child;
        while (child)
        {
            init_item_staticeval(child);
            child = child->next;
        }
    }
}


/********************************************************************
 * STATIC ANALYSIS COMPILER PASS
 ********************************************************************/

/*! \brief
 * Marks a subtree completely dynamic or undoes such a change.
 *
 * \param     sel      Selection subtree to mark.
 * \param[in] bDynamic If TRUE, the \p bStatic flag of the whole
 *   selection subtree is cleared. If FALSE, the flag is restored to
 *   using \ref SEL_DYNAMIC.
 *
 * Does not descend into parameters of methods unless the parameters
 * are evaluated for each atom.
 */
static void
mark_subexpr_dynamic(t_selelem *sel, bool bDynamic)
{
    t_selelem *child;

    sel->cdata->bStatic = (!bDynamic && !(sel->flags & SEL_DYNAMIC));
    child = sel->child;
    while (child)
    {
        if (sel->type != SEL_EXPRESSION || child->type != SEL_SUBEXPRREF
            || (child->u.param->flags & SPAR_ATOMVAL))
        {
            mark_subexpr_dynamic(child, bDynamic);
        }
        child = child->next;
    }
}

/*! \brief
 * Makes an evaluated selection element static.
 *
 * \param     sel   Selection element to make static.
 *
 * The evaluated value becomes the value of the static element.
 * The element type is changed to SEL_CONST and the children are
 * deleted.
 */
static void
make_static(t_selelem *sel)
{
    /* Free the expression data as it is no longer needed */
    _gmx_selelem_free_exprdata(sel);
    /* Make the item static */
    sel->name            = NULL;
    sel->type            = SEL_CONST;
    sel->evaluate        = NULL;
    sel->cdata->evaluate = NULL;
    /* Free the children */
    _gmx_selelem_free_chain(sel->child);
    sel->child           = NULL;
    /* Set the group value.
     * None of the elements for which this function may be called uses
     * the cgrp group, so we can simply overwrite the contents without
     * worrying about memory leaks. */
    if (sel->v.type == GROUP_VALUE)
    {
        gmx_ana_index_set(&sel->u.cgrp, sel->v.u.g->isize, sel->v.u.g->index, NULL, 0);
    }
}

/*! \brief
 * Evaluates a constant expression during analyze_static() and analyze_static2().
 *
 * \param[in]     data Evaluation data.
 * \param[in,out] sel Selection to process.
 * \param[in]     g   The evaluation group.
 * \returns       0 on success, a non-zero error code on error.
 */
static int
process_const(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g)
{
    int  rc;

    rc = 0;
    if (sel->v.type == GROUP_VALUE)
    {
        if (sel->cdata->evaluate)
        {
            rc = sel->cdata->evaluate(data, sel, g);
        }
    }
    /* Other constant expressions do not need evaluation */
    return rc;
}

/*! \brief
 * Sets the parameter value pointer for \ref SEL_SUBEXPRREF params.
 *
 * \param[in,out] sel Selection to process.
 *
 * Copies the value pointer of \p sel to \c sel->u.param if one is present
 * and should receive the value from the compiler
 * (most parameter values are handled during parsing).
 * If \p sel is not of type \ref SEL_SUBEXPRREF, or if \c sel->u.param is NULL,
 * the function does nothing.
 * Also, if the \c sel->u.param does not have \ref SPAR_VARNUM or
 * \ref SPAR_ATOMVAL, the function returns immediately.
 */
static void
store_param_val(t_selelem *sel)
{
    /* Return immediately if there is no parameter. */
    if (sel->type != SEL_SUBEXPRREF || !sel->u.param)
    {
        return;
    }

    /* Or if the value does not need storing. */
    if (!(sel->u.param->flags & (SPAR_VARNUM | SPAR_ATOMVAL)))
    {
        return;
    }

    if (sel->v.type == INT_VALUE || sel->v.type == REAL_VALUE
        || sel->v.type == STR_VALUE)
    {
        _gmx_selvalue_setstore(&sel->u.param->val, sel->v.u.ptr);
    }
}

/*! \brief
 * Handles the initialization of a selection method during analyze_static() pass.
 *
 * \param[in,out] sel Selection element to process.
 * \param[in]     top Topology structure.
 * \param[in]     isize Size of the evaluation group for the element.
 * \returns       0 on success, a non-zero error code on return.
 *
 * Calls sel_initfunc() (and possibly sel_outinitfunc()) to initialize the
 * method.
 * If no \ref SPAR_ATOMVAL parameters are present, multiple initialization
 * is prevented by using \ref SEL_METHODINIT and \ref SEL_OUTINIT flags.
 */
static int
init_method(t_selelem *sel, t_topology *top, int isize)
{
    t_selelem *child;
    bool       bAtomVal;
    int        rc;

    /* Find out whether there are any atom-valued parameters */
    bAtomVal = FALSE;
    child = sel->child;
    while (child)
    {
        if (child->flags & SEL_ATOMVAL)
        {
            bAtomVal = TRUE;
        }
        child = child->next;
    }

    /* Initialize the method */
    if (sel->u.expr.method->init
        && (bAtomVal || !(sel->flags & SEL_METHODINIT)))
    {
        sel->flags |= SEL_METHODINIT;
        /* The allocation flags are cleared first to not to free anything if
         * initialization fails. */
        child = sel->child;
        if (sel->type == SEL_MODIFIER && sel->v.type != NO_VALUE)
        {
            child = child->next;
        }
        while (child)
        {
            child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
            child = child->next;
        }
        rc = sel->u.expr.method->init(top, sel->u.expr.method->nparams,
                sel->u.expr.method->param, sel->u.expr.mdata);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (bAtomVal || !(sel->flags & SEL_OUTINIT))
    {
        sel->flags |= SEL_OUTINIT;
        if (sel->u.expr.method->outinit)
        {
            rc = sel->u.expr.method->outinit(top, &sel->v, sel->u.expr.mdata);
            if (rc != 0)
            {
                return rc;
            }
        }
        else
        {
            alloc_selection_data(sel, isize, TRUE);
            if ((sel->flags & SEL_DYNAMIC)
                && sel->v.type != GROUP_VALUE && sel->v.type != POS_VALUE)
            {
                sel->v.nr = isize;
            }
            /* If the method is char-valued, pre-allocate the strings. */
            if (sel->u.expr.method->flags & SMETH_CHARVAL)
            {
                int  i;

                /* A sanity check */
                if (sel->v.type != STR_VALUE)
                {
                    gmx_bug("internal error");
                    return -1;
                }
                sel->flags |= SEL_ALLOCDATA;
                for (i = 0; i < isize; ++i)
                {
                    if (sel->v.u.s[i] == NULL)
                    {
                        snew(sel->v.u.s[i], 2);
                    }
                }
            }
        }
    }

    return 0;
}

/*! \brief
 * Evaluates the static part of a boolean expression.
 *
 * \param[in]     data Evaluation data.
 * \param[in,out] sel Boolean selection element whose children should be
 *   processed.
 * \param[in]     g   The evaluation group.
 * \returns       0 on success, a non-zero error code on error.
 *
 * reorder_item_static_children() should have been called.
 */
static int
evaluate_boolean_static_part(gmx_sel_evaluate_t *data, t_selelem *sel,
                             gmx_ana_index_t *g)
{
    t_selelem *child, *next;
    int        rc;

    /* Find the last static subexpression */
    child = sel->child;
    while (child->next && child->next->cdata->bStatic)
    {
        child = child->next;
    }
    if (!child->cdata->bStatic)
    {
        return 0;
    }

    /* Evalute the static part if there is more than one expression */
    if (child != sel->child)
    {
        next  = child->next;
        child->next = NULL;
        rc = sel->cdata->evaluate(data, sel, g);
        if (rc != 0)
        {
            return rc;
        }
        /* Replace the subexpressions with the result */
        _gmx_selelem_free_chain(sel->child);
        snew(child, 1);
        child->type       = SEL_CONST;
        child->flags      = SEL_FLAGSSET | SEL_SINGLEVAL | SEL_ALLOCVAL | SEL_ALLOCDATA;
        _gmx_selelem_set_vtype(child, GROUP_VALUE);
        child->evaluate   = NULL;
        _gmx_selvalue_reserve(&child->v, 1);
        gmx_ana_index_copy(child->v.u.g, sel->v.u.g, TRUE);
        init_item_compilerdata(child);
        child->cdata->bStaticEval = sel->cdata->bStaticEval;
        child->next = next;
        sel->child = child;
    }
    else if (child->evaluate)
    {
        rc = child->evaluate(data, child, g);
        if (rc != 0)
        {
            return rc;
        }
    }
    /* Set the evaluation function for the constant element.
     * We never need to evaluate the element again during compilation,
     * but we may need to evaluate the static part again if the
     * expression is not an OR with a static evaluation group.
     * If we reach here with a NOT expression, the NOT expression
     * is also static, and will be made a constant later, so don't waste
     * time copying the group. */
    child->evaluate = NULL;
    if (sel->u.boolt == BOOL_NOT
        || (sel->cdata->bStaticEval && sel->u.boolt == BOOL_OR))
    {
        child->cdata->evaluate = NULL;
    }
    else
    {
        child->cdata->evaluate = &_gmx_sel_evaluate_static;
        /* The cgrp has only been allocated if it originated from an
         * external index group. In that case, we need special handling
         * to preserve the name of the group and to not leak memory.
         * If cgrp has been set in make_static(), it is not allocated,
         * and hence we can overwrite it safely. */
        if (child->u.cgrp.nalloc_index > 0)
        {
            char *name = child->u.cgrp.name;
            gmx_ana_index_copy(&child->u.cgrp, child->v.u.g, FALSE);
            gmx_ana_index_squeeze(&child->u.cgrp);
            child->u.cgrp.name = name;
        }
        else
        {
            gmx_ana_index_copy(&child->u.cgrp, child->v.u.g, TRUE);
        }
    }
    return 0;
}

/*! \brief
 * Evaluates the minimum and maximum groups for a boolean expression.
 *
 * \param[in]  sel  \ref SEL_BOOLEAN element currently being evaluated.
 * \param[in]  g    Group for which \p sel has been evaluated.
 * \param[out] gmin Largest subset of the possible values of \p sel.
 * \param[out] gmax Smallest superset of the possible values of \p sel.
 *
 * This is a helper function for analyze_static() that is called for
 * dynamic \ref SEL_BOOLEAN elements after they have been evaluated.
 * It uses the minimum and maximum groups of the children to calculate
 * the minimum and maximum groups for \p sel, and also updates the static
 * part of \p sel (which is in the first child) if the children give
 * cause for this.
 *
 * This function may allocate some extra memory for \p gmin and \p gmax,
 * but as these groups are freed at the end of analyze_static() (which is
 * reached shortly after this function returns), this should not be a major
 * problem.
 */
static void
evaluate_boolean_minmax_grps(t_selelem *sel, gmx_ana_index_t *g,
                             gmx_ana_index_t *gmin, gmx_ana_index_t *gmax)
{
    t_selelem *child;

    switch (sel->u.boolt)
    {
        case BOOL_NOT:
            gmx_ana_index_reserve(gmin, g->isize);
            gmx_ana_index_reserve(gmax, g->isize);
            gmx_ana_index_difference(gmax, g, sel->child->cdata->gmin);
            gmx_ana_index_difference(gmin, g, sel->child->cdata->gmax);
            break;

        case BOOL_AND:
            gmx_ana_index_copy(gmin, sel->child->cdata->gmin, TRUE);
            gmx_ana_index_copy(gmax, sel->child->cdata->gmax, TRUE);
            child = sel->child->next;
            while (child && gmax->isize > 0)
            {
                gmx_ana_index_intersection(gmin, gmin, child->cdata->gmin);
                gmx_ana_index_intersection(gmax, gmax, child->cdata->gmax);
                child = child->next;
            }
            /* Update the static part if other expressions limit it */
            if (sel->child->cdata->bStatic
                && sel->child->v.u.g->isize > gmax->isize)
            {
                gmx_ana_index_copy(sel->child->v.u.g, gmax, FALSE);
                gmx_ana_index_squeeze(sel->child->v.u.g);
                if (sel->child->u.cgrp.isize > 0)
                {
                    gmx_ana_index_copy(&sel->child->u.cgrp, gmax, FALSE);
                    gmx_ana_index_squeeze(&sel->child->u.cgrp);
                }
            }
            break;

        case BOOL_OR:
            /* We can assume here that the gmin of children do not overlap
             * because of the way _gmx_sel_evaluate_or() works. */
            gmx_ana_index_reserve(gmin, g->isize);
            gmx_ana_index_reserve(gmax, g->isize);
            gmx_ana_index_copy(gmin, sel->child->cdata->gmin, FALSE);
            gmx_ana_index_copy(gmax, sel->child->cdata->gmax, FALSE);
            child = sel->child->next;
            while (child && gmin->isize < g->isize)
            {
                gmx_ana_index_merge(gmin, gmin, child->cdata->gmin);
                gmx_ana_index_union(gmax, gmax, child->cdata->gmax);
                child = child->next;
            }
            /* Update the static part if other expressions have static parts
             * that are not included. */
            if (sel->child->cdata->bStatic
                && sel->child->v.u.g->isize < gmin->isize)
            {
                gmx_ana_index_reserve(sel->child->v.u.g, gmin->isize);
                gmx_ana_index_copy(sel->child->v.u.g, gmin, FALSE);
                if (sel->child->u.cgrp.isize > 0)
                {
                    gmx_ana_index_reserve(&sel->child->u.cgrp, gmin->isize);
                    gmx_ana_index_copy(&sel->child->u.cgrp, gmin, FALSE);
                }
            }
            break;

        case BOOL_XOR: /* Should not be reached */
            gmx_impl("xor expressions not implemented");
            break;
    }
}

/*! \brief
 * Evaluates the static parts of \p sel and analyzes the structure.
 * 
 * \param[in]     data Evaluation data.
 * \param[in,out] sel  Selection currently being evaluated.
 * \param[in]     g    Group for which \p sel should be evaluated.
 * \returns       0 on success, a non-zero error code on error.
 *
 * This function is used as the replacement for the \c t_selelem::evaluate
 * function pointer.
 * It does the single most complex task in the compiler: after all elements
 * have been processed, the \p gmin and \p gmax fields of \p t_compiler_data
 * have been properly initialized, enough memory has been allocated for
 * storing the value of each expression, and the static parts of the 
 * expressions have been evaluated.
 * The above is exactly true only for elements other than subexpressions:
 * another pass is required for subexpressions that are referred to more than
 * once to evaluate the static parts.
 * This second pass is performed by analyze_static2().
 *
 * \see analyze_static2()
 */
static int
analyze_static(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g)
{
    t_selelem       *child, *next;
    gmx_ana_index_t  gmin, gmax;
    bool             bDelayAlloc;
    int              rc;

    gmx_ana_index_clear(&gmin);
    gmx_ana_index_clear(&gmax);
    bDelayAlloc = FALSE;

    if (sel->type != SEL_ROOT && g)
    {
        bDelayAlloc = !alloc_selection_data(sel, g->isize, FALSE);
    }

    /* TODO: This switch is awfully long... */
    rc = 0;
    switch (sel->type)
    {
        case SEL_CONST:
            rc = process_const(data, sel, g);
            break;

        case SEL_EXPRESSION:
        case SEL_MODIFIER:
            rc = _gmx_sel_evaluate_method_params(data, sel, g);
            if (rc != 0)
            {
                return rc;
            }
            rc = init_method(sel, data->top, g->isize);
            if (rc != 0)
            {
                return rc;
            }
            if (!(sel->flags & SEL_DYNAMIC))
            {
                rc = sel->cdata->evaluate(data, sel, g);
                if (rc == 0 && sel->cdata->bStatic)
                {
                    make_static(sel);
                }
            }
            else
            {
                /* Modifiers need to be evaluated even though they process
                 * positions to get the modified output groups from the
                 * maximum possible selections. */
                if (sel->type == SEL_MODIFIER)
                {
                    rc = sel->cdata->evaluate(data, sel, g);
                }
                gmx_ana_index_copy(&gmax, g, TRUE);
            }
            break;

        case SEL_BOOLEAN:
            if (!(sel->flags & SEL_DYNAMIC))
            {
                rc = sel->cdata->evaluate(data, sel, g);
                if (rc == 0 && sel->cdata->bStatic)
                {
                    make_static(sel);
                }
            }
            else
            {
                /* Evalute the static part if there is more than one expression */
                rc = evaluate_boolean_static_part(data, sel, g);
                if (rc != 0)
                {
                    return rc;
                }

                /* Evaluate the selection.
                 * If the type is boolean, we must explicitly handle the
                 * static part evaluated in evaluate_boolean_static_part()
                 * here because g may be larger. */
                if (sel->u.boolt == BOOL_AND && sel->child->type == SEL_CONST)
                {
                    rc = sel->cdata->evaluate(data, sel, sel->child->v.u.g);
                }
                else
                {
                    rc = sel->cdata->evaluate(data, sel, g);
                }
                if (rc != 0)
                {
                    return rc;
                }

                /* Evaluate minimal and maximal selections */
                evaluate_boolean_minmax_grps(sel, g, &gmin, &gmax);
            }
            break;

        case SEL_ROOT:
            rc = sel->cdata->evaluate(data, sel, g);
            break;

        case SEL_SUBEXPR:
            if (sel->u.cgrp.isize == 0)
            {
                gmx_ana_index_reserve(&sel->u.cgrp, g->isize);
                if (bDelayAlloc)
                {
                    /* We need to evaluate the child before we can allocate the
                     * memory. */
                    rc = sel->child->evaluate(data, sel->child, g);
                    if (rc != 0)
                    {
                        return rc;
                    }
                    alloc_selection_data(sel, g->isize, TRUE);
                    /* Do not evaluate the child again */
                    sel->child->evaluate = NULL;
                    rc = sel->cdata->evaluate(data, sel, g);
                    sel->child->evaluate = &analyze_static;
                }
                else
                {
                    alloc_selection_data(sel->child, g->isize, FALSE);
                    rc = sel->cdata->evaluate(data, sel, g);
                }
            }
            else
            {
                int isize = gmx_ana_index_difference_size(g, &sel->u.cgrp);
                if (isize > 0)
                {
                    isize += sel->u.cgrp.isize;
                    gmx_ana_index_reserve(&sel->u.cgrp, isize);
                    if (sel->v.type == GROUP_VALUE || (sel->flags & SEL_ATOMVAL))
                    {
                        alloc_selection_data(sel->child, isize, FALSE);
                        alloc_selection_data(sel,        isize, FALSE);
                    }
                    rc = sel->cdata->evaluate(data, sel, g);
                }
                else
                {
                    rc = sel->cdata->evaluate(data, sel, g);
                }
            }
            break;

        case SEL_SUBEXPRREF:
            /* Evaluate the subexpression if it is not yet evaluated.
             * Can happen when a variable is passed as a parameter or as
             * a selection. */
            if (sel->child->u.cgrp.isize == 0)
            {
                rc = sel->child->evaluate(data, sel->child, g ? g : data->gall);
                if (rc != 0)
                {
                    return rc;
                }
                /* Prevent another evaluation of the child. */
                sel->child->evaluate = NULL;
                alloc_selection_data(sel, sel->child->cdata->gmax->isize, TRUE);
            }
            if (!g)
            {
                alloc_selection_data(sel, sel->child->cdata->gmax->isize, TRUE);
            }
            /* TODO: This is not general enough if/when position references
             * can be evaluated more than once (that is, if there are position
             * methods that do not have SMETH_SINGLEVAL or SMETH_VARNUMVAL). */
            if (sel->v.type == POS_VALUE && !(sel->flags & SEL_OUTINIT))
            {
                gmx_ana_indexmap_copy(&sel->v.u.p->m, &sel->child->child->v.u.p->m, TRUE);
                sel->flags |= SEL_OUTINIT;
            }
            rc = sel->cdata->evaluate(data, sel, g);
            sel->child->evaluate = &analyze_static;
            if (rc != 0)
            {
                return rc;
            }
            /* Store the parameter value if required */
            store_param_val(sel);
            if (!(sel->flags & SEL_DYNAMIC))
            {
                if (sel->cdata->bStatic)
                {
                    make_static(sel);
                }
            }
            else
            {
                if (sel->child->refcount <= 2 || !g)
                {
                    gmx_ana_index_copy(&gmin, sel->child->cdata->gmin, TRUE);
                    gmx_ana_index_copy(&gmax, sel->child->cdata->gmax, TRUE);
                }
                else
                {
                    gmx_ana_index_reserve(&gmin, min(g->isize, sel->child->cdata->gmin->isize));
                    gmx_ana_index_reserve(&gmax, min(g->isize, sel->child->cdata->gmax->isize));
                    gmx_ana_index_intersection(&gmin, sel->child->cdata->gmin, g);
                    gmx_ana_index_intersection(&gmax, sel->child->cdata->gmax, g);
                }
            }
            break;
    }
    /* Exit if there was some problem */
    if (rc != 0)
    {
        return rc;
    }

    /* Update the minimal and maximal evaluation groups */
    if (sel->cdata->bMinMaxAlloc)
    {
        gmx_ana_index_reserve(sel->cdata->gmin, sel->cdata->gmin->isize + gmin.isize);
        gmx_ana_index_reserve(sel->cdata->gmax, sel->cdata->gmax->isize + gmax.isize);
        gmx_ana_index_merge(sel->cdata->gmin, sel->cdata->gmin, &gmin);
        gmx_ana_index_merge(sel->cdata->gmax, sel->cdata->gmax, &gmax);
    }
    /* Replace the result of the evaluation */
    /* This is not necessary for subexpressions or for boolean negations
     * because the evaluation function already has done it properly. */
    if (sel->v.type == GROUP_VALUE && (sel->flags & SEL_DYNAMIC)
        && sel->type != SEL_SUBEXPR
        && !(sel->type == SEL_BOOLEAN && sel->u.boolt == BOOL_NOT))
    {
        if (sel->cdata->bEvalMax)
        {
            gmx_ana_index_copy(sel->v.u.g, &gmax, FALSE);
        }
        else
        {
            gmx_ana_index_copy(sel->v.u.g, &gmin, FALSE);
        }
    }
    gmx_ana_index_deinit(&gmin);
    gmx_ana_index_deinit(&gmax);

    /* Make sure that enough data storage has been allocated */
    /* TODO: Constant expressions could be handled here more intelligently */
    if (sel->type != SEL_ROOT && sel->cdata->bStaticEval)
    {
        alloc_selection_data(sel, sel->cdata->gmax->isize, TRUE);
        /* Make sure that the new value pointer is stored if required */
        store_param_val(sel);
    }

    return 0;
}

/*! \brief
 * Evaluates the static parts of \p sel and analyzes the structure.
 * 
 * \param[in]     data Evaluation data.
 * \param[in,out] sel  Selection currently being evaluated.
 * \param[in]     g   Group for which \p sel should be evaluated.
 * \returns       0 on success, a non-zero error code on error.
 *
 * This function is a simpler version of analyze_static() that is used
 * during a second evaluation round, and can thus use information calculated
 * by analyze_static().
 * It is also used as the replacement for the \c t_selelem::evaluate
 * function pointer.
 * It is used to evaluate the static parts of subexpressions that could not
 * be evaluated during the analyze_static() pass.
 *
 * \see analyze_static()
 */
static int
analyze_static2(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g)
{
    int rc;

    rc = 0;
    switch (sel->type)
    {
        case SEL_CONST:
            rc = process_const(data, sel, g);
            break;

        case SEL_EXPRESSION:
        case SEL_BOOLEAN:
        case SEL_SUBEXPRREF:
            if (sel->cdata->bStatic)
            {
                rc = sel->cdata->evaluate(data, sel, g);
                if (rc == 0)
                {
                    make_static(sel);
                }
            }
            else if (sel->type == SEL_BOOLEAN)
            {
                rc = evaluate_boolean_static_part(data, sel, g);
                if (rc == 0)
                {
                    rc = sel->cdata->evaluate(data, sel, g);
                }
            }
            break;

        case SEL_ROOT:     /* Roots should not be present here */
        case SEL_MODIFIER: /* Modifiers should not be present here */
        case SEL_SUBEXPR:
            rc = sel->cdata->evaluate(data, sel, g);
            break;
    }
    /* Exit if there was some problem */
    if (rc != 0)
    {
        return rc;
    }

    /* Replace the result of the evaluation */
    /* This is not necessary for subexpressions or for boolean negations
     * because the evaluation function already has done it properly. */
    if (sel->v.type == GROUP_VALUE && !sel->cdata->bStatic
        && sel->type != SEL_SUBEXPR
        && !(sel->type == SEL_BOOLEAN && sel->u.boolt == BOOL_NOT))
    {
        if (sel->cdata->bEvalMax)
        {
            gmx_ana_index_copy(sel->v.u.g, sel->cdata->gmax, FALSE);
        }
        else
        {
            gmx_ana_index_copy(sel->v.u.g, sel->cdata->gmin, FALSE);
        }
    }

    return 0;
}


/********************************************************************
 * ROOT ITEM INITIALIZATION COMPILER PASS
 ********************************************************************/

/*! \brief
 * Initializes a \ref SEL_ROOT element.
 *
 * \param     root Root element to initialize.
 * \returns Pointer to the selection element that should replace \p root.
 *   Can be \p root itself or NULL if the selection should be removed.
 *
 * Checks whether it is necessary to evaluate anything through the root
 * element, and either clears the evaluation function or initializes the
 * evaluation group.
 *
 * If the function returns NULL, the memory allocated for \p root is 
 * automatically freed.
 */
static t_selelem *
init_root_item(t_selelem *root)
{
    t_selelem   *expr;
    char        *name;

    /* Process subexpressions */
    if (root->child->type == SEL_SUBEXPR)
    {
        if (root->child->refcount == 1)
        {
            /* Free subexpressions that are no longer used */
            _gmx_selelem_free(root);
            return NULL;
        }
        else if (root->child->v.type == POS_VALUE)
        {
            /* Position values only need to be evaluated once, by the root */
        }
        else if (!root->child->cdata->bStaticEval)
        {
            /* Subexpressions with non-static evaluation group should not be
             * evaluated by the root. */
            root->evaluate = NULL;
            if (root->cdata)
            {
                root->cdata->evaluate = NULL;
            }
        }
    }

    /* Set the evaluation group */
    name = root->u.cgrp.name;
    if (root->evaluate)
    {
        expr = root->child;
        if (expr->type == SEL_SUBEXPR)
        {
            gmx_ana_index_copy(&root->u.cgrp, expr->cdata->gmax, TRUE);
        }
        else
        {
            /* expr should evaluate the positions for a selection */
            if (expr->v.u.p->g)
            {
                _gmx_selelem_set_vtype(root, GROUP_VALUE);
                root->flags  |= (SEL_ALLOCVAL | SEL_ALLOCDATA);
                _gmx_selvalue_reserve(&root->v, 1);
                gmx_ana_index_copy(root->v.u.g, expr->v.u.p->g, TRUE);
            }
            gmx_ana_index_set(&root->u.cgrp, -1, NULL, NULL, 0);
        }
    }
    else
    {
        gmx_ana_index_clear(&root->u.cgrp);
    }
    root->u.cgrp.name = name;
    return root;
}

/*! \brief
 * Initializes the evaluation groups for \ref SEL_ROOT items.
 *
 * \param   root First selection in the whole selection chain.
 * \returns The new first element for the chain.
 *
 * The function also removes static subexpressions that are no longer used.
 */
static t_selelem *
process_roots(t_selelem *root)
{
    t_selelem *item;
    t_selelem *next;

    do
    {
        next = root->next;
        root = init_root_item(root);
        if (!root)
        {
            root = next;
        }
    }
    while (root == next);
    item = root;
    while (item && item->next)
    {
        next = item->next->next;
        item->next = init_root_item(item->next);
        if (!item->next)
        {
            item->next = next;
        }
        else
        {
            item = item->next;
        }
    }
    return root;
}


/********************************************************************
 * SUBEXPRESSION OPTIMIZATION PASS
 ********************************************************************/

/*! \brief
 * Optimizes subexpression evaluation.
 *
 * \param     sel Root of the selection subtree to process.
 *
 * Optimizes away some unnecessary evaluation of subexpressions that are only
 * referenced once.
 */
static void
optimize_item_subexpressions(t_selelem *sel)
{
    t_selelem *child;

    /* Call recursively for all children unless the children have already been processed */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            optimize_item_subexpressions(child);
            child = child->next;
        }
    }

    /* Optimize the evaluation of subexpressions that are used only once */
    if (sel->type == SEL_SUBEXPRREF && sel->cdata->bStaticEval && sel->child->refcount == 2)
    {
        /* The evaluation functions are not always needed, and they do some
         * extra work, but it should be neglible compared to other factors
         * in the evaluation, so set them always for simplicity. */
        sel->evaluate = &_gmx_sel_evaluate_subexprref_pass;
        if (sel->cdata)
        {
            sel->cdata->evaluate = sel->evaluate;
        }
        /* Replace the value of the child */
        _gmx_selelem_free_values(sel->child);
        sel->child->flags            &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
        sel->child->evaluate = &_gmx_sel_evaluate_subexpr_pass;
        if (sel->child->cdata)
        {
            sel->child->cdata->evaluate = sel->child->evaluate;
        }
        /* Replace the value of the grandchild */
        _gmx_selelem_free_values(sel->child->child);
        sel->child->child->flags     &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        _gmx_selvalue_setstore(&sel->child->child->v, sel->v.u.ptr);
    }
}


/********************************************************************
 * COM CALCULATION COMPILER PASS
 ********************************************************************/

/*! \brief
 * Initializes COM/COG calculation for method expressions that require it.
 *
 * \param     sel    Selection subtree to process.
 * \param[in,out] pcc   Position calculation collection to use.
 * \param[in] type   Default position calculation type.
 * \param[in] flags  Flags for default position calculation.
 * \returns   0 on success, a non-zero error code on error.
 *
 * Searches recursively through the selection tree for dynamic
 * \ref SEL_EXPRESSION elements that define the \c gmx_ana_selmethod_t::pupdate
 * function.
 * For each such element found, position calculation is initialized
 * for the maximal evaluation group.
 * The type of the calculation is determined by \p type and \p flags.
 * No calculation is initialized if \p type equals \ref POS_ATOM and
 * the method also defines the \c gmx_ana_selmethod_t::update method.
 */
static int
init_item_comg(t_selelem *sel, gmx_ana_poscalc_coll_t *pcc,
               e_poscalc_t type, int flags)
{
    t_selelem *child;
    int        rc;

    /* Initialize COM calculation for dynamic selections now that we know the maximal evaluation group */
    if (sel->type == SEL_EXPRESSION && sel->u.expr.method
        && sel->u.expr.method->pupdate)
    {
        if (!sel->u.expr.method->update || type != POS_ATOM)
        {
            /* Create a default calculation if one does not yet exist */
            int cflags;
            cflags = 0;
            if (!sel->cdata->bStaticEval)
            {
                cflags |= POS_DYNAMIC;
            }
            if (!sel->u.expr.pc)
            {
                cflags |= flags;
                rc = gmx_ana_poscalc_create(&sel->u.expr.pc, pcc, type, cflags);
                if (rc != 0)
                {
                    return rc;
                }
            }
            else
            {
                gmx_ana_poscalc_set_flags(sel->u.expr.pc, cflags);
            }
            gmx_ana_poscalc_set_maxindex(sel->u.expr.pc, sel->cdata->gmax);
            snew(sel->u.expr.pos, 1);
            gmx_ana_poscalc_init_pos(sel->u.expr.pc, sel->u.expr.pos);
        }
    }

    /* Call recursively for all children unless the children have already been processed */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            rc = init_item_comg(child, pcc, type, flags);
            if (rc != 0)
            {
                return rc;
            }
            child = child->next;
        }
    }
    return 0;
}


/********************************************************************
 * FREE COMPILER DATA PASS
 ********************************************************************/

/*! \brief
 * Frees the allocated compiler data recursively.
 *
 * \param     sel Root of the selection subtree to process.
 *
 * Frees the data allocated for the compilation process.
 */
static void
free_item_compilerdata(t_selelem *sel)
{
    t_selelem *child;

    /* Free compilation data */
    _gmx_selelem_free_compiler_data(sel);

    /* Call recursively for all children unless the children have already been processed */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            free_item_compilerdata(child);
            child = child->next;
        }
    }
}


/********************************************************************
 * INFORMATION UPDATE
 ********************************************************************/

/*! \brief
 * Updates the information about the selection.
 *
 * \param[in]     top   Topology information.
 * \param[in]     ngrps Number of elements in the \p sel array.
 * \param[in,out] sel   Array of selections to update.
 * \param[in]     bMaskOnly TRUE if the positions will always be calculated
 *   for all atoms, i.e., the masses/charges do not change.
 *
 * Initializes total masses and charges.
 */
static void
update_info(t_topology *top, int ngrps, gmx_ana_selection_t *sel[],
            bool bMaskOnly)
{
    int   g, b, i;

    for (g = 0; g < ngrps; ++g)
    {
        sel[g]->g = sel[g]->p.g;
        snew(sel[g]->orgm, sel[g]->p.nr);
        snew(sel[g]->orgq, sel[g]->p.nr);
        for (b = 0; b < sel[g]->p.nr; ++b)
        {
            sel[g]->orgq[b] = 0;
            if (top)
            {
                sel[g]->orgm[b] = 0;
                for (i = sel[g]->p.m.mapb.index[b]; i < sel[g]->p.m.mapb.index[b+1]; ++i)
                {
                    sel[g]->orgm[b] += top->atoms.atom[sel[g]->g->index[i]].m;
                    sel[g]->orgq[b] += top->atoms.atom[sel[g]->g->index[i]].q;
                }
            }
            else
            {
                sel[g]->orgm[b] = 1;
            }
        }
        if (sel[g]->bDynamic && !bMaskOnly)
        {
            snew(sel[g]->m, sel[g]->p.nr);
            snew(sel[g]->q, sel[g]->p.nr);
            for (b = 0; b < sel[g]->p.nr; ++b)
            {
                sel[g]->m[b] = sel[g]->orgm[b];
                sel[g]->q[b] = sel[g]->orgq[b];
            }
        }
        else
        {
            sel[g]->m = sel[g]->orgm;
            sel[g]->q = sel[g]->orgq;
        }
    }
}


/********************************************************************
 * MAIN COMPILATION FUNCTION
 ********************************************************************/

/*!
 * \param[in,out] sc Selection collection to be compiled.
 * \returns       0 on successful compilation, a non-zero error code on error.
 *
 * Before compilation, the selection collection should have been initialized
 * with gmx_ana_selcollection_parse_*().
 * The compiled selection collection can be passed to
 * gmx_ana_selcollection_evaluate() to evaluate the selection for a frame.
 * If an error occurs, \p sc is cleared.
 *
 * The covered fraction information in \p sc is initialized to
 * \ref CFRAC_NONE.
 */
int
gmx_ana_selcollection_compile(gmx_ana_selcollection_t *sc)
{
    gmx_sel_evaluate_t  evaldata;
    t_selelem   *item;
    e_poscalc_t  post;
    int          flags;
    int          rc;

    _gmx_sel_evaluate_init(&evaldata, &sc->gall, sc->top, NULL, NULL);

    /* Clear the symbol table because it is not possible to parse anything
     * after compilation, and variable references in the symbol table can
     * also mess up the compilation and/or become invalid.
     */
    _gmx_selcollection_clear_symtab(sc);

    /* Extract subexpressions into separate roots */
    sc->root = extract_subexpressions(sc->root, &sc->gall);

    /* Initialize the evaluation callbacks and process the tree structure
     * to conform to the expectations of the callback functions. */
    /* Also, initialize and allocate the compiler data structure */
    item = sc->root;
    while (item)
    {
        /* Process boolean expressions */
        optimize_boolean_expressions(item);
        reorder_boolean_static_children(item);
        /* Initialize evaluation */
        if (!init_item_evaluation(item))
        {
            /* FIXME: Clean up the collection */
            return -1;
        }
        /* Initialize the compiler data */
        init_item_compilerdata(item);
        init_item_staticeval(item);
        item = item->next;
    }
    /* Initialize the evaluation index groups */
    initialize_evalgrps(sc);

    /* Evaluate all static parts of the selection and analyze the tree
     * to allocate enough memory to store the value of each dynamic subtree. */
    item = sc->root;
    while (item)
    {
        if (item->child->type == SEL_SUBEXPR && item->child->refcount > 2)
        {
            mark_subexpr_dynamic(item->child, TRUE);
        }
        set_evaluation_function(item, &analyze_static);
        rc = item->evaluate(&evaldata, item, NULL);
        if (rc != 0)
        {
            /* FIXME: Clean up the collection */
            return rc;
        }
        item = item->next;
    }

    /* Do a second pass to evaluate static parts of common subexpressions */
    /* Note that the refcount check skips constant subexpressions completely
     * since they were already evaluated by analyze_static(). */
    item = sc->root;
    while (item)
    {
        if (item->child->type == SEL_SUBEXPR && item->child->refcount > 2)
        {
            mark_subexpr_dynamic(item->child, FALSE);
            item->child->u.cgrp.isize = 0;
            if (item->child->v.type == GROUP_VALUE)
            {
                item->child->child->v.u.g->isize = 0;
            }
            set_evaluation_function(item, &analyze_static2);
            rc = item->evaluate(&evaldata, item->child, item->child->cdata->gmax);
            if (rc != 0)
            {
                /* FIXME: Clean up the collection */
                return rc;
            }
        }
        item = item->next;
    }

    /* Initialize evaluation groups and remove unused subexpressions. */
    sc->root = process_roots(sc->root);

    /* Initialize position calculations for methods, perform some final
     * optimization and free the memory allocated for the compilation. */
    /* By default, use whole residues/molecules. */
    flags = POS_COMPLWHOLE;
    rc = gmx_ana_poscalc_type_from_enum(sc->rpost, &post, &flags);
    if (rc != 0)
    {
        gmx_bug("invalid default reference position type");
        /* FIXME: Clean up the collection */
        return rc;
    }
    item = sc->root;
    while (item)
    {
        optimize_item_subexpressions(item);
        rc = init_item_comg(item, sc->pcc, post, flags);
        if (rc != 0)
        {
            /* FIXME: Clean up the collection */
            return rc;
        }
        free_item_compilerdata(item);
        item = item->next;
    }

    /* Finish up by updating some information */
    update_info(sc->top, sc->nr, sc->sel, sc->bMaskOnly);

    return 0;
}
