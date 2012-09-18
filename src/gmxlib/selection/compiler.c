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
 * At least, the main compilation function leaves the selection collection in
 * a bad state if an error occurs.
 *
 * \todo
 * The memory usage could still be optimized.
 * Use of memory pooling could still be extended, and a lot of redundant
 * gmin/gmax data could be eliminated for complex arithmetic expressions.
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
 *       e.g., "A and B and C". The subexpressions in gmx_boolean expression are
 *       reordered such that static expressions come first without otherwise
 *       altering the relative order of the expressions.
 *    -# The \c t_selelem::evaluate field is set to the correct evaluation
 *       function from evaluate.h.
 *    -# The compiler data structure is allocated for each element, and
 *       the fields are initialized, with the exception of the contents of
 *       \c gmax and \c gmin fields.  In reality, several passes are made
 *       to completely initialize the structure, because some flags are set
 *       recursively based on which elements refer to an element, and these
 *       flags need to be set to initialize other fields.
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
 *     special care needs to be taken with gmx_boolean expressions because they
 *     are short-circuiting. This is done through the
 *     \c SEL_CDATA_EVALMAX flag, which makes dynamic child expressions
 *     of \c BOOL_OR expressions evaluate to empty groups, while subexpressions
 *     of \c BOOL_AND are evaluated to largest possible groups.
 *     Memory is also allocated to store the results of the evaluation.
 *     For each element, analyze_static() calls the actual evaluation function
 *     after the element has been properly initialized.
 *  -# Another evaluation pass is done over subexpressions with more than
 *     one reference to them. These cannot be completely processed during the
 *     first pass, because it is not known whether later references require
 *     additional evaluation of static expressions.
 *  -# Unused subexpressions are removed. For efficiency reasons (and to avoid
 *     some checks), this is actually done several times already earlier in
 *     the compilation process.
 *  -# Most of the processing is now done, and the next pass simply sets the
 *     evaluation group of root elements to the largest selection as determined
 *     in pass 3.  For root elements of subexpressions that should not be
 *     evaluated before they are referred to, the evaluation group/function is
 *     cleared.  At the same time, position calculation data is initialized for
 *     for selection method elements that require it.  Compiler data is also
 *     freed as it is no longer needed.
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
 *     possible group that the selection can evaluate to (only for dynamic
 *     selections), and \ref t_selelem::cgrp has been set to use a NULL group
 *     for evaluation.
 *  -# A subexpression that appears in one or more selections.
 *     Each selection that gives a value for a method parameter is a
 *     potential subexpression, as is any variable value.
 *     Only subexpressions that require evaluation for each frame are left
 *     after the selection is compiled.
 *     Each subexpression appears in the chain before any references to it.
 *     For these elements, \c t_selelem::cgrp has been set to the group
 *     that should be used to evaluate the subexpression.
 *     If \c t_selelem::cgrp is empty, the total evaluation group is not known
 *     in advance or it is more efficient to evaluate the subexpression only
 *     when it is referenced.  If this is the case, \c t_selelem::evaluate is
 *     also NULL.
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
 * evaluated during the current frame.  This is only done if full subexpression
 * evaluation by _gmx_sel_evaluate_subexpr() is needed; the other evaluation
 * functions do not require this memory.
 *
 * \ref SEL_SUBEXPRREF elements are used to describe references to
 * subexpressions. They have always a single child, which is the
 * \ref SEL_SUBEXPR element being referenced.
 *
 * If a subexpression is used only once, the evaluation has been optimized by
 * setting the child of the \ref SEL_SUBEXPR element to evaluate the value of
 * \ref SEL_SUBEXPRREF directly (in the case of memory pooling, this is managed
 * by the evaluation functions).  In such cases, the evaluation routines for the
 * \ref SEL_SUBEXPRREF and \ref SEL_SUBEXPR elements only propagate some status
 * information, but do not unnecessarily copy the values.
 *
 *
 * \subsection selcompiler_tree_gmx_bool Boolean elements
 *
 * \ref SEL_BOOLEAN elements have been merged such that one element
 * may carry out evaluation of more than one operation of the same type.
 * The static parts of the expressions have been evaluated, and are placed
 * in the first child. These are followed by the dynamic expressions, in the
 * order provided by the user.
 *
 *
 * \subsection selcompiler_tree_arith Arithmetic elements
 *
 * Constant and static expressions in \ref SEL_ARITHMETIC elements have been
 * calculated.
 * Currently, no other processing is done.
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
#include "mempool.h"
#include "selcollection.h"
#include "selelem.h"

/*! \internal \brief
 * Compiler flags.
 */
enum
{
    /*! \brief
     * Whether a subexpression needs to evaluated for all atoms.
     *
     * This flag is set for \ref SEL_SUBEXPR elements that are used to
     * evaluate non-atom-valued selection method parameters, as well as
     * those that are used directly as values of selections.
     */
    SEL_CDATA_FULLEVAL    =  1,
    /*! \brief
     * Whether the whole subexpression should be treated as static.
     *
     * This flag is always FALSE if \ref SEL_DYNAMIC is set for the element,
     * but it is also FALSE for static elements within common subexpressions.
     */
    SEL_CDATA_STATIC      =  2,
    /** Whether the subexpression will always be evaluated in the same group. */
    SEL_CDATA_STATICEVAL  =  4,
    /** Whether the compiler evaluation routine should return the maximal selection. */
    SEL_CDATA_EVALMAX     =  8,
    /** Whether memory has been allocated for \p gmin and \p gmax. */
    SEL_CDATA_MINMAXALLOC = 16,
    /** Whether subexpressions use simple pass evaluation functions. */
    SEL_CDATA_SIMPLESUBEXPR = 32,
    /** Whether this expressions is a part of a common subexpression. */
    SEL_CDATA_COMMONSUBEXPR = 64 
};

/*! \internal \brief
 * Internal data structure used by the compiler.
 */
typedef struct t_compiler_data
{
    /** The real evaluation method. */
    sel_evalfunc     evaluate;
    /** Flags for specifying how to treat this element during compilation. */
    int              flags;
    /** Smallest selection that can be selected by the subexpression. */
    gmx_ana_index_t *gmin;
    /** Largest selection that can be selected by the subexpression. */
    gmx_ana_index_t *gmax;
} t_compiler_data;


/********************************************************************
 * COMPILER UTILITY FUNCTIONS
 ********************************************************************/

static void
print_group_info(FILE *fp, const char *name, t_selelem *sel, gmx_ana_index_t *g)
{
    fprintf(fp, " %s=", name);
    if (!g)
    {
        fprintf(fp, "(null)");
    }
    else if (sel->cdata->flags & SEL_CDATA_MINMAXALLOC)
    {
        fprintf(fp, "(%d atoms, %p)", g->isize, (void*)g);
    }
    else if (sel->v.type == GROUP_VALUE && g == sel->v.u.g)
    {
        fprintf(fp, "(static, %p)", (void*)g);
    }
    else
    {
        fprintf(fp, "%p", (void*)g);
    }
}

/*!
 * \param[in] fp      File handle to receive the output.
 * \param[in] sel     Selection element to print.
 * \param[in] level   Indentation level, starting from zero.
 */
void
_gmx_selelem_print_compiler_info(FILE *fp, t_selelem *sel, int level)
{
    if (!sel->cdata)
    {
        return;
    }
    fprintf(fp, "%*c cdata: flg=", level*2+1, ' ');
    if (sel->cdata->flags & SEL_CDATA_FULLEVAL)
    {
        fprintf(fp, "F");
    }
    if (!(sel->cdata->flags & SEL_CDATA_STATIC))
    {
        fprintf(fp, "D");
    }
    if (sel->cdata->flags & SEL_CDATA_STATICEVAL)
    {
        fprintf(fp, "S");
    }
    if (sel->cdata->flags & SEL_CDATA_EVALMAX)
    {
        fprintf(fp, "M");
    }
    if (sel->cdata->flags & SEL_CDATA_MINMAXALLOC)
    {
        fprintf(fp, "A");
    }
    if (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
    {
        fprintf(fp, "Ss");
    }
    if (sel->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
    {
        fprintf(fp, "Sc");
    }
    if (!sel->cdata->flags)
    {
        fprintf(fp, "0");
    }
    fprintf(fp, " eval=");
    _gmx_sel_print_evalfunc_name(fp, sel->cdata->evaluate);
    print_group_info(fp, "gmin", sel, sel->cdata->gmin);
    print_group_info(fp, "gmax", sel, sel->cdata->gmax);
    fprintf(fp, "\n");
}

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
        if (sel->cdata->flags & SEL_CDATA_MINMAXALLOC)
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
 *
 * Allocation of POS_VALUE selection elements is a special case, and is
 * handled by alloc_selection_pos_data().
 */
static gmx_bool
alloc_selection_data(t_selelem *sel, int isize, gmx_bool bChildEval)
{
    int        nalloc;

    if (sel->mempool)
    {
        return TRUE;
    }
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
        nalloc = child->v.nr;
    }
    /* Allocate memory for sel->v.u if needed */
    if (sel->flags & SEL_ALLOCVAL)
    {
        _gmx_selvalue_reserve(&sel->v, nalloc);
    }
    /* Reserve memory inside group structure if SEL_ALLOCDATA is set. */
    if ((sel->flags & SEL_ALLOCDATA) && sel->v.type == GROUP_VALUE)
    {
        gmx_ana_index_reserve(sel->v.u.g, isize);
    }
    return TRUE;
}

/*! \brief
 * Allocates memory for storing the evaluated value of a selection element.
 *
 * \param     sel   Selection element to initialize.
 *
 * Allocation of POS_VALUE selection elements is a special case, and is
 * handled by this function instead of by alloc_selection_data().
 */
static void
alloc_selection_pos_data(t_selelem *sel)
{
    t_selelem *child;
    int        nalloc, isize;

    if (sel->mempool)
    {
        return;
    }

    child = sel;
    if (sel->type == SEL_SUBEXPRREF)
    {
        child = sel->child->child;
    }
    nalloc = child->v.u.p->nr;
    isize = child->v.u.p->m.b.nra;

    /* For positions, we want to allocate just a single structure
     * for nalloc positions. */
    if (sel->flags & SEL_ALLOCVAL)
    {
        _gmx_selvalue_reserve(&sel->v, 1);
    }
    sel->v.nr = 1;
    /* Reserve memory inside position structure if SEL_ALLOCDATA is set. */
    if (sel->flags & SEL_ALLOCDATA)
    {
        gmx_ana_pos_reserve(sel->v.u.p, nalloc, isize);
    }
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
 * SUBEXPRESSION PROCESSING
 ********************************************************************/

/*! \brief
 * Reverses the chain of selection elements starting at \p root.
 *
 * \param   root First selection in the whole selection chain.
 * \returns The new first element for the chain.
 */
static t_selelem *
reverse_selelem_chain(t_selelem *root)
{
    t_selelem *item;
    t_selelem *prev;
    t_selelem *next;

    prev = NULL;
    item = root;
    while (item)
    {
        next = item->next;
        item->next = prev;
        prev = item;
        item = next;
    }
    return prev;
}

/*! \brief
 * Removes subexpressions that don't have any references.
 *
 * \param     root First selection in the whole selection chain.
 * \returns   The new first element for the chain.
 *
 * The elements are processed in reverse order to correctly detect
 * subexpressions only referred to by other subexpressions.
 */
static t_selelem *
remove_unused_subexpressions(t_selelem *root)
{
    t_selelem *item;
    t_selelem *prev;
    t_selelem *next;

    root = reverse_selelem_chain(root);
    while (root->child->type == SEL_SUBEXPR && root->child->refcount == 1)
    {
        next = root->next;
        _gmx_selelem_free(root);
        root = next;
    }
    prev = root;
    item = root->next;
    while (item)
    {
        next = item->next;
        if (item->child->type == SEL_SUBEXPR && item->child->refcount == 1)
        {
            prev->next = next;
            _gmx_selelem_free(item);
        }
        else
        {
            prev = item;
        }
        item = next;
    }
    return reverse_selelem_chain(root);
}

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
 * \param   subexprn Pointer to a subexpression counter.
 * \returns Pointer to a chain of subselections, or NULL if none were found.
 *
 * This function finds recursively all \ref SEL_SUBEXPRREF elements below
 * the given root element and ensures that their children are within
 * \ref SEL_SUBEXPR elements. It also creates a chain of \ref SEL_ROOT elements
 * that contain the subexpression as their children and returns the first
 * of these root elements.
 */
static t_selelem *
extract_item_subselections(t_selelem *sel, int *subexprn)
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
            root = subexpr = extract_item_subselections(child, subexprn);
        }
        else
        {
            subexpr->next = extract_item_subselections(child, subexprn);
        }
        while (subexpr && subexpr->next)
        {
            subexpr = subexpr->next;
        }
        /* The latter check excludes variable references.
         * It also excludes subexpression elements that have already been
         * processed, because they are given a name when they are first
         * encountered.
         * TODO: There should be a more robust mechanism (probably a dedicated
         * flag) for detecting parser-generated subexpressions than relying on
         * a NULL name field. */
        if (child->type == SEL_SUBEXPRREF && (child->child->type != SEL_SUBEXPR
                                              || child->child->name == NULL))
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
            /* Create the subexpression element and/or
             * move the actual subexpression under the created element. */
            if (child->child->type != SEL_SUBEXPR)
            {
                subexpr->child = _gmx_selelem_create(SEL_SUBEXPR);
                _gmx_selelem_set_vtype(subexpr->child, child->v.type);
                subexpr->child->child = child->child;
                child->child          = subexpr->child;
            }
            else
            {
                subexpr->child = child->child;
            }
            create_subexpression_name(subexpr->child, ++*subexprn);
            subexpr->child->refcount++;
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
extract_subexpressions(t_selelem *sel)
{
    t_selelem   *root, *item, *next;
    int          subexprn;

    subexprn = 0;
    root = NULL;
    next = sel;
    while (next)
    {
        item = extract_item_subselections(next, &subexprn);
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
 * BOOLEAN OPERATION REORDERING
 ********************************************************************/

/*! \brief
 * Removes redundant gmx_boolean selection elements.
 *
 * \param  sel Root of the selection subtree to optimize.
 *
 * This function merges similar gmx_boolean operations (e.g., (A or B) or C becomes
 * a single OR operation with three operands).
 */
static void
optimize_gmx_boolean_expressions(t_selelem *sel)
{
    t_selelem *child, *prev;

    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        prev  = NULL;
        child = sel->child;
        while (child)
        {
            optimize_gmx_boolean_expressions(child);
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
 * Reorders children of gmx_boolean expressions such that static selections
 * come first.
 *
 * \param  sel Root of the selection subtree to reorder.
 *
 * The relative order of static expressions does not change.
 * The same is true for the dynamic expressions.
 */
static void
reorder_gmx_boolean_static_children(t_selelem *sel)
{
    t_selelem *child, *prev, *next;

    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            reorder_gmx_boolean_static_children(child);
            child = child->next;
        }
    }

    /* Reorder gmx_boolean expressions such that static selections come first */
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
 * ARITHMETIC EXPRESSION PROCESSING
 ********************************************************************/

/*! \brief
 * Processes arithmetic expressions to simplify and speed up evaluation.
 *
 * \param  sel Root of the selection subtree to process.
 *
 * Currently, this function only converts integer constants to reals
 * within arithmetic expressions.
 */
static gmx_bool
optimize_arithmetic_expressions(t_selelem *sel)
{
    t_selelem  *child;
    gmx_bool        bOk;

    /* Do recursively for children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        child = sel->child;
        while (child)
        {
            bOk = optimize_arithmetic_expressions(child);
            if (!bOk)
            {
                return bOk;
            }
            child = child->next;
        }
    }

    if (sel->type != SEL_ARITHMETIC)
    {
        return TRUE;
    }

    /* Convert integer constants to reals. */
    child = sel->child;
    while (child)
    {
        if (child->v.type == INT_VALUE)
        {
            real  *r;

            if (child->type != SEL_CONST)
            {
                gmx_impl("Non-constant integer expressions not implemented in arithmetic evaluation");
                return FALSE;
            }
            snew(r, 1);
            r[0] = child->v.u.i[0];
            sfree(child->v.u.i);
            child->v.u.r = r;
            child->v.type = REAL_VALUE;
        }
        else if (child->v.type != REAL_VALUE)
        {
            gmx_bug("Internal error");
            return FALSE;
        }
        child = child->next;
    }
    return TRUE;
}


/********************************************************************
 * EVALUATION PREPARATION COMPILER
 ********************************************************************/

/*! \brief
 * Sets the evaluation functions for the selection (sub)tree.
 *
 * \param[in,out] sel Root of the selection subtree to process.
 * \returns       TRUE on success, FALSE if any subexpression fails.
 *
 * This function sets the evaluation function (\c t_selelem::evaluate)
 * for the selection elements.
 */
static gmx_bool
init_item_evalfunc(t_selelem *sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem *child;

        child = sel->child;
        while (child)
        {
            if (!init_item_evalfunc(child))
            {
                return FALSE;
            }
            child = child->next;
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
            if (!(sel->flags & SEL_DYNAMIC) && sel->u.expr.method
                && sel->u.expr.method->init_frame)
            {
                sel->flags |= SEL_INITFRAME;
            }
            sel->evaluate = &_gmx_sel_evaluate_method;
            break;

        case SEL_ARITHMETIC:
            sel->evaluate = &_gmx_sel_evaluate_arithmetic;
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
            sel->evaluate = (sel->refcount == 2
                             ? &_gmx_sel_evaluate_subexpr_simple
                             : &_gmx_sel_evaluate_subexpr);
            break;

        case SEL_SUBEXPRREF:
            sel->name     = sel->child->name;
            sel->evaluate = (sel->child->refcount == 2
                             ? &_gmx_sel_evaluate_subexprref_simple
                             : &_gmx_sel_evaluate_subexprref);
            break;
    }

    return TRUE;
}

/*! \brief
 * Sets the memory pool for selection elements that can use it.
 *
 * \param     sel      Root of the selection subtree to process.
 * \param[in] mempool  Memory pool to use.
 */
static void
setup_memory_pooling(t_selelem *sel, gmx_sel_mempool_t *mempool)
{
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem         *child;

        child = sel->child;
        while (child)
        {
            if ((sel->type == SEL_BOOLEAN && (child->flags & SEL_DYNAMIC))
                || (sel->type == SEL_ARITHMETIC && child->type != SEL_CONST
                    && !(child->flags & SEL_SINGLEVAL))
                || (sel->type == SEL_SUBEXPR && sel->refcount > 2))
            {
                child->mempool = mempool;
                if (child->type == SEL_SUBEXPRREF
                    && child->child->refcount == 2)
                {
                    child->child->child->mempool = mempool;
                }
            }
            setup_memory_pooling(child, mempool);
            child = child->next;
        }
    }
}

/*! \brief
 * Prepares the selection (sub)tree for evaluation.
 *
 * \param[in,out] sel Root of the selection subtree to prepare.
 *
 * It also allocates memory for the \p sel->v.u.g or \p sel->v.u.p
 * structure if required.
 */
static void
init_item_evaloutput(t_selelem *sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem *child;

        child = sel->child;
        while (child)
        {
            init_item_evaloutput(child);
            child = child->next;
        }
    }

    if (sel->type == SEL_SUBEXPR && sel->refcount == 2)
    {
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
        }
    }
    else if (sel->type == SEL_SUBEXPR
             && (sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        sel->evaluate = &_gmx_sel_evaluate_subexpr_staticeval;
        sel->cdata->evaluate = sel->evaluate;
        sel->child->mempool = NULL;
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
        }
    }
    else if (sel->type == SEL_SUBEXPRREF && sel->child->refcount == 2)
    {
        if (sel->v.u.ptr)
        {
            _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
            _gmx_selelem_free_values(sel->child->child);
            sel->child->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
            sel->child->child->flags |= (sel->flags & SEL_ALLOCDATA);
            _gmx_selvalue_setstore(&sel->child->child->v, sel->v.u.ptr);
        }
        else if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_setstore(&sel->v, sel->child->child->v.u.ptr);
        }
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
    }

    /* Make sure that the group/position structure is allocated. */
    if (!sel->v.u.ptr && (sel->flags & SEL_ALLOCVAL))
    {
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_reserve(&sel->v, 1);
            sel->v.nr = 1;
        }
    }
}


/********************************************************************
 * COMPILER DATA INITIALIZATION
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
    sel->cdata->flags = SEL_CDATA_STATICEVAL;
    if (!(sel->flags & SEL_DYNAMIC))
    {
        sel->cdata->flags |= SEL_CDATA_STATIC;
    }
    if (sel->type == SEL_SUBEXPR)
    {
        sel->cdata->flags |= SEL_CDATA_EVALMAX;
    }
    /* Set the full evaluation flag for subexpressions that require it;
     * the subexpression has already been initialized, so we can simply
     * access its compilation flags.*/
    if (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER)
    {
        child = sel->child;
        while (child)
        {
            if (!(child->flags & SEL_ATOMVAL) && child->child)
            {
                child->child->cdata->flags |= SEL_CDATA_FULLEVAL;
            }
            child = child->next;
        }
    }
    else if (sel->type == SEL_ROOT && sel->child->type == SEL_SUBEXPRREF)
    {
        sel->child->child->cdata->flags |= SEL_CDATA_FULLEVAL;
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
        gmx_bool  bEvalMax;

        bEvalMax = (sel->u.boolt == BOOL_AND);
        child = sel->child;
        while (child)
        {
            if (bEvalMax)
            {
                child->cdata->flags |= SEL_CDATA_EVALMAX;
            }
            else if (child->type == SEL_BOOLEAN && child->u.boolt == BOOL_NOT)
            {
                child->child->cdata->flags |= SEL_CDATA_EVALMAX;
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
            child->cdata->flags |= SEL_CDATA_EVALMAX;
            child = child->next;
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
 * reorder_gmx_boolean_static_children() should have been called.
 */
static void
init_item_staticeval(t_selelem *sel)
{
    t_selelem   *child;

    /* Subexpressions with full evaluation should always have bStaticEval,
     * so don't do anything if a reference to them is encountered. */
    if (sel->type == SEL_SUBEXPRREF
        && (sel->child->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        return;
    }

    /* Propagate the bStaticEval flag to children if it is not set */
    if (!(sel->cdata->flags & SEL_CDATA_STATICEVAL))
    {
        child = sel->child;
        while (child)
        {
            if ((sel->type != SEL_EXPRESSION && sel->type != SEL_MODIFIER)
                || (child->flags & SEL_ATOMVAL))
            {
                if (child->cdata->flags & SEL_CDATA_STATICEVAL)
                {
                    child->cdata->flags &= ~SEL_CDATA_STATICEVAL;
                    init_item_staticeval(child);
                }
            }
            child = child->next;
        }
    }
    else /* bStaticEval is set */
    {
        /* For gmx_boolean expressions, any expression after the first dynamic
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
                child->cdata->flags &= ~SEL_CDATA_STATICEVAL;
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

/*! \brief
 * Initializes compiler flags for subexpressions.
 *
 * \param sel Root of the selection subtree to process.
 */
static void
init_item_subexpr_flags(t_selelem *sel)
{
    if (sel->type == SEL_SUBEXPR)
    {
        if (sel->refcount == 2)
        {
            sel->cdata->flags |= SEL_CDATA_SIMPLESUBEXPR;
        }
        else if (!(sel->cdata->flags & SEL_CDATA_FULLEVAL))
        {
            sel->cdata->flags |= SEL_CDATA_COMMONSUBEXPR;
        }
    }
    else if (sel->type == SEL_SUBEXPRREF && sel->child->refcount == 2)
    {
        sel->cdata->flags |= SEL_CDATA_SIMPLESUBEXPR;
    }

    /* Process children, but only follow subexpression references if the
     * common subexpression flag needs to be propagated. */
    if (sel->type != SEL_SUBEXPRREF
        || ((sel->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
            && sel->child->refcount > 2))
    {
        t_selelem *child = sel->child;

        while (child)
        {
            if (!(child->cdata->flags & SEL_CDATA_COMMONSUBEXPR))
            {
                if (sel->type != SEL_EXPRESSION || (child->flags & SEL_ATOMVAL))
                {
                    child->cdata->flags |=
                        (sel->cdata->flags & SEL_CDATA_COMMONSUBEXPR);
                }
                init_item_subexpr_flags(child);
            }
            child = child->next;
        }
    }
}

/*! \brief
 * Initializes the gmin and gmax fields of the compiler data structure.
 *
 * \param sel Root of the selection subtree to process.
 */
static void
init_item_minmax_groups(t_selelem *sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem *child;

        child = sel->child;
        while (child)
        {
            init_item_minmax_groups(child);
            child = child->next;
        }
    }

    /* Initialize the minimum and maximum evaluation groups. */
    if (sel->type != SEL_ROOT && sel->v.type != NO_VALUE)
    {
        if (sel->v.type == GROUP_VALUE
            && (sel->cdata->flags & SEL_CDATA_STATIC))
        {
            sel->cdata->gmin = sel->v.u.g;
            sel->cdata->gmax = sel->v.u.g;
        }
        else if (sel->type == SEL_SUBEXPR
                 && ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                     || (sel->cdata->flags & SEL_CDATA_FULLEVAL)))
        {
            sel->cdata->gmin = sel->child->cdata->gmin;
            sel->cdata->gmax = sel->child->cdata->gmax;
        }
        else
        {
            sel->cdata->flags |= SEL_CDATA_MINMAXALLOC;
            snew(sel->cdata->gmin, 1);
            snew(sel->cdata->gmax, 1);
        }
    }
}


/********************************************************************
 * EVALUATION GROUP INITIALIZATION
 ********************************************************************/

/*! \brief
 * Initializes evaluation groups for root items.
 *
 * \param[in,out] sc   Selection collection data.
 *
 * The evaluation group of each \ref SEL_ROOT element corresponding to a
 * selection in \p sc is set to NULL.  The evaluation grop for \ref SEL_ROOT
 * elements corresponding to subexpressions that need full evaluation is set
 * to \c sc->gall.
 */
static void
initialize_evalgrps(gmx_ana_selcollection_t *sc)
{
    t_selelem   *root;

    root = sc->root;
    while (root)
    {
        if (root->child->type != SEL_SUBEXPR
            || (root->child->v.type != GROUP_VALUE && !(root->flags & SEL_ATOMVAL)))
        {
            gmx_ana_index_set(&root->u.cgrp, -1, 0, root->u.cgrp.name, 0);
        }
        else if (root->child->cdata->flags & SEL_CDATA_FULLEVAL)
        {
            gmx_ana_index_set(&root->u.cgrp, sc->gall.isize, sc->gall.index,
                              root->u.cgrp.name, 0);
        }
        root = root->next;
    }
}


/********************************************************************
 * STATIC ANALYSIS
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
mark_subexpr_dynamic(t_selelem *sel, gmx_bool bDynamic)
{
    t_selelem *child;

    if (!bDynamic && !(sel->flags & SEL_DYNAMIC))
    {
        sel->cdata->flags |= SEL_CDATA_STATIC;
    }
    else
    {
        sel->cdata->flags &= ~SEL_CDATA_STATIC;
    }
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
 * Frees memory for subexpressions that are no longer needed.
 *
 * \param     sel      Selection subtree to check.
 *
 * Checks whether the subtree rooted at \p sel refers to any \ref SEL_SUBEXPR
 * elements that are not referred to by anything else except their own root
 * element. If such elements are found, all memory allocated for them is freed
 * except the actual element. The element is left because otherwise a dangling
 * pointer would be left at the root element, which is not traversed by this
 * function. Later compilation passes remove the stub elements.
 */
static void
release_subexpr_memory(t_selelem *sel)
{
    if (sel->type == SEL_SUBEXPR)
    {
        if (sel->refcount == 2)
        {
            release_subexpr_memory(sel->child);
            sel->name = NULL;
            _gmx_selelem_free_chain(sel->child);
            _gmx_selelem_free_values(sel);
            _gmx_selelem_free_exprdata(sel);
            _gmx_selelem_free_compiler_data(sel);
            sel->child = NULL;
        }
    }
    else
    {
        t_selelem *child;

        child = sel->child;
        while (child)
        {
            release_subexpr_memory(child);
            child = child->next;
        }
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
    /* If this is a subexpression reference and the data is stored in the
     * child, we transfer data ownership before doing anything else. */
    if (sel->type == SEL_SUBEXPRREF
        && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
    {
        if (sel->child->child->flags & SEL_ALLOCDATA)
        {
            sel->flags |= SEL_ALLOCDATA;
            sel->child->child->flags &= ~SEL_ALLOCDATA;
        }
        if (sel->child->child->flags & SEL_ALLOCVAL)
        {
            sel->flags |= SEL_ALLOCVAL;
            sel->v.nalloc = sel->child->child->v.nalloc;
            sel->child->child->flags &= ~SEL_ALLOCVAL;
            sel->child->child->v.nalloc = -1;
        }
    }
    /* When we reach here for parameter elements, the value is already
     * stored in the parent element, so make sure that it is not freed
     * through this element. */
    if (sel->type == SEL_SUBEXPRREF && sel->u.param)
    {
        sel->u.param->val.nalloc = sel->v.nalloc;
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        sel->v.nalloc = -1;
    }
    /* Free the children. */
    release_subexpr_memory(sel);
    _gmx_selelem_free_chain(sel->child);
    sel->child           = NULL;
    /* Free the expression data as it is no longer needed */
    _gmx_selelem_free_exprdata(sel);
    /* Make the item static */
    sel->name            = NULL;
    sel->type            = SEL_CONST;
    sel->evaluate        = NULL;
    sel->cdata->evaluate = NULL;
    /* Set the group value.
     * free_exprdata above frees the cgrp group, so we can just override it. */
    if (sel->v.type == GROUP_VALUE)
    {
        gmx_ana_index_set(&sel->u.cgrp, sel->v.u.g->isize, sel->v.u.g->index, NULL, 0);
    }
}

/*! \brief
 * Evaluates a constant expression during analyze_static().
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
    gmx_bool       bAtomVal;
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
            if (sel->v.type != POS_VALUE && sel->v.type != GROUP_VALUE
                && !(sel->flags & SEL_VARNUMVAL))
            {
                alloc_selection_data(sel, isize, TRUE);
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
        /* Clear the values for dynamic output to avoid valgrind warnings. */
        if ((sel->flags & SEL_DYNAMIC) && sel->v.type == REAL_VALUE)
        {
            int i;

            for (i = 0; i < sel->v.nr; ++i)
            {
                sel->v.u.r[i] = 0.0;
            }
        }
    }

    return 0;
}

/*! \brief
 * Evaluates the static part of a gmx_boolean expression.
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
evaluate_gmx_boolean_static_part(gmx_sel_evaluate_t *data, t_selelem *sel,
                             gmx_ana_index_t *g)
{
    t_selelem *child, *next;
    int        rc;

    /* Find the last static subexpression */
    child = sel->child;
    while (child->next && (child->next->cdata->flags & SEL_CDATA_STATIC))
    {
        child = child->next;
    }
    if (!(child->cdata->flags & SEL_CDATA_STATIC))
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
        init_item_minmax_groups(child);
        child->cdata->flags &= ~SEL_CDATA_STATICEVAL;
        child->cdata->flags |= sel->cdata->flags & SEL_CDATA_STATICEVAL;
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
        || ((sel->cdata->flags & SEL_CDATA_STATICEVAL)
            && sel->u.boolt == BOOL_OR))
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
 * Evaluates the minimum and maximum groups for a gmx_boolean expression.
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
evaluate_gmx_boolean_minmax_grps(t_selelem *sel, gmx_ana_index_t *g,
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
            if ((sel->child->cdata->flags & SEL_CDATA_STATIC)
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
            if ((sel->child->cdata->flags & SEL_CDATA_STATIC)
                && sel->child->v.u.g->isize < gmin->isize)
            {
                gmx_ana_index_reserve(sel->child->v.u.g, gmin->isize);
                gmx_ana_index_copy(sel->child->v.u.g, gmin, FALSE);
                if (sel->child->u.cgrp.nalloc_index > 0)
                {
                    /* Keep the name as in evaluate_boolean_static_part(). */
                    char *name = sel->child->u.cgrp.name;
                    gmx_ana_index_reserve(&sel->child->u.cgrp, gmin->isize);
                    gmx_ana_index_copy(&sel->child->u.cgrp, gmin, FALSE);
                    sel->child->u.cgrp.name = name;
                }
                else
                {
                    sel->child->u.cgrp.isize = sel->child->v.u.g->isize;
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
 * once and whose evaluation group is not known in advance.
 */
static int
analyze_static(gmx_sel_evaluate_t *data, t_selelem *sel, gmx_ana_index_t *g)
{
    t_selelem       *child, *next;
    gmx_bool             bDoMinMax;
    int              rc;

    if (sel->type != SEL_ROOT && g)
    {
        alloc_selection_data(sel, g->isize, FALSE);
    }

    bDoMinMax = (sel->cdata->flags & SEL_CDATA_MINMAXALLOC);
    if (sel->type != SEL_SUBEXPR && bDoMinMax)
    {
        gmx_ana_index_deinit(sel->cdata->gmin);
        gmx_ana_index_deinit(sel->cdata->gmax);
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
            rc = init_method(sel, data->top, g ? g->isize : 0);
            if (rc != 0)
            {
                return rc;
            }
            if (!(sel->flags & SEL_DYNAMIC))
            {
                rc = sel->cdata->evaluate(data, sel, g);
                if (rc == 0 && (sel->cdata->flags & SEL_CDATA_STATIC))
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
                if (bDoMinMax && g)
                {
                    gmx_ana_index_copy(sel->cdata->gmax, g, TRUE);
                }
            }
            break;

        case SEL_BOOLEAN:
            if (!(sel->flags & SEL_DYNAMIC))
            {
                rc = sel->cdata->evaluate(data, sel, g);
                if (rc == 0 && (sel->cdata->flags & SEL_CDATA_STATIC))
                {
                    make_static(sel);
                }
            }
            else
            {
                /* Evalute the static part if there is more than one expression */
                rc = evaluate_gmx_boolean_static_part(data, sel, g);
                if (rc != 0)
                {
                    return rc;
                }

                /* Evaluate the selection.
                 * If the type is gmx_boolean, we must explicitly handle the
                 * static part evaluated in evaluate_gmx_boolean_static_part()
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
                evaluate_gmx_boolean_minmax_grps(sel, g, sel->cdata->gmin,
                                             sel->cdata->gmax);
            }
            break;

        case SEL_ARITHMETIC:
            rc = sel->cdata->evaluate(data, sel, g);
            if (rc != 0)
            {
                return rc;
            }
            if (!(sel->flags & SEL_DYNAMIC))
            {
                if (sel->cdata->flags & SEL_CDATA_STATIC)
                {
                    make_static(sel);
                }
            }
            else if (bDoMinMax)
            {
                gmx_ana_index_copy(sel->cdata->gmax, g, TRUE);
            }
            break;

        case SEL_ROOT:
            rc = sel->cdata->evaluate(data, sel, g);
            break;

        case SEL_SUBEXPR:
            if (sel->cdata->flags & (SEL_CDATA_SIMPLESUBEXPR | SEL_CDATA_FULLEVAL))
            {
                rc = sel->cdata->evaluate(data, sel, g);
                _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
            }
            else if (sel->u.cgrp.isize == 0)
            {
                gmx_ana_index_reserve(&sel->u.cgrp, g->isize);
                rc = sel->cdata->evaluate(data, sel, g);
                if (bDoMinMax)
                {
                    gmx_ana_index_copy(sel->cdata->gmin, sel->child->cdata->gmin, TRUE);
                    gmx_ana_index_copy(sel->cdata->gmax, sel->child->cdata->gmax, TRUE);
                }
            }
            else
            {
                int isize = gmx_ana_index_difference_size(g, &sel->u.cgrp);
                if (isize > 0)
                {
                    isize += sel->u.cgrp.isize;
                    gmx_ana_index_reserve(&sel->u.cgrp, isize);
                    alloc_selection_data(sel, isize, FALSE);
                }
                rc = sel->cdata->evaluate(data, sel, g);
                if (isize > 0 && bDoMinMax)
                {
                    gmx_ana_index_reserve(sel->cdata->gmin,
                                          sel->cdata->gmin->isize
                                          + sel->child->cdata->gmin->isize);
                    gmx_ana_index_reserve(sel->cdata->gmax,
                                          sel->cdata->gmax->isize
                                          + sel->child->cdata->gmax->isize);
                    gmx_ana_index_merge(sel->cdata->gmin, sel->cdata->gmin,
                                        sel->child->cdata->gmin);
                    gmx_ana_index_merge(sel->cdata->gmax, sel->cdata->gmax,
                                        sel->child->cdata->gmax);
                }
            }
            break;

        case SEL_SUBEXPRREF:
            if (!g && !(sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
            {
                /* The subexpression should have been evaluated if g is NULL
                 * (i.e., this is a method parameter or a direct value of a
                 * selection). */
                if (sel->v.type == POS_VALUE)
                {
                    alloc_selection_pos_data(sel);
                }
                else
                {
                    alloc_selection_data(sel, sel->child->cdata->gmax->isize, TRUE);
                }
            }
            rc = sel->cdata->evaluate(data, sel, g);
            if (rc != 0)
            {
                return rc;
            }
            if ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                && (sel->child->child->flags & SEL_ALLOCVAL))
            {
                _gmx_selvalue_setstore(&sel->v, sel->child->child->v.u.ptr);
            }
            /* Store the parameter value if required */
            store_param_val(sel);
            if (!(sel->flags & SEL_DYNAMIC))
            {
                if (sel->cdata->flags & SEL_CDATA_STATIC)
                {
                    make_static(sel);
                }
            }
            else if (bDoMinMax)
            {
                if ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR) || !g)
                {
                    gmx_ana_index_copy(sel->cdata->gmin, sel->child->cdata->gmin, TRUE);
                    gmx_ana_index_copy(sel->cdata->gmax, sel->child->cdata->gmax, TRUE);
                }
                else
                {
                    gmx_ana_index_reserve(sel->cdata->gmin,
                                          min(g->isize, sel->child->cdata->gmin->isize));
                    gmx_ana_index_reserve(sel->cdata->gmax,
                                          min(g->isize, sel->child->cdata->gmax->isize));
                    gmx_ana_index_intersection(sel->cdata->gmin,
                                               sel->child->cdata->gmin, g);
                    gmx_ana_index_intersection(sel->cdata->gmax,
                                               sel->child->cdata->gmax, g);
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
    if (bDoMinMax)
    {
        gmx_ana_index_squeeze(sel->cdata->gmin);
        gmx_ana_index_squeeze(sel->cdata->gmax);
        sfree(sel->cdata->gmin->name);
        sfree(sel->cdata->gmax->name);
        sel->cdata->gmin->name = NULL;
        sel->cdata->gmax->name = NULL;
    }

    /* Replace the result of the evaluation */
    /* This is not necessary for subexpressions or for gmx_boolean negations
     * because the evaluation function already has done it properly. */
    if (sel->v.type == GROUP_VALUE && (sel->flags & SEL_DYNAMIC)
        && sel->type != SEL_SUBEXPR
        && !(sel->type == SEL_BOOLEAN && sel->u.boolt == BOOL_NOT))
    {
        if (sel->cdata->flags & SEL_CDATA_EVALMAX)
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
 * EVALUATION GROUP INITIALIZATION
 ********************************************************************/

/*! \brief
 * Initializes the evaluation group for a \ref SEL_ROOT element.
 *
 * \param     root Root element to initialize.
 * \param[in] gall Group of all atoms.
 *
 * Checks whether it is necessary to evaluate anything through the root
 * element, and either clears the evaluation function or initializes the
 * evaluation group.
 */
static void
init_root_item(t_selelem *root, gmx_ana_index_t *gall)
{
    t_selelem   *expr;
    char        *name;

    expr = root->child;
    /* Subexpressions with non-static evaluation group should not be
     * evaluated by the root, and neither should be single-reference
     * subexpressions that don't evaluate for all atoms. */
    if (expr->type == SEL_SUBEXPR
        && (!(root->child->cdata->flags & SEL_CDATA_STATICEVAL)
            || ((root->child->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                && !(root->child->cdata->flags & SEL_CDATA_FULLEVAL))))
    {
        root->evaluate = NULL;
        if (root->cdata)
        {
            root->cdata->evaluate = NULL;
        }
    }

    /* Set the evaluation group */
    name = root->u.cgrp.name;
    if (root->evaluate)
    {
        /* Non-atom-valued non-group expressions don't care about the group, so
         * don't allocate any memory for it. */
        if ((expr->flags & SEL_VARNUMVAL)
            || ((expr->flags & SEL_SINGLEVAL) && expr->v.type != GROUP_VALUE))
        {
            gmx_ana_index_set(&root->u.cgrp, -1, NULL, NULL, 0);
        }
        else if (expr->cdata->gmax->isize == gall->isize)
        {
            /* Save some memory by only referring to the global group. */
            gmx_ana_index_set(&root->u.cgrp, gall->isize, gall->index, NULL, 0);
        }
        else
        {
            gmx_ana_index_copy(&root->u.cgrp, expr->cdata->gmax, TRUE);
        }
        /* For selections, store the maximum group for
         * gmx_ana_selcollection_evaluate_fin() as the value of the root
         * element (unused otherwise). */
        if (expr->type != SEL_SUBEXPR && expr->v.u.p->g)
        {
            t_selelem *child = expr;

            /* TODO: This code is copied from parsetree.c; it would be better
             * to have this hardcoded only in one place. */
            while (child->type == SEL_MODIFIER)
            {
                child = child->child;
                if (child->type == SEL_SUBEXPRREF)
                {
                    child = child->child->child;
                }
            }
            if (child->type == SEL_SUBEXPRREF)
            {
                child = child->child->child;
            }
            if (child->child->flags & SEL_DYNAMIC)
            {
                _gmx_selelem_set_vtype(root, GROUP_VALUE);
                root->flags  |= (SEL_ALLOCVAL | SEL_ALLOCDATA);
                _gmx_selvalue_reserve(&root->v, 1);
                gmx_ana_index_copy(root->v.u.g, expr->v.u.p->g, TRUE);
            }
        }
    }
    else
    {
        gmx_ana_index_clear(&root->u.cgrp);
    }
    root->u.cgrp.name = name;
}


/********************************************************************
 * FINAL SUBEXPRESSION OPTIMIZATION
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
postprocess_item_subexpressions(t_selelem *sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        t_selelem *child;

        child = sel->child;
        while (child)
        {
            postprocess_item_subexpressions(child);
            child = child->next;
        }
    }

    /* Replace the evaluation function of statically evaluated subexpressions
     * for which the static group was not known in advance. */
    if (sel->type == SEL_SUBEXPR && sel->refcount > 2
        && (sel->cdata->flags & SEL_CDATA_STATICEVAL)
        && !(sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        char *name;

        /* We need to free memory allocated for the group, because it is no
         * longer needed (and would be lost on next call to the evaluation
         * function). But we need to preserve the name. */
        name = sel->u.cgrp.name;
        gmx_ana_index_deinit(&sel->u.cgrp);
        sel->u.cgrp.name = name;

        sel->evaluate = &_gmx_sel_evaluate_subexpr_staticeval;
        if (sel->cdata)
        {
            sel->cdata->evaluate = sel->evaluate;
        }
        _gmx_selelem_free_values(sel->child);
        sel->child->mempool = NULL;
        _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
        sel->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
    }

    /* Adjust memory allocation flags for subexpressions that are used only
     * once.  This is not strictly necessary, but we do it to have the memory
     * managed consistently for all types of subexpressions. */
    if (sel->type == SEL_SUBEXPRREF
        && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
    {
        if (sel->child->child->flags & SEL_ALLOCVAL)
        {
            sel->flags |= SEL_ALLOCVAL;
            sel->flags |= (sel->child->child->flags & SEL_ALLOCDATA);
            sel->v.nalloc = sel->child->child->v.nalloc;
            sel->child->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
            sel->child->child->v.nalloc = -1;
        }
    }

    /* Do the same for subexpressions that are evaluated at once for all atoms. */
    if (sel->type == SEL_SUBEXPR
        && !(sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
        && (sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        sel->flags |= SEL_ALLOCVAL;
        sel->flags |= (sel->child->flags & SEL_ALLOCDATA);
        sel->v.nalloc = sel->child->v.nalloc;
        sel->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        sel->child->v.nalloc = -1;
    }
}


/********************************************************************
 * COM CALCULATION INITIALIZATION
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
            if (!(sel->cdata->flags & SEL_CDATA_STATICEVAL))
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
 * COMPILER DATA FREEING
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
 * MASS AND CHARGE CALCULATION
 ********************************************************************/

/*! \brief
 * Initializes total masses and charges for selections.
 *
 * \param[in]     top   Topology information.
 * \param[in]     ngrps Number of elements in the \p sel array.
 * \param[in,out] sel   Array of selections to update.
 * \param[in]     bMaskOnly TRUE if the positions will always be calculated
 *   for all atoms, i.e., the masses/charges do not change.
 */
static void
calculate_mass_charge(t_topology *top, int ngrps, gmx_ana_selection_t *sel[],
                      gmx_bool bMaskOnly)
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
 * \param[in,out] sc     Selection collection to debug.
 * \param[in]     bDebug If TRUE, later call to gmx_ana_selcollection_compile()
 *     will print out intermediate selection trees.
 */
void
gmx_ana_selcollection_set_compile_debug(gmx_ana_selcollection_t *sc, gmx_bool bDebug)
{
    sc->bDebugCompile = bDebug;
}

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

    rc = _gmx_sel_mempool_create(&sc->mempool);
    if (rc != 0)
    {
        return rc;
    }
    _gmx_sel_evaluate_init(&evaldata, sc->mempool, &sc->gall,
                           sc->top, NULL, NULL);

    /* Clear the symbol table because it is not possible to parse anything
     * after compilation, and variable references in the symbol table can
     * also mess up the compilation and/or become invalid.
     */
    _gmx_selcollection_clear_symtab(sc);

    /* Remove any unused variables. */
    sc->root = remove_unused_subexpressions(sc->root);
    /* Extract subexpressions into separate roots */
    sc->root = extract_subexpressions(sc->root);

    /* Initialize the evaluation callbacks and process the tree structure
     * to conform to the expectations of the callback functions. */
    /* Also, initialize and allocate the compiler data structure */
    item = sc->root;
    while (item)
    {
        /* Process gmx_boolean and arithmetic expressions. */
        optimize_gmx_boolean_expressions(item);
        reorder_gmx_boolean_static_children(item);
        if (!optimize_arithmetic_expressions(item))
        {
            /* FIXME: Clean up the collection */
            return -1;
        }
        /* Initialize evaluation function. */
        if (!init_item_evalfunc(item))
        {
            /* FIXME: Clean up the collection */
            return -1;
        }
        setup_memory_pooling(item, sc->mempool);
        /* Initialize the compiler data */
        init_item_compilerdata(item);
        init_item_staticeval(item);
        item = item->next;
    }
    /* Initialize subexpression flags and evaluation output.
     * Requires compiler flags for the full tree. */
    item = sc->root;
    while (item)
    {
        init_item_subexpr_flags(item);
        init_item_evaloutput(item);
        item = item->next;
    }
    /* Initialize minimum/maximum index groups.
     * Requires evaluation output for the full tree. */
    item = sc->root;
    while (item)
    {
        init_item_minmax_groups(item);
        item = item->next;
    }
    /* Initialize the evaluation index groups */
    initialize_evalgrps(sc);

    if (sc->bDebugCompile)
    {
        fprintf(stderr, "\nTree after initial compiler processing:\n");
        gmx_ana_selcollection_print_tree(stderr, sc, FALSE);
    }

    /* Evaluate all static parts of the selection and analyze the tree
     * to allocate enough memory to store the value of each dynamic subtree. */
    item = sc->root;
    while (item)
    {
        if (item->child->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
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

    /* At this point, static subexpressions no longer have references to them,
     * so they can be removed. */
    sc->root = remove_unused_subexpressions(sc->root);

    if (sc->bDebugCompile)
    {
        fprintf(stderr, "\nTree after first analysis pass:\n");
        gmx_ana_selcollection_print_tree(stderr, sc, FALSE);
    }

    /* Do a second pass to evaluate static parts of common subexpressions */
    item = sc->root;
    while (item)
    {
        if (item->child->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
        {
            gmx_bool bMinMax = item->child->cdata->flags & SEL_CDATA_MINMAXALLOC;

            mark_subexpr_dynamic(item->child, FALSE);
            item->child->u.cgrp.isize = 0;
            /* We won't clear item->child->v.u.g here, because it may
             * be static, and hence actually point to item->child->cdata->gmax,
             * which is used below. We could also check whether this is the
             * case and only clear the group otherwise, but because the value
             * is actually overwritten immediately in the evaluate call, we
             * won't, because similar problems may arise if gmax handling ever
             * changes and the check were not updated.
             * For the same reason, we clear the min/max flag so that the
             * evaluation group doesn't get messed up. */
            set_evaluation_function(item, &analyze_static);
            item->child->cdata->flags &= ~SEL_CDATA_MINMAXALLOC;
            rc = item->evaluate(&evaldata, item->child, item->child->cdata->gmax);
            if (bMinMax)
            {
                item->child->cdata->flags |= SEL_CDATA_MINMAXALLOC;
            }
            if (rc != 0)
            {
                /* FIXME: Clean up the collection */
                return rc;
            }
        }
        item = item->next;
    }

    /* We need a yet another pass of subexpression removal to remove static
     * subexpressions referred to by common dynamic subexpressions. */
    sc->root = remove_unused_subexpressions(sc->root);

    if (sc->bDebugCompile)
    {
        fprintf(stderr, "\nTree after second analysis pass:\n");
        gmx_ana_selcollection_print_tree(stderr, sc, FALSE);
    }

    /* Initialize evaluation groups, position calculations for methods, perform
     * some final optimization, and free the memory allocated for the
     * compilation. */
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
        init_root_item(item, &sc->gall);
        postprocess_item_subexpressions(item);
        rc = init_item_comg(item, sc->pcc, post, flags);
        if (rc != 0)
        {
            /* FIXME: Clean up the collection */
            return rc;
        }
        free_item_compilerdata(item);
        item = item->next;
    }

    /* Allocate memory for the evaluation memory pool. */
    rc = _gmx_sel_mempool_reserve(sc->mempool, 0);
    if (rc != 0)
    {
        return rc;
    }

    /* Finish up by calculating total masses and charges. */
    calculate_mass_charge(sc->top, sc->nr, sc->sel, sc->bMaskOnly);

    return 0;
}
