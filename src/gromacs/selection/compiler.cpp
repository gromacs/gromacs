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
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
/*! \internal
 * \page page_module_selection_compiler Selection compilation
 *
 * The compiler takes the selection element tree from the selection parser
 * (see \ref page_module_selection_parser) as input.
 * The selection parser is quite independent of selection evaluation details,
 * and the compiler processes the tree to conform to what the evaluation
 * functions expect.
 * For better control and optimization possibilities, the compilation is
 * done on all selections simultaneously.
 * Hence, all the selections should be parsed before the compiler can be
 * called.
 *
 * The compiler initializes all fields in gmx::SelectionTreeElement not
 * initialized by the parser: gmx::SelectionTreeElement::v (some fields have
 * already been initialized by the parser),
 * gmx::SelectionTreeElement::evaluate, and gmx::SelectionTreeElement::u
 * (again, some elements have been initialized in the parser).
 * The gmx::SelectionTreeElement::cdata field is used during the compilation to store
 * internal data, but the data is freed when the compiler returns.
 *
 * In addition to initializing the elements, the compiler reorganizes the tree
 * to simplify and optimize evaluation. The compiler also evaluates the static
 * parts of the selection: in the end of the compilation, static parts have
 * been replaced by the result of the evaluation.
 *
 * The compiler is invoked using gmx::SelectionCompiler.
 * The gmx::SelectionCompiler::compile() method does the compilation in several
 * passes over the gmx::SelectionTreeElement tree.
 *  -# Defaults are set for the position type and flags of position calculation
 *     methods that were not explicitly specified in the user input.
 *  -# Subexpressions are extracted: a separate root is created for each
 *     subexpression, and placed before the expression is first used.
 *     Currently, only variables and expressions used to evaluate parameter
 *     values are extracted, but common subexpression could also be detected
 *     here.
 *  -# A second pass (in fact, multiple passes because of interdependencies)
 *     with simple reordering and initialization is done:
 *    -# Boolean expressions are combined such that one element can evaluate,
 *       e.g., "A and B and C". The subexpressions in boolean expression are
 *       reordered such that static expressions come first without otherwise
 *       altering the relative order of the expressions.
 *    -# The compiler data structure is allocated for each element, and
 *       the fields are initialized, with the exception of the contents of
 *       \c gmax and \c gmin fields.  This is the part that needs multiple
 *       passes, because some flags are set recursively based on which elements
 *       refer to an element, and these flags need to be set to initialize
 *       other fields.
 *    -# The gmx::SelectionTreeElement::evaluate field is set to the correct
 *       evaluation function from evaluate.h.
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
 *     in pass 4.  For root elements of subexpressions that should not be
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
 * the compiler or the evaluation functions in evaluate.cpp. All the ones used by
 * currently implemented methods should work, but new combinations might not.
 *
 *
 * \section selcompiler_tree Element tree after compilation
 *
 * After the compilation, the selection element tree is suitable for
 * gmx_ana_selcollection_evaluate().
 * Enough memory has been allocated for gmx::SelectionTreeElement::v
 * (and gmx::SelectionTreeElement::cgrp for \ref SEL_SUBEXPR elements) to allow
 * the selection to be evaluated without allocating any memory.
 *
 *
 * \subsection selcompiler_tree_root Root elements
 *
 * The top level of the tree consists of a chain of \ref SEL_ROOT elements.
 * These are used for two purposes:
 *  -# A selection that should be evaluated.
 *     These elements appear in the same order as the selections in the input.
 *     For these elements, gmx::SelectionTreeElement::v has been set to the
 *     maximum possible group that the selection can evaluate to (only for
 *     dynamic selections), and gmx::SelectionTreeElement::cgrp has been set to
 *     use a NULL group for evaluation.
 *  -# A subexpression that appears in one or more selections.
 *     Each selection that gives a value for a method parameter is a
 *     potential subexpression, as is any variable value.
 *     Only subexpressions that require evaluation for each frame are left
 *     after the selection is compiled.
 *     Each subexpression appears in the chain before any references to it.
 *     For these elements, gmx::SelectionTreeElement::cgrp has been set to the
 *     group that should be used to evaluate the subexpression.
 *     If gmx::SelectionTreeElement::cgrp is empty, the total evaluation group
 *     is not known in advance or it is more efficient to evaluate the
 *     subexpression only when it is referenced.  If this is the case,
 *     gmx::SelectionTreeElement::evaluate is also NULL.
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
 * Several constant elements with a NULL gmx::SelectionTreeElement::evaluate
 * are left for debugging purposes; of these, only the ones for \ref BOOL_OR
 * expressions are used during evaluation.
 *
 * The value is stored in gmx::SelectionTreeElement::v, and for group values
 * with an evaluation function set, also in gmx::SelectionTreeElement::cgrp.
 * For \ref GROUP_VALUE elements, unnecessary atoms (i.e., atoms that
 * could never be selected) have been removed from the value.
 *
 * \ref SEL_CONST elements have no children.
 *
 *
 * \subsection selcompiler_tree_method Method evaluation elements
 *
 * All selection methods that need to be evaluated dynamically are described
 * by a \ref SEL_EXPRESSION element. The gmx::SelectionTreeElement::method and
 * gmx::SelectionTreeElement::mdata fields have already been initialized by the parser,
 * and the compiler only calls the initialization functions in the method
 * data structure to do some additional initialization of these fields at
 * appropriate points. If the gmx::SelectionTreeElement::pc data field has been
 * created by the parser, the compiler initializes the data structure properly
 * once the required positions are known. If the gmx::SelectionTreeElement::pc
 * field is NULL after the parser, but the method provides only
 * sel_updatefunc_pos(), an appropriate position calculation data structure is
 * created.  If gmx::SelectionTreeElement::pc is not NULL,
 * gmx::SelectionTreeElement::pos is also initialized to hold the positions
 * calculated.
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
 * gmx::SelectionTreeElement::cgrp to store the group for which the expression
 * has been evaluated during the current frame.  This is only done if full
 * subexpression evaluation by _gmx_sel_evaluate_subexpr() is needed; the other
 * evaluation functions do not require this memory.
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
 * \subsection selcompiler_tree_bool Boolean elements
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
#include "gmxpre.h"

#include "compiler.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>

#include <algorithm>
#include <memory>
#include <string>

#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "evaluate.h"
#include "keywords.h"
#include "mempool.h"
#include "poscalc.h"
#include "selectioncollection_impl.h"
#include "selelem.h"
#include "selmethod.h"

struct gmx_mtop_t;

using gmx::SelectionLocation;
using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;
using std::min;

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
    SEL_CDATA_FULLEVAL = 1,
    /*! \brief
     * Whether the whole subexpression should be treated as static.
     *
     * This flag is always false if \ref SEL_DYNAMIC is set for the element,
     * but it is also false for static elements within common subexpressions.
     */
    SEL_CDATA_STATIC = 2,
    /** Whether the subexpression will always be evaluated in the same group. */
    SEL_CDATA_STATICEVAL = 4,
    /** Whether the compiler evaluation routine should return the maximal selection. */
    SEL_CDATA_EVALMAX = 8,
    /** Whether memory has been allocated for \p gmin and \p gmax. */
    SEL_CDATA_MINMAXALLOC = 16,
    /** Whether to update \p gmin and \p gmax in static analysis. */
    SEL_CDATA_DOMINMAX = 256,
    /** Whether the subexpression uses simple pass evaluation functions. */
    SEL_CDATA_SIMPLESUBEXPR = 32,
    /*! \brief
     * Whether a static subexpression needs to support multiple evaluations.
     *
     * This flag may only be set on \ref SEL_SUBEXPR elements that also have
     * SEL_CDATA_SIMPLESUBEXPR.
     */
    SEL_CDATA_STATICMULTIEVALSUBEXPR = 64,
    /** Whether this expression is a part of a common subexpression. */
    SEL_CDATA_COMMONSUBEXPR = 128
};

/*! \internal \brief
 * Internal data structure used by the compiler.
 */
typedef struct t_compiler_data
{
    /** The real evaluation method. */
    gmx::sel_evalfunc evaluate;
    /** Number of references to a \ref SEL_SUBEXPR element. */
    int refcount;
    /** Flags for specifying how to treat this element during compilation. */
    int flags;
    /** Smallest selection that can be selected by the subexpression. */
    gmx_ana_index_t* gmin;
    /** Largest selection that can be selected by the subexpression. */
    gmx_ana_index_t* gmax;
} t_compiler_data;


/********************************************************************
 * COMPILER UTILITY FUNCTIONS
 ********************************************************************/

/*! \brief
 * Helper method for printing out debug information about a min/max group.
 */
static void print_group_info(FILE* fp, const char* name, const SelectionTreeElement& sel, gmx_ana_index_t* g)
{
    fprintf(fp, " %s=", name);
    if (!g)
    {
        fprintf(fp, "(null)");
    }
    else if (sel.cdata->flags & SEL_CDATA_MINMAXALLOC)
    {
        fprintf(fp, "(%d atoms, %p)", g->isize, static_cast<void*>(g));
    }
    else if (sel.v.type == GROUP_VALUE && g == sel.v.u.g)
    {
        fprintf(fp, "(static, %p)", static_cast<void*>(g));
    }
    else
    {
        fprintf(fp, "%p", static_cast<void*>(g));
    }
}

/*!
 * \param[in] fp      File handle to receive the output.
 * \param[in] sel     Selection element to print.
 * \param[in] level   Indentation level, starting from zero.
 */
void _gmx_selelem_print_compiler_info(FILE* fp, const SelectionTreeElement& sel, int level)
{
    if (!sel.cdata)
    {
        return;
    }
    fprintf(fp, "%*c cdata: flg=", level * 2 + 1, ' ');
    if (sel.cdata->flags & SEL_CDATA_FULLEVAL)
    {
        fprintf(fp, "F");
    }
    if (!(sel.cdata->flags & SEL_CDATA_STATIC))
    {
        fprintf(fp, "D");
    }
    if (sel.cdata->flags & SEL_CDATA_STATICEVAL)
    {
        fprintf(fp, "S");
    }
    if (sel.cdata->flags & SEL_CDATA_EVALMAX)
    {
        fprintf(fp, "M");
    }
    if (sel.cdata->flags & SEL_CDATA_MINMAXALLOC)
    {
        fprintf(fp, "A");
    }
    if (sel.cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
    {
        fprintf(fp, "Ss");
    }
    if (sel.cdata->flags & SEL_CDATA_STATICMULTIEVALSUBEXPR)
    {
        fprintf(fp, "Sm");
    }
    if (sel.cdata->flags & SEL_CDATA_COMMONSUBEXPR)
    {
        fprintf(fp, "Sc");
    }
    if (!sel.cdata->flags)
    {
        fprintf(fp, "0");
    }
    if (sel.cdata->refcount > 0)
    {
        fprintf(fp, " refc=%d", sel.cdata->refcount);
    }
    fprintf(fp, " eval=");
    _gmx_sel_print_evalfunc_name(fp, sel.cdata->evaluate);
    print_group_info(fp, "gmin", sel, sel.cdata->gmin);
    print_group_info(fp, "gmax", sel, sel.cdata->gmax);
    fprintf(fp, "\n");
}

namespace gmx
{

void SelectionTreeElement::freeCompilerData()
{
    if (cdata)
    {
        evaluate = cdata->evaluate;
        if (cdata->flags & SEL_CDATA_MINMAXALLOC)
        {
            gmx_ana_index_deinit(cdata->gmin);
            gmx_ana_index_deinit(cdata->gmax);
            sfree(cdata->gmin);
            sfree(cdata->gmax);
        }
        sfree(cdata);
    }
    cdata = nullptr;
}

} // namespace gmx

/*! \brief
 * Allocates memory for storing the evaluated value of a selection element.
 *
 * \param     sel   Selection element to initialize
 * \param[in] isize Maximum evaluation group size.
 * \param[in] bChildEval true if children have already been processed.
 *
 * If called more than once, memory is (re)allocated to ensure that the
 * maximum of the \p isize values can be stored.
 *
 * Allocation of POS_VALUE selection elements is a special case, and is
 * handled by alloc_selection_pos_data().
 */
static void alloc_selection_data(const SelectionTreeElementPointer& sel, int isize, bool bChildEval)
{
    int nalloc;

    GMX_RELEASE_ASSERT(sel->v.type != POS_VALUE, "Wrong allocation method called");
    if (sel->mempool)
    {
        return;
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
        // TODO: Consider whether the bChildEval is any longer necessary.
        if (!bChildEval)
        {
            return;
        }
        SelectionTreeElementPointer child = sel;
        if (sel->type == SEL_SUBEXPRREF)
        {
            GMX_RELEASE_ASSERT(sel->child && sel->child->type == SEL_SUBEXPR,
                               "Subexpression expected for subexpression reference");
            child = sel->child->child;
            GMX_RELEASE_ASSERT(child, "Subexpression elements should always have a child element");
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
}

/*! \brief
 * Allocates memory for storing the evaluated value of a selection element.
 *
 * \param     sel   Selection element to initialize.
 *
 * Allocation of POS_VALUE selection elements is a special case, and is
 * handled by this function instead of by alloc_selection_data().
 */
static void alloc_selection_pos_data(const SelectionTreeElementPointer& sel)
{
    int nalloc, isize;

    GMX_RELEASE_ASSERT(sel->v.type == POS_VALUE, "Wrong allocation method called");
    GMX_RELEASE_ASSERT(!(sel->flags & SEL_ATOMVAL), "Per-atom evaluated positions not implemented");
    if (sel->mempool)
    {
        return;
    }

    SelectionTreeElementPointer child = sel;
    if (sel->type == SEL_SUBEXPRREF)
    {
        GMX_RELEASE_ASSERT(sel->child && sel->child->type == SEL_SUBEXPR,
                           "Subexpression expected for subexpression reference");
        child = sel->child->child;
        GMX_RELEASE_ASSERT(child, "Subexpression elements should always have a child element");
    }
    nalloc = child->v.u.p->count();
    isize  = child->v.u.p->m.b.nra;

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
static void set_evaluation_function(const SelectionTreeElementPointer& sel, gmx::sel_evalfunc eval)
{
    sel->evaluate = eval;
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            set_evaluation_function(child, eval);
            child = child->next;
        }
    }
}


/********************************************************************
 * POSITION KEYWORD DEFAULT INITIALIZATION
 ********************************************************************/

/*! \brief
 * Initializes default values for position keyword evaluation.
 *
 * \param[in,out] root       Root of the element tree to initialize.
 * \param[in]     spost      Default output position type.
 * \param[in]     rpost      Default reference position type.
 * \param[in]     sel        Selection that the element evaluates the positions
 *      for, or NULL if the element is an internal element.
 */
static void init_pos_keyword_defaults(SelectionTreeElement*               root,
                                      const char*                         spost,
                                      const char*                         rpost,
                                      const gmx::internal::SelectionData* sel)
{
    /* Selections use largest static group by default, while
     * reference positions use the whole residue/molecule. */
    if (root->type == SEL_EXPRESSION)
    {
        bool bSelection = (sel != nullptr);
        int  flags      = bSelection ? POS_COMPLMAX : POS_COMPLWHOLE;
        if (bSelection)
        {
            if (sel->hasFlag(gmx::efSelection_DynamicMask))
            {
                flags |= POS_MASKONLY;
            }
            if (sel->hasFlag(gmx::efSelection_EvaluateVelocities))
            {
                flags |= POS_VELOCITIES;
            }
            if (sel->hasFlag(gmx::efSelection_EvaluateForces))
            {
                flags |= POS_FORCES;
            }
        }
        _gmx_selelem_set_kwpos_type(root, bSelection ? spost : rpost);
        _gmx_selelem_set_kwpos_flags(root, flags);
    }
    /* Change the defaults once we are no longer processing modifiers */
    if (root->type != SEL_ROOT && root->type != SEL_MODIFIER && root->type != SEL_SUBEXPRREF
        && root->type != SEL_SUBEXPR)
    {
        sel = nullptr;
    }
    /* Recurse into children */
    SelectionTreeElementPointer child = root->child;
    while (child)
    {
        init_pos_keyword_defaults(child.get(), spost, rpost, sel);
        child = child->next;
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
static SelectionTreeElementPointer reverse_selelem_chain(const SelectionTreeElementPointer& root)
{
    SelectionTreeElementPointer prev;
    SelectionTreeElementPointer item = root;
    while (item)
    {
        SelectionTreeElementPointer next = item->next;
        item->next                       = prev;
        prev                             = item;
        item                             = next;
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
static SelectionTreeElementPointer remove_unused_subexpressions(SelectionTreeElementPointer root)
{
    if (!root)
    {
        return SelectionTreeElementPointer();
    }
    root = reverse_selelem_chain(root);
    while (root->child->type == SEL_SUBEXPR && root->child.use_count() == 1)
    {
        // Frees the root element.
        root = root->next;
    }
    SelectionTreeElementPointer prev = root;
    SelectionTreeElementPointer item = root->next;
    while (item)
    {
        SelectionTreeElementPointer next = item->next;
        if (item->child->type == SEL_SUBEXPR && item->child.use_count() == 1)
        {
            // Frees the current item when it goes out of scope.
            prev->next = next;
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
 * Memory is allocated for the name and the name is stored both as the
 * name of the subexpression element and as
 * gmx::SelectionTreeElement::u::cgrp::name; the latter is freed by
 * _gmx_selelem_free().
 */
static void create_subexpression_name(const SelectionTreeElementPointer& sel, int i)
{
    std::string name(gmx::formatString("SubExpr %d", i));
    sel->setName(name);
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
static SelectionTreeElementPointer extract_item_subselections(const SelectionTreeElementPointer& sel,
                                                              int* subexprn)
{
    SelectionTreeElementPointer root;
    SelectionTreeElementPointer subexpr;
    SelectionTreeElementPointer child = sel->child;

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
        /* The latter check excludes variable references. */
        if (child->type == SEL_SUBEXPRREF && child->child->type != SEL_SUBEXPR)
        {
            SelectionLocation location = child->child->location();
            /* Create the root element for the subexpression */
            if (!root)
            {
                root    = std::make_shared<SelectionTreeElement>(SEL_ROOT, location);
                subexpr = root;
            }
            else
            {
                subexpr->next = std::make_shared<SelectionTreeElement>(SEL_ROOT, location);
                subexpr       = subexpr->next;
            }
            /* Create the subexpression element and
             * move the actual subexpression under the created element. */
            subexpr->child = std::make_shared<SelectionTreeElement>(SEL_SUBEXPR, location);
            _gmx_selelem_set_vtype(subexpr->child, child->v.type);
            subexpr->child->child = child->child;
            child->child          = subexpr->child;
            create_subexpression_name(subexpr->child, ++*subexprn);
            /* Set the flags for the created elements */
            subexpr->flags |= (child->flags & SEL_VALFLAGMASK);
            subexpr->child->flags |= (child->flags & SEL_VALFLAGMASK);
        }
        if (child->type == SEL_SUBEXPRREF)
        {
            child->setName(child->child->name());
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
static SelectionTreeElementPointer extract_subexpressions(SelectionTreeElementPointer sel)
{
    SelectionTreeElementPointer root;
    SelectionTreeElementPointer next     = sel;
    int                         subexprn = 0;
    while (next)
    {
        SelectionTreeElementPointer item = extract_item_subselections(next, &subexprn);
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
        sel  = next;
        next = next->next;
    }
    return root;
}


/********************************************************************
 * BOOLEAN OPERATION REORDERING
 ********************************************************************/

/*! \brief
 * Removes redundant boolean selection elements.
 *
 * \param  sel Root of the selection subtree to optimize.
 *
 * This function merges similar boolean operations (e.g., (A or B) or C becomes
 * a single OR operation with three operands).
 */
static void optimize_boolean_expressions(const SelectionTreeElementPointer& sel)
{
    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer prev;
        SelectionTreeElementPointer child = sel->child;
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
                // Discards the two negations.
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
    SelectionTreeElementPointer prev;
    SelectionTreeElementPointer child = sel->child;
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
            // Discards the old child.
            child = prev->next;
        }
        else
        {
            prev  = child;
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
static void reorder_boolean_static_children(const SelectionTreeElementPointer& sel)
{
    /* Do recursively for children */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            reorder_boolean_static_children(child);
            child = child->next;
        }
    }

    /* Reorder boolean expressions such that static selections come first */
    if (sel->type == SEL_BOOLEAN && (sel->flags & SEL_DYNAMIC))
    {
        // Add a dummy head element that precedes the first child.
        SelectionTreeElementPointer dummy(
                new SelectionTreeElement(SEL_BOOLEAN, SelectionLocation::createEmpty()));
        dummy->next                       = sel->child;
        SelectionTreeElementPointer prev  = dummy;
        SelectionTreeElementPointer child = dummy;
        while (child->next)
        {
            /* child is the last handled static expression */
            /* prev is the last handled non-static expression */
            SelectionTreeElementPointer next = prev->next;
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

        sel->child = dummy->next;
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
static void optimize_arithmetic_expressions(const SelectionTreeElementPointer& sel)
{
    /* Do recursively for children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            optimize_arithmetic_expressions(child);
            child = child->next;
        }
    }

    if (sel->type != SEL_ARITHMETIC)
    {
        return;
    }

    /* Convert integer constants to reals. */
    SelectionTreeElementPointer child = sel->child;
    while (child)
    {
        if (child->v.type == INT_VALUE)
        {
            real* r;

            if (child->type != SEL_CONST)
            {
                GMX_THROW(
                        gmx::InconsistentInputError("Non-constant integer expressions not "
                                                    "implemented in arithmetic evaluation"));
            }
            snew(r, 1);
            r[0] = child->v.u.i[0];
            sfree(child->v.u.i);
            child->v.u.r  = r;
            child->v.type = REAL_VALUE;
        }
        else if (child->v.type != REAL_VALUE)
        {
            GMX_THROW(gmx::InternalError("Non-numerical value in arithmetic expression"));
        }
        child = child->next;
    }
}


/********************************************************************
 * EVALUATION PREPARATION COMPILER
 ********************************************************************/

/*! \brief
 * Sets the evaluation functions for the selection (sub)tree.
 *
 * \param[in,out] sel Root of the selection subtree to process.
 *
 * This function sets the evaluation function
 * (gmx::SelectionTreeElement::evaluate) for the selection elements.
 */
static void init_item_evalfunc(const SelectionTreeElementPointer& sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_evalfunc(child);
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
            if (!(sel->flags & SEL_DYNAMIC) && sel->u.expr.method && sel->u.expr.method->init_frame)
            {
                sel->flags |= SEL_INITFRAME;
            }
            sel->evaluate = &_gmx_sel_evaluate_method;
            break;

        case SEL_ARITHMETIC: sel->evaluate = &_gmx_sel_evaluate_arithmetic; break;

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
                case BOOL_OR: sel->evaluate = &_gmx_sel_evaluate_or; break;
                case BOOL_XOR:
                    GMX_THROW(gmx::NotImplementedError("xor expressions not implemented"));
            }
            break;

        case SEL_ROOT: sel->evaluate = &_gmx_sel_evaluate_root; break;

        case SEL_SUBEXPR:
            if ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                && !(sel->cdata->flags & SEL_CDATA_STATICMULTIEVALSUBEXPR))
            {
                sel->evaluate = &_gmx_sel_evaluate_subexpr_simple;
            }
            else
            {
                sel->evaluate = &_gmx_sel_evaluate_subexpr;
            }
            break;

        case SEL_SUBEXPRREF:
            sel->evaluate = ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                                     ? &_gmx_sel_evaluate_subexprref_simple
                                     : &_gmx_sel_evaluate_subexprref);
            break;

        case SEL_GROUPREF: GMX_THROW(gmx::APIError("Unresolved group reference in compilation"));
    }
    sel->cdata->evaluate = sel->evaluate;
}

/*! \brief
 * Sets the memory pool for selection elements that can use it.
 *
 * \param     sel      Root of the selection subtree to process.
 * \param[in] mempool  Memory pool to use.
 */
static void setup_memory_pooling(const SelectionTreeElementPointer& sel, gmx_sel_mempool_t* mempool)
{
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            if ((sel->type == SEL_BOOLEAN && (child->flags & SEL_DYNAMIC))
                || (sel->type == SEL_ARITHMETIC && child->type != SEL_CONST && !(child->flags & SEL_SINGLEVAL))
                || (sel->type == SEL_SUBEXPR && !(sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)))
            {
                child->mempool = mempool;
                if (child->type == SEL_SUBEXPRREF && (child->child->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
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
static void init_item_evaloutput(const SelectionTreeElementPointer& sel)
{
    GMX_ASSERT(!(sel->child == nullptr && (sel->type == SEL_SUBEXPRREF || sel->type == SEL_SUBEXPR)),
               "Subexpression elements should always have a child element");

    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_evaloutput(child);
            child = child->next;
        }
    }

    if (sel->type == SEL_SUBEXPR && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
        && !(sel->cdata->flags & SEL_CDATA_STATICMULTIEVALSUBEXPR))
    {
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
        }
    }
    else if (sel->type == SEL_SUBEXPR && (sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        sel->evaluate        = &_gmx_sel_evaluate_subexpr_staticeval;
        sel->cdata->evaluate = sel->evaluate;
        sel->child->mempool  = nullptr;
        sel->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        if (sel->v.type == GROUP_VALUE || sel->v.type == POS_VALUE)
        {
            _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
        }
    }
    else if (sel->type == SEL_SUBEXPRREF && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
    {
        if (sel->v.u.ptr)
        {
            _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
            sel->child->child->freeValues();
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
static void init_item_compilerdata(const SelectionTreeElementPointer& sel)
{
    /* Allocate the compiler data structure */
    snew(sel->cdata, 1);

    /* Store the real evaluation method because the compiler will replace it */
    sel->cdata->evaluate = sel->evaluate;

    /* This will be computed separately. */
    sel->cdata->refcount = 0;

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
        SelectionTreeElementPointer child = sel->child;
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
        SelectionTreeElementPointer child = sel->child;
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
        bool bEvalMax;

        bEvalMax                          = (sel->u.boolt == BOOL_AND);
        SelectionTreeElementPointer child = sel->child;
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
    else if (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER || sel->type == SEL_SUBEXPR)
    {
        SelectionTreeElementPointer child = sel->child;
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
 * reorder_boolean_static_children() should have been called.
 */
static void init_item_staticeval(const SelectionTreeElementPointer& sel)
{
    /* Subexpressions with full evaluation should always have bStaticEval,
     * so don't do anything if a reference to them is encountered. */
    if (sel->type == SEL_SUBEXPRREF && (sel->child->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        return;
    }

    /* Propagate the bStaticEval flag to children if it is not set */
    if (!(sel->cdata->flags & SEL_CDATA_STATICEVAL))
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            if ((sel->type != SEL_EXPRESSION && sel->type != SEL_MODIFIER) || (child->flags & SEL_ATOMVAL))
            {
                if (child->cdata->flags & SEL_CDATA_STATICEVAL)
                {
                    child->cdata->flags &= ~SEL_CDATA_STATICEVAL;
                    init_item_staticeval(child);
                }
            }
            /* If an expression is evaluated for a dynamic group, then also
             * atom-valued parameters need to be evaluated every time. */
            if ((sel->flags & SEL_DYNAMIC) && (sel->type == SEL_EXPRESSION || sel->type == SEL_MODIFIER)
                && (child->flags & SEL_ATOMVAL))
            {
                child->flags |= SEL_DYNAMIC;
                child->cdata->flags &= ~SEL_CDATA_STATIC;
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
            SelectionTreeElementPointer child = sel->child;
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
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_staticeval(child);
            child = child->next;
        }
    }
}

/*! \brief
 * Compute reference counts for subexpressions.
 *
 * \param sel Root of the selection subtree to process.
 */
static void init_item_subexpr_refcount(const SelectionTreeElementPointer& sel)
{
    // Reset the counter when the subexpression is first encountered.
    if (sel->type == SEL_ROOT && sel->child->type == SEL_SUBEXPR && sel->child->cdata)
    {
        sel->child->cdata->refcount = 0;
    }

    if (sel->type == SEL_SUBEXPRREF)
    {
        ++sel->child->cdata->refcount;
    }
    else
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_subexpr_refcount(child);
            child = child->next;
        }
    }
}

/*! \brief
 * Initializes compiler flags for subexpressions.
 *
 * \param sel Root of the selection subtree to process.
 */
static void init_item_subexpr_flags(const SelectionTreeElementPointer& sel)
{
    if (sel->type == SEL_SUBEXPR)
    {
        if (sel->cdata->refcount == 1)
        {
            sel->cdata->flags |= SEL_CDATA_SIMPLESUBEXPR;
        }
        else if (!(sel->cdata->flags & SEL_CDATA_FULLEVAL))
        {
            sel->cdata->flags |= SEL_CDATA_COMMONSUBEXPR;
        }
    }
    else if (sel->type == SEL_SUBEXPRREF && (sel->child->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
    {
        /* See similar condition in init_item_staticeval(). */
        if ((sel->flags & SEL_ATOMVAL) && (sel->flags & SEL_DYNAMIC) && !(sel->child->flags & SEL_DYNAMIC))
        {
            sel->child->cdata->flags |= SEL_CDATA_STATICMULTIEVALSUBEXPR;
        }
        else
        {
            sel->cdata->flags |= SEL_CDATA_SIMPLESUBEXPR;
        }
    }

    /* Process children, but only follow subexpression references if the
     * common subexpression flag needs to be propagated. */
    if (sel->type != SEL_SUBEXPRREF
        || ((sel->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
            && !(sel->child->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)))
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            if (!(child->cdata->flags & SEL_CDATA_COMMONSUBEXPR))
            {
                if (sel->type != SEL_EXPRESSION || (child->flags & SEL_ATOMVAL))
                {
                    child->cdata->flags |= (sel->cdata->flags & SEL_CDATA_COMMONSUBEXPR);
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
static void init_item_minmax_groups(const SelectionTreeElementPointer& sel)
{
    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_minmax_groups(child);
            child = child->next;
        }
    }

    /* Initialize the minimum and maximum evaluation groups. */
    if (sel->type != SEL_ROOT && sel->v.type != NO_VALUE)
    {
        if (sel->v.type == GROUP_VALUE && (sel->cdata->flags & SEL_CDATA_STATIC))
        {
            sel->cdata->gmin = sel->v.u.g;
            sel->cdata->gmax = sel->v.u.g;
        }
        else if (sel->type == SEL_SUBEXPR
                 && ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                     || (sel->cdata->flags & SEL_CDATA_FULLEVAL)))
        {
            GMX_ASSERT(sel->child, "Subexpression elements should always have a child element");
            sel->cdata->gmin = sel->child->cdata->gmin;
            sel->cdata->gmax = sel->child->cdata->gmax;
        }
        else
        {
            sel->cdata->flags |= SEL_CDATA_MINMAXALLOC | SEL_CDATA_DOMINMAX;
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
 * selection in \p sc is set to NULL.  The evaluation group for \ref SEL_ROOT
 * elements corresponding to subexpressions that need full evaluation is set
 * to \c sc->gall.
 */
static void initialize_evalgrps(gmx_ana_selcollection_t* sc)
{
    SelectionTreeElementPointer root = sc->root;
    while (root)
    {
        GMX_RELEASE_ASSERT(root->child, "Root elements should always have a child");
        if (root->child->type != SEL_SUBEXPR
            || (root->child->v.type != GROUP_VALUE && !(root->flags & SEL_ATOMVAL)))
        {
            gmx_ana_index_set(&root->u.cgrp, -1, nullptr, 0);
        }
        else if (root->child->cdata->flags & SEL_CDATA_FULLEVAL)
        {
            gmx_ana_index_set(&root->u.cgrp, sc->gall.isize, sc->gall.index, 0);
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
 * \param[in] bDynamic If true, the \p bStatic flag of the whole
 *   selection subtree is cleared. If false, the flag is restored to
 *   using \ref SEL_DYNAMIC.
 *
 * Does not descend into parameters of methods unless the parameters
 * are evaluated for each atom.
 */
static void mark_subexpr_dynamic(const SelectionTreeElementPointer& sel, bool bDynamic)
{
    if (!bDynamic && !(sel->flags & SEL_DYNAMIC))
    {
        sel->cdata->flags |= SEL_CDATA_STATIC;
    }
    else
    {
        sel->cdata->flags &= ~SEL_CDATA_STATIC;
    }
    SelectionTreeElementPointer child = sel->child;
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
static void release_subexpr_memory(const SelectionTreeElementPointer& sel)
{
    if (sel->type == SEL_SUBEXPRREF)
    {
        const SelectionTreeElementPointer& subexpr = sel->child;
        if (subexpr.use_count() == 2)
        {
            release_subexpr_memory(subexpr);
            // Free children.
            subexpr->child.reset();
            subexpr->freeValues();
            subexpr->freeExpressionData();
            subexpr->freeCompilerData();
        }
    }
    else
    {
        SelectionTreeElementPointer child = sel->child;
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
static void make_static(const SelectionTreeElementPointer& sel)
{
    /* If this is a subexpression reference and the data is stored in the
     * child, we transfer data ownership before doing anything else. */
    if (sel->type == SEL_SUBEXPRREF && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
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
    /* Free the children. */
    release_subexpr_memory(sel);
    sel->child.reset();
    /* Free the expression data as it is no longer needed */
    sel->freeExpressionData();
    /* Make the item static */
    sel->type            = SEL_CONST;
    sel->evaluate        = nullptr;
    sel->cdata->evaluate = nullptr;
    /* Set the group value.
     * freeExpressionData() frees the cgrp group, so we can just override it.
     * */
    if (sel->v.type == GROUP_VALUE)
    {
        gmx_ana_index_set(&sel->u.cgrp, sel->v.u.g->isize, sel->v.u.g->index, 0);
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
static void process_const(gmx_sel_evaluate_t* data, const SelectionTreeElementPointer& sel, gmx_ana_index_t* g)
{
    if (sel->v.type == GROUP_VALUE)
    {
        if (sel->cdata->evaluate)
        {
            sel->cdata->evaluate(data, sel, g);
        }
    }
    /* Other constant expressions do not need evaluation */
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
static void store_param_val(const SelectionTreeElementPointer& sel)
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

    if (sel->v.type == INT_VALUE || sel->v.type == REAL_VALUE || sel->v.type == STR_VALUE)
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
static void init_method(const SelectionTreeElementPointer& sel, const gmx_mtop_t* top, int isize)
{
    /* Find out whether there are any atom-valued parameters */
    bool                        bAtomVal = false;
    SelectionTreeElementPointer child    = sel->child;
    while (child)
    {
        if (child->flags & SEL_ATOMVAL)
        {
            bAtomVal = true;
        }
        child = child->next;
    }

    /* Initialize the method */
    if (sel->u.expr.method->init && (bAtomVal || !(sel->flags & SEL_METHODINIT)))
    {
        sel->flags |= SEL_METHODINIT;
        sel->u.expr.method->init(
                top, sel->u.expr.method->nparams, sel->u.expr.method->param, sel->u.expr.mdata);
    }
    if (bAtomVal || !(sel->flags & SEL_OUTINIT))
    {
        sel->flags |= SEL_OUTINIT;
        if (sel->u.expr.method->outinit)
        {
            sel->u.expr.method->outinit(top, &sel->v, sel->u.expr.mdata);
            if (sel->v.type != POS_VALUE && sel->v.type != GROUP_VALUE && !(sel->flags & SEL_VARNUMVAL))
            {
                alloc_selection_data(sel, isize, true);
            }
        }
        else
        {
            GMX_RELEASE_ASSERT(sel->v.type != POS_VALUE,
                               "Output initialization must be provided for "
                               "position-valued selection methods");
            GMX_RELEASE_ASSERT(!(sel->flags & SEL_VARNUMVAL),
                               "Output initialization must be provided for "
                               "SMETH_VARNUMVAL selection methods");
            alloc_selection_data(sel, isize, true);
            /* If the method is char-valued, pre-allocate the strings. */
            if (sel->u.expr.method->flags & SMETH_CHARVAL)
            {
                int i;

                /* A sanity check */
                if (sel->v.type != STR_VALUE)
                {
                    GMX_THROW(gmx::InternalError(
                            "Char-valued selection method in non-string element"));
                }
                sel->flags |= SEL_ALLOCDATA;
                for (i = 0; i < isize; ++i)
                {
                    if (sel->v.u.s[i] == nullptr)
                    {
                        snew(sel->v.u.s[i], 2);
                    }
                }
            }
        }
    }
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
static void evaluate_boolean_static_part(gmx_sel_evaluate_t*                data,
                                         const SelectionTreeElementPointer& sel,
                                         gmx_ana_index_t*                   g)
{
    /* Find the last static subexpression */
    SelectionTreeElementPointer child = sel->child;
    while (child->next && (child->next->cdata->flags & SEL_CDATA_STATIC))
    {
        child = child->next;
    }
    if (!(child->cdata->flags & SEL_CDATA_STATIC))
    {
        return;
    }

    /* Evalute the static part if there is more than one expression */
    if (child != sel->child)
    {
        SelectionTreeElementPointer next = child->next;
        child->next.reset();
        sel->cdata->evaluate(data, sel, g);
        /* Replace the subexpressions with the result */
        child = std::make_shared<SelectionTreeElement>(SEL_CONST, SelectionLocation::createEmpty());
        child->flags = SEL_FLAGSSET | SEL_SINGLEVAL | SEL_ALLOCVAL | SEL_ALLOCDATA;
        _gmx_selelem_set_vtype(child, GROUP_VALUE);
        child->evaluate = nullptr;
        _gmx_selvalue_reserve(&child->v, 1);
        gmx_ana_index_copy(child->v.u.g, sel->v.u.g, true);
        init_item_compilerdata(child);
        init_item_minmax_groups(child);
        child->cdata->flags &= ~SEL_CDATA_STATICEVAL;
        child->cdata->flags |= sel->cdata->flags & SEL_CDATA_STATICEVAL;
        child->next = next;
        // Frees the old static subexpressions.
        sel->child = child;
    }
    else if (child->evaluate)
    {
        child->evaluate(data, child, g);
    }
    /* Set the evaluation function for the constant element.
     * We never need to evaluate the element again during compilation,
     * but we may need to evaluate the static part again if the
     * expression is not an OR with a static evaluation group.
     * If we reach here with a NOT expression, the NOT expression
     * is also static, and will be made a constant later, so don't waste
     * time copying the group. */
    child->evaluate = nullptr;
    if (sel->u.boolt == BOOL_NOT || ((sel->cdata->flags & SEL_CDATA_STATICEVAL) && sel->u.boolt == BOOL_OR))
    {
        child->cdata->evaluate = nullptr;
    }
    else
    {
        child->cdata->evaluate = &_gmx_sel_evaluate_static;
        /* The cgrp has only been allocated if it originated from an
         * external index group.
         * If cgrp has been set in make_static(), it is not allocated,
         * and hence we can overwrite it safely. */
        if (child->u.cgrp.nalloc_index > 0)
        {
            gmx_ana_index_copy(&child->u.cgrp, child->v.u.g, false);
            gmx_ana_index_squeeze(&child->u.cgrp);
        }
        else
        {
            gmx_ana_index_copy(&child->u.cgrp, child->v.u.g, true);
        }
    }
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
static void evaluate_boolean_minmax_grps(const SelectionTreeElementPointer& sel,
                                         gmx_ana_index_t*                   g,
                                         gmx_ana_index_t*                   gmin,
                                         gmx_ana_index_t*                   gmax)
{
    SelectionTreeElementPointer child;

    switch (sel->u.boolt)
    {
        case BOOL_NOT:
            GMX_ASSERT(g != nullptr, "Need a valid group");
            gmx_ana_index_reserve(gmin, g->isize);
            gmx_ana_index_reserve(gmax, g->isize);
            gmx_ana_index_difference(gmax, g, sel->child->cdata->gmin);
            gmx_ana_index_difference(gmin, g, sel->child->cdata->gmax);
            break;

        case BOOL_AND:
            gmx_ana_index_copy(gmin, sel->child->cdata->gmin, true);
            gmx_ana_index_copy(gmax, sel->child->cdata->gmax, true);
            child = sel->child->next;
            while (child && gmax->isize > 0)
            {
                gmx_ana_index_intersection(gmin, gmin, child->cdata->gmin);
                gmx_ana_index_intersection(gmax, gmax, child->cdata->gmax);
                child = child->next;
            }
            /* Update the static part if other expressions limit it */
            if ((sel->child->cdata->flags & SEL_CDATA_STATIC) && sel->child->v.u.g->isize > gmax->isize)
            {
                gmx_ana_index_copy(sel->child->v.u.g, gmax, false);
                gmx_ana_index_squeeze(sel->child->v.u.g);
                if (sel->child->u.cgrp.isize > 0)
                {
                    gmx_ana_index_copy(&sel->child->u.cgrp, gmax, false);
                    gmx_ana_index_squeeze(&sel->child->u.cgrp);
                }
            }
            break;

        case BOOL_OR:
            /* We can assume here that the gmin of children do not overlap
             * because of the way _gmx_sel_evaluate_or() works. */
            GMX_ASSERT(g != nullptr, "Need a valid group");
            gmx_ana_index_reserve(gmin, g->isize);
            gmx_ana_index_reserve(gmax, g->isize);
            gmx_ana_index_copy(gmin, sel->child->cdata->gmin, false);
            gmx_ana_index_copy(gmax, sel->child->cdata->gmax, false);
            child = sel->child->next;
            while (child && gmin->isize < g->isize)
            {
                gmx_ana_index_merge(gmin, gmin, child->cdata->gmin);
                gmx_ana_index_union(gmax, gmax, child->cdata->gmax);
                child = child->next;
            }
            /* Update the static part if other expressions have static parts
             * that are not included. */
            if ((sel->child->cdata->flags & SEL_CDATA_STATIC) && sel->child->v.u.g->isize < gmin->isize)
            {
                GMX_RELEASE_ASSERT(sel->child->type == SEL_CONST,
                                   "The first child should have already been evaluated "
                                   "to a constant expression");
                gmx_ana_index_reserve(sel->child->v.u.g, gmin->isize);
                gmx_ana_index_copy(sel->child->v.u.g, gmin, false);
                if (sel->child->u.cgrp.nalloc_index > 0)
                {
                    gmx_ana_index_reserve(&sel->child->u.cgrp, gmin->isize);
                    gmx_ana_index_copy(&sel->child->u.cgrp, gmin, false);
                }
                else
                {
                    GMX_RELEASE_ASSERT(sel->child->u.cgrp.index == sel->child->v.u.g->index,
                                       "If not allocated, the static group should equal the value");
                    sel->child->u.cgrp.isize = sel->child->v.u.g->isize;
                }
            }
            break;

        case BOOL_XOR: /* Should not be reached */
            GMX_THROW(gmx::NotImplementedError("xor expressions not implemented"));
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
 * This function is used as the replacement for the
 * gmx::SelectionTreeElement::evaluate function pointer.
 * It does the single most complex task in the compiler: after all elements
 * have been processed, the \p gmin and \p gmax fields of \p t_compiler_data
 * have been properly initialized, enough memory has been allocated for
 * storing the value of each expression, and the static parts of the
 * expressions have been evaluated.
 * The above is exactly true only for elements other than subexpressions:
 * another pass is required for subexpressions that are referred to more than
 * once and whose evaluation group is not known in advance.
 */
static void analyze_static(gmx_sel_evaluate_t* data, const SelectionTreeElementPointer& sel, gmx_ana_index_t* g)
{
    bool bDoMinMax;

    if (sel->type != SEL_ROOT && g)
    {
        alloc_selection_data(sel, g->isize, false);
    }

    bDoMinMax = ((sel->cdata->flags & SEL_CDATA_DOMINMAX) != 0);
    if (sel->type != SEL_SUBEXPR && bDoMinMax)
    {
        gmx_ana_index_deinit(sel->cdata->gmin);
        gmx_ana_index_deinit(sel->cdata->gmax);
    }

    /* TODO: This switch is awfully long... */
    switch (sel->type)
    {
        case SEL_CONST: process_const(data, sel, g); break;

        case SEL_EXPRESSION:
        case SEL_MODIFIER:
        {
            const int isize = g ? g->isize : 0;
            _gmx_sel_evaluate_method_params(data, sel, g);
            init_method(sel, data->top, isize);
            if (!(sel->flags & SEL_DYNAMIC))
            {
                sel->cdata->evaluate(data, sel, g);
                if (sel->cdata->flags & SEL_CDATA_STATIC)
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
                    sel->cdata->evaluate(data, sel, g);
                }
                else
                {
                    if (sel->v.type != GROUP_VALUE && sel->v.type != POS_VALUE)
                    {
                        sel->v.nr = ((sel->flags & SEL_SINGLEVAL) ? 1 : isize);
                    }
                    // Clear the values for dynamic output to avoid
                    // uninitialized memory.
                    if (sel->v.type == REAL_VALUE)
                    {
                        for (int i = 0; i < sel->v.nr; ++i)
                        {
                            sel->v.u.r[i] = 0.0;
                        }
                    }
                }
                if (bDoMinMax && g)
                {
                    gmx_ana_index_copy(sel->cdata->gmax, g, true);
                }
            }
            break;
        }
        case SEL_BOOLEAN:
            if (!(sel->flags & SEL_DYNAMIC))
            {
                sel->cdata->evaluate(data, sel, g);
                if (sel->cdata->flags & SEL_CDATA_STATIC)
                {
                    make_static(sel);
                }
            }
            else
            {
                /* Evalute the static part if there is more than one expression */
                evaluate_boolean_static_part(data, sel, g);

                /* Evaluate the selection.
                 * If the type is boolean, we must explicitly handle the
                 * static part evaluated in evaluate_boolean_static_part()
                 * here because g may be larger. */
                if (sel->u.boolt == BOOL_AND && sel->child->type == SEL_CONST)
                {
                    sel->cdata->evaluate(data, sel, sel->child->v.u.g);
                }
                else
                {
                    sel->cdata->evaluate(data, sel, g);
                }

                /* Evaluate minimal and maximal selections */
                evaluate_boolean_minmax_grps(sel, g, sel->cdata->gmin, sel->cdata->gmax);
            }
            break;

        case SEL_ARITHMETIC:
            sel->cdata->evaluate(data, sel, g);
            if (!(sel->flags & SEL_DYNAMIC))
            {
                if (sel->cdata->flags & SEL_CDATA_STATIC)
                {
                    make_static(sel);
                }
            }
            else if (bDoMinMax)
            {
                gmx_ana_index_copy(sel->cdata->gmax, g, true);
            }
            break;

        case SEL_ROOT: sel->cdata->evaluate(data, sel, g); break;

        case SEL_SUBEXPR:
            if (((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                 && !(sel->cdata->flags & SEL_CDATA_STATICMULTIEVALSUBEXPR))
                || (sel->cdata->flags & SEL_CDATA_FULLEVAL))
            {
                sel->cdata->evaluate(data, sel, g);
                _gmx_selvalue_setstore(&sel->v, sel->child->v.u.ptr);
            }
            else if (sel->u.cgrp.isize == 0)
            {
                GMX_ASSERT(g, "group cannot be null");
                gmx_ana_index_reserve(&sel->u.cgrp, g->isize);
                sel->cdata->evaluate(data, sel, g);
                if (bDoMinMax)
                {
                    gmx_ana_index_copy(sel->cdata->gmin, sel->child->cdata->gmin, true);
                    gmx_ana_index_copy(sel->cdata->gmax, sel->child->cdata->gmax, true);
                }
            }
            else
            {
                int isize = gmx_ana_index_difference_size(g, &sel->u.cgrp);
                if (isize > 0)
                {
                    isize += sel->u.cgrp.isize;
                    gmx_ana_index_reserve(&sel->u.cgrp, isize);
                    alloc_selection_data(sel, isize, false);
                }
                sel->cdata->evaluate(data, sel, g);
                if (isize > 0 && bDoMinMax)
                {
                    gmx_ana_index_reserve(sel->cdata->gmin,
                                          sel->cdata->gmin->isize + sel->child->cdata->gmin->isize);
                    gmx_ana_index_reserve(sel->cdata->gmax,
                                          sel->cdata->gmax->isize + sel->child->cdata->gmax->isize);
                    gmx_ana_index_merge(sel->cdata->gmin, sel->cdata->gmin, sel->child->cdata->gmin);
                    gmx_ana_index_merge(sel->cdata->gmax, sel->cdata->gmax, sel->child->cdata->gmax);
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
                    gmx_ana_pos_copy(sel->v.u.p, sel->child->v.u.p, true);
                }
                else
                {
                    alloc_selection_data(sel, sel->child->cdata->gmax->isize, true);
                }
            }
            sel->cdata->evaluate(data, sel, g);
            if ((sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR) && (sel->child->child->flags & SEL_ALLOCVAL))
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
                    gmx_ana_index_copy(sel->cdata->gmin, sel->child->cdata->gmin, true);
                    gmx_ana_index_copy(sel->cdata->gmax, sel->child->cdata->gmax, true);
                }
                else
                {
                    gmx_ana_index_reserve(sel->cdata->gmin, min(g->isize, sel->child->cdata->gmin->isize));
                    gmx_ana_index_reserve(sel->cdata->gmax, min(g->isize, sel->child->cdata->gmax->isize));
                    gmx_ana_index_intersection(sel->cdata->gmin, sel->child->cdata->gmin, g);
                    gmx_ana_index_intersection(sel->cdata->gmax, sel->child->cdata->gmax, g);
                }
            }
            break;

        case SEL_GROUPREF: /* Should not be reached */
            GMX_THROW(gmx::APIError("Unresolved group reference in compilation"));
    }

    /* Update the minimal and maximal evaluation groups */
    if (bDoMinMax)
    {
        gmx_ana_index_squeeze(sel->cdata->gmin);
        gmx_ana_index_squeeze(sel->cdata->gmax);
    }

    /* Replace the result of the evaluation */
    /* This is not necessary for subexpressions or for boolean negations
     * because the evaluation function already has done it properly. */
    if (sel->v.type == GROUP_VALUE && (sel->flags & SEL_DYNAMIC) && sel->type != SEL_SUBEXPR
        && !(sel->type == SEL_BOOLEAN && sel->u.boolt == BOOL_NOT))
    {
        if (sel->cdata->flags & SEL_CDATA_EVALMAX)
        {
            gmx_ana_index_copy(sel->v.u.g, sel->cdata->gmax, false);
        }
        else
        {
            gmx_ana_index_copy(sel->v.u.g, sel->cdata->gmin, false);
        }
    }
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
static void init_root_item(const SelectionTreeElementPointer& root, gmx_ana_index_t* gall)
{
    const SelectionTreeElementPointer& expr = root->child;
    /* Subexpressions with non-static evaluation group should not be
     * evaluated by the root, and neither should be single-reference
     * subexpressions that don't evaluate for all atoms. */
    if (expr->type == SEL_SUBEXPR
        && (!(root->child->cdata->flags & SEL_CDATA_STATICEVAL)
            || ((root->child->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
                && !(root->child->cdata->flags & SEL_CDATA_FULLEVAL))))
    {
        root->evaluate = nullptr;
        if (root->cdata)
        {
            root->cdata->evaluate = nullptr;
        }
    }

    /* Set the evaluation group */
    if (root->evaluate)
    {
        /* Non-atom-valued non-group expressions don't care about the group, so
         * don't allocate any memory for it. */
        if ((expr->flags & SEL_VARNUMVAL) || ((expr->flags & SEL_SINGLEVAL) && expr->v.type != GROUP_VALUE))
        {
            gmx_ana_index_set(&root->u.cgrp, -1, nullptr, 0);
        }
        else if (expr->cdata->gmax->isize == gall->isize)
        {
            /* Save some memory by only referring to the global group. */
            gmx_ana_index_set(&root->u.cgrp, gall->isize, gall->index, 0);
        }
        else
        {
            gmx_ana_index_copy(&root->u.cgrp, expr->cdata->gmax, true);
        }
        /* For selections, store the maximum group for
         * gmx_ana_selcollection_evaluate_fin() as the value of the root
         * element (unused otherwise). */
        if (expr->type != SEL_SUBEXPR && expr->v.u.p->m.mapb.a != nullptr)
        {
            SelectionTreeElementPointer child = expr;

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
                gmx_ana_index_t g;
                gmx_ana_index_set(&g, expr->v.u.p->m.mapb.nra, expr->v.u.p->m.mapb.a, 0);
                _gmx_selelem_set_vtype(root, GROUP_VALUE);
                root->flags |= (SEL_ALLOCVAL | SEL_ALLOCDATA);
                _gmx_selvalue_reserve(&root->v, 1);
                gmx_ana_index_copy(root->v.u.g, &g, true);
            }
        }
    }
    else
    {
        gmx_ana_index_clear(&root->u.cgrp);
    }
}


/********************************************************************
 * REQUIRED ATOMS ANALYSIS
 ********************************************************************/

/*! \brief
 * Finds the highest atom index required to evaluate a selection subtree.
 *
 * \param[in]     sel           Root of the selection subtree to process.
 * \param[in,out] requiredAtoms The atoms required to evaluate the subtree.
 *      The existing group is only added to, so multiple calls with
 *      the same parameter will compute the union over several subtrees.
 *
 * For evaluation that starts from a \ref SEL_ROOT element with a fixed group,
 * children will never extend the evaluation group except for method parameter
 * evaluation (which have their own root element), so it is sufficient to check
 * the root.  However, children of \ref SEL_EXPRESSION elements (i.e., the
 * method parameters) may have been independently evaluated to a static group
 * that no longer has a separate root, so those need to be checked as well.
 *
 * Position calculations are not considered here, but are analyzed through the
 * position calculation collection in the main compilation method.
 */
static void init_required_atoms(const SelectionTreeElementPointer& sel, gmx_ana_index_t* requiredAtoms)
{
    // Process children.
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_required_atoms(child, requiredAtoms);
            child = child->next;
        }
    }

    if (sel->type == SEL_ROOT || (sel->type == SEL_CONST && sel->v.type == GROUP_VALUE))
    {
        if (sel->u.cgrp.isize > 0)
        {
            gmx_ana_index_union_unsorted(requiredAtoms, requiredAtoms, &sel->u.cgrp);
        }
    }
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
static void postprocess_item_subexpressions(const SelectionTreeElementPointer& sel)
{
    GMX_ASSERT(!(sel->child == nullptr && (sel->type == SEL_SUBEXPRREF || sel->type == SEL_SUBEXPR)),
               "Subexpression elements should always have a child element");

    /* Process children. */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            postprocess_item_subexpressions(child);
            child = child->next;
        }
    }

    /* Replace the evaluation function of statically evaluated subexpressions
     * for which the static group was not known in advance. */
    if (sel->type == SEL_SUBEXPR && sel->cdata->refcount > 1
        && (sel->cdata->flags & SEL_CDATA_STATICEVAL) && !(sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        /* We need to free memory allocated for the group, because it is no
         * longer needed (and would be lost on next call to the evaluation
         * function). */
        gmx_ana_index_deinit(&sel->u.cgrp);

        sel->evaluate        = &_gmx_sel_evaluate_subexpr_staticeval;
        sel->cdata->evaluate = sel->evaluate;

        sel->child->freeValues();
        sel->child->mempool = nullptr;
        _gmx_selvalue_setstore(&sel->child->v, sel->v.u.ptr);
        sel->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
    }

    /* Adjust memory allocation flags for subexpressions that are used only
     * once.  This is not strictly necessary, but we do it to have the memory
     * managed consistently for all types of subexpressions. */
    if (sel->type == SEL_SUBEXPRREF && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR))
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
    if (sel->type == SEL_SUBEXPR && !(sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
        && (sel->cdata->flags & SEL_CDATA_FULLEVAL))
    {
        sel->flags |= SEL_ALLOCVAL;
        sel->flags |= (sel->child->flags & SEL_ALLOCDATA);
        sel->v.nalloc = sel->child->v.nalloc;
        sel->child->flags &= ~(SEL_ALLOCVAL | SEL_ALLOCDATA);
        sel->child->v.nalloc = -1;
    }

    /* For static subexpressions with a dynamic evaluation group, there is
     * no need to evaluate them again, as the SEL_SUBEXPRREF takes care of
     * everything during evaluation. */
    if (sel->type == SEL_SUBEXPR && (sel->cdata->flags & SEL_CDATA_SIMPLESUBEXPR)
        && (sel->cdata->flags & SEL_CDATA_STATICMULTIEVALSUBEXPR))
    {
        sel->evaluate        = nullptr;
        sel->cdata->evaluate = nullptr;
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
static void init_item_comg(const SelectionTreeElementPointer&  sel,
                           gmx::PositionCalculationCollection* pcc,
                           e_poscalc_t                         type,
                           int                                 flags)
{
    /* Initialize COM calculation for dynamic selections now that we know the maximal evaluation group */
    if (sel->type == SEL_EXPRESSION && sel->u.expr.method && sel->u.expr.method->pupdate)
    {
        if (!sel->u.expr.method->update || type != POS_ATOM)
        {
            /* Create a default calculation if one does not yet exist */
            int cflags = 0;
            if (!(sel->cdata->flags & SEL_CDATA_STATICEVAL))
            {
                cflags |= POS_DYNAMIC;
            }
            if (!sel->u.expr.pc)
            {
                cflags |= flags;
                sel->u.expr.pc = pcc->createCalculation(type, cflags);
            }
            else
            {
                gmx_ana_poscalc_set_flags(sel->u.expr.pc, cflags);
            }
            gmx_ana_poscalc_set_maxindex(sel->u.expr.pc, sel->cdata->gmax);
            sel->u.expr.pos = new gmx_ana_pos_t();
            gmx_ana_poscalc_init_pos(sel->u.expr.pc, sel->u.expr.pos);
        }
    }

    /* Call recursively for all children unless the children have already been processed */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            init_item_comg(child, pcc, type, flags);
            child = child->next;
        }
    }
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
static void free_item_compilerdata(const SelectionTreeElementPointer& sel)
{
    /* Free compilation data */
    sel->freeCompilerData();

    /* Call recursively for all children unless the children have already been processed */
    if (sel->type != SEL_SUBEXPRREF)
    {
        SelectionTreeElementPointer child = sel->child;
        while (child)
        {
            free_item_compilerdata(child);
            child = child->next;
        }
    }
}


/********************************************************************
 * MAIN COMPILATION FUNCTION
 ********************************************************************/

namespace gmx
{

void compileSelection(SelectionCollection* coll)
{
    gmx_ana_selcollection_t*    sc = &coll->impl_->sc_;
    gmx_sel_evaluate_t          evaldata;
    SelectionTreeElementPointer item;
    e_poscalc_t                 post;
    size_t                      i;
    int                         flags;
    bool bDebug = (coll->impl_->debugLevel_ == SelectionCollection::Impl::DebugLevel::Compiled
                   || coll->impl_->debugLevel_ == SelectionCollection::Impl::DebugLevel::Full);

    /* FIXME: Clean up the collection on exceptions */

    sc->mempool = _gmx_sel_mempool_create();
    _gmx_sel_evaluate_init(&evaldata, sc->mempool, &sc->gall, sc->top, nullptr, nullptr);

    /* Clear the symbol table because it is not possible to parse anything
     * after compilation, and variable references in the symbol table can
     * also mess up the compilation and/or become invalid.
     */
    coll->impl_->clearSymbolTable();

    /* Loop through selections and initialize position keyword defaults if no
     * other value has been provided.
     */
    for (i = 0; i < sc->sel.size(); ++i)
    {
        gmx::internal::SelectionData& sel = *sc->sel[i];
        init_pos_keyword_defaults(
                &sel.rootElement(), coll->impl_->spost_.c_str(), coll->impl_->rpost_.c_str(), &sel);
    }

    /* Remove any unused variables. */
    sc->root = remove_unused_subexpressions(sc->root);
    /* Extract subexpressions into separate roots */
    sc->root = extract_subexpressions(sc->root);

    /* Initialize the evaluation callbacks and process the tree structure
     * to conform to the expectations of the callback functions. */
    /* Also, initialize and allocate the compiler data structure */
    // TODO: Processing the tree in reverse root order would be better,
    // as it would make dependency handling easier (all subexpression
    // references would be processed before the actual subexpression) and
    // could remove the need for most of these extra loops.
    item = sc->root;
    while (item)
    {
        /* Process boolean and arithmetic expressions. */
        optimize_boolean_expressions(item);
        reorder_boolean_static_children(item);
        optimize_arithmetic_expressions(item);
        /* Initialize the compiler data */
        init_item_compilerdata(item);
        item = item->next;
    }
    // Initialize the static evaluation compiler flags.
    // Requires the FULLEVAL compiler flag for the whole tree.
    item = sc->root;
    while (item)
    {
        init_item_staticeval(item);
        init_item_subexpr_refcount(item);
        item = item->next;
    }
    /* Initialize subexpression flags.
     * Requires compiler flags for the full tree. */
    item = sc->root;
    while (item)
    {
        init_item_subexpr_flags(item);
        item = item->next;
    }
    /* Initialize evaluation.
     * Requires subexpression flags. */
    item = sc->root;
    while (item)
    {
        init_item_evalfunc(item);
        setup_memory_pooling(item, sc->mempool);
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

    if (bDebug)
    {
        fprintf(stderr, "\nTree after initial compiler processing:\n");
        coll->printTree(stderr, false);
    }

    /* Evaluate all static parts of the selection and analyze the tree
     * to allocate enough memory to store the value of each dynamic subtree. */
    item = sc->root;
    while (item)
    {
        if (item->child->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
        {
            mark_subexpr_dynamic(item->child, true);
        }
        set_evaluation_function(item, &analyze_static);
        item->evaluate(&evaldata, item, nullptr);
        item = item->next;
    }

    /* At this point, static subexpressions no longer have references to them,
     * so they can be removed. */
    sc->root = remove_unused_subexpressions(sc->root);
    // Update the reference counts for consistency (only used for the
    // debugging output below).
    item = sc->root;
    while (item)
    {
        init_item_subexpr_refcount(item);
        item = item->next;
    }

    if (bDebug)
    {
        fprintf(stderr, "\nTree after first analysis pass:\n");
        coll->printTree(stderr, false);
    }

    /* Do a second pass to evaluate static parts of common subexpressions */
    item = sc->root;
    while (item)
    {
        if (item->child->cdata->flags & SEL_CDATA_COMMONSUBEXPR)
        {
            bool bMinMax = (item->child->cdata->flags & SEL_CDATA_DOMINMAX) != 0;

            mark_subexpr_dynamic(item->child, false);
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
            item->child->cdata->flags &= ~SEL_CDATA_DOMINMAX;
            item->evaluate(&evaldata, item->child, item->child->cdata->gmax);
            if (bMinMax)
            {
                item->child->cdata->flags |= SEL_CDATA_DOMINMAX;
            }
        }
        item = item->next;
    }

    /* We need a yet another pass of subexpression removal to remove static
     * subexpressions referred to by common dynamic subexpressions. */
    sc->root = remove_unused_subexpressions(sc->root);
    // Update the reference counts, used by postprocess_item_subexpressions().
    item = sc->root;
    while (item)
    {
        init_item_subexpr_refcount(item);
        item = item->next;
    }

    if (bDebug)
    {
        fprintf(stderr, "\nTree after second analysis pass:\n");
        coll->printTree(stderr, false);
    }

    // Initialize evaluation groups, position calculations for methods, perform
    // some final optimization, and free the memory allocated for the
    // compilation.
    /* By default, use whole residues/molecules. */
    flags = POS_COMPLWHOLE;
    PositionCalculationCollection::typeFromEnum(coll->impl_->rpost_.c_str(), &post, &flags);
    item = sc->root;
    while (item)
    {
        init_root_item(item, &sc->gall);
        postprocess_item_subexpressions(item);
        init_item_comg(item, &sc->pcc, post, flags);
        free_item_compilerdata(item);
        item = item->next;
    }

    // Compute the atoms required for evaluating the selections.
    gmx_ana_index_clear(&coll->impl_->requiredAtoms_);
    gmx_ana_index_reserve(&coll->impl_->requiredAtoms_, sc->gall.isize);
    sc->pcc.getRequiredAtoms(&coll->impl_->requiredAtoms_);
    item = sc->root;
    while (item)
    {
        init_required_atoms(item, &coll->impl_->requiredAtoms_);
        item = item->next;
    }
    gmx_ana_index_squeeze(&coll->impl_->requiredAtoms_);

    /* Allocate memory for the evaluation memory pool. */
    _gmx_sel_mempool_reserve(sc->mempool, 0);

    /* Finish up by calculating total masses and charges. */
    for (i = 0; i < sc->sel.size(); ++i)
    {
        sc->sel[i]->initializeMassesAndCharges(sc->top);
    }
}

} // namespace gmx
