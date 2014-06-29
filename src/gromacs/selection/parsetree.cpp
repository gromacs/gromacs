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
 * \brief
 * Implements functions in parsetree.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
/*! \internal
 * \page page_module_selection_parser Selection parsing
 *
 * The selection parser is implemented in the following files:
 *  - scanner.l:
 *    Tokenizer implemented using Flex, splits the input into tokens
 *    (scanner.c and scanner_flex.h are generated from this file).
 *  - scanner.h, scanner_internal.h, scanner_internal.cpp:
 *    Helper functions for scanner.l and for interfacing between
 *    scanner.l and parser.y. Functions in scanner_internal.h are only
 *    used from scanner.l, while scanner.h is used from the parser.
 *  - symrec.h, symrec.cpp:
 *    Functions used by the tokenizer to handle the symbol table, i.e.,
 *    the recognized keywords. Some basic keywords are hardcoded into
 *    scanner.l, but all method and variable references go through the
 *    symbol table, as do position evaluation keywords.
 *  - parser.y:
 *    Semantic rules for parsing the grammar
 *    (parser.cpp and parser.h are generated from this file by Bison).
 *  - parsetree.h, parsetree.cpp:
 *    Functions called from actions in parser.y to construct the
 *    evaluation elements corresponding to different grammar elements.
 *  - params.c:
 *    Defines a function that processes the parameters of selection
 *    methods and initializes the children of the method element.
 *  - selectioncollection.h, selectioncollection.cpp:
 *    These files define the high-level public interface to the parser
 *    through SelectionCollection::parseFromStdin(),
 *    SelectionCollection::parseFromFile() and
 *    SelectionCollection::parseFromString().
 *
 * The basic control flow in the parser is as follows: when a parser function
 * in SelectionCollection gets called, it performs some
 * initialization, and then calls the _gmx_sel_yyparse() function generated
 * by Bison. This function then calls _gmx_sel_yylex() to repeatedly read
 * tokens from the input (more complex tasks related to token recognition
 * and bookkeeping are done by functions in scanner_internal.cpp) and uses the
 * grammar rules to decide what to do with them. Whenever a grammar rule
 * matches, a corresponding function in parsetree.cpp is called to construct
 * either a temporary representation for the object or a
 * gmx::SelectionTreeElement object
 * (some simple rules are handled internally in parser.y).
 * When a complete selection has been parsed, the functions in parsetree.cpp
 * also take care of updating the ::gmx_ana_selcollection_t structure
 * appropriately.
 *
 * The rest of this page describes the resulting gmx::SelectionTreeElement
 * object tree.
 * Before the selections can be evaluated, this tree needs to be passed to
 * the selection compiler, which is described on a separate page:
 * \ref page_module_selection_compiler
 *
 *
 * \section selparser_tree Element tree constructed by the parser
 *
 * The parser initializes the following fields in all selection elements:
 * gmx::SelectionTreeElement::name, gmx::SelectionTreeElement::type,
 * gmx::SelectionTreeElement::v\c .type,
 * gmx::SelectionTreeElement::flags, gmx::SelectionTreeElement::child, and
 * gmx::SelectionTreeElement::next.
 * Some other fields are also initialized for particular element types as
 * discussed below.
 * Fields that are not initialized are set to zero, NULL, or other similar
 * value.
 *
 *
 * \subsection selparser_tree_root Root elements
 *
 * The parser creates a \ref SEL_ROOT selection element for each variable
 * assignment and each selection. However, there are two exceptions that do
 * not result in a \ref SEL_ROOT element (in these cases, only the symbol
 * table is modified):
 *  - Variable assignments that assign a variable to another variable.
 *  - Variable assignments that assign a non-group constant.
 *  .
 * The \ref SEL_ROOT elements are linked together in a chain in the same order
 * as in the input.
 *
 * The children of the \ref SEL_ROOT elements can be used to distinguish
 * the two types of root elements from each other:
 *  - For variable assignments, the first and only child is always
 *    a \ref SEL_SUBEXPR element.
 *  - For selections, the first child is a \ref SEL_EXPRESSION or a
 *    \ref SEL_MODIFIER element that evaluates the final positions (if the
 *    selection defines a constant position, the child is a \ref SEL_CONST).
 *    The rest of the children are \ref SEL_MODIFIER elements with
 *    \ref NO_VALUE, in the order given by the user.
 *  .
 * The name of the selection/variable is stored in
 * gmx::SelectionTreeElement::cgrp\c .name.
 * It is set to either the name provided by the user or the selection string
 * for selections not explicitly named by the user.
 * \ref SEL_ROOT or \ref SEL_SUBEXPR elements do not appear anywhere else.
 *
 *
 * \subsection selparser_tree_const Constant elements
 *
 * \ref SEL_CONST elements are created for every constant that is required
 * for later evaluation.
 * Currently, \ref SEL_CONST elements can be present for
 *  - selections that consist of a constant position,
 *  - \ref GROUP_VALUE method parameters if provided using external index
 *    groups,
 *  .
 * For group-valued elements, the value is stored in
 * gmx::SelectionTreeElement::cgrp; other types of values are stored in
 * gmx::SelectionTreeElement::v.
 * Constants that appear as parameters for selection methods are not present
 * in the selection tree unless they have \ref GROUP_VALUE.
 * \ref SEL_CONST elements have no children.
 *
 *
 * \subsection selparser_tree_method Method evaluation elements
 *
 * \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements are treated very
 * similarly. The \c gmx_ana_selmethod_t structure corresponding to the
 * evaluation method is in gmx::SelectionTreeElement::method, and the method
 * data in gmx::SelectionTreeElement::mdata has been allocated using
 * sel_datafunc().
 * If a non-standard reference position type was set,
 * gmx::SelectionTreeElement::pc has also been created, but only the type has
 * been set.
 * All children of these elements are of the type \ref SEL_SUBEXPRREF, and
 * each describes a selection that needs to be evaluated to obtain a value
 * for one parameter of the method.
 * No children are present for parameters that were given a constant
 * non-\ref GROUP_VALUE value.
 * The children are sorted in the order in which the parameters appear in the
 * \ref gmx_ana_selmethod_t structure.
 *
 * In addition to actual selection keywords, \ref SEL_EXPRESSION elements
 * are used internally to implement numerical comparisons (e.g., "x < 5")
 * and keyword matching (e.g., "resnr 1 to 3" or "name CA").
 *
 *
 * \subsection selparser_tree_subexpr Subexpression elements
 *
 * \ref SEL_SUBEXPR elements only appear for variables, as described above.
 * gmx::SelectionTreeElement::name points to the name of the variable (from the
 * \ref SEL_ROOT element).
 * The element always has exactly one child, which represents the value of
 * the variable.
 *
 * \ref SEL_SUBEXPRREF elements are used for two purposes:
 *  - Variable references that need to be evaluated (i.e., there is a
 *    \ref SEL_SUBEXPR element for the variable) are represented using
 *    \ref SEL_SUBEXPRREF elements.
 *    In this case, gmx::SelectionTreeElement::param is NULL, and the first and
 *    only child of the element is the \ref SEL_SUBEXPR element of the
 *    variable.
 *    Such references can appear anywhere where the variable value
 *    (the child of the \ref SEL_SUBEXPR element) would be valid.
 *  - Children of \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements are
 *    always of this type. For these elements, gmx::SelectionTreeElement::param
 *    is initialized to point to the parameter that receives the value from
 *    the expression.
 *    Each such element has exactly one child, which can be of any type;
 *    the \ref SEL_SUBEXPR element of a variable is used if the value comes
 *    from a variable, otherwise the child type is not \ref SEL_SUBEXPR.
 *
 *
 * \subsection selparser_tree_bool Boolean elements
 *
 * One \ref SEL_BOOLEAN element is created for each boolean keyword in the
 * input, and the tree structure represents the evaluation order.
 * The gmx::SelectionTreeElement::boolt type gives the type of the operation.
 * Each element has exactly two children (one for \ref BOOL_NOT elements),
 * which are in the order given in the input.
 * The children always have \ref GROUP_VALUE, but different element types
 * are possible.
 *
 *
 * \subsection selparser_tree_arith Arithmetic elements
 *
 * One \ref SEL_ARITHMETIC element is created for each arithmetic operation in
 * the input, and the tree structure represents the evaluation order.
 * The gmx::SelectionTreeElement::optype type gives the name of the operation.
 * Each element has exactly two children (one for unary negation elements),
 * which are in the order given in the input.
 */
#include <stdio.h>
#include <stdarg.h>

#include <boost/exception_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "gromacs/fileio/futil.h"
#include "gromacs/selection/poscalc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selmethod.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/messagestringcollector.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "keywords.h"
#include "parsetree.h"
#include "selectioncollection-impl.h"
#include "selelem.h"
#include "symrec.h"

#include "scanner.h"

using gmx::SelectionParserValue;
using gmx::SelectionParserValueList;
using gmx::SelectionParserValueListPointer;
using gmx::SelectionParserParameter;
using gmx::SelectionParserParameterList;
using gmx::SelectionParserParameterListPointer;
using gmx::SelectionParserValue;
using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;

void
_gmx_selparser_error(yyscan_t scanner, const char *fmt, ...)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    // FIXME: Use an arbitrary length buffer.
    char    buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, 1024, fmt, ap);
    va_end(ap);
    errors->append(buf);
}

bool
_gmx_selparser_handle_exception(yyscan_t scanner, const std::exception &ex)
{
    if (dynamic_cast<const gmx::UserInputError *>(&ex) != NULL)
    {
        // TODO: Consider whether also the non-interactive parser should
        // postpone the exception such that the whole selection can be added as
        // context.
        if (_gmx_sel_is_lexer_interactive(scanner))
        {
            // TODO: Handle exceptions that printing the message may produce.
            gmx::formatExceptionMessageToFile(stderr, ex);
            return true;
        }
    }
    _gmx_sel_lexer_set_exception(scanner, boost::current_exception());
    return false;
}

namespace gmx
{

/********************************************************************
 * SelectionParserValue
 */

SelectionParserValue::SelectionParserValue(e_selvalue_t type)
    : type(type)
{
    memset(&u, 0, sizeof(u));
}

SelectionParserValue::SelectionParserValue(
        const SelectionTreeElementPointer &expr)
    : type(expr->v.type), expr(expr)
{
    memset(&u, 0, sizeof(u));
}

/********************************************************************
 * SelectionParserParameter
 */

SelectionParserParameter::SelectionParserParameter(
        const char                     *name,
        SelectionParserValueListPointer values)
    : name_(name != NULL ? name : ""),
      values_(values ? move(values)
              : SelectionParserValueListPointer(new SelectionParserValueList))
{
}

} // namespace gmx

/*!
 * \param[in,out] sel  Root of the selection element tree to initialize.
 *
 * Propagates the \ref SEL_DYNAMIC flag from the children of \p sel to \p sel
 * (if any child of \p sel is dynamic, \p sel is also marked as such).
 * The \ref SEL_DYNAMIC flag is also set for \ref SEL_EXPRESSION elements with
 * a dynamic method.
 * Also, sets one of the \ref SEL_SINGLEVAL, \ref SEL_ATOMVAL, or
 * \ref SEL_VARNUMVAL flags, either based on the children or on the type of
 * the selection method.
 * If the types of the children conflict, an error is returned.
 *
 * The flags of the children of \p sel are also updated if not done earlier.
 * The flags are initialized only once for any element; if \ref SEL_FLAGSSET
 * is set for an element, the function returns immediately, and the recursive
 * operation does not descend beyond such elements.
 */
void
_gmx_selelem_update_flags(const gmx::SelectionTreeElementPointer &sel)
{
    bool                bUseChildType = false;
    bool                bOnlySingleChildren;

    /* Return if the flags have already been set */
    if (sel->flags & SEL_FLAGSSET)
    {
        return;
    }
    /* Set the flags based on the current element type */
    switch (sel->type)
    {
        case SEL_CONST:
        case SEL_GROUPREF:
            sel->flags   |= SEL_SINGLEVAL;
            bUseChildType = false;
            break;

        case SEL_EXPRESSION:
            if (sel->u.expr.method->flags & SMETH_DYNAMIC)
            {
                sel->flags |= SEL_DYNAMIC;
            }
            if (sel->u.expr.method->flags & SMETH_SINGLEVAL)
            {
                sel->flags |= SEL_SINGLEVAL;
            }
            else if (sel->u.expr.method->flags & SMETH_VARNUMVAL)
            {
                sel->flags |= SEL_VARNUMVAL;
            }
            else
            {
                sel->flags |= SEL_ATOMVAL;
            }
            bUseChildType = false;
            break;

        case SEL_ARITHMETIC:
            sel->flags   |= SEL_ATOMVAL;
            bUseChildType = false;
            break;

        case SEL_MODIFIER:
            if (sel->v.type != NO_VALUE)
            {
                sel->flags |= SEL_VARNUMVAL;
            }
            bUseChildType = false;
            break;

        case SEL_ROOT:
            bUseChildType = false;
            break;

        case SEL_BOOLEAN:
        case SEL_SUBEXPR:
        case SEL_SUBEXPRREF:
            bUseChildType = true;
            break;
    }
    /* Loop through children to propagate their flags upwards */
    bOnlySingleChildren = true;
    SelectionTreeElementPointer child = sel->child;
    while (child)
    {
        /* Update the child */
        _gmx_selelem_update_flags(child);
        /* Propagate the dynamic and unsorted flags */
        sel->flags |= (child->flags & (SEL_DYNAMIC | SEL_UNSORTED));
        /* Propagate the type flag if necessary and check for problems */
        if (bUseChildType)
        {
            if ((sel->flags & SEL_VALTYPEMASK)
                && !(sel->flags & child->flags & SEL_VALTYPEMASK))
            {
                // TODO: Recollect when this is triggered, and whether the type
                // is appropriate.
                GMX_THROW(gmx::InvalidInputError("Invalid combination of selection expressions"));
            }
            sel->flags |= (child->flags & SEL_VALTYPEMASK);
        }
        if (!(child->flags & SEL_SINGLEVAL))
        {
            bOnlySingleChildren = false;
        }

        child = child->next;
    }
    /* For arithmetic expressions consisting only of single values,
     * the result is also a single value. */
    if (sel->type == SEL_ARITHMETIC && bOnlySingleChildren)
    {
        sel->flags = (sel->flags & ~SEL_VALTYPEMASK) | SEL_SINGLEVAL;
    }
    /* For root elements, the type should be propagated here, after the
     * children have been updated. */
    if (sel->type == SEL_ROOT)
    {
        GMX_ASSERT(sel->child, "Root elements should always have a child");
        sel->flags |= (sel->child->flags & SEL_VALTYPEMASK);
    }
    /* Mark that the flags are set */
    sel->flags |= SEL_FLAGSSET;
}

/*!
 * \param[in,out] sel    Selection element to initialize.
 * \param[in]     scanner Scanner data structure.
 *
 * A deep copy of the parameters is made to allow several
 * expressions with the same method to coexist peacefully.
 * Calls sel_datafunc() if one is specified for the method.
 */
void
_gmx_selelem_init_method_params(const gmx::SelectionTreeElementPointer &sel,
                                yyscan_t                                scanner)
{
    int                 nparams;
    gmx_ana_selparam_t *orgparam;
    gmx_ana_selparam_t *param;
    int                 i;
    void               *mdata;

    nparams   = sel->u.expr.method->nparams;
    orgparam  = sel->u.expr.method->param;
    snew(param, nparams);
    memcpy(param, orgparam, nparams*sizeof(gmx_ana_selparam_t));
    for (i = 0; i < nparams; ++i)
    {
        param[i].flags &= ~SPAR_SET;
        _gmx_selvalue_setstore(&param[i].val, NULL);
        if (param[i].flags & SPAR_VARNUM)
        {
            param[i].val.nr = -1;
        }
        /* Duplicate the enum value array if it is given statically */
        if ((param[i].flags & SPAR_ENUMVAL) && orgparam[i].val.u.ptr != NULL)
        {
            int n;

            /* Count the values */
            n = 1;
            while (orgparam[i].val.u.s[n] != NULL)
            {
                ++n;
            }
            _gmx_selvalue_reserve(&param[i].val, n+1);
            memcpy(param[i].val.u.s, orgparam[i].val.u.s,
                   (n+1)*sizeof(param[i].val.u.s[0]));
        }
    }
    mdata = NULL;
    if (sel->u.expr.method->init_data)
    {
        mdata = sel->u.expr.method->init_data(nparams, param);
    }
    if (sel->u.expr.method->set_poscoll)
    {
        gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);

        sel->u.expr.method->set_poscoll(&sc->pcc, mdata);
    }
    /* Store the values */
    sel->u.expr.method->param = param;
    sel->u.expr.mdata         = mdata;
}

/*!
 * \param[in,out] sel    Selection element to initialize.
 * \param[in]     method Selection method to set.
 * \param[in]     scanner Scanner data structure.
 *
 * Makes a copy of \p method and stores it in \p sel->u.expr.method,
 * and calls _gmx_selelem_init_method_params();
 */
void
_gmx_selelem_set_method(const gmx::SelectionTreeElementPointer &sel,
                        gmx_ana_selmethod_t                    *method,
                        yyscan_t                                scanner)
{
    _gmx_selelem_set_vtype(sel, method->type);
    sel->setName(method->name);
    snew(sel->u.expr.method, 1);
    memcpy(sel->u.expr.method, method, sizeof(gmx_ana_selmethod_t));
    _gmx_selelem_init_method_params(sel, scanner);
}

/*! \brief
 * Initializes the reference position calculation for a \ref SEL_EXPRESSION
 * element.
 *
 * \param[in,out] pcc    Position calculation collection to use.
 * \param[in,out] sel    Selection element to initialize.
 * \param[in]     rpost  Reference position type to use (NULL = default).
 * \param[in]     scanner Scanner data structure.
 * \returns       0 on success, a non-zero error code on error.
 */
static void
set_refpos_type(gmx::PositionCalculationCollection *pcc,
                const SelectionTreeElementPointer &sel,
                const char *rpost, yyscan_t scanner)
{
    if (!rpost)
    {
        return;
    }

    if (sel->u.expr.method->pupdate)
    {
        /* By default, use whole residues/molecules. */
        sel->u.expr.pc
            = pcc->createCalculationFromEnum(rpost, POS_COMPLWHOLE);
    }
    else
    {
        // TODO: Should this be treated as a real error?
        _gmx_selparser_error(scanner, "modifier '%s' is not applicable for '%s'",
                             rpost, sel->u.expr.method->name);
    }
}

gmx::SelectionTreeElementPointer
_gmx_sel_init_arithmetic(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         char op, yyscan_t /* scanner */)
{
    SelectionTreeElementPointer sel(new SelectionTreeElement(SEL_ARITHMETIC));
    sel->v.type        = REAL_VALUE;
    switch (op)
    {
        case '+': sel->u.arith.type = ARITH_PLUS; break;
        case '-': sel->u.arith.type = (right ? ARITH_MINUS : ARITH_NEG); break;
        case '*': sel->u.arith.type = ARITH_MULT; break;
        case '/': sel->u.arith.type = ARITH_DIV;  break;
        case '^': sel->u.arith.type = ARITH_EXP;  break;
    }
    char               buf[2];
    buf[0] = op;
    buf[1] = 0;
    sel->setName(buf);
    sel->u.arith.opstr = strdup(buf);
    sel->child         = left;
    sel->child->next   = right;
    return sel;
}

/*!
 * \param[in]  left   Selection element for the left hand side.
 * \param[in]  right  Selection element for the right hand side.
 * \param[in]  cmpop  String representation of the comparison operator.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * comparison expressions.
 */
SelectionTreeElementPointer
_gmx_sel_init_comparison(const gmx::SelectionTreeElementPointer &left,
                         const gmx::SelectionTreeElementPointer &right,
                         const char *cmpop, yyscan_t scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    gmx::MessageStringContext    context(errors, "In comparison initialization");

    SelectionTreeElementPointer  sel(new SelectionTreeElement(SEL_EXPRESSION));
    _gmx_selelem_set_method(sel, &sm_compare, scanner);

    SelectionParserParameterList params;
    const char                  *name;
    // Create the parameter for the left expression.
    name  = left->v.type == INT_VALUE ? "int1" : "real1";
    params.push_back(SelectionParserParameter::createFromExpression(name, left));
    // Create the parameter for the right expression.
    name  = right->v.type == INT_VALUE ? "int2" : "real2";
    params.push_back(SelectionParserParameter::createFromExpression(name, right));
    // Create the parameter for the operator.
    params.push_back(
            SelectionParserParameter::create(
                    "op", SelectionParserValue::createString(cmpop)));
    if (!_gmx_sel_parse_params(params, sel->u.expr.method->nparams,
                               sel->u.expr.method->param, sel, scanner))
    {
        return SelectionTreeElementPointer();
    }

    return sel;
}

/*! \brief
 * Implementation method for keyword expression creation.
 *
 * \param[in]  method Method to use.
 * \param[in]  matchType String matching type (only used if \p method is
 *      a string keyword and \p args is not empty.
 * \param[in]  args   Pointer to the first argument.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * selection methods that do not take parameters.
 */
static SelectionTreeElementPointer
init_keyword_internal(gmx_ana_selmethod_t *method,
                      gmx::SelectionStringMatchType matchType,
                      SelectionParserValueListPointer args,
                      const char *rpost, yyscan_t scanner)
{
    gmx_ana_selcollection_t     *sc = _gmx_sel_lexer_selcollection(scanner);

    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    char  buf[128];
    sprintf(buf, "In keyword '%s'", method->name);
    gmx::MessageStringContext  context(errors, buf);

    if (method->nparams > 0)
    {
        // TODO: Would assert be better?
        GMX_THROW(gmx::InternalError(
                          "Keyword initialization called with non-keyword method"));
    }

    SelectionTreeElementPointer root(new SelectionTreeElement(SEL_EXPRESSION));
    SelectionTreeElementPointer child = root;
    _gmx_selelem_set_method(child, method, scanner);

    /* Initialize the evaluation of keyword matching if values are provided */
    if (args)
    {
        gmx_ana_selmethod_t *kwmethod;
        switch (method->type)
        {
            case INT_VALUE:  kwmethod = &sm_keyword_int;  break;
            case REAL_VALUE: kwmethod = &sm_keyword_real; break;
            case STR_VALUE:  kwmethod = &sm_keyword_str;  break;
            default:
                GMX_THROW(gmx::InternalError(
                                  "Unknown type for keyword selection"));
        }
        /* Initialize the selection element */
        root.reset(new SelectionTreeElement(SEL_EXPRESSION));
        _gmx_selelem_set_method(root, kwmethod, scanner);
        if (method->type == STR_VALUE)
        {
            _gmx_selelem_set_kwstr_match_type(root, matchType);
        }
        SelectionParserParameterList params;
        params.push_back(
                SelectionParserParameter::createFromExpression(NULL, child));
        params.push_back(SelectionParserParameter::create(NULL, move(args)));
        if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                                   root->u.expr.method->param, root, scanner))
        {
            return SelectionTreeElementPointer();
        }
    }
    set_refpos_type(&sc->pcc, child, rpost, scanner);

    return root;
}

/*!
 * \param[in]  method Method to use.
 * \param[in]  args   Pointer to the first argument.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * selection methods that do not take parameters.
 */
SelectionTreeElementPointer
_gmx_sel_init_keyword(gmx_ana_selmethod_t *method,
                      gmx::SelectionParserValueListPointer args,
                      const char *rpost, yyscan_t scanner)
{
    return init_keyword_internal(method, gmx::eStringMatchType_Auto, move(args),
                                 rpost, scanner);
}

/*!
 * \param[in]  method    Method to use.
 * \param[in]  matchType String matching type.
 * \param[in]  args      Pointer to the first argument.
 * \param[in]  rpost     Reference position type to use (NULL = default).
 * \param[in]  scanner   Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * keyword string matching.
 */
SelectionTreeElementPointer
_gmx_sel_init_keyword_strmatch(gmx_ana_selmethod_t *method,
                               gmx::SelectionStringMatchType matchType,
                               gmx::SelectionParserValueListPointer args,
                               const char *rpost, yyscan_t scanner)
{
    GMX_RELEASE_ASSERT(method->type == STR_VALUE,
                       "String keyword method called for a non-string-valued method");
    GMX_RELEASE_ASSERT(args && !args->empty(),
                       "String keyword matching method called without any values");
    return init_keyword_internal(method, matchType, move(args), rpost, scanner);
}

/*!
 * \param[in]  method Method to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * selection methods that take parameters.
 *
 * Part of the behavior of the \c same selection keyword is hardcoded into
 * this function (or rather, into _gmx_selelem_custom_init_same()) to allow the
 * use of any keyword in \c "same KEYWORD as" without requiring special
 * handling somewhere else (or sacrificing the simple syntax).
 */
SelectionTreeElementPointer
_gmx_sel_init_method(gmx_ana_selmethod_t                      *method,
                     gmx::SelectionParserParameterListPointer  params,
                     const char *rpost, yyscan_t scanner)
{
    gmx_ana_selcollection_t     *sc = _gmx_sel_lexer_selcollection(scanner);
    int                          rc;

    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    char  buf[128];
    sprintf(buf, "In keyword '%s'", method->name);
    gmx::MessageStringContext  context(errors, buf);

    _gmx_sel_finish_method(scanner);
    /* The "same" keyword needs some custom massaging of the parameters. */
    rc = _gmx_selelem_custom_init_same(&method, params, scanner);
    if (rc != 0)
    {
        return SelectionTreeElementPointer();
    }
    SelectionTreeElementPointer root(new SelectionTreeElement(SEL_EXPRESSION));
    _gmx_selelem_set_method(root, method, scanner);
    /* Process the parameters */
    if (!_gmx_sel_parse_params(*params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root, scanner))
    {
        return SelectionTreeElementPointer();
    }
    set_refpos_type(&sc->pcc, root, rpost, scanner);

    return root;
}

/*!
 * \param[in]  method Modifier to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  sel    Selection element that the modifier should act on.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * selection modifiers.
 */
SelectionTreeElementPointer
_gmx_sel_init_modifier(gmx_ana_selmethod_t                      *method,
                       gmx::SelectionParserParameterListPointer  params,
                       const gmx::SelectionTreeElementPointer   &sel,
                       yyscan_t                                  scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    char  buf[128];
    sprintf(buf, "In keyword '%s'", method->name);
    gmx::MessageStringContext  context(errors, buf);

    _gmx_sel_finish_method(scanner);
    SelectionTreeElementPointer modifier(new SelectionTreeElement(SEL_MODIFIER));
    _gmx_selelem_set_method(modifier, method, scanner);
    SelectionTreeElementPointer root;
    if (method->type == NO_VALUE)
    {
        SelectionTreeElementPointer child = sel;
        while (child->next)
        {
            child = child->next;
        }
        child->next = modifier;
        root        = sel;
    }
    else
    {
        params->push_front(
                SelectionParserParameter::createFromExpression(NULL, sel));
        root = modifier;
    }
    /* Process the parameters */
    if (!_gmx_sel_parse_params(*params, modifier->u.expr.method->nparams,
                               modifier->u.expr.method->param, modifier, scanner))
    {
        return SelectionTreeElementPointer();
    }

    return root;
}

/*!
 * \param[in]  expr    Input selection element for the position calculation.
 * \param[in]  type    Reference position type or NULL for default.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a gmx::SelectionTreeElement object for
 * evaluation of reference positions.
 */
SelectionTreeElementPointer
_gmx_sel_init_position(const gmx::SelectionTreeElementPointer &expr,
                       const char *type, yyscan_t scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    char  buf[128];
    sprintf(buf, "In position evaluation");
    gmx::MessageStringContext   context(errors, buf);

    SelectionTreeElementPointer root(new SelectionTreeElement(SEL_EXPRESSION));
    _gmx_selelem_set_method(root, &sm_keyword_pos, scanner);
    _gmx_selelem_set_kwpos_type(root.get(), type);
    /* Create the parameters for the parameter parser. */
    SelectionParserParameterList params;
    params.push_back(SelectionParserParameter::createFromExpression(NULL, expr));
    /* Parse the parameters. */
    if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root, scanner))
    {
        return SelectionTreeElementPointer();
    }

    return root;
}

/*!
 * \param[in] x,y,z  Coordinates for the position.
 * \returns   The creates selection element.
 */
SelectionTreeElementPointer
_gmx_sel_init_const_position(real x, real y, real z)
{
    rvec                        pos;

    SelectionTreeElementPointer sel(new SelectionTreeElement(SEL_CONST));
    _gmx_selelem_set_vtype(sel, POS_VALUE);
    _gmx_selvalue_reserve(&sel->v, 1);
    pos[XX] = x;
    pos[YY] = y;
    pos[ZZ] = z;
    gmx_ana_pos_init_const(sel->v.u.p, pos);
    return sel;
}

/*!
 * \param[in] name  Name of an index group to search for.
 * \param[in] scanner Scanner data structure.
 * \returns   The created selection element.
 *
 * See gmx_ana_indexgrps_find() for information on how \p name is matched
 * against the index groups.
 */
SelectionTreeElementPointer
_gmx_sel_init_group_by_name(const char *name, yyscan_t scanner)
{

    SelectionTreeElementPointer sel(new SelectionTreeElement(SEL_GROUPREF));
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    sel->setName(gmx::formatString("group \"%s\"", name));
    sel->u.gref.name = strdup(name);
    sel->u.gref.id   = -1;

    if (_gmx_sel_lexer_has_groups_set(scanner))
    {
        gmx_ana_indexgrps_t *grps = _gmx_sel_lexer_indexgrps(scanner);
        sel->resolveIndexGroupReference(grps);
    }

    return sel;
}

/*!
 * \param[in] id    Zero-based index number of the group to extract.
 * \param[in] scanner Scanner data structure.
 * \returns   The created selection element.
 */
SelectionTreeElementPointer
_gmx_sel_init_group_by_id(int id, yyscan_t scanner)
{
    SelectionTreeElementPointer sel(new SelectionTreeElement(SEL_GROUPREF));
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    sel->setName(gmx::formatString("group %d", id));
    sel->u.gref.name = NULL;
    sel->u.gref.id   = id;

    if (_gmx_sel_lexer_has_groups_set(scanner))
    {
        gmx_ana_indexgrps_t *grps = _gmx_sel_lexer_indexgrps(scanner);
        sel->resolveIndexGroupReference(grps);
    }

    return sel;
}

/*!
 * \param[in,out] sel  Value of the variable.
 * \returns       The created selection element that references \p sel.
 *
 * The reference count of \p sel is updated, but no other modifications are
 * made.
 */
SelectionTreeElementPointer
_gmx_sel_init_variable_ref(const gmx::SelectionTreeElementPointer &sel)
{
    SelectionTreeElementPointer ref;

    if (sel->v.type == POS_VALUE && sel->type == SEL_CONST)
    {
        ref = sel;
    }
    else
    {
        ref.reset(new SelectionTreeElement(SEL_SUBEXPRREF));
        _gmx_selelem_set_vtype(ref, sel->v.type);
        ref->setName(sel->name());
        ref->child = sel;
    }
    return ref;
}

/*!
 * \param[in]  name     Name for the selection
 *     (if NULL, a default name is constructed).
 * \param[in]  sel      The selection element that evaluates the selection.
 * \param      scanner  Scanner data structure.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root (\ref SEL_ROOT)
 * gmx::SelectionTreeElement objects for selections.
 */
SelectionTreeElementPointer
_gmx_sel_init_selection(const char                             *name,
                        const gmx::SelectionTreeElementPointer &sel,
                        yyscan_t                                scanner)
{
    if (sel->v.type != POS_VALUE)
    {
        /* FIXME: Better handling of this error */
        GMX_THROW(gmx::InternalError(
                          "Each selection must evaluate to a position"));
    }

    SelectionTreeElementPointer root(new SelectionTreeElement(SEL_ROOT));
    root->child = sel;
    if (name)
    {
        root->setName(name);
    }
    /* Update the flags */
    _gmx_selelem_update_flags(root);
    gmx::ExceptionInitializer errors("Invalid index group reference(s)");
    root->checkUnsortedAtoms(true, &errors);
    if (errors.hasNestedExceptions())
    {
        GMX_THROW(gmx::InconsistentInputError(errors));
    }

    root->fillNameIfMissing(_gmx_sel_lexer_pselstr(scanner));

    /* Print out some information if the parser is interactive */
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        fprintf(stderr, "Selection '%s' parsed\n",
                _gmx_sel_lexer_pselstr(scanner));
    }

    return root;
}


/*!
 * \param[in]  name     Name of the variable.
 * \param[in]  expr     The selection element that evaluates the variable.
 * \param      scanner  Scanner data structure.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root gmx::SelectionTreeElement objects
 * for variable assignments. A \ref SEL_ROOT element and a \ref SEL_SUBEXPR
 * element are both created.
 */
SelectionTreeElementPointer
_gmx_sel_assign_variable(const char                             *name,
                         const gmx::SelectionTreeElementPointer &expr,
                         yyscan_t                                scanner)
{
    gmx_ana_selcollection_t     *sc      = _gmx_sel_lexer_selcollection(scanner);
    const char                  *pselstr = _gmx_sel_lexer_pselstr(scanner);
    SelectionTreeElementPointer  root;

    _gmx_selelem_update_flags(expr);
    /* Check if this is a constant non-group value */
    if (expr->type == SEL_CONST && expr->v.type != GROUP_VALUE)
    {
        /* If so, just assign the constant value to the variable */
        sc->symtab->addVariable(name, expr);
    }
    /* Check if we are assigning a variable to another variable */
    else if (expr->type == SEL_SUBEXPRREF)
    {
        /* If so, make a simple alias */
        sc->symtab->addVariable(name, expr->child);
    }
    else
    {
        /* Create the root element */
        root.reset(new SelectionTreeElement(SEL_ROOT));
        root->setName(name);
        /* Create the subexpression element */
        root->child.reset(new SelectionTreeElement(SEL_SUBEXPR));
        root->child->setName(name);
        _gmx_selelem_set_vtype(root->child, expr->v.type);
        root->child->child  = expr;
        /* Update flags */
        _gmx_selelem_update_flags(root);
        gmx::ExceptionInitializer errors("Invalid index group reference(s)");
        root->checkUnsortedAtoms(true, &errors);
        if (errors.hasNestedExceptions())
        {
            GMX_THROW(gmx::InconsistentInputError(errors));
        }
        /* Add the variable to the symbol table */
        sc->symtab->addVariable(name, root->child);
    }
    srenew(sc->varstrs, sc->nvars + 1);
    sc->varstrs[sc->nvars] = strdup(pselstr);
    ++sc->nvars;
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        fprintf(stderr, "Variable '%s' parsed\n", pselstr);
    }
    return root;
}

/*!
 * \param         sel   Selection to append (can be NULL, in which
 *   case nothing is done).
 * \param         last  Last selection, or NULL if not present or not known.
 * \param         scanner  Scanner data structure.
 * \returns       The last selection after the append.
 *
 * Appends \p sel after the last root element, and returns either \p sel
 * (if it was non-NULL) or the last element (if \p sel was NULL).
 */
SelectionTreeElementPointer
_gmx_sel_append_selection(const gmx::SelectionTreeElementPointer &sel,
                          gmx::SelectionTreeElementPointer        last,
                          yyscan_t                                scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);

    /* Append sel after last, or the last element of sc if last is NULL */
    if (last)
    {
        last->next = sel;
    }
    else
    {
        if (sc->root)
        {
            last = sc->root;
            while (last->next)
            {
                last = last->next;
            }
            last->next = sel;
        }
        else
        {
            sc->root = sel;
        }
    }
    /* Initialize a selection object if necessary */
    if (sel)
    {
        last = sel;
        /* Add the new selection to the collection if it is not a variable. */
        if (sel->child->type != SEL_SUBEXPR)
        {
            gmx::SelectionDataPointer selPtr(
                    new gmx::internal::SelectionData(
                            sel.get(), _gmx_sel_lexer_pselstr(scanner)));
            sc->sel.push_back(gmx::move(selPtr));
        }
    }
    /* Clear the selection string now that we've saved it */
    _gmx_sel_lexer_clear_pselstr(scanner);
    return last;
}

/*!
 * \param[in] scanner Scanner data structure.
 * \returns   true if the parser should finish, false if parsing should
 *   continue.
 *
 * This function is called always after _gmx_sel_append_selection() to
 * check whether a sufficient number of selections has already been provided.
 * This is used to terminate interactive parsers when the correct number of
 * selections has been provided.
 */
bool
_gmx_sel_parser_should_finish(yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    return (int)sc->sel.size() == _gmx_sel_lexer_exp_selcount(scanner);
}
