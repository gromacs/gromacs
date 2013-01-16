/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \brief Implementation of functions of in parsetree.h.
 */
/*! \internal
 * \page selparser Selection parsing
 *
 * The selection parser is implemented in the following files:
 *  - scanner.l:
 *    Tokenizer implemented using Flex, splits the input into tokens
 *    (scanner.c and scanner_flex.h are generated from this file).
 *  - scanner.h, scanner_internal.h, scanner_internal.c:
 *    Helper functions for scanner.l and for interfacing between
 *    scanner.l and parser.y. Functions in scanner_internal.h are only
 *    used from scanner.l, while scanner.h is used from the parser.
 *  - symrec.h, symrec.c:
 *    Functions used by the tokenizer to handle the symbol table, i.e.,
 *    the recognized keywords. Some basic keywords are hardcoded into
 *    scanner.l, but all method and variable references go through the
 *    symbol table, as do position evaluation keywords.
 *  - parser.y:
 *    Semantic rules for parsing the grammar
 *    (parser.c and parser.h are generated from this file by Bison).
 *  - parsetree.h, parsetree.c:
 *    Functions called from actions in parser.y to construct the
 *    evaluation elements corresponding to different grammar elements.
 *    parsetree.c also defines the external interface of the parser,
 *    i.e., the \c gmx_ana_selcollection_parse_*() functions declared
 *    in selection.h.
 *  - params.c:
 *    Defines a function that processes the parameters of selection
 *    methods and initializes the children of the method element.
 *
 * The basic control flow in the parser is as follows: when a
 * fmx_ana_selcollection_parse_*() gets called, it performs some
 * initialization, and then calls the _gmx_sel_yyparse() function generated
 * by Bison. This function then calls _gmx_sel_yylex() to repeatedly read
 * tokens from the input (more complex tasks related to token recognition
 * and bookkeeping are done by functions in scanner_internal.c) and uses the
 * grammar rules to decide what to do with them. Whenever a grammar rule
 * matches, a corresponding function in parsetree.c is called to construct
 * either a temporary representation for the object or a \c t_selelem object
 * (some simple rules are handled internally in parser.y).
 * When a complete selection has been parsed, the functions in parsetree.c
 * also take care of updating the \c gmx_ana_selcollection_t structure
 * appropriately.
 *
 * The rest of this page describes the resulting \c t_selelem object tree.
 * Before the selections can be evaluated, this tree needs to be passed to
 * the selection compiler, which is described on a separate page:
 * \ref selcompiler
 *
 *
 * \section selparser_tree Element tree constructed by the parser
 *
 * The parser initializes the following fields in all selection elements:
 * \c t_selelem::name, \c t_selelem::type, \c t_selelem::v\c .type,
 * \c t_selelem::flags, \c t_selelem::child, \c t_selelem::next, and
 * \c t_selelem::refcount.
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
 * The name of the selection/variable is stored in \c t_selelem::cgrp\c .name.
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
 * For group-valued elements, the value is stored in \c t_selelem::cgrp;
 * other types of values are stored in \c t_selelem::v.
 * Constants that appear as parameters for selection methods are not present
 * in the selection tree unless they have \ref GROUP_VALUE.
 * \ref SEL_CONST elements have no children.
 *
 *
 * \subsection selparser_tree_method Method evaluation elements
 *
 * \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements are treated very
 * similarly. The \c gmx_ana_selmethod_t structure corresponding to the
 * evaluation method is in \c t_selelem::method, and the method data in
 * \c t_selelem::mdata has been allocated using sel_datafunc().
 * If a non-standard reference position type was set, \c t_selelem::pc has
 * also been created, but only the type has been set.
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
 * \c t_selelem::name points to the name of the variable (from the
 * \ref SEL_ROOT element).
 * The element always has exactly one child, which represents the value of
 * the variable.
 * \ref SEL_SUBEXPR element is the only element type that can have
 * \c t_selelem::refcount different from 1.
 *
 * \ref SEL_SUBEXPRREF elements are used for two purposes:
 *  - Variable references that need to be evaluated (i.e., there is a
 *    \ref SEL_SUBEXPR element for the variable) are represented using
 *    \ref SEL_SUBEXPRREF elements.
 *    In this case, \c t_selelem::param is NULL, and the first and only
 *    child of the element is the \ref SEL_SUBEXPR element of the variable.
 *    Such references can appear anywhere where the variable value
 *    (the child of the \ref SEL_SUBEXPR element) would be valid.
 *  - Children of \ref SEL_EXPRESSION and \ref SEL_MODIFIER elements are
 *    always of this type. For these elements, \c t_selelem::param is
 *    initialized to point to the parameter that receives the value from
 *    the expression.
 *    Each such element has exactly one child, which can be of any type;
 *    the \ref SEL_SUBEXPR element of a variable is used if the value comes
 *    from a variable, otherwise the child type is not \ref SEL_SUBEXPR.
 *
 *
 * \subsection selparser_tree_gmx_bool Boolean elements
 *
 * One \ref SEL_BOOLEAN element is created for each gmx_boolean keyword in the
 * input, and the tree structure represents the evaluation order.
 * The \c t_selelem::boolt type gives the type of the operation.
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
 * The \c t_selelem::optype type gives the name of the operation.
 * Each element has exactly two children (one for unary negation elements),
 * which are in the order given in the input.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdarg.h>

#include <futil.h>
#include <smalloc.h>
#include <string2.h>
#include <gmx_fatal.h>
#include <assert.h>

#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>

#include "keywords.h"
#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "selhelp.h"
#include "symrec.h"

#include "scanner.h"

/* In parser.y */
/*! \brief
 * Parser function generated by Bison.
 */
int
_gmx_sel_yybparse(void *scanner);

/*!
 * It is a simple wrapper for fprintf(stderr, ...).
 */
void
_gmx_selparser_error(const char *fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "selection parser: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
}

/*!
 * \param[in] type  Type for the new value.
 * \returns   Pointer to the newly allocated value.
 */
t_selexpr_value *
_gmx_selexpr_create_value(e_selvalue_t type)
{
    t_selexpr_value *value;
    snew(value, 1);
    value->type  = type;
    value->bExpr = FALSE;
    value->next  = NULL;
    return value;
}

/*!
 * \param[in] expr  Expression for the value.
 * \returns   Pointer to the newly allocated value.
 */
t_selexpr_value *
_gmx_selexpr_create_value_expr(t_selelem *expr)
{
    t_selexpr_value *value;
    snew(value, 1);
    value->type   = expr->v.type;
    value->bExpr  = TRUE;
    value->u.expr = expr;
    value->next   = NULL;
    return value;
}

/*!
 * \param[in] name Name for the new parameter.
 * \returns   Pointer to the newly allocated parameter.
 *
 * No copy of \p name is made.
 */
t_selexpr_param *
_gmx_selexpr_create_param(char *name)
{
    t_selexpr_param *param;
    snew(param, 1);
    param->name = name;
    param->next = NULL;
    return param;
}

/*!
 * \param value Pointer to the beginning of the value list to free.
 *
 * The expressions referenced by the values are also freed
 * (to prevent this, set the expression to NULL before calling the function).
 */
void
_gmx_selexpr_free_values(t_selexpr_value *value)
{
    t_selexpr_value *old;

    while (value)
    {
        if (value->bExpr)
        {
            if (value->u.expr)
            {
                _gmx_selelem_free(value->u.expr);
            }
        }
        else if (value->type == STR_VALUE)
        {
            sfree(value->u.s);
        }
        old   = value;
        value = value->next;
        sfree(old);
    }
}

/*!
 * \param param Pointer the the beginning of the parameter list to free.
 *
 * The values of the parameters are freed with free_selexpr_values().
 */
void
_gmx_selexpr_free_params(t_selexpr_param *param)
{
    t_selexpr_param *old;

    while (param)
    {
        _gmx_selexpr_free_values(param->value);
        old   = param;
        param = param->next;
        sfree(old->name);
        sfree(old);
    }
}

/*!
 * \param[in,out] sel  Root of the selection element tree to initialize.
 * \returns       0 on success, an error code on error.
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
int
_gmx_selelem_update_flags(t_selelem *sel)
{
    t_selelem              *child;
    int                     rc;
    gmx_bool                bUseChildType = FALSE;
    gmx_bool                bOnlySingleChildren;

    /* Return if the flags have already been set */
    if (sel->flags & SEL_FLAGSSET)
    {
        return 0;
    }
    /* Set the flags based on the current element type */
    switch (sel->type)
    {
        case SEL_CONST:
            sel->flags   |= SEL_SINGLEVAL;
            bUseChildType = FALSE;
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
            bUseChildType = FALSE;
            break;

        case SEL_ARITHMETIC:
            sel->flags   |= SEL_ATOMVAL;
            bUseChildType = FALSE;
            break;

        case SEL_MODIFIER:
            if (sel->v.type != NO_VALUE)
            {
                sel->flags |= SEL_VARNUMVAL;
            }
            bUseChildType = FALSE;
            break;

        case SEL_ROOT:
            bUseChildType = FALSE;
            break;

        case SEL_BOOLEAN:
        case SEL_SUBEXPR:
        case SEL_SUBEXPRREF:
            bUseChildType = TRUE;
            break;
    }
    /* Loop through children to propagate their flags upwards */
    bOnlySingleChildren = TRUE;
    child               = sel->child;
    while (child)
    {
        /* Update the child */
        rc = _gmx_selelem_update_flags(child);
        if (rc != 0)
        {
            return rc;
        }
        /* Propagate the dynamic flag */
        sel->flags |= (child->flags & SEL_DYNAMIC);
        /* Propagate the type flag if necessary and check for problems */
        if (bUseChildType)
        {
            if ((sel->flags & SEL_VALTYPEMASK)
                && !(sel->flags & child->flags & SEL_VALTYPEMASK))
            {
                _gmx_selparser_error("invalid combination of selection expressions");
                return EINVAL;
            }
            sel->flags |= (child->flags & SEL_VALTYPEMASK);
        }
        if (!(child->flags & SEL_SINGLEVAL))
        {
            bOnlySingleChildren = FALSE;
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
        assert(sel->child);
        sel->flags |= (sel->child->flags & SEL_VALTYPEMASK);
    }
    /* Mark that the flags are set */
    sel->flags |= SEL_FLAGSSET;
    return 0;
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
_gmx_selelem_init_method_params(t_selelem *sel, yyscan_t scanner)
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
        if (mdata == NULL)
        {
            gmx_fatal(FARGS, "Method data initialization failed");
        }
    }
    if (sel->u.expr.method->set_poscoll)
    {
        gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);

        sel->u.expr.method->set_poscoll(sc->pcc, mdata);
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
_gmx_selelem_set_method(t_selelem *sel, gmx_ana_selmethod_t *method,
                        yyscan_t scanner)
{
    int      i;

    _gmx_selelem_set_vtype(sel, method->type);
    sel->name   = method->name;
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
 * \returns       0 on success, a non-zero error code on error.
 */
static int
set_refpos_type(gmx_ana_poscalc_coll_t *pcc, t_selelem *sel, const char *rpost)
{
    int  rc;

    if (!rpost)
    {
        return 0;
    }

    rc = 0;
    if (sel->u.expr.method->pupdate)
    {
        /* By default, use whole residues/molecules. */
        rc = gmx_ana_poscalc_create_enum(&sel->u.expr.pc, pcc, rpost,
                                         POS_COMPLWHOLE);
    }
    else
    {
        _gmx_selparser_error("warning: '%s' modifier for '%s' ignored",
                             rpost, sel->u.expr.method->name);
    }
    return rc;
}

/*!
 * \param[in]  left    Selection element for the left hand side.
 * \param[in]  right   Selection element for the right hand side.
 * \param[in]  op      String representation of the operator.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * arithmetic expressions.
 */
t_selelem *
_gmx_sel_init_arithmetic(t_selelem *left, t_selelem *right, char op,
                         yyscan_t scanner)
{
    t_selelem         *sel;
    char               buf[2];

    buf[0]             = op;
    buf[1]             = 0;
    sel                = _gmx_selelem_create(SEL_ARITHMETIC);
    sel->v.type        = REAL_VALUE;
    switch (op)
    {
        case '+': sel->u.arith.type = ARITH_PLUS; break;
        case '-': sel->u.arith.type = (right ? ARITH_MINUS : ARITH_NEG); break;
        case '*': sel->u.arith.type = ARITH_MULT; break;
        case '/': sel->u.arith.type = ARITH_DIV;  break;
        case '^': sel->u.arith.type = ARITH_EXP;  break;
    }
    sel->u.arith.opstr = strdup(buf);
    sel->name          = sel->u.arith.opstr;
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
 * This function handles the creation of a \c t_selelem object for
 * comparison expressions.
 */
t_selelem *
_gmx_sel_init_comparison(t_selelem *left, t_selelem *right, char *cmpop,
                         yyscan_t scanner)
{
    t_selelem         *sel;
    t_selexpr_param   *params, *param;
    const char        *name;
    int                rc;

    sel = _gmx_selelem_create(SEL_EXPRESSION);
    _gmx_selelem_set_method(sel, &sm_compare, scanner);
    /* Create the parameter for the left expression */
    name               = left->v.type == INT_VALUE ? "int1" : "real1";
    params             = param     = _gmx_selexpr_create_param(strdup(name));
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value_expr(left);
    /* Create the parameter for the right expression */
    name               = right->v.type == INT_VALUE ? "int2" : "real2";
    param              = _gmx_selexpr_create_param(strdup(name));
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value_expr(right);
    params->next       = param;
    /* Create the parameter for the operator */
    param              = _gmx_selexpr_create_param(strdup("op"));
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value(STR_VALUE);
    param->value->u.s  = cmpop;
    params->next->next = param;
    if (!_gmx_sel_parse_params(params, sel->u.expr.method->nparams,
                               sel->u.expr.method->param, sel, scanner))
    {
        _gmx_selparser_error("error in comparison initialization");
        _gmx_selelem_free(sel);
        return NULL;
    }

    return sel;
}

/*!
 * \param[in]  method Method to use.
 * \param[in]  args   Pointer to the first argument.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection methods that do not take parameters.
 */
t_selelem *
_gmx_sel_init_keyword(gmx_ana_selmethod_t *method, t_selexpr_value *args,
                      const char *rpost, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root, *child;
    t_selexpr_param         *params, *param;
    t_selexpr_value         *arg;
    int                      nargs;
    int                      rc;

    if (method->nparams > 0)
    {
        gmx_bug("internal error");
        return NULL;
    }

    root  = _gmx_selelem_create(SEL_EXPRESSION);
    child = root;
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
                _gmx_selparser_error("unknown type for keyword selection");
                _gmx_selexpr_free_values(args);
                goto on_error;
        }
        /* Count the arguments */
        nargs = 0;
        arg   = args;
        while (arg)
        {
            ++nargs;
            arg = arg->next;
        }
        /* Initialize the selection element */
        root = _gmx_selelem_create(SEL_EXPRESSION);
        _gmx_selelem_set_method(root, kwmethod, scanner);
        params         = param = _gmx_selexpr_create_param(NULL);
        param->nval    = 1;
        param->value   = _gmx_selexpr_create_value_expr(child);
        param          = _gmx_selexpr_create_param(NULL);
        param->nval    = nargs;
        param->value   = args;
        params->next   = param;
        if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                                   root->u.expr.method->param, root, scanner))
        {
            _gmx_selparser_error("error in keyword selection initialization");
            goto on_error;
        }
    }
    rc = set_refpos_type(sc->pcc, child, rpost);
    if (rc != 0)
    {
        goto on_error;
    }

    return root;

/* On error, free all memory and return NULL. */
on_error:
    _gmx_selelem_free(root);
    return NULL;
}

/*!
 * \param[in]  method Method to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection methods that take parameters.
 *
 * Part of the behavior of the \c same selection keyword is hardcoded into
 * this function (or rather, into _gmx_selelem_custom_init_same()) to allow the
 * use of any keyword in \c "same KEYWORD as" without requiring special
 * handling somewhere else (or sacrificing the simple syntax).
 */
t_selelem *
_gmx_sel_init_method(gmx_ana_selmethod_t *method, t_selexpr_param *params,
                     const char *rpost, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root;
    int                      rc;

    _gmx_sel_finish_method(scanner);
    /* The "same" keyword needs some custom massaging of the parameters. */
    rc = _gmx_selelem_custom_init_same(&method, params, scanner);
    if (rc != 0)
    {
        _gmx_selexpr_free_params(params);
        return NULL;
    }
    root = _gmx_selelem_create(SEL_EXPRESSION);
    _gmx_selelem_set_method(root, method, scanner);
    /* Process the parameters */
    if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root, scanner))
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    rc = set_refpos_type(sc->pcc, root, rpost);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }

    return root;
}

/*!
 * \param[in]  method Modifier to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  sel    Selection element that the modifier should act on.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection modifiers.
 */
t_selelem *
_gmx_sel_init_modifier(gmx_ana_selmethod_t *method, t_selexpr_param *params,
                       t_selelem *sel, yyscan_t scanner)
{
    t_selelem         *root;
    t_selelem         *mod;
    t_selexpr_param   *vparam;
    int                i;

    _gmx_sel_finish_method(scanner);
    mod = _gmx_selelem_create(SEL_MODIFIER);
    _gmx_selelem_set_method(mod, method, scanner);
    if (method->type == NO_VALUE)
    {
        t_selelem *child;

        child = sel;
        while (child->next)
        {
            child = child->next;
        }
        child->next = mod;
        root        = sel;
    }
    else
    {
        vparam        = _gmx_selexpr_create_param(NULL);
        vparam->nval  = 1;
        vparam->value = _gmx_selexpr_create_value_expr(sel);
        vparam->next  = params;
        params        = vparam;
        root          = mod;
    }
    /* Process the parameters */
    if (!_gmx_sel_parse_params(params, mod->u.expr.method->nparams,
                               mod->u.expr.method->param, mod, scanner))
    {
        _gmx_selelem_free(mod);
        return NULL;
    }

    return root;
}

/*!
 * \param[in]  expr    Input selection element for the position calculation.
 * \param[in]  type    Reference position type or NULL for default.
 * \param[in]  scanner Scanner data structure.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * evaluation of reference positions.
 */
t_selelem *
_gmx_sel_init_position(t_selelem *expr, const char *type, yyscan_t scanner)
{
    t_selelem       *root;
    t_selexpr_param *params;

    root = _gmx_selelem_create(SEL_EXPRESSION);
    _gmx_selelem_set_method(root, &sm_keyword_pos, scanner);
    _gmx_selelem_set_kwpos_type(root, type);
    /* Create the parameters for the parameter parser. */
    params        = _gmx_selexpr_create_param(NULL);
    params->nval  = 1;
    params->value = _gmx_selexpr_create_value_expr(expr);
    /* Parse the parameters. */
    if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root, scanner))
    {
        _gmx_selelem_free(root);
        return NULL;
    }

    return root;
}

/*!
 * \param[in] x,y,z  Coordinates for the position.
 * \returns   The creates selection element.
 */
t_selelem *
_gmx_sel_init_const_position(real x, real y, real z)
{
    t_selelem *sel;
    rvec       pos;

    sel = _gmx_selelem_create(SEL_CONST);
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
 * \returns   The created constant selection element, or NULL if no matching
 *     index group found.
 *
 * See gmx_ana_indexgrps_find() for information on how \p name is matched
 * against the index groups.
 */
t_selelem *
_gmx_sel_init_group_by_name(const char *name, yyscan_t scanner)
{
    gmx_ana_indexgrps_t *grps = _gmx_sel_lexer_indexgrps(scanner);
    t_selelem           *sel;

    if (!grps)
    {
        return NULL;
    }
    sel = _gmx_selelem_create(SEL_CONST);
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    /* FIXME: The constness should not be cast away */
    if (!gmx_ana_indexgrps_find(&sel->u.cgrp, grps, (char *)name))
    {
        _gmx_selelem_free(sel);
        return NULL;
    }
    sel->name = sel->u.cgrp.name;
    return sel;
}

/*!
 * \param[in] id    Zero-based index number of the group to extract.
 * \param[in] scanner Scanner data structure.
 * \returns   The created constant selection element, or NULL if no matching
 *     index group found.
 */
t_selelem *
_gmx_sel_init_group_by_id(int id, yyscan_t scanner)
{
    gmx_ana_indexgrps_t *grps = _gmx_sel_lexer_indexgrps(scanner);
    t_selelem           *sel;

    if (!grps)
    {
        return NULL;
    }
    sel = _gmx_selelem_create(SEL_CONST);
    _gmx_selelem_set_vtype(sel, GROUP_VALUE);
    if (!gmx_ana_indexgrps_extract(&sel->u.cgrp, grps, id))
    {
        _gmx_selelem_free(sel);
        return NULL;
    }
    sel->name = sel->u.cgrp.name;
    return sel;
}

/*!
 * \param[in,out] sel  Value of the variable.
 * \returns       The created selection element that references \p sel.
 *
 * The reference count of \p sel is updated, but no other modifications are
 * made.
 */
t_selelem *
_gmx_sel_init_variable_ref(t_selelem *sel)
{
    t_selelem *ref;

    if (sel->v.type == POS_VALUE && sel->type == SEL_CONST)
    {
        ref = sel;
    }
    else
    {
        ref = _gmx_selelem_create(SEL_SUBEXPRREF);
        _gmx_selelem_set_vtype(ref, sel->v.type);
        ref->name  = sel->name;
        ref->child = sel;
    }
    sel->refcount++;
    return ref;
}

/*! \brief
 * Initializes default values for position keyword evaluation.
 *
 * \param[in,out] root       Root of the element tree to initialize.
 * \param[in]     sc         Selection collection to use defaults from.
 * \param[in]     bSelection Whether the element evaluates the positions for a
 *   selection.
 */
static void
init_pos_keyword_defaults(t_selelem *root, gmx_ana_selcollection_t *sc, gmx_bool bSelection)
{
    t_selelem               *child;
    int                      flags;

    /* Selections use largest static group by default, while
     * reference positions use the whole residue/molecule. */
    if (root->type == SEL_EXPRESSION)
    {
        flags = bSelection ? POS_COMPLMAX : POS_COMPLWHOLE;
        if (bSelection && sc->bMaskOnly)
        {
            flags |= POS_MASKONLY;
        }
        if (bSelection && sc->bVelocities)
        {
            flags |= POS_VELOCITIES;
        }
        if (bSelection && sc->bForces)
        {
            flags |= POS_FORCES;
        }
        _gmx_selelem_set_kwpos_type(root, bSelection ? sc->spost : sc->rpost);
        _gmx_selelem_set_kwpos_flags(root, flags);
    }
    /* Change the defaults once we are no longer processing modifiers */
    if (root->type != SEL_ROOT && root->type != SEL_MODIFIER
        && root->type != SEL_SUBEXPRREF && root->type != SEL_SUBEXPR)
    {
        bSelection = FALSE;
    }
    /* Recurse into children */
    child = root->child;
    while (child)
    {
        init_pos_keyword_defaults(child, sc, bSelection);
        child = child->next;
    }
}

/*!
 * \param[in]  name     Name for the selection
 *     (if NULL, a default name is constructed).
 * \param[in]  sel      The selection element that evaluates the selection.
 * \param      scanner  Scanner data structure.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root (\ref SEL_ROOT) \c t_selelem
 * objects for selections.
 */
t_selelem *
_gmx_sel_init_selection(char *name, t_selelem *sel, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root;
    int                      rc;

    if (sel->v.type != POS_VALUE)
    {
        gmx_bug("each selection must evaluate to a position");
        /* FIXME: Better handling of this error */
        sfree(name);
        return NULL;
    }

    root        = _gmx_selelem_create(SEL_ROOT);
    root->child = sel;
    /* Assign the name (this is done here to free it automatically in the case
     * of an error below). */
    if (name)
    {
        root->name = root->u.cgrp.name = name;
    }
    /* Update the flags */
    rc = _gmx_selelem_update_flags(root);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    /* Initialize defaults for position keywords */
    init_pos_keyword_defaults(sel, sc, TRUE);

    /* If there is no name provided by the user, check whether the actual
     * selection given was from an external group, and if so, use the name
     * of the external group. */
    if (!root->name)
    {
        t_selelem *child = root->child;
        while (child->type == SEL_MODIFIER)
        {
            if (!child->child || child->child->type != SEL_SUBEXPRREF
                || !child->child->child)
            {
                break;
            }
            child = child->child->child;
        }
        if (child->type == SEL_EXPRESSION
            && child->child && child->child->type == SEL_SUBEXPRREF
            && child->child->child
            && child->child->child->type == SEL_CONST
            && child->child->child->v.type == GROUP_VALUE)
        {
            root->name = root->u.cgrp.name =
                    strdup(child->child->child->u.cgrp.name);
        }
    }
    /* If there still is no name, use the selection string */
    if (!root->name)
    {
        root->name = root->u.cgrp.name
                   = strdup(_gmx_sel_lexer_pselstr(scanner));
    }

    /* Print out some information if the parser is interactive */
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        fprintf(stderr, "Selection '%s' parsed\n",
                _gmx_sel_lexer_pselstr(scanner));
    }

    return root;
}


/*!
 * \param[in]  name     Name of the variable (should not be freed after this
 *   function).
 * \param[in]  expr     The selection element that evaluates the variable.
 * \param      scanner  Scanner data structure.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root \c t_selelem objects for
 * variable assignments. A \ref SEL_ROOT element and a \ref SEL_SUBEXPR
 * element are both created.
 */
t_selelem *
_gmx_sel_assign_variable(char *name, t_selelem *expr, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc      = _gmx_sel_lexer_selcollection(scanner);
    const char              *pselstr = _gmx_sel_lexer_pselstr(scanner);
    t_selelem               *root    = NULL;
    int                      rc;

    rc = _gmx_selelem_update_flags(expr);
    if (rc != 0)
    {
        sfree(name);
        _gmx_selelem_free(expr);
        return NULL;
    }
    /* Check if this is a constant non-group value */
    if (expr->type == SEL_CONST && expr->v.type != GROUP_VALUE)
    {
        /* If so, just assign the constant value to the variable */
        if (!_gmx_sel_add_var_symbol(sc->symtab, name, expr))
        {
            _gmx_selelem_free(expr);
            sfree(name);
            return NULL;
        }
        _gmx_selelem_free(expr);
        sfree(name);
        goto finish;
    }
    /* Check if we are assigning a variable to another variable */
    if (expr->type == SEL_SUBEXPRREF)
    {
        /* If so, make a simple alias */
        if (!_gmx_sel_add_var_symbol(sc->symtab, name, expr->child))
        {
            _gmx_selelem_free(expr);
            sfree(name);
            return NULL;
        }
        _gmx_selelem_free(expr);
        sfree(name);
        goto finish;
    }
    /* Create the root element */
    root                = _gmx_selelem_create(SEL_ROOT);
    root->name          = name;
    root->u.cgrp.name   = name;
    /* Create the subexpression element */
    root->child = _gmx_selelem_create(SEL_SUBEXPR);
    _gmx_selelem_set_vtype(root->child, expr->v.type);
    root->child->name   = name;
    root->child->child  = expr;
    /* Update flags */
    rc = _gmx_selelem_update_flags(root);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }
    /* Add the variable to the symbol table */
    if (!_gmx_sel_add_var_symbol(sc->symtab, name, root->child))
    {
        _gmx_selelem_free(root);
        return NULL;
    }
finish:
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
t_selelem *
_gmx_sel_append_selection(t_selelem *sel, t_selelem *last, yyscan_t scanner)
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
            int        i;

            sc->nr++;
            srenew(sc->sel, sc->nr);
            i = sc->nr - 1;
            snew(sc->sel[i], 1);
            sc->sel[i]->name   = strdup(sel->name);
            sc->sel[i]->selstr = strdup(_gmx_sel_lexer_pselstr(scanner));

            if (sel->child->type == SEL_CONST)
            {
                gmx_ana_pos_copy(&sc->sel[i]->p, sel->child->v.u.p, TRUE);
                sc->sel[i]->bDynamic = FALSE;
            }
            else
            {
                t_selelem *child;

                child             = sel->child;
                child->flags     &= ~SEL_ALLOCVAL;
                _gmx_selvalue_setstore(&child->v, &sc->sel[i]->p);
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
                sc->sel[i]->bDynamic = (child->child->flags & SEL_DYNAMIC);
            }
            /* The group will be set after compilation */
            sc->sel[i]->g        = NULL;
            sc->sel[i]->selelem  = sel;
            gmx_ana_selection_init_coverfrac(sc->sel[i], CFRAC_NONE);
        }
    }
    /* Clear the selection string now that we've saved it */
    _gmx_sel_lexer_clear_pselstr(scanner);
    return last;
}

/*!
 * \param[in] scanner Scanner data structure.
 * \returns   TRUE if the parser should finish, FALSE if parsing should
 *   continue.
 *
 * This function is called always after _gmx_sel_append_selection() to
 * check whether a sufficient number of selections has already been provided.
 * This is used to terminate interactive parsers when the correct number of
 * selections has been provided.
 */
gmx_bool
_gmx_sel_parser_should_finish(yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    return sc->nr == _gmx_sel_lexer_exp_selcount(scanner);
}

/*!
 * \param[in] scanner Scanner data structure.
 */
void
_gmx_sel_handle_empty_cmd(yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc   = _gmx_sel_lexer_selcollection(scanner);
    gmx_ana_indexgrps_t     *grps = _gmx_sel_lexer_indexgrps(scanner);
    int                      i;

    if (!_gmx_sel_is_lexer_interactive(scanner))
    {
        return;
    }

    if (grps)
    {
        fprintf(stderr, "Available index groups:\n");
        gmx_ana_indexgrps_print(_gmx_sel_lexer_indexgrps(scanner), 0);
    }
    if (sc->nvars > 0 || sc->nr > 0)
    {
        fprintf(stderr, "Currently provided selections:\n");
        for (i = 0; i < sc->nvars; ++i)
        {
            fprintf(stderr, "     %s\n", sc->varstrs[i]);
        }
        for (i = 0; i < sc->nr; ++i)
        {
            fprintf(stderr, " %2d. %s\n", i+1, sc->sel[i]->selstr);
        }
    }
}

/*!
 * \param[in] topic   Topic for which help was requested, or NULL for general
 *                    help.
 * \param[in] scanner Scanner data structure.
 *
 * \p topic is freed by this function.
 */
void
_gmx_sel_handle_help_cmd(char *topic, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);

    _gmx_sel_print_help(sc, topic);
    if (topic)
    {
        sfree(topic);
    }
}

/*! \brief
 * Internal helper function used by gmx_ana_selcollection_parse_*() to do the actual work.
 *
 * \param[in]     maxnr   Maximum number of selections to parse
 *   (if -1, parse as many as provided by the user).
 * \param[in,out] scanner Scanner data structure.
 * \returns       0 on success, -1 on error.
 */
static int
run_parser(int maxnr, yyscan_t scanner)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    gmx_bool                 bOk;
    int                      nr;

    nr  = sc->nr;
    bOk = !_gmx_sel_yybparse(scanner);
    _gmx_sel_free_lexer(scanner);
    nr = sc->nr - nr;
    if (maxnr > 0 && nr != maxnr)
    {
        return -1;
    }
    return bOk ? 0 : -1;
}

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     nr    Number of selections to parse
 *   (if -1, parse as many as provided by the user).
 * \param[in]     grps  External index groups (can be NULL).
 * \param[in]     bInteractive Whether the parser should behave interactively.
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_stdin(gmx_ana_selcollection_t *sc, int nr,
                                  gmx_ana_indexgrps_t *grps, gmx_bool bInteractive)
{
    yyscan_t scanner;
    int      rc;

    rc = _gmx_sel_init_lexer(&scanner, sc, bInteractive, nr, grps);
    if (rc != 0)
    {
        return rc;
    }
    /* We don't set the lexer input here, which causes it to use a special
     * internal implementation for reading from stdin. */
    return run_parser(nr, scanner);
}

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     fnm   Name of the file to parse selections from.
 * \param[in]     grps  External index groups (can be NULL).
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_file(gmx_ana_selcollection_t *sc, const char *fnm,
                                 gmx_ana_indexgrps_t *grps)
{
    yyscan_t scanner;
    FILE    *fp;
    int      rc;

    rc = _gmx_sel_init_lexer(&scanner, sc, FALSE, -1, grps);
    if (rc != 0)
    {
        return rc;
    }
    fp = ffopen(fnm, "r");
    _gmx_sel_set_lex_input_file(scanner, fp);
    rc = run_parser(-1, scanner);
    ffclose(fp);
    return rc;
}

/*!
 * \param[in,out] sc    Selection collection to use for output.
 * \param[in]     str   String to parse selections from.
 * \param[in]     grps  External index groups (can be NULL).
 * \returns       0 on success, -1 on error.
 *
 * The number of selections parsed can be accessed with
 * gmx_ana_selcollection_get_count() (note that if you call the parser
 * multiple times, this function returns the total count).
 */
int
gmx_ana_selcollection_parse_str(gmx_ana_selcollection_t *sc, const char *str,
                                gmx_ana_indexgrps_t *grps)
{
    yyscan_t scanner;
    int      rc;

    rc = _gmx_sel_init_lexer(&scanner, sc, FALSE, -1, grps);
    if (rc != 0)
    {
        return rc;
    }
    _gmx_sel_set_lex_input_str(scanner, str);
    return run_parser(-1, scanner);
}
