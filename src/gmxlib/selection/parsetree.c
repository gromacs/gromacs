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
 * \brief Implementation of functions of in parsetree.h.
 */
/*! \internal
 * \page selparser Selection parsing
 *
 * \todo
 * Write some more details of how the parser works.
 *
 *
 * \section selparser_tree Element tree constructed by the parser
 *
 * The parser initializes the \c t_selelem::name (for most elements),
 * \c t_selelem::type, and \c t_selelem::v\c .type fields, as well as the
 * \c t_selelem::child, \c t_selelem::next, and \c t_selelem::refcount fields.
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
 * The name of the selection/variable is stored in \c t_selelem::cgrp\c .name
 * (this is NULL if an explicit name has not been provided for a selection).
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
 * \subsection selparser_tree_bool Boolean elements
 *
 * One \ref SEL_BOOLEAN element is created for each boolean keyword in the
 * input, and the tree structure represents the evaluation order.
 * The \c t_selelem::boolt type gives the type of the operation.
 * Each element has exactly two children (one for \ref BOOL_NOT elements),
 * which are in the order given in the input.
 * The children always have \ref GROUP_VALUE, but different element types
 * are possible.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include <futil.h>
#include <smalloc.h>
#include <string2.h>

#include <poscalc.h>
#include <selection.h>
#include <selmethod.h>

#include "keywords.h"
#include "parsetree.h"
#include "selcollection.h"
#include "selelem.h"
#include "symrec.h"

#include "scanner.h"

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
_gmx_selexpr_create_param(const char *name)
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
        old = value;
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
        old = param;
        param = param->next;
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
    t_selelem          *child;
    int                 rc;
    bool                bUseChildType;

    /* Return if the flags have already been set */
    if (sel->flags & SEL_FLAGSSET)
    {
        return 0;
    }
    /* Set the flags based on the current element type */
    switch (sel->type)
    {
        case SEL_CONST:
            sel->flags |= SEL_SINGLEVAL;
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

        default:
            bUseChildType = TRUE;
            break;
    }
    /* Loop through children to propagate their flags upwards */
    child = sel->child;
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

        child = child->next;
    }
    /* Mark that the flags are set */
    sel->flags |= SEL_FLAGSSET;
    /* For root elements, the type should be propagated here, after the
     * children have been updated. */
    if (sel->type == SEL_ROOT)
    {
        sel->flags |= (sel->child->flags & SEL_VALTYPEMASK);
    }
    return 0;
}

/*! \brief
 * Initializes the method parameter data of \ref SEL_EXPRESSION and
 * \ref SEL_MODIFIER elements.
 *
 * \param[in]     sc     Selection collection.
 * \param[in,out] sel    Selection element to initialize.
 *
 * A deep copy of the parameters is made to allow several
 * expressions with the same method to coexist peacefully.
 * Calls sel_datafunc() if one is specified for the method.
 */
static void
init_method_params(gmx_ana_selcollection_t *sc, t_selelem *sel)
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
        _gmx_selvalue_clear(&param[i].val);
        if (param[i].flags & SPAR_VARNUM)
        {
            param[i].val.nr = -1;
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
        sel->u.expr.method->set_poscoll(sc->pcc, mdata);
    }
    /* Store the values */
    sel->u.expr.method->param = param;
    sel->u.expr.mdata         = mdata;
}

/*! \brief
 * Initializes the method for a \ref SEL_EXPRESSION selection element.
 *
 * \param[in]     sc     Selection collection.
 * \param[in,out] sel    Selection element to initialize.
 * \param[in]     method Selection method to set.
 *
 * Makes a copy of \p method and stores it in \p sel->u.expr.method,
 * and calls init_method_params();
 */
static void
set_method(gmx_ana_selcollection_t *sc, t_selelem *sel,
           gmx_ana_selmethod_t *method)
{
    int      i;

    _gmx_selelem_set_vtype(sel, method->type);
    sel->name   = method->name;
    snew(sel->u.expr.method, 1);
    memcpy(sel->u.expr.method, method, sizeof(gmx_ana_selmethod_t));
    init_method_params(sc, sel);
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
        _gmx_selparser_error("warning: '%d' modifier for '%s' ignored",
                             rpost, sel->u.expr.method->name);
    }
    return rc;
}

/*!
 * \param[in]  sc     Selection collection.
 * \param[in]  left   Selection element for the left hand side.
 * \param[in]  right  Selection element for the right hand side.
 * \param[in]  cmpop  String representation of the comparison operator.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * comparison expressions.
 */
t_selelem *
_gmx_sel_init_comparison(gmx_ana_selcollection_t *sc,
                         t_selelem *left, t_selelem *right, char *cmpop)
{
    t_selelem         *sel;
    t_selexpr_param   *params, *param;
    int                rc;

    sel = _gmx_selelem_create(SEL_EXPRESSION);
    set_method(sc, sel, &sm_compare);
    /* Create the parameter for the left expression */
    params = param     = _gmx_selexpr_create_param(left->v.type == INT_VALUE ? "int1" : "real1");
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value_expr(left);
    /* Create the parameter for the right expression */
    param              = _gmx_selexpr_create_param(right->v.type == INT_VALUE ? "int2" : "real2");
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value_expr(right);
    params->next       = param;
    /* Create the parameter for the operator */
    param              = _gmx_selexpr_create_param("op");
    param->nval        = 1;
    param->value       = _gmx_selexpr_create_value(STR_VALUE);
    param->value->u.s  = cmpop;
    params->next->next = param;
    if (!_gmx_sel_parse_params(params, sel->u.expr.method->nparams,
                               sel->u.expr.method->param, sel))
    {
        _gmx_selparser_error("error in comparison initialization");
        _gmx_selelem_free(sel);
        return NULL;
    }

    return sel;
}

/*!
 * \param[in]  sc     Selection collection (used for position evaluation).
 * \param[in]  method Method to use.
 * \param[in]  nargs  Number of arguments for keyword matching.
 * \param[in]  args   Pointer to the first argument.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection methods that do not take parameters.
 */
t_selelem *
_gmx_sel_init_keyword(gmx_ana_selcollection_t *sc,
                      gmx_ana_selmethod_t *method, int nargs,
                      t_selexpr_value *args, const char *rpost)
{
    t_selelem         *root, *child;
    t_selexpr_param   *params, *param;
    int                rc;

    if (method->nparams > 0)
    {
        gmx_bug("internal error");
        return NULL;
    }

    root = _gmx_selelem_create(SEL_EXPRESSION);
    child = root;
    set_method(sc, child, method);

    /* Initialize the evaluation of keyword matching if values are provided */
    if (nargs > 0)
    {
        gmx_ana_selmethod_t *kwmethod;
        switch (method->type)
        {
            case INT_VALUE: kwmethod = &sm_keyword_int; break;
            case STR_VALUE: kwmethod = &sm_keyword_str; break;
            default:
                _gmx_selparser_error("unknown type for keyword selection");
                goto on_error;
        }
        root = _gmx_selelem_create(SEL_EXPRESSION);
        set_method(sc, root, kwmethod);
        params = param = _gmx_selexpr_create_param(NULL);
        param->nval    = 1;
        param->value   = _gmx_selexpr_create_value_expr(child);
        param          = _gmx_selexpr_create_param(NULL);
        param->nval    = nargs;
        param->value   = args;
        params->next   = param;
        if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                                   root->u.expr.method->param, root))
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
 * \param[in]  sc      Selection collection.
 * \param[in]  expr    Input selection element for the position calculation.
 * \param[in]  type    Reference position type or NULL for default.
 * \param[in]  bSelPos Whether the element evaluates the positions for a
 *   selection.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * evaluation of reference positions.
 */
t_selelem *
_gmx_sel_init_position(gmx_ana_selcollection_t *sc, t_selelem *expr,
                       const char *type, bool bSelPos)
{
    t_selelem       *root;
    t_selexpr_param *params;
    int              flags;

    root = _gmx_selelem_create(SEL_EXPRESSION);
    set_method(sc, root, &sm_keyword_pos);
    /* Selections use largest static group by default, while
     * reference positions use the whole residue/molecule. */
    flags = bSelPos ? POS_COMPLMAX : POS_COMPLWHOLE;
    if (bSelPos && sc->bMaskOnly)
    {
        flags |= POS_MASKONLY;
    }
    /* FIXME: It would be better not to have the string here hardcoded. */
    if (type[0] != 'a')
    {
        root->u.expr.method->flags |= SMETH_REQTOP;
    }
    _gmx_selelem_set_kwpos_type(type, flags, root->u.expr.mdata);
    /* Create the parameters for the parameter parser. */
    params        = _gmx_selexpr_create_param(NULL);
    params->nval  = 1;
    params->value = _gmx_selexpr_create_value_expr(expr);
    /* Parse the parameters. */
    if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root))
    {
        _gmx_selelem_free(root);
        return NULL;
    }

    return root;
}

/*!
 * \param[in]  sc     Selection collection (used for position evaluation).
 * \param[in]  method Method to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  rpost  Reference position type to use (NULL = default).
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection methods that take parameters.
 */
t_selelem *
_gmx_sel_init_method(gmx_ana_selcollection_t *sc, gmx_ana_selmethod_t *method,
                     t_selexpr_param *params, const char *rpost)
{
    t_selelem       *root;
    int              rc;

    root = _gmx_selelem_create(SEL_EXPRESSION);
    set_method(sc, root, method);
    /* Process the parameters */
    if (!_gmx_sel_parse_params(params, root->u.expr.method->nparams,
                               root->u.expr.method->param, root))
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
 * \param[in]  sc     Selection collection.
 * \param[in]  method Modifier to use for initialization.
 * \param[in]  params Pointer to the first parameter.
 * \param[in]  sel    Selection element that the modifier should act on.
 * \returns    The created selection element.
 *
 * This function handles the creation of a \c t_selelem object for
 * selection modifiers.
 */
t_selelem *
_gmx_sel_init_modifier(gmx_ana_selcollection_t *sc, gmx_ana_selmethod_t *method,
                       t_selexpr_param *params, t_selelem *sel)
{
    t_selelem         *root;
    t_selelem         *mod;
    t_selexpr_param   *vparam;
    int                i;

    mod = _gmx_selelem_create(SEL_MODIFIER);
    set_method(sc, mod, method);
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
                               mod->u.expr.method->param, mod))
    {
        if (mod->child != sel)
        {
            _gmx_selelem_free(sel);
        }
        _gmx_selelem_free(mod);
        return NULL;
    }

    return root;
}

/*!
 * \param      scanner  Scanner data structure.
 * \param[in]  sel      The selection element that evaluates the selection.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root (\ref SEL_ROOT) \c t_selelem
 * objects for selections.
 */
t_selelem *
_gmx_sel_init_selection(gmx_sel_lexer_t *scanner, t_selelem *sel)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root;
    int                      rc;

    if (sel->v.type != POS_VALUE)
    {
        gmx_bug("each selection must evaluate to a position");
        /* FIXME: Better handling of this error */
        return NULL;
    }

    root = _gmx_selelem_create(SEL_ROOT);
    root->child = sel;
    /* Update the flags */
    rc = _gmx_selelem_update_flags(root);
    if (rc != 0)
    {
        _gmx_selelem_free(root);
        return NULL;
    }

    /* Print out some information if the parser is interactive */
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        /* TODO: It would be nice to print the whole selection here */
        fprintf(stderr, "Selection parsed\n");
    }

    return root;
}


/*!
 * \param      scanner  Scanner data structure.
 * \param[in]  name     Name of the variable (should not be freed after this
 *   function).
 * \param[in]  expr     The selection element that evaluates the variable.
 * \returns    The created root selection element.
 *
 * This function handles the creation of root \c t_selelem objects for
 * variable assignments. A \ref SEL_ROOT element and a \ref SEL_SUBEXPR
 * element are both created.
 */
t_selelem *
_gmx_sel_assign_variable(gmx_sel_lexer_t *scanner, char *name, t_selelem *expr)
{
    gmx_ana_selcollection_t *sc = _gmx_sel_lexer_selcollection(scanner);
    t_selelem               *root;
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
        if (_gmx_sel_is_lexer_interactive(scanner))
        {
            fprintf(stderr, "Variable '%s' parsed\n", name);
        }
        sfree(name);
        return NULL;
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
        if (_gmx_sel_is_lexer_interactive(scanner))
        {
            fprintf(stderr, "Variable '%s' parsed\n", name);
        }
        sfree(name);
        return NULL;
    }
    /* Create the root element */
    root = _gmx_selelem_create(SEL_ROOT);
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
    if (_gmx_sel_is_lexer_interactive(scanner))
    {
        fprintf(stderr, "Variable '%s' parsed\n", name);
    }
    return root;
}

/*!
 * \param[in,out] sc    Selection collection to append to.
 * \param         sel   Selection to append to \p sc (can be NULL, in which
 *   case nothing is done).
 * \param         last  Last selection in \p sc, or NULL if not present or not
 *   known.
 * \returns       The last selection in \p sc after the append.
 *
 * Appends \p sel after the last root element in \p sc, and returns either
 * \p sel (if it was non-NULL) or the last element in \p sc (if \p sel was
 * NULL).
 */
t_selelem *
_gmx_sel_append_selection(gmx_ana_selcollection_t *sc, t_selelem *sel,
                          t_selelem *last)
{
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

            if (sel->child->type == SEL_CONST)
            {
                gmx_ana_pos_copy(&sc->sel[i]->p, sel->child->v.u.p, TRUE);
                sc->sel[i]->bDynamic = FALSE;
            }
            else
            {
                t_selelem *child;

                child = sel->child;
                child->flags     &= ~SEL_ALLOCVAL;
                _gmx_selvalue_setstore(&child->v, &sc->sel[i]->p);
                /* We should also skip any modifiers to determine the dynamic
                 * status. */
                while (child->type == SEL_MODIFIER)
                {
                    child = child->child;
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
    return last;
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
                                  gmx_ana_indexgrps_t *grps, bool bInteractive)
{
    gmx_sel_lexer_t *scanner;

    _gmx_sel_init_lexer(&scanner, sc, bInteractive);
    _gmx_sel_set_lex_input_file(scanner, stdin);
    return _gmx_sel_run_parser(scanner, sc, grps, nr);
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
    gmx_sel_lexer_t *scanner;
    FILE *fp;
    int   rc;

    _gmx_sel_init_lexer(&scanner, sc, FALSE);
    fp = ffopen(fnm, "r");
    _gmx_sel_set_lex_input_file(scanner, fp);
    rc = _gmx_sel_run_parser(scanner, sc, grps, -1);
    fclose(fp);
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
    gmx_sel_lexer_t *scanner;

    _gmx_sel_init_lexer(&scanner, sc, FALSE);
    _gmx_sel_set_lex_input_str(scanner, str);
    return _gmx_sel_run_parser(scanner, sc, grps, -1);
}
