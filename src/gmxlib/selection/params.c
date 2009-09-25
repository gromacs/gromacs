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
 * \brief
 * Implementation of functions in selparam.h.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <smalloc.h>
#include <string2.h>
#include <vec.h>

#include <position.h>
#include <selmethod.h>
#include <selparam.h>

#include "parsetree.h"
#include "position.h"
#include "selelem.h"

/*!
 * \param[in] name   Name of the parameter to search.
 * \param[in] nparam Number of parameters in the \p param array.
 * \param[in] param  Parameter array to search.
 * \returns   Pointer to the parameter in the \p param
 *   or NULL if no parameter with name \p name was found.
 *
 * The comparison is case-sensitive.
 */
gmx_ana_selparam_t *
gmx_ana_selparam_find(const char *name, int nparam, gmx_ana_selparam_t *param)
{
    int                i;

    if (nparam == 0)
    {
        return NULL;
    }
    /* Find the first non-null parameter */
    i = 0;
    while (i < nparam && param[i].name == NULL)
    {
        ++i;
    }
    /* Process the special case of a NULL parameter */
    if (name == NULL)
    {
        return (i == 0) ? NULL : &param[i-1];
    }
    for ( ; i < nparam; ++i)
    {
        if (!strcmp(param[i].name, name))
        {
            return &param[i];
        }
        /* Check for 'no' prefix on boolean parameters */
        if (param[i].val.type == NO_VALUE
            && strlen(name) > 2 && name[0] == 'n' && name[1] == 'o'
            && !strcmp(param[i].name, name+2))
        {
            return &param[i];
        }
    }
    return NULL;
}

/*! \brief
 * Comparison function for sorting integer ranges.
 * 
 * \param[in] a Pointer to the first range.
 * \param[in] b Pointer to the second range.
 * \returns   -1, 0, or 1 depending on the relative order of \p a and \p b.
 *
 * The ranges are primarily sorted based on their starting point, and
 * secondarily based on length (longer ranges come first).
 */
static int
cmp_int_range(const void *a, const void *b)
{
    if (((int *)a)[0] < ((int *)b)[0])
    {
        return -1;
    }
    if (((int *)a)[0] > ((int *)b)[0])
    {
        return 1;
    }
    if (((int *)a)[1] > ((int *)b)[1])
    {
        return -1;
    }
    return 0;
}

/*! \brief
 * Parses the values for a parameter that takes integer ranges.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 */
static bool
parse_values_range(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param)
{
    t_selexpr_value    *value;
    int                *data;
    int                 i, j, n;

    param->flags &= ~SPAR_DYNAMIC;
    if (param->val.type != INT_VALUE)
    {
        gmx_bug("internal error");
        return FALSE;
    }
    snew(data, nval*2);
    value = values;
    i = 0;
    while (value)
    {
        if (value->bExpr)
        {
            _gmx_selparser_error("expressions not supported within range parameters");
            return FALSE;
        }
        if (value->type != INT_VALUE)
        {
            gmx_bug("internal error");
            return FALSE;
        }
        /* Make sure the input range is in increasing order */
        if (value->u.i.i1 > value->u.i.i2)
        {
            int tmp       = value->u.i.i1;
            value->u.i.i1 = value->u.i.i2;
            value->u.i.i2 = tmp;
        }
        /* Check if the new range overlaps or extends the previous one */
        if (i > 0 && value->u.i.i1 <= data[i-1]+1 && value->u.i.i2 >= data[i-2]-1)
        {
            data[i-2] = min(data[i-2], value->u.i.i1);
            data[i-1] = max(data[i-1], value->u.i.i2);
        }
        else
        {
            data[i++] = value->u.i.i1;
            data[i++] = value->u.i.i2;
        }
        value = value->next;
    }
    n = i/2;
    /* Sort the ranges and merge consequent ones */
    qsort(data, n, 2*sizeof(int), &cmp_int_range);
    for (i = j = 2; i < 2*n; i += 2)
    {
        if (data[j-1]+1 >= data[i])
        {
            if (data[i+1] > data[j-1])
            {
                data[j-1] = data[i+1];
            }
        }
        else
        {
            data[j]   = data[i];
            data[j+1] = data[i+1];
            j += 2;
        }
    }
    n = j/2;
    /* Store the values */
    if (param->flags & SPAR_VARNUM)
    {
        srenew(data, j);
        param->val.nr  = n;
        _gmx_selvalue_setstore_alloc(&param->val, data, j);
    }
    else
    {
        if (n != param->val.nr)
        {
            _gmx_selparser_error("the value of parameter '%s' should consist of exactly one range",
                                 param->name);       
            sfree(data);
            return FALSE;
        }
        memcpy(param->val.u.i, data, 2*n*sizeof(int));
        sfree(data);
    }
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;

    return TRUE;
}

/*! \brief
 * Parses the values for a parameter that takes a variable number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_varnum(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param)
{
    t_selexpr_value    *value;
    int                 i, j;

    param->flags &= ~SPAR_DYNAMIC;
    /* Update nval if there are integer ranges. */
    if (param->val.type == INT_VALUE)
    {
        value = values;
        while (value)
        {
            if (value->type == INT_VALUE && !value->bExpr)
            {
                nval += abs(value->u.i.i2 - value->u.i.i1);
            }
            value = value->next;
        }
    }

    /* Check that the value type is actually implemented */
    if (param->val.type != INT_VALUE && param->val.type != REAL_VALUE
        && param->val.type != STR_VALUE && param->val.type != POS_VALUE)
    {
        gmx_bug("internal error");
        return FALSE;
    }

    /* Reserve appropriate amount of memory */
    if (param->val.type == POS_VALUE)
    {
        gmx_ana_pos_reserve(param->val.u.p, nval, 0);
        gmx_ana_indexmap_init(&param->val.u.p->m, NULL, NULL, INDEX_UNKNOWN);
    }
    else
    {
        _gmx_selvalue_reserve(&param->val, nval);
    }

    value = values;
    i     = 0;
    while (value)
    {
        if (value->bExpr)
        {
            _gmx_selparser_error("expressions not supported within value lists");
            return FALSE;
        }
        if (value->type != param->val.type)
        {
            gmx_bug("internal error");
            return FALSE;
        }
        switch (param->val.type)
        {
            case INT_VALUE:
                if (value->u.i.i1 <= value->u.i.i2)
                {
                    for (j = value->u.i.i1; j <= value->u.i.i2; ++j)
                    {
                        param->val.u.i[i++] = j;
                    }
                }
                else
                {
                    for (j = value->u.i.i1; j >= value->u.i.i2; --j)
                    {
                        param->val.u.i[i++] = j;
                    }
                }
                break;
            case REAL_VALUE: param->val.u.r[i++] = value->u.r;         break;
            case STR_VALUE:  param->val.u.s[i++] = strdup(value->u.s); break;
            case POS_VALUE:  copy_rvec(value->u.x, param->val.u.p->x[i++]); break;
            default: /* Should not be reached */
                gmx_bug("internal error");
                return FALSE;
        }
        value = value->next;
    }
    param->val.nr = i;
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;

    return TRUE;
}

/*! \brief
 * Adds a new subexpression reference to a selection element.
 *
 * \param[in,out] root  Root element to which the subexpression is added.
 * \param[in]     param Parameter for which this expression is a value.
 * \param[in]     expr  Expression to add.
 * \returns       The created child element.
 *
 * Creates a new \ref SEL_SUBEXPRREF element and adds it into the child
 * list of \p root.
 * If \p expr is already a \ref SEL_SUBEXPRREF, it is used as it is.
 * \ref SEL_ALLOCVAL is cleared for the returned element.
 */
static t_selelem *
add_child(t_selelem *root, gmx_ana_selparam_t *param, t_selelem *expr)
{
    gmx_ana_selparam_t *ps;
    int                 n;
    t_selelem          *child;
    int                 rc;

    if (root->type != SEL_EXPRESSION && root->type != SEL_MODIFIER)
    {
        gmx_bug("unsupported root element for selection parameter parser");
        return NULL;
    }
    ps = root->u.expr.method->param;
    n  = param - ps;
    /* Create a subexpression reference element if necessary */
    if (expr->type == SEL_SUBEXPRREF)
    {
        child = expr;
    }
    else
    {
        child = _gmx_selelem_create(SEL_SUBEXPRREF);
        if (!child)
        {
            return NULL;
        }
        _gmx_selelem_set_vtype(child, expr->v.type);
        child->child  = expr;
    }
    /* Put the child element in the correct place */
    if (!root->child || n < root->child->u.param - ps)
    {
        child->next = root->child;
        root->child = child;
    }
    else
    {
        t_selelem *prev;

        prev = root->child;
        while (prev->next && prev->next->u.param - ps >= n)
        {
            prev = prev->next;
        }
        child->next = prev->next;
        prev->next  = child;
    }
    /* Setup the child element */
    child->flags &= ~SEL_ALLOCVAL;
    child->u.param = param;
    if (child->v.type != param->val.type)
    {
        _gmx_selparser_error("invalid expression value for parameter '%s'",
                             param->name);
        goto on_error;
    }
    rc = _gmx_selelem_update_flags(child);
    if (rc != 0)
    {
        goto on_error;
    }
    if ((child->flags & SEL_DYNAMIC) && !(param->flags & SPAR_DYNAMIC))
    {
        _gmx_selparser_error("parameter '%s' does not support dynamic values",
                             param->name);
        goto on_error;
    }
    if (!(child->flags & SEL_DYNAMIC))
    {
        param->flags &= ~SPAR_DYNAMIC;
    }
    return child;

on_error:
    if (child != expr)
    {
        _gmx_selelem_free(child);
    }
    return NULL;
}

/*! \brief
 * Parses an expression value for a parameter that takes a variable number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 */
static bool
parse_values_varnum_expr(int nval, t_selexpr_value *values,
                         gmx_ana_selparam_t *param, t_selelem *root)
{
    t_selexpr_value    *value;
    t_selelem          *child;
    t_selelem          *expr;

    if (nval != 1 || !values->bExpr)
    {
        gmx_bug("internal error");
        return FALSE;
    }

    value = values;
    child = add_child(root, param, value->u.expr);
    value->u.expr = NULL;
    if (!child)
    {
        return FALSE;
    }

    /* Process single-valued expressions */
    /* TODO: We should also handle SEL_SINGLEVAL expressions here */
    if (child->v.type == POS_VALUE || child->v.type == GROUP_VALUE)
    {
        /* Set the value storage */
        _gmx_selvalue_setstore(&child->v, param->val.u.ptr);
        param->val.nr = 1;
        if (param->nvalptr)
        {
            *param->nvalptr = param->val.nr;
        }
        param->nvalptr = NULL;
        return TRUE;
    }

    if (!(child->flags & SEL_VARNUMVAL))
    {
        _gmx_selparser_error("invalid expression value for parameter '%s'",
                             param->name);
        return FALSE;
    }

    param->val.nr   = -1;
    *param->nvalptr = param->val.nr;
    /* Rest of the initialization is done during compilation in
     * init_method(). */

    return TRUE;
}

/*! \brief
 * Initializes the storage of an expression value.
 *
 * \param[in,out] sel   Selection element that evaluates the value.
 * \param[in]     param Parameter to receive the value.
 * \param[in]     i     The value of \p sel evaluates the value \p i for
 *   \p param.
 *
 * Initializes the data pointer of \p sel such that the result is stored
 * as the value \p i of \p param.
 * This function is used internally by parse_values_std().
 */
static bool
set_expr_value_store(t_selelem *sel, gmx_ana_selparam_t *param, int i)
{
    if (sel->v.type != GROUP_VALUE && !(sel->flags & SEL_SINGLEVAL))
    {
        _gmx_selparser_error("invalid expression value for parameter '%s'",
                             param->name);
        return FALSE;
    }
    switch (sel->v.type)
    {
        case INT_VALUE:   sel->v.u.i = &param->val.u.i[i]; break;
        case REAL_VALUE:  sel->v.u.r = &param->val.u.r[i]; break;
        case STR_VALUE:   sel->v.u.s = &param->val.u.s[i]; break;
        case POS_VALUE:   sel->v.u.p = &param->val.u.p[i]; break;
        case GROUP_VALUE: sel->v.u.g = &param->val.u.g[i]; break;
        default: /* Error */
            gmx_bug("internal error");
            return FALSE;
    }
    sel->v.nalloc = -1;
    return TRUE;
}

/*! \brief
 * Parses the values for a parameter that takes a constant number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_std(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param,
                 t_selelem *root)
{
    t_selexpr_value   *value;
    t_selelem         *child;
    int                i, j;
    bool               bDynamic;

    /* Handle atom-valued parameters */
    if (param->flags & SPAR_ATOMVAL)
    {
        if (nval > 1)
        {
            _gmx_selparser_error("extra values for parameter '%s'", param->name);
            return FALSE;
        }
        value = values;
        if (value->bExpr)
        {
            child = add_child(root, param, value->u.expr);
            value->u.expr = NULL;
            if (!child)
            {
                return FALSE;
            }
            if (child->v.type != GROUP_VALUE && (child->flags & SEL_ATOMVAL))
            {
                /* Rest of the initialization is done during compilation in
                 * init_method(). */
                /* TODO: Positions are not correctly handled */
                param->val.nr = -1;
                if (param->nvalptr)
                {
                    *param->nvalptr = -1;
                }
                return TRUE;
            }
            param->flags  &= ~SPAR_ATOMVAL;
            param->val.nr  = 1;
            if (param->nvalptr)
            {
                *param->nvalptr = 1;
            }
            param->nvalptr = NULL;
            if (param->val.type == INT_VALUE || param->val.type == REAL_VALUE
                || param->val.type == STR_VALUE)
            {
                _gmx_selvalue_reserve(&param->val, 1);
            }
            return set_expr_value_store(child, param, 0);
        }
        /* If we reach here, proceed with normal parameter handling */
        param->val.nr = 1;
        if (param->val.type == INT_VALUE || param->val.type == REAL_VALUE
            || param->val.type == STR_VALUE)
        {
            _gmx_selvalue_reserve(&param->val, 1);
        }
        param->flags &= ~SPAR_ATOMVAL;
        param->flags &= ~SPAR_DYNAMIC;
    }

    value = values;
    i = 0;
    bDynamic = FALSE;
    while (value && i < param->val.nr)
    {
        if (value->type != param->val.type)
        {
            _gmx_selparser_error("incorrect value for parameter '%s' skipped", param->name);
            value = value->next;
            continue;
        }
        if (value->bExpr)
        {
            child = add_child(root, param, value->u.expr);
            /* Clear the expression from the value once it is stored */
            value->u.expr = NULL;
            /* Check that the expression is valid */
            if (!child)
            {
                return FALSE;
            }
            if (!set_expr_value_store(child, param, i))
            {
                return FALSE;
            }
            if (child->flags & SEL_DYNAMIC)
            {
                bDynamic = TRUE;
            }
        }
        else
        {
            /* Value is not an expression */
            switch (value->type)
            {
                case INT_VALUE:
                    if (value->u.i.i1 <= value->u.i.i2)
                    {
                        for (j = value->u.i.i1; j <= value->u.i.i2 && i < param->val.nr; ++j)
                        {
                            param->val.u.i[i++] = j;
                        }
                        if (j != value->u.i.i2 + 1)
                        {
                            _gmx_selparser_error("extra values for parameter '%s' skipped", param->name);
                        }
                    }
                    else
                    {
                        for (j = value->u.i.i1; j >= value->u.i.i2 && i < param->val.nr; --j)
                        {
                            param->val.u.i[i++] = j;
                        }
                        if (j != value->u.i.i2 - 1)
                        {
                            _gmx_selparser_error("extra values for parameter '%s' skipped", param->name);
                        }
                    }
                    --i;
                    break;
                case REAL_VALUE:
                    param->val.u.r[i] = value->u.r;
                    break;
                case STR_VALUE:
                    param->val.u.s[i] = strdup(value->u.s);
                    break;
                case POS_VALUE:
                    gmx_ana_pos_init_const(&param->val.u.p[i], value->u.x);
                    break;
                case NO_VALUE:
                case GROUP_VALUE:
                    gmx_bug("internal error");
                    return FALSE;
            }
        }
        ++i;
        value = value->next;
    }
    if (value)
    {
        _gmx_selparser_error("extra values for parameter '%s'", param->name);
        return FALSE;
    }
    if (i < param->val.nr)
    {
        _gmx_selparser_error("not enough values for parameter '%s'", param->name);
        return FALSE;
    }
    if (!bDynamic)
    {
        param->flags &= ~SPAR_DYNAMIC;
    }
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;

    return TRUE;
}

/*! \brief
 * Parses the values for a boolean parameter.
 *
 * \param[in] name   Name by which the parameter was given.
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 */
static bool
parse_values_bool(const char *name, int nval, t_selexpr_value *values, gmx_ana_selparam_t *param)
{
    bool bSetNo;
    int  len;

    if (param->val.type != NO_VALUE)
    {
        gmx_bug("internal error");
        return FALSE;
    }
    if (nval > 1 || (values && values->type != INT_VALUE))
    {
        /* We should be here only if the parser has screwed up */
        _gmx_selparser_error("boolean parameter '%s' should have at most one integer value", param->name);
        return FALSE;
    }

    bSetNo = FALSE;
    /* Check if the parameter name is given with a 'no' prefix */
    len = strlen(name);
    if (len > 2 && name[0] == 'n' && name[1] == 'o'
        && strncmp(name+2, param->name, len-2) == 0)
    {
        bSetNo = TRUE;
    }
    if (bSetNo && nval > 0)
    {
        /* We should be here only if the parser has screwed up */
        _gmx_selparser_error("boolean parameter 'no%s' should not have a value", param->name);
        return FALSE;
    }
    if (values && values->u.i.i1 == 0)
    {
        bSetNo = TRUE;
    }

    *param->val.u.b = bSetNo ? FALSE : TRUE;
    return TRUE;
}

/*! \brief
 * Parses the values for an enumeration parameter.
 *
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \returns   TRUE if the values were parsed successfully, FALSE otherwise.
 */
static bool
parse_values_enum(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param)
{
    int  i, len, match;

    if (nval != 1)
    {
        _gmx_selparser_error("a single value is required for parameter '%s'", param->name);
        return FALSE;
    }
    if (values->type != STR_VALUE || param->val.type != STR_VALUE)
    {
        gmx_bug("internal error");
        return FALSE;
    }
    if (values->bExpr)
    {
        _gmx_selparser_error("expression value for enumerated parameter '%s' not supported", param->name);
        return FALSE;
    }

    len = strlen(values->u.s);
    i = 1;
    match = 0;
    while (param->val.u.s[i] != NULL)
    {
        if (strncmp(values->u.s, param->val.u.s[i], len) == 0)
        {
            /* Check if there is a duplicate match */
            if (match > 0)
            {
                _gmx_selparser_error("ambiguous value for parameter '%s'", param->name);
                return FALSE;
            }
            match = i;
        }
        ++i;
    }
    if (match == 0)
    {
        _gmx_selparser_error("invalid value for parameter '%s'", param->name);
        return FALSE;
    }
    param->val.u.s[0] = param->val.u.s[match];
    return TRUE;
}

/*! \brief
 * Replaces constant expressions with their values.
 *
 * \param[in,out] values First element in the value list to process.
 */
static void
convert_const_values(t_selexpr_value *values)
{
    t_selexpr_value *val;

    val = values;
    while (val)
    {
        if (val->bExpr && val->u.expr->v.type != GROUP_VALUE &&
            val->u.expr->type == SEL_CONST)
        {
            t_selelem *expr = val->u.expr;
            val->bExpr = FALSE;
            switch (expr->v.type)
            {
                case INT_VALUE:
                    val->u.i.i1 = val->u.i.i2 = expr->v.u.i[0];
                    break;
                case REAL_VALUE:
                    val->u.r = expr->v.u.r[0];
                    break;
                case POS_VALUE:
                    copy_rvec(expr->v.u.p->x[0], val->u.x);
                    break;
                default:
                    gmx_bug("internal error");
                    break;
            }
            _gmx_selelem_free(expr);
        }
        val = val->next;
    }
}

/*!
 * \param     pparams List of parameters from the selection parser.
 * \param[in] nparam  Number of parameters in \p params.
 * \param     params  Array of parameters to parse.
 * \param     root    Selection element to which child expressions are added.
 * \returns   TRUE if the parameters were parsed successfully, FALSE otherwise.
 *
 * Initializes the \p params array based on the parameters in \p pparams.
 * See the documentation of \c gmx_ana_selparam_t for different options
 * available for parsing.
 *
 * The list \p pparams and any associated values are freed after the parameters
 * have been processed, no matter is there was an error or not.
 */
bool
_gmx_sel_parse_params(t_selexpr_param *pparams, int nparam, gmx_ana_selparam_t *params,
                      t_selelem *root)
{
    t_selexpr_param    *pparam;
    gmx_ana_selparam_t *oparam;
    bool                bOk, rc;
    int                 i;

    /* Check that the value pointers of SPAR_VARNUM parameters are NULL and
     * that they are not NULL for other parameters */
    bOk = TRUE;
    for (i = 0; i < nparam; ++i)
    {
        if (params[i].val.type != POS_VALUE && (params[i].flags & (SPAR_VARNUM | SPAR_ATOMVAL)))
        {
            if (params[i].val.u.ptr != NULL)
            {
                _gmx_selparser_error("warning: value pointer of parameter '%s' is not NULL\n"
                                     "         although it should be for SPAR_VARNUM and SPAR_ATOMVAL parameters\n",
                                     params[i].name);
            }
            if ((params[i].flags & SPAR_VARNUM)
                && (params[i].flags & SPAR_DYNAMIC) && !params[i].nvalptr)
            {
                _gmx_selparser_error("error: nvalptr of parameter '%s' is NULL\n"
                                     "       but both SPAR_VARNUM and SPAR_DYNAMIC are specified\n",
                                     params[i].name);
                bOk = FALSE;
            }
        }
        else
        {
            if (params[i].val.u.ptr == NULL)
            {
                _gmx_selparser_error("error: value pointer of parameter '%s' is NULL\n",
                                     params[i].name);
                bOk = FALSE;
            }
        }
    }
    if (!bOk)
    {
        return FALSE;
    }
    /* Parse the parameters */
    pparam = pparams;
    i      = 0;
    while (pparam)
    {
        /* Find the parameter and make some checks */
        if (pparam->name != NULL)
        {
            i = -1;
            oparam = gmx_ana_selparam_find(pparam->name, nparam, params);
        }
        else if (i >= 0)
        {
            oparam = &params[i];
            if (oparam->name != NULL)
            {
                oparam = NULL;
                _gmx_selparser_error("too many NULL parameters provided");
                bOk = FALSE;
                pparam = pparam->next;
                continue;
            }
            ++i;
        }
        else
        {
            _gmx_selparser_error("all NULL parameters should appear in the beginning of the list");
            bOk = FALSE;
            pparam = pparam->next;
            continue;
        }
        if (!oparam)
        {
            _gmx_selparser_error("unknown parameter '%s' skipped", pparam->name);
            bOk = FALSE;
            pparam = pparam->next;
            continue;
        }
        if (oparam->flags & SPAR_SET)
        {
            _gmx_selparser_error("parameter '%s' set multiple times, extra values skipped", pparam->name);
            bOk = FALSE;
            pparam = pparam->next;
            continue;
        }
        oparam->flags |= SPAR_SET;
        /* Process the values for the parameter */
        convert_const_values(pparam->value);
        if (oparam->val.type == NO_VALUE)
        {
            rc = parse_values_bool(pparam->name, pparam->nval, pparam->value, oparam);
        }
        else if (oparam->flags & SPAR_RANGES)
        {
            rc = parse_values_range(pparam->nval, pparam->value, oparam);
        }
        else if (oparam->flags & SPAR_VARNUM)
        {
            if (pparam->nval == 1 && pparam->value->bExpr)
            {
                rc = parse_values_varnum_expr(pparam->nval, pparam->value, oparam, root);
            }
            else
            {
                rc = parse_values_varnum(pparam->nval, pparam->value, oparam);
            }
        }
        else if (oparam->flags & SPAR_ENUMVAL)
        {
            rc = parse_values_enum(pparam->nval, pparam->value, oparam);
        }
        else
        {
            rc = parse_values_std(pparam->nval, pparam->value, oparam, root);
        }
        if (!rc)
        {
            bOk = FALSE;
        }
        /* Advance to the next parameter */
        pparam = pparam->next;
    }
    /* Check that all required parameters are present */
    for (i = 0; i < nparam; ++i)
    {
        if (!(params[i].flags & SPAR_OPTIONAL) && !(params[i].flags & SPAR_SET))
        {
            _gmx_selparser_error("required parameter '%s' not specified", params[i].name);
            bOk = FALSE;
        }
    }
    
    _gmx_selexpr_free_params(pparams);
    return bOk;
}
