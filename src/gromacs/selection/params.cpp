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
 * Implements functions in selparam.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>

#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/string2.h"
#include "gromacs/legacyheaders/vec.h"

#include "gromacs/selection/position.h"
#include "gromacs/selection/selmethod.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/utility/errorcodes.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/messagestringcollector.h"
#include "gromacs/utility/stringutil.h"

#include "parsetree.h"
#include "position.h"
#include "scanner.h"
#include "selelem.h"

using std::min;
using std::max;

using gmx::SelectionTreeElement;
using gmx::SelectionTreeElementPointer;

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
 * Does a type conversion on a \c t_selexpr_value.
 *
 * \param[in,out] value    Value to convert.
 * \param[in]     type     Type to convert to.
 * \param[in]     scanner  Scanner data structure.
 * \returns       0 on success, a non-zero value on error.
 */
static int
convert_value(t_selexpr_value *value, e_selvalue_t type, void *scanner)
{
    if (value->type == type || type == NO_VALUE)
    {
        return 0;
    }
    if (value->hasExpressionValue())
    {
        /* Conversion from atom selection to position using default
         * reference positions. */
        if (value->type == GROUP_VALUE && type == POS_VALUE)
        {
            value->expr =
                _gmx_sel_init_position(value->expr, NULL, scanner);
            // FIXME: Use exceptions
            if (!value->expr)
            {
                return -1;
            }
            value->type = type;
            return 0;
        }
        return -1;
    }
    else
    {
        /* Integers to floating point are easy */
        if (value->type == INT_VALUE && type == REAL_VALUE)
        {
            real r1 = (real)value->u.i.i1;
            real r2 = (real)value->u.i.i2;
            value->u.r.r1 = r1;
            value->u.r.r2 = r2;
            value->type = type;
            return 0;
        }
        /* Reals that are integer-valued can also be converted */
        if (value->type == REAL_VALUE && type == INT_VALUE
            && gmx_within_tol(value->u.r.r1, (int)value->u.r.r1, GMX_REAL_EPS)
            && gmx_within_tol(value->u.r.r2, (int)value->u.r.r2, GMX_REAL_EPS))
        {
            int i1 = (int)value->u.r.r1;
            int i2 = (int)value->u.r.r2;
            value->u.i.i1 = i1;
            value->u.i.i2 = i2;
            value->type = type;
            return 0;
        }
    }
    return -1;
}

/*! \brief
 * Does a type conversion on a list of values.
 *
 * \param[in,out] values   Values to convert.
 * \param[in]     type     Type to convert to.
 * \param[in]     scanner  Scanner data structure.
 * \returns       0 on success, a non-zero value on error.
 */
static int
convert_values(t_selexpr_value *values, e_selvalue_t type, void *scanner)
{
    t_selexpr_value *value;
    int              rc, rc1;

    rc = 0;
    value = values;
    while (value)
    {
        rc1 = convert_value(value, type, scanner);
        if (rc1 != 0 && rc == 0)
        {
            rc = rc1;
        }
        value = value->next;
    }
    /* FIXME: More informative error messages */
    return rc;
}

/*! \brief
 * Adds a child element for a parameter, keeping the parameter order.
 *
 * \param[in,out] root  Root element to which the child is added.
 * \param[in]     child Child to add.
 * \param[in]     param Parameter for which this child is a value.
 *
 * Puts \p child in the child list of \p root such that the list remains
 * in the same order as the corresponding parameters.
 */
static void
place_child(const SelectionTreeElementPointer &root,
            const SelectionTreeElementPointer &child,
            gmx_ana_selparam_t *param)
{
    gmx_ana_selparam_t *ps;
    int                 n;

    ps = root->u.expr.method->param;
    n  = param - ps;
    /* Put the child element in the correct place */
    if (!root->child || n < root->child->u.param - ps)
    {
        child->next = root->child;
        root->child = child;
    }
    else
    {
        SelectionTreeElementPointer prev = root->child;
        while (prev->next && prev->next->u.param - ps >= n)
        {
            prev = prev->next;
        }
        child->next = prev->next;
        prev->next  = child;
    }
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
 * Comparison function for sorting real ranges.
 *
 * \param[in] a Pointer to the first range.
 * \param[in] b Pointer to the second range.
 * \returns   -1, 0, or 1 depending on the relative order of \p a and \p b.
 *
 * The ranges are primarily sorted based on their starting point, and
 * secondarily based on length (longer ranges come first).
 */
static int
cmp_real_range(const void *a, const void *b)
{
    if (((real *)a)[0] < ((real *)b)[0])
    {
        return -1;
    }
    if (((real *)a)[0] > ((real *)b)[0])
    {
        return 1;
    }
    if (((real *)a)[1] > ((real *)b)[1])
    {
        return -1;
    }
    return 0;
}

/*! \brief
 * Parses the values for a parameter that takes integer or real ranges.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_range(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param,
                   void *scanner)
{
    t_selexpr_value    *value;
    int                *idata;
    real               *rdata;
    int                 i, j, n;

    param->flags &= ~SPAR_DYNAMIC;
    if (param->val.type != INT_VALUE && param->val.type != REAL_VALUE)
    {
        GMX_ERROR_NORET(gmx::eeInternalError, "Invalid range parameter type");
        return false;
    }
    idata = NULL;
    rdata = NULL;
    if (param->val.type == INT_VALUE)
    {
        snew(idata, nval*2);
    }
    else
    {
        snew(rdata, nval*2);
    }
    value = values;
    i = 0;
    while (value)
    {
        if (value->hasExpressionValue())
        {
            _gmx_selparser_error(scanner, "expressions not supported within range parameters");
            return false;
        }
        if (value->type != param->val.type)
        {
            GMX_ERROR_NORET(gmx::eeInternalError, "Invalid range value type");
            return false;
        }
        if (param->val.type == INT_VALUE)
        {
            /* Make sure the input range is in increasing order */
            if (value->u.i.i1 > value->u.i.i2)
            {
                int tmp       = value->u.i.i1;
                value->u.i.i1 = value->u.i.i2;
                value->u.i.i2 = tmp;
            }
            /* Check if the new range overlaps or extends the previous one */
            if (i > 0 && value->u.i.i1 <= idata[i-1]+1 && value->u.i.i2 >= idata[i-2]-1)
            {
                idata[i-2] = min(idata[i-2], value->u.i.i1);
                idata[i-1] = max(idata[i-1], value->u.i.i2);
            }
            else
            {
                idata[i++] = value->u.i.i1;
                idata[i++] = value->u.i.i2;
            }
        }
        else
        {
            /* Make sure the input range is in increasing order */
            if (value->u.r.r1 > value->u.r.r2)
            {
                real tmp      = value->u.r.r1;
                value->u.r.r1 = value->u.r.r2;
                value->u.r.r2 = tmp;
            }
            /* Check if the new range overlaps or extends the previous one */
            if (i > 0 && value->u.r.r1 <= rdata[i-1] && value->u.r.r2 >= rdata[i-2])
            {
                rdata[i-2] = min(rdata[i-2], value->u.r.r1);
                rdata[i-1] = max(rdata[i-1], value->u.r.r2);
            }
            else
            {
                rdata[i++] = value->u.r.r1;
                rdata[i++] = value->u.r.r2;
            }
        }
        value = value->next;
    }
    n = i/2;
    /* Sort the ranges and merge consequent ones */
    if (param->val.type == INT_VALUE)
    {
        qsort(idata, n, 2*sizeof(int), &cmp_int_range);
        for (i = j = 2; i < 2*n; i += 2)
        {
            if (idata[j-1]+1 >= idata[i])
            {
                if (idata[i+1] > idata[j-1])
                {
                    idata[j-1] = idata[i+1];
                }
            }
            else
            {
                idata[j]   = idata[i];
                idata[j+1] = idata[i+1];
                j += 2;
            }
        }
    }
    else
    {
        qsort(rdata, n, 2*sizeof(real), &cmp_real_range);
        for (i = j = 2; i < 2*n; i += 2)
        {
            if (rdata[j-1]+1 >= rdata[i])
            {
                if (rdata[i+1] > rdata[j-1])
                {
                    rdata[j-1] = rdata[i+1];
                }
            }
            else
            {
                rdata[j]   = rdata[i];
                rdata[j+1] = rdata[i+1];
                j += 2;
            }
        }
    }
    n = j/2;
    /* Store the values */
    if (param->flags & SPAR_VARNUM)
    {
        param->val.nr  = n;
        if (param->val.type == INT_VALUE)
        {
            srenew(idata, j);
            _gmx_selvalue_setstore_alloc(&param->val, idata, j);
        }
        else
        {
            srenew(rdata, j);
            _gmx_selvalue_setstore_alloc(&param->val, rdata, j);
        }
    }
    else
    {
        if (n != param->val.nr)
        {
            _gmx_selparser_error(scanner, "the value should consist of exactly one range");
            sfree(idata);
            sfree(rdata);
            return false;
        }
        if (param->val.type == INT_VALUE)
        {
            memcpy(param->val.u.i, idata, 2*n*sizeof(int));
            sfree(idata);
        }
        else
        {
            memcpy(param->val.u.r, rdata, 2*n*sizeof(real));
            sfree(rdata);
        }
    }
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;

    return true;
}

/*! \brief
 * Parses the values for a parameter that takes a variable number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_varnum(int nval, t_selexpr_value *values,
                    gmx_ana_selparam_t *param,
                    const SelectionTreeElementPointer &root,
                    void *scanner)
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
            if (value->type == INT_VALUE && !value->hasExpressionValue())
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
        GMX_ERROR_NORET(gmx::eeInternalError,
                        "Variable-count value type not implemented");
        return false;
    }

    /* Reserve appropriate amount of memory */
    if (param->val.type == POS_VALUE)
    {
        gmx_ana_pos_reserve(param->val.u.p, nval, 0);
        gmx_ana_pos_set_nr(param->val.u.p, nval);
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
        if (value->hasExpressionValue())
        {
            _gmx_selparser_error(scanner, "expressions not supported within value lists");
            return false;
        }
        if (value->type != param->val.type)
        {
            GMX_ERROR_NORET(gmx::eeInternalError, "Invalid value type");
            return false;
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
            case REAL_VALUE:
                if (value->u.r.r1 != value->u.r.r2)
                {
                    _gmx_selparser_error(scanner, "real ranges not supported");
                    return false;
                }
                param->val.u.r[i++] = value->u.r.r1;
                break;
            case STR_VALUE:  param->val.u.s[i++] = strdup(value->u.s); break;
            case POS_VALUE:  copy_rvec(value->u.x, param->val.u.p->x[i++]); break;
            default: /* Should not be reached */
                GMX_ERROR_NORET(gmx::eeInternalError, "Invalid value type");
                return false;
        }
        value = value->next;
    }
    param->val.nr = i;
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;
    /* Create a dummy child element to store the string values.
     * This element is responsible for freeing the values, but carries no
     * other function. */
    if (param->val.type == STR_VALUE)
    {
        SelectionTreeElementPointer child(new SelectionTreeElement(SEL_CONST));
        _gmx_selelem_set_vtype(child, STR_VALUE);
        child->setName(param->name);
        child->flags &= ~SEL_ALLOCVAL;
        child->flags |= SEL_FLAGSSET | SEL_VARNUMVAL | SEL_ALLOCDATA;
        child->v.nr = param->val.nr;
        _gmx_selvalue_setstore(&child->v, param->val.u.s);
        /* Because the child is not group-valued, the u union is not used
         * for anything, so we can abuse it by storing the parameter value
         * as place_child() expects, but this is really ugly... */
        child->u.param = param;
        place_child(root, child, param);
    }

    return true;
}

/*! \brief
 * Adds a new subexpression reference to a selection element.
 *
 * \param[in,out] root  Root element to which the subexpression is added.
 * \param[in]     param Parameter for which this expression is a value.
 * \param[in]     expr  Expression to add.
 * \param[in]     scanner Scanner data structure.
 * \returns       The created child element.
 *
 * Creates a new \ref SEL_SUBEXPRREF element and adds it into the child
 * list of \p root.
 * If \p expr is already a \ref SEL_SUBEXPRREF, it is used as it is.
 * \ref SEL_ALLOCVAL is cleared for the returned element.
 */
static SelectionTreeElementPointer
add_child(const SelectionTreeElementPointer &root, gmx_ana_selparam_t *param,
          const SelectionTreeElementPointer &expr, void *scanner)
{
    GMX_RELEASE_ASSERT(root->type == SEL_EXPRESSION || root->type == SEL_MODIFIER,
                       "Unsupported root element for selection parameter parser");
    SelectionTreeElementPointer child;
    /* Create a subexpression reference element if necessary */
    if (expr->type == SEL_SUBEXPRREF)
    {
        child = expr;
    }
    else
    {
        child.reset(new SelectionTreeElement(SEL_SUBEXPRREF));
        _gmx_selelem_set_vtype(child, expr->v.type);
        child->child  = expr;
    }
    /* Setup the child element */
    child->flags &= ~SEL_ALLOCVAL;
    child->u.param = param;
    if (child->v.type != param->val.type)
    {
        _gmx_selparser_error(scanner, "invalid expression value");
        // FIXME: Use exceptions.
        return SelectionTreeElementPointer();
    }
    _gmx_selelem_update_flags(child, scanner);
    if ((child->flags & SEL_DYNAMIC) && !(param->flags & SPAR_DYNAMIC))
    {
        _gmx_selparser_error(scanner, "dynamic values not supported");
        // FIXME: Use exceptions.
        return SelectionTreeElementPointer();
    }
    if (!(child->flags & SEL_DYNAMIC))
    {
        param->flags &= ~SPAR_DYNAMIC;
    }
    /* Put the child element in the correct place */
    place_child(root, child, param);
    return child;
}

/*! \brief
 * Parses an expression value for a parameter that takes a variable number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_varnum_expr(int nval, t_selexpr_value *values,
                         gmx_ana_selparam_t *param,
                         const SelectionTreeElementPointer &root,
                         void *scanner)
{
    if (nval != 1 || !values->hasExpressionValue())
    {
        GMX_ERROR_NORET(gmx::eeInternalError, "Invalid expression value");
        return false;
    }

    t_selexpr_value *value = values;
    SelectionTreeElementPointer child
        = add_child(root, param, value->expr, scanner);
    if (!child)
    {
        return false;
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
        return true;
    }

    if (!(child->flags & SEL_VARNUMVAL))
    {
        _gmx_selparser_error(scanner, "invalid expression value");
        return false;
    }

    child->flags   |= SEL_ALLOCVAL;
    param->val.nr   = -1;
    *param->nvalptr = param->val.nr;
    /* Rest of the initialization is done during compilation in
     * init_method(). */

    return true;
}

/*! \brief
 * Initializes the storage of an expression value.
 *
 * \param[in,out] sel   Selection element that evaluates the value.
 * \param[in]     param Parameter to receive the value.
 * \param[in]     i     The value of \p sel evaluates the value \p i for
 *   \p param.
 * \param[in]     scanner Scanner data structure.
 *
 * Initializes the data pointer of \p sel such that the result is stored
 * as the value \p i of \p param.
 * This function is used internally by parse_values_std().
 */
static bool
set_expr_value_store(const SelectionTreeElementPointer &sel,
                     gmx_ana_selparam_t *param, int i, void *scanner)
{
    if (sel->v.type != GROUP_VALUE && !(sel->flags & SEL_SINGLEVAL))
    {
        _gmx_selparser_error(scanner, "invalid expression value");
        return false;
    }
    switch (sel->v.type)
    {
        case INT_VALUE:   sel->v.u.i = &param->val.u.i[i]; break;
        case REAL_VALUE:  sel->v.u.r = &param->val.u.r[i]; break;
        case STR_VALUE:   sel->v.u.s = &param->val.u.s[i]; break;
        case POS_VALUE:   sel->v.u.p = &param->val.u.p[i]; break;
        case GROUP_VALUE: sel->v.u.g = &param->val.u.g[i]; break;
        default: /* Error */
            GMX_ERROR_NORET(gmx::eeInternalError, "Invalid value type");
            return false;
    }
    sel->v.nr = 1;
    sel->v.nalloc = -1;
    return true;
}

/*! \brief
 * Parses the values for a parameter that takes a constant number of values.
 * 
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_std(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param,
                 const SelectionTreeElementPointer &root, void *scanner)
{
    t_selexpr_value   *value;
    int                i, j;
    bool               bDynamic;

    /* Handle atom-valued parameters */
    if (param->flags & SPAR_ATOMVAL)
    {
        if (nval > 1)
        {
            _gmx_selparser_error(scanner, "more than one value not supported");
            return false;
        }
        value = values;
        if (value->hasExpressionValue())
        {
            SelectionTreeElementPointer child
                = add_child(root, param, value->expr, scanner);
            if (!child)
            {
                return false;
            }
            child->flags |= SEL_ALLOCVAL;
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
                return true;
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
            return set_expr_value_store(child, param, 0, scanner);
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
    bDynamic = false;
    while (value && i < param->val.nr)
    {
        if (value->type != param->val.type)
        {
            _gmx_selparser_error(scanner, "incorrect value skipped");
            value = value->next;
            continue;
        }
        if (value->hasExpressionValue())
        {
            SelectionTreeElementPointer child
                = add_child(root, param, value->expr, scanner);
            /* Check that the expression is valid */
            if (!child)
            {
                return false;
            }
            if (!set_expr_value_store(child, param, i, scanner))
            {
                return false;
            }
            if (child->flags & SEL_DYNAMIC)
            {
                bDynamic = true;
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
                            _gmx_selparser_error(scanner, "extra values skipped");
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
                            _gmx_selparser_error(scanner, "extra values skipped");
                        }
                    }
                    --i;
                    break;
                case REAL_VALUE:
                    if (value->u.r.r1 != value->u.r.r2)
                    {
                        _gmx_selparser_error(scanner, "real ranges not supported");
                        return false;
                    }
                    param->val.u.r[i] = value->u.r.r1;
                    break;
                case STR_VALUE:
                    param->val.u.s[i] = strdup(value->u.s);
                    break;
                case POS_VALUE:
                    gmx_ana_pos_init_const(&param->val.u.p[i], value->u.x);
                    break;
                case NO_VALUE:
                case GROUP_VALUE:
                    GMX_ERROR_NORET(gmx::eeInternalError,
                                    "Invalid non-expression value");
                    return false;
            }
        }
        ++i;
        value = value->next;
    }
    if (value != NULL)
    {
        _gmx_selparser_error(scanner, "extra values'");
        return false;
    }
    if (i < param->val.nr)
    {
        _gmx_selparser_error(scanner, "not enough values");
        return false;
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

    return true;
}

/*! \brief
 * Parses the values for a boolean parameter.
 *
 * \param[in] name   Name by which the parameter was given.
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_bool(const char *name, int nval, t_selexpr_value *values,
                  gmx_ana_selparam_t *param, void *scanner)
{
    bool bSetNo;
    int  len;

    if (param->val.type != NO_VALUE)
    {
        GMX_ERROR_NORET(gmx::eeInternalError, "Invalid boolean parameter");
        return false;
    }
    if (nval > 1 || (values && values->type != INT_VALUE))
    {
        _gmx_selparser_error(scanner, "parameter takes only a yes/no/on/off/0/1 value");
        return false;
    }

    bSetNo = false;
    /* Check if the parameter name is given with a 'no' prefix */
    len = name != NULL ? strlen(name) : 0;
    if (len > 2 && name[0] == 'n' && name[1] == 'o'
        && strncmp(name+2, param->name, len-2) == 0)
    {
        bSetNo = true;
    }
    if (bSetNo && nval > 0)
    {
        _gmx_selparser_error(scanner, "parameter 'no%s' should not have a value",
                             param->name);
        return false;
    }
    if (values && values->u.i.i1 == 0)
    {
        bSetNo = true;
    }

    *param->val.u.b = bSetNo ? false : true;
    return true;
}

/*! \brief
 * Parses the values for an enumeration parameter.
 *
 * \param[in] nval   Number of values in \p values.
 * \param[in] values Pointer to the list of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_enum(int nval, t_selexpr_value *values, gmx_ana_selparam_t *param,
                  void *scanner)
{
    int  i, len, match;

    if (nval != 1)
    {
        _gmx_selparser_error(scanner, "a single value is required");
        return false;
    }
    if (values->type != STR_VALUE || param->val.type != STR_VALUE)
    {
        GMX_ERROR_NORET(gmx::eeInternalError, "Invalid enum parameter");
        return false;
    }
    if (values->hasExpressionValue())
    {
        _gmx_selparser_error(scanner, "expression value for enumerated parameter not supported");
        return false;
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
                _gmx_selparser_error(scanner, "ambiguous value");
                return false;
            }
            match = i;
        }
        ++i;
    }
    if (match == 0)
    {
        _gmx_selparser_error(scanner, "invalid value");
        return false;
    }
    param->val.u.s[0] = param->val.u.s[match];
    return true;
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
        if (val->hasExpressionValue() && val->expr->v.type != GROUP_VALUE &&
            val->expr->type == SEL_CONST)
        {
            SelectionTreeElementPointer expr = val->expr;
            val->expr.reset();
            switch (expr->v.type)
            {
                case INT_VALUE:
                    val->u.i.i1 = val->u.i.i2 = expr->v.u.i[0];
                    break;
                case REAL_VALUE:
                    val->u.r.r1 = val->u.r.r2 = expr->v.u.r[0];
                    break;
                case STR_VALUE:
                    val->u.s = expr->v.u.s[0];
                    break;
                case POS_VALUE:
                    copy_rvec(expr->v.u.p->x[0], val->u.x);
                    break;
                default:
                    GMX_ERROR_NORET(gmx::eeInternalError,
                                    "Unsupported value type");
                    break;
            }
        }
        val = val->next;
    }
}

/*!
 * \param     pparams List of parameters from the selection parser.
 * \param[in] nparam  Number of parameters in \p params.
 * \param     params  Array of parameters to parse.
 * \param     root    Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the parameters were parsed successfully, false otherwise.
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
                      const SelectionTreeElementPointer &root, void *scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    t_selexpr_param    *pparam;
    gmx_ana_selparam_t *oparam;
    bool                bOk, rc;
    int                 i;

    /* Check that the value pointers of SPAR_VARNUM parameters are NULL and
     * that they are not NULL for other parameters */
    bOk = true;
    for (i = 0; i < nparam; ++i)
    {
        std::string contextStr = gmx::formatString("In parameter '%s'", params[i].name);
        gmx::MessageStringContext  context(errors, contextStr);
        if (params[i].val.type != POS_VALUE && (params[i].flags & (SPAR_VARNUM | SPAR_ATOMVAL)))
        {
            if (params[i].val.u.ptr != NULL)
            {
                _gmx_selparser_error(scanner, "value pointer is not NULL "
                                     "although it should be for SPAR_VARNUM "
                                     "and SPAR_ATOMVAL parameters");
            }
            if ((params[i].flags & SPAR_VARNUM)
                && (params[i].flags & SPAR_DYNAMIC) && !params[i].nvalptr)
            {
                _gmx_selparser_error(scanner, "nvalptr is NULL but both "
                                     "SPAR_VARNUM and SPAR_DYNAMIC are specified");
                bOk = false;
            }
        }
        else
        {
            if (params[i].val.u.ptr == NULL)
            {
                _gmx_selparser_error(scanner, "value pointer is NULL");
                bOk = false;
            }
        }
    }
    if (!bOk)
    {
        _gmx_selexpr_free_params(pparams);
        return false;
    }
    /* Parse the parameters */
    pparam = pparams;
    i      = 0;
    while (pparam)
    {
        std::string contextStr;
        /* Find the parameter and make some checks */
        if (pparam->name != NULL)
        {
            contextStr = gmx::formatString("In parameter '%s'", pparam->name);
            i = -1;
            oparam = gmx_ana_selparam_find(pparam->name, nparam, params);
        }
        else if (i >= 0)
        {
            contextStr = gmx::formatString("In value %d", i + 1);
            oparam = &params[i];
            if (oparam->name != NULL)
            {
                oparam = NULL;
                _gmx_selparser_error(scanner, "too many NULL parameters provided");
                bOk = false;
                pparam = pparam->next;
                continue;
            }
            ++i;
        }
        else
        {
            _gmx_selparser_error(scanner, "all NULL parameters should appear in the beginning of the list");
            bOk = false;
            pparam = pparam->next;
            continue;
        }
        gmx::MessageStringContext  context(errors, contextStr);
        if (!oparam)
        {
            _gmx_selparser_error(scanner, "unknown parameter skipped");
            bOk = false;
            goto next_param;
        }
        if (oparam->val.type != NO_VALUE && pparam->value == NULL)
        {
            GMX_ASSERT(pparam->nval == 0, "Inconsistent values and value count");
            _gmx_selparser_error(scanner, "no value provided");
            bOk = false;
            goto next_param;
        }
        if (oparam->flags & SPAR_SET)
        {
            _gmx_selparser_error(scanner, "parameter set multiple times, extra values skipped");
            bOk = false;
            goto next_param;
        }
        oparam->flags |= SPAR_SET;
        /* Process the values for the parameter */
        convert_const_values(pparam->value);
        if (convert_values(pparam->value, oparam->val.type, scanner) != 0)
        {
            _gmx_selparser_error(scanner, "invalid value");
            bOk = false;
            goto next_param;
        }
        if (oparam->val.type == NO_VALUE)
        {
            rc = parse_values_bool(pparam->name, pparam->nval, pparam->value, oparam, scanner);
        }
        else if (oparam->flags & SPAR_RANGES)
        {
            rc = parse_values_range(pparam->nval, pparam->value, oparam, scanner);
        }
        else if (oparam->flags & SPAR_VARNUM)
        {
            if (pparam->nval == 1 && pparam->value->hasExpressionValue())
            {
                rc = parse_values_varnum_expr(pparam->nval, pparam->value, oparam, root, scanner);
            }
            else
            {
                rc = parse_values_varnum(pparam->nval, pparam->value, oparam, root, scanner);
            }
        }
        else if (oparam->flags & SPAR_ENUMVAL)
        {
            rc = parse_values_enum(pparam->nval, pparam->value, oparam, scanner);
        }
        else
        {
            rc = parse_values_std(pparam->nval, pparam->value, oparam, root, scanner);
        }
        if (!rc)
        {
            bOk = false;
        }
        /* Advance to the next parameter */
next_param:
        pparam = pparam->next;
    }
    /* Check that all required parameters are present */
    for (i = 0; i < nparam; ++i)
    {
        if (!(params[i].flags & SPAR_OPTIONAL) && !(params[i].flags & SPAR_SET))
        {
            _gmx_selparser_error(scanner, "required parameter '%s' not specified", params[i].name);
            bOk = false;
        }
    }

    _gmx_selexpr_free_params(pparams);
    return bOk;
}
