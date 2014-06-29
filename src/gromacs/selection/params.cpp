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
 * Implements functions in selparam.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include <algorithm>
#include <string>

#include "gromacs/legacyheaders/vec.h"

#include "gromacs/selection/position.h"
#include "gromacs/selection/selmethod.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/messagestringcollector.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "parsetree.h"
#include "position.h"
#include "scanner.h"
#include "selelem.h"

using gmx::SelectionParserValue;
using gmx::SelectionParserValueList;
using gmx::SelectionParserParameter;
using gmx::SelectionParserParameterList;
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
    for (; i < nparam; ++i)
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
 * Does a type conversion on a SelectionParserValue.
 *
 * \param[in,out] value    Value to convert.
 * \param[in]     type     Type to convert to.
 * \param[in]     scanner  Scanner data structure.
 * \returns       0 on success, a non-zero value on error.
 */
static int
convert_value(SelectionParserValue *value, e_selvalue_t type, void *scanner)
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
            SelectionTreeElementPointer expr =
                _gmx_sel_init_position(value->expr, NULL, scanner);
            // FIXME: Use exceptions
            if (!expr)
            {
                return -1;
            }
            *value = SelectionParserValue::createExpr(expr);
            return 0;
        }
        return -1;
    }
    else
    {
        /* Integers to floating point are easy */
        if (value->type == INT_VALUE && type == REAL_VALUE)
        {
            *value = SelectionParserValue::createRealRange(value->u.i.i1,
                                                           value->u.i.i2);
            return 0;
        }
        /* Reals that are integer-valued can also be converted */
        if (value->type == REAL_VALUE && type == INT_VALUE)
        {
            int i1 = static_cast<int>(value->u.r.r1);
            int i2 = static_cast<int>(value->u.r.r2);
            if (gmx_within_tol(value->u.r.r1, i1, GMX_REAL_EPS)
                && gmx_within_tol(value->u.r.r2, i2, GMX_REAL_EPS))
            {
                *value = SelectionParserValue::createIntegerRange(i1, i2);
                return 0;
            }
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
convert_values(SelectionParserValueList *values, e_selvalue_t type, void *scanner)
{
    int rc = 0;
    SelectionParserValueList::iterator value;
    for (value = values->begin(); value != values->end(); ++value)
    {
        int rc1 = convert_value(&*value, type, scanner);
        if (rc1 != 0 && rc == 0)
        {
            rc = rc1;
        }
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
            gmx_ana_selparam_t                *param)
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
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_range(const SelectionParserValueList &values,
                   gmx_ana_selparam_t *param, void *scanner)
{
    int                *idata;
    real               *rdata;
    int                 i, j, n;

    param->flags &= ~SPAR_DYNAMIC;
    GMX_RELEASE_ASSERT(param->val.type == INT_VALUE || param->val.type == REAL_VALUE,
                       "Invalid range parameter type");
    idata = NULL;
    rdata = NULL;
    if (param->val.type == INT_VALUE)
    {
        snew(idata, values.size()*2);
    }
    else
    {
        snew(rdata, values.size()*2);
    }
    i = 0;
    SelectionParserValueList::const_iterator value;
    for (value = values.begin(); value != values.end(); ++value)
    {
        if (value->hasExpressionValue())
        {
            _gmx_selparser_error(scanner, "expressions not supported within range parameters");
            return false;
        }
        GMX_RELEASE_ASSERT(value->type == param->val.type,
                           "Invalid range value type (should have been caught earlier)");
        if (param->val.type == INT_VALUE)
        {
            int i1 = std::min(value->u.i.i1, value->u.i.i2);
            int i2 = std::max(value->u.i.i1, value->u.i.i2);
            /* Check if the new range overlaps or extends the previous one */
            if (i > 0 && i1 <= idata[i-1]+1 && i2 >= idata[i-2]-1)
            {
                idata[i-2] = std::min(idata[i-2], i1);
                idata[i-1] = std::max(idata[i-1], i2);
            }
            else
            {
                idata[i++] = i1;
                idata[i++] = i2;
            }
        }
        else
        {
            real r1 = std::min(value->u.r.r1, value->u.r.r2);
            real r2 = std::max(value->u.r.r1, value->u.r.r2);
            /* Check if the new range overlaps or extends the previous one */
            if (i > 0 && r1 <= rdata[i-1] && r2 >= rdata[i-2])
            {
                rdata[i-2] = std::min(rdata[i-2], r1);
                rdata[i-1] = std::max(rdata[i-1], r2);
            }
            else
            {
                rdata[i++] = r1;
                rdata[i++] = r2;
            }
        }
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
                j         += 2;
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
                j         += 2;
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
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_varnum(const SelectionParserValueList    &values,
                    gmx_ana_selparam_t                *param,
                    const SelectionTreeElementPointer &root,
                    void                              *scanner)
{
    int                 i, j;

    param->flags &= ~SPAR_DYNAMIC;
    /* Compute number of values, considering also integer ranges. */
    size_t valueCount = values.size();
    if (param->val.type == INT_VALUE)
    {
        SelectionParserValueList::const_iterator value;
        for (value = values.begin(); value != values.end(); ++value)
        {
            if (value->type == INT_VALUE && !value->hasExpressionValue())
            {
                valueCount += abs(value->u.i.i2 - value->u.i.i1);
            }
        }
    }

    /* Check that the value type is actually implemented */
    if (param->val.type != INT_VALUE && param->val.type != REAL_VALUE
        && param->val.type != STR_VALUE && param->val.type != POS_VALUE)
    {
        GMX_THROW(gmx::InternalError("Variable-count value type not implemented"));
    }

    /* Reserve appropriate amount of memory */
    if (param->val.type == POS_VALUE)
    {
        gmx_ana_pos_reserve(param->val.u.p, valueCount, 0);
        gmx_ana_indexmap_init(&param->val.u.p->m, NULL, NULL, INDEX_UNKNOWN);
        gmx_ana_pos_set_nr(param->val.u.p, valueCount);
    }
    else
    {
        _gmx_selvalue_reserve(&param->val, valueCount);
    }

    i     = 0;
    SelectionParserValueList::const_iterator value;
    for (value = values.begin(); value != values.end(); ++value)
    {
        if (value->hasExpressionValue())
        {
            _gmx_selparser_error(scanner, "expressions not supported within value lists");
            return false;
        }
        GMX_RELEASE_ASSERT(value->type == param->val.type,
                           "Invalid value type (should have been caught earlier)");
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
            case STR_VALUE:
                param->val.u.s[i++] = strdup(value->stringValue().c_str());
                break;
            case POS_VALUE:  copy_rvec(value->u.x, param->val.u.p->x[i++]); break;
            default: /* Should not be reached */
                GMX_THROW(gmx::InternalError("Variable-count value type not implemented"));
        }
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
        child->v.nr   = param->val.nr;
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
    child->flags  &= ~SEL_ALLOCVAL;
    child->u.param = param;
    if (child->v.type != param->val.type)
    {
        _gmx_selparser_error(scanner, "invalid expression value");
        // FIXME: Use exceptions.
        return SelectionTreeElementPointer();
    }
    _gmx_selelem_update_flags(child);
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
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_varnum_expr(const SelectionParserValueList    &values,
                         gmx_ana_selparam_t                *param,
                         const SelectionTreeElementPointer &root,
                         void                              *scanner)
{
    GMX_RELEASE_ASSERT(values.size() == 1 && values.front().hasExpressionValue(),
                       "Called with an invalid type of value");

    SelectionTreeElementPointer child
        = add_child(root, param, values.front().expr, scanner);
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
            GMX_THROW(gmx::InternalError("Invalid value type"));
    }
    sel->v.nr     = 1;
    sel->v.nalloc = -1;
    return true;
}

/*! \brief
 * Parses the values for a parameter that takes a constant number of values.
 *
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static bool
parse_values_std(const SelectionParserValueList &values,
                 gmx_ana_selparam_t *param,
                 const SelectionTreeElementPointer &root, void *scanner)
{
    int                i, j;
    bool               bDynamic;

    /* Handle atom-valued parameters */
    if (param->flags & SPAR_ATOMVAL)
    {
        if (values.size() > 1)
        {
            _gmx_selparser_error(scanner, "more than one value not supported");
            return false;
        }
        if (values.front().hasExpressionValue())
        {
            SelectionTreeElementPointer child
                = add_child(root, param, values.front().expr, scanner);
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

    i        = 0;
    bDynamic = false;
    SelectionParserValueList::const_iterator value;
    for (value = values.begin(); value != values.end() && i < param->val.nr; ++value)
    {
        if (value->type != param->val.type)
        {
            _gmx_selparser_error(scanner, "incorrect value skipped");
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
                    param->val.u.s[i] = strdup(value->stringValue().c_str());
                    break;
                case POS_VALUE:
                    gmx_ana_pos_init_const(&param->val.u.p[i], value->u.x);
                    break;
                case NO_VALUE:
                case GROUP_VALUE:
                    GMX_THROW(gmx::InternalError("Invalid non-expression value type"));
            }
        }
        ++i;
    }
    if (value != values.end())
    {
        _gmx_selparser_error(scanner, "extra values");
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
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_bool(const std::string &name,
                  const SelectionParserValueList &values,
                  gmx_ana_selparam_t *param, void *scanner)
{
    GMX_ASSERT(param->val.type == NO_VALUE,
               "Boolean parser called for non-boolean parameter");
    if (values.size() > 1 || (!values.empty() && values.front().type != INT_VALUE))
    {
        _gmx_selparser_error(scanner, "parameter takes only a yes/no/on/off/0/1 value");
        return false;
    }

    bool bSetNo = false;
    /* Check if the parameter name is given with a 'no' prefix */
    if (name.length() > 2 && name[0] == 'n' && name[1] == 'o'
        && name.compare(2, name.length() - 2, param->name) == 0)
    {
        bSetNo = true;
    }
    if (bSetNo && !values.empty())
    {
        _gmx_selparser_error(scanner, "parameter 'no%s' should not have a value",
                             param->name);
        return false;
    }
    if (!values.empty() && values.front().u.i.i1 == 0)
    {
        bSetNo = true;
    }

    *param->val.u.b = bSetNo ? false : true;
    return true;
}

/*! \brief
 * Parses the values for an enumeration parameter.
 *
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static bool
parse_values_enum(const SelectionParserValueList &values,
                  gmx_ana_selparam_t             *param,
                  void                           *scanner)
{
    GMX_ASSERT(param->val.type == STR_VALUE,
               "Enum parser called for non-string parameter");
    if (values.size() != 1)
    {
        _gmx_selparser_error(scanner, "a single value is required");
        return false;
    }
    const SelectionParserValue &value = values.front();
    GMX_RELEASE_ASSERT(value.type == param->val.type,
                       "Invalid value type (should have been caught earlier)");
    if (value.hasExpressionValue())
    {
        _gmx_selparser_error(scanner, "expression value for enumerated parameter not supported");
        return false;
    }

    const std::string &svalue = value.stringValue();
    int                i      = 1;
    int                match  = 0;
    while (param->val.u.s[i] != NULL)
    {
        if (gmx::startsWith(param->val.u.s[i], svalue))
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
convert_const_values(SelectionParserValueList *values)
{
    SelectionParserValueList::iterator value;
    for (value = values->begin(); value != values->end(); ++value)
    {
        if (value->hasExpressionValue() && value->expr->v.type != GROUP_VALUE &&
            value->expr->type == SEL_CONST)
        {
            SelectionTreeElementPointer expr = value->expr;
            switch (expr->v.type)
            {
                case INT_VALUE:
                    *value = SelectionParserValue::createInteger(expr->v.u.i[0]);
                    break;
                case REAL_VALUE:
                    *value = SelectionParserValue::createReal(expr->v.u.r[0]);
                    break;
                case STR_VALUE:
                    *value = SelectionParserValue::createString(expr->v.u.s[0]);
                    break;
                case POS_VALUE:
                    *value = SelectionParserValue::createPosition(expr->v.u.p->x[0]);
                    break;
                default:
                    GMX_THROW(gmx::InternalError(
                                      "Unsupported constant expression value type"));
            }
        }
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
_gmx_sel_parse_params(const gmx::SelectionParserParameterList &pparams,
                      int nparam, gmx_ana_selparam_t *params,
                      const gmx::SelectionTreeElementPointer &root,
                      void *scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    gmx_ana_selparam_t          *oparam;
    bool                         bOk, rc;
    int                          i;

    /* Check that the value pointers of SPAR_VARNUM parameters are NULL and
     * that they are not NULL for other parameters */
    bOk = true;
    for (i = 0; i < nparam; ++i)
    {
        std::string                contextStr = gmx::formatString("In parameter '%s'", params[i].name);
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
        return false;
    }
    /* Parse the parameters */
    i = 0;
    SelectionParserParameterList::const_iterator pparam;
    for (pparam = pparams.begin(); pparam != pparams.end(); ++pparam)
    {
        std::string contextStr;
        /* Find the parameter and make some checks */
        if (!pparam->name().empty())
        {
            contextStr = gmx::formatString("In parameter '%s'", pparam->name().c_str());
            i          = -1;
            oparam     = gmx_ana_selparam_find(pparam->name().c_str(), nparam, params);
        }
        else if (i >= 0)
        {
            contextStr = gmx::formatString("In value %d", i + 1);
            oparam     = &params[i];
            if (oparam->name != NULL)
            {
                oparam = NULL;
                _gmx_selparser_error(scanner, "too many NULL parameters provided");
                bOk = false;
                continue;
            }
            ++i;
        }
        else
        {
            _gmx_selparser_error(scanner, "all NULL parameters should appear in the beginning of the list");
            bOk = false;
            continue;
        }
        gmx::MessageStringContext  context(errors, contextStr);
        if (!oparam)
        {
            _gmx_selparser_error(scanner, "unknown parameter skipped");
            bOk = false;
            continue;
        }
        if (oparam->val.type != NO_VALUE && pparam->values().empty())
        {
            _gmx_selparser_error(scanner, "no value provided");
            bOk = false;
            continue;
        }
        if (oparam->flags & SPAR_SET)
        {
            _gmx_selparser_error(scanner, "parameter set multiple times, extra values skipped");
            bOk = false;
            continue;
        }
        oparam->flags |= SPAR_SET;
        /* Process the values for the parameter */
        convert_const_values(pparam->values_.get());
        if (convert_values(pparam->values_.get(), oparam->val.type, scanner) != 0)
        {
            _gmx_selparser_error(scanner, "invalid value");
            bOk = false;
            continue;
        }
        if (oparam->val.type == NO_VALUE)
        {
            rc = parse_values_bool(pparam->name(), pparam->values(), oparam, scanner);
        }
        else if (oparam->flags & SPAR_RANGES)
        {
            rc = parse_values_range(pparam->values(), oparam, scanner);
        }
        else if (oparam->flags & SPAR_VARNUM)
        {
            if (pparam->values().size() == 1
                && pparam->values().front().hasExpressionValue())
            {
                rc = parse_values_varnum_expr(pparam->values(), oparam, root, scanner);
            }
            else
            {
                rc = parse_values_varnum(pparam->values(), oparam, root, scanner);
            }
        }
        else if (oparam->flags & SPAR_ENUMVAL)
        {
            rc = parse_values_enum(pparam->values(), oparam, scanner);
        }
        else
        {
            rc = parse_values_std(pparam->values(), oparam, root, scanner);
        }
        if (!rc)
        {
            bOk = false;
        }
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

    return bOk;
}
