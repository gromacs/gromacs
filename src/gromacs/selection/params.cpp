/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include <algorithm>
#include <string>

#include "gromacs/math/vec.h"
#include "gromacs/selection/position.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/scoped_cptr.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "parsetree.h"
#include "scanner.h"
#include "selelem.h"
#include "selmethod.h"
#include "selparam.h"

using namespace gmx;

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
 * \param[in]     errors   Errors will be reported into this as nested exceptions.
 * \param[in]     scanner  Scanner data structure.
 */
static void
convert_value(SelectionParserValue *value, e_selvalue_t type,
              ExceptionInitializer *errors, void *scanner)
{
    if (value->type == type || type == NO_VALUE)
    {
        return;
    }
    if (value->hasExpressionValue())
    {
        /* Conversion from atom selection to position using default
         * reference positions. */
        if (value->type == GROUP_VALUE && type == POS_VALUE)
        {
            try
            {
                SelectionTreeElementPointer expr =
                    _gmx_sel_init_position(value->expr, NULL, scanner);
                *value = SelectionParserValue::createExpr(expr);
            }
            catch (UserInputError &ex)
            {
                std::string text(_gmx_sel_lexer_get_text(scanner, value->location()));
                std::string context(formatString("In '%s'", text.c_str()));
                ex.prependContext(context);
                errors->addCurrentExceptionAsNested();
            }
            return;
        }
    }
    else
    {
        /* Integers to floating point are easy */
        if (value->type == INT_VALUE && type == REAL_VALUE)
        {
            *value = SelectionParserValue::createRealRange(value->u.i.i1,
                                                           value->u.i.i2,
                                                           value->location());
            return;
        }
        /* Reals that are integer-valued can also be converted */
        if (value->type == REAL_VALUE && type == INT_VALUE)
        {
            int i1 = static_cast<int>(value->u.r.r1);
            int i2 = static_cast<int>(value->u.r.r2);
            if (gmx_within_tol(value->u.r.r1, i1, GMX_REAL_EPS)
                && gmx_within_tol(value->u.r.r2, i2, GMX_REAL_EPS))
            {
                *value = SelectionParserValue::createIntegerRange(i1, i2, value->location());
                return;
            }
        }
    }
    std::string text(_gmx_sel_lexer_get_text(scanner, value->location()));
    std::string message(
            formatString("Expression '%s' evaluates to a type is not valid in this context",
                         text.c_str()));
    InvalidInputError ex(message);
    errors->addNested(ex);
}

/*! \brief
 * Does a type conversion on a list of values.
 *
 * \param[in,out] values   Values to convert.
 * \param[in]     type     Type to convert to.
 * \param[in]     scanner  Scanner data structure.
 */
static void
convert_values(SelectionParserValueList *values, e_selvalue_t type, void *scanner)
{
    ExceptionInitializer               errors("");
    SelectionParserValueList::iterator value;
    for (value = values->begin(); value != values->end(); ++value)
    {
        convert_value(&*value, type, &errors, scanner);
    }
    if (errors.hasNestedExceptions())
    {
        GMX_THROW(InvalidInputError(errors));
    }
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
 */
static void
parse_values_range(const SelectionParserValueList &values,
                   gmx_ana_selparam_t *param, void *scanner)
{
    int                 i, j, n;

    param->flags &= ~SPAR_DYNAMIC;
    GMX_RELEASE_ASSERT(param->val.type == INT_VALUE || param->val.type == REAL_VALUE,
                       "Invalid range parameter type");
    int                *idata = NULL;
    real               *rdata = NULL;
    scoped_guard_sfree  dataGuard;
    if (param->val.type == INT_VALUE)
    {
        snew(idata, values.size()*2);
        dataGuard.reset(idata);
    }
    else
    {
        snew(rdata, values.size()*2);
        dataGuard.reset(rdata);
    }
    i = 0;
    SelectionParserValueList::const_iterator value;
    for (value = values.begin(); value != values.end(); ++value)
    {
        GMX_RELEASE_ASSERT(value->type == param->val.type,
                           "Invalid range value type (should have been caught earlier)");
        if (value->hasExpressionValue())
        {
            std::string       text(_gmx_sel_lexer_get_text(scanner, value->location()));
            std::string       message("Only simple values or 'A to B' ranges are "
                                      "supported in this context");
            InvalidInputError ex(message);
            ex.prependContext(formatString("Invalid expression '%s'", text.c_str()));
            GMX_THROW(ex);
        }
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
            if (rdata[j-1] >= rdata[i])
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
        dataGuard.release();
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
            GMX_ASSERT(n == 1,
                       "Range parameters with a fixed count > 1 do not make sense");
            GMX_THROW(InvalidInputError("Only one value or 'A to B' range is "
                                        "supported in this context"));
        }
        if (param->val.type == INT_VALUE)
        {
            memcpy(param->val.u.i, idata, 2*n*sizeof(int));
        }
        else
        {
            memcpy(param->val.u.r, rdata, 2*n*sizeof(real));
        }
    }
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;
}

/*! \brief
 * Parses the values for a parameter that takes a variable number of values.
 *
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static void
parse_values_varnum(const SelectionParserValueList    &values,
                    gmx_ana_selparam_t                *param,
                    const SelectionTreeElementPointer &root,
                    void                              *scanner)
{
    int                 i, j;

    param->flags &= ~SPAR_DYNAMIC;
    /* Compute number of values, considering also integer ranges. */
    int valueCount = static_cast<int>(values.size());
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
        GMX_THROW(InternalError("Variable-count value type not implemented"));
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
    /* Create a dummy child element to store the string values.
     * This element is responsible for freeing the values, but carries no
     * other function. */
    if (param->val.type == STR_VALUE)
    {
        SelectionTreeElementPointer child(
                new SelectionTreeElement(SEL_CONST, SelectionLocation::createEmpty()));
        _gmx_selelem_set_vtype(child, STR_VALUE);
        child->setName(param->name);
        child->flags &= ~SEL_ALLOCVAL;
        child->flags |= SEL_FLAGSSET | SEL_VARNUMVAL | SEL_ALLOCDATA;
        child->v.nr   = valueCount;
        _gmx_selvalue_setstore(&child->v, param->val.u.s);
        /* Because the child is not group-valued, the u union is not used
         * for anything, so we can abuse it by storing the parameter value
         * as place_child() expects, but this is really ugly... */
        child->u.param = param;
        place_child(root, child, param);
    }
    param->val.nr = valueCount;

    i     = 0;
    SelectionParserValueList::const_iterator value;
    for (value = values.begin(); value != values.end(); ++value)
    {
        GMX_RELEASE_ASSERT(value->type == param->val.type,
                           "Invalid value type (should have been caught earlier)");
        if (value->hasExpressionValue())
        {
            std::string       text(_gmx_sel_lexer_get_text(scanner, value->location()));
            std::string       message("Selection expressions are not supported in this "
                                      "context when multiple values are provided");
            InvalidInputError ex(message);
            ex.prependContext(formatString("Invalid expression '%s'", text.c_str()));
            GMX_THROW(ex);
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
                    std::string text(_gmx_sel_lexer_get_text(scanner, value->location()));
                    std::string message
                        = formatString("Real range ('%s') is not supported in this context",
                                       text.c_str());
                    InvalidInputError ex(message);
                    GMX_THROW(ex);
                }
                param->val.u.r[i++] = value->u.r.r1;
                break;
            case STR_VALUE:
                param->val.u.s[i++] = gmx_strdup(value->stringValue().c_str());
                break;
            case POS_VALUE:  copy_rvec(value->u.x, param->val.u.p->x[i++]); break;
            default: /* Should not be reached */
                GMX_RELEASE_ASSERT(false, "Variable-count value type not implemented");
        }
    }
    GMX_RELEASE_ASSERT(i == valueCount,
                       "Inconsistent value count wrt. the actual value population");
    if (param->nvalptr)
    {
        *param->nvalptr = param->val.nr;
    }
    param->nvalptr = NULL;
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
        // TODO: Initialize such that it includes the parameter.
        child.reset(new SelectionTreeElement(SEL_SUBEXPRREF, expr->location()));
        _gmx_selelem_set_vtype(child, expr->v.type);
        child->child  = expr;
    }
    /* Setup the child element */
    child->flags  &= ~SEL_ALLOCVAL;
    child->u.param = param;
    if (child->v.type != param->val.type)
    {
        // TODO: It would be nice to say what is the expected type.
        std::string text(_gmx_sel_lexer_get_text(scanner, expr->location()));
        std::string message
            = formatString("Expression '%s' is not valid in this context "
                           "(produces the wrong type of values)",
                           text.c_str());
        GMX_THROW(InvalidInputError(message));
    }
    _gmx_selelem_update_flags(child);
    if ((child->flags & SEL_DYNAMIC) && !(param->flags & SPAR_DYNAMIC))
    {
        std::string text(_gmx_sel_lexer_get_text(scanner, expr->location()));
        std::string message
            = formatString("Expression '%s' is dynamic, which is not "
                           "valid in this context",
                           text.c_str());
        GMX_THROW(InvalidInputError(message));
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
 */
static void
parse_values_varnum_expr(const SelectionParserValueList    &values,
                         gmx_ana_selparam_t                *param,
                         const SelectionTreeElementPointer &root,
                         void                              *scanner)
{
    GMX_RELEASE_ASSERT(values.size() == 1 && values.front().hasExpressionValue(),
                       "Called with an invalid type of value");

    SelectionTreeElementPointer child
        = add_child(root, param, values.front().expr, scanner);

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
        return;
    }

    if (!(child->flags & SEL_VARNUMVAL))
    {
        std::string text(_gmx_sel_lexer_get_text(scanner, values.front().location()));
        std::string message
            = formatString("Expression '%s' is invalid in this context",
                           text.c_str());
        GMX_THROW(InvalidInputError(message));
    }

    child->flags   |= SEL_ALLOCVAL;
    param->val.nr   = -1;
    *param->nvalptr = param->val.nr;
    /* Rest of the initialization is done during compilation in
     * init_method(). */
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
static void
set_expr_value_store(const SelectionTreeElementPointer &sel,
                     gmx_ana_selparam_t *param, int i, void *scanner)
{
    if (sel->v.type != GROUP_VALUE && !(sel->flags & SEL_SINGLEVAL))
    {
        std::string text(_gmx_sel_lexer_get_text(scanner, sel->location()));
        std::string message
            = formatString("Expression '%s' is invalid in this context",
                           text.c_str());
        GMX_THROW(InvalidInputError(message));
    }
    switch (sel->v.type)
    {
        case INT_VALUE:   sel->v.u.i = &param->val.u.i[i]; break;
        case REAL_VALUE:  sel->v.u.r = &param->val.u.r[i]; break;
        case STR_VALUE:   sel->v.u.s = &param->val.u.s[i]; break;
        case POS_VALUE:   sel->v.u.p = &param->val.u.p[i]; break;
        case GROUP_VALUE: sel->v.u.g = &param->val.u.g[i]; break;
        default: /* Error */
            GMX_THROW(InternalError("Invalid value type"));
    }
    sel->v.nr     = 1;
    sel->v.nalloc = -1;
}

/*! \brief
 * Parses the values for a parameter that takes a constant number of values.
 *
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param     root   Selection element to which child expressions are added.
 * \param[in] scanner Scanner data structure.
 *
 * For integer ranges, the sequence of numbers from the first to second value
 * is stored, each as a separate value.
 */
static void
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
            GMX_THROW(InvalidInputError(
                              "Only a single value or a single expression is "
                              "supported in this context"));
        }
        if (values.front().hasExpressionValue())
        {
            SelectionTreeElementPointer child
                          = add_child(root, param, values.front().expr, scanner);
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
                return;
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
            set_expr_value_store(child, param, 0, scanner);
            return;
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
        GMX_RELEASE_ASSERT(value->type == param->val.type,
                           "Invalid value type (should have been caught earlier)");
        if (value->hasExpressionValue())
        {
            SelectionTreeElementPointer child
                = add_child(root, param, value->expr, scanner);
            set_expr_value_store(child, param, i, scanner);
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
                {
                    bool bTooManyValues;
                    if (value->u.i.i1 <= value->u.i.i2)
                    {
                        for (j = value->u.i.i1; j <= value->u.i.i2 && i < param->val.nr; ++j)
                        {
                            param->val.u.i[i++] = j;
                        }
                        bTooManyValues = (j != value->u.i.i2 + 1);
                    }
                    else
                    {
                        for (j = value->u.i.i1; j >= value->u.i.i2 && i < param->val.nr; --j)
                        {
                            param->val.u.i[i++] = j;
                        }
                        bTooManyValues = (j != value->u.i.i2 - 1);
                    }
                    if (bTooManyValues)
                    {
                        std::string text(_gmx_sel_lexer_get_text(scanner, value->location()));
                        std::string message
                            = formatString("Range ('%s') produces more values than is "
                                           "accepted in this context",
                                           text.c_str());
                        GMX_THROW(InvalidInputError(message));
                    }
                    --i;
                    break;
                }
                case REAL_VALUE:
                    if (value->u.r.r1 != value->u.r.r2)
                    {
                        std::string text(_gmx_sel_lexer_get_text(scanner, value->location()));
                        std::string message
                            = formatString("Real range ('%s') is not supported in this context",
                                           text.c_str());
                        GMX_THROW(InvalidInputError(message));
                    }
                    param->val.u.r[i] = value->u.r.r1;
                    break;
                case STR_VALUE:
                    param->val.u.s[i] = gmx_strdup(value->stringValue().c_str());
                    break;
                case POS_VALUE:
                    gmx_ana_pos_init_const(&param->val.u.p[i], value->u.x);
                    break;
                case NO_VALUE:
                case GROUP_VALUE:
                    GMX_THROW(InternalError("Invalid non-expression value type"));
            }
        }
        ++i;
    }
    if (value != values.end())
    {
        std::string message
            = formatString("Too many values provided, expected %d",
                           param->val.nr);
        GMX_THROW(InvalidInputError(message));
    }
    if (i < param->val.nr)
    {
        std::string message
            = formatString("Too few values provided, expected %d",
                           param->val.nr);
        GMX_THROW(InvalidInputError(message));
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
}

/*! \brief
 * Parses the values for a boolean parameter.
 *
 * \param[in] name   Name by which the parameter was given.
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 */
static void
parse_values_bool(const std::string &name,
                  const SelectionParserValueList &values,
                  gmx_ana_selparam_t *param, void *scanner)
{
    GMX_UNUSED_VALUE(scanner);
    GMX_ASSERT(param->val.type == NO_VALUE,
               "Boolean parser called for non-boolean parameter");
    if (values.size() > 1 || (!values.empty() && values.front().type != INT_VALUE))
    {
        std::string message
            = formatString("'%s' only accepts yes/no/on/off/0/1 (and empty) as a value",
                           param->name);
        GMX_THROW(InvalidInputError(message));
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
        std::string message
            = formatString("'no%s' cannot be followed by any value",
                           param->name);
        GMX_THROW(InvalidInputError(message));
    }
    if (!values.empty() && values.front().u.i.i1 == 0)
    {
        bSetNo = true;
    }

    *param->val.u.b = bSetNo ? false : true;
}

/*! \brief
 * Parses the values for an enumeration parameter.
 *
 * \param[in] values List of values.
 * \param     param  Parameter to parse.
 * \param[in] scanner Scanner data structure.
 * \returns   true if the values were parsed successfully, false otherwise.
 */
static void
parse_values_enum(const SelectionParserValueList &values,
                  gmx_ana_selparam_t             *param,
                  void                           *scanner)
{
    GMX_ASSERT(param->val.type == STR_VALUE,
               "Enum parser called for non-string parameter");
    if (values.size() != 1)
    {
        GMX_THROW(InvalidInputError(
                          "Only a single string value is supported in this context"));
    }
    const SelectionParserValue &value = values.front();
    GMX_RELEASE_ASSERT(value.type == param->val.type,
                       "Invalid value type (should have been caught earlier)");
    if (value.hasExpressionValue())
    {
        std::string text(_gmx_sel_lexer_get_text(scanner, value.location()));
        std::string message
            = formatString("Expression ('%s') is not supported in this context",
                           text.c_str());
        GMX_THROW(InvalidInputError(message));
    }

    const std::string &svalue = value.stringValue();
    int                i      = 1;
    int                match  = 0;
    while (param->val.u.s[i] != NULL)
    {
        if (startsWith(param->val.u.s[i], svalue))
        {
            /* Check if there is a duplicate match */
            if (match > 0)
            {
                std::string message
                    = formatString("Value '%s' is ambiguous", svalue.c_str());
                GMX_THROW(InvalidInputError(message));
            }
            match = i;
        }
        ++i;
    }
    if (match == 0)
    {
        std::string message
            = formatString("Value '%s' is not recognized", svalue.c_str());
        GMX_THROW(InvalidInputError(message));
    }
    param->val.u.s[0] = param->val.u.s[match];
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
            SelectionTreeElementPointer expr     = value->expr;
            const SelectionLocation    &location = value->location();
            switch (expr->v.type)
            {
                case INT_VALUE:
                    *value = SelectionParserValue::createInteger(expr->v.u.i[0], location);
                    break;
                case REAL_VALUE:
                    *value = SelectionParserValue::createReal(expr->v.u.r[0], location);
                    break;
                case STR_VALUE:
                    *value = SelectionParserValue::createString(expr->v.u.s[0], location);
                    break;
                case POS_VALUE:
                    *value = SelectionParserValue::createPosition(expr->v.u.p->x[0], location);
                    break;
                default:
                    GMX_RELEASE_ASSERT(false,
                                       "Unsupported constant expression value type");
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
 *
 * Initializes the \p params array based on the parameters in \p pparams.
 * See the documentation of \c gmx_ana_selparam_t for different options
 * available for parsing.
 *
 * The list \p pparams and any associated values are freed after the parameters
 * have been processed, no matter is there was an error or not.
 */
void
_gmx_sel_parse_params(const gmx::SelectionParserParameterList &pparams,
                      int nparam, gmx_ana_selparam_t *params,
                      const gmx::SelectionTreeElementPointer &root,
                      void *scanner)
{
    ExceptionInitializer errors("");
    /* Check that the value pointers of SPAR_VARNUM parameters are NULL and
     * that they are not NULL for other parameters */
    for (int i = 0; i < nparam; ++i)
    {
        if (params[i].val.type != POS_VALUE
            && (params[i].flags & (SPAR_VARNUM | SPAR_ATOMVAL)))
        {
            GMX_RELEASE_ASSERT(params[i].val.u.ptr == NULL,
                               "value pointer is not NULL "
                               "although it should be for SPAR_VARNUM "
                               "and SPAR_ATOMVAL parameters");
            GMX_RELEASE_ASSERT(!((params[i].flags & SPAR_VARNUM)
                                 && (params[i].flags & SPAR_DYNAMIC))
                               || params[i].nvalptr != NULL,
                               "nvalptr is NULL but both "
                               "SPAR_VARNUM and SPAR_DYNAMIC are specified");
        }
        else
        {
            GMX_RELEASE_ASSERT(params[i].val.u.ptr != NULL,
                               "value pointer is NULL");
        }
    }
    /* Parse the parameters */
    int nullParamIndex = 0;
    SelectionParserParameterList::const_iterator pparam;
    for (pparam = pparams.begin(); pparam != pparams.end(); ++pparam)
    {
        try
        {
            // Always assigned afterwards, but cppcheck does not see that.
            gmx_ana_selparam_t *oparam = NULL;
            /* Find the parameter and make some checks */
            if (!pparam->name().empty())
            {
                nullParamIndex = -1;
                oparam
                    = gmx_ana_selparam_find(pparam->name().c_str(), nparam, params);
                GMX_RELEASE_ASSERT(oparam != NULL, "Inconsistent selection parameter");
            }
            else if (nullParamIndex >= 0)
            {
                oparam = &params[nullParamIndex];
                if (oparam->name != NULL)
                {
                    std::string text(_gmx_sel_lexer_get_text(scanner, pparam->location()));
                    std::string message
                        = formatString("Unexpected '%s'", text.c_str());
                    GMX_THROW(InvalidInputError(message));
                }
                ++nullParamIndex;
            }
            else
            {
                GMX_RELEASE_ASSERT(false, "All NULL parameters should appear in "
                                   "the beginning of the list");
            }
            if (oparam->flags & SPAR_SET)
            {
                std::string message
                    = formatString("'%s' appears multiple times",
                                   pparam->name().c_str());
                GMX_THROW(InvalidInputError(message));
            }
            oparam->flags |= SPAR_SET;
            if (oparam->val.type != NO_VALUE && pparam->values().empty())
            {
                std::string text;
                if (pparam->name().empty())
                {
                    text = root->name();
                }
                else
                {
                    text = _gmx_sel_lexer_get_text(scanner, pparam->location());
                }
                std::string message
                    = formatString("'%s' should be followed by a value/expression",
                                   text.c_str());
                GMX_THROW(InvalidInputError(message));
            }
            /* Process the values for the parameter */
            convert_const_values(pparam->values_.get());
            convert_values(pparam->values_.get(), oparam->val.type, scanner);
            if (oparam->val.type == NO_VALUE)
            {
                parse_values_bool(pparam->name(), pparam->values(), oparam, scanner);
            }
            else if (oparam->flags & SPAR_RANGES)
            {
                parse_values_range(pparam->values(), oparam, scanner);
            }
            else if (oparam->flags & SPAR_VARNUM)
            {
                if (pparam->values().size() == 1
                    && pparam->values().front().hasExpressionValue())
                {
                    parse_values_varnum_expr(pparam->values(), oparam, root, scanner);
                }
                else
                {
                    parse_values_varnum(pparam->values(), oparam, root, scanner);
                }
            }
            else if (oparam->flags & SPAR_ENUMVAL)
            {
                parse_values_enum(pparam->values(), oparam, scanner);
            }
            else
            {
                parse_values_std(pparam->values(), oparam, root, scanner);
            }
        }
        catch (UserInputError &ex)
        {
            if (!pparam->name().empty())
            {
                std::string text(_gmx_sel_lexer_get_text(scanner, pparam->location()));
                ex.prependContext(formatString("In '%s'", text.c_str()));
            }
            errors.addCurrentExceptionAsNested();
        }
    }
    /* Check that all required parameters are present */
    for (int i = 0; i < nparam; ++i)
    {
        if (!(params[i].flags & SPAR_OPTIONAL) && !(params[i].flags & SPAR_SET))
        {
            std::string message;
            if (params[i].name == NULL)
            {
                message = formatString("'%s' should be followed by a value/expression",
                                       root->name().c_str());
            }
            else
            {
                message = formatString("'%s' is missing", params[i].name);
            }
            InvalidInputError ex(message);
            errors.addNested(ex);
        }
    }
    if (errors.hasNestedExceptions())
    {
        GMX_THROW(InvalidInputError(errors));
    }
}
