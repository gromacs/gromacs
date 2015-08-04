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
 * Implements internal selection method for comparison expressions.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <cmath>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/utilities.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "selmethod.h"

/** Defines the comparison operator for comparison expressions. */
typedef enum
{
    CMP_INVALID,        /**< Indicates an error */
    CMP_LESS,           /**< '<' */
    CMP_LEQ,            /**< '<=' */
    CMP_GTR,            /**< '>' */
    CMP_GEQ,            /**< '>=' */
    CMP_EQUAL,          /**< '==' */
    CMP_NEQ             /**< '!=' */
} e_comparison_t;

/** The operand has a single value. */
#define CMP_SINGLEVAL  1
/** The operand value is dynamic. */
#define CMP_DYNAMICVAL 2
/** The value is real. */
#define CMP_REALVAL    4
/** The integer array is allocated. */
#define CMP_ALLOCINT   16
/** The real array is allocated. */
#define CMP_ALLOCREAL  32

/*! \internal \brief
 * Data structure for comparison expression operand values.
 */
typedef struct
{
    /** Flags that describe the type of the operand. */
    int             flags;
    /** (Array of) integer value(s). */
    int            *i;
    /** (Array of) real value(s). */
    real           *r;
} t_compare_value;

/*! \internal \brief
 * Data structure for comparison expression evaluation.
 */
typedef struct
{
    /** Comparison operator as a string. */
    char            *cmpop;
    /** Comparison operator type. */
    e_comparison_t   cmpt;
    /** Left value. */
    t_compare_value  left;
    /** Right value. */
    t_compare_value  right;
} t_methoddata_compare;

/*! \brief
 * Allocates data for comparison expression evaluation.
 *
 * \param[in]     npar  Not used (should be 5).
 * \param[in,out] param Method parameters (should point to a copy of
 *   \ref smparams_compare).
 * \returns       Pointer to the allocated data (\c t_methoddata_compare).
 *
 * Allocates memory for a \c t_methoddata_compare structure.
 */
static void *
init_data_compare(int npar, gmx_ana_selparam_t *param);
/*! \brief
 * Initializes data for comparison expression evaluation.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 5).
 * \param[in] param Method parameters (should point to \ref smparams_compare).
 * \param[in] data  Should point to a \c t_methoddata_compare.
 * \returns   0 if the input data is valid, -1 on error.
 */
static void
init_compare(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the memory allocated for comparison expression evaluation. */
static void
free_data_compare(void *data);
/** Evaluates comparison expressions. */
static void
evaluate_compare(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                 gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/** Parameters for comparison expression evaluation. */
static gmx_ana_selparam_t smparams_compare[] = {
    {"int1",  {INT_VALUE,  -1, {NULL}}, NULL,
     SPAR_OPTIONAL | SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"real1", {REAL_VALUE, -1, {NULL}}, NULL,
     SPAR_OPTIONAL | SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"op",    {STR_VALUE,   1, {NULL}}, NULL, 0},
    {"int2",  {INT_VALUE,  -1, {NULL}}, NULL,
     SPAR_OPTIONAL | SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"real2", {REAL_VALUE, -1, {NULL}}, NULL,
     SPAR_OPTIONAL | SPAR_DYNAMIC | SPAR_ATOMVAL},
};

/** \internal Selection method data for comparison expression evaluation. */
gmx_ana_selmethod_t sm_compare = {
    "cmp", GROUP_VALUE, SMETH_SINGLEVAL,
    asize(smparams_compare), smparams_compare,
    &init_data_compare,
    NULL,
    &init_compare,
    NULL,
    &free_data_compare,
    NULL,
    &evaluate_compare,
    NULL,
    {NULL, NULL, 0, NULL},
};

/*! \brief
 * Returns a \c e_comparison_t value corresponding to an operator.
 *
 * \param[in] str  String to process.
 * \returns   The comparison type corresponding to the first one or two
 *   characters of \p str.
 *
 * \p str can contain any number of characters; only the first two
 * are used.
 * If the beginning of \p str does not match any of the recognized types,
 * \ref CMP_INVALID is returned.
 */
static e_comparison_t
comparison_type(char *str)
{
    switch (str[0])
    {
        case '<': return (str[1] == '=') ? CMP_LEQ   : CMP_LESS;
        case '>': return (str[1] == '=') ? CMP_GEQ   : CMP_GTR;
        case '=': return (str[1] == '=') ? CMP_EQUAL : CMP_INVALID;
        case '!': return (str[1] == '=') ? CMP_NEQ   : CMP_INVALID;
    }
    return CMP_INVALID;
}

/*! \brief
 * Returns a string corresponding to a \c e_comparison_t value.
 *
 * \param[in] cmpt  Comparison type to convert.
 * \returns   Pointer to a string that corresponds to \p cmpt.
 *
 * The return value points to a string constant and should not be \p free'd.
 *
 * The function returns NULL if \p cmpt is not one of the valid values.
 */
static const char *
comparison_type_str(e_comparison_t cmpt)
{
    const char *p = NULL;
    switch (cmpt)
    {
        case CMP_INVALID: p = "INVALID"; break;
        case CMP_LESS:    p = "<";       break;
        case CMP_LEQ:     p = "<=";      break;
        case CMP_GTR:     p = ">";       break;
        case CMP_GEQ:     p = ">=";      break;
        case CMP_EQUAL:   p = "==";      break;
        case CMP_NEQ:     p = "!=";      break;
            // No default clause so we intentionally get compiler errors
            // if new selection choices are added later.
    }
    return p;
}

/*!
 * \param[in] fp    File to receive the output.
 * \param[in] data  Should point to a \c t_methoddata_compare.
 */
void
_gmx_selelem_print_compare_info(FILE *fp, void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;

    fprintf(fp, " \"");
    /* Print the left value */
    if ((d->left.flags & CMP_SINGLEVAL) && !(d->left.flags & CMP_DYNAMICVAL))
    {
        if (d->left.flags & CMP_REALVAL)
        {
            fprintf(fp, "%f ", d->left.r[0]);
        }
        else
        {
            fprintf(fp, "%d ", d->left.i[0]);
        }
    }
    /* Print the operator */
    if (d->cmpt != CMP_INVALID)
    {
        fprintf(fp, "%s", comparison_type_str(d->cmpt));
    }
    else
    {
        fprintf(fp, "%s", d->cmpop);
    }
    /* Print the right value */
    if ((d->right.flags & CMP_SINGLEVAL) && !(d->right.flags & CMP_DYNAMICVAL))
    {
        if (d->right.flags & CMP_REALVAL)
        {
            fprintf(fp, " %f", d->right.r[0]);
        }
        else
        {
            fprintf(fp, " %d", d->right.i[0]);
        }
    }
    fprintf(fp, "\"");
}

static void *
init_data_compare(int /* npar */, gmx_ana_selparam_t *param)
{
    t_methoddata_compare *data;

    snew(data, 1);
    param[2].val.u.s = &data->cmpop;
    return data;
}

/*! \brief
 * Reverses a comparison operator.
 *
 * \param[in] type  Comparison operator to reverse.
 * \returns   The correct comparison operator that equals \p type when the
 *   left and right sides are interchanged.
 */
static e_comparison_t
reverse_comparison_type(e_comparison_t type)
{
    switch (type)
    {
        case CMP_LESS: return CMP_GTR;
        case CMP_LEQ:  return CMP_GEQ;
        case CMP_GTR:  return CMP_LESS;
        case CMP_GEQ:  return CMP_LEQ;
        default:       break;
    }
    return type;
}

/*! \brief
 * Initializes the value storage for comparison expression.
 *
 * \param[out] val   Value structure to initialize.
 * \param[in]  param Parameters to use for initialization.
 * \returns    The number of values provided for the value, 0 on error.
 */
static int
init_comparison_value(t_compare_value *val, gmx_ana_selparam_t param[2])
{
    int  n;

    val->flags = 0;
    if (param[0].flags & SPAR_SET)
    {
        val->flags |=  (param[0].flags & SPAR_DYNAMIC) ? CMP_DYNAMICVAL : 0;
        val->flags |= !(param[0].flags & SPAR_ATOMVAL) ? CMP_SINGLEVAL  : 0;
        n           = param[0].val.nr;
        val->i      = param[0].val.u.i;
    }
    else if (param[1].flags & SPAR_SET)
    {
        val->flags |=  (param[1].flags & SPAR_DYNAMIC) ? CMP_DYNAMICVAL : 0;
        val->flags |= !(param[1].flags & SPAR_ATOMVAL) ? CMP_SINGLEVAL  : 0;
        val->flags |= CMP_REALVAL;
        n           = param[1].val.nr;
        val->r      = param[1].val.u.r;
    }
    else
    {
        n           = 0;
        val->i      = NULL;
        val->r      = NULL;
    }
    return n;
}

/*! \brief
 * Converts an integer value to floating point.
 *
 * \param[in]     n   Number of values in the \p val->u array.
 * \param[in,out] val Value to convert.
 */
static void
convert_int_real(int n, t_compare_value *val)
{
    int   i;
    real *rv;

    snew(rv, n);
    for (i = 0; i < n; ++i)
    {
        rv[i] = (real)val->i[i];
    }
    /* Free the previous value if one is present. */
    sfree(val->r);
    val->r      = rv;
    val->flags |= CMP_REALVAL | CMP_ALLOCREAL;
}

/*! \brief
 * Converts a floating point value to integer.
 *
 * \param[in]     n      Number of values in the \p val->u array.
 * \param[in,out] val    Value to convert.
 * \param[in]     cmpt   Comparison operator type.
 * \param[in]     bRight true if \p val appears on the right hand size of
 *   \p cmpt.
 * \returns       0 on success, EINVAL on error.
 *
 * The values are rounded such that the same comparison operator can be used.
 */
static void
convert_real_int(int n, t_compare_value *val, e_comparison_t cmpt, bool bRight)
{
    int   i;
    int  *iv;

    if (!bRight)
    {
        cmpt = reverse_comparison_type(cmpt);
    }
    snew(iv, n);
    /* Round according to the comparison type */
    for (i = 0; i < n; ++i)
    {
        switch (cmpt)
        {
            case CMP_LESS:
            case CMP_GEQ:
                iv[i] = static_cast<int>(std::ceil(val->r[i]));
                break;
            case CMP_GTR:
            case CMP_LEQ:
                iv[i] = static_cast<int>(std::floor(val->r[i]));
                break;
            case CMP_EQUAL:
            case CMP_NEQ:
                sfree(iv);
                /* TODO: Implement, although it isn't typically very useful.
                 * Implementation is only a matter of proper initialization,
                 * the evaluation function can already handle this case with
                 * proper preparations. */
                GMX_THROW(gmx::NotImplementedError("Equality comparison between dynamic integer and static real expressions not implemented"));
            case CMP_INVALID: /* Should not be reached */
                sfree(iv);
                GMX_THROW(gmx::InternalError("Invalid comparison type"));
        }
    }
    /* Free the previous value if one is present. */
    sfree(val->i);
    val->i      = iv;
    val->flags &= ~CMP_REALVAL;
    val->flags |= CMP_ALLOCINT;
}

static void
init_compare(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;
    int                   n1, n2;

    /* Store the values */
    n1 = init_comparison_value(&d->left, &param[0]);
    n2 = init_comparison_value(&d->right, &param[3]);
    /* Store the comparison type */
    d->cmpt = comparison_type(d->cmpop);
    if (d->cmpt == CMP_INVALID)
    {
        GMX_THROW(gmx::InternalError("Invalid comparison type"));
    }
    /* Convert the values to the same type */
    /* TODO: Currently, there are no dynamic integer-valued selection methods,
     * which means that only the branches with convert_int_real() will ever be
     * taken. It should be considered whether it is necessary to support these
     * other cases at all.
     */
    if ((d->left.flags & CMP_REALVAL) && !(d->right.flags & CMP_REALVAL))
    {
        if (d->left.flags & d->right.flags & CMP_DYNAMICVAL)
        {
            /* Nothing can be done */
        }
        else if (!(d->right.flags & CMP_DYNAMICVAL))
        {
            convert_int_real(n2, &d->right);
        }
        else /* d->left is static */
        {
            convert_real_int(n1, &d->left, d->cmpt, false);
        }
    }
    else if (!(d->left.flags & CMP_REALVAL) && (d->right.flags & CMP_REALVAL))
    {
        if (d->left.flags & d->right.flags & CMP_DYNAMICVAL)
        {
            /* Reverse the sides to place the integer on the right */
            int    flags;
            d->left.r      = d->right.r;
            d->right.r     = NULL;
            d->right.i     = d->left.i;
            d->left.i      = NULL;
            flags          = d->left.flags;
            d->left.flags  = d->right.flags;
            d->right.flags = flags;
            d->cmpt        = reverse_comparison_type(d->cmpt);
        }
        else if (!(d->left.flags & CMP_DYNAMICVAL))
        {
            convert_int_real(n1, &d->left);
        }
        else /* d->right is static */
        {
            convert_real_int(n2, &d->right, d->cmpt, true);
        }
    }
}

/*!
 * \param data Data to free (should point to a \c t_methoddata_compare).
 *
 * Frees the memory allocated for \c t_methoddata_compare.
 */
static void
free_data_compare(void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;

    sfree(d->cmpop);
    if (d->left.flags & CMP_ALLOCINT)
    {
        sfree(d->left.i);
    }
    if (d->left.flags & CMP_ALLOCREAL)
    {
        sfree(d->left.r);
    }
    if (d->right.flags & CMP_ALLOCINT)
    {
        sfree(d->right.i);
    }
    if (d->right.flags & CMP_ALLOCREAL)
    {
        sfree(d->right.r);
    }
    sfree(d);
}

/*! \brief
 * Implementation for evaluate_compare() for integer values.
 *
 * \param[in]  top   Not used.
 * \param[in]  fr    Not used.
 * \param[in]  pbc   Not used.
 * \param[in]  g     Evaluation index group.
 * \param[out] out   Output data structure (\p out->u.g is used).
 * \param[in]  data  Should point to a \c t_methoddata_compare.
 */
static void
evaluate_compare_int(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;
    int                   i, i1, i2, ig;
    int                   a, b;
    bool                  bAccept;

    GMX_UNUSED_VALUE(top);
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(pbc);
    for (i = i1 = i2 = ig = 0; i < g->isize; ++i)
    {
        a       = d->left.i[i1];
        b       = d->right.i[i2];
        bAccept = false;
        switch (d->cmpt)
        {
            case CMP_INVALID: break;
            case CMP_LESS:    bAccept = a <  b; break;
            case CMP_LEQ:     bAccept = a <= b; break;
            case CMP_GTR:     bAccept = a >  b; break;
            case CMP_GEQ:     bAccept = a >= b; break;
            case CMP_EQUAL:   bAccept = a == b; break;
            case CMP_NEQ:     bAccept = a != b; break;
        }
        if (bAccept)
        {
            out->u.g->index[ig++] = g->index[i];
        }
        if (!(d->left.flags & CMP_SINGLEVAL))
        {
            ++i1;
        }
        if (!(d->right.flags & CMP_SINGLEVAL))
        {
            ++i2;
        }
    }
    out->u.g->isize = ig;
}

/*! \brief
 * Implementation for evaluate_compare() if either value is non-integer.
 *
 * \param[in]  top   Not used.
 * \param[in]  fr    Not used.
 * \param[in]  pbc   Not used.
 * \param[in]  g     Evaluation index group.
 * \param[out] out   Output data structure (\p out->u.g is used).
 * \param[in]  data  Should point to a \c t_methoddata_compare.
 *
 * Left value is assumed to be real-valued; right value can be either.
 * This is ensured by the initialization method.
 */
static void
evaluate_compare_real(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                      gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;
    int                   i, i1, i2, ig;
    real                  a, b;
    bool                  bAccept;

    GMX_UNUSED_VALUE(top);
    GMX_UNUSED_VALUE(fr);
    GMX_UNUSED_VALUE(pbc);
    for (i = i1 = i2 = ig = 0; i < g->isize; ++i)
    {
        a       = d->left.r[i1];
        b       = (d->right.flags & CMP_REALVAL) ? d->right.r[i2] : d->right.i[i2];
        bAccept = false;
        switch (d->cmpt)
        {
            case CMP_INVALID: break;
            case CMP_LESS:    bAccept = a <  b; break;
            case CMP_LEQ:     bAccept = a <= b; break;
            case CMP_GTR:     bAccept = a >  b; break;
            case CMP_GEQ:     bAccept = a >= b; break;
            case CMP_EQUAL:   bAccept =  gmx_within_tol(a, b, GMX_REAL_EPS); break;
            case CMP_NEQ:     bAccept = !gmx_within_tol(a, b, GMX_REAL_EPS); break;
        }
        if (bAccept)
        {
            out->u.g->index[ig++] = g->index[i];
        }
        if (!(d->left.flags & CMP_SINGLEVAL))
        {
            ++i1;
        }
        if (!(d->right.flags & CMP_SINGLEVAL))
        {
            ++i2;
        }
    }
    out->u.g->isize = ig;
}

/*!
 * \param[in]  top   Not used.
 * \param[in]  fr    Not used.
 * \param[in]  pbc   Not used.
 * \param[in]  g     Evaluation index group.
 * \param[out] out   Output data structure (\p out->u.g is used).
 * \param[in]  data  Should point to a \c t_methoddata_compare.
 */
static void
evaluate_compare(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                 gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_compare *d = (t_methoddata_compare *)data;

    if (!((d->left.flags | d->right.flags) & CMP_REALVAL))
    {
        evaluate_compare_int(top, fr, pbc, g, out, data);
    }
    else
    {
        evaluate_compare_real(top, fr, pbc, g, out, data);
    }
}
