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
 * Implements the \p same selection method.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <stdlib.h>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "keywords.h"
#include "parsetree.h"
#include "selelem.h"
#include "selmethod.h"

/*! \internal
 * \brief
 * Data structure for the \p same selection method.
 *
 * To avoid duplicate initialization code, the same data structure is used
 * for matching both integer and string keywords; hence the unions.
 *
 * \ingroup module_selection
 */
typedef struct
{
    /** Value for each atom to match. */
    union
    {
        int                 *i;
        char               **s;
        void                *ptr;
    }                        val;
    /*! \brief
     * Number of values in the \p as array.
     *
     * For string values, this is actually the number of values in the
     * \p as_s_sorted array.
     */
    int                      nas;
    /** Values to match against. */
    union
    {
        int                 *i;
        char               **s;
        void                *ptr;
    }                        as;
    /*! \brief
     * Separate array for sorted \p as.s array.
     *
     * The array of strings returned as the output value of a parameter should
     * not be messed with to avoid memory corruption (the pointers in the array
     * may be reused for several evaluations), so we keep our own copy for
     * modifications.
     */
    char                   **as_s_sorted;
    /** Whether simple matching can be used. */
    bool                     bSorted;
} t_methoddata_same;

/*! \brief
 * Allocates data for the \p same selection method.
 *
 * \param[in]     npar  Not used (should be 2).
 * \param[in,out] param Method parameters (should point to a copy of
 *      ::smparams_same_int or ::smparams_same_str).
 * \returns Pointer to the allocated data (\ref t_methoddata_same).
 */
static void *
init_data_same(int npar, gmx_ana_selparam_t *param);
/*! \brief
 * Initializes the \p same selection method.
 *
 * \param   top   Not used.
 * \param   npar  Not used (should be 2).
 * \param   param Initialized method parameters (should point to a copy of
 *      ::smparams_same_int or ::smparams_same_str).
 * \param   data  Pointer to \ref t_methoddata_same to initialize.
 * \returns 0 on success, -1 on failure.
 */
static void
init_same(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the data allocated for the \p same selection method. */
static void
free_data_same(void *data);
/*! \brief
 * Initializes the evaluation of the \p same selection method for a frame.
 *
 * \param[in]  top  Not used.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \ref t_methoddata_same.
 *
 * Sorts the \c data->as.i array and removes identical values for faster and
 * simpler lookup.
 */
static void
init_frame_same_int(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/** Evaluates the \p same selection method. */
static void
evaluate_same_int(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/*! \brief
 * Initializes the evaluation of the \p same selection method for a frame.
 *
 * \param[in]  top  Not used.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \ref t_methoddata_same.
 *
 * Sorts the \c data->as.s array and removes identical values for faster and
 * simpler lookup.
 */
static void
init_frame_same_str(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/** Evaluates the \p same selection method. */
static void
evaluate_same_str(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/** Parameters for the \p same selection method. */
static gmx_ana_selparam_t smparams_same_int[] = {
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"as", {INT_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
};

/** Parameters for the \p same selection method. */
static gmx_ana_selparam_t smparams_same_str[] = {
    {NULL, {STR_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"as", {STR_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
};

/** Help text for the \p same selection method. */
static const char *const help_same[] = {
    "::",
    "",
    "  same KEYWORD as ATOM_EXPR",
    "",

    "The keyword [TT]same[tt] can be used to select all atoms for which",
    "the given [TT]KEYWORD[tt] matches any of the atoms in [TT]ATOM_EXPR[tt].",
    "Keywords that evaluate to integer or string values are supported.",
};

/*! \internal \brief Selection method data for the \p same method. */
gmx_ana_selmethod_t sm_same = {
    "same", GROUP_VALUE, 0,
    asize(smparams_same_int), smparams_same_int,
    &init_data_same,
    NULL,
    &init_same,
    NULL,
    &free_data_same,
    &init_frame_same_int,
    &evaluate_same_int,
    NULL,
    {"same KEYWORD as ATOM_EXPR",
     "Extending selections", asize(help_same), help_same},
};

/*! \brief
 * Selection method data for the \p same method.
 *
 * This selection method is used for matching string keywords. The parser
 * never sees this method; _gmx_selelem_custom_init_same() replaces sm_same
 * with this method in cases where it is required.
 */
static gmx_ana_selmethod_t sm_same_str = {
    "same", GROUP_VALUE, SMETH_SINGLEVAL,
    asize(smparams_same_str), smparams_same_str,
    &init_data_same,
    NULL,
    &init_same,
    NULL,
    &free_data_same,
    &init_frame_same_str,
    &evaluate_same_str,
    NULL,
    {NULL, NULL, 0, NULL},
};

static void *
init_data_same(int /* npar */, gmx_ana_selparam_t *param)
{
    t_methoddata_same *data;

    snew(data, 1);
    data->as_s_sorted = NULL;
    param[1].nvalptr  = &data->nas;
    return data;
}

/*!
 * \param[in,out] method  The method to initialize.
 * \param[in,out] params  Pointer to the first parameter.
 * \param[in]     scanner Scanner data structure.
 *
 * If \p *method is not a \c same method, this function returns
 * immediately.
 */
void
_gmx_selelem_custom_init_same(gmx_ana_selmethod_t                           **method,
                              const gmx::SelectionParserParameterListPointer &params,
                              void                                           *scanner)
{

    /* Do nothing if this is not a same method. */
    if (!*method || (*method)->name != sm_same.name || params->empty())
    {
        return;
    }

    const gmx::SelectionParserValueList &kwvalues = params->front().values();
    if (kwvalues.size() != 1 || !kwvalues.front().hasExpressionValue()
        || kwvalues.front().expr->type != SEL_EXPRESSION)
    {
        GMX_THROW(gmx::InvalidInputError("'same' should be followed by a single keyword"));
    }
    gmx_ana_selmethod_t *kwmethod = kwvalues.front().expr->u.expr.method;
    if (kwmethod->type == STR_VALUE)
    {
        *method = &sm_same_str;
    }

    /* We do custom processing for the "as" parameter. */
    gmx::SelectionParserParameterList::iterator asparam = ++params->begin();
    if (asparam != params->end() && asparam->name() == sm_same.param[1].name)
    {
        const gmx::SelectionParserValueList &asvalues = asparam->values();
        if (asvalues.size() != 1 || !asvalues.front().hasExpressionValue())
        {
            // TODO: Think about providing more informative context.
            GMX_THROW(gmx::InvalidInputError("'same ... as' should be followed by a single expression"));
        }
        const gmx::SelectionTreeElementPointer &child = asvalues.front().expr;
        /* Create a second keyword evaluation element for the keyword given as
         * the first parameter, evaluating the keyword in the group given by the
         * second parameter. */
        gmx::SelectionTreeElementPointer kwelem
            = _gmx_sel_init_keyword_evaluator(kwmethod, child, scanner);
        /* Replace the second parameter with one with a value from \p kwelem. */
        std::string pname = asparam->name();
        *asparam = gmx::SelectionParserParameter::createFromExpression(pname, kwelem);
    }
}

static void
init_same(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;

    d->val.ptr = param[0].val.u.ptr;
    d->as.ptr  = param[1].val.u.ptr;
    if (param[1].val.type == STR_VALUE)
    {
        snew(d->as_s_sorted, d->nas);
    }
    if (!(param[0].flags & SPAR_ATOMVAL))
    {
        GMX_THROW(gmx::InvalidInputError(
                          "The 'same' selection keyword combined with a "
                          "non-keyword does not make sense"));
    }
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_same).
 */
static void
free_data_same(void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;

    sfree(d->as_s_sorted);
    sfree(d);
}

/*! \brief
 * Helper function for comparison of two integers.
 */
static int
cmp_int(const void *a, const void *b)
{
    if (*(int *)a < *(int *)b)
    {
        return -1;
    }
    if (*(int *)a > *(int *)b)
    {
        return 1;
    }
    return 0;
}

static void
init_frame_same_int(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;
    int                i, j;

    /* Collapse adjacent values, and check whether the array is sorted. */
    d->bSorted = true;
    if (d->nas == 0)
    {
        return;
    }
    for (i = 1, j = 0; i < d->nas; ++i)
    {
        if (d->as.i[i] != d->as.i[j])
        {
            if (d->as.i[i] < d->as.i[j])
            {
                d->bSorted = false;
            }
            ++j;
            d->as.i[j] = d->as.i[i];
        }
    }
    d->nas = j + 1;

    if (!d->bSorted)
    {
        qsort(d->as.i, d->nas, sizeof(d->as.i[0]), &cmp_int);
        /* More identical values may become adjacent after sorting. */
        for (i = 1, j = 0; i < d->nas; ++i)
        {
            if (d->as.i[i] != d->as.i[j])
            {
                ++j;
                d->as.i[j] = d->as.i[i];
            }
        }
        d->nas = j + 1;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_same.
 *
 * Calculates which values in \c data->val.i can be found in \c data->as.i
 * (assumed sorted), and writes the corresponding atoms to output.
 * If \c data->val is sorted, uses a linear scan of both arrays, otherwise a
 * binary search of \c data->as is performed for each block of values in
 * \c data->val.
 */
static void
evaluate_same_int(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_same     *d = (t_methoddata_same *)data;
    int                    i, j;

    out->u.g->isize = 0;
    i               = j = 0;
    while (j < g->isize)
    {
        if (d->bSorted)
        {
            /* If we are sorted, we can do a simple linear scan. */
            while (i < d->nas && d->as.i[i] < d->val.i[j])
            {
                ++i;
            }
        }
        else
        {
            /* If not, we must do a binary search of all the values. */
            int i1, i2;

            i1 = 0;
            i2 = d->nas;
            while (i2 - i1 > 1)
            {
                int itry = (i1 + i2) / 2;
                if (d->as.i[itry] <= d->val.i[j])
                {
                    i1 = itry;
                }
                else
                {
                    i2 = itry;
                }
            }
            i = (d->as.i[i1] == d->val.i[j] ? i1 : d->nas);
        }
        /* Check whether the value was found in the as list. */
        if (i == d->nas || d->as.i[i] != d->val.i[j])
        {
            /* If not, skip all atoms with the same value. */
            int tmpval = d->val.i[j];
            ++j;
            while (j < g->isize && d->val.i[j] == tmpval)
            {
                ++j;
            }
        }
        else
        {
            /* Copy all the atoms with this value to the output. */
            while (j < g->isize && d->val.i[j] == d->as.i[i])
            {
                out->u.g->index[out->u.g->isize++] = g->index[j];
                ++j;
            }
        }
        if (j < g->isize && d->val.i[j] < d->val.i[j - 1])
        {
            d->bSorted = false;
        }
    }
}

/*! \brief
 * Helper function for comparison of two strings.
 */
static int
cmp_str(const void *a, const void *b)
{
    return strcmp(*(char **)a, *(char **)b);
}

static void
init_frame_same_str(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;
    int                i, j;

    /* Collapse adjacent values.
     * For strings, it's unlikely that the values would be sorted originally,
     * so set bSorted always to false. */
    d->bSorted        = false;
    if (d->nas == 0)
    {
        return;
    }
    d->as_s_sorted[0] = d->as.s[0];
    for (i = 1, j = 0; i < d->nas; ++i)
    {
        if (strcmp(d->as.s[i], d->as_s_sorted[j]) != 0)
        {
            ++j;
            d->as_s_sorted[j] = d->as.s[i];
        }
    }
    d->nas = j + 1;

    qsort(d->as_s_sorted, d->nas, sizeof(d->as_s_sorted[0]), &cmp_str);
    /* More identical values may become adjacent after sorting. */
    for (i = 1, j = 0; i < d->nas; ++i)
    {
        if (strcmp(d->as_s_sorted[i], d->as_s_sorted[j]) != 0)
        {
            ++j;
            d->as_s_sorted[j] = d->as_s_sorted[i];
        }
    }
    d->nas = j + 1;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_same.
 *
 * Calculates which strings in \c data->val.s can be found in \c data->as.s
 * (assumed sorted), and writes the corresponding atoms to output.
 * A binary search of \c data->as is performed for each block of values in
 * \c data->val.
 */
static void
evaluate_same_str(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_same     *d = (t_methoddata_same *)data;
    int                    j;

    out->u.g->isize = 0;
    j               = 0;
    while (j < g->isize)
    {
        /* Do a binary search of the strings. */
        void *ptr;
        ptr = bsearch(&d->val.s[j], d->as_s_sorted, d->nas,
                      sizeof(d->as_s_sorted[0]), &cmp_str);
        /* Check whether the value was found in the as list. */
        if (ptr == NULL)
        {
            /* If not, skip all atoms with the same value. */
            const char *tmpval = d->val.s[j];
            ++j;
            while (j < g->isize && strcmp(d->val.s[j], tmpval) == 0)
            {
                ++j;
            }
        }
        else
        {
            const char *tmpval = d->val.s[j];
            /* Copy all the atoms with this value to the output. */
            while (j < g->isize && strcmp(d->val.s[j], tmpval) == 0)
            {
                out->u.g->index[out->u.g->isize++] = g->index[j];
                ++j;
            }
        }
    }
}
