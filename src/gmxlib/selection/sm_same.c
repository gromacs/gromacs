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
 * \brief Implementation of the \p same selection method.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

#include <macros.h>
#include <smalloc.h>
#include <string2.h>

#include <selmethod.h>

#include "keywords.h"
#include "parsetree.h"
#include "selelem.h"

/*! \internal \brief
 * Data structure for the \p same selection method.
 */
typedef struct
{
    /** Value for each atom to match. */
    int                     *val;
    /** Number of values in the \p as array. */
    int                      nas;
    /** Values to match against. */
    int                     *as;
    /** Whether simple matching can be used. */
    bool                     bSorted;
} t_methoddata_same;

/** Allocates data for the \p same selection method. */
static void *
init_data_same(int npar, gmx_ana_selparam_t *param);
/** Initializes the \p same selection method. */
static int
init_same(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the data allocated for the \p same selection method. */
static void
free_data_same(void *data);
/** Initializes the evaluation of the \p same selection method for a frame. */
static int
init_frame_same(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/** Evaluates the \p same selection method. */
static int
evaluate_same(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/** Parameters for the \p same selection method. */
static gmx_ana_selparam_t smparams_same[] = {
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_ATOMVAL},
    {"as", {INT_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
};

/** Help text for the \p same selection method. */
static const char *help_same[] = {
    "EXTENDING SELECTIONS[PAR]",

    "[TT]same KEYWORD as ATOM_EXPR[tt][PAR]",

    "The keyword [TT]same[tt] can be used to select all atoms for which",
    "the given [TT]KEYWORD[tt] matches any of the atoms in [TT]ATOM_EXPR[tt].",
    "Currently, only keywords that evaluate to integer values are supported.",
};

/*! \internal \brief Selection method data for the \p same method. */
gmx_ana_selmethod_t sm_same = {
    "same", GROUP_VALUE, 0,
    asize(smparams_same), smparams_same,
    &init_data_same,
    NULL,
    &init_same,
    NULL,
    &free_data_same,
    &init_frame_same,
    &evaluate_same,
    NULL,
    {"same KEYWORD as ATOM_EXPR", asize(help_same), help_same},
};

/*!
 * \param[in]     npar  Not used (should be 2).
 * \param[in,out] param Method parameters (should point to 
 *   \ref smparams_same).
 * \returns Pointer to the allocated data (\ref t_methoddata_same).
 */
static void *
init_data_same(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_same *data;

    snew(data, 1);
    param[1].nvalptr = &data->nas;
    return data;
}

/*!
 * \param[in]     method  Selection method to query.
 * \returns   TRUE if \p method is a \c same method.
 */
bool
_gmx_selelem_is_method_same(gmx_ana_selmethod_t *method)
{
    return (method && method->name == sm_same.name);
}

/*!
 * \param[in,out] params  Pointer to the first parameter.
 * \param[in]     scanner Scanner data structure.
 * \returns       0 on success, a non-zero error code on error.
 */
int
_gmx_selelem_custom_init_same(t_selexpr_param *params,
                              void *scanner)
{
    gmx_ana_selmethod_t *method;
    t_selelem           *kwelem;
    t_selexpr_param     *param;
    char                *pname;
    int                  rc;

    if (params->nval != 1 || !params->value->bExpr
        || params->value->u.expr->type != SEL_EXPRESSION)
    {
        _gmx_selparser_error("error: 'same' should be followed by a single keyword");
        return -1;
    }
    method = params->value->u.expr->u.expr.method;

    /* We do custom processing with the second parameter, so remove it from
     * the params list, but save the name for later. */
    param        = params->next;
    params->next = NULL;
    pname        = param->name;
    param->name  = NULL;
    /* Create a second keyword evaluation element for the keyword given as
     * the first parameter, evaluating the keyword in the group given by the
     * second parameter. */
    rc = _gmx_sel_init_keyword_evaluator(&kwelem, method, param, scanner);
    if (rc != 0)
    {
        sfree(pname);
        return rc;
    }
    /* Replace the second parameter with one with a value from \p kwelem. */
    param        = _gmx_selexpr_create_param(pname);
    param->nval  = 1;
    param->value = _gmx_selexpr_create_value_expr(kwelem);
    params->next = param;
    return 0;
}

/*!
 * \param   top   Not used.
 * \param   npar  Not used (should be 2).
 * \param   param Initialized method parameters (should point to a copy of
 *      \ref smparams_same).
 * \param   data  Pointer to \ref t_methoddata_same to initialize.
 * \returns 0 on success, -1 on failure.
 */
static int
init_same(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;

    d->val = param[0].val.u.i;
    d->as  = param[1].val.u.i;
    if (!(param[0].flags & SPAR_ATOMVAL))
    {
        fprintf(stderr, "ERROR: the same selection keyword combined with a "
                        "non-keyword does not make sense\n");
        return -1;
    }
    return 0;
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_same).
 */
static void
free_data_same(void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;

    sfree(d->val);
    sfree(d->as);
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

/*!
 * \param[in]  top  Not used.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \ref t_methoddata_same.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Sorts the \c data->as array and removes identical values for faster and
 * simpler lookup.
 */
static int
init_frame_same(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;
    int                i, j;

    /* Collapse adjacent values, and check whether the array is sorted. */
    d->bSorted = TRUE;
    for (i = 1, j = 0; i < d->nas; ++i)
    {
        if (d->as[i] != d->as[j])
        {
            if (d->as[i] < d->as[j])
            {
                d->bSorted = FALSE;
            }
            ++j;
            d->as[j] = d->as[i];
        }
    }
    d->nas = j + 1;

    if (!d->bSorted)
    {
        qsort(d->as, d->nas, sizeof(d->as[0]), &cmp_int);
        /* More identical values may become adjacent after sorting. */
        for (i = 1, j = 0; i < d->nas; ++i)
        {
            if (d->as[i] != d->as[j])
            {
                ++j;
                d->as[j] = d->as[i];
            }
        }
        d->nas = j + 1;
    }
    return 0;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_same.
 *
 * Calculates which values in \c data->val can be found in \c data->as
 * (assumed sorted), and writes the corresponding atoms to output.
 * If \c data->val is sorted, uses a linear scan of both arrays, otherwise a
 * binary search of \c data->as is performed for each block of values in
 * \c data->val.
 */
static int
evaluate_same(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_same *d = (t_methoddata_same *)data;
    int                    i, j;

    out->u.g->isize = 0;
    i = j = 0;
    while (j < g->isize)
    {
        if (d->bSorted)
        {
            /* If we are sorted, we can do a simple linear scan. */
            while (i < d->nas && d->as[i] < d->val[j]) ++i;
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
                if (d->as[itry] <= d->val[j])
                {
                    i1 = itry;
                }
                else
                {
                    i2 = itry;
                }
            }
            i = (d->as[i1] == d->val[j] ? i1 : d->nas);
        }
        /* Check whether the value was found in the as list. */
        if (i == d->nas || d->as[i] != d->val[j])
        {
            /* If not, skip all atoms with the same value. */
            int tmpval = d->val[j];
            ++j;
            while (j < g->isize && d->val[j] == tmpval) ++j;
        }
        else
        {
            /* Copy all the atoms with this value to the output. */
            while (j < g->isize && d->val[j] == d->as[i])
            {
                out->u.g->index[out->u.g->isize++] = g->index[j];
                ++j;
            }
        }
        if (j < g->isize && d->val[j] < d->val[j - 1])
        {
            d->bSorted = FALSE;
        }
    }
    return 0;
}

