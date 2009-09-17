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
 * \brief Implementations of internal selection methods for integer and
 * string keyword evaluation.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <fnmatch.h>
#include <regex.h>
#include <string.h>

#include <macros.h>
#include <smalloc.h>

#include <selmethod.h>

/** Allocates data for integer keyword evaluation. */
static void *
init_data_kwint(int npar, gmx_ana_selparam_t *param);
/** Allocates data for string keyword evaluation. */
static void *
init_data_kwstr(int npar, gmx_ana_selparam_t *param);
/** Initializes data for integer keyword evaluation. */
static int
init_kwint(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Initializes data for string keyword evaluation. */
static int
init_kwstr(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the memory allocated for integer keyword evaluation. */
static void
free_data_kwint(void *data);
/** Frees the memory allocated for string keyword evaluation. */
static void
free_data_kwstr(void *data);
/** Evaluates integer selection keywords. */
static int
evaluate_keyword_int(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates string selection keywords. */
static int
evaluate_keyword_str(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/*! \internal \brief
 * Data structure for integer keyword expression evaluation.
 */
typedef struct t_methoddata_kwint
{
    /** Array of values for the keyword. */
    int               *v;
    /** Number of ranges in the \p r array. */
    int                n;
    /*! \brief
     * Array of sorted integer ranges to match against.
     *
     * Each range is made of two integers, giving the endpoints (inclusive).
     * This field stores the pointer to the ranges allocated by the
     * parameter parser; see \ref SPAR_RANGES for more information.
     */
    int               *r;
} t_methoddata_kwint;

/*! \internal \brief
 * Data structure for string keyword expression evaluation.
 */
typedef struct t_methoddata_kwstr
{
    /** Array of values for the keyword. */
    char             **v;
    /** Number of elements in the \p val array. */
    int                n;
    /** Array of strings/regular expressions to match against. */
    struct t_methoddata_kwstr_match {
        /** TRUE if the expression is a regular expression, FALSE otherwise. */
        bool           bRegExp;
        /** The value to match against. */
        union {
            /** Compiled regular expression if \p bRegExp is TRUE. */
            regex_t    r;
            /** The string if \p bRegExp is FALSE; */
            char      *s;
        }              u;
    }                 *m;
} t_methoddata_kwstr;

/** Parameters for integer keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_int[] = {
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_ATOMVAL},
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_RANGES | SPAR_VARNUM},
};

/** Parameters for string keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_str[] = {
    {NULL, {STR_VALUE, -1, {NULL}}, NULL, SPAR_ATOMVAL},
    {NULL, {STR_VALUE, -1, {NULL}}, NULL, SPAR_VARNUM},
};

/** \internal Selection method data for integer keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_int = {
    "kw_int", GROUP_VALUE, SMETH_SINGLEVAL,
    asize(smparams_keyword_int), smparams_keyword_int,
    &init_data_kwint,
     NULL,
    &init_kwint,
     NULL,
    &free_data_kwint,
     NULL,
    &evaluate_keyword_int,
     NULL,
    {NULL, 0, NULL},
};

/** \internal Selection method data for string keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_str = {
    "kw_str", GROUP_VALUE, SMETH_SINGLEVAL,
    asize(smparams_keyword_str), smparams_keyword_str,
    &init_data_kwstr,
     NULL,
    &init_kwstr,
     NULL,
    &free_data_kwstr,
     NULL,
    &evaluate_keyword_str,
     NULL,
    {NULL, 0, NULL},
};


/********************************************************************
 * INTEGER KEYWORD EVALUATION
 ********************************************************************/

/*!
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns   Pointer to the allocated data (\ref t_methoddata_kwint).
 *
 * Allocates memory for a \ref t_methoddata_kwint structure.
 */
static void *
init_data_kwint(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_kwint *data;

    snew(data, 1);
    return data;
}

/*!
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_int).
 * \param[in] data  Should point to \ref t_methoddata_kwint.
 * \returns   0 (the initialization always succeeds).
 */
static int
init_kwint(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kwint *d = (t_methoddata_kwint *)data;

    d->v = param[0].val.u.i;
    d->n = param[1].val.nr;
    d->r = param[1].val.u.i;
    return 0;
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_kwint).
 *
 * Frees the memory allocated for t_methoddata_kwint::r.
 */
static void
free_data_kwint(void *data)
{
    t_methoddata_kwint *d = (t_methoddata_kwint *)data;

    sfree(d->v);
    sfree(d->r);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwint.
 *
 * Does a binary search to find which atoms match the ranges in the
 * \c t_methoddata_kwint structure for this selection.
 * Matching atoms are stored in \p out->u.g.
 */
static int
evaluate_keyword_int(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kwint *d = (t_methoddata_kwint *)data;
    int                 n, i, j, jmin, jmax;
    int                 val;

    out->u.g->isize = 0;
    n    = d->n;
    for (i = 0; i < g->isize; ++i)
    {
        val = d->v[i];
        if (d->r[0] > val || d->r[2*n-1] < val)
        {
            continue;
        }
        jmin = 0;
        jmax = n;
        while (jmax - jmin > 1)
        {
            j = jmin + (jmax - jmin) / 2;
            if (val < d->r[2*j])
            {
                jmax = j;
            }
            else
            {
                jmin = j;
                if (val <= d->r[2*j+1])
                {
                    break;
                }
                /* ++jmin;*/
            }
        }
        if (val <= d->r[2*jmin+1])
        {
            out->u.g->index[out->u.g->isize++] = g->index[i];
        }
    }
    return 0;
}


/********************************************************************
 * STRING KEYWORD EVALUATION
 ********************************************************************/

/*!
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns Pointer to the allocated data (\ref t_methoddata_kwstr).
 *
 * Allocates memory for a \ref t_methoddata_kwstr structure.
 */
static void *
init_data_kwstr(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_kwstr *data;

    snew(data, 1);
    return data;
}

/*!
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_str).
 * \param[in] data  Should point to \ref t_methoddata_kwstr.
 * \returns   0 (the initialization always succeeds).
 */
static int
init_kwstr(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)data;
    char               *buf;
    char               *s;
    int                 i;
    size_t              j;
    bool                bRegExp;

    d->v   = param[0].val.u.s;
    d->n   = param[1].val.nr;
    /* Return if this is not the first time */
    if (d->m)
    {
        return 0;
    }
    snew(d->m, d->n);
    for (i = 0; i < d->n; ++i)
    {
        s = param[1].val.u.s[i];
        bRegExp = FALSE;
        for (j = 0; j < strlen(s); ++j)
        {
            if (ispunct(s[j]) && s[j] != '?' && s[j] != '*')
            {
                bRegExp = TRUE;
                break;
            }
        }
        if (bRegExp)
        {
            snew(buf, strlen(s) + 3);
            sprintf(buf, "^%s$", s);
            if (regcomp(&d->m[i].u.r, buf, REG_EXTENDED | REG_NOSUB))
            {
                bRegExp = FALSE;
                fprintf(stderr, "warning: will match '%s' as a simple string\n", s);
            }
            else
            {
                sfree(s);
            }
            sfree(buf);
        }
        if (!bRegExp)
        {
            d->m[i].u.s = s;
        }
        d->m[i].bRegExp = bRegExp;
    }
    sfree(param[1].val.u.s);
    return 0;
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_kwstr).
 *
 * Frees the memory allocated for t_methoddata_kwstr::val.
 */
static void
free_data_kwstr(void *data)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)data;
    int                 i;

    sfree(d->v);
    for (i = 0; i < d->n; ++i)
    {
        if (d->m[i].bRegExp)
        {
            regfree(&d->m[i].u.r);
        }
        else
        {
            sfree(d->m[i].u.s);
        }
    }
    sfree(d->m);
}


/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwstr.
 *
 * Does a linear search to find which atoms match the strings in the
 * \c t_methoddata_kwstr structure for this selection.
 * Wildcards are allowed in the strings.
 * Matching atoms are stored in \p out->u.g.
 */
static int
evaluate_keyword_str(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)data;
    int                 i, j;
    bool                bFound;

    out->u.g->isize = 0;
    for (i = 0; i < g->isize; ++i)
    {
        bFound = FALSE;
        for (j = 0; j < d->n && !bFound; ++j)
        {
            if (d->m[j].bRegExp)
            {
                if (!regexec(&d->m[j].u.r, d->v[i], 0, NULL, 0))
                {
                    bFound = TRUE;
                }
            }
            else
            {
                if (!fnmatch(d->m[j].u.s, d->v[i], 0))
                {
                    bFound = TRUE;
                }
            }
        }
        if (bFound)
        {
            out->u.g->index[out->u.g->isize++] = g->index[i];
        }
    }
    return 0;
}
