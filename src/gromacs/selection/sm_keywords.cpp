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
 * Implements internal selection methods for numeric and string keyword
 * evaluation.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_selection
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>  /*old Mac needs types before regex.h*/
#endif
#ifdef HAVE_REGEX_H
#include <regex.h>
#define USE_REGEX
#endif

#include <string>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/string2.h"

#include "gromacs/selection/selmethod.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/messagestringcollector.h"
#include "gromacs/utility/stringutil.h"

#include "keywords.h"
#include "parsetree.h"
#include "scanner.h"
#include "selelem.h"

/** Allocates data for integer keyword evaluation. */
static void *
init_data_kwint(int npar, gmx_ana_selparam_t *param);
/** Allocates data for real keyword evaluation. */
static void *
init_data_kwreal(int npar, gmx_ana_selparam_t *param);
/** Allocates data for string keyword evaluation. */
static void *
init_data_kwstr(int npar, gmx_ana_selparam_t *param);
/** Initializes data for integer keyword evaluation. */
static void
init_kwint(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Initializes data for real keyword evaluation. */
static void
init_kwreal(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Initializes data for string keyword evaluation. */
static void
init_kwstr(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the memory allocated for string keyword evaluation. */
static void
free_data_kwstr(void *data);
/** Evaluates integer selection keywords. */
static void
evaluate_keyword_int(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates real selection keywords. */
static void
evaluate_keyword_real(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                      gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates string selection keywords. */
static void
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
 * Data structure for real keyword expression evaluation.
 */
typedef struct t_methoddata_kwreal
{
    /** Array of values for the keyword. */
    real              *v;
    /** Number of ranges in the \p r array. */
    int                n;
    /*! \brief
     * Array of sorted ranges to match against.
     *
     * Each range is made of two values, giving the endpoints (inclusive).
     * This field stores the pointer to the ranges allocated by the
     * parameter parser; see \ref SPAR_RANGES for more information.
     */
    real              *r;
} t_methoddata_kwreal;

/*! \internal \brief
 * Data structure for string keyword expression evaluation.
 */
typedef struct t_methoddata_kwstr
{
    /** Matching type for the strings. */
    gmx::SelectionStringMatchType       matchType;
    /** Array of values for the keyword. */
    char             **v;
    /** Number of elements in the \p val array. */
    int                n;
    /*! \internal \brief
     * Array of strings/regular expressions to match against.
     */
    struct t_methoddata_kwstr_match {
        /** true if the expression is a regular expression, false otherwise. */
        bool           bRegExp;
        /** The value to match against. */
        union {
#ifdef USE_REGEX
            /** Compiled regular expression if \p bRegExp is true. */
            regex_t    r;
#endif
            /** The string if \p bRegExp is false; */
            const char *s;
        }              u;
    }                 *m;
    /**< Array of strings/regular expressions to match against.*/
} t_methoddata_kwstr;

/** Parameters for integer keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_int[] = {
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_ATOMVAL},
    {NULL, {INT_VALUE, -1, {NULL}}, NULL, SPAR_RANGES | SPAR_VARNUM},
};

/** Parameters for real keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_real[] = {
    {NULL, {REAL_VALUE, -1, {NULL}}, NULL, SPAR_ATOMVAL | SPAR_DYNAMIC},
    {NULL, {REAL_VALUE, -1, {NULL}}, NULL, SPAR_RANGES | SPAR_VARNUM},
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
     NULL,
     NULL,
    &evaluate_keyword_int,
     NULL,
    {NULL, 0, NULL},
};

/** \internal Selection method data for real keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_real = {
    "kw_real", GROUP_VALUE, SMETH_SINGLEVAL,
    asize(smparams_keyword_real), smparams_keyword_real,
    &init_data_kwreal,
     NULL,
    &init_kwreal,
     NULL,
     NULL,
     NULL,
    &evaluate_keyword_real,
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

/** Initializes keyword evaluation for an arbitrary group. */
static void
init_kweval(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Initializes output for keyword evaluation in an arbitrary group. */
static void
init_output_kweval(t_topology *top, gmx_ana_selvalue_t *out, void *data);
/** Frees the data allocated for keyword evaluation in an arbitrary group. */
static void
free_data_kweval(void *data);
/** Initializes frame evaluation for keyword evaluation in an arbitrary group. */
static void
init_frame_kweval(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/** Evaluates keywords in an arbitrary group. */
static void
evaluate_kweval(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/*! \internal \brief
 * Data structure for keyword evaluation in arbitrary groups.
 */
typedef struct
{
    /** Wrapped keyword method for evaluating the values. */
    gmx_ana_selmethod_t  *kwmethod;
    /** Method data for \p kwmethod. */
    void                 *kwmdata;
    /** Group in which \p kwmethod should be evaluated. */
    gmx_ana_index_t       g;
} t_methoddata_kweval;

/** Parameters for keyword evaluation in an arbitrary group. */
static gmx_ana_selparam_t smparams_kweval[] = {
    {NULL,   {GROUP_VALUE, 1, {NULL}}, NULL, SPAR_DYNAMIC},
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
 */
static void
init_kwint(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kwint *d = (t_methoddata_kwint *)data;

    d->v = param[0].val.u.i;
    d->n = param[1].val.nr;
    d->r = param[1].val.u.i;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwint.
 *
 * Does a binary search to find which atoms match the ranges in the
 * \c t_methoddata_kwint structure for this selection.
 * Matching atoms are stored in \p out->u.g.
 */
static void
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
}


/********************************************************************
 * REAL KEYWORD EVALUATION
 ********************************************************************/

/*!
 * \param[in] npar  Not used.
 * \param     param Not used.
 * \returns   Pointer to the allocated data (\ref t_methoddata_kwreal).
 *
 * Allocates memory for a \ref t_methoddata_kwreal structure.
 */
static void *
init_data_kwreal(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_kwreal *data;

    snew(data, 1);
    return data;
}

/*!
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_real).
 * \param[in] data  Should point to \ref t_methoddata_kwreal.
 * \returns   0 (the initialization always succeeds).
 */
static void
init_kwreal(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kwreal *d = (t_methoddata_kwreal *)data;

    d->v = param[0].val.u.r;
    d->n = param[1].val.nr;
    d->r = param[1].val.u.r;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kwreal.
 *
 * Does a binary search to find which atoms match the ranges in the
 * \c t_methoddata_kwreal structure for this selection.
 * Matching atoms are stored in \p out->u.g.
 */
static void
evaluate_keyword_real(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kwreal *d = (t_methoddata_kwreal *)data;
    int                  n, i, j, jmin, jmax;
    real                 val;

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
    data->matchType = gmx::eStringMatchType_Auto;
    return data;
}

/*!
 * \param[in,out] sel   Selection element to initialize.
 * \param[in]     matchType  Method to use to match string values.
 *
 * Sets the string matching method for string keyword matching.
 */
void
_gmx_selelem_set_kwstr_match_type(const gmx::SelectionTreeElementPointer &sel,
                                  gmx::SelectionStringMatchType matchType)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)sel->u.expr.mdata;

    if (sel->type != SEL_EXPRESSION || !sel->u.expr.method
        || sel->u.expr.method->name != sm_keyword_str.name)
    {
        return;
    }
    d->matchType = matchType;
}

/*!
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_keyword_str).
 * \param[in] data  Should point to \ref t_methoddata_kwstr.
 */
static void
init_kwstr(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)data;
    int                 i;

    d->v   = param[0].val.u.s;
    d->n   = param[1].val.nr;
    /* Return if this is not the first time */
    if (d->m)
    {
        return;
    }
    snew(d->m, d->n);
    for (i = 0; i < d->n; ++i)
    {
        const char *s = param[1].val.u.s[i];
        bool bRegExp = (d->matchType == gmx::eStringMatchType_RegularExpression);
        if (d->matchType == gmx::eStringMatchType_Auto)
        {
            for (size_t j = 0; j < strlen(s); ++j)
            {
                if (ispunct(s[j]) && s[j] != '?' && s[j] != '*')
                {
                    bRegExp = true;
                    break;
                }
            }
        }
        if (bRegExp)
        {
#ifdef USE_REGEX
            std::string buf(gmx::formatString("^%s$", s));
            if (regcomp(&d->m[i].u.r, buf.c_str(), REG_EXTENDED | REG_NOSUB) != 0)
            {
                GMX_THROW(gmx::InvalidInputError(gmx::formatString(
                                "Error in regular expression \"%s\"", s)));
            }
#else
            GMX_THROW(gmx::InvalidInputError(gmx::formatString(
                            "No regular expression support, cannot match \"%s\"", s)));
#endif
        }
        else
        {
            d->m[i].u.s = s;
        }
        d->m[i].bRegExp = bRegExp;
    }
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

    for (i = 0; i < d->n; ++i)
    {
        if (d->m[i].bRegExp)
        {
#ifdef USE_REGEX
            /* This branch should only be taken if regular expressions
             * are available, but the ifdef is still needed. */
            regfree(&d->m[i].u.r);
#endif
        }
    }
    sfree(d->m);
    sfree(d);
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
static void
evaluate_keyword_str(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kwstr *d = (t_methoddata_kwstr *)data;
    int                 i, j;
    bool                bFound;

    out->u.g->isize = 0;
    for (i = 0; i < g->isize; ++i)
    {
        bFound = false;
        for (j = 0; j < d->n && !bFound; ++j)
        {
            if (d->m[j].bRegExp)
            {
#ifdef USE_REGEX
                /* This branch should only be taken if regular expressions
                 * are available, but the ifdef is still needed. */
                bFound = (regexec(&d->m[j].u.r, d->v[i], 0, NULL, 0) == 0);
#endif
            }
            else if (d->matchType == gmx::eStringMatchType_Exact)
            {
                bFound = (strcmp(d->m[j].u.s, d->v[i]) == 0);
            }
            else
            {
                bFound = (gmx_wcmatch(d->m[j].u.s, d->v[i]) == 0);
            }
        }
        if (bFound)
        {
            out->u.g->index[out->u.g->isize++] = g->index[i];
        }
    }
}


/********************************************************************
 * KEYWORD EVALUATION FOR ARBITRARY GROUPS
 ********************************************************************/

/*!
 * \param[in] top   Not used.
 * \param[in] npar  Not used.
 * \param[in] param Not used.
 * \param[in] data  Should point to \ref t_methoddata_kweval.
 * \returns   0 on success, a non-zero error code on return.
 *
 * Calls the initialization method of the wrapped keyword.
 */
static void
init_kweval(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_kweval *d = (t_methoddata_kweval *)data;

    d->kwmethod->init(top, 0, NULL, d->kwmdata);
}

/*!
 * \param[in]     top   Not used.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_kweval.
 * \returns       0 for success.
 */
static void
init_output_kweval(t_topology *top, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kweval *d = (t_methoddata_kweval *)data;

    out->nr = d->g.isize;
    _gmx_selvalue_reserve(out, out->nr);
}

/*!
 * \param data Data to free (should point to a \c t_methoddata_kweval).
 *
 * Frees the memory allocated for all the members of \c t_methoddata_kweval.
 */
static void
free_data_kweval(void *data)
{
    t_methoddata_kweval *d = (t_methoddata_kweval *)data;

    _gmx_selelem_free_method(d->kwmethod, d->kwmdata);
    sfree(d);
}

/*!
 * \param[in]  top  Topology.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \ref t_methoddata_kweval.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Creates a lookup structure that enables fast queries of whether a point
 * is within the solid angle or not.
 */
static void
init_frame_kweval(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data)
{
    t_methoddata_kweval *d = (t_methoddata_kweval *)data;

    d->kwmethod->init_frame(top, fr, pbc, d->kwmdata);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_kweval.
 *
 * Calls the evaluation function of the wrapped keyword with the given
 * parameters, with the exception of using \c t_methoddata_kweval::g for the
 * evaluation group.
 */
static void
evaluate_kweval(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_kweval *d = (t_methoddata_kweval *)data;

    d->kwmethod->update(top, fr, pbc, &d->g, out, d->kwmdata);
}

/*!
 * \param[in]   method  Keyword selection method to evaluate.
 * \param[in]   params  Parameter that gives the group to evaluate \p method in.
 * \param[in]   scanner Scanner data structure.
 * \returns     Pointer to the created selection element (NULL on error).
 *
 * Creates a \ref SEL_EXPRESSION selection element that evaluates the keyword
 * method given by \p method in the group given by \p param.
 *
 * The name of \p param should be empty.
 */
gmx::SelectionTreeElementPointer
_gmx_sel_init_keyword_evaluator(gmx_ana_selmethod_t *method,
                                const gmx::SelectionParserParameterList &params,
                                void *scanner)
{
    gmx::MessageStringCollector *errors = _gmx_sel_lexer_error_reporter(scanner);
    char  buf[1024];
    sprintf(buf, "In evaluation of '%s'", method->name);
    gmx::MessageStringContext   context(errors, buf);

    if ((method->flags & (SMETH_SINGLEVAL | SMETH_VARNUMVAL))
        || method->outinit || method->pupdate)
    {
        GMX_THROW(gmx::InternalError(
                "Unsupported keyword method for arbitrary group evaluation"));
    }

    gmx::SelectionTreeElementPointer sel(
            new gmx::SelectionTreeElement(SEL_EXPRESSION));
    _gmx_selelem_set_method(sel, method, scanner);

    t_methoddata_kweval  *data;
    snew(data, 1);
    data->kwmethod = sel->u.expr.method;
    data->kwmdata  = sel->u.expr.mdata;
    gmx_ana_index_clear(&data->g);

    snew(sel->u.expr.method, 1);
    memcpy(sel->u.expr.method, data->kwmethod, sizeof(gmx_ana_selmethod_t));
    sel->u.expr.method->flags       |= SMETH_VARNUMVAL;
    sel->u.expr.method->init_data    = NULL;
    sel->u.expr.method->set_poscoll  = NULL;
    sel->u.expr.method->init         = method->init ? &init_kweval : NULL;
    sel->u.expr.method->outinit      = &init_output_kweval;
    sel->u.expr.method->free         = &free_data_kweval;
    sel->u.expr.method->init_frame   = method->init_frame ? &init_frame_kweval : NULL;
    sel->u.expr.method->update       = &evaluate_kweval;
    sel->u.expr.method->pupdate      = NULL;
    sel->u.expr.method->nparams      = asize(smparams_kweval);
    sel->u.expr.method->param        = smparams_kweval;
    _gmx_selelem_init_method_params(sel, scanner);
    sel->u.expr.mdata = data;

    sel->u.expr.method->param[0].val.u.g = &data->g;

    if (!_gmx_sel_parse_params(params, sel->u.expr.method->nparams,
                               sel->u.expr.method->param, sel, scanner))
    {
        return gmx::SelectionTreeElementPointer();
    }
    return sel;
}
