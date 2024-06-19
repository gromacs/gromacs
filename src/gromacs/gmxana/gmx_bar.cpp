/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdlib/energyoutput.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/trajectory/energyframe.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t;


/* Structure for the names of lambda vector components */
typedef struct lambda_components_t
{
    char** names; /* Array of strings with names for the lambda vector
                     components */
    int N;        /* The number of components */
    int Nalloc;   /* The number of allocated components */
} lambda_components_t;

/* Structure for a lambda vector or a dhdl derivative direction */
typedef struct lambda_vec_t
{
    double* val;                   /* The lambda vector component values. Only valid if
                                      dhdl == -1 */
    int dhdl;                      /* The coordinate index for the derivative described by this
                                      structure, or -1 */
    const lambda_components_t* lc; /* the associated lambda_components
                                      structure */
    int index;                     /* The state number (init-lambda-state) of this lambda
                                      vector, if known. If not, it is set to -1 */
} lambda_vec_t;

/* the dhdl.xvg data from a simulation */
typedef struct xvg_t
{
    const char*   filename;
    int           nset;         /* number of lambdas, including dhdl */
    int*          np;           /* number of data points (du or hists) per lambda */
    double        temp;         /* temperature */
    lambda_vec_t* lambda;       /* the lambdas (of first index for y). */
    double*       t;            /* the times (of second index for y) */
    double**      y;            /* the dU values. y[0] holds the derivative, while
                                   further ones contain the energy differences between
                                   the native lambda and the 'foreign' lambdas. */
    lambda_vec_t native_lambda; /* the native lambda */

} xvg_t;


typedef struct hist_t
{
    unsigned int* bin[2]; /* the (forward + reverse) histogram values */
    double        dx[2];  /* the histogram spacing. The reverse
                             dx is the negative of the forward dx.*/
    int64_t x0[2];        /* the (forward + reverse) histogram start
                                     point(s) as int */

    int     nbin[2]; /* the (forward+reverse) number of bins */
    int64_t sum;     /* the total number of counts. Must be
                                the same for forward + reverse.  */
    int nhist;       /* number of hist datas (forward or reverse) */

    double start_time, delta_time; /* start time, end time of histogram */
} hist_t;


/* an aggregate of samples for partial free energy calculation */
typedef struct samples_t
{
    lambda_vec_t* native_lambda;  /* pointer to native lambda vector */
    lambda_vec_t* foreign_lambda; /* pointer to foreign lambda vector */
    double        temp;           /* the temperature */
    gmx_bool      derivative;     /* whether this sample is a derivative */

    /* The samples come either as either delta U lists: */
    int     ndu;                    /* the number of delta U samples */
    double* du;                     /* the delta u's */
    double* t;                      /* the times associated with those samples, or: */
    double  start_time, delta_time; /*start time and delta time for linear time*/

    /* or as histograms: */
    hist_t* hist; /* a histogram */

    /* allocation data: (not NULL for data 'owned' by this struct) */
    double* du_alloc;  /* allocated delta u arrays  */
    size_t  ndu_alloc; /* pre-allocated sizes */

    int64_t     ntot;     /* total number of samples */
    const char* filename; /* the file name this sample comes from */
} samples_t;

/* a sample range (start to end for du-style data, or boolean
    for both du-style data and histograms */
typedef struct sample_range_t
{
    int      start, end; /* start and end index for du style data */
    gmx_bool use;        /* whether to use this sample */

    samples_t* s; /* the samples this range belongs to */
} sample_range_t;


/* a collection of samples for a partial free energy calculation
    (i.e. the collection of samples from one native lambda to one
    foreign lambda) */
typedef struct sample_coll_t
{
    lambda_vec_t* native_lambda;  /* these should be the same for all samples
                                     in the histogram */
    lambda_vec_t* foreign_lambda; /* collection */
    double        temp;           /* the temperature */

    int             nsamples;       /* the number of samples */
    samples_t**     s;              /* the samples themselves */
    sample_range_t* r;              /* the sample ranges */
    int             nsamples_alloc; /* number of allocated samples */

    int64_t ntot; /* total number of samples in the ranges of
                             this collection */

    struct sample_coll_t *next, *prev; /* next and previous in the list */
} sample_coll_t;

/* all the samples associated with a lambda point */
typedef struct lambda_data_t
{
    lambda_vec_t* lambda; /* the native lambda (at start time if dynamic) */
    double        temp;   /* temperature */

    sample_coll_t* sc; /* the samples */

    sample_coll_t sc_head; /*the pre-allocated list head for the linked list.*/

    struct lambda_data_t *next, *prev; /* the next and prev in the list */
} lambda_data_t;

/* Top-level data structure of simulation data */
typedef struct sim_data_t
{
    lambda_data_t* lb;      /* a lambda data linked list */
    lambda_data_t  lb_head; /* The head element of the linked list */

    lambda_components_t lc; /* the allowed components of the lambda
                               vectors */
} sim_data_t;

/* Top-level data structure with calculated values. */
typedef struct
{
    sample_coll_t *a, *b; /* the simulation data */

    double dg;     /* the free energy difference */
    double dg_err; /* the free energy difference */

    double dg_disc_err;      /* discretization error */
    double dg_histrange_err; /* histogram range error */

    double sa;     /* relative entropy of b in state a */
    double sa_err; /* error in sa */
    double sb;     /* relative entropy of a in state b */
    double sb_err; /* error in sb */

    double dg_stddev;     /* expected dg stddev per sample */
    double dg_stddev_err; /* error in dg_stddev */
} barres_t;


/* Initialize a lambda_components structure */
static void lambda_components_init(lambda_components_t* lc)
{
    lc->N      = 0;
    lc->Nalloc = 2;
    snew(lc->names, lc->Nalloc);
}

/* Add a component to a lambda_components structure */
static void lambda_components_add(lambda_components_t* lc, const char* name, size_t name_length)
{
    while (lc->N + 1 > lc->Nalloc)
    {
        lc->Nalloc = (lc->Nalloc == 0) ? 2 : 2 * lc->Nalloc;
        srenew(lc->names, lc->Nalloc);
    }
    snew(lc->names[lc->N], name_length + 1);
    // GCC 12.1 has a false positive about the missing \0. But it is already there, nothing to worry about.
    GCC_DIAGNOSTIC_IGNORE("-Wstringop-truncation")
    std::strncpy(lc->names[lc->N], name, name_length);
    GCC_DIAGNOSTIC_RESET
    lc->N++;
}

/* check whether a component with index 'index' matches the given name, or
   is also NULL. Returns TRUE if this is the case.
   the string name does not need to end */
static gmx_bool lambda_components_check(const lambda_components_t* lc, int index, const char* name, size_t name_length)
{
    size_t len;
    if (!lc || index >= lc->N)
    {
        return FALSE;
    }
    if ((name == nullptr) && (lc->names[index] == nullptr))
    {
        return TRUE;
    }
    if (((name != nullptr) && (lc->names[index] == nullptr))
        || ((name == nullptr) && (lc->names[index] != nullptr)))
    {
        return FALSE;
    }
    GMX_RELEASE_ASSERT(
            (name != nullptr) || (name_length == 0),
            "If name is empty, the length of the substring to examine within it must be zero");
    len = std::strlen(lc->names[index]);
    if (len != name_length)
    {
        return FALSE;
    }
    if (name_length == 0)
    {
        // Everything matches a zero-length substring. This branch is
        // needed because name could in principle be nullptr.
        return TRUE;
    }
    return std::strncmp(lc->names[index], name, name_length) == 0;
}

/* Find the index of a given lambda component name, or -1 if not found */
static int lambda_components_find(const lambda_components_t* lc, const char* name, size_t name_length)
{
    int i;

    for (i = 0; i < lc->N; i++)
    {
        if (std::strncmp(lc->names[i], name, name_length) == 0)
        {
            return i;
        }
    }
    return -1;
}


/* initialize a lambda vector */
static void lambda_vec_init(lambda_vec_t* lv, const lambda_components_t* lc)
{
    snew(lv->val, lc->N);
    lv->index = -1;
    lv->dhdl  = -1;
    lv->lc    = lc;
}

static void lambda_vec_copy(lambda_vec_t* lv, const lambda_vec_t* orig)
{
    int i;

    lambda_vec_init(lv, orig->lc);
    lv->dhdl  = orig->dhdl;
    lv->index = orig->index;
    for (i = 0; i < lv->lc->N; i++)
    {
        lv->val[i] = orig->val[i];
    }
}

/* write a lambda vec to a preallocated string */
static void lambda_vec_print(const lambda_vec_t* lv, char* str, gmx_bool named)
{
    int i;

    str[0] = 0; /* reset the string */
    if (lv->dhdl < 0)
    {
        if (named)
        {
            str += sprintf(str, "delta H to ");
        }
        if (lv->lc->N > 1)
        {
            str += sprintf(str, "(");
        }
        for (i = 0; i < lv->lc->N; i++)
        {
            str += sprintf(str, "%g", lv->val[i]);
            if (i < lv->lc->N - 1)
            {
                str += sprintf(str, ", ");
            }
        }
        if (lv->lc->N > 1)
        {
            sprintf(str, ")");
        }
    }
    else
    {
        /* this lambda vector describes a derivative */
        str += sprintf(str, "dH/dl");
        if (std::strlen(lv->lc->names[lv->dhdl]) > 0)
        {
            sprintf(str, " (%s)", lv->lc->names[lv->dhdl]);
        }
    }
}

/* write a shortened version of the lambda vec to a preallocated string */
static void lambda_vec_print_short(const lambda_vec_t* lv, char* str)
{
    if (lv->index >= 0)
    {
        sprintf(str, "%6d", lv->index);
    }
    else
    {
        if (lv->dhdl < 0)
        {
            sprintf(str, "%6.3f", lv->val[0]);
        }
        else
        {
            sprintf(str, "dH/dl[%d]", lv->dhdl);
        }
    }
}

/* write an intermediate version of two lambda vecs to a preallocated string */
static void lambda_vec_print_intermediate(const lambda_vec_t* a, const lambda_vec_t* b, char* str)
{
    str[0] = 0;
    if ((a->index >= 0) && (b->index >= 0))
    {
        sprintf(str, "%6.3f", (a->index + b->index) / 2.0);
    }
    else
    {
        if ((a->dhdl < 0) && (b->dhdl < 0))
        {
            sprintf(str, "%6.3f", (a->val[0] + b->val[0]) / 2.0);
        }
    }
}


/* calculate and return the absolute difference in lambda vectors: c = |a-b|.
   a and b must describe non-derivative lambda points */
static double lambda_vec_abs_diff(const lambda_vec_t* a, const lambda_vec_t* b)
{
    int    i;
    double ret = 0.;

    if ((a->dhdl > 0) || (b->dhdl > 0))
    {
        gmx_fatal(
                FARGS,
                "Trying to calculate the difference between derivatives instead of lambda points");
    }
    if (a->lc != b->lc)
    {
        gmx_fatal(FARGS, "Trying to calculate the difference lambdas with differing basis set");
    }
    for (i = 0; i < a->lc->N; i++)
    {
        double df = a->val[i] - b->val[i];
        ret += df * df;
    }
    return std::sqrt(ret);
}


/* check whether two lambda vectors are the same */
static gmx_bool lambda_vec_same(const lambda_vec_t* a, const lambda_vec_t* b)
{
    int i;

    if (a->lc != b->lc)
    {
        return FALSE;
    }
    if (a->dhdl < 0)
    {
        for (i = 0; i < a->lc->N; i++)
        {
            if (!gmx_within_tol(a->val[i], b->val[i], 10 * GMX_REAL_EPS))
            {
                return FALSE;
            }
        }
        return TRUE;
    }
    else
    {
        /* they're derivatives, so we check whether the indices match */
        return (a->dhdl == b->dhdl);
    }
}

/* Compare the sort order of two foreign lambda vectors

    returns 1 if a is 'bigger' than b,
    returns 0 if they're the same,
    returns -1 if a is 'smaller' than b.*/
static int lambda_vec_cmp_foreign(const lambda_vec_t* a, const lambda_vec_t* b)
{
    int      i;
    double   norm_a = 0, norm_b = 0;
    gmx_bool different = FALSE;

    if (a->lc != b->lc)
    {
        gmx_fatal(FARGS, "Can't compare lambdas with differing basis sets");
    }
    /* if either one has an index we sort based on that */
    if ((a->index >= 0) || (b->index >= 0))
    {
        if (a->index == b->index)
        {
            return 0;
        }
        return (a->index > b->index) ? 1 : -1;
    }
    if (a->dhdl >= 0 || b->dhdl >= 0)
    {
        /* lambda vectors that are derivatives always sort higher than those
           without derivatives */
        if ((a->dhdl >= 0) != (b->dhdl >= 0))
        {
            return (a->dhdl >= 0) ? 1 : -1;
        }
        return 0;
    }

    /* neither has an index, so we can only sort on the lambda components,
       which is only valid if there is one component */
    for (i = 0; i < a->lc->N; i++)
    {
        if (!gmx_within_tol(a->val[i], b->val[i], 10 * GMX_REAL_EPS))
        {
            different = TRUE;
        }
        norm_a += a->val[i] * a->val[i];
        norm_b += b->val[i] * b->val[i];
    }
    if (!different)
    {
        return 0;
    }
    return (norm_a > norm_b) ? 1 : -1;
}

/* Compare the sort order of two native lambda vectors

    returns 1 if a is 'bigger' than b,
    returns 0 if they're the same,
    returns -1 if a is 'smaller' than b.*/
static int lambda_vec_cmp_native(const lambda_vec_t* a, const lambda_vec_t* b)
{
    if (a->lc != b->lc)
    {
        gmx_fatal(FARGS, "Can't compare lambdas with differing basis sets");
    }
    /* if either one has an index we sort based on that */
    if ((a->index >= 0) || (b->index >= 0))
    {
        if (a->index == b->index)
        {
            return 0;
        }
        return (a->index > b->index) ? 1 : -1;
    }
    /* neither has an index, so we can only sort on the lambda components,
       which is only valid if there is one component */
    if (a->lc->N > 1)
    {
        gmx_fatal(FARGS, "Can't compare lambdas with no index and > 1 component");
    }
    if (a->dhdl >= 0 || b->dhdl >= 0)
    {
        gmx_fatal(FARGS, "Can't compare native lambdas that are derivatives");
    }
    if (gmx_within_tol(a->val[0], b->val[0], 10 * GMX_REAL_EPS))
    {
        return 0;
    }
    return a->val[0] > b->val[0] ? 1 : -1;
}


static void hist_init(hist_t* h, int nhist, int* nbin)
{
    int i;
    if (nhist > 2)
    {
        gmx_fatal(FARGS, "histogram with more than two sets of data!");
    }
    for (i = 0; i < nhist; i++)
    {
        snew(h->bin[i], nbin[i]);
        h->x0[i]      = 0;
        h->nbin[i]    = nbin[i];
        h->start_time = h->delta_time = 0;
        h->dx[i]                      = 0;
    }
    h->sum   = 0;
    h->nhist = nhist;
}

static void xvg_init(xvg_t* ba)
{
    ba->filename = nullptr;
    ba->nset     = 0;
    ba->np       = nullptr;
    ba->y        = nullptr;
}

static void samples_init(samples_t*    s,
                         lambda_vec_t* native_lambda,
                         lambda_vec_t* foreign_lambda,
                         double        temp,
                         gmx_bool      derivative,
                         const char*   filename)
{
    s->native_lambda  = native_lambda;
    s->foreign_lambda = foreign_lambda;
    s->temp           = temp;
    s->derivative     = derivative;

    s->ndu        = 0;
    s->du         = nullptr;
    s->t          = nullptr;
    s->start_time = s->delta_time = 0;
    s->hist                       = nullptr;
    s->du_alloc                   = nullptr;
    s->ndu_alloc                  = 0;

    s->ntot     = 0;
    s->filename = filename;
}

static void sample_range_init(sample_range_t* r, samples_t* s)
{
    r->start = 0;
    r->end   = s->ndu;
    r->use   = TRUE;
    r->s     = nullptr;
}

static void sample_coll_init(sample_coll_t* sc, lambda_vec_t* native_lambda, lambda_vec_t* foreign_lambda, double temp)
{
    sc->native_lambda  = native_lambda;
    sc->foreign_lambda = foreign_lambda;
    sc->temp           = temp;

    sc->nsamples       = 0;
    sc->s              = nullptr;
    sc->r              = nullptr;
    sc->nsamples_alloc = 0;

    sc->ntot = 0;
    sc->next = sc->prev = nullptr;
}

static void sample_coll_destroy(sample_coll_t* sc)
{
    /* don't free the samples themselves */
    sfree(sc->r);
    sfree(sc->s);
}


static void lambda_data_init(lambda_data_t* l, lambda_vec_t* native_lambda, double temp)
{
    l->lambda = native_lambda;
    l->temp   = temp;

    l->next = nullptr;
    l->prev = nullptr;

    l->sc = &(l->sc_head);

    sample_coll_init(l->sc, native_lambda, nullptr, 0.);
    l->sc->next = l->sc;
    l->sc->prev = l->sc;
}

static void barres_init(barres_t* br)
{
    br->dg            = 0;
    br->dg_err        = 0;
    br->sa            = 0;
    br->sa_err        = 0;
    br->sb            = 0;
    br->sb_err        = 0;
    br->dg_stddev     = 0;
    br->dg_stddev_err = 0;

    br->a = nullptr;
    br->b = nullptr;
}


/* calculate the total number of samples in a sample collection */
static void sample_coll_calc_ntot(sample_coll_t* sc)
{
    int i;

    sc->ntot = 0;
    for (i = 0; i < sc->nsamples; i++)
    {
        if (sc->r[i].use)
        {
            if (sc->s[i]->hist)
            {
                sc->ntot += sc->s[i]->ntot;
            }
            else
            {
                sc->ntot += sc->r[i].end - sc->r[i].start;
            }
        }
    }
}


/* find the barsamples_t associated with a lambda that corresponds to
   a specific foreign lambda */
static sample_coll_t* lambda_data_find_sample_coll(lambda_data_t* l, lambda_vec_t* foreign_lambda)
{
    sample_coll_t* sc = l->sc->next;

    while (sc != l->sc)
    {
        if (lambda_vec_same(sc->foreign_lambda, foreign_lambda))
        {
            return sc;
        }
        sc = sc->next;
    }

    return nullptr;
}

/* insert li into an ordered list of lambda_colls */
static void lambda_data_insert_sample_coll(lambda_data_t* l, sample_coll_t* sc)
{
    sample_coll_t* scn = l->sc->next;
    while ((scn != l->sc))
    {
        if (lambda_vec_cmp_foreign(scn->foreign_lambda, sc->foreign_lambda) > 0)
        {
            break;
        }
        scn = scn->next;
    }
    /* now insert it before the found scn */
    sc->next        = scn;
    sc->prev        = scn->prev;
    scn->prev->next = sc;
    scn->prev       = sc;
}

/* insert li into an ordered list of lambdas */
static void lambda_data_insert_lambda(lambda_data_t* head, lambda_data_t* li)
{
    lambda_data_t* lc = head->next;
    while (lc != head)
    {
        if (lambda_vec_cmp_native(lc->lambda, li->lambda) > 0)
        {
            break;
        }
        lc = lc->next;
    }
    /* now insert ourselves before the found lc */
    li->next       = lc;
    li->prev       = lc->prev;
    lc->prev->next = li;
    lc->prev       = li;
}

/* insert a sample and a sample_range into a sample_coll. The
    samples are stored as a pointer, the range is copied. */
static void sample_coll_insert_sample(sample_coll_t* sc, samples_t* s, sample_range_t* r)
{
    /* first check if it belongs here */
    GMX_ASSERT(sc->next->s, "Next not properly initialized!");
    if (sc->temp != s->temp)
    {
        gmx_fatal(FARGS,
                  "Temperatures in files %s and %s are not the same!",
                  s->filename,
                  sc->next->s[0]->filename);
    }
    if (!lambda_vec_same(sc->native_lambda, s->native_lambda))
    {
        gmx_fatal(FARGS,
                  "Native lambda in files %s and %s are not the same (and they should be)!",
                  s->filename,
                  sc->next->s[0]->filename);
    }
    if (!lambda_vec_same(sc->foreign_lambda, s->foreign_lambda))
    {
        gmx_fatal(FARGS,
                  "Foreign lambda in files %s and %s are not the same (and they should be)!",
                  s->filename,
                  sc->next->s[0]->filename);
    }

    /* check if there's room */
    if ((sc->nsamples + 1) > sc->nsamples_alloc)
    {
        sc->nsamples_alloc = std::max(2 * sc->nsamples_alloc, 2);
        srenew(sc->s, sc->nsamples_alloc);
        srenew(sc->r, sc->nsamples_alloc);
    }
    sc->s[sc->nsamples] = s;
    sc->r[sc->nsamples] = *r;
    sc->nsamples++;

    sample_coll_calc_ntot(sc);
}

/* insert a sample into a lambda_list, creating the right sample_coll if
   neccesary */
static void lambda_data_list_insert_sample(lambda_data_t* head, samples_t* s)
{
    gmx_bool       found = FALSE;
    sample_coll_t* sc;
    sample_range_t r;

    lambda_data_t* l = head->next;

    /* first search for the right lambda_data_t */
    while (l != head)
    {
        if (lambda_vec_same(l->lambda, s->native_lambda))
        {
            found = TRUE;
            break;
        }
        l = l->next;
    }

    if (!found)
    {
        snew(l, 1);                                     /* allocate a new one */
        lambda_data_init(l, s->native_lambda, s->temp); /* initialize it */
        lambda_data_insert_lambda(head, l);             /* add it to the list */
    }

    /* now look for a sample collection */
    sc = lambda_data_find_sample_coll(l, s->foreign_lambda);
    if (!sc)
    {
        snew(sc, 1); /* allocate a new one */
        sample_coll_init(sc, s->native_lambda, s->foreign_lambda, s->temp);
        lambda_data_insert_sample_coll(l, sc);
    }

    /* now insert the samples into the sample coll */
    sample_range_init(&r, s);
    sample_coll_insert_sample(sc, s, &r);
}


/* make a histogram out of a sample collection */
static void sample_coll_make_hist(sample_coll_t* sc, std::vector<int>* bin, double* dx, double* xmin, int nbin_default)
{
    int      i, j, k;
    gmx_bool dx_set   = FALSE;
    gmx_bool xmin_set = FALSE;

    gmx_bool xmax_set      = FALSE;
    gmx_bool xmax_set_hard = FALSE; /* whether the xmax is bounded by the
                                       limits of a histogram */
    double xmax = -1;

    /* first determine dx and xmin; try the histograms */
    for (i = 0; i < sc->nsamples; i++)
    {
        if (sc->s[i]->hist)
        {
            hist_t* hist = sc->s[i]->hist;
            for (k = 0; k < hist->nhist; k++)
            {
                double hdx      = hist->dx[k];
                double xmax_now = (hist->x0[k] + hist->nbin[k]) * hdx;

                /* we use the biggest dx*/
                if ((!dx_set) || hist->dx[0] > *dx)
                {
                    dx_set = TRUE;
                    *dx    = hist->dx[0];
                }
                if ((!xmin_set) || (hist->x0[k] * hdx) < *xmin)
                {
                    xmin_set = TRUE;
                    *xmin    = (hist->x0[k] * hdx);
                }

                if ((!xmax_set) || (xmax_now > xmax && !xmax_set_hard))
                {
                    xmax_set = TRUE;
                    xmax     = xmax_now;
                    if (hist->bin[k][hist->nbin[k] - 1] != 0)
                    {
                        xmax_set_hard = TRUE;
                    }
                }
                if (hist->bin[k][hist->nbin[k] - 1] != 0 && (xmax_now < xmax))
                {
                    xmax_set_hard = TRUE;
                    xmax          = xmax_now;
                }
            }
        }
    }
    /* and the delta us */
    for (i = 0; i < sc->nsamples; i++)
    {
        if (sc->s[i]->ndu > 0)
        {
            /* determine min and max */
            int    starti  = sc->r[i].start;
            int    endi    = sc->r[i].end;
            double du_xmin = sc->s[i]->du[starti];
            double du_xmax = sc->s[i]->du[starti];
            for (j = starti + 1; j < endi; j++)
            {
                if (sc->s[i]->du[j] < du_xmin)
                {
                    du_xmin = sc->s[i]->du[j];
                }
                if (sc->s[i]->du[j] > du_xmax)
                {
                    du_xmax = sc->s[i]->du[j];
                }
            }

            /* and now change the limits */
            if ((!xmin_set) || (du_xmin < *xmin))
            {
                xmin_set = TRUE;
                *xmin    = du_xmin;
            }
            if ((!xmax_set) || ((du_xmax > xmax) && !xmax_set_hard))
            {
                xmax_set = TRUE;
                xmax     = du_xmax;
            }
        }
    }

    if (!xmax_set || !xmin_set)
    {
        bin->clear();
        return;
    }


    if (!dx_set)
    {
        bin->resize(nbin_default);
        *dx = (xmax - (*xmin)) / ((bin->size()) - 2); /* -2 because we want the last bin to
                                                   be 0, and we count from 0 */
    }
    else
    {
        bin->resize(static_cast<int>((xmax - (*xmin)) / (*dx)));
    }

    /* reset the histogram */
    std::fill(bin->begin(), bin->end(), 0);

    /* now add the actual data */
    for (i = 0; i < sc->nsamples; i++)
    {
        if (sc->s[i]->hist)
        {
            hist_t* hist = sc->s[i]->hist;
            for (k = 0; k < hist->nhist; k++)
            {
                double hdx       = hist->dx[k];
                double xmin_hist = hist->x0[k] * hdx;
                for (j = 0; j < hist->nbin[k]; j++)
                {
                    /* calculate the bin corresponding to the middle of the
                       original bin */
                    double x     = hdx * (j + 0.5) + xmin_hist;
                    int    binnr = static_cast<int>((x - (*xmin)) / (*dx));

                    if (binnr >= gmx::ssize(*bin) || binnr < 0)
                    {
                        binnr = (bin->size()) - 1;
                    }

                    (*bin)[binnr] += hist->bin[k][j];
                }
            }
        }
        else
        {
            int starti = sc->r[i].start;
            int endi   = sc->r[i].end;
            for (j = starti; j < endi; j++)
            {
                int binnr = static_cast<int>((sc->s[i]->du[j] - (*xmin)) / (*dx));
                if (binnr >= gmx::ssize(*bin) || binnr < 0)
                {
                    binnr = (bin->size()) - 1;
                }

                (*bin)[binnr]++;
            }
        }
    }
}

/* write a collection of histograms to a file */
static void sim_data_histogram(sim_data_t* sd, const char* filename, int nbin_default, const gmx_output_env_t* oenv)
{
    char                     label_x[STRLEN];
    const char *             dhdl = "dH/d\\lambda", *deltag = "\\DeltaH", *lambda = "\\lambda";
    const char*              title   = "N(\\DeltaH)";
    const char*              label_y = "Samples";
    FILE*                    fp;
    lambda_data_t*           bl;
    std::vector<std::string> setnames;
    gmx_bool                 first_set = FALSE;
    /* histogram data: */
    std::vector<int> hist;
    double           dx      = 0;
    double           minval  = 0;
    lambda_data_t*   bl_head = sd->lb;

    printf("\nWriting histogram to %s\n", filename);
    sprintf(label_x, "\\DeltaH (%s)", unit_energy);

    fp = xvgropen_type(filename, title, label_x, label_y, exvggtXNY, oenv);

    /* first get all the set names */
    bl = bl_head->next;
    /* iterate over all lambdas */
    while (bl != bl_head)
    {
        sample_coll_t* sc = bl->sc->next;

        /* iterate over all samples */
        while (sc != bl->sc)
        {
            char buf[STRLEN], buf2[STRLEN];

            if (sc->foreign_lambda->dhdl < 0)
            {
                lambda_vec_print(sc->native_lambda, buf, FALSE);
                lambda_vec_print(sc->foreign_lambda, buf2, FALSE);
                setnames.emplace_back(gmx::formatString(
                        "N(%s(%s=%s) | %s=%s)", deltag, lambda, buf2, lambda, buf));
            }
            else
            {
                lambda_vec_print(sc->native_lambda, buf, FALSE);
                setnames.emplace_back(gmx::formatString("N(%s | %s=%s)", dhdl, lambda, buf));
            }
            sc = sc->next;
        }

        bl = bl->next;
    }
    xvgrLegend(fp, setnames, oenv);


    /* now make the histograms */
    bl = bl_head->next;
    /* iterate over all lambdas */
    while (bl != bl_head)
    {
        sample_coll_t* sc = bl->sc->next;

        /* iterate over all samples */
        while (sc != bl->sc)
        {
            if (!first_set)
            {
                xvgrNewDataset(fp, 0, {}, oenv);
            }

            sample_coll_make_hist(sc, &hist, &dx, &minval, nbin_default);

            for (gmx::Index i = 0; i < gmx::ssize(hist); i++)
            {
                double xmin = i * dx + minval;
                double xmax = (i + 1) * dx + minval;

                fprintf(fp, "%g %d\n%g %d\n", xmin, hist[i], xmax, hist[i]);
            }

            first_set = FALSE;
            sc        = sc->next;
        }

        bl = bl->next;
    }

    xvgrclose(fp);
}

static int snprint_lambda_vec(char* str, int sz, const char* label, lambda_vec_t* lambda)
{
    int n = 0;

    n += snprintf(str + n, sz - n, "lambda vector [%s]: ", label);
    if (lambda->index >= 0)
    {
        n += snprintf(str + n, sz - n, " init-lambda-state=%d", lambda->index);
    }
    if (lambda->dhdl >= 0)
    {
        n += snprintf(str + n, sz - n, " dhdl index=%d", lambda->dhdl);
    }
    else
    {
        int i;
        for (i = 0; i < lambda->lc->N; i++)
        {
            n += snprintf(str + n, sz - n, " (%s) l=%g", lambda->lc->names[i], lambda->val[i]);
        }
    }
    return n;
}

/* create a collection (array) of barres_t object given a ordered linked list
   of barlamda_t sample collections */
static barres_t* barres_list_create(sim_data_t* sd, int* nres, gmx_bool use_dhdl)
{
    lambda_data_t* bl;
    int            nlambda = 0;
    barres_t*      res;
    gmx_bool       dhdl    = FALSE;
    gmx_bool       first   = TRUE;
    lambda_data_t* bl_head = sd->lb;

    /* first count the lambdas */
    bl = bl_head->next;
    while (bl != bl_head)
    {
        nlambda++;
        bl = bl->next;
    }
    snew(res, nlambda - 1);

    /* next put the right samples in the res */
    *nres = 0;
    bl    = bl_head->next->next; /* we start with the second one. */
    while (bl != bl_head)
    {
        sample_coll_t *sc, *scprev;
        barres_t*      br = &(res[*nres]);
        /* there is always a previous one. we search for that as a foreign
           lambda: */
        scprev = lambda_data_find_sample_coll(bl->prev, bl->lambda);
        sc     = lambda_data_find_sample_coll(bl, bl->prev->lambda);

        barres_init(br);

        if (use_dhdl)
        {
            /* we use dhdl */

            scprev = lambda_data_find_sample_coll(bl->prev, bl->prev->lambda);
            sc     = lambda_data_find_sample_coll(bl, bl->lambda);

            if (first)
            {
                printf("\nWARNING: Using the derivative data (dH/dlambda) to extrapolate delta H "
                       "values.\nThis will only work if the Hamiltonian is linear in lambda.\n");
                dhdl = TRUE;
            }
            if (!dhdl)
            {
                gmx_fatal(FARGS,
                          "Some dhdl files contain only one value (dH/dl), while others \ncontain "
                          "multiple values (dH/dl and/or Delta H), will not proceed \nbecause of "
                          "possible inconsistencies.\n");
            }
        }
        else if (!scprev && !sc)
        {
            char descX[STRLEN], descY[STRLEN];
            snprint_lambda_vec(descX, STRLEN, "X", bl->prev->lambda);
            snprint_lambda_vec(descY, STRLEN, "Y", bl->lambda);

            gmx_fatal(FARGS,
                      "There is no path between the states X & Y below that is covered by foreign "
                      "lambdas:\ncannot proceed with BAR.\nUse thermodynamic integration of dH/dl "
                      "by calculating the averages of dH/dl\nwith gmx analyze and integrating "
                      "them.\nAlternatively, use the -extp option if (and only if) the "
                      "Hamiltonian\ndepends linearly on lambda, which is NOT normally the "
                      "case.\n\n%s\n%s\n",
                      descX,
                      descY);
        }

        /* normal delta H */
        if (!scprev)
        {
            char descX[STRLEN], descY[STRLEN];
            snprint_lambda_vec(descX, STRLEN, "X", bl->lambda);
            snprint_lambda_vec(descY, STRLEN, "Y", bl->prev->lambda);
            gmx_fatal(FARGS,
                      "Could not find a set for foreign lambda (state X below)\nin the files for "
                      "main lambda (state Y below)\n\n%s\n%s\n",
                      descX,
                      descY);
        }
        if (!sc)
        {
            char descX[STRLEN], descY[STRLEN];
            snprint_lambda_vec(descX, STRLEN, "X", bl->prev->lambda);
            snprint_lambda_vec(descY, STRLEN, "Y", bl->lambda);
            gmx_fatal(FARGS,
                      "Could not find a set for foreign lambda (state X below)\nin the files for "
                      "main lambda (state Y below)\n\n%s\n%s\n",
                      descX,
                      descY);
        }
        br->a = scprev;
        br->b = sc;

        first = FALSE;
        (*nres)++;
        bl = bl->next;
    }
    return res;
}

/* estimate the maximum discretization error */
static double barres_list_max_disc_err(barres_t* res, int nres)
{
    int    i, j;
    double disc_err = 0.;
    double delta_lambda;

    for (i = 0; i < nres; i++)
    {
        barres_t* br = &(res[i]);

        delta_lambda = lambda_vec_abs_diff(br->b->native_lambda, br->a->native_lambda);

        for (j = 0; j < br->a->nsamples; j++)
        {
            if (br->a->s[j]->hist)
            {
                double Wfac = 1.;
                if (br->a->s[j]->derivative)
                {
                    Wfac = delta_lambda;
                }

                disc_err = std::max(disc_err, Wfac * br->a->s[j]->hist->dx[0]);
            }
        }
        for (j = 0; j < br->b->nsamples; j++)
        {
            if (br->b->s[j]->hist)
            {
                double Wfac = 1.;
                if (br->b->s[j]->derivative)
                {
                    Wfac = delta_lambda;
                }
                disc_err = std::max(disc_err, Wfac * br->b->s[j]->hist->dx[0]);
            }
        }
    }
    return disc_err;
}


/* impose start and end times on a sample collection, updating sample_ranges */
static void sample_coll_impose_times(sample_coll_t* sc, double begin_t, double end_t)
{
    int i;
    for (i = 0; i < sc->nsamples; i++)
    {
        samples_t*      s = sc->s[i];
        sample_range_t* r = &(sc->r[i]);
        if (s->hist)
        {
            double end_time = s->hist->delta_time * s->hist->sum + s->hist->start_time;
            if (s->hist->start_time < begin_t || end_time > end_t)
            {
                r->use = FALSE;
            }
        }
        else
        {
            if (!s->t)
            {
                double end_time;
                if (s->start_time < begin_t)
                {
                    r->start = static_cast<int>((begin_t - s->start_time) / s->delta_time);
                }
                end_time = s->delta_time * s->ndu + s->start_time;
                if (end_time > end_t)
                {
                    r->end = static_cast<int>((end_t - s->start_time) / s->delta_time);
                }
            }
            else
            {
                int j;
                for (j = 0; j < s->ndu; j++)
                {
                    if (s->t[j] < begin_t)
                    {
                        r->start = j;
                    }

                    if (s->t[j] >= end_t)
                    {
                        r->end = j;
                        break;
                    }
                }
            }
            if (r->start > r->end)
            {
                r->use = FALSE;
            }
        }
    }
    sample_coll_calc_ntot(sc);
}

static void sim_data_impose_times(sim_data_t* sd, double begin, double end)
{
    double         first_t, last_t;
    double         begin_t, end_t;
    lambda_data_t* lc;
    lambda_data_t* head = sd->lb;
    int            j;

    if (begin <= 0 && end < 0)
    {
        return;
    }

    /* first determine the global start and end times */
    first_t = -1;
    last_t  = -1;
    lc      = head->next;
    while (lc != head)
    {
        sample_coll_t* sc = lc->sc->next;
        while (sc != lc->sc)
        {
            for (j = 0; j < sc->nsamples; j++)
            {
                double start_t, end_t;

                start_t = sc->s[j]->start_time;
                end_t   = sc->s[j]->start_time;
                if (sc->s[j]->hist)
                {
                    end_t += sc->s[j]->delta_time * sc->s[j]->hist->sum;
                }
                else
                {
                    if (sc->s[j]->t)
                    {
                        end_t = sc->s[j]->t[sc->s[j]->ndu - 1];
                    }
                    else
                    {
                        end_t += sc->s[j]->delta_time * sc->s[j]->ndu;
                    }
                }

                if (start_t < first_t || first_t < 0)
                {
                    first_t = start_t;
                }
                if (end_t > last_t)
                {
                    last_t = end_t;
                }
            }
            sc = sc->next;
        }
        lc = lc->next;
    }

    /* calculate the actual times */
    if (begin > 0)
    {
        begin_t = begin;
    }
    else
    {
        begin_t = first_t;
    }

    if (end > 0)
    {
        end_t = end;
    }
    else
    {
        end_t = last_t;
    }
    printf("\n   Samples in time interval: %.3f - %.3f\n", first_t, last_t);

    if (begin_t > end_t)
    {
        return;
    }
    printf("Removing samples outside of: %.3f - %.3f\n", begin_t, end_t);

    /* then impose them */
    lc = head->next;
    while (lc != head)
    {
        sample_coll_t* sc = lc->sc->next;
        while (sc != lc->sc)
        {
            sample_coll_impose_times(sc, begin_t, end_t);
            sc = sc->next;
        }
        lc = lc->next;
    }
}


/* create subsample i out of ni from an existing sample_coll */
static gmx_bool sample_coll_create_subsample(sample_coll_t* sc, sample_coll_t* sc_orig, int i, int ni)
{
    int j;

    int64_t ntot_start;
    int64_t ntot_end;
    int64_t ntot_so_far;

    *sc = *sc_orig; /* just copy all fields */

    /* allocate proprietary memory */
    snew(sc->s, sc_orig->nsamples);
    snew(sc->r, sc_orig->nsamples);

    /* copy the samples */
    for (j = 0; j < sc_orig->nsamples; j++)
    {
        sc->s[j] = sc_orig->s[j];
        sc->r[j] = sc_orig->r[j]; /* copy the ranges too */
    }

    /* now fix start and end fields */
    /* the casts avoid possible overflows */
    ntot_start = static_cast<int64_t>(sc_orig->ntot * static_cast<double>(i) / static_cast<double>(ni));
    ntot_end = static_cast<int64_t>(sc_orig->ntot * static_cast<double>(i + 1) / static_cast<double>(ni));
    ntot_so_far = 0;
    for (j = 0; j < sc->nsamples; j++)
    {
        int64_t ntot_add;
        int64_t new_start, new_end;

        if (sc->r[j].use)
        {
            if (sc->s[j]->hist)
            {
                ntot_add = sc->s[j]->hist->sum;
            }
            else
            {
                ntot_add = sc->r[j].end - sc->r[j].start;
            }
        }
        else
        {
            ntot_add = 0;
        }

        if (!sc->s[j]->hist)
        {
            if (ntot_so_far < ntot_start)
            {
                /* adjust starting point */
                new_start = sc->r[j].start + (ntot_start - ntot_so_far);
            }
            else
            {
                new_start = sc->r[j].start;
            }
            /* adjust end point */
            new_end = sc->r[j].start + (ntot_end - ntot_so_far);
            if (new_end > sc->r[j].end)
            {
                new_end = sc->r[j].end;
            }

            /* check if we're in range at all */
            if ((new_end < new_start) || (new_start > sc->r[j].end))
            {
                new_start = 0;
                new_end   = 0;
            }
            /* and write the new range */
            GMX_RELEASE_ASSERT(new_start <= std::numeric_limits<int>::max(),
                               "Value of 'new_start' too large for int converstion");
            GMX_RELEASE_ASSERT(new_end <= std::numeric_limits<int>::max(),
                               "Value of 'new_end' too large for int converstion");
            sc->r[j].start = static_cast<int>(new_start);
            sc->r[j].end   = static_cast<int>(new_end);
        }
        else
        {
            if (sc->r[j].use)
            {
                double overlap;
                double ntot_start_norm, ntot_end_norm;
                /* calculate the amount of overlap of the
                   desired range (ntot_start -- ntot_end) onto
                   the histogram range (ntot_so_far -- ntot_so_far+ntot_add)*/

                /* first calculate normalized bounds
                   (where 0 is the start of the hist range, and 1 the end) */
                ntot_start_norm = (ntot_start - ntot_so_far) / static_cast<double>(ntot_add);
                ntot_end_norm   = (ntot_end - ntot_so_far) / static_cast<double>(ntot_add);

                /* now fix the boundaries */
                ntot_start_norm = std::min(1.0, std::max(0.0, ntot_start_norm));
                ntot_end_norm   = std::max(0.0, std::min(1.0, ntot_end_norm));

                /* and calculate the overlap */
                overlap = ntot_end_norm - ntot_start_norm;

                if (overlap > 0.95) /* we allow for 5% slack */
                {
                    sc->r[j].use = TRUE;
                }
                else if (overlap < 0.05)
                {
                    sc->r[j].use = FALSE;
                }
                else
                {
                    return FALSE;
                }
            }
        }
        ntot_so_far += ntot_add;
    }
    sample_coll_calc_ntot(sc);

    return TRUE;
}

/* calculate minimum and maximum work values in sample collection */
static void sample_coll_min_max(sample_coll_t* sc, double Wfac, double* Wmin, double* Wmax)
{
    int i, j;

    *Wmin = std::numeric_limits<float>::max();
    *Wmax = -std::numeric_limits<float>::max();

    for (i = 0; i < sc->nsamples; i++)
    {
        samples_t*      s = sc->s[i];
        sample_range_t* r = &(sc->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for (j = r->start; j < r->end; j++)
                {
                    *Wmin = std::min(*Wmin, s->du[j] * Wfac);
                    *Wmax = std::max(*Wmax, s->du[j] * Wfac);
                }
            }
            else
            {
                int    hd = 0; /* determine the histogram direction: */
                double dx;
                if ((s->hist->nhist > 1) && (Wfac < 0))
                {
                    hd = 1;
                }
                dx = s->hist->dx[hd];

                for (j = s->hist->nbin[hd] - 1; j >= 0; j--)
                {
                    *Wmin = std::min(*Wmin, Wfac * (s->hist->x0[hd]) * dx);
                    *Wmax = std::max(*Wmax, Wfac * (s->hist->x0[hd]) * dx);
                    /* look for the highest value bin with values */
                    if (s->hist->bin[hd][j] > 0)
                    {
                        *Wmin = std::min(*Wmin, Wfac * (j + s->hist->x0[hd] + 1) * dx);
                        *Wmax = std::max(*Wmax, Wfac * (j + s->hist->x0[hd] + 1) * dx);
                        break;
                    }
                }
            }
        }
    }
}

/* Initialize a sim_data structure */
static void sim_data_init(sim_data_t* sd)
{
    /* make linked list */
    sd->lb       = &(sd->lb_head);
    sd->lb->next = sd->lb;
    sd->lb->prev = sd->lb;

    lambda_components_init(&(sd->lc));
}


static double calc_bar_sum(int n, const double* W, double Wfac, double sbMmDG)
{
    int    i;
    double sum;

    sum = 0;

    for (i = 0; i < n; i++)
    {
        sum += 1. / (1. + std::exp(Wfac * W[i] + sbMmDG));
    }

    return sum;
}

/* calculate the BAR average given a histogram

    if type== 0, calculate the best estimate for the average,
    if type==-1, calculate the minimum possible value given the histogram
    if type== 1, calculate the maximum possible value given the histogram */
static double calc_bar_sum_hist(const hist_t* hist, double Wfac, double sbMmDG, int type)
{
    double sum = 0.;
    int    i;
    int    maxbin;
    /* normalization factor multiplied with bin width and
       number of samples (we normalize through M): */
    double normdx = 1.;
    int    hd     = 0; /* determine the histogram direction: */
    double dx;

    if ((hist->nhist > 1) && (Wfac < 0))
    {
        hd = 1;
    }
    dx     = hist->dx[hd];
    maxbin = hist->nbin[hd] - 1;
    if (type == 1)
    {
        maxbin = hist->nbin[hd]; /* we also add whatever was out of range */
    }

    for (i = 0; i < maxbin; i++)
    {
        double x    = Wfac * ((i + hist->x0[hd]) + 0.5) * dx; /* bin middle */
        double pxdx = hist->bin[0][i] * normdx;               /* p(x)dx */

        sum += pxdx / (1. + std::exp(x + sbMmDG));
    }

    return sum;
}

static double calc_bar_lowlevel(sample_coll_t* ca, sample_coll_t* cb, double temp, double tol, int type)
{
    double kT, beta, M;
    int    i;
    double Wfac1, Wfac2, Wmin, Wmax;
    double DG0, DG1, DG2, dDG1;
    double n1, n2; /* numbers of samples as doubles */

    kT   = gmx::c_boltz * temp;
    beta = 1 / kT;

    /* count the numbers of samples */
    n1 = ca->ntot;
    n2 = cb->ntot;

    M = std::log(n1 / n2);

    /*if (!lambda_vec_same(ca->native_lambda, ca->foreign_lambda))*/
    if (ca->foreign_lambda->dhdl < 0)
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a
           multiplication factor.  */
        /*double delta_lambda=cb->native_lambda-ca->native_lambda;*/
        double delta_lambda = lambda_vec_abs_diff(cb->native_lambda, ca->native_lambda);
        if (cb->native_lambda->lc->N > 1)
        {
            gmx_fatal(FARGS, "Can't (yet) do multi-component dhdl interpolation");
        }

        Wfac1 = beta * delta_lambda;
        Wfac2 = -beta * delta_lambda;
    }

    if (beta < 1)
    {
        /* We print the output both in kT and kJ/mol.
         * Here we determine DG in kT, so when beta < 1
         * the precision has to be increased.
         */
        tol *= beta;
    }

    /* Calculate minimum and maximum work to give an initial estimate of
     * delta G  as their average.
     */
    {
        double Wmin1, Wmin2, Wmax1, Wmax2;
        sample_coll_min_max(ca, Wfac1, &Wmin1, &Wmax1);
        sample_coll_min_max(cb, Wfac2, &Wmin2, &Wmax2);

        Wmin = std::min(Wmin1, Wmin2);
        Wmax = std::max(Wmax1, Wmax2);
    }

    DG0 = Wmin;
    DG2 = Wmax;

    if (debug)
    {
        fprintf(debug, "DG %9.5f %9.5f\n", DG0, DG2);
    }
    /* We approximate by bisection: given our initial estimates
       we keep checking whether the halfway point is greater or
       smaller than what we get out of the BAR averages.

       For the comparison we can use twice the tolerance. */
    while (DG2 - DG0 > 2 * tol)
    {
        DG1 = 0.5 * (DG0 + DG2);

        /* calculate the BAR averages */
        dDG1 = 0.;

        for (i = 0; i < ca->nsamples; i++)
        {
            samples_t*      s = ca->s[i];
            sample_range_t* r = &(ca->r[i]);
            if (r->use)
            {
                if (s->hist)
                {
                    dDG1 += calc_bar_sum_hist(s->hist, Wfac1, (M - DG1), type);
                }
                else
                {
                    dDG1 += calc_bar_sum(r->end - r->start, s->du + r->start, Wfac1, (M - DG1));
                }
            }
        }
        for (i = 0; i < cb->nsamples; i++)
        {
            samples_t*      s = cb->s[i];
            sample_range_t* r = &(cb->r[i]);
            if (r->use)
            {
                if (s->hist)
                {
                    dDG1 -= calc_bar_sum_hist(s->hist, Wfac2, -(M - DG1), type);
                }
                else
                {
                    dDG1 -= calc_bar_sum(r->end - r->start, s->du + r->start, Wfac2, -(M - DG1));
                }
            }
        }

        if (dDG1 < 0)
        {
            DG0 = DG1;
        }
        else
        {
            DG2 = DG1;
        }
        if (debug)
        {
            fprintf(debug, "DG %9.5f %9.5f\n", DG0, DG2);
        }
    }

    return 0.5 * (DG0 + DG2);
}

static void calc_rel_entropy(sample_coll_t* ca, sample_coll_t* cb, double temp, double dg, double* sa, double* sb)
{
    int    i, j;
    double W_ab = 0.;
    double W_ba = 0.;
    double kT, beta;
    double Wfac1, Wfac2;
    double n1, n2;

    kT   = gmx::c_boltz * temp;
    beta = 1 / kT;

    /* count the numbers of samples */
    n1 = ca->ntot;
    n2 = cb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    /*if (!lambda_vec_same(ca->native_lambda, ca->foreign_lambda))*/
    if (ca->foreign_lambda->dhdl < 0)
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a
           multiplication factor.  */
        double delta_lambda = lambda_vec_abs_diff(cb->native_lambda, ca->native_lambda);
        Wfac1               = beta * delta_lambda;
        Wfac2               = -beta * delta_lambda;
    }

    /* first calculate the average work in both directions */
    for (i = 0; i < ca->nsamples; i++)
    {
        samples_t*      s = ca->s[i];
        sample_range_t* r = &(ca->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for (j = r->start; j < r->end; j++)
                {
                    W_ab += Wfac1 * s->du[j];
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int    hd = 0; /* histogram direction */
                if ((s->hist->nhist > 1) && (Wfac1 < 0))
                {
                    hd = 1;
                }
                dx = s->hist->dx[hd];

                for (j = 0; j < s->hist->nbin[0]; j++)
                {
                    double x    = Wfac1 * ((j + s->hist->x0[0]) + 0.5) * dx; /*bin ctr*/
                    double pxdx = s->hist->bin[0][j] * normdx;               /* p(x)dx */
                    W_ab += pxdx * x;
                }
            }
        }
    }
    W_ab /= n1;

    for (i = 0; i < cb->nsamples; i++)
    {
        samples_t*      s = cb->s[i];
        sample_range_t* r = &(cb->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for (j = r->start; j < r->end; j++)
                {
                    W_ba += Wfac1 * s->du[j];
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int    hd = 0; /* histogram direction */
                if ((s->hist->nhist > 1) && (Wfac2 < 0))
                {
                    hd = 1;
                }
                dx = s->hist->dx[hd];

                for (j = 0; j < s->hist->nbin[0]; j++)
                {
                    double x    = Wfac1 * ((j + s->hist->x0[0]) + 0.5) * dx; /*bin ctr*/
                    double pxdx = s->hist->bin[0][j] * normdx;               /* p(x)dx */
                    W_ba += pxdx * x;
                }
            }
        }
    }
    W_ba /= n2;

    /* then calculate the relative entropies */
    *sa = (W_ab - dg);
    *sb = (W_ba + dg);
}

static void calc_dg_stddev(sample_coll_t* ca, sample_coll_t* cb, double temp, double dg, double* stddev)
{
    int    i, j;
    double M;
    double sigmafact = 0.;
    double kT, beta;
    double Wfac1, Wfac2;
    double n1, n2;

    kT   = gmx::c_boltz * temp;
    beta = 1 / kT;

    /* count the numbers of samples */
    n1 = ca->ntot;
    n2 = cb->ntot;

    /* to ensure the work values are the same as during the delta_G */
    /*if (!lambda_vec_same(ca->native_lambda, ca->foreign_lambda))*/
    if (ca->foreign_lambda->dhdl < 0)
    {
        /* this is the case when the delta U were calculated directly
           (i.e. we're not scaling dhdl) */
        Wfac1 = beta;
        Wfac2 = beta;
    }
    else
    {
        /* we're using dhdl, so delta_lambda needs to be a
           multiplication factor.  */
        double delta_lambda = lambda_vec_abs_diff(cb->native_lambda, ca->native_lambda);
        Wfac1               = beta * delta_lambda;
        Wfac2               = -beta * delta_lambda;
    }

    M = std::log(n1 / n2);


    /* calculate average in both directions */
    for (i = 0; i < ca->nsamples; i++)
    {
        samples_t*      s = ca->s[i];
        sample_range_t* r = &(ca->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for (j = r->start; j < r->end; j++)
                {
                    sigmafact += 1. / (2. + 2. * std::cosh((M + Wfac1 * s->du[j] - dg)));
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int    hd = 0; /* histogram direction */
                if ((s->hist->nhist > 1) && (Wfac1 < 0))
                {
                    hd = 1;
                }
                dx = s->hist->dx[hd];

                for (j = 0; j < s->hist->nbin[0]; j++)
                {
                    double x    = Wfac1 * ((j + s->hist->x0[0]) + 0.5) * dx; /*bin ctr*/
                    double pxdx = s->hist->bin[0][j] * normdx;               /* p(x)dx */

                    sigmafact += pxdx / (2. + 2. * std::cosh((M + x - dg)));
                }
            }
        }
    }
    for (i = 0; i < cb->nsamples; i++)
    {
        samples_t*      s = cb->s[i];
        sample_range_t* r = &(cb->r[i]);
        if (r->use)
        {
            if (!s->hist)
            {
                for (j = r->start; j < r->end; j++)
                {
                    sigmafact += 1. / (2. + 2. * std::cosh((M - Wfac2 * s->du[j] - dg)));
                }
            }
            else
            {
                /* normalization factor multiplied with bin width and
                   number of samples (we normalize through M): */
                double normdx = 1.;
                double dx;
                int    hd = 0; /* histogram direction */
                if ((s->hist->nhist > 1) && (Wfac2 < 0))
                {
                    hd = 1;
                }
                dx = s->hist->dx[hd];

                for (j = 0; j < s->hist->nbin[0]; j++)
                {
                    double x    = Wfac2 * ((j + s->hist->x0[0]) + 0.5) * dx; /*bin ctr*/
                    double pxdx = s->hist->bin[0][j] * normdx;               /* p(x)dx */

                    sigmafact += pxdx / (2. + 2. * std::cosh((M - x - dg)));
                }
            }
        }
    }

    sigmafact /= (n1 + n2);


    /* Eq. 10 from
       Shirts, Bair, Hooker & Pande, Phys. Rev. Lett 91, 140601 (2003): */
    *stddev = std::sqrt(((1.0 / sigmafact) - ((n1 + n2) / n1 + (n1 + n2) / n2)));
}


static void calc_bar(barres_t* br, double tol, int npee_min, int npee_max, gmx_bool* bEE, double* partsum)
{
    int    npee, p;
    double dg_sig2, sa_sig2, sb_sig2, stddev_sig2; /* intermediate variance values
                                                      for calculated quantities */
    double   temp = br->a->temp;
    int      i;
    double   dg_min, dg_max;
    gmx_bool have_hist = FALSE;

    br->dg = calc_bar_lowlevel(br->a, br->b, temp, tol, 0);

    br->dg_disc_err      = 0.;
    br->dg_histrange_err = 0.;

    /* check if there are histograms */
    for (i = 0; i < br->a->nsamples; i++)
    {
        if (br->a->r[i].use && br->a->s[i]->hist)
        {
            have_hist = TRUE;
            break;
        }
    }
    if (!have_hist)
    {
        for (i = 0; i < br->b->nsamples; i++)
        {
            if (br->b->r[i].use && br->b->s[i]->hist)
            {
                have_hist = TRUE;
                break;
            }
        }
    }

    /* calculate histogram-specific errors */
    if (have_hist)
    {
        dg_min = calc_bar_lowlevel(br->a, br->b, temp, tol, -1);
        dg_max = calc_bar_lowlevel(br->a, br->b, temp, tol, 1);

        if (std::abs(dg_max - dg_min) > GMX_REAL_EPS * 10)
        {
            /* the histogram range  error is the biggest of the differences
               between the best estimate and the extremes */
            br->dg_histrange_err = std::abs(dg_max - dg_min);
        }
        br->dg_disc_err = 0.;
        for (i = 0; i < br->a->nsamples; i++)
        {
            if (br->a->s[i]->hist)
            {
                br->dg_disc_err = std::max(br->dg_disc_err, br->a->s[i]->hist->dx[0]);
            }
        }
        for (i = 0; i < br->b->nsamples; i++)
        {
            if (br->b->s[i]->hist)
            {
                br->dg_disc_err = std::max(br->dg_disc_err, br->b->s[i]->hist->dx[0]);
            }
        }
    }
    calc_rel_entropy(br->a, br->b, temp, br->dg, &(br->sa), &(br->sb));

    calc_dg_stddev(br->a, br->b, temp, br->dg, &(br->dg_stddev));

    dg_sig2     = 0;
    sa_sig2     = 0;
    sb_sig2     = 0;
    stddev_sig2 = 0;

    *bEE = TRUE;
    {
        sample_coll_t ca, cb;

        /* initialize the samples */
        sample_coll_init(&ca, br->a->native_lambda, br->a->foreign_lambda, br->a->temp);
        sample_coll_init(&cb, br->b->native_lambda, br->b->foreign_lambda, br->b->temp);

        for (npee = npee_min; npee <= npee_max; npee++)
        {
            double dgs      = 0;
            double dgs2     = 0;
            double dsa      = 0;
            double dsb      = 0;
            double dsa2     = 0;
            double dsb2     = 0;
            double dstddev  = 0;
            double dstddev2 = 0;


            for (p = 0; p < npee; p++)
            {
                double   dgp;
                double   stddevc;
                double   sac, sbc;
                gmx_bool cac, cbc;

                cac = sample_coll_create_subsample(&ca, br->a, p, npee);
                cbc = sample_coll_create_subsample(&cb, br->b, p, npee);

                if (!cac || !cbc)
                {
                    printf("WARNING: histogram number incompatible with block number for "
                           "averaging: can't do error estimate\n");
                    *bEE = FALSE;
                    if (cac)
                    {
                        sample_coll_destroy(&ca);
                    }
                    if (cbc)
                    {
                        sample_coll_destroy(&cb);
                    }
                    return;
                }

                dgp = calc_bar_lowlevel(&ca, &cb, temp, tol, 0);
                dgs += dgp;
                dgs2 += dgp * dgp;

                partsum[npee * (npee_max + 1) + p] += dgp;

                calc_rel_entropy(&ca, &cb, temp, dgp, &sac, &sbc);
                dsa += sac;
                dsa2 += sac * sac;
                dsb += sbc;
                dsb2 += sbc * sbc;
                calc_dg_stddev(&ca, &cb, temp, dgp, &stddevc);

                dstddev += stddevc;
                dstddev2 += stddevc * stddevc;

                sample_coll_destroy(&ca);
                sample_coll_destroy(&cb);
            }
            dgs /= npee;
            dgs2 /= npee;
            dg_sig2 += (dgs2 - dgs * dgs) / (npee - 1);

            dsa /= npee;
            dsa2 /= npee;
            dsb /= npee;
            dsb2 /= npee;
            sa_sig2 += (dsa2 - dsa * dsa) / (npee - 1);
            sb_sig2 += (dsb2 - dsb * dsb) / (npee - 1);

            dstddev /= npee;
            dstddev2 /= npee;
            stddev_sig2 += (dstddev2 - dstddev * dstddev) / (npee - 1);
        }
        br->dg_err        = std::sqrt(dg_sig2 / (npee_max - npee_min + 1));
        br->sa_err        = std::sqrt(sa_sig2 / (npee_max - npee_min + 1));
        br->sb_err        = std::sqrt(sb_sig2 / (npee_max - npee_min + 1));
        br->dg_stddev_err = std::sqrt(stddev_sig2 / (npee_max - npee_min + 1));
    }
}


static double bar_err(int nbmin, int nbmax, const double* partsum)
{
    int    nb, b;
    double svar, s, s2, dg;

    svar = 0;
    for (nb = nbmin; nb <= nbmax; nb++)
    {
        s  = 0;
        s2 = 0;
        for (b = 0; b < nb; b++)
        {
            dg = partsum[nb * (nbmax + 1) + b];
            s += dg;
            s2 += dg * dg;
        }
        s /= nb;
        s2 /= nb;
        svar += (s2 - s * s) / (nb - 1);
    }

    return std::sqrt(svar / (nbmax + 1 - nbmin));
}


/* Seek the end of an identifier (consecutive non-spaces), followed by
   an optional number of spaces or '='-signs. Returns a pointer to the
   first non-space value found after that. Returns NULL if the string
   ends before that.
 */
static const char* find_value(const char* str)
{
    gmx_bool name_end_found = FALSE;

    /* if the string is a NULL pointer, return a NULL pointer. */
    if (str == nullptr)
    {
        return nullptr;
    }
    while (*str != '\0')
    {
        /* first find the end of the name */
        if (!name_end_found)
        {
            if (std::isspace(*str) || (*str == '='))
            {
                name_end_found = TRUE;
            }
        }
        else
        {
            if (!(std::isspace(*str) || (*str == '=')))
            {
                return str;
            }
        }
        str++;
    }
    return nullptr;
}


/* read a vector-notation description of a lambda vector */
static gmx_bool read_lambda_compvec(const char*                str,
                                    lambda_vec_t*              lv,
                                    const lambda_components_t* lc_in,
                                    lambda_components_t*       lc_out,
                                    const char**               end,
                                    const char*                fn)
{
    gmx_bool initialize_lc = FALSE;  /* whether to initialize the lambda
                                        components, or to check them */
    gmx_bool start_reached = FALSE;  /* whether the start of component names
                                        has been reached */
    gmx_bool    vector    = FALSE;   /* whether there are multiple components */
    int         n         = 0;       /* current component number */
    const char* val_start = nullptr; /* start of the component name, or NULL
                                        if not in a value */
    char*    strtod_end;
    gmx_bool OK = TRUE;

    if (end)
    {
        *end = str;
    }


    if (lc_out && lc_out->N == 0)
    {
        initialize_lc = TRUE;
    }

    if (lc_in == nullptr)
    {
        lc_in = lc_out;
    }

    while (true)
    {
        if (!start_reached)
        {
            if (std::isalnum(*str))
            {
                vector        = FALSE;
                start_reached = TRUE;
                val_start     = str;
            }
            else if (*str == '(')
            {
                vector        = TRUE;
                start_reached = TRUE;
            }
            else if (!std::isspace(*str))
            {
                gmx_fatal(FARGS, "Error in lambda components in %s", fn);
            }
        }
        else
        {
            if (val_start)
            {
                if (std::isspace(*str) || *str == ')' || *str == ',' || *str == '\0')
                {
                    /* end of value */
                    if (lv == nullptr)
                    {
                        if (initialize_lc)
                        {
                            lambda_components_add(lc_out, val_start, (str - val_start));
                        }
                        else
                        {
                            if (!lambda_components_check(lc_out, n, val_start, (str - val_start)))
                            {
                                return FALSE;
                            }
                        }
                    }
                    else
                    {
                        /* add a vector component to lv */
                        lv->val[n] = strtod(val_start, &strtod_end);
                        if (val_start == strtod_end)
                        {
                            gmx_fatal(FARGS, "Error reading lambda vector in %s", fn);
                        }
                    }
                    /* reset for the next identifier */
                    val_start = nullptr;
                    n++;
                    if (!vector)
                    {
                        return OK;
                    }
                }
            }
            else if (std::isalnum(*str))
            {
                val_start = str;
            }
            if (*str == ')')
            {
                str++;
                if (end)
                {
                    *end = str;
                }
                if (!vector)
                {
                    gmx_fatal(FARGS, "Error in lambda components in %s", fn);
                }
                else
                {
                    GMX_RELEASE_ASSERT(lc_in != nullptr, "Internal inconsistency? lc_in==NULL");
                    if (n == lc_in->N)
                    {
                        return OK;
                    }
                    else if (lv == nullptr)
                    {
                        return FALSE;
                    }
                    else
                    {
                        gmx_fatal(FARGS, "Incomplete lambda vector data in %s", fn);
                        return FALSE;
                    }
                }
            }
        }
        if (*str == '\0')
        {
            break;
        }
        str++;
        if (end)
        {
            *end = str;
        }
    }
    if (vector)
    {
        gmx_fatal(FARGS, "Incomplete lambda components data in %s", fn);
        return FALSE;
    }
    return OK;
}

/* read and check the component names from a string */
static gmx_bool read_lambda_components(const char* str, lambda_components_t* lc, const char** end, const char* fn)
{
    return read_lambda_compvec(str, nullptr, nullptr, lc, end, fn);
}

/* read an initialized lambda vector from a string */
static gmx_bool read_lambda_vector(const char* str, lambda_vec_t* lv, const char** end, const char* fn)
{
    return read_lambda_compvec(str, lv, lv->lc, nullptr, end, fn);
}


/* deduce lambda value from legend.
    fn = the file name
    legend = the legend string
    ba = the xvg data
    lam = the initialized lambda vector
   returns whether to use the data in this set.
 */
static gmx_bool legend2lambda(const char* fn, const char* legend, lambda_vec_t* lam)
{
    const char *ptr = nullptr, *ptr2 = nullptr;
    gmx_bool    ok    = FALSE;
    gmx_bool    bdhdl = FALSE;
    const char* tostr = " to ";

    if (legend == nullptr)
    {
        gmx_fatal(FARGS, "There is no legend in file '%s', can not deduce lambda", fn);
    }

    /* look for the last 'to': */
    ptr2 = legend;
    do
    {
        ptr2 = std::strstr(ptr2, tostr);
        if (ptr2 != nullptr)
        {
            ptr = ptr2;
            ptr2++;
        }
    } while (ptr2 != nullptr && *ptr2 != '\0');

    if (ptr)
    {
        ptr += std::strlen(tostr) - 1; /* and advance past that 'to' */
    }
    else
    {
        /* look for the = sign */
        ptr = std::strrchr(legend, '=');
        if (!ptr)
        {
            /* otherwise look for the last space */
            ptr = std::strrchr(legend, ' ');
        }
    }

    if (std::strstr(legend, "dH"))
    {
        ok    = TRUE;
        bdhdl = TRUE;
    }
    else if (std::strchr(legend, 'D') != nullptr && std::strchr(legend, 'H') != nullptr)
    {
        ok    = TRUE;
        bdhdl = FALSE;
    }
    else /*if (std::strstr(legend, "pV"))*/
    {
        return FALSE;
    }
    if (!ptr)
    {
        ok = FALSE;
    }

    if (!ok)
    {
        gmx_fatal(FARGS, "There is no proper lambda legend in file '%s', can not deduce lambda", fn);
    }
    if (!bdhdl)
    {
        ptr = find_value(ptr);
        if (!ptr || !read_lambda_vector(ptr, lam, nullptr, fn))
        {
            gmx_fatal(FARGS, "lambda vector '%s' %s faulty", legend, fn);
        }
    }
    else
    {
        int         dhdl_index;
        const char* end;

        ptr = std::strrchr(legend, '=');
        end = ptr;
        if (ptr)
        {
            /* there must be a component name */
            ptr--;
            if (ptr < legend)
            {
                gmx_fatal(FARGS, "dhdl legend '%s' %s faulty", legend, fn);
            }
            /* now backtrack to the start of the identifier */
            while (isspace(*ptr))
            {
                end = ptr;
                ptr--;
                if (ptr < legend)
                {
                    gmx_fatal(FARGS, "dhdl legend '%s' %s faulty", legend, fn);
                }
            }
            while (!std::isspace(*ptr))
            {
                ptr--;
                if (ptr < legend)
                {
                    gmx_fatal(FARGS, "dhdl legend '%s' %s faulty", legend, fn);
                }
            }
            ptr++;
            dhdl_index = lambda_components_find(lam->lc, ptr, (end - ptr));
            if (dhdl_index < 0)
            {
                char buf[STRLEN];
                std::strncpy(buf, ptr, (end - ptr));
                buf[(end - ptr)] = '\0';
                gmx_fatal(FARGS, "Did not find lambda component for '%s' in %s", buf, fn);
            }
        }
        else
        {
            if (lam->lc->N > 1)
            {
                gmx_fatal(FARGS, "dhdl without component name with >1 lambda component in %s", fn);
            }
            dhdl_index = 0;
        }
        lam->dhdl = dhdl_index;
    }
    return TRUE;
}

static gmx_bool subtitle2lambda(const char* subtitle, xvg_t* ba, const char* fn, lambda_components_t* lc)
{
    gmx_bool    bFound;
    const char* ptr;
    char*       end;
    double      native_lambda;

    bFound = FALSE;

    /* first check for a state string */
    ptr = std::strstr(subtitle, "state");
    if (ptr)
    {
        int         index = -1;
        const char* val_end;

        /* the new 4.6 style lambda vectors */
        ptr = find_value(ptr);
        if (ptr)
        {
            index = std::strtol(ptr, &end, 10);
            if (ptr == end)
            {
                gmx_fatal(FARGS, "Incomplete state data in %s", fn);
                return FALSE;
            }
            ptr = end;
        }
        else
        {
            gmx_fatal(FARGS, "Incomplete state data in %s", fn);
            return FALSE;
        }
        /* now find the lambda vector component names */
        while (*ptr != '(' && !std::isalnum(*ptr))
        {
            ptr++;
            if (*ptr == '\0')
            {
                gmx_fatal(FARGS, "Incomplete lambda vector component data in %s", fn);
                return FALSE;
            }
        }
        val_end = ptr;
        if (!read_lambda_components(ptr, lc, &val_end, fn))
        {
            gmx_fatal(FARGS, "lambda vector components in %s don't match those previously read", fn);
        }
        ptr = find_value(val_end);
        if (!ptr)
        {
            gmx_fatal(FARGS, "Incomplete state data in %s", fn);
            return FALSE;
        }
        lambda_vec_init(&(ba->native_lambda), lc);
        if (!read_lambda_vector(ptr, &(ba->native_lambda), nullptr, fn))
        {
            gmx_fatal(FARGS, "lambda vector in %s faulty", fn);
        }
        ba->native_lambda.index = index;
        bFound                  = TRUE;
    }
    else
    {
        /* compatibility mode: check for lambda in other ways. */
        /* plain text lambda string */
        ptr = std::strstr(subtitle, "lambda");
        if (ptr == nullptr)
        {
            /* xmgrace formatted lambda string */
            ptr = std::strstr(subtitle, "\\xl\\f{}");
        }
        if (ptr == nullptr)
        {
            /* xmgr formatted lambda string */
            ptr = std::strstr(subtitle, "\\8l\\4");
        }
        if (ptr != nullptr)
        {
            ptr = std::strstr(ptr, "=");
        }
        if (ptr != nullptr)
        {
            bFound = (sscanf(ptr + 1, "%lf", &(native_lambda)) == 1);
            /* add the lambda component name as an empty string */
            if (lc->N > 0)
            {
                if (!lambda_components_check(lc, 0, "", 0))
                {
                    gmx_fatal(FARGS,
                              "lambda vector components in %s don't match those previously read",
                              fn);
                }
            }
            else
            {
                lambda_components_add(lc, "", 0);
            }
            lambda_vec_init(&(ba->native_lambda), lc);
            ba->native_lambda.val[0] = native_lambda;
        }
    }

    return bFound;
}

static void read_bar_xvg_lowlevel(const char* fn, const real* temp, xvg_t* ba, lambda_components_t* lc)
{
    int      i;
    char *   subtitle, **legend, *ptr;
    int      np;
    gmx_bool native_lambda_read = FALSE;
    char     buf[STRLEN];

    xvg_init(ba);

    ba->filename = fn;

    np = read_xvg_legend(fn, &ba->y, &ba->nset, &subtitle, &legend);
    if (!ba->y)
    {
        gmx_fatal(FARGS, "File %s contains no usable data.", fn);
    }
    /* Reorder the data */
    ba->t = ba->y[0];
    for (i = 1; i < ba->nset; i++)
    {
        ba->y[i - 1] = ba->y[i];
    }
    ba->nset--;

    snew(ba->np, ba->nset);
    for (i = 0; i < ba->nset; i++)
    {
        ba->np[i] = np;
    }

    ba->temp = -1;
    if (subtitle != nullptr)
    {
        /* try to extract temperature */
        ptr = std::strstr(subtitle, "T =");
        if (ptr != nullptr)
        {
            ptr += 3;
            if (sscanf(ptr, "%lf", &ba->temp) == 1)
            {
                if (ba->temp <= 0)
                {
                    gmx_fatal(FARGS, "Found temperature of %f in file '%s'", ba->temp, fn);
                }
            }
        }
    }
    if (ba->temp < 0)
    {
        if (*temp <= 0)
        {
            gmx_fatal(FARGS,
                      "Did not find a temperature in the subtitle in file '%s', use the -temp "
                      "option of [TT]gmx bar[tt]",
                      fn);
        }
        ba->temp = *temp;
    }

    /* Try to deduce lambda from the subtitle */
    if (subtitle)
    {
        if (subtitle2lambda(subtitle, ba, fn, lc))
        {
            native_lambda_read = TRUE;
        }
    }

    if (!native_lambda_read)
    {
        gmx_fatal(FARGS, "File %s contains multiple sets but no indication of the native lambda", fn);
    }

    snew(ba->lambda, ba->nset);
    if (legend == nullptr)
    {
        /* Check if we have a single set, no legend, nset=1 means t and dH/dl */
        if (ba->nset == 1)
        {
            ba->lambda[0] = ba->native_lambda;
        }
        else
        {
            gmx_fatal(FARGS,
                      "File %s contains multiple sets but no legends, can not determine the lambda "
                      "values",
                      fn);
        }
    }
    else
    {
        for (i = 0; i < ba->nset;)
        {
            /* Read lambda from the legend */
            lambda_vec_init(&(ba->lambda[i]), lc);
            lambda_vec_copy(&(ba->lambda[i]), &(ba->native_lambda));
            gmx_bool use = legend2lambda(fn, legend[i], &(ba->lambda[i]));
            if (use)
            {
                lambda_vec_print(&(ba->lambda[i]), buf, FALSE);
                i++;
            }
            else
            {
                int j;
                printf("%s: Ignoring set '%s'.\n", fn, legend[i]);
                for (j = i + 1; j < ba->nset; j++)
                {
                    ba->y[j - 1]  = ba->y[j];
                    legend[j - 1] = legend[j];
                }
                ba->nset--;
            }
        }
    }

    if (legend != nullptr)
    {
        for (i = 0; i < ba->nset - 1; i++)
        {
            sfree(legend[i]);
        }
        sfree(legend);
    }
}

static void read_bar_xvg(const char* fn, real* temp, sim_data_t* sd)
{
    xvg_t*     barsim;
    samples_t* s;
    int        i;

    snew(barsim, 1);

    read_bar_xvg_lowlevel(fn, temp, barsim, &(sd->lc));

    if (barsim->nset < 1)
    {
        gmx_fatal(FARGS, "File '%s' contains fewer than two columns", fn);
    }

    if (!gmx_within_tol(*temp, barsim->temp, GMX_FLOAT_EPS) && (*temp > 0))
    {
        gmx_fatal(FARGS, "Temperature in file %s different from earlier files or setting\n", fn);
    }
    *temp = barsim->temp;

    /* now create a series of samples_t */
    snew(s, barsim->nset);
    for (i = 0; i < barsim->nset; i++)
    {
        samples_init(s + i,
                     &(barsim->native_lambda),
                     &(barsim->lambda[i]),
                     barsim->temp,
                     lambda_vec_same(&(barsim->native_lambda), &(barsim->lambda[i])),
                     fn);
        s[i].du  = barsim->y[i];
        s[i].ndu = barsim->np[i];
        s[i].t   = barsim->t;

        lambda_data_list_insert_sample(sd->lb, s + i);
    }
    {
        char buf[STRLEN];

        lambda_vec_print(s[0].native_lambda, buf, FALSE);
        printf("%s: %.1f - %.1f; lambda = %s\n    dH/dl & foreign lambdas:\n",
               fn,
               s[0].t[0],
               s[0].t[s[0].ndu - 1],
               buf);
        for (i = 0; i < barsim->nset; i++)
        {
            lambda_vec_print(s[i].foreign_lambda, buf, TRUE);
            printf("        %s (%d pts)\n", buf, s[i].ndu);
        }
    }
    printf("\n\n");
}

static void read_edr_rawdh_block(samples_t**   smp,
                                 int*          ndu,
                                 t_enxblock*   blk,
                                 double        start_time,
                                 double        delta_time,
                                 lambda_vec_t* native_lambda,
                                 double        temp,
                                 double*       last_t,
                                 const char*   filename)
{
    int           i, j;
    lambda_vec_t* foreign_lambda;
    int           type;
    samples_t*    s; /* convenience pointer */
    int           startj;

    /* check the block types etc. */
    if ((blk->nsub < 3) || (blk->sub[0].type != XdrDataType::Int)
        || (blk->sub[1].type != XdrDataType::Double)
        || ((blk->sub[2].type != XdrDataType::Float) && (blk->sub[2].type != XdrDataType::Double))
        || (blk->sub[0].nr < 1) || (blk->sub[1].nr < 1))
    {
        gmx_fatal(FARGS, "Unexpected/corrupted block data in file %s around time %f.", filename, start_time);
    }

    snew(foreign_lambda, 1);
    lambda_vec_init(foreign_lambda, native_lambda->lc);
    lambda_vec_copy(foreign_lambda, native_lambda);
    type = blk->sub[0].ival[0];
    if (type == dhbtDH)
    {
        for (i = 0; i < native_lambda->lc->N; i++)
        {
            foreign_lambda->val[i] = blk->sub[1].dval[i];
        }
    }
    else
    {
        if (blk->sub[0].nr > 1)
        {
            foreign_lambda->dhdl = blk->sub[0].ival[1];
        }
        else
        {
            foreign_lambda->dhdl = 0;
        }
    }

    if (!*smp)
    {
        /* initialize the samples structure if it's empty. */
        snew(*smp, 1);
        samples_init(*smp, native_lambda, foreign_lambda, temp, type == dhbtDHDL, filename);
        (*smp)->start_time = start_time;
        (*smp)->delta_time = delta_time;
    }

    /* set convenience pointer */
    s = *smp;

    /* now double check */
    if (!lambda_vec_same(s->foreign_lambda, foreign_lambda))
    {
        char buf[STRLEN], buf2[STRLEN];
        lambda_vec_print(foreign_lambda, buf, FALSE);
        lambda_vec_print(s->foreign_lambda, buf2, FALSE);
        fprintf(stderr, "Got foreign lambda=%s, expected: %s\n", buf, buf2);
        gmx_fatal(FARGS, "Corrupted data in file %s around t=%f.", filename, start_time);
    }

    /* make room for the data */
    if (gmx::Index(s->ndu_alloc) < s->ndu + blk->sub[2].nr)
    {
        s->ndu_alloc += (s->ndu_alloc < static_cast<size_t>(blk->sub[2].nr)) ? blk->sub[2].nr * 2
                                                                             : s->ndu_alloc;
        srenew(s->du_alloc, s->ndu_alloc);
        s->du = s->du_alloc;
    }
    startj = s->ndu;
    s->ndu += blk->sub[2].nr;
    s->ntot += blk->sub[2].nr;
    *ndu = blk->sub[2].nr;

    /* and copy the data*/
    for (j = 0; j < blk->sub[2].nr; j++)
    {
        if (blk->sub[2].type == XdrDataType::Float)
        {
            s->du[startj + j] = blk->sub[2].fval[j];
        }
        else
        {
            s->du[startj + j] = blk->sub[2].dval[j];
        }
    }
    if (start_time + blk->sub[2].nr * delta_time > *last_t)
    {
        *last_t = start_time + blk->sub[2].nr * delta_time;
    }
}

static samples_t* read_edr_hist_block(int*          nsamples,
                                      t_enxblock*   blk,
                                      double        start_time,
                                      double        delta_time,
                                      lambda_vec_t* native_lambda,
                                      double        temp,
                                      double*       last_t,
                                      const char*   filename)
{
    int           i, j;
    samples_t*    s;
    int           nhist;
    lambda_vec_t* foreign_lambda;
    int           type;
    int           nbins[2];

    /* check the block types etc. */
    if ((blk->nsub < 2) || (blk->sub[0].type != XdrDataType::Double)
        || (blk->sub[1].type != XdrDataType::Int64) || (blk->sub[0].nr < 2) || (blk->sub[1].nr < 2))
    {
        gmx_fatal(FARGS, "Unexpected/corrupted block data in file %s around time %f", filename, start_time);
    }

    nhist = blk->nsub - 2;
    if (nhist == 0)
    {
        return nullptr;
    }
    if (nhist > 2)
    {
        gmx_fatal(FARGS, "Unexpected/corrupted block data in file %s around time %f", filename, start_time);
    }

    snew(s, 1);
    *nsamples = 1;

    snew(foreign_lambda, 1);
    lambda_vec_init(foreign_lambda, native_lambda->lc);
    lambda_vec_copy(foreign_lambda, native_lambda);
    type = static_cast<int>(blk->sub[1].lval[1]);
    if (type == dhbtDH)
    {
        double old_foreign_lambda;

        old_foreign_lambda = blk->sub[0].dval[0];
        if (old_foreign_lambda >= 0)
        {
            foreign_lambda->val[0] = old_foreign_lambda;
            if (foreign_lambda->lc->N > 1)
            {
                gmx_fatal(FARGS, "Single-component lambda in multi-component file %s", filename);
            }
        }
        else
        {
            for (i = 0; i < native_lambda->lc->N; i++)
            {
                foreign_lambda->val[i] = blk->sub[0].dval[i + 2];
            }
        }
    }
    else
    {
        if (foreign_lambda->lc->N > 1)
        {
            if (blk->sub[1].nr < 3 + nhist)
            {
                gmx_fatal(FARGS, "Missing derivative coord in multi-component file %s", filename);
            }
            foreign_lambda->dhdl = blk->sub[1].lval[2 + nhist];
        }
        else
        {
            foreign_lambda->dhdl = 0;
        }
    }

    samples_init(s, native_lambda, foreign_lambda, temp, type == dhbtDHDL, filename);
    snew(s->hist, 1);

    for (i = 0; i < nhist; i++)
    {
        nbins[i] = blk->sub[i + 2].nr;
    }

    hist_init(s->hist, nhist, nbins);

    for (i = 0; i < nhist; i++)
    {
        s->hist->x0[i] = blk->sub[1].lval[2 + i];
        s->hist->dx[i] = blk->sub[0].dval[1];
        if (i == 1)
        {
            s->hist->dx[i] = -s->hist->dx[i];
        }
    }

    s->hist->start_time = start_time;
    s->hist->delta_time = delta_time;
    s->start_time       = start_time;
    s->delta_time       = delta_time;

    for (i = 0; i < nhist; i++)
    {
        int64_t sum = 0;

        for (j = 0; j < s->hist->nbin[i]; j++)
        {
            int binv = static_cast<int>(blk->sub[i + 2].ival[j]);

            s->hist->bin[i][j] = binv;
            sum += binv;
        }
        if (i == 0)
        {
            s->ntot      = sum;
            s->hist->sum = sum;
        }
        else
        {
            if (s->ntot != sum)
            {
                gmx_fatal(FARGS, "Histogram counts don't match in %s", filename);
            }
        }
    }

    if (start_time + s->hist->sum * delta_time > *last_t)
    {
        *last_t = start_time + s->hist->sum * delta_time;
    }
    return s;
}


static void read_barsim_edr(const char* fn, real* temp, sim_data_t* sd)
{
    int            i, j;
    ener_file_t    fp;
    t_enxframe*    fr;
    int            nre;
    gmx_enxnm_t*   enm           = nullptr;
    double         first_t       = -1;
    double         last_t        = -1;
    samples_t**    samples_rawdh = nullptr; /* contains samples for raw delta_h  */
    int*           nhists        = nullptr; /* array to keep count & print at end */
    int*           npts          = nullptr; /* array to keep count & print at end */
    lambda_vec_t** lambdas       = nullptr; /* array to keep count & print at end */
    lambda_vec_t*  native_lambda;
    int            nsamples = 0;
    lambda_vec_t   start_lambda;

    fp = open_enx(fn, "r");
    do_enxnms(fp, &nre, &enm);
    snew(fr, 1);

    snew(native_lambda, 1);
    start_lambda.lc  = nullptr;
    start_lambda.val = nullptr;

    while (do_enx(fp, fr))
    {
        /* count the data blocks */
        int nblocks_raw  = 0;
        int nblocks_hist = 0;
        int nlam         = 0;
        int k;
        /* DHCOLL block information: */
        double start_time = 0, delta_time = 0, old_start_lambda = 0, delta_lambda = 0;
        double rtemp = 0;

        /* count the blocks and handle collection information: */
        for (i = 0; i < fr->nblock; i++)
        {
            if (fr->block[i].id == enxDHHIST)
            {
                nblocks_hist++;
            }
            if (fr->block[i].id == enxDH)
            {
                nblocks_raw++;
            }
            if (fr->block[i].id == enxDHCOLL)
            {
                nlam++;
                if ((fr->block[i].nsub < 1) || (fr->block[i].sub[0].type != XdrDataType::Double)
                    || (fr->block[i].sub[0].nr < 5))
                {
                    gmx_fatal(FARGS, "Unexpected block data in file %s", fn);
                }

                /* read the data from the DHCOLL block */
                rtemp            = fr->block[i].sub[0].dval[0];
                start_time       = fr->block[i].sub[0].dval[1];
                delta_time       = fr->block[i].sub[0].dval[2];
                old_start_lambda = fr->block[i].sub[0].dval[3];
                delta_lambda     = fr->block[i].sub[0].dval[4];

                if (delta_lambda != 0)
                {
                    gmx_fatal(FARGS, "Lambda values not constant in %s: can't apply BAR method", fn);
                }
                if ((*temp != rtemp) && (*temp > 0))
                {
                    gmx_fatal(FARGS,
                              "Temperature in file %s different from earlier files or setting\n",
                              fn);
                }
                *temp = rtemp;

                if (old_start_lambda >= 0)
                {
                    if (sd->lc.N > 0)
                    {
                        if (!lambda_components_check(&(sd->lc), 0, "", 0))
                        {
                            gmx_fatal(FARGS,
                                      "lambda vector components in %s don't match those previously "
                                      "read",
                                      fn);
                        }
                    }
                    else
                    {
                        lambda_components_add(&(sd->lc), "", 0);
                    }
                    if (!start_lambda.lc)
                    {
                        lambda_vec_init(&start_lambda, &(sd->lc));
                    }
                    start_lambda.val[0] = old_start_lambda;
                }
                else
                {
                    /* read lambda vector */
                    int      n_lambda_vec;
                    gmx_bool check = (sd->lc.N > 0);
                    if (fr->block[i].nsub < 2)
                    {
                        gmx_fatal(FARGS, "No lambda vector, but start_lambda=%f\n", old_start_lambda);
                    }
                    n_lambda_vec = fr->block[i].sub[1].ival[1];
                    for (j = 0; j < n_lambda_vec; j++)
                    {
                        const char* name =
                                enumValueToStringSingular(static_cast<FreeEnergyPerturbationCouplingType>(
                                        fr->block[i].sub[1].ival[1 + j]));
                        if (check)
                        {
                            /* check the components */
                            lambda_components_check(&(sd->lc), j, name, std::strlen(name));
                        }
                        else
                        {
                            lambda_components_add(&(sd->lc), name, std::strlen(name));
                        }
                    }
                    lambda_vec_init(&start_lambda, &(sd->lc));
                    start_lambda.index = fr->block[i].sub[1].ival[0];
                    for (j = 0; j < n_lambda_vec; j++)
                    {
                        start_lambda.val[j] = fr->block[i].sub[0].dval[5 + j];
                    }
                }
                if (first_t < 0)
                {
                    first_t = start_time;
                }
            }
        }

        if (nlam != 1)
        {
            gmx_fatal(FARGS, "Did not find delta H information in file %s", fn);
        }
        if (nblocks_raw > 0 && nblocks_hist > 0)
        {
            gmx_fatal(FARGS, "Can't handle both raw delta U data and histograms in the same file %s", fn);
        }

        if (nsamples == 0)
        {
            /* this is the first round; allocate the associated data
               structures */
            /*native_lambda=start_lambda;*/
            lambda_vec_init(native_lambda, &(sd->lc));
            lambda_vec_copy(native_lambda, &start_lambda);
            nsamples = nblocks_raw + nblocks_hist;
            snew(nhists, nsamples);
            snew(npts, nsamples);
            snew(lambdas, nsamples);
            snew(samples_rawdh, nsamples);
            for (i = 0; i < nsamples; i++)
            {
                nhists[i]        = 0;
                npts[i]          = 0;
                lambdas[i]       = nullptr;
                samples_rawdh[i] = nullptr; /* init to NULL so we know which
                                               ones contain values */
            }
        }
        else
        {
            // nsamples > 0 means this is NOT the first iteration

            /* check the native lambda */
            if (!lambda_vec_same(&start_lambda, native_lambda))
            {
                gmx_fatal(FARGS,
                          "Native lambda not constant in file %s: started at %f, and becomes %f at "
                          "time %f",
                          fn,
                          native_lambda->val[0],
                          start_lambda.val[0],
                          start_time);
            }
            /* check the number of samples against the previous number */
            if (((nblocks_raw + nblocks_hist) != nsamples) || (nlam != 1))
            {
                gmx_fatal(FARGS,
                          "Unexpected block count in %s: was %d, now %d\n",
                          fn,
                          nsamples + 1,
                          nblocks_raw + nblocks_hist + nlam);
            }
            /* check whether last iterations's end time matches with
               the currrent start time */
            if ((std::abs(last_t - start_time) > 2 * delta_time) && last_t >= 0)
            {
                /* it didn't. We need to store our samples and reallocate */
                for (i = 0; i < nsamples; i++)
                {
                    if (samples_rawdh[i])
                    {
                        /* insert it into the existing list */
                        lambda_data_list_insert_sample(sd->lb, samples_rawdh[i]);
                        /* and make sure we'll allocate a new one this time
                           around */
                        samples_rawdh[i] = nullptr;
                    }
                }
            }
        }

        /* and read them */
        k = 0; /* counter for the lambdas, etc. arrays */
        for (i = 0; i < fr->nblock; i++)
        {
            if (fr->block[i].id == enxDH)
            {
                int type = (fr->block[i].sub[0].ival[0]);
                if (type == dhbtDH || type == dhbtDHDL)
                {
                    int ndu;
                    read_edr_rawdh_block(&(samples_rawdh[k]),
                                         &ndu,
                                         &(fr->block[i]),
                                         start_time,
                                         delta_time,
                                         native_lambda,
                                         rtemp,
                                         &last_t,
                                         fn);
                    npts[k] += ndu;
                    if (samples_rawdh[k])
                    {
                        lambdas[k] = samples_rawdh[k]->foreign_lambda;
                    }
                    k++;
                }
            }
            else if (fr->block[i].id == enxDHHIST)
            {
                int type = static_cast<int>(fr->block[i].sub[1].lval[1]);
                if (type == dhbtDH || type == dhbtDHDL)
                {
                    int        j;
                    int        nb = 0;
                    samples_t* s; /* this is where the data will go */
                    s = read_edr_hist_block(
                            &nb, &(fr->block[i]), start_time, delta_time, native_lambda, rtemp, &last_t, fn);
                    nhists[k] += nb;
                    if (nb > 0)
                    {
                        lambdas[k] = s->foreign_lambda;
                    }
                    k++;
                    /* and insert the new sample immediately */
                    for (j = 0; j < nb; j++)
                    {
                        lambda_data_list_insert_sample(sd->lb, s + j);
                    }
                }
            }
        }
    }
    /* Now store all our extant sample collections */
    for (i = 0; i < nsamples; i++)
    {
        if (samples_rawdh[i])
        {
            /* insert it into the existing list */
            lambda_data_list_insert_sample(sd->lb, samples_rawdh[i]);
        }
    }


    {
        char buf[STRLEN];
        printf("\n");
        lambda_vec_print(native_lambda, buf, FALSE);
        printf("%s: %.1f - %.1f; lambda = %s\n    foreign lambdas:\n", fn, first_t, last_t, buf);
        for (i = 0; i < nsamples; i++)
        {
            if (lambdas[i])
            {
                lambda_vec_print(lambdas[i], buf, TRUE);
                if (nhists[i] > 0)
                {
                    printf("        %s (%d hists)\n", buf, nhists[i]);
                }
                else
                {
                    printf("        %s (%d pts)\n", buf, npts[i]);
                }
            }
        }
    }
    printf("\n\n");
    sfree(npts);
    sfree(nhists);
    sfree(lambdas);
}


int gmx_bar(int argc, char* argv[])
{
    static const char* desc[] = {
        "[THISMODULE] calculates free energy difference estimates through ",
        "Bennett's acceptance ratio method (BAR). It also automatically",
        "adds series of individual free energies obtained with BAR into",
        "a combined free energy estimate.[PAR]",

        "Every individual BAR free energy difference relies on two ",
        "simulations at different states: say state A and state B, as",
        "controlled by a parameter, [GRK]lambda[grk] (see the [REF].mdp[ref] parameter",
        "[TT]init_lambda[tt]). The BAR method calculates a ratio of weighted",
        "average of the Hamiltonian difference of state B given state A and",
        "vice versa.",
        "The energy differences to the other state must be calculated",
        "explicitly during the simulation. This can be done with",
        "the [REF].mdp[ref] option [TT]foreign_lambda[tt].[PAR]",

        "Input option [TT]-f[tt] expects multiple [TT]dhdl.xvg[tt] files. ",
        "Two types of input files are supported:",
        "",
        " * Files with more than one [IT]y[it]-value. ",
        "   The files should have columns ",
        "   with dH/d[GRK]lambda[grk] and [GRK]Delta[grk][GRK]lambda[grk]. ",
        "   The [GRK]lambda[grk] values are inferred ",
        "   from the legends: [GRK]lambda[grk] of the simulation from the legend of ",
        "   dH/d[GRK]lambda[grk] and the foreign [GRK]lambda[grk] values from the ",
        "   legends of Delta H",
        " * Files with only one [IT]y[it]-value. Using the",
        "   [TT]-extp[tt] option for these files, it is assumed",
        "   that the [IT]y[it]-value is dH/d[GRK]lambda[grk] and that the ",
        "   Hamiltonian depends linearly on [GRK]lambda[grk]. ",
        "   The [GRK]lambda[grk] value of the simulation is inferred from the ",
        "   subtitle (if present), otherwise from a number in the subdirectory ",
        "   in the file name.",
        "",

        "The [GRK]lambda[grk] of the simulation is parsed from ",
        "[TT]dhdl.xvg[tt] file's legend containing the string 'dH', the ",
        "foreign [GRK]lambda[grk] values from the legend containing the ",
        "capitalized letters 'D' and 'H'. The temperature is parsed from ",
        "the legend line containing 'T ='.[PAR]",

        "The input option [TT]-g[tt] expects multiple [REF].edr[ref] files. ",
        "These can contain either lists of energy differences (see the ",
        "[REF].mdp[ref] option [TT]separate_dhdl_file[tt]), or a series of ",
        "histograms (see the [REF].mdp[ref] options [TT]dh_hist_size[tt] and ",
        "[TT]dh_hist_spacing[tt]).",
        "The temperature and [GRK]lambda[grk] ",
        "values are automatically deduced from the [TT]ener.edr[tt] file.[PAR]",

        "In addition to the [REF].mdp[ref] option [TT]foreign_lambda[tt], ",
        "the energy difference can also be extrapolated from the ",
        "dH/d[GRK]lambda[grk] values. This is done with the[TT]-extp[tt]",
        "option, which assumes that the system's Hamiltonian depends linearly",
        "on [GRK]lambda[grk], which is not normally the case.[PAR]",

        "The free energy estimates are determined using BAR with bisection, ",
        "with the precision of the output set with [TT]-prec[tt]. ",
        "An error estimate taking into account time correlations ",
        "is made by splitting the data into blocks and determining ",
        "the free energy differences over those blocks and assuming ",
        "the blocks are independent. ",
        "The final error estimate is determined from the average variance ",
        "over 5 blocks. A range of block numbers for error estimation can ",
        "be provided with the options [TT]-nbmin[tt] and [TT]-nbmax[tt].[PAR]",

        "[THISMODULE] tries to aggregate samples with the same 'native' and ",
        "'foreign' [GRK]lambda[grk] values, but always assumes independent ",
        "samples. [BB]Note[bb] that when aggregating energy ",
        "differences/derivatives with different sampling intervals, this is ",
        "almost certainly not correct. Usually subsequent energies are ",
        "correlated and different time intervals mean different degrees ",
        "of correlation between samples.[PAR]",

        "The results are split in two parts: the last part contains the final ",
        "results in kJ/mol, together with the error estimate for each part ",
        "and the total. The first part contains detailed free energy ",
        "difference estimates and phase space overlap measures in units of ",
        "kT (together with their computed error estimate). The printed ",
        "values are:",
        "",
        " * lam_A: the [GRK]lambda[grk] values for point A.",
        " * lam_B: the [GRK]lambda[grk] values for point B.",
        " *    DG: the free energy estimate.",
        " *   s_A: an estimate of the relative entropy of B in A.",
        " *   s_B: an estimate of the relative entropy of A in B.",
        " * stdev: an estimate expected per-sample standard deviation.",
        "",

        "The relative entropy of both states in each other's ensemble can be ",
        "interpreted as a measure of phase space overlap: ",
        "the relative entropy s_A of the work samples of lambda_B in the ",
        "ensemble of lambda_A (and vice versa for s_B), is a ",
        "measure of the 'distance' between Boltzmann distributions of ",
        "the two states, that goes to zero for identical distributions. See ",
        "Wu & Kofke, J. Chem. Phys. 123 084109 (2005) for more information.",
        "[PAR]",
        "The estimate of the expected per-sample standard deviation, as given ",
        "in Bennett's original BAR paper: Bennett, J. Comp. Phys. 22, p 245 (1976).",
        "Eq. 10 therein gives an estimate of the quality of sampling (not directly",
        "of the actual statistical error, because it assumes independent samples).[PAR]",

        "To get a visual estimate of the phase space overlap, use the ",
        "[TT]-oh[tt] option to write series of histograms, together with the ",
        "[TT]-nbin[tt] option.[PAR]"
    };
    static real begin = 0, end = -1, temp = -1;
    int         nd = 2, nbmin = 5, nbmax = 5;
    int         nbin     = 100;
    gmx_bool    use_dhdl = FALSE;
    t_pargs     pa[]     = {
        { "-b", FALSE, etREAL, { &begin }, "Begin time for BAR" },
        { "-e", FALSE, etREAL, { &end }, "End time for BAR" },
        { "-temp", FALSE, etREAL, { &temp }, "Temperature (K)" },
        { "-prec", FALSE, etINT, { &nd }, "The number of digits after the decimal point" },
        { "-nbmin", FALSE, etINT, { &nbmin }, "Minimum number of blocks for error estimation" },
        { "-nbmax", FALSE, etINT, { &nbmax }, "Maximum number of blocks for error estimation" },
        { "-nbin", FALSE, etINT, { &nbin }, "Number of bins for histogram output" },
        { "-extp",
          FALSE,
          etBOOL,
          { &use_dhdl },
          "Whether to linearly extrapolate dH/dl values to use as energies" }
    };

    t_filenm fnm[] = { { efXVG, "-f", "dhdl", ffOPTRDMULT },
                       { efEDR, "-g", "ener", ffOPTRDMULT },
                       { efXVG, "-o", "bar", ffOPTWR },
                       { efXVG, "-oi", "barint", ffOPTWR },
                       { efXVG, "-oh", "histogram", ffOPTWR } };
#define NFILE asize(fnm)

    int        f;
    int        nfile_tot; /* total number of input files */
    sim_data_t sim_data;  /* the simulation data */
    barres_t*  results;   /* the results */
    int        nresults;  /* number of results in results array */

    double*           partsum;
    double            prec, dg_tot;
    FILE *            fpb, *fpi;
    char              dgformat[20], xvg2format[STRLEN], xvg3format[STRLEN];
    char              buf[STRLEN], buf2[STRLEN];
    char              ktformat[STRLEN], sktformat[STRLEN];
    char              kteformat[STRLEN], skteformat[STRLEN];
    gmx_output_env_t* oenv;
    double            kT;
    gmx_bool          result_OK = TRUE, bEE = TRUE;

    gmx_bool disc_err          = FALSE;
    double   sum_disc_err      = 0.; /* discretization error */
    gmx_bool histrange_err     = FALSE;
    double   sum_histrange_err = 0.; /* histogram range error */
    double   stat_err          = 0.; /* statistical error */

    if (!parse_common_args(
                &argc, argv, PCA_CAN_VIEW, NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }

    gmx::ArrayRef<const std::string> xvgFiles = opt2fnsIfOptionSet("-f", NFILE, fnm);
    gmx::ArrayRef<const std::string> edrFiles = opt2fnsIfOptionSet("-g", NFILE, fnm);

    sim_data_init(&sim_data);
#if 0
    /* make linked list */
    lb = &lambda_head;
    lambda_data_init(lb, 0, 0);
    lb->next = lb;
    lb->prev = lb;
#endif


    nfile_tot = xvgFiles.size() + edrFiles.size();

    if (nfile_tot == 0)
    {
        gmx_fatal(FARGS, "No input files!");
    }

    if (nd < 0)
    {
        gmx_fatal(FARGS, "Can not have negative number of digits");
    }
    prec = std::pow(10.0, static_cast<double>(-nd));

    snew(partsum, (nbmax + 1) * (nbmax + 1));

    /* read in all files. First xvg files */
    for (const std::string& filenm : xvgFiles)
    {
        read_bar_xvg(filenm.c_str(), &temp, &sim_data);
    }
    /* then .edr files */
    for (const std::string& filenm : edrFiles)
    {
        read_barsim_edr(filenm.c_str(), &temp, &sim_data);
    }

    /* fix the times to allow for equilibration */
    sim_data_impose_times(&sim_data, begin, end);

    if (opt2bSet("-oh", NFILE, fnm))
    {
        sim_data_histogram(&sim_data, opt2fn("-oh", NFILE, fnm), nbin, oenv);
    }

    /* assemble the output structures from the lambdas */
    results = barres_list_create(&sim_data, &nresults, use_dhdl);

    sum_disc_err = barres_list_max_disc_err(results, nresults);

    if (nresults == 0)
    {
        printf("\nNo results to calculate.\n");
        return 0;
    }

    if (sum_disc_err > prec)
    {
        prec = sum_disc_err;
        nd   = static_cast<int>(std::ceil(-std::log10(prec)));
        printf("WARNING: setting the precision to %g because that is the minimum\n         "
               "reasonable number, given the expected discretization error.\n",
               prec);
    }


    /*sprintf(lamformat,"%%6.3f");*/
    sprintf(dgformat, "%%%d.%df", 3 + nd, nd);
    /* the format strings of the results in kT */
    sprintf(ktformat, "%%%d.%df", 5 + nd, nd);
    sprintf(sktformat, "%%%ds", 6 + nd);
    /* the format strings of the errors in kT */
    sprintf(kteformat, "%%%d.%df", 3 + nd, nd);
    sprintf(skteformat, "%%%ds", 4 + nd);
    sprintf(xvg2format, "%s %s\n", "%s", dgformat);
    sprintf(xvg3format, "%s %s %s\n", "%s", dgformat, dgformat);


    fpb = nullptr;
    if (opt2bSet("-o", NFILE, fnm))
    {
        sprintf(buf, "%s (%s)", "\\DeltaG", "kT");
        fpb = xvgropen_type(
                opt2fn("-o", NFILE, fnm), "Free energy differences", "\\lambda", buf, exvggtXYDY, oenv);
    }

    fpi = nullptr;
    if (opt2bSet("-oi", NFILE, fnm))
    {
        sprintf(buf, "%s (%s)", "\\DeltaG", "kT");
        fpi = xvgropen(opt2fn("-oi", NFILE, fnm), "Free energy integral", "\\lambda", buf, oenv);
    }


    if (nbmin > nbmax)
    {
        nbmin = nbmax;
    }

    /* first calculate results */
    bEE      = TRUE;
    disc_err = FALSE;
    for (f = 0; f < nresults; f++)
    {
        /* Determine the free energy difference with a factor of 10
         * more accuracy than requested for printing.
         */
        calc_bar(&(results[f]), 0.1 * prec, nbmin, nbmax, &bEE, partsum);

        if (results[f].dg_disc_err > prec / 10.)
        {
            disc_err = TRUE;
        }
        if (results[f].dg_histrange_err > prec / 10.)
        {
            histrange_err = TRUE;
        }
    }

    /* print results in kT */
    kT = gmx::c_boltz * temp;

    printf("\nTemperature: %g K\n", temp);

    printf("\nDetailed results in kT (see help for explanation):\n\n");
    printf("%6s ", " lam_A");
    printf("%6s ", " lam_B");
    printf(sktformat, "DG ");
    if (bEE)
    {
        printf(skteformat, "+/- ");
    }
    if (disc_err)
    {
        printf(skteformat, "disc ");
    }
    if (histrange_err)
    {
        printf(skteformat, "range ");
    }
    printf(sktformat, "s_A ");
    if (bEE)
    {
        printf(skteformat, "+/- ");
    }
    printf(sktformat, "s_B ");
    if (bEE)
    {
        printf(skteformat, "+/- ");
    }
    printf(sktformat, "stdev ");
    if (bEE)
    {
        printf(skteformat, "+/- ");
    }
    printf("\n");
    for (f = 0; f < nresults; f++)
    {
        lambda_vec_print_short(results[f].a->native_lambda, buf);
        printf("%s ", buf);
        lambda_vec_print_short(results[f].b->native_lambda, buf);
        printf("%s ", buf);
        printf(ktformat, results[f].dg);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].dg_err);
            printf(" ");
        }
        if (disc_err)
        {
            printf(kteformat, results[f].dg_disc_err);
            printf(" ");
        }
        if (histrange_err)
        {
            printf(kteformat, results[f].dg_histrange_err);
            printf(" ");
        }
        printf(ktformat, results[f].sa);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].sa_err);
            printf(" ");
        }
        printf(ktformat, results[f].sb);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].sb_err);
            printf(" ");
        }
        printf(ktformat, results[f].dg_stddev);
        printf(" ");
        if (bEE)
        {
            printf(kteformat, results[f].dg_stddev_err);
        }
        printf("\n");

        /* Check for negative relative entropy with a 95% certainty. */
        if (results[f].sa < -2 * results[f].sa_err || results[f].sb < -2 * results[f].sb_err)
        {
            result_OK = FALSE;
        }
    }

    if (!result_OK)
    {
        printf("\nWARNING: Some of these results violate the Second Law of "
               "Thermodynamics: \n"
               "         This is can be the result of severe undersampling, or "
               "(more likely)\n"
               "         there is something wrong with the simulations.\n");
    }


    /* final results in kJ/mol */
    printf("\n\nFinal results in kJ/mol:\n\n");
    dg_tot = 0;
    for (f = 0; f < nresults; f++)
    {

        if (fpi != nullptr)
        {
            lambda_vec_print_short(results[f].a->native_lambda, buf);
            fprintf(fpi, xvg2format, buf, dg_tot);
        }


        if (fpb != nullptr)
        {
            lambda_vec_print_intermediate(results[f].a->native_lambda, results[f].b->native_lambda, buf);

            fprintf(fpb, xvg3format, buf, results[f].dg, results[f].dg_err);
        }

        printf("point ");
        lambda_vec_print_short(results[f].a->native_lambda, buf);
        lambda_vec_print_short(results[f].b->native_lambda, buf2);
        printf("%s - %s", buf, buf2);
        printf(",   DG ");

        printf(dgformat, results[f].dg * kT);
        if (bEE)
        {
            printf(" +/- ");
            printf(dgformat, results[f].dg_err * kT);
        }
        if (histrange_err)
        {
            printf(" (max. range err. = ");
            printf(dgformat, results[f].dg_histrange_err * kT);
            printf(")");
            sum_histrange_err += results[f].dg_histrange_err * kT;
        }

        printf("\n");
        dg_tot += results[f].dg;
    }
    printf("\n");
    printf("total ");
    lambda_vec_print_short(results[0].a->native_lambda, buf);
    lambda_vec_print_short(results[nresults - 1].b->native_lambda, buf2);
    printf("%s - %s", buf, buf2);
    printf(",   DG ");

    printf(dgformat, dg_tot * kT);
    if (bEE)
    {
        stat_err = bar_err(nbmin, nbmax, partsum) * kT;
        printf(" +/- ");
        printf(dgformat, std::max(std::max(stat_err, sum_disc_err), sum_histrange_err));
    }
    printf("\n");
    if (disc_err)
    {
        printf("\nmaximum discretization error = ");
        printf(dgformat, sum_disc_err);
        if (bEE && stat_err < sum_disc_err)
        {
            printf("WARNING: discretization error (%g) is larger than statistical error.\n       "
                   "Decrease histogram spacing for more accurate results\n",
                   stat_err);
        }
    }
    if (histrange_err)
    {
        printf("\nmaximum histogram range error = ");
        printf(dgformat, sum_histrange_err);
        if (bEE && stat_err < sum_histrange_err)
        {
            printf("WARNING: histogram range error (%g) is larger than statistical error.\n       "
                   "Increase histogram range for more accurate results\n",
                   stat_err);
        }
    }
    printf("\n");


    if (fpi != nullptr)
    {
        lambda_vec_print_short(results[nresults - 1].b->native_lambda, buf);
        fprintf(fpi, xvg2format, buf, dg_tot);
        xvgrclose(fpi);
    }
    if (fpb != nullptr)
    {
        xvgrclose(fpb);
    }

    do_view(oenv, opt2fn_null("-o", NFILE, fnm), "-xydy");
    do_view(oenv, opt2fn_null("-oi", NFILE, fnm), "-xydy");

    return 0;
}
