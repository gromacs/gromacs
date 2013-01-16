/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * \brief Implementation of distance-based selection methods.
 *
 * This file implements the \p distance, \p mindistance and \p within
 * selection methods.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <macros.h>
#include <pbc.h>
#include <smalloc.h>
#include <vec.h>

#include <nbsearch.h>
#include <position.h>
#include <selmethod.h>

/*! \internal \brief
 * Data structure for distance-based selection method.
 *
 * The same data structure is used by all the distance-based methods.
 */
typedef struct
{
    /** Cutoff distance. */
    real                cutoff;
    /** Positions of the reference points. */
    gmx_ana_pos_t       p;
    /** Neighborhood search data. */
    gmx_ana_nbsearch_t *nb;
} t_methoddata_distance;

/** Allocates data for distance-based selection methods. */
static void *
init_data_common(int npar, gmx_ana_selparam_t *param);
/** Initializes a distance-based selection method. */
static int
init_common(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the data allocated for a distance-based selection method. */
static void
free_data_common(void *data);
/** Initializes the evaluation of a distance-based within selection method for a frame. */
static int
init_frame_common(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data);
/** Evaluates the \p distance selection method. */
static int
evaluate_distance(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p within selection method. */
static int
evaluate_within(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);

/** Parameters for the \p distance selection method. */
static gmx_ana_selparam_t smparams_distance[] = {
    {"cutoff", {REAL_VALUE, 1, {NULL}}, NULL, SPAR_OPTIONAL},
    {"from",   {POS_VALUE,  1, {NULL}}, NULL, SPAR_DYNAMIC},
};

/** Parameters for the \p mindistance selection method. */
static gmx_ana_selparam_t smparams_mindistance[] = {
    {"cutoff", {REAL_VALUE, 1, {NULL}}, NULL, SPAR_OPTIONAL},
    {"from",   {POS_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
};

/** Parameters for the \p within selection method. */
static gmx_ana_selparam_t smparams_within[] = {
    {NULL, {REAL_VALUE,  1, {NULL}}, NULL, 0},
    {"of", {POS_VALUE,  -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
};

/** Help text for the distance selection methods. */
static const char *help_distance[] = {
    "DISTANCE-BASED SELECTION KEYWORDS[PAR]",

    "[TT]distance from POS [cutoff REAL][tt][BR]",
    "[TT]mindistance from POS_EXPR [cutoff REAL][tt][BR]",
    "[TT]within REAL of POS_EXPR[tt][PAR]",

    "[TT]distance[tt] and [TT]mindistance[tt] calculate the distance from the",
    "given position(s), the only difference being in that [TT]distance[tt]",
    "only accepts a single position, while any number of positions can be",
    "given for [TT]mindistance[tt], which then calculates the distance to the",
    "closest position.",
    "[TT]within[tt] directly selects atoms that are within [TT]REAL[tt] of",
    "[TT]POS_EXPR[tt].[PAR]",

    "For the first two keywords, it is possible to specify a cutoff to speed",
    "up the evaluation: all distances above the specified cutoff are",
    "returned as equal to the cutoff.",
    "Currently, this does nothing, but in the future, it allows the use of",
    "grid-based neighborhood search techniques.",
};

/** \internal Selection method data for the \p distance method. */
gmx_ana_selmethod_t sm_distance = {
    "distance", REAL_VALUE, SMETH_DYNAMIC,
    asize(smparams_distance), smparams_distance,
    &init_data_common,
    NULL,
    &init_common,
    NULL,
    &free_data_common,
    &init_frame_common,
    NULL,
    &evaluate_distance,
    {"distance from POS [cutoff REAL]", asize(help_distance), help_distance},
};

/** \internal Selection method data for the \p distance method. */
gmx_ana_selmethod_t sm_mindistance = {
    "mindistance", REAL_VALUE, SMETH_DYNAMIC,
    asize(smparams_mindistance), smparams_mindistance,
    &init_data_common,
    NULL,
    &init_common,
    NULL,
    &free_data_common,
    &init_frame_common,
    NULL,
    &evaluate_distance,
    {"mindistance from POS_EXPR [cutoff REAL]", asize(help_distance), help_distance},
};

/** \internal Selection method data for the \p within method. */
gmx_ana_selmethod_t sm_within = {
    "within", GROUP_VALUE, SMETH_DYNAMIC,
    asize(smparams_within), smparams_within,
    &init_data_common,
    NULL,
    &init_common,
    NULL,
    &free_data_common,
    &init_frame_common,
    NULL,
    &evaluate_within,
    {"within REAL of POS_EXPR", asize(help_distance), help_distance},
};

/*!
 * \param[in]     npar  Not used (should be 2).
 * \param[in,out] param Method parameters (should point to one of the distance
 *   parameter arrays).
 * \returns       Pointer to the allocated data (\c t_methoddata_distance).
 *
 * Allocates memory for a \c t_methoddata_distance structure and
 * initializes the parameter as follows:
 *  - the first parameter defines the value for
 *    \c t_methoddata_distance::cutoff.
 *  - the second parameter defines the reference positions and the value is
 *    stored in \c t_methoddata_distance::p.
 */
static void *
init_data_common(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_distance *data;

    snew(data, 1);
    data->cutoff     = -1;
    param[0].val.u.r = &data->cutoff;
    param[1].val.u.p = &data->p;
    return data;
}

/*!
 * \param   top   Not used.
 * \param   npar  Not used (should be 2).
 * \param   param Method parameters (should point to one of the distance
 *   parameter arrays).
 * \param   data  Pointer to \c t_methoddata_distance to initialize.
 * \returns 0 on success, a non-zero error code on failure.
 *
 * Initializes the neighborhood search data structure
 * (\c t_methoddata_distance::nb).
 * Also checks that the cutoff is valid.
 */
static int
init_common(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_distance *d = (t_methoddata_distance *)data;

    if ((param[0].flags & SPAR_SET) && d->cutoff <= 0)
    {
        fprintf(stderr, "error: distance cutoff should be > 0");
        return -1;
    }
    return gmx_ana_nbsearch_create(&d->nb, d->cutoff, d->p.nr);
}

/*!
 * \param data Data to free (should point to a \c t_methoddata_distance).
 *
 * Frees the memory allocated for \c t_methoddata_distance::xref and
 * \c t_methoddata_distance::nb.
 */
static void
free_data_common(void *data)
{
    t_methoddata_distance *d = (t_methoddata_distance *)data;

    if (d->nb)
    {
        gmx_ana_nbsearch_free(d->nb);
    }
}

/*!
 * \param[in]  top  Not used.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \c t_methoddata_distance.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Initializes the neighborhood search for the current frame.
 */
static int
init_frame_common(t_topology *top, t_trxframe *fr, t_pbc *pbc, void *data)
{
    t_methoddata_distance *d = (t_methoddata_distance *)data;

    return gmx_ana_nbsearch_pos_init(d->nb, pbc, &d->p);
}

/*!
 * See sel_updatefunc_pos() for description of the parameters.
 * \p data should point to a \c t_methoddata_distance.
 *
 * Calculates the distance of each position from \c t_methoddata_distance::p
 * and puts them in \p out->u.r.
 */
static int
evaluate_distance(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_distance *d = (t_methoddata_distance *)data;
    int                    b, i;
    real                   n;

    out->nr = pos->g->isize;
    for (b = 0; b < pos->nr; ++b)
    {
        n = gmx_ana_nbsearch_pos_mindist(d->nb, pos, b);
        for (i = pos->m.mapb.index[b]; i < pos->m.mapb.index[b+1]; ++i)
        {
            out->u.r[i] = n;
        }
    }
    return 0;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_distance.
 *
 * Finds the atoms that are closer than the defined cutoff to
 * \c t_methoddata_distance::xref and puts them in \p out.g.
 */
static int
evaluate_within(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_distance *d = (t_methoddata_distance *)data;
    int                    b;

    out->u.g->isize = 0;
    for (b = 0; b < pos->nr; ++b)
    {
        if (gmx_ana_nbsearch_pos_is_within(d->nb, pos, b))
        {
            gmx_ana_pos_append(NULL, out->u.g, pos, b, 0);
        }
    }
    return 0;
}
