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
 * Implements distance-based selection methods.
 *
 * This file implements the \p distance, \p mindistance and \p within
 * selection methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/position.h"
#include "gromacs/utility/exceptions.h"

#include "selmethod.h"

/*! \internal
 * \brief
 * Data structure for distance-based selection method.
 *
 * The same data structure is used by all the distance-based methods.
 *
 * \ingroup module_selection
 */
struct t_methoddata_distance
{
    t_methoddata_distance() : cutoff(-1.0)
    {
    }

    /** Cutoff distance. */
    real                             cutoff;
    /** Positions of the reference points. */
    gmx_ana_pos_t                    p;
    /** Neighborhood search data. */
    gmx::AnalysisNeighborhood        nb;
    /** Neighborhood search for an invididual frame. */
    gmx::AnalysisNeighborhoodSearch  nbsearch;
};

/*! \brief
 * Allocates data for distance-based selection methods.
 *
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
init_data_common(int npar, gmx_ana_selparam_t *param);
/*! \brief
 * Initializes a distance-based selection method.
 *
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
static void
init_common(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Frees the data allocated for a distance-based selection method. */
static void
free_data_common(void *data);
/*! \brief
 * Initializes the evaluation of a distance-based within selection method for a
 * frame.
 *
 * \param[in]  top  Not used.
 * \param[in]  fr   Current frame.
 * \param[in]  pbc  PBC structure.
 * \param      data Should point to a \c t_methoddata_distance.
 * \returns    0 on success, a non-zero error code on error.
 *
 * Initializes the neighborhood search for the current frame.
 */
static void
init_frame_common(t_topology *top, t_trxframe * fr, t_pbc *pbc, void *data);
/** Evaluates the \p distance selection method. */
static void
evaluate_distance(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p within selection method. */
static void
evaluate_within(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
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

//! Help title for distance selection methods.
static const char        helptitle_distance[] = "Selecting based on distance";
//! Help text for distance selection methods.
static const char *const help_distance[] = {
    "::",
    "",
    "  distance from POS [cutoff REAL]",
    "  mindistance from POS_EXPR [cutoff REAL]",
    "  within REAL of POS_EXPR",
    "",
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
};

/** Selection method data for the \p distance method. */
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
    {"distance from POS [cutoff REAL]",
     helptitle_distance, asize(help_distance), help_distance},
};

/** Selection method data for the \p distance method. */
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
    {"mindistance from POS_EXPR [cutoff REAL]",
     helptitle_distance, asize(help_distance), help_distance},
};

/** Selection method data for the \p within method. */
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
    {"within REAL of POS_EXPR",
     helptitle_distance, asize(help_distance), help_distance},
};

static void *
init_data_common(int /* npar */, gmx_ana_selparam_t *param)
{
    t_methoddata_distance *data = new t_methoddata_distance();
    param[0].val.u.r = &data->cutoff;
    param[1].val.u.p = &data->p;
    return data;
}

static void
init_common(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_distance *d = static_cast<t_methoddata_distance *>(data);

    if ((param[0].flags & SPAR_SET) && d->cutoff <= 0)
    {
        GMX_THROW(gmx::InvalidInputError("Distance cutoff should be > 0"));
    }
    d->nb.setCutoff(d->cutoff);
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
    delete static_cast<t_methoddata_distance *>(data);
}

static void
init_frame_common(t_topology * /* top */, t_trxframe * /* fr */, t_pbc *pbc, void *data)
{
    t_methoddata_distance *d = static_cast<t_methoddata_distance *>(data);

    d->nbsearch.reset();
    gmx::AnalysisNeighborhoodPositions pos(d->p.x, d->p.count());
    d->nbsearch = d->nb.initSearch(pbc, pos);
}

/*!
 * See sel_updatefunc_pos() for description of the parameters.
 * \p data should point to a \c t_methoddata_distance.
 *
 * Calculates the distance of each position from \c t_methoddata_distance::p
 * and puts them in \p out->u.r.
 */
static void
evaluate_distance(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_distance *d = static_cast<t_methoddata_distance *>(data);

    out->nr = pos->count();
    for (int i = 0; i < pos->count(); ++i)
    {
        out->u.r[i] = d->nbsearch.minimumDistance(pos->x[i]);
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_distance.
 *
 * Finds the atoms that are closer than the defined cutoff to
 * \c t_methoddata_distance::xref and puts them in \p out.g.
 */
static void
evaluate_within(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_distance *d = static_cast<t_methoddata_distance *>(data);

    out->u.g->isize = 0;
    for (int b = 0; b < pos->count(); ++b)
    {
        if (d->nbsearch.isWithin(pos->x[b]))
        {
            gmx_ana_pos_add_to_group(out->u.g, pos, b);
        }
    }
}
