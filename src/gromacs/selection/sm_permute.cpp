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
 * Implements the \p permute selection modifier.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/math/vec.h"
#include "gromacs/selection/position.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "selmethod.h"

/*! \internal \brief
 * Data structure for the \p permute selection modifier.
 */
typedef struct
{
    /** Positions to permute. */
    gmx_ana_pos_t    p;
    /** Number of elements in the permutation. */
    int              n;
    /** Array describing the permutation. */
    int             *perm;
    /** Array that has the permutation reversed. */
    int             *rperm;
} t_methoddata_permute;

/*! \brief
 * Allocates data for the \p permute selection modifier.
 *
 * \param[in]     npar  Not used (should be 2).
 * \param[in,out] param Method parameters (should point to a copy of
 *   \ref smparams_permute).
 * \returns Pointer to the allocated data (\p t_methoddata_permute).
 *
 * Allocates memory for a \p t_methoddata_permute structure.
 */
static void *
init_data_permute(int npar, gmx_ana_selparam_t *param);
/*! \brief
 * Initializes data for the \p permute selection modifier.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2).
 * \param[in] param Method parameters (should point to \ref smparams_permute).
 * \param[in] data  Should point to a \p t_methoddata_permute.
 * \returns   0 if the input permutation is valid, -1 on error.
 */
static void
init_permute(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/*! \brief
 * Initializes output for the \p permute selection modifier.
 *
 * \param[in]     top   Topology data structure.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_permute.
 */
static void
init_output_permute(t_topology *top, gmx_ana_selvalue_t *out, void *data);
/** Frees the memory allocated for the \p permute selection modifier. */
static void
free_data_permute(void *data);
static void
/*! \brief
 * Evaluates the \p permute selection modifier.
 *
 * \param[in]  top   Not used.
 * \param[in]  fr    Not used.
 * \param[in]  pbc   Not used.
 * \param[in]  p     Positions to permute (should point to \p data->p).
 * \param[out] out   Output data structure (\p out->u.p is used).
 * \param[in]  data  Should point to a \p t_methoddata_permute.
 * \returns    0 if \p p could be permuted, -1 on error.
 *
 * Returns -1 if the size of \p p is not divisible by the number of
 * elements in the permutation.
 */
evaluate_permute(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                 gmx_ana_pos_t *p, gmx_ana_selvalue_t *out, void *data);

/** Parameters for the \p permute selection modifier. */
static gmx_ana_selparam_t smparams_permute[] = {
    {NULL,       {POS_VALUE, -1, {NULL}}, NULL, SPAR_DYNAMIC | SPAR_VARNUM},
    {NULL,       {INT_VALUE, -1, {NULL}}, NULL, SPAR_VARNUM},
};

/** Help text for the \p permute selection modifier. */
static const char *const help_permute[] = {
    "::",
    "",
    "  permute P1 ... PN",
    "",
    "By default, all selections are evaluated such that the atom indices are",
    "returned in ascending order. This can be changed by appending",
    "[TT]permute P1 P2 ... PN[tt] to an expression.",
    "The [TT]Pi[tt] should form a permutation of the numbers 1 to N.",
    "This keyword permutes each N-position block in the selection such that",
    "the i'th position in the block becomes Pi'th.",
    "Note that it is the positions that are permuted, not individual atoms.",
    "A fatal error occurs if the size of the selection is not a multiple of n.",
    "It is only possible to permute the whole selection expression, not any",
    "subexpressions, i.e., the [TT]permute[tt] keyword should appear last in",
    "a selection.",
};

/** Selection method data for the \p permute modifier. */
gmx_ana_selmethod_t sm_permute = {
    "permute", POS_VALUE, SMETH_MODIFIER,
    asize(smparams_permute), smparams_permute,
    &init_data_permute,
    NULL,
    &init_permute,
    &init_output_permute,
    &free_data_permute,
    NULL,
    NULL,
    &evaluate_permute,
    {"POSEXPR permute P1 ... PN",
     "Permuting selections", asize(help_permute), help_permute},
};

static void *
init_data_permute(int /* npar */, gmx_ana_selparam_t *param)
{
    t_methoddata_permute *data = new t_methoddata_permute();
    data->n          = 0;
    data->perm       = NULL;
    data->rperm      = NULL;
    param[0].val.u.p = &data->p;
    return data;
}

static void
init_permute(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_permute *d = (t_methoddata_permute *)data;
    int                   i;

    d->n    = param[1].val.nr;
    d->perm = param[1].val.u.i;
    if (d->p.count() % d->n != 0)
    {
        GMX_THROW(gmx::InconsistentInputError(
                          gmx::formatString("The number of positions to be permuted is not divisible by %d", d->n)));
    }
    snew(d->rperm, d->n);
    for (i = 0; i < d->n; ++i)
    {
        d->rperm[i] = -1;
    }
    for (i = 0; i < d->n; ++i)
    {
        d->perm[i]--;
        if (d->perm[i] < 0 || d->perm[i] >= d->n)
        {
            GMX_THROW(gmx::InvalidInputError("Invalid permutation"));
        }
        if (d->rperm[d->perm[i]] >= 0)
        {
            GMX_THROW(gmx::InvalidInputError("Invalid permutation"));
        }
        d->rperm[d->perm[i]] = i;
    }
}

static void
init_output_permute(t_topology * /* top */, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_permute *d = (t_methoddata_permute *)data;
    int                   i, j, b;

    out->u.p->m.type = d->p.m.type;
    gmx_ana_pos_reserve_for_append(out->u.p, d->p.count(), d->p.m.b.nra,
                                   d->p.v != NULL, d->p.f != NULL);
    gmx_ana_pos_empty_init(out->u.p);
    for (i = 0; i < d->p.count(); i += d->n)
    {
        for (j = 0; j < d->n; ++j)
        {
            b = i + d->rperm[j];
            gmx_ana_pos_append_init(out->u.p, &d->p, b);
        }
    }
}

/*!
 * \param data Data to free (should point to a \p t_methoddata_permute).
 *
 * Frees the memory allocated for \c t_methoddata_permute.
 */
static void
free_data_permute(void *data)
{
    t_methoddata_permute *d = (t_methoddata_permute *)data;

    sfree(d->rperm);
    delete d;
}

static void
evaluate_permute(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                 gmx_ana_pos_t * /*p*/, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_permute *d = (t_methoddata_permute *)data;
    int                   i, j, b;
    int                   refid;

    if (d->p.count() % d->n != 0)
    {
        GMX_THROW(gmx::InconsistentInputError(
                          gmx::formatString("The number of positions to be permuted is not divisible by %d", d->n)));
    }
    gmx_ana_pos_empty(out->u.p);
    for (i = 0; i < d->p.count(); i += d->n)
    {
        for (j = 0; j < d->n; ++j)
        {
            b     = i + d->rperm[j];
            refid = d->p.m.refid[b];
            if (refid != -1)
            {
                /* De-permute the reference ID */
                refid = refid - (refid % d->n) + d->perm[refid % d->n];
            }
            gmx_ana_pos_append(out->u.p, &d->p, b, refid);
        }
    }
    gmx_ana_pos_append_finish(out->u.p);
}
