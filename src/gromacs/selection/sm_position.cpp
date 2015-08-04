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
 * Implements position evaluation selection methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"

#include "keywords.h"
#include "poscalc.h"
#include "selelem.h"
#include "selmethod.h"

/*! \internal \brief
 * Data structure for position keyword evaluation.
 */
typedef struct
{
    /** Position calculation collection to use. */
    gmx::PositionCalculationCollection *pcc;
    /** Index group for which the center should be evaluated. */
    gmx_ana_index_t                     g;
    /** Position evaluation data structure. */
    gmx_ana_poscalc_t                  *pc;
    /** true if periodic boundary conditions should be used. */
    bool                                bPBC;
    /** Type of positions to calculate. */
    char                               *type;
    /** Flags for the position calculation. */
    int                                 flags;
} t_methoddata_pos;

/** Allocates data for position evaluation selection methods. */
static void *
init_data_pos(int npar, gmx_ana_selparam_t *param);
/** Sets the position calculation collection for position evaluation selection methods. */
static void
set_poscoll_pos(gmx::PositionCalculationCollection *pcc, void *data);
/*! \brief
 * Initializes position evaluation keywords.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used.
 * \param[in] param Not used.
 * \param[in,out] data  Should point to \c t_methoddata_pos.
 * \returns       0 on success, a non-zero error code on error.
 *
 * The \c t_methoddata_pos::type field should have been initialized
 * externally using _gmx_selelem_set_kwpos_type().
 */
static void
init_kwpos(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/*! \brief
 * Initializes the \p cog selection method.
 *
 * \param[in]     top   Topology data structure.
 * \param[in]     npar  Not used.
 * \param[in]     param Not used.
 * \param[in,out] data  Should point to \c t_methoddata_pos.
 * \returns       0 on success, a non-zero error code on error.
 */
static void
init_cog(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/*! \brief
 * Initializes the \p cog selection method.
 *
 * \param[in]     top   Topology data structure.
 * \param[in]     npar  Not used.
 * \param[in]     param Not used.
 * \param[in,out] data  Should point to \c t_methoddata_pos.
 * \returns       0 on success, a non-zero error code on error.
 */
static void
init_com(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/*! \brief
 * Initializes output for position evaluation selection methods.
 *
 * \param[in]     top   Topology data structure.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_pos.
 * \returns       0 for success.
 */
static void
init_output_pos(t_topology *top, gmx_ana_selvalue_t *out, void *data);
/** Frees the data allocated for position evaluation selection methods. */
static void
free_data_pos(void *data);
/** Evaluates position evaluation selection methods. */
static void
evaluate_pos(t_topology * /* top */, t_trxframe *fr, t_pbc *pbc, gmx_ana_index_t * /* g */, gmx_ana_selvalue_t *out, void *data);

/** Parameters for position keyword evaluation. */
static gmx_ana_selparam_t smparams_keyword_pos[] = {
    {NULL,   {GROUP_VALUE, 1, {NULL}}, NULL, SPAR_DYNAMIC},
};

/** Parameters for the \p cog and \p com selection methods. */
static gmx_ana_selparam_t smparams_com[] = {
    {"of",   {GROUP_VALUE, 1, {NULL}}, NULL, SPAR_DYNAMIC},
    {"pbc",  {NO_VALUE,    0, {NULL}}, NULL, 0},
};

/** Selection method data for position keyword evaluation. */
gmx_ana_selmethod_t sm_keyword_pos = {
    "kw_pos", POS_VALUE, SMETH_DYNAMIC | SMETH_VARNUMVAL | SMETH_ALLOW_UNSORTED,
    asize(smparams_keyword_pos), smparams_keyword_pos,
    &init_data_pos,
    &set_poscoll_pos,
    &init_kwpos,
    &init_output_pos,
    &free_data_pos,
    NULL,
    &evaluate_pos,
    NULL,
    {NULL, NULL, 0, NULL},
};

/** Selection method data for the \p cog method. */
gmx_ana_selmethod_t sm_cog = {
    "cog", POS_VALUE, SMETH_DYNAMIC | SMETH_SINGLEVAL,
    asize(smparams_com), smparams_com,
    &init_data_pos,
    &set_poscoll_pos,
    &init_cog,
    &init_output_pos,
    &free_data_pos,
    NULL,
    &evaluate_pos,
    NULL,
    {"cog of ATOM_EXPR [pbc]", NULL, 0, NULL},
};

/** Selection method data for the \p com method. */
gmx_ana_selmethod_t sm_com = {
    "com", POS_VALUE, SMETH_REQTOP | SMETH_DYNAMIC | SMETH_SINGLEVAL,
    asize(smparams_com), smparams_com,
    &init_data_pos,
    &set_poscoll_pos,
    &init_com,
    &init_output_pos,
    &free_data_pos,
    NULL,
    &evaluate_pos,
    NULL,
    {"com of ATOM_EXPR [pbc]", NULL, 0, NULL},
};

/*!
 * \param[in]     npar  Should be 1 or 2.
 * \param[in,out] param Method parameters (should point to
 *   \ref smparams_keyword_pos or \ref smparams_com).
 * \returns       Pointer to the allocated data (\c t_methoddata_pos).
 *
 * Allocates memory for a \c t_methoddata_pos structure and initializes
 * the first parameter to define the value for \c t_methoddata_pos::g.
 * If a second parameter is present, it is used for setting the
 * \c t_methoddata_pos::bPBC flag.
 */
static void *
init_data_pos(int npar, gmx_ana_selparam_t *param)
{
    t_methoddata_pos *data;

    snew(data, 1);
    param[0].val.u.g = &data->g;
    if (npar > 1)
    {
        param[1].val.u.b = &data->bPBC;
    }
    data->pc       = NULL;
    data->bPBC     = false;
    data->type     = NULL;
    data->flags    = -1;
    return data;
}

/*!
 * \param[in]     pcc   Position calculation collection to use.
 * \param[in,out] data  Should point to \c t_methoddata_pos.
 */
static void
set_poscoll_pos(gmx::PositionCalculationCollection *pcc, void *data)
{
    ((t_methoddata_pos *)data)->pcc = pcc;
}

bool
_gmx_selelem_is_default_kwpos(const gmx::SelectionTreeElement &sel)
{
    if (sel.type != SEL_EXPRESSION || !sel.u.expr.method
        || sel.u.expr.method->name != sm_keyword_pos.name)
    {
        return false;
    }

    t_methoddata_pos *d = static_cast<t_methoddata_pos *>(sel.u.expr.mdata);
    return d->type == NULL;
}

/*!
 * \param[in,out] sel   Selection element to initialize.
 * \param[in]     type  One of the enum values acceptable for
 *     PositionCalculationCollection::typeFromEnum().
 *
 * Initializes the reference position type for position evaluation.
 * If called multiple times, the first setting takes effect, and later calls
 * are neglected.
 */
void
_gmx_selelem_set_kwpos_type(gmx::SelectionTreeElement *sel, const char *type)
{
    t_methoddata_pos *d = (t_methoddata_pos *)sel->u.expr.mdata;

    if (sel->type != SEL_EXPRESSION || !sel->u.expr.method
        || sel->u.expr.method->name != sm_keyword_pos.name)
    {
        return;
    }
    if (!d->type && type)
    {
        d->type  = gmx_strdup(type);
        /* FIXME: It would be better not to have the string here hardcoded. */
        if (type[0] != 'a')
        {
            sel->u.expr.method->flags |= SMETH_REQTOP;
        }
    }
}

/*!
 * \param[in,out] sel   Selection element to initialize.
 * \param[in]     flags Default completion flags
 *     (see PositionCalculationCollection::typeFromEnum()).
 *
 * Initializes the flags for position evaluation.
 * If called multiple times, the first setting takes effect, and later calls
 * are neglected.
 */
void
_gmx_selelem_set_kwpos_flags(gmx::SelectionTreeElement *sel, int flags)
{
    t_methoddata_pos *d = (t_methoddata_pos *)sel->u.expr.mdata;

    if (sel->type != SEL_EXPRESSION || !sel->u.expr.method
        || sel->u.expr.method->name != sm_keyword_pos.name)
    {
        return;
    }
    if (d->flags == -1)
    {
        d->flags = flags;
    }
}

static void
init_kwpos(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    if (!(param[0].flags & SPAR_DYNAMIC))
    {
        d->flags &= ~(POS_DYNAMIC | POS_MASKONLY);
    }
    else if (!(d->flags & POS_MASKONLY))
    {
        d->flags |= POS_DYNAMIC;
    }
    d->pc = d->pcc->createCalculationFromEnum(d->type, d->flags);
    gmx_ana_poscalc_set_maxindex(d->pc, &d->g);
}

static void
init_cog(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    d->flags = (param[0].flags & SPAR_DYNAMIC) ? POS_DYNAMIC : 0;
    d->pc    = d->pcc->createCalculation(d->bPBC ? POS_ALL_PBC : POS_ALL, d->flags);
    gmx_ana_poscalc_set_maxindex(d->pc, &d->g);
}

static void
init_com(t_topology * /* top */, int /* npar */, gmx_ana_selparam_t *param, void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    d->flags  = (param[0].flags & SPAR_DYNAMIC) ? POS_DYNAMIC : 0;
    d->flags |= POS_MASS;
    d->pc     = d->pcc->createCalculation(d->bPBC ? POS_ALL_PBC : POS_ALL, d->flags);
    gmx_ana_poscalc_set_maxindex(d->pc, &d->g);
}

static void
init_output_pos(t_topology * /* top */, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    gmx_ana_poscalc_init_pos(d->pc, out->u.p);
}

/*!
 * \param data Data to free (should point to a \c t_methoddata_pos).
 *
 * Frees the memory allocated for \c t_methoddata_pos::g and
 * \c t_methoddata_pos::pc.
 */
static void
free_data_pos(void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    sfree(d->type);
    gmx_ana_poscalc_free(d->pc);
    sfree(d);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_pos.
 *
 * Calculates the positions using \c t_methoddata_pos::pc for the index group
 * in \c t_methoddata_pos::g and stores the results in \p out->u.p.
 */
static void
evaluate_pos(t_topology * /* top */, t_trxframe *fr, t_pbc *pbc,
             gmx_ana_index_t * /* g */, gmx_ana_selvalue_t *out, void *data)
{
    t_methoddata_pos *d = (t_methoddata_pos *)data;

    gmx_ana_poscalc_update(d->pc, out->u.p, &d->g, fr, pbc);
}
