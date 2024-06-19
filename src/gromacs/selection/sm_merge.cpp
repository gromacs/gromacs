/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements the \p merge and \p plus selection modifiers.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include "gromacs/math/vec.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/topology/block.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/smalloc.h"

#include "selmethod.h"
#include "selmethod_impl.h"

struct gmx_mtop_t;

/*! \internal \brief
 * Data structure for the merging selection modifiers.
 */
typedef struct methoddata_merge
{
    /** Input positions. */
    gmx_ana_pos_t p1;
    /** Other input positions. */
    gmx_ana_pos_t p2;
    /** Stride for merging (\c stride values from \c p1 for each in \c p2). */
    int stride;
} t_methoddata_merge;

/** Allocates data for the merging selection modifiers. */
static void* init_data_merge(int npar, gmx_ana_selparam_t* param);
/*! \brief
 * Initializes data for the merging selection modifiers.
 *
 * \param[in] top   Not used.
 * \param[in] npar  Not used (should be 2 or 3).
 * \param[in] param Method parameters (should point to \ref smparams_merge).
 * \param[in] data  Should point to a \p t_methoddata_merge.
 * \returns   0 if everything is successful, -1 on error.
 */
static void init_merge(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Initializes output for the \p merge selection modifier. */
static void init_output_merge(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data);
/** Initializes output for the \p plus selection modifier. */
static void init_output_plus(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data);
/** Frees the memory allocated for the merging selection modifiers. */
static void free_data_merge(void* data);
/*! \brief
 * Evaluates the \p merge selection modifier.
 *
 * \param[in]  context Not used.
 * \param[in]  p     Positions to merge (should point to \p data->p1).
 * \param[out] out   Output data structure (\p out->u.p is used).
 * \param[in]  data  Should point to a \p t_methoddata_merge.
 */
static void evaluate_merge(const gmx::SelMethodEvalContext& context,
                           gmx_ana_pos_t*                   p,
                           gmx_ana_selvalue_t*              out,
                           void*                            data);
/*! \brief
 * Evaluates the \p plus selection modifier.
 *
 * \param[in]  context Not used.
 * \param[in]  p     Positions to merge (should point to \p data->p1).
 * \param[out] out   Output data structure (\p out->u.p is used).
 * \param[in]  data  Should point to a \p t_methoddata_merge.
 */
static void evaluate_plus(const gmx::SelMethodEvalContext& context,
                          gmx_ana_pos_t*                   p,
                          gmx_ana_selvalue_t*              out,
                          void*                            data);

/** Parameters for the merging selection modifiers. */
static gmx_ana_selparam_t smparams_merge[] = {
    { nullptr, { POS_VALUE, -1, { nullptr } }, nullptr, SPAR_DYNAMIC | SPAR_VARNUM },
    { nullptr, { POS_VALUE, -1, { nullptr } }, nullptr, SPAR_DYNAMIC | SPAR_VARNUM },
    { "stride", { INT_VALUE, 1, { nullptr } }, nullptr, SPAR_OPTIONAL },
};

//! Help title for the merging selection modifiers.
static const char helptitle_merge[] = "Merging selections";
//! Help text for the merging selection modifiers.
static const char* const help_merge[] = {
    "::",
    "",
    "  POSEXPR merge POSEXPR [stride INT]",
    "  POSEXPR merge POSEXPR [merge POSEXPR ...]",
    "  POSEXPR plus POSEXPR [plus POSEXPR ...]",
    "",
    "Basic selection keywords can only create selections where each atom",
    "occurs at most once. The [TT]merge[tt] and [TT]plus[tt] selection",
    "keywords can be used to work around this limitation. Both create",
    "a selection that contains the positions from all the given position",
    "expressions, even if they contain duplicates.",
    "The difference between the two is that [TT]merge[tt] expects two or more",
    "selections with the same number of positions, and the output contains",
    "the input positions selected from each expression in turn, i.e.,",
    "the output is like A1 B1 A2 B2 and so on. It is also possible to merge",
    "selections of unequal size as long as the size of the first is a",
    "multiple of the second one. The [TT]stride[tt] parameter can be used",
    "to explicitly provide this multiplicity.",
    "[TT]plus[tt] simply concatenates the positions after each other, and",
    "can work also with selections of different sizes.",
    "These keywords are valid only at the selection level, not in any",
    "subexpressions.",
};

/** Selection method data for the \p plus modifier. */
gmx_ana_selmethod_t sm_merge = {
    "merge",
    POS_VALUE,
    SMETH_MODIFIER,
    asize(smparams_merge),
    smparams_merge,
    &init_data_merge,
    nullptr,
    &init_merge,
    &init_output_merge,
    &free_data_merge,
    nullptr,
    nullptr,
    &evaluate_merge,
    { "merge POSEXPR", helptitle_merge, asize(help_merge), help_merge },
};

/** Selection method data for the \p plus modifier. */
gmx_ana_selmethod_t sm_plus = {
    "plus",
    POS_VALUE,
    SMETH_MODIFIER,
    asize(smparams_merge) - 1,
    smparams_merge,
    &init_data_merge,
    nullptr,
    &init_merge,
    &init_output_plus,
    &free_data_merge,
    nullptr,
    nullptr,
    &evaluate_plus,
    { "plus POSEXPR", helptitle_merge, asize(help_merge), help_merge },
};

/*!
 * \param[in]     npar  Should be 2 for \c plus and 3 for \c merge.
 * \param[in,out] param Method parameters (should point to a copy of
 *   \ref smparams_merge).
 * \returns Pointer to the allocated data (\p t_methoddata_merge).
 *
 * Allocates memory for a \p t_methoddata_merge structure.
 */
static void* init_data_merge(int npar, gmx_ana_selparam_t* param)
{
    t_methoddata_merge* data = new t_methoddata_merge();
    data->stride             = 0;
    param[0].val.u.p         = &data->p1;
    param[1].val.u.p         = &data->p2;
    if (npar > 2)
    {
        param[2].val.u.i = &data->stride;
    }
    return data;
}

static void init_merge(const gmx_mtop_t* /* top */, int /* npar */, gmx_ana_selparam_t* /* param */, void* data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);

    if (d->stride < 0)
    {
        GMX_THROW(gmx::InvalidInputError("Stride for merging should be positive"));
    }
    /* If no stride given, deduce it from the input sizes */
    if (d->stride == 0)
    {
        d->stride = d->p1.count() / d->p2.count();
    }
    if (d->p1.count() != d->stride * d->p2.count())
    {
        GMX_THROW(gmx::InconsistentInputError(
                "The number of positions to be merged are not compatible"));
    }
}

/*! \brief
 * Does common initialization to all merging modifiers.
 *
 * \param[in]     top   Topology data structure.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_merge.
 */
static void init_output_common(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);

    GMX_UNUSED_VALUE(top);
    if (d->p1.m.type != d->p2.m.type)
    {
        /* TODO: Maybe we could pick something else here? */
        out->u.p->m.type = INDEX_UNKNOWN;
    }
    else
    {
        out->u.p->m.type = d->p1.m.type;
    }
    gmx_ana_pos_reserve_for_append(out->u.p,
                                   d->p1.count() + d->p2.count(),
                                   d->p1.m.b.nra + d->p2.m.b.nra,
                                   d->p1.v != nullptr,
                                   d->p1.f != nullptr);
    gmx_ana_pos_empty_init(out->u.p);
}

/*!
 * \param[in]     top   Topology data structure.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_merge.
 */
static void init_output_merge(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);
    int                 i, j;

    init_output_common(top, out, data);
    for (i = 0; i < d->p2.count(); ++i)
    {
        for (j = 0; j < d->stride; ++j)
        {
            gmx_ana_pos_append_init(out->u.p, &d->p1, d->stride * i + j);
        }
        gmx_ana_pos_append_init(out->u.p, &d->p2, i);
    }
}

/*!
 * \param[in]     top   Topology data structure.
 * \param[in,out] out   Pointer to output data structure.
 * \param[in,out] data  Should point to \c t_methoddata_merge.
 */
static void init_output_plus(const gmx_mtop_t* top, gmx_ana_selvalue_t* out, void* data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);
    int                 i;

    init_output_common(top, out, data);
    for (i = 0; i < d->p1.count(); ++i)
    {
        gmx_ana_pos_append_init(out->u.p, &d->p1, i);
    }
    for (i = 0; i < d->p2.count(); ++i)
    {
        gmx_ana_pos_append_init(out->u.p, &d->p2, i);
    }
}

/*!
 * \param data Data to free (should point to a \p t_methoddata_merge).
 *
 * Frees the memory allocated for \c t_methoddata_merge.
 */
static void free_data_merge(void* data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);
    delete d;
}

static void evaluate_merge(const gmx::SelMethodEvalContext& /*context*/,
                           gmx_ana_pos_t* /* p */,
                           gmx_ana_selvalue_t* out,
                           void*               data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);
    int                 i, j;
    int                 refid;

    if (d->p1.count() != d->stride * d->p2.count())
    {
        GMX_THROW(gmx::InconsistentInputError(
                "The number of positions to be merged are not compatible"));
    }
    gmx_ana_pos_empty(out->u.p);
    for (i = 0; i < d->p2.count(); ++i)
    {
        for (j = 0; j < d->stride; ++j)
        {
            refid = d->p1.m.refid[d->stride * i + j];
            if (refid != -1)
            {
                refid = (d->stride + 1) * (refid / d->stride) + (refid % d->stride);
            }
            gmx_ana_pos_append(out->u.p, &d->p1, d->stride * i + j, refid);
        }
        refid = (d->stride + 1) * d->p2.m.refid[i] + d->stride;
        gmx_ana_pos_append(out->u.p, &d->p2, i, refid);
    }
    gmx_ana_pos_append_finish(out->u.p);
}

static void evaluate_plus(const gmx::SelMethodEvalContext& /*context*/,
                          gmx_ana_pos_t* /* p */,
                          gmx_ana_selvalue_t* out,
                          void*               data)
{
    t_methoddata_merge* d = static_cast<t_methoddata_merge*>(data);
    int                 i;
    int                 refid;

    gmx_ana_pos_empty(out->u.p);
    for (i = 0; i < d->p1.count(); ++i)
    {
        refid = d->p1.m.refid[i];
        gmx_ana_pos_append(out->u.p, &d->p1, i, refid);
    }
    for (i = 0; i < d->p2.count(); ++i)
    {
        refid = d->p2.m.refid[i];
        if (refid != -1)
        {
            refid += d->p1.m.b.nr;
        }
        gmx_ana_pos_append(out->u.p, &d->p2, i, refid);
    }
    gmx_ana_pos_append_finish(out->u.p);
}
