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
/*! \internal
 * \page page_module_selection_insolidangle Selection method: insolidangle
 *
 * This method selects a subset of particles that are located in a solid
 * angle defined by a center and a set of points.
 * The solid angle is constructed as a union of small cones whose axis
 * goes through the center and a point.
 * So there's such a cone for each position, and a
 * point is in the solid angle if it lies within any of these cones.
 * The width of the cones can be adjusted.
 *
 * The method is implemented by partitioning the surface of the unit sphere
 * into bins using the polar coordinates \f$(\theta, \phi)\f$.
 * The partitioning is always uniform in the zenith angle \f$\theta\f$,
 * while the partitioning in the azimuthal angle \f$\phi\f$ varies.
 * For each reference point, the unit vector from the center to the point
 * is constructed, and it is stored in all the bins that overlap with the
 * cone defined by the point.
 * Bins that are completely covered by a single cone are marked as such.
 * Checking whether a point is in the solid angle is then straightforward
 * with this data structure: one finds the bin that corresponds to the point,
 * and checks whether the bin is completely covered. If it is not, one
 * additionally needs to check whether it is within the specified cutoff of
 * any of the stored points.
 *
 * The above construction gives quite a lot of flexibility for constructing
 * the bins without modifying the rest of the code.
 * The current (quite inefficient) implementation is discussed below, but
 * it should be optimized to get the most out of the code.
 *
 * The current way of constructing the bins constructs the boundaries
 * statically: the bin size in the zenith direction is set to approximately
 * half the angle cutoff, and the bins in the azimuthal direction have
 * sizes such that the shortest edge of the bin is approximately equal to
 * half the angle cutoff (for the regions close to the poles, a single bin
 * is used).
 * Each reference point is then added to the bins as follows:
 *  -# Find the zenith angle range that is spanned by the cone centered at the
 *     point (this is simple addition/subtraction).
 *  -# Calculate the maximal span of the cone in the azimuthal direction using
 *     the formula
 *     \f[\sin \Delta \phi_{max} = \frac{\sin \alpha}{\sin \theta}\f]
 *     (a sine formula in spherical coordinates),
 *     where \f$\alpha\f$ is the width of the cone and \f$\theta\f$ is the
 *     zenith angle of the cone center.
 *     Similarly, the zenith angle at which this extent is achieved is
 *     calculated using
 *     \f[\cos \theta_{max} = \frac{\cos \theta}{\cos \alpha}\f]
 *     (Pythagoras's theorem in spherical coordinates).
 *  -# For each zenith angle bin that is at least partially covered by the
 *     cone, calculate the span of the cone at the edges using
 *     \f[\sin^2 \frac{\Delta \phi}{2} = \frac{\sin^2 \frac{\alpha}{2} - \sin^2 \frac{\theta -
 * \theta'}{2}}{\sin \theta \sin \theta'}\f] (distance in spherical geometry), where \f$\theta'\f$
 * is the zenith angle of the bin edge. Treat zenith angle bins that are completely covered by the
 * cone (in the case that the cone is centered close to the pole) as a special case.
 *  -# Using the values calculated above, loop through the azimuthal bins that
 *     are partially or completely covered by the cone and update them.
 *
 * The total solid angle (for covered fraction calculations) is estimated by
 * taking the total area of completely covered bins plus
 * half the area of partially covered bins.
 * The second one is an approximation, but should give reasonable estimates
 * for the averages as well as in cases where the bin size is small.
 */
/*! \internal \file
 * \brief
 * Implements the \ref sm_insolidangle "insolidangle" selection method.
 *
 * \todo
 * The implementation could be optimized quite a bit.
 *
 * \todo
 * Move the covered fraction stuff somewhere else and make it more generic
 * (along the lines it is handled in selection.h and trajana.h in the old C
 * API).
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <cmath>

#include <algorithm>
#include <memory>

#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

#include "selelem.h"
#include "selmethod.h"
#include "selmethod_impl.h"

struct gmx_mtop_t;

using std::max;
using std::min;

/*! \internal
 * \brief
 * Internal data structure for the \p insolidangle selection method.
 *
 * \see \c t_partition
 *
 * \ingroup module_selection
 */
typedef struct
{
    /** Left edge of the partition. */
    real left;
    /** Bin index corresponding to this partition. */
    int bin;
} t_partition_item;

/*! \internal
 * \brief
 * Internal data structure for the \p insolidangle selection method.
 *
 * Describes the surface partitioning within one slice along the zenith angle.
 * The slice from azimuthal angle \p p[i].left to \p p[i+1].left belongs to
 * bin \p p[i].bin.
 *
 * \ingroup module_selection
 */
typedef struct partition
{
    /** Number of partition items (\p p contains \p n+1 items). */
    int n;
    /** Array of partition edges and corresponding bins. */
    t_partition_item* p;
} t_partition;

/*! \internal
 * \brief
 * Internal data structure for the \p insolidangle selection method.
 *
 * Contains the reference points that partially cover a certain region on the
 * surface of the unit sphere.
 * If \p n is -1, the whole region described by the bin is covered.
 *
 * \ingroup module_selection
 */
typedef struct spheresurfacebin
{
    /** Number of points in the array \p x, -1 if whole bin covered. */
    int n;
    /** Number of elements allocated for \p x. */
    int n_alloc;
    /** Array of points that partially cover the bin. */
    rvec* x;
} t_spheresurfacebin;

/*! \internal
 * \brief
 * Data structure for the \p insolidangle selection method.
 *
 * All angle values are in the units of radians.
 *
 * \ingroup module_selection
 */
typedef struct methoddata_insolidangle
{
    /** Center of the solid angle. */
    gmx_ana_pos_t center;
    /** Positions that span the solid angle. */
    gmx_ana_pos_t span;
    /** Cutoff angle. */
    real angcut;
    /** Estimate of the covered fraction. */
    real cfrac;

    /** Cutoff for the cosine (equals cos(angcut)). */
    real distccut;
    /** Bin size to be used as the target bin size when constructing the bins. */
    real targetbinsize;

    /** Number of bins in the \p tbin array. */
    int ntbins;
    /** Size of one bin in the zenith angle direction. */
    real tbinsize;
    /** Array of zenith angle slices. */
    t_partition* tbin;
    /** Number of elements allocated for the \p bin array. */
    int maxbins;
    /** Number of elements used in the \p bin array. */
    int nbins;
    /** Array of individual bins. */
    t_spheresurfacebin* bin;
} t_methoddata_insolidangle;

/*! \brief
 * Allocates data for the \p insolidangle selection method.
 *
 * \param[in]     npar  Not used (should be 3).
 * \param[in,out] param Method parameters (should point to
 *   \ref smparams_insolidangle).
 * \returns Pointer to the allocated data (\ref t_methoddata_insolidangle).
 *
 * Allocates memory for a \ref t_methoddata_insolidangle structure and
 * initializes the parameter as follows:
 *  - \p center defines the value for t_methoddata_insolidangle::center.
 *  - \p span   defines the value for t_methoddata_insolidangle::span.
 *  - \p cutoff defines the value for t_methoddata_insolidangle::angcut.
 */
static void* init_data_insolidangle(int npar, gmx_ana_selparam_t* param);
/*! \brief
 * Initializes the \p insolidangle selection method.
 *
 * \param   top  Not used.
 * \param   npar Not used.
 * \param   param Not used.
 * \param   data Pointer to \ref t_methoddata_insolidangle to initialize.
 * \returns 0 on success, -1 on failure.
 *
 * Converts t_methoddata_insolidangle::angcut to radians and allocates
 * and allocates memory for the bins used during the evaluation.
 */
static void init_insolidangle(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Frees the data allocated for the \p insolidangle selection method. */
static void free_data_insolidangle(void* data);
/*! \brief
 * Initializes the evaluation of the \p insolidangle selection method for a frame.
 *
 * \param[in]  context Evaluation context.
 * \param      data    Should point to a \ref t_methoddata_insolidangle.
 *
 * Creates a lookup structure that enables fast queries of whether a point
 * is within the solid angle or not.
 */
static void init_frame_insolidangle(const gmx::SelMethodEvalContext& context, void* data);
/** Internal helper function for evaluate_insolidangle(). */
static bool accept_insolidangle(rvec x, const t_pbc* pbc, void* data);
/** Evaluates the \p insolidangle selection method. */
static void evaluate_insolidangle(const gmx::SelMethodEvalContext& context,
                                  gmx_ana_pos_t*                   pos,
                                  gmx_ana_selvalue_t*              out,
                                  void*                            data);

/** Calculates the distance between unit vectors. */
static real sph_distc(rvec x1, rvec x2);
/** Does a binary search on a \p t_partition to find a bin for a value. */
static int find_partition_bin(t_partition* p, real value);
/** Finds a bin that corresponds to a location on the unit sphere surface. */
static int find_surface_bin(t_methoddata_insolidangle* surf, rvec x);
/** Clears/initializes the bins on the unit sphere surface. */
static void clear_surface_points(t_methoddata_insolidangle* surf);
/** Frees memory allocated for storing the reference points in the surface bins. */
static void free_surface_points(t_methoddata_insolidangle* surf);
/** Adds a reference point to a given bin. */
static void add_surface_point(t_methoddata_insolidangle* surf, int tbin, int pbin, rvec x);
/** Marks a bin as completely covered. */
static void mark_surface_covered(t_methoddata_insolidangle* surf, int tbin, int pbin);
/** Helper function for store_surface_point() to update a single zenith angle bin. */
static void update_surface_bin(t_methoddata_insolidangle* surf,
                               int                        tbin,
                               real                       phi,
                               real                       pdelta1,
                               real                       pdelta2,
                               real                       pdeltamax,
                               rvec                       x);
/** Adds a single reference point and updates the surface bins. */
static void store_surface_point(t_methoddata_insolidangle* surf, rvec x);
/*! \brief
 * Optimizes the surface bins for faster searching.
 *
 * \param[in,out] surf Surface data structure.
 *
 * Currently, this function does nothing.
 */
static void optimize_surface_points(t_methoddata_insolidangle* surf);
/** Estimates the area covered by the reference cones. */
static real estimate_covered_fraction(t_methoddata_insolidangle* surf);
/** Checks whether a point lies within a solid angle. */
static bool is_surface_covered(t_methoddata_insolidangle* surf, rvec x);

/** Parameters for the \p insolidangle selection method. */
static gmx_ana_selparam_t smparams_insolidangle[] = {
    { "center", { POS_VALUE, 1, { nullptr } }, nullptr, SPAR_DYNAMIC },
    { "span", { POS_VALUE, -1, { nullptr } }, nullptr, SPAR_DYNAMIC | SPAR_VARNUM },
    { "cutoff", { REAL_VALUE, 1, { nullptr } }, nullptr, SPAR_OPTIONAL },
};

/** Help text for the \p insolidangle selection method. */
static const char* const help_insolidangle[] = {
    "::",
    "",
    "  insolidangle center POS span POS_EXPR [cutoff REAL]",
    "",
    "This keyword selects atoms that are within [TT]REAL[tt] degrees",
    "(default=5) of any position in [TT]POS_EXPR[tt] as seen from [TT]POS[tt]",
    "a position expression that evaluates to a single position), i.e., atoms",
    "in the solid angle spanned by the positions in [TT]POS_EXPR[tt] and",
    "centered at [TT]POS[tt].[PAR]",

    "Technically, the solid angle is constructed as a union of small cones",
    "whose tip is at [TT]POS[tt] and the axis goes through a point in",
    "[TT]POS_EXPR[tt]. There is such a cone for each position in",
    "[TT]POS_EXPR[tt], and point is in the solid angle if it lies within any",
    "of these cones. The cutoff determines the width of the cones.",
};

/** Selection method data for the \p insolidangle method. */
gmx_ana_selmethod_t sm_insolidangle = {
    "insolidangle",
    GROUP_VALUE,
    SMETH_DYNAMIC,
    asize(smparams_insolidangle),
    smparams_insolidangle,
    &init_data_insolidangle,
    nullptr,
    &init_insolidangle,
    nullptr,
    &free_data_insolidangle,
    &init_frame_insolidangle,
    nullptr,
    &evaluate_insolidangle,
    { "insolidangle center POS span POS_EXPR [cutoff REAL]",
      "Selecting atoms in a solid angle",
      asize(help_insolidangle),
      help_insolidangle },
};

static void* init_data_insolidangle(int /* npar */, gmx_ana_selparam_t* param)
{
    t_methoddata_insolidangle* data = new t_methoddata_insolidangle();
    data->angcut                    = 5.0;
    data->cfrac                     = 0.0;

    data->distccut      = 0.0;
    data->targetbinsize = 0.0;

    data->ntbins   = 0;
    data->tbinsize = 0.0;
    data->tbin     = nullptr;
    data->maxbins  = 0;
    data->nbins    = 0;
    data->bin      = nullptr;

    param[0].val.u.p = &data->center;
    param[1].val.u.p = &data->span;
    param[2].val.u.r = &data->angcut;
    return data;
}

static void init_insolidangle(const gmx_mtop_t* /* top */,
                              int /* npar */,
                              gmx_ana_selparam_t* /* param */,
                              void* data)
{
    t_methoddata_insolidangle* surf = static_cast<t_methoddata_insolidangle*>(data);
    int                        i, c;

    if (surf->angcut <= 0)
    {
        GMX_THROW(gmx::InvalidInputError("Angle cutoff should be > 0"));
    }

    surf->angcut *= gmx::c_deg2Rad;

    surf->distccut      = -std::cos(surf->angcut);
    surf->targetbinsize = surf->angcut / 2;
    surf->ntbins        = static_cast<int>(M_PI / surf->targetbinsize);
    surf->tbinsize      = (180.0 / surf->ntbins) * gmx::c_deg2Rad;

    snew(surf->tbin, static_cast<int>(M_PI / surf->tbinsize) + 1);
    surf->maxbins = 0;
    for (i = 0; i < surf->ntbins; ++i)
    {
        c = static_cast<int>(std::max(std::sin(surf->tbinsize * i), std::sin(surf->tbinsize * (i + 1)))
                             * M_2PI / surf->targetbinsize)
            + 1;
        snew(surf->tbin[i].p, c + 1);
        surf->maxbins += c;
    }
    surf->nbins = 0;
    snew(surf->bin, surf->maxbins);
}

/*!
 * \param data Data to free (should point to a \ref t_methoddata_insolidangle).
 *
 * Frees the memory allocated for \c t_methoddata_insolidangle::center and
 * \c t_methoddata_insolidangle::span, as well as the memory for the internal
 * bin structure.
 */
static void free_data_insolidangle(void* data)
{
    t_methoddata_insolidangle* d = static_cast<t_methoddata_insolidangle*>(data);
    int                        i;

    if (d->tbin)
    {
        for (i = 0; i < d->ntbins; ++i)
        {
            sfree(d->tbin[i].p);
        }
        sfree(d->tbin);
    }
    free_surface_points(d);
    sfree(d->bin);
    delete d;
}

static void init_frame_insolidangle(const gmx::SelMethodEvalContext& context, void* data)
{
    t_methoddata_insolidangle* d = static_cast<t_methoddata_insolidangle*>(data);
    rvec                       dx;
    int                        i;

    free_surface_points(d);
    clear_surface_points(d);
    for (i = 0; i < d->span.count(); ++i)
    {
        if (context.pbc_)
        {
            pbc_dx(context.pbc_, d->span.x[i], d->center.x[0], dx);
        }
        else
        {
            rvec_sub(d->span.x[i], d->center.x[0], dx);
        }
        unitv(dx, dx);
        store_surface_point(d, dx);
    }
    optimize_surface_points(d);
    d->cfrac = -1;
}

/*!
 * \param[in] x    Test point.
 * \param[in] pbc  PBC data (if NULL, no PBC are used).
 * \param[in] data Pointer to a \c t_methoddata_insolidangle data structure.
 * \returns   true if \p x is within the solid angle, false otherwise.
 */
static bool accept_insolidangle(rvec x, const t_pbc* pbc, void* data)
{
    t_methoddata_insolidangle* d = static_cast<t_methoddata_insolidangle*>(data);
    rvec                       dx;

    if (pbc)
    {
        pbc_dx(pbc, x, d->center.x[0], dx);
    }
    else
    {
        rvec_sub(x, d->center.x[0], dx);
    }
    unitv(dx, dx);
    return is_surface_covered(d, dx);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data should point to a \c t_methoddata_insolidangle.
 *
 * Calculates which atoms in \p g are within the solid angle spanned by
 * \c t_methoddata_insolidangle::span and centered at
 * \c t_methoddata_insolidangle::center, and stores the result in \p out->u.g.
 */
static void evaluate_insolidangle(const gmx::SelMethodEvalContext& context,
                                  gmx_ana_pos_t*                   pos,
                                  gmx_ana_selvalue_t*              out,
                                  void*                            data)
{
    out->u.g->isize = 0;
    for (int b = 0; b < pos->count(); ++b)
    {
        if (accept_insolidangle(pos->x[b], context.pbc_, data))
        {
            gmx_ana_pos_add_to_group(out->u.g, pos, b);
        }
    }
}

/*!
 * \param[in] sel Selection element to query.
 * \returns   true if the covered fraction can be estimated for \p sel with
 *   _gmx_selelem_estimate_coverfrac(), false otherwise.
 */
bool _gmx_selelem_can_estimate_cover(const gmx::SelectionTreeElement& sel)
{
    if (sel.type == SEL_BOOLEAN && sel.u.boolt == BOOL_OR)
    {
        return false;
    }
    bool                             bFound    = false;
    bool                             bDynFound = false;
    gmx::SelectionTreeElementPointer child     = sel.child;
    while (child)
    {
        if (child->type == SEL_EXPRESSION)
        {
            if (child->u.expr.method->name == sm_insolidangle.name)
            {
                if (bFound || bDynFound)
                {
                    return false;
                }
                bFound = true;
            }
            else if (child->u.expr.method && (child->u.expr.method->flags & SMETH_DYNAMIC))
            {
                if (bFound)
                {
                    return false;
                }
                bDynFound = true;
            }
        }
        else if (!_gmx_selelem_can_estimate_cover(*child))
        {
            return false;
        }
        child = child->next;
    }
    return true;
}

/*!
 * \param[in] sel Selection for which the fraction should be calculated.
 * \returns Fraction of angles covered by the selection (between zero and one).
 *
 * The return value is undefined if _gmx_selelem_can_estimate_cover() returns
 * false.
 * Should be called after gmx_ana_evaluate_selections() has been called for the
 * frame.
 */
real _gmx_selelem_estimate_coverfrac(const gmx::SelectionTreeElement& sel)
{
    real cfrac;

    if (sel.type == SEL_EXPRESSION && sel.u.expr.method->name == sm_insolidangle.name)
    {
        t_methoddata_insolidangle* d = static_cast<t_methoddata_insolidangle*>(sel.u.expr.mdata);
        if (d->cfrac < 0)
        {
            d->cfrac = estimate_covered_fraction(d);
        }
        return d->cfrac;
    }
    if (sel.type == SEL_BOOLEAN && sel.u.boolt == BOOL_NOT)
    {
        cfrac = _gmx_selelem_estimate_coverfrac(*sel.child);
        if (cfrac < 1.0)
        {
            return 1 - cfrac;
        }
        return 1;
    }

    /* Here, we assume that the selection is simple enough */
    gmx::SelectionTreeElementPointer child = sel.child;
    while (child)
    {
        cfrac = _gmx_selelem_estimate_coverfrac(*child);
        if (cfrac < 1.0)
        {
            return cfrac;
        }
        child = child->next;
    }
    return 1.0;
}

/*!
 * \param[in] x1  Unit vector 1.
 * \param[in] x2  Unit vector 2.
 * \returns   Minus the dot product of \p x1 and \p x2.
 *
 * This function is used internally to calculate the distance between the
 * unit vectors \p x1 and \p x2 to find out whether \p x2 is within the
 * cone centered at \p x1. Currently, the cosine of the angle is used
 * for efficiency, and the minus is there to make it behave like a normal
 * distance (larger values mean longer distances).
 */
static real sph_distc(rvec x1, rvec x2)
{
    return -iprod(x1, x2);
}

/*!
 * \param[in] p     Partition to search.
 * \param[in] value Value to search for.
 * \returns   The partition index in \p p that contains \p value.
 *
 * If \p value is outside the range of \p p, the first/last index is returned.
 * Otherwise, the return value \c i satisfies \c p->p[i].left<=value and
 * \c p->p[i+1].left>value
 */
static int find_partition_bin(t_partition* p, real value)
{
    int pmin, pmax, pbin;

    /* Binary search the partition */
    pmin = 0;
    pmax = p->n;
    while (pmax > pmin + 1)
    {
        pbin = pmin + (pmax - pmin) / 2;
        if (p->p[pbin].left <= value)
        {
            pmin = pbin;
        }
        else
        {
            pmax = pbin;
        }
    }
    pbin = pmin;
    return pbin;
}

/*!
 * \param[in] surf  Surface data structure to search.
 * \param[in] x     Unit vector to find.
 * \returns   The bin index that contains \p x.
 *
 * The return value is an index to the \p surf->bin array.
 */
static int find_surface_bin(t_methoddata_insolidangle* surf, rvec x)
{
    real theta, phi;
    int  tbin, pbin;

    theta = std::acos(x[ZZ]);
    phi   = std::atan2(x[YY], x[XX]);
    tbin  = static_cast<int>(std::floor(theta / surf->tbinsize));
    if (tbin >= surf->ntbins)
    {
        tbin = surf->ntbins - 1;
    }
    pbin = find_partition_bin(&surf->tbin[tbin], phi);
    return surf->tbin[tbin].p[pbin].bin;
}

/*!
 * \param[in,out] surf Surface data structure.
 *
 * Clears the reference points from the bins and (re)initializes the edges
 * of the azimuthal bins.
 */
static void clear_surface_points(t_methoddata_insolidangle* surf)
{
    int i, j, c;

    surf->nbins = 0;
    for (i = 0; i < surf->ntbins; ++i)
    {
        c = static_cast<int>(std::min(std::sin(surf->tbinsize * i), std::sin(surf->tbinsize * (i + 1)))
                             * M_2PI / surf->targetbinsize)
            + 1;
        if (c <= 0)
        {
            c = 1;
        }
        surf->tbin[i].n = c;
        for (j = 0; j < c; ++j)
        {
            surf->tbin[i].p[j].left  = -M_PI + j * M_2PI / c - 0.0001;
            surf->tbin[i].p[j].bin   = surf->nbins;
            surf->bin[surf->nbins].n = 0;
            surf->nbins++;
        }
        surf->tbin[i].p[c].left = M_PI + 0.0001;
        surf->tbin[i].p[c].bin  = -1;
    }
}

/*!
 * \param[in,out] surf Surface data structure.
 */
static void free_surface_points(t_methoddata_insolidangle* surf)
{
    int i;

    for (i = 0; i < surf->nbins; ++i)
    {
        if (surf->bin[i].x)
        {
            sfree(surf->bin[i].x);
        }
        surf->bin[i].n_alloc = 0;
        surf->bin[i].x       = nullptr;
    }
}

/*!
 * \param[in,out] surf Surface data structure.
 * \param[in]     tbin Bin number in the zenith angle direction.
 * \param[in]     pbin Bin number in the azimuthal angle direction.
 * \param[in]     x    Point to store.
 */
static void add_surface_point(t_methoddata_insolidangle* surf, int tbin, int pbin, rvec x)
{
    int bin;

    bin = surf->tbin[tbin].p[pbin].bin;
    /* Return if bin is already completely covered */
    if (surf->bin[bin].n == -1)
    {
        return;
    }
    /* Allocate more space if necessary */
    if (surf->bin[bin].n == surf->bin[bin].n_alloc)
    {
        surf->bin[bin].n_alloc += 10;
        srenew(surf->bin[bin].x, surf->bin[bin].n_alloc);
    }
    /* Add the point to the bin */
    copy_rvec(x, surf->bin[bin].x[surf->bin[bin].n]);
    ++surf->bin[bin].n;
}

/*!
 * \param[in,out] surf Surface data structure.
 * \param[in]     tbin Bin number in the zenith angle direction.
 * \param[in]     pbin Bin number in the azimuthal angle direction.
 */
static void mark_surface_covered(t_methoddata_insolidangle* surf, int tbin, int pbin)
{
    int bin;

    bin              = surf->tbin[tbin].p[pbin].bin;
    surf->bin[bin].n = -1;
}

/*!
 * \param[in,out] surf      Surface data structure.
 * \param[in]     tbin      Bin number in the zenith angle direction.
 * \param[in]     phi       Azimuthal angle of \p x.
 * \param[in]     pdelta1   Width of the cone at the lower edge of \p tbin.
 * \param[in]     pdelta2   Width of the cone at the uppper edge of \p tbin.
 * \param[in]     pdeltamax Max. width of the cone inside \p tbin.
 * \param[in]     x         Point to store (should have unit length).
 */
static void update_surface_bin(t_methoddata_insolidangle* surf,
                               int                        tbin,
                               real                       phi,
                               real                       pdelta1,
                               real                       pdelta2,
                               real                       pdeltamax,
                               rvec                       x)
{
    real pdelta, phi1, phi2;
    int  pbin1, pbin2, pbiniter, pbin;

    /* Find the edges of the bins affected */
    pdelta = max(max(pdelta1, pdelta2), pdeltamax);
    phi1   = phi - pdelta;
    if (phi1 >= -M_PI)
    {
        pbin  = find_partition_bin(&surf->tbin[tbin], phi1);
        pbin1 = pbin;
    }
    else
    {
        pbin  = find_partition_bin(&surf->tbin[tbin], phi1 + M_2PI);
        pbin1 = pbin - surf->tbin[tbin].n;
    }
    phi2 = phi + pdelta;
    if (phi2 <= M_PI)
    {
        pbin2 = find_partition_bin(&surf->tbin[tbin], phi2);
    }
    else
    {
        pbin2 = find_partition_bin(&surf->tbin[tbin], phi2 - M_2PI);
        pbin2 += surf->tbin[tbin].n;
    }
    ++pbin2;
    if (pbin2 - pbin1 > surf->tbin[tbin].n)
    {
        pbin2 = pbin1 + surf->tbin[tbin].n;
    }
    /* Find the edges of completely covered region */
    pdelta = min(pdelta1, pdelta2);
    phi1   = phi - pdelta;
    if (phi1 < -M_PI)
    {
        phi1 += M_2PI;
    }
    phi2 = phi + pdelta;
    /* Loop over all affected bins */
    for (pbiniter = pbin1; pbiniter != pbin2; ++pbiniter, ++pbin)
    {
        /* Wrap bin around if end reached */
        if (pbin == surf->tbin[tbin].n)
        {
            pbin = 0;
            phi1 -= M_2PI;
            phi2 -= M_2PI;
        }
        /* Check if bin is completely covered and update */
        if (surf->tbin[tbin].p[pbin].left >= phi1 && surf->tbin[tbin].p[pbin + 1].left <= phi2)
        {
            mark_surface_covered(surf, tbin, pbin);
        }
        else
        {
            add_surface_point(surf, tbin, pbin, x);
        }
    }
}

/*!
 * \param[in,out] surf Surface data structure.
 * \param[in]     x    Point to store (should have unit length).
 *
 * Finds all the bins covered by the cone centered at \p x and calls
 * update_surface_bin() to update them.
 */
static void store_surface_point(t_methoddata_insolidangle* surf, rvec x)
{
    real theta, phi;
    real pdeltamax, tmax;
    real theta1, theta2, pdelta1, pdelta2;
    int  tbin;

    theta = std::acos(x[ZZ]);
    phi   = std::atan2(x[YY], x[XX]);
    /* Find the maximum extent in the phi direction */
    if (theta <= surf->angcut)
    {
        pdeltamax = M_PI;
        tmax      = 0;
    }
    else if (theta >= M_PI - surf->angcut)
    {
        pdeltamax = M_PI;
        tmax      = M_PI;
    }
    else
    {
        pdeltamax = std::asin(std::sin(surf->angcut) / std::sin(theta));
        tmax      = std::acos(std::cos(theta) / std::cos(surf->angcut));
    }
    /* Find the first affected bin */
    tbin   = max(static_cast<int>(std::floor((theta - surf->angcut) / surf->tbinsize)), 0);
    theta1 = tbin * surf->tbinsize;
    if (theta1 < theta - surf->angcut)
    {
        pdelta1 = 0;
    }
    else
    {
        pdelta1 = M_PI;
    }
    /* Loop through all affected bins */
    while (tbin < std::ceil((theta + surf->angcut) / surf->tbinsize) && tbin < surf->ntbins)
    {
        /* Calculate the next boundaries */
        theta2 = (tbin + 1) * surf->tbinsize;
        if (theta2 > theta + surf->angcut)
        {
            /* The circle is completely outside the cone */
            pdelta2 = 0;
        }
        else if (theta2 <= -(theta - surf->angcut) || theta2 >= M_2PI - (theta + surf->angcut)
                 || tbin == surf->ntbins - 1)
        {
            /* The circle is completely inside the cone, or we are in the
             * 360 degree bin covering the pole. */
            pdelta2 = M_PI;
        }
        else
        {
            /* TODO: This formula is numerically unstable if theta is very
             * close to the pole.  In practice, it probably does not matter
             * much, but it would be nicer to adjust the theta bin boundaries
             * such that the case above catches this instead of falling through
             * here. */
            pdelta2 = 2
                      * std::asin(std::sqrt((gmx::square(std::sin(surf->angcut / 2))
                                             - gmx::square(std::sin((theta2 - theta) / 2)))
                                            / (std::sin(theta) * std::sin(theta2))));
        }
        /* Update the bin */
        if (tmax >= theta1 && tmax <= theta2)
        {
            update_surface_bin(surf, tbin, phi, pdelta1, pdelta2, pdeltamax, x);
        }
        else
        {
            update_surface_bin(surf, tbin, phi, pdelta1, pdelta2, 0, x);
        }
        /* Next bin */
        theta1  = theta2;
        pdelta1 = pdelta2;
        ++tbin;
    }
}

static void optimize_surface_points(t_methoddata_insolidangle* /* surf */)
{
    /* TODO: Implement */
}

/*!
 * \param[in] surf Surface data structure.
 * \returns   An estimate for the area covered by the reference points.
 */
static real estimate_covered_fraction(t_methoddata_insolidangle* surf)
{
    int  t, p, n;
    real cfrac, tfrac, pfrac;

    cfrac = 0.0;
    for (t = 0; t < surf->ntbins; ++t)
    {
        tfrac = std::cos(t * surf->tbinsize) - std::cos((t + 1) * surf->tbinsize);
        for (p = 0; p < surf->tbin[t].n; ++p)
        {
            pfrac = surf->tbin[t].p[p + 1].left - surf->tbin[t].p[p].left;
            n     = surf->bin[surf->tbin[t].p[p].bin].n;
            if (n == -1) /* Bin completely covered */
            {
                cfrac += tfrac * pfrac;
            }
            else if (n > 0) /* Bin partially covered */
            {
                cfrac += tfrac * pfrac / 2; /* A rough estimate */
            }
        }
    }
    return cfrac / (4 * M_PI);
}

/*!
 * \param[in] surf  Surface data structure to search.
 * \param[in] x     Unit vector to check.
 * \returns   true if \p x is within the solid angle, false otherwise.
 */
static bool is_surface_covered(t_methoddata_insolidangle* surf, rvec x)
{
    int bin, i;

    bin = find_surface_bin(surf, x);
    /* Check for completely covered bin */
    if (surf->bin[bin].n == -1)
    {
        return true;
    }
    /* Check each point that partially covers the bin */
    for (i = 0; i < surf->bin[bin].n; ++i)
    {
        if (sph_distc(x, surf->bin[bin].x[i]) < surf->distccut)
        {
            return true;
        }
    }
    return false;
}
