/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2013- The GROMACS Authors
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
 * Implements functions in swapcoords.h.
 *
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \ingroup module_swap
 */
#include "gmxpre.h"

#include "swapcoords.h"

#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdrunutility/handlerestart.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdrunoptions.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/mdtypes/swaphistory.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class ForceProviders;
class IMDOutputProvider;
class IMdpOptionProvider;
struct MDModulesNotifiers;
} // namespace gmx

static const std::string SwS      = { "SWAP:" }; /**< For output that comes from the swap module */
static const std::string SwSEmpty = { "     " }; /**< Placeholder for multi-line output */
static constexpr gmx::EnumerationArray<Compartment, const char*> CompStr = { "A", "B" }; /**< Compartment name */
static constexpr gmx::EnumerationArray<SwapType, const char*> SwapStr = { "", "X-", "Y-", "Z-" }; /**< Name for the swap types. */
static const char* const DimStr[DIM + 1] = { "X", "Y", "Z", nullptr }; /**< Name for the swap dimension. */

/** Keep track of through which channel the ions have passed */
enum class ChannelHistory : int
{
    None,
    Ch0,
    Ch1,
    Count
};
static constexpr gmx::EnumerationArray<ChannelHistory, const char*> ChannelString = { "none",
                                                                                      "channel0",
                                                                                      "channel1" }; /**< Name for the channels */

/*! \brief Domain identifier.
 *
 * Keeps track of from which compartment the ions came before passing the
 * channel.
 */
enum class Domain : int
{
    Notset,
    A,
    B,
    Count
};
static constexpr gmx::EnumerationArray<Domain, const char*> DomainString = { "not_assigned",
                                                                             "Domain_A",
                                                                             "Domain_B" }; /**< Name for the domains */

namespace gmx
{

extern template LocalAtomSet LocalAtomSetManager::add<void, void>(ArrayRef<const int> globalAtomIndex);

/*! \internal
 * \brief Implement Computational Electrophysiology swapping.
 */
class SwapCoordinates final : public IMDModule
{
    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return nullptr; }
    IMDOutputProvider*  outputProvider() override { return nullptr; }
    void                initForceProviders(ForceProviders* /* forceProviders */) override {}
    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
};

std::unique_ptr<IMDModule> createSwapCoordinatesModule()
{
    return std::make_unique<SwapCoordinates>();
}


} // namespace gmx

/*! \internal \brief
 * Structure containing compartment-specific data.
 */
typedef struct swap_compartment
{
    int nMol;        /**< Number of ion or water molecules detected
                          in this compartment.                          */
    int  nMolBefore; /**< Number of molecules before swapping.          */
    int  nMolReq;    /**< Requested number of molecules in compartment. */
    real nMolAv;     /**< Time-averaged number of molecules matching
                          the compartment conditions.                   */
    int*  nMolPast;  /**< Past molecule counts for time-averaging.      */
    int*  ind;       /**< Indices to collective array of atoms.         */
    real* dist;      /**< Distance of atom to bulk layer, which is
                          normally the center layer of the compartment  */
    int nalloc;      /**< Allocation size for ind array.                */
    int inflow_net;  /**< Net inflow of ions into this compartment.     */
} t_compartment;


/*! \internal \brief
 * This structure contains data needed for the groups involved in swapping:
 * split group 0, split group 1, solvent group, ion groups.
 */
typedef struct swap_group
{
    /*!\brief Construct a swap group given the managed swap atoms.
     *
     * \param[in] atomset Managed indices of atoms that are part of the swap group.
     */
    swap_group(const gmx::LocalAtomSet& atomset);
    char* molname = nullptr;        /**< Name of the group or ion type                         */
    int   apm     = 0;              /**< Number of atoms in each molecule                      */
    gmx::LocalAtomSet atomSet;      /**< The atom indices in the swap group                    */
    rvec*             xc = nullptr; /**< Collective array of group atom positions (size nat)   */
    ivec*   xc_shifts    = nullptr; /**< Current (collective) shifts (size nat)                */
    ivec*   xc_eshifts   = nullptr; /**< Extra shifts since last DD step (size nat)            */
    rvec*   xc_old       = nullptr; /**< Old (collective) positions (size nat)                 */
    real    q            = 0.;      /**< Total charge of one molecule of this group            */
    real*   m            = nullptr; /**< Masses (can be omitted, size apm)                     */
    Domain* comp_from    = nullptr; /**< (Collective) Stores from which compartment this
                                              molecule has come. This way we keep track of
                                              through which channel an ion permeates
                                              (size nMol = nat/apm)                                 */
    Domain* comp_now = nullptr;     /**< In which compartment this ion is now (size nMol)      */
    ChannelHistory* channel_label = nullptr; /**< Which channel was passed at last by this ion?
                                               (size nMol) */
    rvec center; /**< Center of the group; COM if masses are used           */
    gmx::EnumerationArray<Compartment, t_compartment> comp; /**< Distribution of particles of this
                                       group across the two compartments */
    gmx::EnumerationArray<Compartment, real> vacancy; /**< How many molecules need to be swapped in? */
    gmx::EnumerationArray<Channel, int>      fluxfromAtoB; /**< Net flux of ions per channel */
    gmx::EnumerationArray<Channel, int> nCyl; /**< Number of ions residing in a channel         */
    int nCylBoth = 0; /**< Ions assigned to cyl0 and cyl1. Not good.             */
} t_swapgrp;

t_swapgrp::swap_group(const gmx::LocalAtomSet& atomset) : atomSet{ atomset }
{
    center[0] = 0;
    center[1] = 0;
    center[2] = 0;
    for (auto compartment : keysOf(comp))
    {
        comp[compartment]    = {};
        vacancy[compartment] = 0;
    }
    for (auto channel : keysOf(fluxfromAtoB))
    {
        fluxfromAtoB[channel] = 0;
        nCyl[channel]         = 0;
    }
}

/*! \internal \brief
 * Main (private) data structure for the position swapping protocol.
 */
struct t_swap
{
    int                    swapdim;  /**< One of XX, YY, ZZ                               */
    t_pbc*                 pbc;      /**< Needed to make molecules whole.                 */
    FILE*                  fpout;    /**< Output file.                                    */
    int                    ngrp;     /**< Number of t_swapgrp groups                      */
    std::vector<t_swapgrp> group;    /**< Separate groups for channels, solvent, ions     */
    int                    fluxleak; /**< Flux not going through any of the channels.     */
    real                   deltaQ;   /**< The charge imbalance between the compartments.  */
};


/*! \brief Check whether point is in channel.
 *
 * A channel is a cylinder defined by a disc
 * with radius r around its center c. The thickness of the cylinder is
 * d_up - d_down.
 *
 * \code
 *               ^  d_up
 *               |
 *     r         |
 *     <---------c--------->
 *               |
 *               v  d_down
 *
 * \endcode
 *
 * \param[in] point    The position (xyz) under consideration.
 * \param[in] center   The center of the cylinder.
 * \param[in] d_up     The upper extension of the cylinder.
 * \param[in] d_down   The lower extension.
 * \param[in] r_cyl2   Cylinder radius squared.
 * \param[in] pbc      Structure with info about periodic boundary conditions.
 * \param[in] normal   The membrane normal direction is typically 3, i.e. z, but can be x or y also.
 *
 * \returns   Whether the point is inside the defined cylindric channel.
 */
static gmx_bool is_in_channel(rvec point, rvec center, real d_up, real d_down, real r_cyl2, t_pbc* pbc, int normal)
{
    rvec dr;
    int  plane1, plane2; /* Directions tangential to membrane */


    plane1 = (normal + 1) % 3; /* typically 0, i.e. XX */
    plane2 = (normal + 2) % 3; /* typically 1, i.e. YY */

    /* Get the distance vector dr between the point and the center of the cylinder */
    pbc_dx(pbc, point, center, dr); /* This puts center in the origin */

    /* Check vertical direction */
    if ((dr[normal] > d_up) || (dr[normal] < -d_down))
    {
        return FALSE;
    }

    /* Check radial direction */
    if ((dr[plane1] * dr[plane1] + dr[plane2] * dr[plane2]) > r_cyl2)
    {
        return FALSE;
    }

    /* All check passed, this point is in the cylinder */
    return TRUE;
}


/*! \brief Prints output to CompEL output file.
 *
 * Prints to swap output file how many ions are in each compartment,
 * where the centers of the split groups are, and how many ions of each type
 * passed the channels.
 */
static void print_ionlist(t_swap* s, double time, const char comment[])
{
    // Output time
    fprintf(s->fpout, "%12.5e", time);

    // Output number of molecules and difference to reference counts for each
    // compartment and ion type
    for (auto iComp : gmx::EnumerationWrapper<Compartment>{})
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            t_compartment* comp = &s->group[ig].comp[iComp];

            fprintf(s->fpout, "%10d%10.1f%10d", comp->nMol, comp->nMolAv - comp->nMolReq, comp->inflow_net);
        }
    }

    // Output center of split groups
    fprintf(s->fpout,
            "%10g%10g",
            s->group[static_cast<int>(SwapGroupSplittingType::Split0)].center[s->swapdim],
            s->group[static_cast<int>(SwapGroupSplittingType::Split1)].center[s->swapdim]);

    // Output ion flux for each channel and ion type
    for (auto iChan : gmx::EnumerationWrapper<Channel>{})
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            t_swapgrp* g = &s->group[ig];
            fprintf(s->fpout, "%10d", g->fluxfromAtoB[iChan]);
        }
    }

    /* Output the number of molecules that leaked from A to B */
    fprintf(s->fpout, "%10d", s->fluxleak);

    fprintf(s->fpout, "%s\n", comment);
}


/*! \brief Get the center of a group of nat atoms.
 *
 * Since with PBC an atom group might not be whole, use the first atom as the
 * reference atom and determine the center with respect to this reference.
 */
static void get_molecule_center(rvec x[], int nat, const real* weights, rvec center, t_pbc* pbc)
{
    int  i;
    rvec weightedPBCimage;
    real wi, wsum;
    rvec reference, correctPBCimage, dx;


    /* Use the first atom as the reference and put other atoms near that one */
    /* This does not work for large molecules that span > half of the box! */
    copy_rvec(x[0], reference);

    /* Calculate either the weighted center or simply the center of geometry */
    wsum = 0.0;
    clear_rvec(center);
    for (i = 0; i < nat; i++)
    {
        /* PBC distance between position and reference */
        pbc_dx(pbc, x[i], reference, dx);

        /* Add PBC distance to reference */
        rvec_add(reference, dx, correctPBCimage);

        /* Take weight into account */
        if (nullptr == weights)
        {
            wi = 1.0;
        }
        else
        {
            wi = weights[i];
        }
        wsum += wi;
        svmul(wi, correctPBCimage, weightedPBCimage);

        /* Add to center */
        rvec_inc(center, weightedPBCimage);
    }

    /* Normalize */
    svmul(1.0 / wsum, center, center);
}


/*! \brief Return TRUE if position x of ion (or water) is found in the compartment,
 * i.e. between w1 and w2.
 *
 * One can define and additional offset "b" if one wants to exchange ions/water
 * to or from a plane not directly in the middle of w1 and w2. The offset can be
 * in  ]-1.0, ..., +1.0 [.
 * A bulkOffset of 0.0 means 'no offset', so the swap-layer is directly in the
 * middle between w1 and w2. Offsets -1.0 < b <  0.0 will yield swaps nearer to w1,
 * whereas offsets  0.0 < 0 < +1.0 will yield swaps nearer to w2.
 *
 * \code
 *
 * ||--------------+-------------|-------------+------------------------||
 *                w1  ? ? ? ? ? ? ? ? ? ? ?   w2
 * ||--------------+-------------|----b--------+------------------------||
 *                -1            0.0           +1
 *
 * \endcode
 *
 * \param[in]  w1               Position of 'wall' atom 1.
 * \param[in]  w2               Position of 'wall' atom 2.
 * \param[in]  x                Position of the ion or the water molecule under consideration.
 * \param[in]  l                Length of the box, from || to || in the sketch.
 * \param[in]  bulkOffset       Where is the bulk layer "b" to be found between w1 and w2?
 * \param[out] distance_from_b  Distance of x to the bulk layer "b".
 *
 * \returns TRUE if x is between w1 and w2.
 *
 * Also computes the distance of x to the compartment center (the layer that is
 * normally situated in the middle of w1 and w2 that would be considered as having
 * the bulk concentration of ions).
 */
static gmx_bool compartment_contains_atom(real w1, real w2, real x, real l, real bulkOffset, real* distance_from_b)
{
    real m, l_2;
    real width;


    /* First set the origin in the middle of w1 and w2 */
    m = 0.5 * (w1 + w2);
    w1 -= m;
    w2 -= m;
    x -= m;
    width = w2 - w1;

    /* Now choose the PBC image of x that is closest to the origin: */
    l_2 = 0.5 * l;
    while (x > l_2)
    {
        x -= l;
    }
    while (x <= -l_2)
    {
        x += l;
    }

    *distance_from_b = std::fabs(x - bulkOffset * 0.5_real * width);

    /* Return TRUE if we now are in area "????" */
    return (x >= w1) && (x < w2);
}


/*! \brief Updates the time-averaged number of ions in a compartment. */
static void update_time_window(t_compartment* comp, int values, int replace)
{
    real average;
    int  i;


    /* Put in the new value */
    if (replace >= 0)
    {
        comp->nMolPast[replace] = comp->nMol;
    }

    /* Compute the new time-average */
    average = 0.0;
    for (i = 0; i < values; i++)
    {
        average += comp->nMolPast[i];
    }
    average /= values;
    comp->nMolAv = average;
}


/*! \brief Add the atom with collective index ci to the atom list in compartment 'comp'.
 *
 * \param[in]     ci        Index of this ion in the collective xc array.
 * \param[inout]  comp      Compartment to add this atom to.
 * \param[in]     distance  Shortest distance of this atom to the bulk layer,
 *                          from which ion/water pairs are selected for swapping.
 */
static void add_to_list(int ci, t_compartment* comp, real distance)
{
    int nr = comp->nMol;

    if (nr >= comp->nalloc)
    {
        comp->nalloc = over_alloc_dd(nr + 1);
        srenew(comp->ind, comp->nalloc);
        srenew(comp->dist, comp->nalloc);
    }
    comp->ind[nr]  = ci;
    comp->dist[nr] = distance;
    comp->nMol++;
}


/*! \brief Determine the compartment boundaries from the channel centers. */
static void get_compartment_boundaries(Compartment c, t_swap* s, const matrix box, real* left, real* right)
{
    real pos0, pos1;
    real leftpos, rightpos, leftpos_orig;


    if (c >= Compartment::Count)
    {
        gmx_fatal(FARGS, "Compartment out of range");
    }

    pos0 = s->group[static_cast<int>(SwapGroupSplittingType::Split0)].center[s->swapdim];
    pos1 = s->group[static_cast<int>(SwapGroupSplittingType::Split1)].center[s->swapdim];

    if (pos0 < pos1)
    {
        leftpos  = pos0;
        rightpos = pos1;
    }
    else
    {
        leftpos  = pos1;
        rightpos = pos0;
    }

    /* This gets us the other compartment: */
    if (c == Compartment::B)
    {
        leftpos_orig = leftpos;
        leftpos      = rightpos;
        rightpos     = leftpos_orig + box[s->swapdim][s->swapdim];
    }

    *left  = leftpos;
    *right = rightpos;
}


/*! \brief Determine the per-channel ion flux.
 *
 * To determine the flux through the individual channels, we
 * remember the compartment and channel history of each ion. An ion can be
 * either in channel0 or channel1, or in the remaining volume of compartment
 * A or B.
 *
 * \code
 *    +-----------------+
 *    |                 | B
 *    |                 | B compartment
 *    ||||||||||0|||||||| bilayer with channel 0
 *    |                 | A
 *    |                 | A
 *    |                 | A compartment
 *    |                 | A
 *    |||||1||||||||||||| bilayer with channel 1
 *    |                 | B
 *    |                 | B compartment
 *    +-----------------+
 *
 * \endcode
 */
static void detect_flux_per_channel(t_swapgrp*          g,
                                    int                 iAtom,
                                    Compartment         comp,
                                    rvec                atomPosition,
                                    Domain*             comp_now,
                                    Domain*             comp_from,
                                    ChannelHistory*     channel_label,
                                    const t_swapcoords* sc,
                                    t_swap*             s,
                                    real                cyl0_r2,
                                    real                cyl1_r2,
                                    int64_t             step,
                                    gmx_bool            bRerun,
                                    FILE*               fpout)
{
    int      sd, chan_nr;
    gmx_bool in_cyl0, in_cyl1;
    char     buf[STRLEN];


    sd = s->swapdim;

    /* Check whether ion is inside any of the channels */
    in_cyl0 = is_in_channel(atomPosition,
                            s->group[static_cast<int>(SwapGroupSplittingType::Split0)].center,
                            sc->cyl0u,
                            sc->cyl0l,
                            cyl0_r2,
                            s->pbc,
                            sd);
    in_cyl1 = is_in_channel(atomPosition,
                            s->group[static_cast<int>(SwapGroupSplittingType::Split1)].center,
                            sc->cyl1u,
                            sc->cyl1l,
                            cyl1_r2,
                            s->pbc,
                            sd);

    if (in_cyl0 && in_cyl1)
    {
        /* Ion appears to be in both channels. Something is severely wrong! */
        g->nCylBoth++;
        *comp_now      = Domain::Notset;
        *comp_from     = Domain::Notset;
        *channel_label = ChannelHistory::None;
    }
    else if (in_cyl0)
    {
        /* Ion is in channel 0 now */
        *channel_label = ChannelHistory::Ch0;
        *comp_now      = Domain::Notset;
        g->nCyl[Channel::Zero]++;
    }
    else if (in_cyl1)
    {
        /* Ion is in channel 1 now */
        *channel_label = ChannelHistory::Ch1;
        *comp_now      = Domain::Notset;
        g->nCyl[Channel::One]++;
    }
    else
    {
        /* Ion is not in any of the channels, so it must be in domain A or B */
        if (Compartment::A == comp)
        {
            *comp_now = Domain::A;
        }
        else
        {
            *comp_now = Domain::B;
        }
    }

    /* Only take action, if ion is now in domain A or B, and was before
     * in the other domain!
     */
    if (Domain::Notset == *comp_from)
    {
        /* Maybe we can set the domain now */
        *comp_from = *comp_now; /* Could still be Domain::Notset, though */
    }
    else if ((*comp_now != Domain::Notset) /* if in channel */
             && (*comp_from != *comp_now))
    {
        /* Obviously the ion changed its domain.
         * Count this for the channel through which it has passed. */
        switch (*channel_label)
        {
            case ChannelHistory::None:
                ++s->fluxleak;

                fprintf(stderr,
                        " %s Warning! Step %s, ion %d moved from %s to %s\n",
                        SwS.c_str(),
                        gmx_step_str(step, buf),
                        iAtom,
                        DomainString[*comp_from],
                        DomainString[*comp_now]);
                if (bRerun)
                {
                    fprintf(stderr, ", possibly due to a swap in the original simulation.\n");
                }
                else
                {
                    fprintf(stderr,
                            "but did not pass cyl0 or cyl1 as defined in the .mdp file.\n"
                            "Do you have an ion somewhere within the membrane?\n");
                    /* Write this info to the CompEL output file: */
                    fprintf(s->fpout,
                            " # Warning: step %s, ion %d moved from %s to %s (probably through the "
                            "membrane)\n",
                            gmx_step_str(step, buf),
                            iAtom,
                            DomainString[*comp_from],
                            DomainString[*comp_now]);
                }
                break;
            case ChannelHistory::Ch0:
            case ChannelHistory::Ch1:
                if (*channel_label == ChannelHistory::Ch0)
                {
                    chan_nr = 0;
                }
                else
                {
                    chan_nr = 1;
                }

                if (Domain::A == *comp_from)
                {
                    g->fluxfromAtoB[chan_nr]++;
                }
                else
                {
                    g->fluxfromAtoB[chan_nr]--;
                }
                fprintf(fpout, "# Atom nr. %d finished passing %s.\n", iAtom, ChannelString[*channel_label]);
                break;
            default:
                gmx_fatal(FARGS, "%s Unknown channel history entry for ion type '%s'\n", SwS.c_str(), g->molname);
        }

        /* This ion has moved to the _other_ compartment ... */
        *comp_from = *comp_now;
        /* ... and it did not pass any channel yet */
        *channel_label = ChannelHistory::None;
    }
}


/*! \brief Determines which ions or solvent molecules are in compartment A and B */
static void sortMoleculesIntoCompartments(t_swapgrp*          g,
                                          t_commrec*          cr,
                                          const t_swapcoords* sc,
                                          t_swap*             s,
                                          const matrix        box,
                                          int64_t             step,
                                          FILE*               fpout,
                                          gmx_bool            bRerun,
                                          gmx_bool            bIsSolvent)
{
    gmx::EnumerationArray<Compartment, int> nMolNotInComp; /* consistency check */
    real                                    cyl0_r2 = sc->cyl0r * sc->cyl0r;
    real                                    cyl1_r2 = sc->cyl1r * sc->cyl1r;

    /* Get us a counter that cycles in the range of [0 ... sc->nAverage[ */
    int replace = (step / sc->nstswap) % sc->nAverage;

    for (auto comp : gmx::EnumerationWrapper<Compartment>{})
    {
        real left, right;

        /* Get lists of atoms that match criteria for this compartment */
        get_compartment_boundaries(comp, s, box, &left, &right);

        /* First clear the ion molecule lists */
        g->comp[comp].nMol  = 0;
        nMolNotInComp[comp] = 0; /* consistency check */

        /* Loop over the molecules and atoms of this group */
        for (int iMol = 0, iAtom = 0; iAtom < static_cast<int>(g->atomSet.numAtomsGlobal());
             iAtom += g->apm, iMol++)
        {
            real dist;
            int  sd = s->swapdim;

            /* Is this first atom of the molecule in the compartment that we look at? */
            if (compartment_contains_atom(
                        left, right, g->xc[iAtom][sd], box[sd][sd], sc->bulkOffset[comp], &dist))
            {
                /* Add the first atom of this molecule to the list of molecules in this compartment */
                add_to_list(iAtom, &g->comp[comp], dist);

                /* Main also checks for ion groups through which channel each ion has passed */
                if (MAIN(cr) && (g->comp_now != nullptr) && !bIsSolvent)
                {
                    int globalAtomNr = g->atomSet.globalIndex()[iAtom] + 1; /* PDB index starts at 1 ... */
                    detect_flux_per_channel(g,
                                            globalAtomNr,
                                            comp,
                                            g->xc[iAtom],
                                            &g->comp_now[iMol],
                                            &g->comp_from[iMol],
                                            &g->channel_label[iMol],
                                            sc,
                                            s,
                                            cyl0_r2,
                                            cyl1_r2,
                                            step,
                                            bRerun,
                                            fpout);
                }
            }
            else
            {
                nMolNotInComp[comp]++;
            }
        }
        /* Correct the time-averaged number of ions in the compartment */
        if (!bIsSolvent)
        {
            update_time_window(&g->comp[comp], sc->nAverage, replace);
        }
    }

    /* Flux detection warnings */
    if (MAIN(cr) && !bIsSolvent)
    {
        if (g->nCylBoth > 0)
        {
            fprintf(stderr,
                    "\n"
                    "%s Warning: %d atoms were detected as being in both channels! Probably your "
                    "split\n"
                    "%s          cylinder is way too large, or one compartment has collapsed (step "
                    "%" PRId64 ")\n",
                    SwS.c_str(),
                    g->nCylBoth,
                    SwS.c_str(),
                    step);

            fprintf(s->fpout, "Warning: %d atoms were assigned to both channels!\n", g->nCylBoth);

            g->nCylBoth = 0;
        }
    }

    if (bIsSolvent && nullptr != fpout)
    {
        fprintf(fpout,
                "# Solv. molecules in comp.%s: %d   comp.%s: %d\n",
                CompStr[Compartment::A],
                g->comp[Compartment::A].nMol,
                CompStr[Compartment::B],
                g->comp[Compartment::B].nMol);
    }

    /* Consistency checks */
    const auto numMolecules = static_cast<int>(g->atomSet.numAtomsGlobal() / g->apm);
    if (nMolNotInComp[Compartment::A] + nMolNotInComp[Compartment::B] != numMolecules)
    {
        fprintf(stderr,
                "%s Warning: Inconsistency while assigning '%s' molecules to compartments. !inA: "
                "%d, !inB: %d, total molecules %d\n",
                SwS.c_str(),
                g->molname,
                nMolNotInComp[Compartment::A],
                nMolNotInComp[Compartment::B],
                numMolecules);
    }

    int sum = g->comp[Compartment::A].nMol + g->comp[Compartment::B].nMol;
    if (sum != numMolecules)
    {
        fprintf(stderr,
                "%s Warning: %d molecules are in group '%s', but altogether %d have been assigned "
                "to the compartments.\n",
                SwS.c_str(),
                numMolecules,
                g->molname,
                sum);
    }
}


/*! \brief Find out how many group atoms are in the compartments initially */
static void get_initial_ioncounts(const t_inputrec* ir,
                                  t_swap*           s,
                                  const rvec        x[], /* the initial positions */
                                  const matrix      box,
                                  t_commrec*        cr,
                                  gmx_bool          bRerun)
{
    t_swapcoords* sc;
    t_swapgrp*    g;

    sc = ir->swap;

    /* Loop over the user-defined (ion) groups */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g = &s->group[ig];

        /* Copy the initial positions of the atoms in the group
         * to the collective array so that we can compartmentalize */
        for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
        {
            int ind = g->atomSet.globalIndex()[i];
            copy_rvec(x[ind], g->xc[i]);
        }

        /* Set up the compartments and get lists of atoms in each compartment */
        sortMoleculesIntoCompartments(g, cr, sc, s, box, 0, s->fpout, bRerun, FALSE);

        /* Set initial molecule counts if requested (as signaled by "-1" value) */
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            int requested = sc->grp[ig].nmolReq[ic];
            if (requested < 0)
            {
                g->comp[ic].nMolReq = g->comp[ic].nMol;
            }
            else
            {
                g->comp[ic].nMolReq = requested;
            }
        }

        /* Check whether the number of requested molecules adds up to the total number */
        int req = g->comp[Compartment::A].nMolReq + g->comp[Compartment::B].nMolReq;
        int tot = g->comp[Compartment::A].nMol + g->comp[Compartment::B].nMol;

        if ((req != tot))
        {
            gmx_fatal(FARGS,
                      "Mismatch of the number of %s ions summed over both compartments.\n"
                      "You requested a total of %d ions (%d in A and %d in B),\n"
                      "but there are a total of %d ions of this type in the system.\n",
                      g->molname,
                      req,
                      g->comp[Compartment::A].nMolReq,
                      g->comp[Compartment::B].nMolReq,
                      tot);
        }

        /* Initialize time-averaging:
         * Write initial concentrations to all time bins to start with */
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            g->comp[ic].nMolAv = g->comp[ic].nMol;
            for (int i = 0; i < sc->nAverage; i++)
            {
                g->comp[ic].nMolPast[i] = g->comp[ic].nMol;
            }
        }
    }
}


/*! \brief Copy history of ion counts from checkpoint file.
 *
 * When called, the checkpoint file has already been read in. Here we copy
 * over the values from .cpt file to the swap data structure.
 */
static void get_initial_ioncounts_from_cpt(const t_inputrec* ir,
                                           t_swap*           s,
                                           swaphistory_t*    swapstate,
                                           t_commrec*        cr,
                                           gmx_bool          bVerbose)
{
    t_swapcoords*    sc;
    t_swapgrp*       g;
    swapstateIons_t* gs;

    sc = ir->swap;

    if (MAIN(cr))
    {
        /* Copy the past values from the checkpoint values that have been read in already */
        if (bVerbose)
        {
            fprintf(stderr, "%s Copying values from checkpoint\n", SwS.c_str());
        }

        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            g  = &s->group[ig];
            gs = &swapstate->ionType[ig - static_cast<int>(SwapGroupSplittingType::Count)];

            for (auto ic : gmx::EnumerationWrapper<Compartment>{})
            {
                g->comp[ic].nMolReq    = gs->nMolReq[ic];
                g->comp[ic].inflow_net = gs->inflow_net[ic];

                if (bVerbose)
                {
                    fprintf(stderr,
                            "%s ... Influx netto: %d   Requested: %d   Past values: ",
                            SwS.c_str(),
                            g->comp[ic].inflow_net,
                            g->comp[ic].nMolReq);
                }

                for (int j = 0; j < sc->nAverage; j++)
                {
                    g->comp[ic].nMolPast[j] = gs->nMolPast[ic][j];
                    if (bVerbose)
                    {
                        fprintf(stderr, "%d ", g->comp[ic].nMolPast[j]);
                    }
                }
                if (bVerbose)
                {
                    fprintf(stderr, "\n");
                }
            }
        }
    }
}


/*! \brief The main lets all others know about the initial ion counts. */
static void bc_initial_concentrations(t_commrec* cr, t_swapcoords* swap, t_swap* s)
{
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        t_swapgrp* g = &s->group[ig];

        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            gmx_bcast(sizeof(g->comp[ic].nMolReq), &(g->comp[ic].nMolReq), cr->mpi_comm_mygroup);
            gmx_bcast(sizeof(g->comp[ic].nMol), &(g->comp[ic].nMol), cr->mpi_comm_mygroup);
            gmx_bcast(swap->nAverage * sizeof(g->comp[ic].nMolPast[0]), g->comp[ic].nMolPast, cr->mpi_comm_mygroup);
        }
    }
}


/*! \brief Ensure that each atom belongs to at most one of the swap groups. */
static void check_swap_groups(t_swap* s, int nat, gmx_bool bVerbose)
{
    int* nGroup = nullptr; /* This array counts for each atom in the MD system to
                              how many swap groups it belongs (should be 0 or 1!) */
    int ind       = -1;
    int nMultiple = 0; /* Number of atoms belonging to multiple groups */


    if (bVerbose)
    {
        fprintf(stderr, "%s Making sure each atom belongs to at most one of the swap groups.\n", SwS.c_str());
    }

    /* Add one to the group count of atoms belonging to a swap group: */
    snew(nGroup, nat);
    for (int i = 0; i < s->ngrp; i++)
    {
        t_swapgrp* g = &s->group[i];
        for (size_t j = 0; j < g->atomSet.numAtomsGlobal(); j++)
        {
            /* Get the global index of this atom of this group: */
            ind = g->atomSet.globalIndex()[j];
            nGroup[ind]++;
        }
    }
    /* Make sure each atom belongs to at most one of the groups: */
    for (int i = 0; i < nat; i++)
    {
        if (nGroup[i] > 1)
        {
            nMultiple++;
        }
    }
    sfree(nGroup);

    if (nMultiple)
    {
        gmx_fatal(FARGS,
                  "%s Cannot perform swapping since %d atom%s allocated to more than one swap "
                  "index group.\n"
                  "%s Each atom must be allocated to at most one of the split groups, the swap "
                  "groups, or the solvent.\n"
                  "%s Check the .mdp file settings regarding the swap index groups or the index "
                  "groups themselves.\n",
                  SwS.c_str(),
                  nMultiple,
                  (1 == nMultiple) ? " is" : "s are",
                  SwSEmpty.c_str(),
                  SwSEmpty.c_str());
    }
}


/*! \brief Get the number of atoms per molecule for this group.
 *
 * Also ensure that all the molecules in this group have this number of atoms.
 */
static int get_group_apm_check(int igroup, t_swap* s, gmx_bool bVerbose, const gmx_mtop_t& mtop)
{
    t_swapgrp* g   = &s->group[igroup];
    const int* ind = s->group[igroup].atomSet.globalIndex().data();
    int        nat = s->group[igroup].atomSet.numAtomsGlobal();

    /* Determine the number of solvent atoms per solvent molecule from the
     * first solvent atom: */
    int molb = 0;
    mtopGetMolblockIndex(mtop, ind[0], &molb, nullptr, nullptr);
    int apm = mtop.moleculeBlockIndices[molb].numAtomsPerMolecule;

    if (bVerbose)
    {
        fprintf(stderr,
                "%s Checking whether all %s molecules consist of %d atom%s\n",
                SwS.c_str(),
                g->molname,
                apm,
                apm > 1 ? "s" : "");
    }

    /* Check whether this is also true for all other solvent atoms */
    for (int i = 1; i < nat; i++)
    {
        mtopGetMolblockIndex(mtop, ind[i], &molb, nullptr, nullptr);
        if (apm != mtop.moleculeBlockIndices[molb].numAtomsPerMolecule)
        {
            gmx_fatal(FARGS, "Not all molecules of swap group %d consist of %d atoms.", igroup, apm);
        }
    }

    // TODO: check whether charges and masses of each molecule are identical!
    return apm;
}


/*! \brief Print the legend to the swap output file.
 *
 * Also print the initial values of ion counts and position of split groups.
 */
static void print_ionlist_legend(const t_inputrec* ir, t_swap* s, const gmx_output_env_t* oenv)
{
    std::vector<std::string> legend;

    // Number of molecules and difference to reference counts for each
    // compartment and ion type
    for (auto ic : gmx::EnumerationWrapper<Compartment>{})
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            t_swapGroup* g = &ir->swap->grp[ig];
            real         q = s->group[ig].q;

            legend.emplace_back(gmx::formatString(
                    "%s %s ions (charge %s%g)", CompStr[ic], g->molname, q > 0 ? "+" : "", q));
            legend.emplace_back(gmx::formatString(
                    "%s av. mismatch to %d %s ions", CompStr[ic], s->group[ig].comp[ic].nMolReq, g->molname));

            legend.emplace_back(gmx::formatString("%s net %s ion influx", CompStr[ic], g->molname));
        }
    }

    // Center of split groups
    legend.emplace_back(gmx::formatString(
            "%scenter of %s of split group 0",
            SwapStr[ir->eSwapCoords],
            (nullptr != s->group[static_cast<int>(SwapGroupSplittingType::Split0)].m) ? "mass" : "geometry"));
    legend.emplace_back(gmx::formatString(
            "%scenter of %s of split group 1",
            SwapStr[ir->eSwapCoords],
            (nullptr != s->group[static_cast<int>(SwapGroupSplittingType::Split1)].m) ? "mass" : "geometry"));

    // Ion flux for each channel and ion type
    for (auto ic : gmx::EnumerationWrapper<Channel>{})
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            t_swapGroup* g = &ir->swap->grp[ig];
            legend.emplace_back(gmx::formatString(
                    "A->ch%d->B %s permeations", static_cast<int>(ic), g->molname));
        }
    }

    // Number of molecules that leaked from A to B
    legend.emplace_back("leakage");

    xvgrLegend(s->fpout, legend, oenv);

    fprintf(s->fpout,
            "# Instantaneous ion counts and time-averaged differences to requested numbers\n");

    // We add a simple text legend helping to identify the columns with xvgr legend strings
    fprintf(s->fpout, "#  time (ps)");
    for (int i = 0; i < gmx::ssize(legend); i++)
    {
        fprintf(s->fpout, "%10s", gmx::formatString("s%d", i).c_str());
    }
    fprintf(s->fpout, "\n");
    fflush(s->fpout);
}


/*! \brief Initialize channel ion flux detection routine.
 *
 * Initialize arrays that keep track of where the ions come from and where
 * they go.
 */
static void detect_flux_per_channel_init(t_swap* s, swaphistory_t* swapstate, const bool isRestart)
{
    t_swapgrp*       g;
    swapstateIons_t* gs;

    /* All these flux detection routines run on the main only */
    if (swapstate == nullptr)
    {
        return;
    }

    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g  = &s->group[ig];
        gs = &swapstate->ionType[ig - static_cast<int>(SwapGroupSplittingType::Count)];

        /******************************************************/
        /* Channel and domain history for the individual ions */
        /******************************************************/
        if (isRestart) /* set the pointers right */
        {
            g->comp_from     = gs->comp_from;
            g->channel_label = gs->channel_label;
        }
        else /* allocate memory for molecule counts */
        {
            snew(g->comp_from, g->atomSet.numAtomsGlobal() / g->apm);
            gs->comp_from = g->comp_from;
            snew(g->channel_label, g->atomSet.numAtomsGlobal() / g->apm);
            gs->channel_label = g->channel_label;
        }
        snew(g->comp_now, g->atomSet.numAtomsGlobal() / g->apm);

        /* Initialize the channel and domain history counters */
        for (size_t i = 0; i < g->atomSet.numAtomsGlobal() / g->apm; i++)
        {
            g->comp_now[i] = Domain::Notset;
            if (!isRestart)
            {
                g->comp_from[i]     = Domain::Notset;
                g->channel_label[i] = ChannelHistory::None;
            }
        }

        /************************************/
        /* Channel fluxes for both channels */
        /************************************/
        g->nCyl[Channel::Zero] = 0;
        g->nCyl[Channel::One]  = 0;
        g->nCylBoth            = 0;
    }

    if (isRestart)
    {
        fprintf(stderr, "%s Copying channel fluxes from checkpoint file data\n", SwS.c_str());
    }


    // Loop over ion types (and both channels)
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g  = &s->group[ig];
        gs = &swapstate->ionType[ig - static_cast<int>(SwapGroupSplittingType::Count)];

        for (auto ic : gmx::EnumerationWrapper<Channel>{})
        {
            fprintf(stderr,
                    "%s Channel %d flux history for ion type %s (charge %g): ",
                    SwS.c_str(),
                    static_cast<int>(ic),
                    g->molname,
                    g->q);
            if (isRestart)
            {
                g->fluxfromAtoB[ic] = gs->fluxfromAtoB[ic];
            }
            else
            {
                g->fluxfromAtoB[ic] = 0;
            }

            fprintf(stderr, "%d molecule%s", g->fluxfromAtoB[ic], g->fluxfromAtoB[ic] == 1 ? "" : "s");
            fprintf(stderr, "\n");
        }
    }

    /* Set pointers for checkpoint writing */
    swapstate->fluxleak_p = &s->fluxleak;
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g  = &s->group[ig];
        gs = &swapstate->ionType[ig - static_cast<int>(SwapGroupSplittingType::Count)];

        for (auto ic : gmx::EnumerationWrapper<Channel>{})
        {
            gs->fluxfromAtoB_p[ic] = &g->fluxfromAtoB[ic];
        }
    }
}


/*! \brief Outputs the initial structure to PDB file for debugging reasons.
 *
 * Output the starting structure so that in case of multimeric channels
 * the user can check whether we have the correct PBC image for all atoms.
 * If this is not correct, the ion counts per channel will be very likely
 * wrong.
 */
static void outputStartStructureIfWanted(const gmx_mtop_t& mtop, rvec* x, PbcType pbcType, const matrix box)
{
    char* env = getenv("GMX_COMPELDUMP");

    if (env != nullptr)
    {
        fprintf(stderr,
                "\n%s Found env.var. GMX_COMPELDUMP, will output CompEL starting structure made "
                "whole.\n"
                "%s In case of multimeric channels, please check whether they have the correct PBC "
                "representation.\n",
                SwS.c_str(),
                SwSEmpty.c_str());

        write_sto_conf_mtop(
                "CompELAssumedWholeConfiguration.pdb", *mtop.name, mtop, x, nullptr, pbcType, box);
    }
}


/*! \brief Initialize the swapstate structure, used for checkpoint writing.
 *
 * The swapstate struct stores the information we need to make the channels
 * whole again after restarts from a checkpoint file. Here we do the following:
 * a) If we did not start from .cpt, we prepare the struct for proper .cpt writing,
 * b) if we did start from .cpt, we copy over the last whole structures from .cpt,
 * c) in any case, for subsequent checkpoint writing, we set the pointers in
 * swapstate to the x_old arrays, which contain the correct PBC representation of
 * multimeric channels at the last time step.
 */
static void init_swapstate(swaphistory_t*    swapstate,
                           t_swapcoords*     sc,
                           t_swap*           s,
                           const gmx_mtop_t& mtop,
                           const rvec*       x, /* the initial positions */
                           const matrix      box,
                           const t_inputrec* ir)
{
    rvec*      x_pbc = nullptr; /* positions of the whole MD system with molecules made whole */
    t_swapgrp* g;


    /* We always need the last whole positions such that
     * in the next time step we can make the channels whole again in PBC */
    if (swapstate->bFromCpt)
    {
        /* Copy the last whole positions of each channel from .cpt */
        g = &(s->group[static_cast<int>(SwapGroupSplittingType::Split0)]);
        for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
        {
            copy_rvec(swapstate->xc_old_whole[Channel::Zero][i], g->xc_old[i]);
        }
        g = &(s->group[static_cast<int>(SwapGroupSplittingType::Split1)]);
        for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
        {
            copy_rvec(swapstate->xc_old_whole[Channel::One][i], g->xc_old[i]);
        }
    }
    else
    {
        swapstate->eSwapCoords = ir->eSwapCoords;

        /* Set the number of ion types and allocate memory for checkpointing */
        swapstate->nIonTypes = s->ngrp - static_cast<int>(SwapGroupSplittingType::Count);
        snew(swapstate->ionType, swapstate->nIonTypes);

        /* Store the total number of ions of each type in the swapstateIons
         * structure that is accessible during checkpoint writing */
        for (int ii = 0; ii < swapstate->nIonTypes; ii++)
        {
            swapstateIons_t* gs = &swapstate->ionType[ii];
            gs->nMol            = sc->grp[ii + static_cast<int>(SwapGroupSplittingType::Count)].nat;
        }

        /* Extract the initial split group positions. */

        /* Remove pbc, make molecule whole. */
        snew(x_pbc, mtop.natoms);
        copy_rvecn(x, x_pbc, 0, mtop.natoms);

        /* This can only make individual molecules whole, not multimers */
        do_pbc_mtop(ir->pbcType, box, &mtop, x_pbc);

        /* Output the starting structure? */
        outputStartStructureIfWanted(mtop, x_pbc, ir->pbcType, box);

        /* If this is the first run (i.e. no checkpoint present) we assume
         * that the starting positions give us the correct PBC representation */
        for (int ig = static_cast<int>(SwapGroupSplittingType::Split0);
             ig <= static_cast<int>(SwapGroupSplittingType::Split1);
             ig++)
        {
            g = &(s->group[ig]);
            for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
            {
                copy_rvec(x_pbc[g->atomSet.globalIndex()[i]], g->xc_old[i]);
            }
        }
        sfree(x_pbc);

        /* Prepare swapstate arrays for later checkpoint writing */
        swapstate->nat[Channel::Zero] =
                s->group[static_cast<int>(SwapGroupSplittingType::Split0)].atomSet.numAtomsGlobal();
        swapstate->nat[Channel::One] =
                s->group[static_cast<int>(SwapGroupSplittingType::Split1)].atomSet.numAtomsGlobal();
    }

    /* For subsequent checkpoint writing, set the swapstate pointers to the xc_old
     * arrays that get updated at every swapping step */
    swapstate->xc_old_whole_p[Channel::Zero] =
            &s->group[static_cast<int>(SwapGroupSplittingType::Split0)].xc_old;
    swapstate->xc_old_whole_p[Channel::One] =
            &s->group[static_cast<int>(SwapGroupSplittingType::Split1)].xc_old;
}

/*! \brief Determine the total charge imbalance resulting from the swap groups */
static real getRequestedChargeImbalance(t_swap* s)
{
    int                                      ig;
    real                                     DeltaQ = 0.0;
    t_swapgrp*                               g;
    real                                     particle_charge;
    gmx::EnumerationArray<Compartment, real> particle_number;

    for (ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g = &s->group[ig];

        particle_charge                 = g->q;
        particle_number[Compartment::A] = g->comp[Compartment::A].nMolReq;
        particle_number[Compartment::B] = g->comp[Compartment::B].nMolReq;

        DeltaQ += particle_charge * (particle_number[Compartment::A] - particle_number[Compartment::B]);
    }

    return DeltaQ;
}


/*! \brief Sorts anions and cations into two separate groups
 *
 * This routine should be called for the 'anions' and 'cations' group,
 * of which the indices were lumped together in the older version of the code.
 */
static void copyIndicesToGroup(const int* indIons, int nIons, t_swapGroup* g, t_commrec* cr)
{
    g->nat = nIons;

    /* If explicit ion counts were requested in the .mdp file
     * (by setting positive values for the number of ions),
     * we can make an additional consistency check here */
    if ((g->nmolReq[Compartment::A] < 0) && (g->nmolReq[Compartment::B] < 0))
    {
        if (g->nat != (g->nmolReq[Compartment::A] + g->nmolReq[Compartment::B]))
        {
            gmx_fatal_collective(FARGS,
                                 cr->mpi_comm_mysim,
                                 MAIN(cr),
                                 "%s Inconsistency while importing swap-related data from an old "
                                 "input file version.\n"
                                 "%s The requested ion counts in compartments A (%d) and B (%d)\n"
                                 "%s do not add up to the number of ions (%d) of this type for the "
                                 "group '%s'.\n",
                                 SwS.c_str(),
                                 SwSEmpty.c_str(),
                                 g->nmolReq[Compartment::A],
                                 g->nmolReq[Compartment::B],
                                 SwSEmpty.c_str(),
                                 g->nat,
                                 g->molname);
        }
    }

    srenew(g->ind, g->nat);
    for (int i = 0; i < g->nat; i++)
    {
        g->ind[i] = indIons[i];
    }
}


/*! \brief Converts old .tpr file CompEL contents to new data layout.
 *
 *  If we have read an old .tpr file (tpxv <= tpxv_CompElPolyatomicIonsAndMultipleIonTypes),
 * anions and cations are stored together in group #3. In the new
 * format we store each ion type in a separate group.
 * The 'classic' groups are:
 * #0 split group 0  - OK
 * #1 split group 1  - OK
 * #2 solvent        - OK
 * #3 anions         - contains also cations, needs to be converted
 * #4 cations        - empty before conversion
 *
 */
static void convertOldToNewGroupFormat(t_swapcoords* sc, const gmx_mtop_t& mtop, gmx_bool bVerbose, t_commrec* cr)
{
    t_swapGroup* g = &sc->grp[3];

    /* Loop through the atom indices of group #3 (anions) and put all indices
     * that belong to cations into the cation group.
     */
    int  nAnions    = 0;
    int  nCations   = 0;
    int* indAnions  = nullptr;
    int* indCations = nullptr;
    snew(indAnions, g->nat);
    snew(indCations, g->nat);

    int molb = 0;
    for (int i = 0; i < g->nat; i++)
    {
        const t_atom& atom = mtopGetAtomParameters(mtop, g->ind[i], &molb);
        if (atom.q < 0)
        {
            // This is an anion, add it to the list of anions
            indAnions[nAnions++] = g->ind[i];
        }
        else
        {
            // This is a cation, add it to the list of cations
            indCations[nCations++] = g->ind[i];
        }
    }

    if (bVerbose)
    {
        fprintf(stdout,
                "%s Sorted %d ions into separate groups of %d anions and %d cations.\n",
                SwS.c_str(),
                g->nat,
                nAnions,
                nCations);
    }


    /* Now we have the correct lists of anions and cations.
     * Copy it to the right groups.
     */
    copyIndicesToGroup(indAnions, nAnions, g, cr);
    g = &sc->grp[4];
    copyIndicesToGroup(indCations, nCations, g, cr);
    sfree(indAnions);
    sfree(indCations);
}


/*! \brief Returns TRUE if we started from an old .tpr
 *
 * Then we need to re-sort anions and cations into separate groups */
static gmx_bool bConvertFromOldTpr(t_swapcoords* sc)
{
    // If the last group has no atoms it means we need to convert!
    return (sc->ngrp >= 5) && (0 == sc->grp[4].nat);
}


t_swap* init_swapcoords(FILE*                       fplog,
                        const t_inputrec*           ir,
                        const char*                 fn,
                        const gmx_mtop_t&           mtop,
                        const t_state*              globalState,
                        ObservablesHistory*         oh,
                        t_commrec*                  cr,
                        gmx::LocalAtomSetManager*   atomSets,
                        const gmx_output_env_t*     oenv,
                        const gmx::MdrunOptions&    mdrunOptions,
                        const gmx::StartingBehavior startingBehavior)
{
    t_swapgrp*       g;
    swapstateIons_t* gs;
    swaphistory_t*   swapstate = nullptr;

    if ((PAR(cr)) && !haveDDAtomOrdering(*cr))
    {
        gmx_fatal(FARGS, "Position swapping is only implemented for domain decomposition!");
    }

    auto* sc = ir->swap;
    auto* s  = new t_swap();

    if (mdrunOptions.rerun)
    {
        if (PAR(cr))
        {
            gmx_fatal(FARGS,
                      "%s This module does not support reruns in parallel\nPlease request a serial "
                      "run with -nt 1 / -np 1\n",
                      SwS.c_str());
        }

        fprintf(stderr, "%s Rerun - using every available frame\n", SwS.c_str());
        sc->nstswap  = 1;
        sc->nAverage = 1; /* averaging makes no sense for reruns */
    }

    if (MAIN(cr) && startingBehavior == gmx::StartingBehavior::NewSimulation)
    {
        fprintf(fplog, "\nInitializing ion/water position exchanges\n");
        please_cite(fplog, "Kutzner2011b");
    }

    switch (ir->eSwapCoords)
    {
        case SwapType::X: s->swapdim = XX; break;
        case SwapType::Y: s->swapdim = YY; break;
        case SwapType::Z: s->swapdim = ZZ; break;
        default: s->swapdim = -1; break;
    }

    const gmx_bool bVerbose = mdrunOptions.verbose;

    // For compatibility with old .tpr files
    if (bConvertFromOldTpr(sc))
    {
        convertOldToNewGroupFormat(sc, mtop, bVerbose && MAIN(cr), cr);
    }

    /* Copy some data and pointers to the group structures for convenience */
    /* Number of atoms in the group */
    s->ngrp = sc->ngrp;
    for (int i = 0; i < s->ngrp; i++)
    {
        s->group.emplace_back(atomSets->add(
                gmx::ArrayRef<const int>(sc->grp[i].ind, sc->grp[i].ind + sc->grp[i].nat)));
        s->group[i].molname = sc->grp[i].molname;
    }

    /* Check for overlapping atoms */
    check_swap_groups(s, mtop.natoms, bVerbose && MAIN(cr));

    /* Allocate space for the collective arrays for all groups */
    /* For the collective position array */
    for (int i = 0; i < s->ngrp; i++)
    {
        g = &s->group[i];
        snew(g->xc, g->atomSet.numAtomsGlobal());

        /* For the split groups (the channels) we need some extra memory to
         * be able to make the molecules whole even if they span more than
         * half of the box size. */
        if ((i == static_cast<int>(SwapGroupSplittingType::Split0))
            || (i == static_cast<int>(SwapGroupSplittingType::Split1)))
        {
            snew(g->xc_shifts, g->atomSet.numAtomsGlobal());
            snew(g->xc_eshifts, g->atomSet.numAtomsGlobal());
            snew(g->xc_old, g->atomSet.numAtomsGlobal());
        }
    }

    if (MAIN(cr))
    {
        if (oh->swapHistory == nullptr)
        {
            oh->swapHistory = std::make_unique<swaphistory_t>(swaphistory_t{});
        }
        swapstate = oh->swapHistory.get();

        init_swapstate(swapstate, sc, s, mtop, globalState->x.rvec_array(), globalState->box, ir);
    }

    /* After init_swapstate we have a set of (old) whole positions for our
     * channels. Now transfer that to all nodes */
    if (PAR(cr))
    {
        for (int ig = static_cast<int>(SwapGroupSplittingType::Split0);
             ig <= static_cast<int>(SwapGroupSplittingType::Split1);
             ig++)
        {
            g = &(s->group[ig]);
            gmx_bcast((g->atomSet.numAtomsGlobal()) * sizeof((g->xc_old)[0]), g->xc_old, cr->mpi_comm_mygroup);
        }
    }

    /* Make sure that all molecules in the solvent and ion groups contain the
     * same number of atoms each */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Solvent); ig < s->ngrp; ig++)
    {
        real charge;

        g      = &(s->group[ig]);
        g->apm = get_group_apm_check(ig, s, MAIN(cr) && bVerbose, mtop);

        /* Since all molecules of a group are equal, we only need enough space
         * to determine properties of a single molecule at at time */
        snew(g->m, g->apm); /* For the center of mass */
        charge   = 0;       /* To determine the charge imbalance */
        int molb = 0;
        for (int j = 0; j < g->apm; j++)
        {
            const t_atom& atom = mtopGetAtomParameters(mtop, g->atomSet.globalIndex()[j], &molb);
            g->m[j]            = atom.m;
            charge += atom.q;
        }
        /* Total charge of one molecule of this group: */
        g->q = charge;
    }


    /* Need mass-weighted center of split group? */
    for (int j = static_cast<int>(SwapGroupSplittingType::Split0);
         j <= static_cast<int>(SwapGroupSplittingType::Split1);
         j++)
    {
        g = &(s->group[j]);
        if (sc->massw_split[j])
        {
            /* Save the split group masses if mass-weighting is requested */
            snew(g->m, g->atomSet.numAtomsGlobal());
            int molb = 0;
            for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
            {
                g->m[i] = mtopGetAtomMass(mtop, g->atomSet.globalIndex()[i], &molb);
            }
        }
    }

    /* Make a t_pbc struct on all nodes so that the molecules
     * chosen for an exchange can be made whole. */
    snew(s->pbc, 1);

    bool restartWithAppending = (startingBehavior == gmx::StartingBehavior::RestartWithAppending);
    if (MAIN(cr))
    {
        if (bVerbose)
        {
            fprintf(stderr,
                    "%s Opening output file %s%s\n",
                    SwS.c_str(),
                    fn,
                    restartWithAppending ? " for appending" : "");
        }

        s->fpout = gmx_fio_fopen(fn, restartWithAppending ? "a" : "w");

        if (!restartWithAppending)
        {
            xvgr_header(s->fpout, "Molecule counts", "Time (ps)", "counts", exvggtXNY, oenv);

            for (int ig = 0; ig < s->ngrp; ig++)
            {
                auto enumValue = static_cast<SwapGroupSplittingType>(ig);
                g              = &(s->group[ig]);
                fprintf(s->fpout,
                        "# %s group '%s' contains %d atom%s",
                        ig < static_cast<int>(SwapGroupSplittingType::Count) ? enumValueToString(enumValue)
                                                                             : "Ion",
                        g->molname,
                        static_cast<int>(g->atomSet.numAtomsGlobal()),
                        (g->atomSet.numAtomsGlobal() > 1) ? "s" : "");
                if (!(SwapGroupSplittingType::Split0 == enumValue
                      || SwapGroupSplittingType::Split1 == enumValue))
                {
                    fprintf(s->fpout,
                            " with %d atom%s in each molecule of charge %g",
                            g->apm,
                            (g->apm > 1) ? "s" : "",
                            g->q);
                }
                fprintf(s->fpout, ".\n");
            }

            fprintf(s->fpout, "#\n# Initial positions of split groups:\n");
        }

        for (int j = static_cast<int>(SwapGroupSplittingType::Split0);
             j <= static_cast<int>(SwapGroupSplittingType::Split1);
             j++)
        {
            auto enumValue = static_cast<SwapGroupSplittingType>(j);
            g              = &(s->group[j]);
            for (size_t i = 0; i < g->atomSet.numAtomsGlobal(); i++)
            {
                copy_rvec(globalState->x[sc->grp[j].ind[i]], g->xc[i]);
            }
            /* xc has the correct PBC representation for the two channels, so we do
             * not need to correct for that */
            get_center(g->xc, g->m, g->atomSet.numAtomsGlobal(), g->center);
            if (!restartWithAppending)
            {
                fprintf(s->fpout,
                        "# %s group %s-center %5f nm\n",
                        enumValueToString(enumValue),
                        DimStr[s->swapdim],
                        g->center[s->swapdim]);
            }
        }

        if (!restartWithAppending)
        {
            if ((0 != sc->bulkOffset[Compartment::A]) || (0 != sc->bulkOffset[Compartment::B]))
            {
                fprintf(s->fpout, "#\n");
                fprintf(s->fpout,
                        "# You provided an offset for the position of the bulk layer(s).\n");
                fprintf(s->fpout,
                        "# That means the layers to/from which ions and water molecules are "
                        "swapped\n");
                fprintf(s->fpout,
                        "# are not midway (= at 0.0) between the compartment-defining layers (at "
                        "+/- 1.0).\n");
                fprintf(s->fpout, "# bulk-offsetA = %g\n", sc->bulkOffset[Compartment::A]);
                fprintf(s->fpout, "# bulk-offsetB = %g\n", sc->bulkOffset[Compartment::B]);
            }

            fprintf(s->fpout, "#\n");
            fprintf(s->fpout,
                    "# Split0 cylinder radius %f nm, up %f nm, down %f nm\n",
                    sc->cyl0r,
                    sc->cyl0u,
                    sc->cyl0l);
            fprintf(s->fpout,
                    "# Split1 cylinder radius %f nm, up %f nm, down %f nm\n",
                    sc->cyl1r,
                    sc->cyl1u,
                    sc->cyl1l);

            fprintf(s->fpout, "#\n");
            if (!mdrunOptions.rerun)
            {
                fprintf(s->fpout,
                        "# Coupling constant (number of swap attempt steps to average over): %d  "
                        "(translates to %f ps).\n",
                        sc->nAverage,
                        sc->nAverage * sc->nstswap * ir->delta_t);
                fprintf(s->fpout, "# Threshold is %f\n", sc->threshold);
                fprintf(s->fpout, "#\n");
                fprintf(s->fpout,
                        "# Remarks about which atoms passed which channel use global atoms numbers "
                        "starting at one.\n");
            }
        }
    }
    else
    {
        s->fpout = nullptr;
    }

    /* Allocate memory to remember the past particle counts for time averaging */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        g = &(s->group[ig]);
        for (auto ic : keysOf(g->comp))
        {
            snew(g->comp[ic].nMolPast, sc->nAverage);
        }
    }

    /* Get the initial particle concentrations and let the other nodes know */
    if (MAIN(cr))
    {
        if (startingBehavior != gmx::StartingBehavior::NewSimulation)
        {
            get_initial_ioncounts_from_cpt(ir, s, swapstate, cr, bVerbose);
        }
        else
        {
            fprintf(stderr, "%s Determining initial numbers of ions per compartment.\n", SwS.c_str());
            get_initial_ioncounts(
                    ir, s, globalState->x.rvec_array(), globalState->box, cr, mdrunOptions.rerun);
        }

        /* Prepare (further) checkpoint writes ... */
        if (startingBehavior != gmx::StartingBehavior::NewSimulation)
        {
            /* Consistency check */
            if (swapstate->nAverage != sc->nAverage)
            {
                gmx_fatal(FARGS,
                          "%s Ion count averaging steps mismatch! checkpoint: %d, tpr: %d",
                          SwS.c_str(),
                          swapstate->nAverage,
                          sc->nAverage);
            }
        }
        else
        {
            swapstate->nAverage = sc->nAverage;
        }
        fprintf(stderr, "%s Setting pointers for checkpoint writing\n", SwS.c_str());
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
            {
                g  = &s->group[ig];
                gs = &swapstate->ionType[ig - static_cast<int>(SwapGroupSplittingType::Count)];

                gs->nMolReq_p[ic]    = &(g->comp[ic].nMolReq);
                gs->nMolPast_p[ic]   = &(g->comp[ic].nMolPast[0]);
                gs->inflow_net_p[ic] = &(g->comp[ic].inflow_net);
            }
        }

        /* Determine the total charge imbalance */
        s->deltaQ = getRequestedChargeImbalance(s);

        if (bVerbose)
        {
            fprintf(stderr, "%s Requested charge imbalance is Q(A) - Q(B) = %g e.\n", SwS.c_str(), s->deltaQ);
        }
        if (!restartWithAppending)
        {
            fprintf(s->fpout, "# Requested charge imbalance is Q(A)-Q(B) = %g e.\n", s->deltaQ);
        }
    }

    if (PAR(cr))
    {
        bc_initial_concentrations(cr, ir->swap, s);
    }

    /* Update the time-averaged number of molecules for all groups and compartments */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < sc->ngrp; ig++)
    {
        g = &s->group[ig];
        for (auto ic : keysOf(g->comp))
        {
            update_time_window(&g->comp[ic], sc->nAverage, -1);
        }
    }

    /* Initialize arrays that keep track of through which channel the ions go */
    detect_flux_per_channel_init(s, swapstate, startingBehavior != gmx::StartingBehavior::NewSimulation);

    /* We need to print the legend if we open this file for the first time. */
    if (MAIN(cr) && !restartWithAppending)
    {
        print_ionlist_legend(ir, s, oenv);
    }
    return s;
}


void finish_swapcoords(t_swap* s)
{
    if (s == nullptr)
    {
        return;
    }
    if (s->fpout)
    {
        // Close the swap output file
        gmx_fio_fclose(s->fpout);
    }
}

/*! \brief Do we need to swap a molecule in any of the ion groups with a water molecule at this step?
 *
 * From the requested and average molecule counts we determine whether a swap is needed
 * at this time step.
 */
static gmx_bool need_swap(const t_swapcoords* sc, t_swap* s)
{
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < sc->ngrp; ig++)
    {
        t_swapgrp* g = &s->group[ig];

        for (auto ic : keysOf(g->comp))
        {
            if (g->comp[ic].nMolReq - g->comp[ic].nMolAv >= sc->threshold)
            {
                return TRUE;
            }
        }
    }
    return FALSE;
}


/*! \brief Return the index of an atom or molecule suitable for swapping.
 *
 * Returns the index of an atom that is far off the compartment boundaries,
 * that is near to the bulk layer to/from which the swaps take place.
 * Other atoms of the molecule (if any) will directly follow the returned index.
 *
 * \param[in] comp    Structure containing compartment-specific data.
 * \param[in] molname Name of the molecule.
 *
 * \returns Index of the first atom of the molecule chosen for a position exchange.
 */
static int get_index_of_distant_atom(t_compartment* comp, const char molname[])
{
    int  ibest = -1;
    real d     = GMX_REAL_MAX;


    /* comp->nat contains the original number of atoms in this compartment
     * prior to doing any swaps. Some of these atoms may already have been
     * swapped out, but then they are marked with a distance of GMX_REAL_MAX
     */
    for (int iMol = 0; iMol < comp->nMolBefore; iMol++)
    {
        if (comp->dist[iMol] < d)
        {
            ibest = iMol;
            d     = comp->dist[ibest];
        }
    }

    if (ibest < 0)
    {
        gmx_fatal(FARGS,
                  "Could not get index of %s atom. Compartment contains %d %s molecules before "
                  "swaps.",
                  molname,
                  comp->nMolBefore,
                  molname);
    }

    /* Set the distance of this index to infinity such that it won't get selected again in
     * this time step
     */
    comp->dist[ibest] = GMX_REAL_MAX;

    return comp->ind[ibest];
}


/*! \brief Swaps centers of mass and makes molecule whole if broken */
static void translate_positions(rvec* x, int apm, rvec old_com, rvec new_com, t_pbc* pbc)
{
    int  i;
    rvec reference, dx, correctPBCimage;


    /* Use the first atom as the reference for PBC */
    copy_rvec(x[0], reference);

    for (i = 0; i < apm; i++)
    {
        /* PBC distance between position and reference */
        pbc_dx(pbc, x[i], reference, dx);

        /* Add PBC distance to reference */
        rvec_add(reference, dx, correctPBCimage);

        /* Subtract old_com from correct image and add new_com */
        rvec_dec(correctPBCimage, old_com);
        rvec_inc(correctPBCimage, new_com);

        copy_rvec(correctPBCimage, x[i]);
    }
}


/*! \brief Write back the modified local positions from the collective array to the official positions. */
static void apply_modified_positions(swap_group* g, rvec x[])
{
    auto collectiveIndex = g->atomSet.collectiveIndex().begin();
    for (const auto localIndex : g->atomSet.localIndex())
    {
        /* Copy the possibly modified position */
        copy_rvec(g->xc[*collectiveIndex], x[localIndex]);
        ++collectiveIndex;
    }
}


gmx_bool do_swapcoords(t_commrec*        cr,
                       int64_t           step,
                       double            t,
                       const t_inputrec* ir,
                       t_swap*           s,
                       gmx_wallcycle*    wcycle,
                       rvec              x[],
                       matrix            box,
                       gmx_bool          bVerbose,
                       gmx_bool          bRerun)
{
    const t_swapcoords* sc    = ir->swap;
    gmx_bool            bSwap = FALSE;
    t_swapgrp*          gsol;
    int                 isol, iion;
    rvec                com_solvent, com_particle; /* solvent and swap molecule's center of mass */


    wallcycle_start(wcycle, WallCycleCounter::Swap);

    set_pbc(s->pbc, ir->pbcType, box);

    /* Assemble the positions of the split groups, i.e. the channels.
     * Here we also pass a shifts array to communicate_group_positions(), so that it can make
     * the molecules whole even in cases where they span more than half of the box in
     * any dimension */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Split0);
         ig <= static_cast<int>(SwapGroupSplittingType::Split1);
         ig++)
    {
        t_swapgrp* g = &(s->group[ig]);
        communicate_group_positions(cr,
                                    g->xc,
                                    g->xc_shifts,
                                    g->xc_eshifts,
                                    TRUE,
                                    x,
                                    g->atomSet.numAtomsGlobal(),
                                    g->atomSet.numAtomsLocal(),
                                    g->atomSet.localIndex().data(),
                                    g->atomSet.collectiveIndex().data(),
                                    g->xc_old,
                                    box);

        get_center(g->xc, g->m, g->atomSet.numAtomsGlobal(), g->center); /* center of split groups == channels */
    }

    /* Assemble the positions of the ions (ig = 3, 4, ...). These molecules should
     * be small and we can always make them whole with a simple distance check.
     * Therefore we pass NULL as third argument. */
    for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
    {
        t_swapgrp* g = &(s->group[ig]);
        communicate_group_positions(cr,
                                    g->xc,
                                    nullptr,
                                    nullptr,
                                    FALSE,
                                    x,
                                    g->atomSet.numAtomsGlobal(),
                                    g->atomSet.numAtomsLocal(),
                                    g->atomSet.localIndex().data(),
                                    g->atomSet.collectiveIndex().data(),
                                    nullptr,
                                    nullptr);

        /* Determine how many ions of this type each compartment contains */
        sortMoleculesIntoCompartments(g, cr, sc, s, box, step, s->fpout, bRerun, FALSE);
    }

    /* Output how many ions are in the compartments */
    if (MAIN(cr))
    {
        print_ionlist(s, t, "");
    }

    /* If we are doing a rerun, we are finished here, since we cannot perform
     * swaps anyway */
    if (bRerun)
    {
        return FALSE;
    }

    /* Do we have to perform a swap? */
    bSwap = need_swap(sc, s);
    if (bSwap)
    {
        /* Since we here know that we have to perform ion/water position exchanges,
         * we now assemble the solvent positions */
        t_swapgrp* g = &(s->group[static_cast<int>(SwapGroupSplittingType::Solvent)]);
        communicate_group_positions(cr,
                                    g->xc,
                                    nullptr,
                                    nullptr,
                                    FALSE,
                                    x,
                                    g->atomSet.numAtomsGlobal(),
                                    g->atomSet.numAtomsLocal(),
                                    g->atomSet.localIndex().data(),
                                    g->atomSet.collectiveIndex().data(),
                                    nullptr,
                                    nullptr);

        /* Determine how many molecules of solvent each compartment contains */
        sortMoleculesIntoCompartments(g, cr, sc, s, box, step, s->fpout, bRerun, TRUE);

        /* Save number of solvent molecules per compartment prior to any swaps */
        g->comp[Compartment::A].nMolBefore = g->comp[Compartment::A].nMol;
        g->comp[Compartment::B].nMolBefore = g->comp[Compartment::B].nMol;

        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            t_swapgrp* g = &(s->group[ig]);

            for (auto ic : keysOf(g->comp))
            {
                /* Determine in which compartment ions are missing and where they are too many */
                g->vacancy[ic] = g->comp[ic].nMolReq - g->comp[ic].nMolAv;

                /* Save number of ions per compartment prior to swaps */
                g->comp[ic].nMolBefore = g->comp[ic].nMol;
            }
        }

        /* Now actually perform the particle exchanges, one swap group after another */
        gsol = &s->group[static_cast<int>(SwapGroupSplittingType::Solvent)];
        for (int ig = static_cast<int>(SwapGroupSplittingType::Count); ig < s->ngrp; ig++)
        {
            int        nswaps = 0;
            t_swapgrp* g      = &s->group[ig];
            for (auto thisC : gmx::EnumerationWrapper<Compartment>{})
            {
                /* Index to the other compartment */
                auto otherC = thisC == Compartment::A ? Compartment::B : Compartment::A;

                while (g->vacancy[thisC] >= sc->threshold)
                {
                    /* Swap in an ion */

                    /* Get the xc-index of the first atom of a solvent molecule of this compartment */
                    isol = get_index_of_distant_atom(&gsol->comp[thisC], gsol->molname);

                    /* Get the xc-index of a particle from the other compartment */
                    iion = get_index_of_distant_atom(&g->comp[otherC], g->molname);

                    get_molecule_center(&gsol->xc[isol], gsol->apm, gsol->m, com_solvent, s->pbc);
                    get_molecule_center(&g->xc[iion], g->apm, g->m, com_particle, s->pbc);

                    /* Subtract solvent molecule's center of mass and add swap particle's center of mass */
                    translate_positions(&gsol->xc[isol], gsol->apm, com_solvent, com_particle, s->pbc);
                    /* Similarly for the swap particle, subtract com_particle and add com_solvent */
                    translate_positions(&g->xc[iion], g->apm, com_particle, com_solvent, s->pbc);

                    /* Keep track of the changes */
                    g->vacancy[thisC]--;
                    g->vacancy[otherC]++;
                    g->comp[thisC].nMol++;
                    g->comp[otherC].nMol--;
                    g->comp[thisC].inflow_net++;
                    g->comp[otherC].inflow_net--;
                    /* Correct the past time window to still get the right averages from now on */
                    g->comp[thisC].nMolAv++;
                    g->comp[otherC].nMolAv--;
                    for (int j = 0; j < sc->nAverage; j++)
                    {
                        g->comp[thisC].nMolPast[j]++;
                        g->comp[otherC].nMolPast[j]--;
                    }
                    /* Clear ion history */
                    if (MAIN(cr))
                    {
                        int iMol               = iion / g->apm;
                        g->channel_label[iMol] = ChannelHistory::None;
                        g->comp_from[iMol]     = Domain::Notset;
                    }
                    /* That was the swap */
                    nswaps++;
                }
            }

            if (nswaps && bVerbose)
            {
                fprintf(stderr,
                        "%s Performed %d swap%s in step %" PRId64 " for iontype %s.\n",
                        SwS.c_str(),
                        nswaps,
                        nswaps > 1 ? "s" : "",
                        step,
                        g->molname);
            }
        }

        if (s->fpout != nullptr)
        {
            print_ionlist(s, t, "  # after swap");
        }

        /* For the solvent and user-defined swap groups, each rank writes back its
         * (possibly modified) local positions to the official position array. */
        for (int ig = static_cast<int>(SwapGroupSplittingType::Solvent); ig < s->ngrp; ig++)
        {
            t_swapgrp* g = &s->group[ig];
            apply_modified_positions(g, x);
        }

    } /* end of if(bSwap) */

    wallcycle_stop(wcycle, WallCycleCounter::Swap);

    return bSwap;
}
