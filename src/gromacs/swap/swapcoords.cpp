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
#include <numeric>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/localatomset.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdrunutility/handlerestart.h"
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
#include "gromacs/utility/mpicomm.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/vec.h"

namespace gmx
{
class ForceProviders;
class IMDOutputProvider;
class IMdpOptionProvider;
struct MDModulesNotifiers;
} // namespace gmx

//! For output that comes from the swap module
static const std::string SwS = { "SWAP:" };
//! Placeholder for multi-line output
static const std::string SwSEmpty = { "     " };
//! Compartment name
static constexpr gmx::EnumerationArray<Compartment, const char*> CompStr = { "A", "B" };
//! Name for the swap types.
static constexpr gmx::EnumerationArray<SwapType, const char*> SwapStr = { "", "X-", "Y-", "Z-" };
//! Name for the swap dimension.
static const char* const DimStr[DIM + 1] = { "X", "Y", "Z", nullptr };

//! Keep track of through which channel the ions have passed
enum class ChannelHistory : int
{
    None,
    Ch0,
    Ch1,
    Count
};
//! Names for the channels
static constexpr gmx::EnumerationArray<ChannelHistory, const char*> ChannelString = { "none",
                                                                                      "channel0",
                                                                                      "channel1" };

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
//! Names for the domains
static constexpr gmx::EnumerationArray<Domain, const char*> DomainString = { "not_assigned",
                                                                             "Domain_A",
                                                                             "Domain_B" };

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
    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /* notifiers */) override {}
    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /* notifiers */) override {}
};

std::unique_ptr<IMDModule> SwapCoordinatesModuleInfo::create()
{
    return std::make_unique<SwapCoordinates>();
}


} // namespace gmx

/*! \internal \brief
 * Structure containing compartment-specific data.
 */
typedef struct swap_compartment
{
    //! Number of ion or water molecules detected in this compartment.
    int nMol = 0;
    //! Number of molecules before swapping.
    int nMolBefore = 0;
    //! Requested number of molecules in compartment.
    int nMolReq = 0;
    //! Time-averaged number of molecules matching the compartment conditions.
    real nMolAv = 0._real;
    //! Past molecule counts for time-averaging.
    std::vector<int> nMolPast;
    //! Indices to collective array of atoms.
    std::vector<int> ind;
    //! Distance of atom to bulk layer, which is normally the center layer of the compartment
    std::vector<real> dist;
    //! Net inflow of ions into this compartment.
    int inflow_net = 0;
} t_compartment;


/*! \internal \brief
 * This structure contains data needed for the groups involved in swapping:
 * split group 0, split group 1, solvent group, ion groups.
 */
struct t_swapgrp
{
    /*!\brief Construct a swap group given the managed swap atoms.
     *
     * \param[in] atomset Managed indices of atoms that are part of the swap group.
     */
    t_swapgrp(const gmx::LocalAtomSet& atomset);
    //! Name of the group or ion type
    std::string molname;
    //! Number of atoms in each molecule
    int apm = 0;
    //! The atom indices in the swap group
    gmx::LocalAtomSet atomSet;
    //! Collective array of group atom positions (size nat)
    std::vector<gmx::RVec> xc;
    //! Current (collective) shifts (size nat)
    std::vector<gmx::IVec> xc_shifts;
    //! Extra shifts since last DD step (size nat)
    std::vector<gmx::IVec> xc_eshifts;
    //! Old (collective) positions (size nat)
    std::vector<gmx::RVec> xc_old;
    //! Total charge of one molecule of this group
    real q = 0.;
    //! Masses (not relevant for split groups, size apm)
    std::vector<real> m;
    /*! \brief (Collective) Stores from which compartment this
     * molecule has come. This way we keep track of through which
     * channel an ion permeates (size nMol = nat/apm)
     *
     * Sometimes this is filled during checkpoint reading, and
     * sometimes during tpr reading, so we use std::shared_ptr so
     * there will be only one vector and it can be allocated from
     * either place.
     */
    std::shared_ptr<std::vector<Domain>> comp_from;
    //! In which compartment this ion is now (size nMol)
    std::vector<Domain> comp_now;
    /*! \brief Which channel was passed at last by this ion? (size nMol)
     *
     * Sometimes this is filled during checkpoint reading, and
     * sometimes during tpr reading, so we use std::shared_ptr so
     * there will be only one vector and it can be allocated from
     * either place.
     */
    std::shared_ptr<std::vector<ChannelHistory>> channel_label;
    //! Center of the group; COM if masses are used
    gmx::RVec center = { 0, 0, 0 };
    //! Distribution of particles of this group across the two compartments
    gmx::EnumerationArray<Compartment, t_compartment> comp;
    //! How many molecules need to be swapped in?
    gmx::EnumerationArray<Compartment, real> vacancy = { { 0 } };
    //! Net flux of ions per channel
    gmx::EnumerationArray<Channel, int> fluxfromAtoB = { { 0 } };
    //! Number of ions residing in a channel
    gmx::EnumerationArray<Channel, int> nCyl = { { 0 } };
    //! Ions assigned to cyl0 and cyl1. Not good.
    int nCylBoth = 0;
};

t_swapgrp::t_swapgrp(const gmx::LocalAtomSet& atomset) : atomSet{ atomset } {}

/*! \internal \brief
 * Main (private) data structure for the position swapping protocol.
 */
class SwapCoords::Impl
{
public:
    ~Impl()
    {
        if (fpout)
        {
            gmx_fio_fclose(fpout);
        }
    }

    //! One of XX, YY, ZZ
    int swapdim = -1;
    //! Needed to make molecules whole before exchanges.
    t_pbc pbc;
    //! Output file.
    FILE* fpout = nullptr;
    //! Groups for channels, solvent and ions
    std::vector<t_swapgrp> groups;
    //! Flux not going through any of the channels.
    int fluxleak = 0;
    //! The charge imbalance between the compartments.
    real deltaQ = 0._real;
    //! Three required groups: channel0, channel1, solvent
    static constexpr int sc_numRequiredGroups = static_cast<int>(SwapGroupSplittingType::Count);
    //! Return an ArrayRef to the split groups (a.k.a. channel groups)
    gmx::ArrayRef<t_swapgrp> splitGroups()
    {
        return { groups.data() + static_cast<int>(SwapGroupSplittingType::Split0),
                 groups.data() + static_cast<int>(SwapGroupSplittingType::Split1) + 1 };
    }
    //! Return an ArrayRef to the ion groups
    gmx::ArrayRef<t_swapgrp> ionGroups()
    {
        return { groups.data() + sc_numRequiredGroups, groups.data() + groups.size() };
    }
    //! Return an ArrayRef-to-const to the ion groups
    gmx::ArrayRef<const t_swapgrp> ionGroups() const
    {
        return { groups.data() + sc_numRequiredGroups, groups.data() + groups.size() };
    }
    //! Return an ArrayRef to the solvent and ion groups
    gmx::ArrayRef<t_swapgrp> solventAndIonGroups()
    {
        return { groups.data() + static_cast<int>(SwapGroupSplittingType::Solvent),
                 groups.data() + groups.size() };
    }
    //! Return an ArrayRef-to-const to the solvent and ion groups
    gmx::ArrayRef<const t_swapgrp> solventAndIonGroups() const
    {
        return { groups.data() + static_cast<int>(SwapGroupSplittingType::Solvent),
                 groups.data() + groups.size() };
    }
    //! Return a reference to the split or solvent group \c groupType
    t_swapgrp& requiredGroup(SwapGroupSplittingType groupType)
    {
        return groups[static_cast<int>(groupType)];
    }
    //! Return a const reference to the split or solvent group \c groupType
    const t_swapgrp& requiredGroup(SwapGroupSplittingType groupType) const
    {
        return groups[static_cast<int>(groupType)];
    }
    //! Return a reference to the ion group \c index (i.e within the set of ion groups)
    t_swapgrp& ionGroup(size_t index) { return groups[index + sc_numRequiredGroups]; }
    //! Return a const reference to the ion group \c index (i.e within the set of ion groups)
    const t_swapgrp& ionGroup(size_t index) const { return groups[index + sc_numRequiredGroups]; }
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
static gmx_bool is_in_channel(const gmx::RVec point,
                              const gmx::RVec center,
                              real            d_up,
                              real            d_down,
                              real            r_cyl2,
                              t_pbc*          pbc,
                              int             normal)
{
    gmx::RVec dr;
    int       plane1, plane2; /* Directions tangential to membrane */


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
static void print_ionlist(SwapCoords::Impl* s, double time, const char comment[])
{
    // Output time
    fprintf(s->fpout, "%12.5e", time);

    // Output number of molecules and difference to reference counts for each
    // compartment and ion type
    for (auto iComp : gmx::EnumerationWrapper<Compartment>{})
    {
        for (const t_swapgrp& group : s->ionGroups())
        {
            const t_compartment& comp = group.comp[iComp];
            fprintf(s->fpout, "%10d%10.1f%10d", comp.nMol, comp.nMolAv - comp.nMolReq, comp.inflow_net);
        }
    }

    // Output center of split groups
    fprintf(s->fpout,
            "%10g%10g",
            s->requiredGroup(SwapGroupSplittingType::Split0).center[s->swapdim],
            s->requiredGroup(SwapGroupSplittingType::Split1).center[s->swapdim]);

    // Output ion flux for each channel and ion type
    for (auto iChan : gmx::EnumerationWrapper<Channel>{})
    {
        for (const t_swapgrp& group : s->ionGroups())
        {
            fprintf(s->fpout, "%10d", group.fluxfromAtoB[iChan]);
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
static gmx::RVec get_molecule_center(gmx::ArrayRef<const gmx::RVec> x,
                                     gmx::ArrayRef<const real>      weights,
                                     t_pbc*                         pbc)
{
    /* Use the first atom as the reference and put other atoms near that one */
    /* This does not work for large molecules that span > half of the box! */
    const gmx::RVec reference = x[0];

    /* Calculate either the weighted center or simply the center of geometry */
    real      wsum   = 0.0;
    gmx::RVec center = { 0, 0, 0 };
    for (size_t i = 0; i != x.size(); i++)
    {
        /* PBC distance between position and reference */
        gmx::RVec dx;
        pbc_dx(pbc, x[i], reference, dx);

        /* Add PBC distance to reference */
        const gmx::RVec correctPBCimage = reference + dx;

        /* Take weight into account */
        const real wi = weights.empty() ? 1.0 : weights[i];
        wsum += wi;
        const gmx::RVec weightedPBCimage = wi * correctPBCimage;

        /* Add to center */
        center += weightedPBCimage;
    }

    // Normalize and return
    return center / wsum;
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
 * This atom is expected to be the first atom of a molecule.
 *
 * \param[in]     ci        Index of this ion in the collective xc array.
 * \param[inout]  comp      Compartment to add this atom to.
 * \param[in]     distance  Shortest distance of this atom to the bulk layer,
 *                          from which ion/water pairs are selected for swapping.
 */
static void add_to_list(int ci, t_compartment* comp, real distance)
{
    comp->ind.push_back(ci);
    comp->dist.push_back(distance);
    comp->nMol++;
}


/*! \brief Determine the compartment boundaries from the channel centers. */
static void get_compartment_boundaries(Compartment c, SwapCoords::Impl* s, const matrix box, real* left, real* right)
{
    real pos0, pos1;
    real leftpos, rightpos, leftpos_orig;


    if (c >= Compartment::Count)
    {
        gmx_fatal(FARGS, "Compartment out of range");
    }

    pos0 = s->requiredGroup(SwapGroupSplittingType::Split0).center[s->swapdim];
    pos1 = s->requiredGroup(SwapGroupSplittingType::Split1).center[s->swapdim];

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
static void detect_flux_per_channel(t_swapgrp*          group,
                                    int                 iAtom,
                                    Compartment         comp,
                                    const gmx::RVec     atomPosition,
                                    Domain*             comp_now,
                                    Domain*             comp_from,
                                    ChannelHistory*     channel_label,
                                    const t_swapcoords& sc,
                                    SwapCoords::Impl*   s,
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
                            s->requiredGroup(SwapGroupSplittingType::Split0).center,
                            sc.cyl0u,
                            sc.cyl0l,
                            cyl0_r2,
                            &s->pbc,
                            sd);
    in_cyl1 = is_in_channel(atomPosition,
                            s->requiredGroup(SwapGroupSplittingType::Split1).center,
                            sc.cyl1u,
                            sc.cyl1l,
                            cyl1_r2,
                            &s->pbc,
                            sd);

    if (in_cyl0 && in_cyl1)
    {
        /* Ion appears to be in both channels. Something is severely wrong! */
        group->nCylBoth++;
        *comp_now      = Domain::Notset;
        *comp_from     = Domain::Notset;
        *channel_label = ChannelHistory::None;
    }
    else if (in_cyl0)
    {
        /* Ion is in channel 0 now */
        *channel_label = ChannelHistory::Ch0;
        *comp_now      = Domain::Notset;
        group->nCyl[Channel::Zero]++;
    }
    else if (in_cyl1)
    {
        /* Ion is in channel 1 now */
        *channel_label = ChannelHistory::Ch1;
        *comp_now      = Domain::Notset;
        group->nCyl[Channel::One]++;
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
                    group->fluxfromAtoB[chan_nr]++;
                }
                else
                {
                    group->fluxfromAtoB[chan_nr]--;
                }
                fprintf(fpout, "# Atom nr. %d finished passing %s.\n", iAtom, ChannelString[*channel_label]);
                break;
            default:
                gmx_fatal(FARGS,
                          "%s Unknown channel history entry for ion type '%s'\n",
                          SwS.c_str(),
                          group->molname.c_str());
        }

        /* This ion has moved to the _other_ compartment ... */
        *comp_from = *comp_now;
        /* ... and it did not pass any channel yet */
        *channel_label = ChannelHistory::None;
    }
}


/*! \brief Determines which ions or solvent molecules are in compartment A and B */
static void sortMoleculesIntoCompartments(t_swapgrp*          group,
                                          const gmx::MpiComm& mpiComm,
                                          const t_swapcoords& sc,
                                          SwapCoords::Impl*   s,
                                          const matrix        box,
                                          int64_t             step,
                                          FILE*               fpout,
                                          gmx_bool            bRerun,
                                          gmx_bool            bIsSolvent)
{
    gmx::EnumerationArray<Compartment, int> nMolNotInComp; /* consistency check */
    real                                    cyl0_r2 = sc.cyl0r * sc.cyl0r;
    real                                    cyl1_r2 = sc.cyl1r * sc.cyl1r;

    /* Get us a counter that cycles in the range of [0 ... sc.nAverage[ */
    int replace = (step / sc.nstswap) % sc.nAverage;

    for (auto comp : gmx::EnumerationWrapper<Compartment>{})
    {
        real left, right;

        /* Get lists of atoms that match criteria for this compartment */
        get_compartment_boundaries(comp, s, box, &left, &right);

        /* First clear the ion molecule lists */
        group->comp[comp].nMol = 0;
        nMolNotInComp[comp]    = 0; /* consistency check */
        group->comp[comp].ind.clear();
        group->comp[comp].dist.clear();

        /* Loop over the molecules and atoms of this group */
        for (int iMol = 0, iAtom = 0; iAtom < static_cast<int>(group->atomSet.numAtomsGlobal());
             iAtom += group->apm, iMol++)
        {
            real dist;
            int  sd = s->swapdim;

            /* Is this first atom of the molecule in the compartment that we look at? */
            if (compartment_contains_atom(
                        left, right, group->xc[iAtom][sd], box[sd][sd], sc.bulkOffset[comp], &dist))
            {
                /* Add the first atom of this molecule to the list of molecules in this compartment */
                add_to_list(iAtom, &group->comp[comp], dist);

                /* Main also checks for ion groups through which channel each ion has passed */
                if (mpiComm.isMainRank() && !group->comp_now.empty() && !bIsSolvent)
                {
                    int globalAtomNr =
                            group->atomSet.globalIndex()[iAtom] + 1; /* PDB index starts at 1 ... */
                    detect_flux_per_channel(group,
                                            globalAtomNr,
                                            comp,
                                            group->xc[iAtom],
                                            &group->comp_now[iMol],
                                            &(*group->comp_from)[iMol],
                                            &(*group->channel_label)[iMol],
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
            update_time_window(&group->comp[comp], sc.nAverage, replace);
        }
    }

    /* Flux detection warnings */
    if (mpiComm.isMainRank() && !bIsSolvent)
    {
        if (group->nCylBoth > 0)
        {
            fprintf(stderr,
                    "\n"
                    "%s Warning: %d atoms were detected as being in both channels! Probably your "
                    "split\n"
                    "%s          cylinder is way too large, or one compartment has collapsed (step "
                    "%" PRId64 ")\n",
                    SwS.c_str(),
                    group->nCylBoth,
                    SwS.c_str(),
                    step);

            fprintf(s->fpout, "Warning: %d atoms were assigned to both channels!\n", group->nCylBoth);

            group->nCylBoth = 0;
        }
    }

    if (bIsSolvent && nullptr != fpout)
    {
        fprintf(fpout,
                "# Solv. molecules in comp.%s: %d   comp.%s: %d\n",
                CompStr[Compartment::A],
                group->comp[Compartment::A].nMol,
                CompStr[Compartment::B],
                group->comp[Compartment::B].nMol);
    }

    /* Consistency checks */
    const auto numMolecules = static_cast<int>(group->atomSet.numAtomsGlobal() / group->apm);
    if (nMolNotInComp[Compartment::A] + nMolNotInComp[Compartment::B] != numMolecules)
    {
        fprintf(stderr,
                "%s Warning: Inconsistency while assigning '%s' molecules to compartments. !inA: "
                "%d, !inB: %d, total molecules %d\n",
                SwS.c_str(),
                group->molname.c_str(),
                nMolNotInComp[Compartment::A],
                nMolNotInComp[Compartment::B],
                numMolecules);
    }

    int sum = group->comp[Compartment::A].nMol + group->comp[Compartment::B].nMol;
    if (sum != numMolecules)
    {
        fprintf(stderr,
                "%s Warning: %d molecules are in group '%s', but altogether %d have been assigned "
                "to the compartments.\n",
                SwS.c_str(),
                numMolecules,
                group->molname.c_str(),
                sum);
    }
}


/*! \brief Find out how many group atoms are in the compartments initially */
static void get_initial_ioncounts(const t_swapcoords&            sc,
                                  SwapCoords::Impl*              s,
                                  gmx::ArrayRef<const gmx::RVec> x, /* the initial positions */
                                  const matrix                   box,
                                  const gmx::MpiComm&            mpiComm,
                                  gmx_bool                       bRerun)
{
    /* Loop over the user-defined (ion) groups */
    for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
    {
        t_swapgrp& group = s->ionGroup(ig);
        /* Copy the initial positions of the atoms in the group
         * to the collective array so that we can compartmentalize */
        for (size_t i = 0; i < group.atomSet.numAtomsGlobal(); i++)
        {
            int ind     = group.atomSet.globalIndex()[i];
            group.xc[i] = x[ind];
        }

        /* Set up the compartments and get lists of atoms in each compartment */
        sortMoleculesIntoCompartments(&group, mpiComm, sc, s, box, 0, s->fpout, bRerun, FALSE);

        /* Set initial molecule counts if requested (as signaled by "-1" value) */
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            int requested = sc.ionGroup(ig).nmolReq[ic];
            if (requested < 0)
            {
                group.comp[ic].nMolReq = group.comp[ic].nMol;
            }
            else
            {
                group.comp[ic].nMolReq = requested;
            }
        }

        /* Check whether the number of requested molecules adds up to the total number */
        int req = group.comp[Compartment::A].nMolReq + group.comp[Compartment::B].nMolReq;
        int tot = group.comp[Compartment::A].nMol + group.comp[Compartment::B].nMol;

        if ((req != tot))
        {
            gmx_fatal(FARGS,
                      "Mismatch of the number of %s ions summed over both compartments.\n"
                      "You requested a total of %d ions (%d in A and %d in B),\n"
                      "but there are a total of %d ions of this type in the system.\n",
                      group.molname.c_str(),
                      req,
                      group.comp[Compartment::A].nMolReq,
                      group.comp[Compartment::B].nMolReq,
                      tot);
        }

        /* Initialize time-averaging:
         * Write initial concentrations to all time bins to start with */
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            group.comp[ic].nMolAv = group.comp[ic].nMol;
            for (int i = 0; i < sc.nAverage; i++)
            {
                group.comp[ic].nMolPast[i] = group.comp[ic].nMol;
            }
        }
    }
}


/*! \brief Copy history of ion counts from checkpoint file.
 *
 * When called, the checkpoint file has already been read in. Here we copy
 * over the values from .cpt file to the swap data structure.
 */
static void get_initial_ioncounts_from_cpt(const t_swapcoords& sc,
                                           SwapCoords::Impl*   s,
                                           swaphistory_t*      swapstate,
                                           const gmx::MpiComm& mpiComm,
                                           gmx_bool            bVerbose)
{
    if (mpiComm.isMainRank())
    {
        /* Copy the past values from the checkpoint values that have been read in already */
        if (bVerbose)
        {
            fprintf(stderr, "%s Copying values from checkpoint\n", SwS.c_str());
        }

        for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
        {
            t_swapgrp&             group     = s->ionGroup(ig);
            const swapstateIons_t& swapGroup = swapstate->ionType[ig];

            for (auto ic : gmx::EnumerationWrapper<Compartment>{})
            {
                group.comp[ic].nMolReq    = swapGroup.nMolReq[ic];
                group.comp[ic].inflow_net = swapGroup.inflow_net[ic];

                if (bVerbose)
                {
                    fprintf(stderr,
                            "%s ... Influx netto: %d   Requested: %d   Past values: ",
                            SwS.c_str(),
                            group.comp[ic].inflow_net,
                            group.comp[ic].nMolReq);
                }

                for (int j = 0; j < sc.nAverage; j++)
                {
                    group.comp[ic].nMolPast[j] = swapGroup.nMolPast[ic][j];
                    if (bVerbose)
                    {
                        fprintf(stderr, "%d ", group.comp[ic].nMolPast[j]);
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
static void bc_initial_concentrations(const gmx::MpiComm& mpiComm, t_swapcoords* swap, SwapCoords::Impl* s)
{
    for (t_swapgrp& group : s->ionGroups())
    {
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            gmx_bcast(sizeof(group.comp[ic].nMolReq), &(group.comp[ic].nMolReq), mpiComm.comm());
            gmx_bcast(sizeof(group.comp[ic].nMol), &(group.comp[ic].nMol), mpiComm.comm());
            gmx_bcast(swap->nAverage * sizeof(group.comp[ic].nMolPast[0]),
                      group.comp[ic].nMolPast.data(),
                      mpiComm.comm());
        }
    }
}


/*! \brief Ensure that each atom belongs to at most one of the swap groups. */
static void check_swap_groups(SwapCoords::Impl* s, int nat, gmx_bool bVerbose)
{
    if (bVerbose)
    {
        fprintf(stderr, "%s Making sure each atom belongs to at most one of the swap groups.\n", SwS.c_str());
    }

    // This array counts for each atom in the MD system to
    // how many swap groups it belongs (should be 0 or 1!)
    std::vector<int> numGroupsContainingAtom(nat);
    for (const t_swapgrp& group : s->groups)
    {
        for (size_t j = 0; j < group.atomSet.numAtomsGlobal(); j++)
        {
            // Get the global index of this atom of this group
            const int ind = group.atomSet.globalIndex()[j];
            numGroupsContainingAtom[ind]++;
        }
    }
    // Make sure each atom belongs to at most one of the groups
    const int numAtomsInMultipleGroups =
            std::count_if(numGroupsContainingAtom.begin(),
                          numGroupsContainingAtom.end(),
                          [](const int numGroups) { return numGroups > 1; });
    if (numAtomsInMultipleGroups)
    {
        gmx_fatal(FARGS,
                  "%s Cannot perform swapping since %d atom%s allocated to more than one swap "
                  "index group.\n"
                  "%s Each atom must be allocated to at most one of the split groups, the swap "
                  "groups, or the solvent.\n"
                  "%s Check the .mdp file settings regarding the swap index groups or the index "
                  "groups themselves.\n",
                  SwS.c_str(),
                  numAtomsInMultipleGroups,
                  (1 == numAtomsInMultipleGroups) ? " is" : "s are",
                  SwSEmpty.c_str(),
                  SwSEmpty.c_str());
    }
}


/*! \brief Get the number of atoms per molecule for this group.
 *
 * Also ensure that all the molecules in this group have this number of atoms.
 */
static int get_group_apm_check(const t_swapgrp& group, gmx_bool bVerbose, const gmx_mtop_t& mtop)
{
    const int* ind = group.atomSet.globalIndex().data();
    int        nat = group.atomSet.numAtomsGlobal();

    /* Determine the number of solvent atoms per solvent molecule from the
     * first solvent atom: */
    int molb = 0;
    mtopGetMolblockIndex(mtop, ind[0], &molb, nullptr, nullptr);
    const int apm = mtop.moleculeBlockIndices[molb].numAtomsPerMolecule;

    if (bVerbose)
    {
        fprintf(stderr,
                "%s Checking whether all %s molecules consist of %d atom%s\n",
                SwS.c_str(),
                group.molname.c_str(),
                apm,
                apm > 1 ? "s" : "");
    }

    /* Check whether this is also true for all other solvent atoms */
    for (int i = 1; i < nat; i++)
    {
        mtopGetMolblockIndex(mtop, ind[i], &molb, nullptr, nullptr);
        if (apm != mtop.moleculeBlockIndices[molb].numAtomsPerMolecule)
        {
            gmx_fatal(
                    FARGS,
                    "Not all molecules of swap group containing %s molecules consist of %d atoms.",
                    group.molname.c_str(),
                    apm);
        }
    }

    // TODO: check whether charges and masses of each molecule are identical!
    return apm;
}


/*! \brief Print the legend to the swap output file.
 *
 * Also print the initial values of ion counts and position of split groups.
 */
static void print_ionlist_legend(const t_inputrec* ir, SwapCoords::Impl* s, const gmx_output_env_t* oenv)
{
    std::vector<std::string> legend;

    // Number of molecules and difference to reference counts for each
    // compartment and ion type
    for (auto ic : gmx::EnumerationWrapper<Compartment>{})
    {
        for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
        {
            const t_swapgrp& group   = s->ionGroup(ig);
            const char*      molname = group.molname.c_str();
            const real       q       = group.q;
            legend.emplace_back(gmx::formatString(
                    "%s %s ions (charge %s%g)", CompStr[ic], molname, q > 0 ? "+" : "", q));
            legend.emplace_back(gmx::formatString(
                    "%s av. mismatch to %d %s ions", CompStr[ic], group.comp[ic].nMolReq, molname));

            legend.emplace_back(gmx::formatString("%s net %s ion influx", CompStr[ic], molname));
        }
    }

    // Center of split groups
    legend.emplace_back(gmx::formatString(
            "%scenter of %s of split group 0",
            SwapStr[ir->eSwapCoords],
            !s->requiredGroup(SwapGroupSplittingType::Split0).m.empty() ? "mass" : "geometry"));
    legend.emplace_back(gmx::formatString(
            "%scenter of %s of split group 1",
            SwapStr[ir->eSwapCoords],
            !s->requiredGroup(SwapGroupSplittingType::Split1).m.empty() ? "mass" : "geometry"));

    // Ion flux for each channel and ion type
    for (auto ic : gmx::EnumerationWrapper<Channel>{})
    {
        const int label = static_cast<int>(ic);
        for (const t_swapGroup& group : ir->swap->ionGroups())
        {
            legend.emplace_back(
                    gmx::formatString("A->ch%d->B %s permeations", label, group.molname.c_str()));
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
    std::fflush(s->fpout);
}


/*! \brief Initialize channel ion flux detection routine.
 *
 * Initialize arrays that keep track of where the ions come from and where
 * they go.
 */
static void detect_flux_per_channel_init(SwapCoords::Impl* s, swaphistory_t* swapstate, const bool isRestart)
{
    /* All these flux detection routines run on the main only */
    if (swapstate == nullptr)
    {
        return;
    }

    for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
    {
        t_swapgrp&       group     = s->ionGroup(ig);
        swapstateIons_t& swapGroup = swapstate->ionType[ig];

        /******************************************************/
        /* Channel and domain history for the individual ions */
        /******************************************************/
        const int numMolecules = group.atomSet.numAtomsGlobal() / group.apm;
        if (isRestart)
        {
            // Copy the shared pointers, i.e. refer to the same contents
            group.comp_from     = swapGroup.comp_from;
            group.channel_label = swapGroup.channel_label;
        }
        else /* allocate memory for molecule counts */
        {
            group.channel_label = std::make_shared<std::vector<ChannelHistory>>(numMolecules);
            group.comp_from     = std::make_shared<std::vector<Domain>>(numMolecules);
            // Copy the shared pointers, i.e. refer to the same contents
            swapGroup.comp_from     = group.comp_from;
            swapGroup.channel_label = group.channel_label;
        }
        group.comp_now.resize(numMolecules);

        /* Initialize the channel and domain history counters */
        for (int i = 0; i < numMolecules; i++)
        {
            group.comp_now[i] = Domain::Notset;
            if (!isRestart)
            {
                (*group.comp_from)[i]     = Domain::Notset;
                (*group.channel_label)[i] = ChannelHistory::None;
            }
        }

        /************************************/
        /* Channel fluxes for both channels */
        /************************************/
        group.nCyl[Channel::Zero] = 0;
        group.nCyl[Channel::One]  = 0;
        group.nCylBoth            = 0;
    }

    if (isRestart)
    {
        fprintf(stderr, "%s Copying channel fluxes from checkpoint file data\n", SwS.c_str());
    }


    // Loop over ion types (and both channels)
    for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
    {
        t_swapgrp&             group     = s->ionGroup(ig);
        const swapstateIons_t& swapGroup = swapstate->ionType[ig];

        for (auto ic : gmx::EnumerationWrapper<Channel>{})
        {
            fprintf(stderr,
                    "%s Channel %d flux history for ion type %s (charge %g): ",
                    SwS.c_str(),
                    static_cast<int>(ic),
                    group.molname.c_str(),
                    group.q);
            if (isRestart)
            {
                group.fluxfromAtoB[ic] = swapGroup.fluxfromAtoB[ic];
            }
            else
            {
                group.fluxfromAtoB[ic] = 0;
            }

            fprintf(stderr, "%d molecule%s", group.fluxfromAtoB[ic], group.fluxfromAtoB[ic] == 1 ? "" : "s");
            fprintf(stderr, "\n");
        }
    }

    /* Set pointers for checkpoint writing */
    swapstate->fluxleak_p = &s->fluxleak;
    for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
    {
        t_swapgrp&       group     = s->ionGroup(ig);
        swapstateIons_t& swapGroup = swapstate->ionType[ig];

        for (auto ic : gmx::EnumerationWrapper<Channel>{})
        {
            swapGroup.fluxfromAtoB_p[ic] = &group.fluxfromAtoB[ic];
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
    char* env = std::getenv("GMX_COMPELDUMP");

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
static void init_swapstate(swaphistory_t*                 swapstate,
                           SwapCoords::Impl*              s,
                           const gmx_mtop_t&              mtop,
                           gmx::ArrayRef<const gmx::RVec> x, /* the initial positions */
                           const matrix                   box,
                           const t_inputrec*              ir)
{
    /* We always need the last whole positions such that
     * in the next time step we can make the channels whole again in PBC */
    if (swapstate->bFromCpt)
    {
        /* Copy the last whole positions of each channel from .cpt */
        std::copy(swapstate->xc_old_whole[Channel::Zero].begin(),
                  swapstate->xc_old_whole[Channel::Zero].end(),
                  s->requiredGroup(SwapGroupSplittingType::Split0).xc_old.begin());
        std::copy(swapstate->xc_old_whole[Channel::One].begin(),
                  swapstate->xc_old_whole[Channel::One].end(),
                  s->requiredGroup(SwapGroupSplittingType::Split1).xc_old.begin());
    }
    else
    {
        swapstate->eSwapCoords = ir->eSwapCoords;

        /* Set the number of ion types and allocate memory for checkpointing */
        swapstate->ionType.resize(s->ionGroups().size());

        /* Extract the initial split group positions. */

        /* Remove pbc, make molecule whole. */
        std::vector<gmx::RVec> x_pbc(x.begin(), x.end());

        /* This can only make individual molecules whole, not multimers */
        do_pbc_mtop(ir->pbcType, box, &mtop, as_rvec_array(x_pbc.data()));

        /* Output the starting structure? */
        outputStartStructureIfWanted(mtop, as_rvec_array(x_pbc.data()), ir->pbcType, box);

        /* If this is the first run (i.e. no checkpoint present) we assume
         * that the starting positions give us the correct PBC representation */
        for (t_swapgrp& group : s->splitGroups())
        {
            gmx::ArrayRef<const int> globalIndex = group.atomSet.globalIndex();
            for (size_t i = 0; i < group.atomSet.numAtomsGlobal(); i++)
            {
                group.xc_old[i] = x_pbc[globalIndex[i]];
            }
        }

        /* Prepare swapstate arrays for later checkpoint writing */
        swapstate->xc_old_whole[Channel::Zero].resize(
                s->requiredGroup(SwapGroupSplittingType::Split0).atomSet.numAtomsGlobal(), { 0, 0, 0 });
        swapstate->xc_old_whole[Channel::One].resize(
                s->requiredGroup(SwapGroupSplittingType::Split1).atomSet.numAtomsGlobal(), { 0, 0, 0 });
    }

    /* For subsequent checkpoint writing, set the swapstate pointers to the xc_old
     * arrays that get updated at every swapping step */
    swapstate->xc_old_whole_p[Channel::Zero] = &s->requiredGroup(SwapGroupSplittingType::Split0).xc_old;
    swapstate->xc_old_whole_p[Channel::One] = &s->requiredGroup(SwapGroupSplittingType::Split1).xc_old;
}

/*! \brief Determine the total charge imbalance resulting from the swap groups */
static real getRequestedChargeImbalance(SwapCoords::Impl* s)
{
    return std::accumulate(
            s->ionGroups().begin(),
            s->ionGroups().end(),
            0._real,
            [](real subTotal, const t_swapgrp& group) -> real
            {
                return subTotal
                       + group.q * (group.comp[Compartment::A].nMolReq - group.comp[Compartment::B].nMolReq);
            });
}


/*! \brief Sorts anions and cations into two separate groups
 *
 * This routine should be called for the 'anions' and 'cations' group,
 * of which the indices were lumped together in the older version of the code.
 */
static void copyIndicesToGroup(gmx::ArrayRef<const int> indIons, t_swapGroup* group, const gmx::MpiComm& mpiComm)
{
    /* If explicit ion counts were requested in the .mdp file
     * (by setting positive values for the number of ions),
     * we can make an additional consistency check here */
    if ((group->nmolReq[Compartment::A] > 0) && (group->nmolReq[Compartment::B] > 0))
    {
        if (indIons.ssize() != (group->nmolReq[Compartment::A] + group->nmolReq[Compartment::B]))
        {
            gmx_fatal_collective(
                    FARGS,
                    mpiComm.comm(),
                    mpiComm.isMainRank(),
                    "%s Inconsistency while importing swap-related data from an old "
                    "input file version.\n"
                    "%s The requested ion counts in compartments A (%d) and B (%d)\n"
                    "%s do not add up to the number of ions (%zu) of this type for the "
                    "group '%s'.\n",
                    SwS.c_str(),
                    SwSEmpty.c_str(),
                    group->nmolReq[Compartment::A],
                    group->nmolReq[Compartment::B],
                    SwSEmpty.c_str(),
                    indIons.size(),
                    group->molname.c_str());
        }
    }

    group->ind.assign(indIons.begin(), indIons.end());
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
static void convertOldToNewGroupFormat(t_swapcoords*       sc,
                                       const gmx_mtop_t&   mtop,
                                       gmx_bool            bVerbose,
                                       const gmx::MpiComm& mpiComm)
{
    t_swapGroup& group = sc->ionGroup(0);

    /* Loop through the atom indices of group #3 (anions) and put all indices
     * that belong to cations into the cation group.
     */
    std::vector<int> indAnions;
    std::vector<int> indCations;
    indAnions.reserve(group.ind.size());
    indCations.reserve(group.ind.size());

    int molb = 0;
    for (size_t i = 0; i < group.ind.size(); i++)
    {
        const t_atom& atom = mtopGetAtomParameters(mtop, group.ind[i], &molb);
        if (atom.q < 0)
        {
            // This is an anion, add it to the list of anions
            indAnions.push_back(group.ind[i]);
        }
        else
        {
            // This is a cation, add it to the list of cations
            indCations.push_back(group.ind[i]);
        }
    }

    if (bVerbose)
    {
        fprintf(stdout,
                "%s Sorted %zu ions into separate groups of %zu anions and %zu cations.\n",
                SwS.c_str(),
                group.ind.size(),
                indAnions.size(),
                indCations.size());
    }


    /* Now we have the correct lists of anions and cations.
     * Copy it to the right groups.
     */
    copyIndicesToGroup(indAnions, &group, mpiComm);
    group = sc->ionGroup(1);
    copyIndicesToGroup(indCations, &group, mpiComm);
}


/*! \brief Returns TRUE if we started from an old .tpr
 *
 * Then we need to re-sort anions and cations into separate groups */
static gmx_bool bConvertFromOldTpr(t_swapcoords* sc)
{
    // If the second group has no atoms it means we need to convert!
    return (sc->ionGroups().size() >= 2) && sc->ionGroup(1).ind.empty();
}

std::unique_ptr<SwapCoords> init_swapcoords(FILE*                       fplog,
                                            const t_inputrec*           ir,
                                            const char*                 fn,
                                            const gmx_mtop_t&           mtop,
                                            const t_state*              globalState,
                                            ObservablesHistory*         oh,
                                            const gmx::MpiComm&         mpiComm,
                                            const gmx_domdec_t*         dd,
                                            gmx::LocalAtomSetManager*   atomSets,
                                            const gmx_output_env_t*     oenv,
                                            const gmx::MdrunOptions&    mdrunOptions,
                                            const gmx::StartingBehavior startingBehavior)
{
    swaphistory_t* swapstate = nullptr;

    if (mpiComm.size() > 1 && dd == nullptr)
    {
        gmx_fatal(FARGS, "Position swapping is only implemented for domain decomposition!");
    }

    t_swapcoords*               scModifiable = ir->swap.get();
    std::unique_ptr<SwapCoords> swapCoords   = std::make_unique<SwapCoords>();
    SwapCoords::Impl*           s            = swapCoords->impl_.get();
    // TODO is this a good idea?
    // set_pbc(&s->pbc, ir->pbcType, globalState->box);

    if (mdrunOptions.rerun)
    {
        if (mpiComm.size() > 1)
        {
            gmx_fatal(FARGS,
                      "%s This module does not support reruns in parallel\nPlease request a serial "
                      "run with -nt 1 / -np 1\n",
                      SwS.c_str());
        }

        fprintf(stderr, "%s Rerun - using every available frame\n", SwS.c_str());
        scModifiable->nstswap  = 1;
        scModifiable->nAverage = 1; /* averaging makes no sense for reruns */
    }

    if (mpiComm.isMainRank() && startingBehavior == gmx::StartingBehavior::NewSimulation)
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
    if (bConvertFromOldTpr(scModifiable))
    {
        convertOldToNewGroupFormat(scModifiable, mtop, bVerbose && mpiComm.isMainRank(), mpiComm);
    }
    // Now sc can be const
    const t_swapcoords& sc = *ir->swap;

    /* Copy some data and pointers to the group structures for convenience */
    s->groups.reserve(sc.groups.size());
    for (const t_swapGroup& swapGroup : sc.groups)
    {
        s->groups.emplace_back(atomSets->add(swapGroup.ind));
        s->groups.back().molname = swapGroup.molname;
    }

    /* Check for overlapping atoms */
    check_swap_groups(s, mtop.natoms, bVerbose && mpiComm.isMainRank());

    /* Allocate space for the collective arrays for all groups */
    for (int ig = 0; ig != gmx::ssize(s->groups); ++ig)
    {
        t_swapgrp& group = s->groups[ig];
        // Allocate space for the collective arrays for all groups
        const size_t numAtomsGlobal = group.atomSet.numAtomsGlobal();
        // Allocate for the collective position array
        group.xc.resize(numAtomsGlobal, { 0, 0, 0 });
    }
    for (t_swapgrp& group : s->splitGroups())
    {
        // For the split groups (the channels) we need some extra memory to
        // be able to make the molecules whole even if they span more than
        // half of the box size.
        const size_t numAtomsGlobal = group.atomSet.numAtomsGlobal();
        group.xc_shifts.resize(numAtomsGlobal, { 0, 0, 0 });
        group.xc_eshifts.resize(numAtomsGlobal, { 0, 0, 0 });
        group.xc_old.resize(numAtomsGlobal, { 0, 0, 0 });
    }

    if (mpiComm.isMainRank())
    {
        if (oh->swapHistory == nullptr)
        {
            oh->swapHistory = std::make_unique<swaphistory_t>(swaphistory_t{});
        }
        swapstate = oh->swapHistory.get();

        init_swapstate(swapstate, s, mtop, globalState->x, globalState->box, ir);
    }

    /* After init_swapstate we have a set of (old) whole positions for our
     * channels. Now transfer that to all nodes */
    if (mpiComm.size() > 1)
    {
        for (t_swapgrp& group : s->splitGroups())
        {
            gmx_bcast(group.xc_old.size() * sizeof(group.xc_old[0]), group.xc_old.data(), mpiComm.comm());
        }
    }

    /* Make sure that all molecules in the solvent and ion groups contain the
     * same number of atoms each */
    for (t_swapgrp& group : s->solventAndIonGroups())
    {
        group.apm = get_group_apm_check(group, mpiComm.isMainRank() && bVerbose, mtop);

        /* Since all molecules of a group are equal, we only need enough space
         * to determine properties of a single molecule at at time */
        group.m.resize(group.apm); /* For the center of mass */
        real charge = 0;           /* To determine the total charge */
        int  molb   = 0;
        for (int j = 0; j < group.apm; j++)
        {
            const t_atom& atom = mtopGetAtomParameters(mtop, group.atomSet.globalIndex()[j], &molb);
            group.m[j]         = atom.m;
            charge += atom.q;
        }
        /* Total charge of one molecule of this group: */
        group.q = charge;
    };

    for (size_t ig = 0; ig != s->splitGroups().size(); ++ig)
    {
        t_swapgrp& group = s->splitGroups()[ig];
        if (sc.massw_split[ig])
        {
            // Save the split-group masses if mass-weighting is requested
            group.m.resize(group.atomSet.numAtomsGlobal());
            int molb = 0;
            for (size_t i = 0; i < group.atomSet.numAtomsGlobal(); i++)
            {
                group.m[i] = mtopGetAtomMass(mtop, group.atomSet.globalIndex()[i], &molb);
            }
        }
    }

    bool restartWithAppending = (startingBehavior == gmx::StartingBehavior::RestartWithAppending);
    if (mpiComm.isMainRank())
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

            FILE* fpout              = s->fpout;
            auto  reportGroupDetails = [fpout](const t_swapgrp& group,
                                              const bool       reportMoleculeDetails,
                                              const char*      name) -> void
            {
                fprintf(fpout,
                        "# %s group '%s' contains %zu atom%s",
                        name,
                        group.molname.c_str(),
                        group.atomSet.numAtomsGlobal(),
                        (group.atomSet.numAtomsGlobal() > 1) ? "s" : "");
                if (reportMoleculeDetails)
                {
                    fprintf(fpout,
                            " with %d atom%s in each molecule of charge %g",
                            group.apm,
                            (group.apm > 1) ? "s" : "",
                            group.q);
                }
                fprintf(fpout, ".\n");
            };
            reportGroupDetails(s->requiredGroup(SwapGroupSplittingType::Split0),
                               false,
                               enumValueToString(SwapGroupSplittingType::Split0));
            reportGroupDetails(s->requiredGroup(SwapGroupSplittingType::Split1),
                               false,
                               enumValueToString(SwapGroupSplittingType::Split1));
            reportGroupDetails(s->requiredGroup(SwapGroupSplittingType::Solvent),
                               true,
                               enumValueToString(SwapGroupSplittingType::Solvent));
            for (const t_swapgrp& group : s->ionGroups())
            {
                reportGroupDetails(group, true, "Ion");
            }

            fprintf(s->fpout, "#\n# Initial positions of split groups:\n");
        }

        for (const auto groupType : { SwapGroupSplittingType::Split0, SwapGroupSplittingType::Split1 })
        {
            t_swapgrp& group = s->requiredGroup(groupType);

            const gmx::ArrayRef<const int> indices(sc.requiredGroup(groupType).ind);
            for (size_t i = 0; i < group.atomSet.numAtomsGlobal(); i++)
            {
                group.xc[i] = globalState->x[indices[i]];
            }
            // xc has the correct PBC representation for the two channels, so we do
            // not need to correct for that
            get_center(as_rvec_array(group.xc.data()), group.m.data(), group.atomSet.numAtomsGlobal(), group.center);

            if (!restartWithAppending)
            {
                fprintf(s->fpout,
                        "# %s group %s-center %5f nm\n",
                        enumValueToString(groupType),
                        DimStr[s->swapdim],
                        group.center[s->swapdim]);
            }
        }

        if (!restartWithAppending)
        {
            if ((0 != sc.bulkOffset[Compartment::A]) || (0 != sc.bulkOffset[Compartment::B]))
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
                fprintf(s->fpout, "# bulk-offsetA = %g\n", sc.bulkOffset[Compartment::A]);
                fprintf(s->fpout, "# bulk-offsetB = %g\n", sc.bulkOffset[Compartment::B]);
            }

            fprintf(s->fpout, "#\n");
            fprintf(s->fpout, "# Split0 cylinder radius %f nm, up %f nm, down %f nm\n", sc.cyl0r, sc.cyl0u, sc.cyl0l);
            fprintf(s->fpout, "# Split1 cylinder radius %f nm, up %f nm, down %f nm\n", sc.cyl1r, sc.cyl1u, sc.cyl1l);

            fprintf(s->fpout, "#\n");
            if (!mdrunOptions.rerun)
            {
                fprintf(s->fpout,
                        "# Coupling constant (number of swap attempt steps to average over): %d  "
                        "(translates to %f ps).\n",
                        sc.nAverage,
                        sc.nAverage * sc.nstswap * ir->delta_t);
                fprintf(s->fpout, "# Threshold is %f\n", sc.threshold);
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
    for (t_swapgrp& group : s->ionGroups())
    {
        for (auto ic : keysOf(group.comp))
        {
            group.comp[ic].nMolPast.resize(sc.nAverage);
        }
    }

    /* Get the initial particle concentrations and let the other nodes know */
    if (mpiComm.isMainRank())
    {
        if (startingBehavior != gmx::StartingBehavior::NewSimulation)
        {
            get_initial_ioncounts_from_cpt(sc, s, swapstate, mpiComm, bVerbose);
        }
        else
        {
            fprintf(stderr, "%s Determining initial numbers of ions per compartment.\n", SwS.c_str());
            get_initial_ioncounts(sc, s, globalState->x, globalState->box, mpiComm, mdrunOptions.rerun);
        }

        /* Prepare (further) checkpoint writes ... */
        if (startingBehavior != gmx::StartingBehavior::NewSimulation)
        {
            /* Consistency check */
            if (swapstate->nAverage != sc.nAverage)
            {
                gmx_fatal(FARGS,
                          "%s Ion count averaging steps mismatch! checkpoint: %d, tpr: %d",
                          SwS.c_str(),
                          swapstate->nAverage,
                          sc.nAverage);
            }
        }
        else
        {
            swapstate->nAverage = sc.nAverage;
        }
        fprintf(stderr, "%s Setting pointers for checkpoint writing\n", SwS.c_str());
        for (auto ic : gmx::EnumerationWrapper<Compartment>{})
        {
            for (size_t ig = 0; ig != s->ionGroups().size(); ++ig)
            {
                t_swapgrp&       group     = s->ionGroup(ig);
                swapstateIons_t& swapGroup = swapstate->ionType[ig];

                swapGroup.nMolReq_p[ic]    = &(group.comp[ic].nMolReq);
                swapGroup.nMolPast_p[ic]   = &(group.comp[ic].nMolPast[0]);
                swapGroup.inflow_net_p[ic] = &(group.comp[ic].inflow_net);
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

    if (mpiComm.size() > 1)
    {
        bc_initial_concentrations(mpiComm, ir->swap.get(), s);
    }

    /* Update the time-averaged number of molecules for all groups and compartments */
    for (t_swapgrp& group : s->ionGroups())
    {
        for (auto ic : keysOf(group.comp))
        {
            update_time_window(&group.comp[ic], sc.nAverage, -1);
        }
    }

    /* Initialize arrays that keep track of through which channel the ions go */
    detect_flux_per_channel_init(s, swapstate, startingBehavior != gmx::StartingBehavior::NewSimulation);

    /* We need to print the legend if we open this file for the first time. */
    if (mpiComm.isMainRank() && !restartWithAppending)
    {
        print_ionlist_legend(ir, s, oenv);
    }
    return swapCoords;
}

SwapCoords::SwapCoords() : impl_(std::make_unique<Impl>()) {}

SwapCoords::~SwapCoords() = default;

/*! \brief Do we need to swap a molecule in any of the ion groups with a water molecule at this step?
 *
 * From the requested and average molecule counts we determine whether a swap is needed
 * at this time step.
 */
static gmx_bool need_swap(const t_swapcoords& sc, SwapCoords::Impl* s)
{
    for (const t_swapgrp& group : s->ionGroups())
    {
        for (auto ic : keysOf(group.comp))
        {
            if (group.comp[ic].nMolReq - group.comp[ic].nMolAv >= sc.threshold)
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
static int get_index_of_distant_atom(t_compartment* comp, const std::string& molname)
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
                  molname.c_str(),
                  comp->nMolBefore,
                  molname.c_str());
    }

    /* Set the distance of this index to infinity such that it won't get selected again in
     * this time step
     */
    comp->dist[ibest] = GMX_REAL_MAX;

    return comp->ind[ibest];
}


/*! \brief Swaps centers of mass and makes molecule whole if broken */
static void translate_positions(gmx::ArrayRef<gmx::RVec> moleculeX,
                                const gmx::RVec          old_com,
                                const gmx::RVec          new_com,
                                t_pbc*                   pbc)
{
    /* Use the first atom as the reference for PBC */
    const gmx::RVec reference = moleculeX[0];

    for (gmx::RVec& x : moleculeX)
    {
        /* PBC distance between position and reference */
        gmx::RVec dx;
        pbc_dx(pbc, x, reference, dx);

        /* Add PBC distance to reference */
        gmx::RVec correctPBCimage = reference + dx;

        /* Subtract old_com from correct image and add new_com */
        correctPBCimage -= old_com;
        correctPBCimage += new_com;

        x = correctPBCimage;
    }
}


/*! \brief Write back the modified local positions from the collective array to the official positions. */
static void apply_modified_positions(const t_swapgrp& group, gmx::ArrayRef<gmx::RVec> x)
{
    auto collectiveIndex = group.atomSet.collectiveIndex().begin();
    for (const auto localIndex : group.atomSet.localIndex())
    {
        /* Copy the possibly modified position */
        x[localIndex] = group.xc[*collectiveIndex];
        ++collectiveIndex;
    }
}


gmx_bool do_swapcoords(const gmx::MpiComm&      mpiComm,
                       int64_t                  step,
                       double                   t,
                       const t_inputrec*        ir,
                       SwapCoords*              swapCoords,
                       gmx_wallcycle*           wcycle,
                       gmx::ArrayRef<gmx::RVec> x,
                       matrix                   box,
                       gmx_bool                 bVerbose,
                       gmx_bool                 bRerun)
{
    SwapCoords::Impl*   s     = swapCoords->impl_.get();
    const t_swapcoords& sc    = *ir->swap;
    gmx_bool            bSwap = FALSE;


    wallcycle_start(wcycle, WallCycleCounter::Swap);

    set_pbc(&s->pbc, ir->pbcType, box);

    /* Assemble the positions of the split groups, i.e. the channels.
     * Here we also pass a shifts array to communicate_group_positions(), so that it can make
     * the molecules whole even in cases where they span more than half of the box in
     * any dimension */
    for (t_swapgrp& group : s->splitGroups())
    {
        communicate_group_positions(mpiComm,
                                    as_rvec_array(group.xc.data()),
                                    as_ivec_array(group.xc_shifts.data()),
                                    as_ivec_array(group.xc_eshifts.data()),
                                    TRUE,
                                    as_rvec_array(x.data()),
                                    group.atomSet.numAtomsGlobal(),
                                    group.atomSet.numAtomsLocal(),
                                    group.atomSet.localIndex().data(),
                                    group.atomSet.collectiveIndex().data(),
                                    as_rvec_array(group.xc_old.data()),
                                    box);

        get_center(as_rvec_array(group.xc.data()),
                   group.m.data(),
                   group.atomSet.numAtomsGlobal(),
                   group.center); /* center of split groups == channels */
    };

    /* Assemble the positions of the ions (ig = 3, 4, ...). These molecules should
     * be small and we can always make them whole with a simple distance check.
     * Therefore we pass NULL as third argument to communicate_group_positions. */
    for (t_swapgrp& group : s->ionGroups())
    {
        communicate_group_positions(mpiComm,
                                    as_rvec_array(group.xc.data()),
                                    nullptr,
                                    nullptr,
                                    FALSE,
                                    as_rvec_array(x.data()),
                                    group.atomSet.numAtomsGlobal(),
                                    group.atomSet.numAtomsLocal(),
                                    group.atomSet.localIndex().data(),
                                    group.atomSet.collectiveIndex().data(),
                                    nullptr,
                                    nullptr);

        /* Determine how many ions of this type each compartment contains */
        sortMoleculesIntoCompartments(&group, mpiComm, sc, s, box, step, s->fpout, bRerun, FALSE);
    }

    /* Output how many ions are in the compartments */
    if (mpiComm.isMainRank())
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
        t_swapgrp& solventGroup = s->requiredGroup(SwapGroupSplittingType::Solvent);
        communicate_group_positions(mpiComm,
                                    as_rvec_array(solventGroup.xc.data()),
                                    nullptr,
                                    nullptr,
                                    FALSE,
                                    as_rvec_array(x.data()),
                                    solventGroup.atomSet.numAtomsGlobal(),
                                    solventGroup.atomSet.numAtomsLocal(),
                                    solventGroup.atomSet.localIndex().data(),
                                    solventGroup.atomSet.collectiveIndex().data(),
                                    nullptr,
                                    nullptr);

        /* Determine how many molecules of solvent each compartment contains */
        sortMoleculesIntoCompartments(&solventGroup, mpiComm, sc, s, box, step, s->fpout, bRerun, TRUE);

        /* Save number of solvent molecules per compartment prior to any swaps */
        solventGroup.comp[Compartment::A].nMolBefore = solventGroup.comp[Compartment::A].nMol;
        solventGroup.comp[Compartment::B].nMolBefore = solventGroup.comp[Compartment::B].nMol;

        for (t_swapgrp& group : s->ionGroups())
        {
            for (auto ic : keysOf(group.comp))
            {
                /* Determine in which compartment ions are missing and where they are too many */
                group.vacancy[ic] = group.comp[ic].nMolReq - group.comp[ic].nMolAv;

                /* Save number of ions per compartment prior to swaps */
                group.comp[ic].nMolBefore = group.comp[ic].nMol;
            }
        }

        /* Now actually perform the particle exchanges, one swap group after another */
        for (t_swapgrp& group : s->ionGroups())
        {
            int nswaps = 0;
            for (auto thisC : gmx::EnumerationWrapper<Compartment>{})
            {
                /* Index to the other compartment */
                auto otherC = thisC == Compartment::A ? Compartment::B : Compartment::A;

                while (group.vacancy[thisC] >= sc.threshold)
                {
                    /* Swap in an ion */

                    /* Get the xc-index of the first atom of a solvent molecule of this compartment */
                    int isol = get_index_of_distant_atom(&solventGroup.comp[thisC], solventGroup.molname);
                    gmx::ArrayRef<gmx::RVec> solventMoleculePositions =
                            gmx::ArrayRef<gmx::RVec>(solventGroup.xc).subArray(isol, solventGroup.m.size());

                    /* Get the xc-index of an ion from the other compartment */
                    int iion = get_index_of_distant_atom(&group.comp[otherC], group.molname);
                    gmx::ArrayRef<gmx::RVec> ionPositions =
                            gmx::ArrayRef<gmx::RVec>(group.xc).subArray(iion, group.m.size());

                    /* solvent and swap molecule's center of mass */
                    const gmx::RVec com_solvent =
                            get_molecule_center(solventMoleculePositions, solventGroup.m, &s->pbc);
                    const gmx::RVec com_ion = get_molecule_center(ionPositions, group.m, &s->pbc);

                    /* Subtract solvent molecule's center of mass and add swap particle's center of mass */
                    translate_positions(solventMoleculePositions, com_solvent, com_ion, &s->pbc);
                    /* Similarly for the swap particle, subtract com_particle and add com_solvent */
                    translate_positions(ionPositions, com_ion, com_solvent, &s->pbc);

                    /* Keep track of the changes */
                    group.vacancy[thisC]--;
                    group.vacancy[otherC]++;
                    group.comp[thisC].nMol++;
                    group.comp[otherC].nMol--;
                    group.comp[thisC].inflow_net++;
                    group.comp[otherC].inflow_net--;
                    /* Correct the past time window to still get the right averages from now on */
                    group.comp[thisC].nMolAv++;
                    group.comp[otherC].nMolAv--;
                    for (int j = 0; j < sc.nAverage; j++)
                    {
                        group.comp[thisC].nMolPast[j]++;
                        group.comp[otherC].nMolPast[j]--;
                    }
                    /* Clear ion history */
                    if (mpiComm.isMainRank())
                    {
                        int iMol                     = iion / group.apm;
                        (*group.channel_label)[iMol] = ChannelHistory::None;
                        (*group.comp_from)[iMol]     = Domain::Notset;
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
                        group.molname.c_str());
            }
        }

        if (s->fpout != nullptr)
        {
            print_ionlist(s, t, "  # after swap");
        }

        /* For the solvent and user-defined swap groups, each rank writes back its
         * (possibly modified) local positions to the official position array. */
        for (t_swapgrp& group : s->solventAndIonGroups())
        {
            apply_modified_positions(group, x);
        }

    } /* end of if(bSwap) */

    wallcycle_stop(wcycle, WallCycleCounter::Swap);

    return bSwap;
}
