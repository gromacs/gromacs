/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2021, by the GROMACS development team, led by
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
 *
 * \brief This file defines functions used by the domdec module
 * while managing the construction, use and error checking for
 * topologies local to a DD rank.
 *
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "gromacs/domdec/localtopologychecker.h"

#include <optional>
#include <string>
#include <vector>

#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/reversetopology.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/textwriter.h"

#include "dump.h"

namespace gmx
{

/*! \brief Checks whether interactions have been assigned for one function type
 *
 * Loops over a list of interactions in the local topology of one function type
 * and flags each of the interactions as assigned in the global \p isAssigned list.
 * Exits with an inconsistency error when an interaction is assigned more than once.
 */
static void flagInteractionsForType(const int              ftype,
                                    const InteractionList& il,
                                    const reverse_ilist_t& ril,
                                    const Range<int>&      atomRange,
                                    const int              numAtomsPerMolecule,
                                    ArrayRef<const int>    globalAtomIndices,
                                    ArrayRef<int>          isAssigned)
{
    const int nril_mol = ril.index[numAtomsPerMolecule];
    const int nral     = NRAL(ftype);

    for (int i = 0; i < il.size(); i += 1 + nral)
    {
        // ia[0] is the interaction type, ia[1, ...] the atom indices
        const int* ia = il.iatoms.data() + i;
        // Extract the global atom index of the first atom in this interaction
        const int a0 = globalAtomIndices[ia[1]];
        /* Check if this interaction is in
         * the currently checked molblock.
         */
        if (atomRange.isInRange(a0))
        {
            // The molecule index in the list of this molecule type
            const int moleculeIndex = (a0 - atomRange.begin()) / numAtomsPerMolecule;
            const int atomOffset = (a0 - atomRange.begin()) - moleculeIndex * numAtomsPerMolecule;
            const int globalAtomStartInMolecule = atomRange.begin() + moleculeIndex * numAtomsPerMolecule;
            int       j_mol                     = ril.index[atomOffset];
            bool found                          = false;
            while (j_mol < ril.index[atomOffset + 1] && !found)
            {
                const int j       = moleculeIndex * nril_mol + j_mol;
                const int ftype_j = ril.il[j_mol];
                /* Here we need to check if this interaction has
                 * not already been assigned, since we could have
                 * multiply defined interactions.
                 */
                if (ftype == ftype_j && ia[0] == ril.il[j_mol + 1] && isAssigned[j] == 0)
                {
                    /* Check the atoms */
                    found = true;
                    for (int a = 0; a < nral; a++)
                    {
                        if (globalAtomIndices[ia[1 + a]]
                            != globalAtomStartInMolecule + ril.il[j_mol + 2 + a])
                        {
                            found = false;
                        }
                    }
                    if (found)
                    {
                        isAssigned[j] = 1;
                    }
                }
                j_mol += 2 + nral_rt(ftype_j);
            }
            if (!found)
            {
                gmx_incons("Some interactions seem to be assigned multiple times");
            }
        }
    }
}

/*! \brief Help print error output when interactions are missing in a molblock
 *
 * \note This function needs to be called on all ranks (contains a global summation)
 */
static std::string printMissingInteractionsMolblock(const t_commrec*         cr,
                                                    const gmx_reverse_top_t& rt,
                                                    const char*              moltypename,
                                                    const reverse_ilist_t&   ril,
                                                    const Range<int>&        atomRange,
                                                    const int                numAtomsPerMolecule,
                                                    const int                numMolecules,
                                                    const InteractionDefinitions& idef)
{
    const int          nril_mol = ril.index[numAtomsPerMolecule];
    std::vector<int>   isAssigned(numMolecules * nril_mol, 0);
    StringOutputStream stream;
    TextWriter         log(&stream);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rt.options()))
        {
            flagInteractionsForType(
                    ftype, idef.il[ftype], ril, atomRange, numAtomsPerMolecule, cr->dd->globalAtomIndices, isAssigned);
        }
    }

    gmx_sumi(isAssigned.size(), isAssigned.data(), cr);

    const int numMissingToPrint = 10;
    int       i                 = 0;
    for (int mol = 0; mol < numMolecules; mol++)
    {
        int j_mol = 0;
        while (j_mol < nril_mol)
        {
            int ftype = ril.il[j_mol];
            int nral  = NRAL(ftype);
            int j     = mol * nril_mol + j_mol;
            if (isAssigned[j] == 0 && !(interaction_function[ftype].flags & IF_VSITE))
            {
                if (DDMASTER(cr->dd))
                {
                    if (i == 0)
                    {
                        log.writeLineFormatted("Molecule type '%s'", moltypename);
                        log.writeLineFormatted(
                                "the first %d missing interactions, except for exclusions:",
                                numMissingToPrint);
                    }
                    log.writeStringFormatted("%20s atoms", interaction_function[ftype].longname);
                    int a = 0;
                    for (; a < nral; a++)
                    {
                        log.writeStringFormatted("%5d", ril.il[j_mol + 2 + a] + 1);
                    }
                    while (a < 4)
                    {
                        log.writeString("     ");
                        a++;
                    }
                    log.writeString(" global");
                    for (int a = 0; a < nral; a++)
                    {
                        log.writeStringFormatted("%6d",
                                                 atomRange.begin() + mol * numAtomsPerMolecule
                                                         + ril.il[j_mol + 2 + a] + 1);
                    }
                    log.ensureLineBreak();
                }
                i++;
                if (i >= numMissingToPrint)
                {
                    break;
                }
            }
            j_mol += 2 + nral_rt(ftype);
        }
    }

    return stream.toString();
}

/*! \brief Help print error output when interactions are missing */
static void printMissingInteractionsAtoms(const MDLogger&               mdlog,
                                          const t_commrec*              cr,
                                          const gmx_mtop_t&             mtop,
                                          const InteractionDefinitions& idef)
{
    const gmx_reverse_top_t& rt = *cr->dd->reverse_top;

    /* Print the atoms in the missing interactions per molblock */
    int a_end = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& moltype = mtop.moltype[molb.type];
        const int            a_start = a_end;
        a_end                        = a_start + molb.nmol * moltype.atoms.nr;
        const Range<int> atomRange(a_start, a_end);

        auto warning = printMissingInteractionsMolblock(
                cr,
                rt,
                *(moltype.name),
                cr->dd->reverse_top->interactionListForMoleculeType(molb.type),
                atomRange,
                moltype.atoms.nr,
                molb.nmol,
                idef);

        GMX_LOG(mdlog.warning).appendText(warning);
    }
}

/*! \brief Print error output when interactions are missing */
[[noreturn]] static void dd_print_missing_interactions(const MDLogger&  mdlog,
                                                       const t_commrec* cr,
                                                       const int numBondedInteractionsOverAllDomains,
                                                       const int expectedNumGlobalBondedInteractions,
                                                       const gmx_mtop_t&     top_global,
                                                       const gmx_localtop_t* top_local,
                                                       ArrayRef<const RVec>  x,
                                                       const matrix          box)
{
    int           cl[F_NRE];
    gmx_domdec_t* dd = cr->dd;

    GMX_LOG(mdlog.warning)
            .appendText(
                    "Not all bonded interactions have been properly assigned to the domain "
                    "decomposition cells");

    const int ndiff_tot = numBondedInteractionsOverAllDomains - expectedNumGlobalBondedInteractions;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        const int nral = NRAL(ftype);
        cl[ftype]      = top_local->idef.il[ftype].size() / (1 + nral);
    }

    gmx_sumi(F_NRE, cl, cr);

    if (DDMASTER(dd))
    {
        GMX_LOG(mdlog.warning).appendText("A list of missing interactions:");
        int rest_global = expectedNumGlobalBondedInteractions;
        int rest        = numBondedInteractionsOverAllDomains;
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            /* In the reverse and local top all constraints are merged
             * into F_CONSTR. So in the if statement we skip F_CONSTRNC
             * and add these constraints when doing F_CONSTR.
             */
            if (dd_check_ftype(ftype, dd->reverse_top->options()) && ftype != F_CONSTRNC)
            {
                int n = gmx_mtop_ftype_count(top_global, ftype);
                if (ftype == F_CONSTR)
                {
                    n += gmx_mtop_ftype_count(top_global, F_CONSTRNC);
                }
                int ndiff = cl[ftype] - n;
                if (ndiff != 0)
                {
                    GMX_LOG(mdlog.warning)
                            .appendTextFormatted("%20s of %6d missing %6d",
                                                 interaction_function[ftype].longname,
                                                 n,
                                                 -ndiff);
                }
                rest_global -= n;
                rest -= cl[ftype];
            }
        }

        int ndiff = rest - rest_global;
        if (ndiff != 0)
        {
            GMX_LOG(mdlog.warning).appendTextFormatted("%20s of %6d missing %6d", "exclusions", rest_global, -ndiff);
        }
    }

    printMissingInteractionsAtoms(mdlog, cr, top_global, top_local->idef);
    write_dd_pdb("dd_dump_err", 0, "dump", top_global, cr, -1, as_rvec_array(x.data()), box);

    std::string errorMessage;

    if (ndiff_tot > 0)
    {
        errorMessage =
                "One or more interactions were assigned to multiple domains of the domain "
                "decompostion. Please report this bug.";
    }
    else
    {
        errorMessage = formatString(
                "%d of the %d bonded interactions could not be calculated because some atoms "
                "involved moved further apart than the multi-body cut-off distance (%g nm) or the "
                "two-body cut-off distance (%g nm), see option -rdd, for pairs and tabulated bonds "
                "also see option -ddcheck",
                -ndiff_tot,
                expectedNumGlobalBondedInteractions,
                dd_cutoff_multibody(dd),
                dd_cutoff_twobody(dd));
    }
    gmx_fatal_collective(FARGS, cr->mpi_comm_mygroup, MASTER(cr), "%s", errorMessage.c_str());
}

/*! \brief Data to help check local topology construction
 *
 * Partitioning could incorrectly miss a bonded interaction.
 * However, checking for that requires a global communication
 * stage, which does not otherwise happen during partitioning. So,
 * for performance, we do that alongside the first global energy
 * reduction after a new DD is made. These variables handle
 * whether the check happens, its input for this domain, output
 * across all domains, and the expected value it should match. */
class LocalTopologyChecker::Impl
{
public:
    //! Constructor
    Impl(const MDLogger& mdlog, const t_commrec* cr, const gmx_mtop_t& mtop, bool useUpdateGroups);

    //! Objects used when reporting that interactions are missing
    //! {
    //! Logger
    const MDLogger& mdlog_;
    //! Communication record
    const t_commrec* cr_;
    //! Global system topology
    const gmx_mtop_t& mtop_;
    //! }

    /*! \brief Number of bonded interactions found in the local
     * topology for this domain. */
    int numBondedInteractionsToReduce_ = 0;
    /*! \brief Whether to check at the next global communication
     * stage the total number of bonded interactions found.
     *
     * Cleared after that number is found. */
    bool shouldCheckNumberOfBondedInteractions_ = false;
    /*! \brief The total number of bonded interactions found in
     * the local topology across all domains.
     *
     * Only has a value after reduction across all ranks, which is
     * removed when it is again time to check after a new
     * partition. */
    std::optional<int> numBondedInteractionsOverAllDomains_;
    //! The number of bonded interactions computed from the full system topology
    int expectedNumGlobalBondedInteractions_ = 0;
};


/*! \brief Compute the total bonded interaction count
 *
 * \param[in] mtop              The global system topology
 * \param[in] useUpdateGroups   Whether update groups are in use
 *
 * When using domain decomposition without update groups,
 * constraint-type interactions can be split across domains, and so we
 * do not consider them in this correctness check. Otherwise, we
 * include them.
 */
static int computeExpectedNumGlobalBondedInteractions(const gmx_mtop_t& mtop, const bool useUpdateGroups)
{
    int expectedNumGlobalBondedInteractions = gmx_mtop_interaction_count(mtop, IF_BOND);
    if (useUpdateGroups)
    {
        expectedNumGlobalBondedInteractions += gmx_mtop_interaction_count(mtop, IF_CONSTRAINT);
    }
    return expectedNumGlobalBondedInteractions;
}

LocalTopologyChecker::Impl::Impl(const MDLogger&   mdlog,
                                 const t_commrec*  cr,
                                 const gmx_mtop_t& mtop,
                                 const bool        useUpdateGroups) :
    mdlog_(mdlog),
    cr_(cr),
    mtop_(mtop),
    expectedNumGlobalBondedInteractions_(computeExpectedNumGlobalBondedInteractions(mtop, useUpdateGroups))
{
}

LocalTopologyChecker::LocalTopologyChecker(const MDLogger&   mdlog,
                                           const t_commrec*  cr,
                                           const gmx_mtop_t& mtop,
                                           const bool        useUpdateGroups) :
    impl_(std::make_unique<Impl>(mdlog, cr, mtop, useUpdateGroups))
{
}

LocalTopologyChecker::~LocalTopologyChecker() = default;

LocalTopologyChecker::LocalTopologyChecker(LocalTopologyChecker&&) noexcept = default;

LocalTopologyChecker& LocalTopologyChecker::operator=(LocalTopologyChecker&& other) noexcept
{
    impl_ = std::move(other.impl_);
    return *this;
}

void LocalTopologyChecker::scheduleCheckOfLocalTopology(const int numBondedInteractionsToReduce)
{
    impl_->numBondedInteractionsToReduce_ = numBondedInteractionsToReduce;
    // Note that it's possible for this to still be true from the last
    // time it was set, e.g. if repartitioning was triggered before
    // global communication that would have acted on the true
    // value. This could happen for example when replica exchange took
    // place soon after a partition.
    impl_->shouldCheckNumberOfBondedInteractions_ = true;
    // Clear the old global value, which is now invalid
    impl_->numBondedInteractionsOverAllDomains_.reset();
}

bool LocalTopologyChecker::shouldCheckNumberOfBondedInteractions() const
{
    return impl_->shouldCheckNumberOfBondedInteractions_;
}

int LocalTopologyChecker::numBondedInteractions() const
{
    return impl_->numBondedInteractionsToReduce_;
}

void LocalTopologyChecker::setNumberOfBondedInteractionsOverAllDomains(const int newValue)
{
    GMX_RELEASE_ASSERT(!impl_->numBondedInteractionsOverAllDomains_.has_value(),
                       "Cannot set number of bonded interactions because it is already set");
    impl_->numBondedInteractionsOverAllDomains_.emplace(newValue);
}

void LocalTopologyChecker::checkNumberOfBondedInteractions(const gmx_localtop_t* top_local,
                                                           ArrayRef<const RVec>  x,
                                                           const matrix          box)
{
    if (impl_->shouldCheckNumberOfBondedInteractions_)
    {
        GMX_RELEASE_ASSERT(impl_->numBondedInteractionsOverAllDomains_.has_value(),
                           "The check for the total number of bonded interactions requires the "
                           "value to have been reduced across all domains");
        if (impl_->numBondedInteractionsOverAllDomains_.value() != impl_->expectedNumGlobalBondedInteractions_)
        {
            dd_print_missing_interactions(impl_->mdlog_,
                                          impl_->cr_,
                                          impl_->numBondedInteractionsOverAllDomains_.value(),
                                          impl_->expectedNumGlobalBondedInteractions_,
                                          impl_->mtop_,
                                          top_local,
                                          x,
                                          box); // Does not return
        }
        // Now that the value is set and the check complete, future
        // global communication should not compute the value until
        // after the next partitioning.
        impl_->shouldCheckNumberOfBondedInteractions_ = false;
    }
}

} // namespace gmx
