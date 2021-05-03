/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006 - 2014, The GROMACS development team.
 * Copyright (c) 2015,2016,2017,2018,2019,2020,2021, by the GROMACS development team, led by
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
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include <cassert>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <memory>
#include <optional>
#include <string>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/domdec/options.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/mshift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringstream.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"

#include "domdec_constraints.h"
#include "domdec_internal.h"
#include "domdec_vsite.h"
#include "dump.h"

using gmx::ArrayRef;
using gmx::DDBondedChecking;
using gmx::ListOfLists;
using gmx::RVec;

/*! \brief The number of integer item in the local state, used for broadcasting of the state */
#define NITEM_DD_INIT_LOCAL_STATE 5

struct reverse_ilist_t
{
    std::vector<int> index;              /* Index for each atom into il          */
    std::vector<int> il;                 /* ftype|type|a0|...|an|ftype|...       */
    int              numAtomsInMolecule; /* The number of atoms in this molecule */
};

struct MolblockIndices
{
    int a_start;
    int a_end;
    int natoms_mol;
    int type;
};

/*! \brief Struct for thread local work data for local topology generation */
struct thread_work_t
{
    /*! \brief Constructor
     *
     * \param[in] ffparams  The interaction parameters, the lifetime of the created object should not exceed the lifetime of the passed parameters
     */
    thread_work_t(const gmx_ffparams_t& ffparams) : idef(ffparams) {}

    InteractionDefinitions         idef;               /**< Partial local topology */
    std::unique_ptr<gmx::VsitePbc> vsitePbc = nullptr; /**< vsite PBC structure */
    int              numBondedInteractions  = 0; /**< The number of bonded interactions observed */
    ListOfLists<int> excl;                       /**< List of exclusions */
};

/*! \brief Options for setting up gmx_reverse_top_t */
struct ReverseTopOptions
{
    //! Constructor, constraints and settles are not including with a single argument
    ReverseTopOptions(DDBondedChecking ddBondedChecking,
                      bool             includeConstraints = false,
                      bool             includeSettles     = false) :
        ddBondedChecking(ddBondedChecking),
        includeConstraints(includeConstraints),
        includeSettles(includeSettles)
    {
    }

    //! \brief For which bonded interactions to check assignments
    const DDBondedChecking ddBondedChecking;
    //! \brief Whether constraints are stored in this reverse top
    const bool includeConstraints;
    //! \brief Whether settles are stored in this reverse top
    const bool includeSettles;
};

/*! \brief Struct for the reverse topology: links bonded interactions to atoms */
struct gmx_reverse_top_t::Impl
{
    //! Constructs a reverse topology from \p mtop
    Impl(const gmx_mtop_t& mtop, bool useFreeEnergy, const ReverseTopOptions& reverseTopOptions);

    //! @cond Doxygen_Suppress
    //! Options for the setup of this reverse topology
    const ReverseTopOptions options;
    //! Are there interaction of type F_POSRES and/or F_FBPOSRES
    bool havePositionRestraints;
    //! \brief Are there bondeds/exclusions between atoms?
    bool bInterAtomicInteractions = false;
    //! \brief Reverse ilist for all moltypes
    std::vector<reverse_ilist_t> ril_mt;
    //! \brief The size of ril_mt[?].index summed over all entries
    int ril_mt_tot_size = 0;
    //! \brief The sorting state of bondeds for free energy
    int ilsort = ilsortUNKNOWN;
    //! \brief molblock to global atom index for quick lookup of molblocks on atom index
    std::vector<MolblockIndices> mbi;

    //! \brief Do we have intermolecular interactions?
    bool bIntermolecularInteractions = false;
    //! \brief Intermolecular reverse ilist
    reverse_ilist_t ril_intermol;

    /*! \brief Data to help check reverse topology construction
     *
     * Partitioning could incorrectly miss a bonded interaction.
     * However, checking for that requires a global communication
     * stage, which does not otherwise happen during partitioning. So,
     * for performance, we do that alongside the first global energy
     * reduction after a new DD is made. These variables handle
     * whether the check happens, its input for this domain, output
     * across all domains, and the expected value it should match. */
    /*! \{ */
    /*! \brief Number of bonded interactions found in the reverse
     * topology for this domain. */
    int numBondedInteractions = 0;
    /*! \brief Whether to check at the next global communication
     * stage the total number of bonded interactions found.
     *
     * Cleared after that number is found. */
    bool shouldCheckNumberOfBondedInteractions = false;
    /*! \brief The total number of bonded interactions found in
     * the reverse topology across all domains.
     *
     * Only has a value after reduction across all ranks, which is
     * removed when it is again time to check after a new
     * partition. */
    std::optional<int> numBondedInteractionsOverAllDomains;
    //! The number of bonded interactions computed from the full topology
    int expectedNumGlobalBondedInteractions = 0;
    /*! \} */

    /* Work data structures for multi-threading */
    //! \brief Thread work array for local topology generation
    std::vector<thread_work_t> th_work;
    //! @endcond
};


/*! \brief Returns the number of atom entries for il in gmx_reverse_top_t */
static int nral_rt(int ftype)
{
    int nral = NRAL(ftype);
    if (interaction_function[ftype].flags & IF_VSITE)
    {
        /* With vsites the reverse topology contains an extra entry
         * for storing if constructing atoms are vsites.
         */
        nral += 1;
    }

    return nral;
}

/*! \brief Return whether interactions of type \p ftype need to be assigned exactly once */
static gmx_bool dd_check_ftype(const int ftype, const ReverseTopOptions rtOptions)
{
    return ((((interaction_function[ftype].flags & IF_BOND) != 0U)
             && ((interaction_function[ftype].flags & IF_VSITE) == 0U)
             && ((rtOptions.ddBondedChecking == DDBondedChecking::All)
                 || ((interaction_function[ftype].flags & IF_LIMZERO) == 0U)))
            || (rtOptions.includeConstraints && (ftype == F_CONSTR || ftype == F_CONSTRNC))
            || (rtOptions.includeSettles && ftype == F_SETTLE));
}

/*! \brief Checks whether interactions have been assigned for one function type
 *
 * Loops over a list of interactions in the local topology of one function type
 * and flags each of the interactions as assigned in the global \p isAssigned list.
 * Exits with an inconsistency error when an interaction is assigned more than once.
 */
static void flagInteractionsForType(const int                ftype,
                                    const InteractionList&   il,
                                    const reverse_ilist_t&   ril,
                                    const gmx::Range<int>&   atomRange,
                                    const int                numAtomsPerMolecule,
                                    gmx::ArrayRef<const int> globalAtomIndices,
                                    gmx::ArrayRef<int>       isAssigned)
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
static std::string printMissingInteractionsMolblock(t_commrec*               cr,
                                                    const gmx_reverse_top_t& rt,
                                                    const char*              moltypename,
                                                    const reverse_ilist_t&   ril,
                                                    const gmx::Range<int>&   atomRange,
                                                    const int                numAtomsPerMolecule,
                                                    const int                numMolecules,
                                                    const InteractionDefinitions& idef)
{
    const int               nril_mol = ril.index[numAtomsPerMolecule];
    std::vector<int>        isAssigned(numMolecules * nril_mol, 0);
    gmx::StringOutputStream stream;
    gmx::TextWriter         log(&stream);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rt.impl_->options))
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
static void printMissingInteractionsAtoms(const gmx::MDLogger&          mdlog,
                                          t_commrec*                    cr,
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
        const gmx::Range<int> atomRange(a_start, a_end);

        auto warning = printMissingInteractionsMolblock(
                cr, rt, *(moltype.name), rt.impl_->ril_mt[molb.type], atomRange, moltype.atoms.nr, molb.nmol, idef);

        GMX_LOG(mdlog.warning).appendText(warning);
    }
}

void dd_print_missing_interactions(const gmx::MDLogger&  mdlog,
                                   t_commrec*            cr,
                                   int                   numBondedInteractionsOverAllDomains,
                                   const gmx_mtop_t&     top_global,
                                   const gmx_localtop_t* top_local,
                                   gmx::ArrayRef<const gmx::RVec> x,
                                   const matrix                   box)
{
    int           cl[F_NRE];
    gmx_domdec_t* dd = cr->dd;

    GMX_LOG(mdlog.warning)
            .appendText(
                    "Not all bonded interactions have been properly assigned to the domain "
                    "decomposition cells");

    const int ndiff_tot = numBondedInteractionsOverAllDomains
                          - dd->reverse_top->impl_->expectedNumGlobalBondedInteractions;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        const int nral = NRAL(ftype);
        cl[ftype]      = top_local->idef.il[ftype].size() / (1 + nral);
    }

    gmx_sumi(F_NRE, cl, cr);

    if (DDMASTER(dd))
    {
        GMX_LOG(mdlog.warning).appendText("A list of missing interactions:");
        int rest_global = dd->reverse_top->impl_->expectedNumGlobalBondedInteractions;
        int rest        = numBondedInteractionsOverAllDomains;
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            /* In the reverse and local top all constraints are merged
             * into F_CONSTR. So in the if statement we skip F_CONSTRNC
             * and add these constraints when doing F_CONSTR.
             */
            if (dd_check_ftype(ftype, dd->reverse_top->impl_->options) && ftype != F_CONSTRNC)
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
        errorMessage = gmx::formatString(
                "%d of the %d bonded interactions could not be calculated because some atoms "
                "involved moved further apart than the multi-body cut-off distance (%g nm) or the "
                "two-body cut-off distance (%g nm), see option -rdd, for pairs and tabulated bonds "
                "also see option -ddcheck",
                -ndiff_tot,
                dd->reverse_top->impl_->expectedNumGlobalBondedInteractions,
                dd_cutoff_multibody(dd),
                dd_cutoff_twobody(dd));
    }
    gmx_fatal_collective(FARGS, cr->mpi_comm_mygroup, MASTER(cr), "%s", errorMessage.c_str());
}

/*! \brief Return global topology molecule information for global atom index \p i_gl */
static void global_atomnr_to_moltype_ind(ArrayRef<const MolblockIndices> molblockIndices,
                                         int                             i_gl,
                                         int*                            mb,
                                         int*                            mt,
                                         int*                            mol,
                                         int*                            i_mol)
{
    const MolblockIndices* mbi   = molblockIndices.data();
    int                    start = 0;
    int                    end   = molblockIndices.size(); /* exclusive */
    int                    mid   = 0;

    /* binary search for molblock_ind */
    while (TRUE)
    {
        mid = (start + end) / 2;
        if (i_gl >= mbi[mid].a_end)
        {
            start = mid + 1;
        }
        else if (i_gl < mbi[mid].a_start)
        {
            end = mid;
        }
        else
        {
            break;
        }
    }

    *mb = mid;
    mbi += mid;

    *mt    = mbi->type;
    *mol   = (i_gl - mbi->a_start) / mbi->natoms_mol;
    *i_mol = (i_gl - mbi->a_start) - (*mol) * mbi->natoms_mol;
}

/*! \brief Returns the maximum number of exclusions per atom */
static int getMaxNumExclusionsPerAtom(const ListOfLists<int>& excls)
{
    int maxNumExcls = 0;
    for (gmx::index at = 0; at < excls.ssize(); at++)
    {
        const auto list     = excls[at];
        const int  numExcls = list.ssize();

        GMX_RELEASE_ASSERT(numExcls != 1 || list[0] == at,
                           "With 1 exclusion we expect a self-exclusion");

        maxNumExcls = std::max(maxNumExcls, numExcls);
    }

    return maxNumExcls;
}

//! Options for linking atoms in make_reverse_ilist
enum class AtomLinkRule
{
    FirstAtom,        //!< Link all interactions to the first atom in the atom list
    AllAtomsInBondeds //!< Link bonded interactions to all atoms involved, don't link vsites
};

/*! \brief Run the reverse ilist generation and store it in r_il when \p bAssign = TRUE */
static int low_make_reverse_ilist(const InteractionLists&  il_mt,
                                  const t_atom*            atom,
                                  int*                     count,
                                  const ReverseTopOptions& rtOptions,
                                  gmx::ArrayRef<const int> r_index,
                                  gmx::ArrayRef<int>       r_il,
                                  const AtomLinkRule       atomLinkRule,
                                  const bool               assignReverseIlist)
{
    const bool             includeConstraints = rtOptions.includeConstraints;
    const bool             includeSettles     = rtOptions.includeSettles;
    const DDBondedChecking ddBondedChecking   = rtOptions.ddBondedChecking;

    int nint = 0;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if ((interaction_function[ftype].flags & (IF_BOND | IF_VSITE))
            || (includeConstraints && (ftype == F_CONSTR || ftype == F_CONSTRNC))
            || (includeSettles && ftype == F_SETTLE))
        {
            const bool  isVSite = ((interaction_function[ftype].flags & IF_VSITE) != 0U);
            const int   nral    = NRAL(ftype);
            const auto& il      = il_mt[ftype];
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                const int* ia = il.iatoms.data() + i;
                // Virtual sites should not be linked for bonded interactions
                const int nlink = (atomLinkRule == AtomLinkRule::FirstAtom) ? 1 : (isVSite ? 0 : nral);
                for (int link = 0; link < nlink; link++)
                {
                    const int a = ia[1 + link];
                    if (assignReverseIlist)
                    {
                        GMX_ASSERT(!r_il.empty(), "with bAssign not allowed to be empty");
                        GMX_ASSERT(!r_index.empty(), "with bAssign not allowed to be empty");
                        r_il[r_index[a] + count[a]]     = (ftype == F_CONSTRNC ? F_CONSTR : ftype);
                        r_il[r_index[a] + count[a] + 1] = ia[0];
                        for (int j = 1; j < 1 + nral; j++)
                        {
                            /* Store the molecular atom number */
                            r_il[r_index[a] + count[a] + 1 + j] = ia[j];
                        }
                    }
                    if (interaction_function[ftype].flags & IF_VSITE)
                    {
                        if (assignReverseIlist)
                        {
                            /* Add an entry to iatoms for storing
                             * which of the constructing atoms are
                             * vsites again.
                             */
                            r_il[r_index[a] + count[a] + 2 + nral] = 0;
                            for (int j = 2; j < 1 + nral; j++)
                            {
                                if (atom[ia[j]].ptype == ParticleType::VSite)
                                {
                                    r_il[r_index[a] + count[a] + 2 + nral] |= (2 << j);
                                }
                            }
                        }
                    }
                    else
                    {
                        /* We do not count vsites since they are always
                         * uniquely assigned and can be assigned
                         * to multiple nodes with recursive vsites.
                         */
                        if (ddBondedChecking == DDBondedChecking::All
                            || !(interaction_function[ftype].flags & IF_LIMZERO))
                        {
                            nint++;
                        }
                    }
                    count[a] += 2 + nral_rt(ftype);
                }
            }
        }
    }

    return nint;
}

/*! \brief Make the reverse ilist: a list of bonded interactions linked to atoms */
static int make_reverse_ilist(const InteractionLists&  ilist,
                              const t_atoms*           atoms,
                              const ReverseTopOptions& rtOptions,
                              const AtomLinkRule       atomLinkRule,
                              reverse_ilist_t*         ril_mt)
{
    /* Count the interactions */
    const int        nat_mt = atoms->nr;
    std::vector<int> count(nat_mt);
    low_make_reverse_ilist(ilist, atoms->atom, count.data(), rtOptions, {}, {}, atomLinkRule, FALSE);

    ril_mt->index.push_back(0);
    for (int i = 0; i < nat_mt; i++)
    {
        ril_mt->index.push_back(ril_mt->index[i] + count[i]);
        count[i] = 0;
    }
    ril_mt->il.resize(ril_mt->index[nat_mt]);

    /* Store the interactions */
    int nint_mt = low_make_reverse_ilist(
            ilist, atoms->atom, count.data(), rtOptions, ril_mt->index, ril_mt->il, atomLinkRule, TRUE);

    ril_mt->numAtomsInMolecule = atoms->nr;

    return nint_mt;
}

gmx_reverse_top_t::gmx_reverse_top_t(const gmx_mtop_t&        mtop,
                                     bool                     useFreeEnergy,
                                     const ReverseTopOptions& reverseTopOptions) :
    impl_(std::make_unique<Impl>(mtop, useFreeEnergy, reverseTopOptions))
{
}

gmx_reverse_top_t::~gmx_reverse_top_t() {}

/*! \brief Generate the reverse topology */
gmx_reverse_top_t::Impl::Impl(const gmx_mtop_t&        mtop,
                              const bool               useFreeEnergy,
                              const ReverseTopOptions& reverseTopOptions) :
    options(reverseTopOptions),
    havePositionRestraints(
            gmx_mtop_ftype_count(mtop, F_POSRES) + gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0),
    bInterAtomicInteractions(mtop.bIntermolecularInteractions)
{
    bInterAtomicInteractions = mtop.bIntermolecularInteractions;
    ril_mt.resize(mtop.moltype.size());
    ril_mt_tot_size = 0;
    std::vector<int> nint_mt;
    for (size_t mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const gmx_moltype_t& molt = mtop.moltype[mt];
        if (molt.atoms.nr > 1)
        {
            bInterAtomicInteractions = true;
        }

        /* Make the atom to interaction list for this molecule type */
        int numberOfInteractions = make_reverse_ilist(
                molt.ilist, &molt.atoms, options, AtomLinkRule::FirstAtom, &ril_mt[mt]);
        nint_mt.push_back(numberOfInteractions);

        ril_mt_tot_size += ril_mt[mt].index[molt.atoms.nr];
    }
    if (debug)
    {
        fprintf(debug, "The total size of the atom to interaction index is %d integers\n", ril_mt_tot_size);
    }

    expectedNumGlobalBondedInteractions = 0;
    for (const gmx_molblock_t& molblock : mtop.molblock)
    {
        expectedNumGlobalBondedInteractions += molblock.nmol * nint_mt[molblock.type];
    }

    /* Make an intermolecular reverse top, if necessary */
    bIntermolecularInteractions = mtop.bIntermolecularInteractions;
    if (bIntermolecularInteractions)
    {
        t_atoms atoms_global;

        atoms_global.nr   = mtop.natoms;
        atoms_global.atom = nullptr; /* Only used with virtual sites */

        GMX_RELEASE_ASSERT(mtop.intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        expectedNumGlobalBondedInteractions += make_reverse_ilist(
                *mtop.intermolecular_ilist, &atoms_global, options, AtomLinkRule::FirstAtom, &ril_intermol);
    }

    if (useFreeEnergy && gmx_mtop_bondeds_free_energy(&mtop))
    {
        ilsort = ilsortFE_UNSORTED;
    }
    else
    {
        ilsort = ilsortNO_FE;
    }

    /* Make a molblock index for fast searching */
    int i = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb           = mtop.molblock[mb];
        const int             numAtomsPerMol = mtop.moltype[molb.type].atoms.nr;
        MolblockIndices       mbiMolblock;
        mbiMolblock.a_start = i;
        i += molb.nmol * numAtomsPerMol;
        mbiMolblock.a_end      = i;
        mbiMolblock.natoms_mol = numAtomsPerMol;
        mbiMolblock.type       = molb.type;
        mbi.push_back(mbiMolblock);
    }

    for (int th = 0; th < gmx_omp_nthreads_get(ModuleMultiThread::Domdec); th++)
    {
        th_work.emplace_back(mtop.ffparams);
    }
}

void dd_make_reverse_top(FILE*                           fplog,
                         gmx_domdec_t*                   dd,
                         const gmx_mtop_t&               mtop,
                         const gmx::VirtualSitesHandler* vsite,
                         const t_inputrec&               inputrec,
                         const DDBondedChecking          ddBondedChecking)
{
    if (fplog)
    {
        fprintf(fplog, "\nLinking all bonded interactions to atoms\n");
    }

    /* If normal and/or settle constraints act only within charge groups,
     * we can store them in the reverse top and simply assign them to domains.
     * Otherwise we need to assign them to multiple domains and set up
     * the parallel version constraint algorithm(s).
     */
    GMX_RELEASE_ASSERT(ddBondedChecking == DDBondedChecking::ExcludeZeroLimit
                               || ddBondedChecking == DDBondedChecking::All,
                       "Invalid enum value for mdrun -ddcheck");
    const ReverseTopOptions rtOptions(ddBondedChecking,
                                      !dd->comm->systemInfo.mayHaveSplitConstraints,
                                      !dd->comm->systemInfo.mayHaveSplitSettles);

    dd->reverse_top = std::make_unique<gmx_reverse_top_t>(
            mtop, inputrec.efep != FreeEnergyPerturbationType::No, rtOptions);

    dd->haveExclusions = false;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const int maxNumExclusionsPerAtom = getMaxNumExclusionsPerAtom(mtop.moltype[molb.type].excls);
        // We checked above that max 1 exclusion means only self exclusions
        if (maxNumExclusionsPerAtom > 1)
        {
            dd->haveExclusions = true;
        }
    }

    const int numInterUpdategroupVirtualSites =
            (vsite == nullptr ? 0 : vsite->numInterUpdategroupVirtualSites());
    if (numInterUpdategroupVirtualSites > 0)
    {
        if (fplog)
        {
            fprintf(fplog,
                    "There are %d inter update-group virtual sites,\n"
                    "will perform an extra communication step for selected coordinates and "
                    "forces\n",
                    numInterUpdategroupVirtualSites);
        }
        init_domdec_vsites(dd, numInterUpdategroupVirtualSites);
    }

    if (dd->comm->systemInfo.mayHaveSplitConstraints || dd->comm->systemInfo.mayHaveSplitSettles)
    {
        init_domdec_constraints(dd, mtop);
    }
    if (fplog)
    {
        fprintf(fplog, "\n");
    }
}

/*! \brief Store a vsite interaction at the end of \p il
 *
 * This routine is very similar to add_ifunc, but vsites interactions
 * have more work to do than other kinds of interactions, and the
 * complex way nral (and thus vector contents) depends on ftype
 * confuses static analysis tools unless we fuse the vsite
 * atom-indexing organization code with the ifunc-adding code, so that
 * they can see that nral is the same value. */
static inline void add_ifunc_for_vsites(t_iatom*           tiatoms,
                                        const gmx_ga2la_t& ga2la,
                                        int                nral,
                                        gmx_bool           bHomeA,
                                        int                a,
                                        int                a_gl,
                                        int                a_mol,
                                        const t_iatom*     iatoms,
                                        InteractionList*   il)
{
    /* Copy the type */
    tiatoms[0] = iatoms[0];

    if (bHomeA)
    {
        /* We know the local index of the first atom */
        tiatoms[1] = a;
    }
    else
    {
        /* Convert later in make_local_vsites */
        tiatoms[1] = -a_gl - 1;
    }

    GMX_ASSERT(nral >= 2 && nral <= 5, "Invalid nral for vsites");
    for (int k = 2; k < 1 + nral; k++)
    {
        int ak_gl = a_gl + iatoms[k] - a_mol;
        if (const int* homeIndex = ga2la.findHome(ak_gl))
        {
            tiatoms[k] = *homeIndex;
        }
        else
        {
            /* Copy the global index, convert later in make_local_vsites */
            tiatoms[k] = -(ak_gl + 1);
        }
        // Note that ga2la_get_home always sets the third parameter if
        // it returns TRUE
    }
    il->push_back(tiatoms[0], nral, tiatoms + 1);
}

/*! \brief Store a position restraint in idef and iatoms, complex because the parameters are different for each entry */
static void add_posres(int                     mol,
                       int                     a_mol,
                       int                     numAtomsInMolecule,
                       const gmx_molblock_t&   molb,
                       gmx::ArrayRef<int>      iatoms,
                       const t_iparams*        ip_in,
                       InteractionDefinitions* idef)
{
    /* This position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    const int n = idef->il[F_POSRES].size() / 2;

    /* Get the position restraint coordinates from the molblock */
    const int a_molb = mol * numAtomsInMolecule + a_mol;
    GMX_ASSERT(a_molb < ssize(molb.posres_xA),
               "We need a sufficient number of position restraint coordinates");
    /* Copy the force constants */
    t_iparams ip        = ip_in[iatoms[0]];
    ip.posres.pos0A[XX] = molb.posres_xA[a_molb][XX];
    ip.posres.pos0A[YY] = molb.posres_xA[a_molb][YY];
    ip.posres.pos0A[ZZ] = molb.posres_xA[a_molb][ZZ];
    if (!molb.posres_xB.empty())
    {
        ip.posres.pos0B[XX] = molb.posres_xB[a_molb][XX];
        ip.posres.pos0B[YY] = molb.posres_xB[a_molb][YY];
        ip.posres.pos0B[ZZ] = molb.posres_xB[a_molb][ZZ];
    }
    else
    {
        ip.posres.pos0B[XX] = ip.posres.pos0A[XX];
        ip.posres.pos0B[YY] = ip.posres.pos0A[YY];
        ip.posres.pos0B[ZZ] = ip.posres.pos0A[ZZ];
    }
    /* Set the parameter index for idef->iparams_posres */
    iatoms[0] = n;
    idef->iparams_posres.push_back(ip);

    GMX_ASSERT(int(idef->iparams_posres.size()) == n + 1,
               "The index of the parameter type should match n");
}

/*! \brief Store a flat-bottomed position restraint in idef and iatoms, complex because the parameters are different for each entry */
static void add_fbposres(int                     mol,
                         int                     a_mol,
                         int                     numAtomsInMolecule,
                         const gmx_molblock_t&   molb,
                         gmx::ArrayRef<int>      iatoms,
                         const t_iparams*        ip_in,
                         InteractionDefinitions* idef)
{
    /* This flat-bottom position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    const int n = idef->il[F_FBPOSRES].size() / 2;

    /* Get the position restraint coordinats from the molblock */
    const int a_molb = mol * numAtomsInMolecule + a_mol;
    GMX_ASSERT(a_molb < ssize(molb.posres_xA),
               "We need a sufficient number of position restraint coordinates");
    /* Copy the force constants */
    t_iparams ip = ip_in[iatoms[0]];
    /* Take reference positions from A position of normal posres */
    ip.fbposres.pos0[XX] = molb.posres_xA[a_molb][XX];
    ip.fbposres.pos0[YY] = molb.posres_xA[a_molb][YY];
    ip.fbposres.pos0[ZZ] = molb.posres_xA[a_molb][ZZ];

    /* Note: no B-type for flat-bottom posres */

    /* Set the parameter index for idef->iparams_fbposres */
    iatoms[0] = n;
    idef->iparams_fbposres.push_back(ip);

    GMX_ASSERT(int(idef->iparams_fbposres.size()) == n + 1,
               "The index of the parameter type should match n");
}

/*! \brief Store a virtual site interaction, complex because of PBC and recursion */
static void add_vsite(const gmx_ga2la_t&       ga2la,
                      gmx::ArrayRef<const int> index,
                      gmx::ArrayRef<const int> rtil,
                      int                      ftype,
                      int                      nral,
                      gmx_bool                 bHomeA,
                      int                      a,
                      int                      a_gl,
                      int                      a_mol,
                      const t_iatom*           iatoms,
                      InteractionDefinitions*  idef)
{
    t_iatom tiatoms[1 + MAXATOMLIST];

    /* Add this interaction to the local topology */
    add_ifunc_for_vsites(tiatoms, ga2la, nral, bHomeA, a, a_gl, a_mol, iatoms, &idef->il[ftype]);

    if (iatoms[1 + nral])
    {
        /* Check for recursion */
        for (int k = 2; k < 1 + nral; k++)
        {
            if ((iatoms[1 + nral] & (2 << k)) && (tiatoms[k] < 0))
            {
                /* This construction atoms is a vsite and not a home atom */
                if (gmx_debug_at)
                {
                    fprintf(debug,
                            "Constructing atom %d of vsite atom %d is a vsite and non-home\n",
                            iatoms[k] + 1,
                            a_mol + 1);
                }
                /* Find the vsite construction */

                /* Check all interactions assigned to this atom */
                int j = index[iatoms[k]];
                while (j < index[iatoms[k] + 1])
                {
                    int ftype_r = rtil[j++];
                    int nral_r  = NRAL(ftype_r);
                    if (interaction_function[ftype_r].flags & IF_VSITE)
                    {
                        /* Add this vsite (recursion) */
                        add_vsite(ga2la,
                                  index,
                                  rtil,
                                  ftype_r,
                                  nral_r,
                                  FALSE,
                                  -1,
                                  a_gl + iatoms[k] - iatoms[1],
                                  iatoms[k],
                                  rtil.data() + j,
                                  idef);
                    }
                    j += 1 + nral_rt(ftype_r);
                }
            }
        }
    }
}

/*! \brief Returns the squared distance between atoms \p i and \p j */
static real dd_dist2(const t_pbc* pbc_null, ArrayRef<const RVec> coordinates, const int i, const int j)
{
    rvec dx;

    if (pbc_null)
    {
        pbc_dx_aiuc(pbc_null, coordinates[i], coordinates[j], dx);
    }
    else
    {
        rvec_sub(coordinates[i], coordinates[j], dx);
    }

    return norm2(dx);
}

/*! \brief Append t_idef structures 1 to nsrc in src to *dest */
static void combine_idef(InteractionDefinitions* dest, gmx::ArrayRef<const thread_work_t> src)
{
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        int n = 0;
        for (gmx::index s = 1; s < src.ssize(); s++)
        {
            n += src[s].idef.il[ftype].size();
        }
        if (n > 0)
        {
            for (gmx::index s = 1; s < src.ssize(); s++)
            {
                dest->il[ftype].append(src[s].idef.il[ftype]);
            }

            /* Position restraints need an additional treatment */
            if (ftype == F_POSRES || ftype == F_FBPOSRES)
            {
                int                     nposres = dest->il[ftype].size() / 2;
                std::vector<t_iparams>& iparams_dest =
                        (ftype == F_POSRES ? dest->iparams_posres : dest->iparams_fbposres);

                /* Set nposres to the number of original position restraints in dest */
                for (gmx::index s = 1; s < src.ssize(); s++)
                {
                    nposres -= src[s].idef.il[ftype].size() / 2;
                }

                for (gmx::index s = 1; s < src.ssize(); s++)
                {
                    const std::vector<t_iparams>& iparams_src =
                            (ftype == F_POSRES ? src[s].idef.iparams_posres : src[s].idef.iparams_fbposres);
                    iparams_dest.insert(iparams_dest.end(), iparams_src.begin(), iparams_src.end());

                    /* Correct the indices into iparams_posres */
                    for (int i = 0; i < src[s].idef.il[ftype].size() / 2; i++)
                    {
                        /* Correct the index into iparams_posres */
                        dest->il[ftype].iatoms[nposres * 2] = nposres;
                        nposres++;
                    }
                }
                GMX_RELEASE_ASSERT(
                        int(iparams_dest.size()) == nposres,
                        "The number of parameters should match the number of restraints");
            }
        }
    }
}

//! Options for assigning interactions for atoms
enum class InteractionConnectivity
{
    Intramolecular, //!< Only intra-molecular interactions
    Intermolecular  //!< Only inter-molecular interactions
};

/*! \brief Determine whether the local domain has responsibility for
 * any of the bonded interactions for local atom \p atomIndex
 * and assign those to the local domain.
 *
 * \returns The total number of bonded interactions for this atom for
 * which this domain is responsible.
 */
template<InteractionConnectivity interactionConnectivity>
static inline int assignInteractionsForAtom(const int                 atomIndex,
                                            const int                 globalAtomIndex,
                                            const int                 atomIndexInMolecule,
                                            gmx::ArrayRef<const int>  index,
                                            gmx::ArrayRef<const int>  rtil,
                                            const int                 ind_start,
                                            const int                 ind_end,
                                            const gmx_ga2la_t&        ga2la,
                                            const gmx_domdec_zones_t* zones,
                                            const bool                checkDistanceMultiBody,
                                            const ivec                rcheck,
                                            const bool                checkDistanceTwoBody,
                                            const real                cutoffSquared,
                                            const t_pbc*              pbc_null,
                                            ArrayRef<const RVec>      coordinates,
                                            InteractionDefinitions*   idef,
                                            const int                 iz,
                                            const DDBondedChecking    ddBondedChecking)
{
    gmx::ArrayRef<const DDPairInteractionRanges> iZones = zones->iZones;

    int numBondedInteractions = 0;

    int j = ind_start;
    while (j < ind_end)
    {
        t_iatom tiatoms[1 + MAXATOMLIST];

        const int ftype  = rtil[j++];
        auto      iatoms = gmx::constArrayRefFromArray(rtil.data() + j, rtil.size() - j);
        const int nral   = NRAL(ftype);
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            GMX_ASSERT(interactionConnectivity == InteractionConnectivity::Intramolecular,
                       "All vsites should be intramolecular");

            /* The vsite construction goes where the vsite itself is */
            if (iz == 0)
            {
                add_vsite(ga2la,
                          index,
                          rtil,
                          ftype,
                          nral,
                          TRUE,
                          atomIndex,
                          globalAtomIndex,
                          atomIndexInMolecule,
                          iatoms.data(),
                          idef);
            }
        }
        else
        {
            bool bUse = false;

            /* Copy the type */
            tiatoms[0] = iatoms[0];

            if (nral == 1)
            {
                GMX_ASSERT(interactionConnectivity == InteractionConnectivity::Intramolecular,
                           "All interactions that involve a single atom are intramolecular");

                /* Assign single-body interactions to the home zone.
                 * Position restraints are not handled here, but separately.
                 */
                if (iz == 0 && !(ftype == F_POSRES || ftype == F_FBPOSRES))
                {
                    bUse       = true;
                    tiatoms[1] = atomIndex;
                }
            }
            else if (nral == 2)
            {
                /* This is a two-body interaction, we can assign
                 * analogous to the non-bonded assignments.
                 */
                const int k_gl = (interactionConnectivity == InteractionConnectivity::Intramolecular)
                                         ?
                                         /* Get the global index using the offset in the molecule */
                                         (globalAtomIndex + iatoms[2] - atomIndexInMolecule)
                                         : iatoms[2];
                if (const auto* entry = ga2la.find(k_gl))
                {
                    int kz = entry->cell;
                    if (kz >= zones->n)
                    {
                        kz -= zones->n;
                    }
                    /* Check zone interaction assignments */
                    bUse = ((iz < iZones.ssize() && iz <= kz && iZones[iz].jZoneRange.isInRange(kz))
                            || (kz < iZones.ssize() && iz > kz && iZones[kz].jZoneRange.isInRange(iz)));
                    if (bUse)
                    {
                        GMX_ASSERT(ftype != F_CONSTR || (iz == 0 && kz == 0),
                                   "Constraint assigned here should only involve home atoms");

                        tiatoms[1] = atomIndex;
                        tiatoms[2] = entry->la;
                        /* If necessary check the cgcm distance */
                        if (checkDistanceTwoBody
                            && dd_dist2(pbc_null, coordinates, tiatoms[1], tiatoms[2]) >= cutoffSquared)
                        {
                            bUse = false;
                        }
                    }
                }
                else
                {
                    bUse = false;
                }
            }
            else
            {
                /* Assign this multi-body bonded interaction to
                 * the local domain if we have all the atoms involved
                 * (local or communicated) and the minimum zone shift
                 * in each dimension is zero, for dimensions
                 * with 2 DD cells an extra check may be necessary.
                 */
                ivec k_zero, k_plus;
                bUse = true;
                clear_ivec(k_zero);
                clear_ivec(k_plus);
                for (int k = 1; k <= nral && bUse; k++)
                {
                    const int k_gl =
                            (interactionConnectivity == InteractionConnectivity::Intramolecular)
                                    ?
                                    /* Get the global index using the offset in the molecule */
                                    (globalAtomIndex + iatoms[k] - atomIndexInMolecule)
                                    : iatoms[k];
                    const auto* entry = ga2la.find(k_gl);
                    if (entry == nullptr || entry->cell >= zones->n)
                    {
                        /* We do not have this atom of this interaction
                         * locally, or it comes from more than one cell
                         * away.
                         */
                        bUse = FALSE;
                    }
                    else
                    {
                        tiatoms[k] = entry->la;
                        for (int d = 0; d < DIM; d++)
                        {
                            if (zones->shift[entry->cell][d] == 0)
                            {
                                k_zero[d] = k;
                            }
                            else
                            {
                                k_plus[d] = k;
                            }
                        }
                    }
                }
                bUse = (bUse && (k_zero[XX] != 0) && (k_zero[YY] != 0) && (k_zero[ZZ] != 0));
                if (checkDistanceMultiBody)
                {
                    for (int d = 0; (d < DIM && bUse); d++)
                    {
                        /* Check if the distance falls within
                         * the cut-off to avoid possible multiple
                         * assignments of bonded interactions.
                         */
                        if (rcheck[d] && k_plus[d]
                            && dd_dist2(pbc_null, coordinates, tiatoms[k_zero[d]], tiatoms[k_plus[d]])
                                       >= cutoffSquared)
                        {
                            bUse = FALSE;
                        }
                    }
                }
            }
            if (bUse)
            {
                /* Add this interaction to the local topology */
                idef->il[ftype].push_back(tiatoms[0], nral, tiatoms + 1);
                /* Sum so we can check in global_stat
                 * if we have everything.
                 */
                if (ddBondedChecking == DDBondedChecking::All
                    || !(interaction_function[ftype].flags & IF_LIMZERO))
                {
                    numBondedInteractions++;
                }
            }
        }
        j += 1 + nral_rt(ftype);
    }

    return numBondedInteractions;
}

/*! \brief Determine whether the local domain has responsibility for
 * any of the position restraints for local atom \p atomIndex
 * and assign those to the local domain.
 *
 * \returns The total number of bonded interactions for this atom for
 * which this domain is responsible.
 */
static inline int assignPositionRestraintsForAtom(const int                atomIndex,
                                                  const int                mol,
                                                  const int                atomIndexInMolecule,
                                                  const int                numAtomsInMolecule,
                                                  gmx::ArrayRef<const int> rtil,
                                                  const int                ind_start,
                                                  const int                ind_end,
                                                  const gmx_molblock_t&    molb,
                                                  const t_iparams*         ip_in,
                                                  InteractionDefinitions*  idef)
{
    constexpr int nral = 1;

    int numBondedInteractions = 0;

    int j = ind_start;
    while (j < ind_end)
    {
        const int ftype  = rtil[j++];
        auto      iatoms = gmx::constArrayRefFromArray(rtil.data() + j, rtil.size() - j);

        if (ftype == F_POSRES || ftype == F_FBPOSRES)
        {
            std::array<int, 1 + nral> tiatoms = { iatoms[0], atomIndex };
            if (ftype == F_POSRES)
            {
                add_posres(mol, atomIndexInMolecule, numAtomsInMolecule, molb, tiatoms, ip_in, idef);
            }
            else
            {
                add_fbposres(mol, atomIndexInMolecule, numAtomsInMolecule, molb, tiatoms, ip_in, idef);
            }
            idef->il[ftype].push_back(tiatoms[0], nral, tiatoms.data() + 1);
            numBondedInteractions++;
        }
        j += 1 + nral_rt(ftype);
    }

    return numBondedInteractions;
}

/*! \brief This function looks up and assigns bonded interactions for zone iz.
 *
 * With thread parallelizing each thread acts on a different atom range:
 * at_start to at_end.
 */
static int make_bondeds_zone(gmx_reverse_top_t*                 rt,
                             ArrayRef<const int>                globalAtomIndices,
                             const gmx_ga2la_t&                 ga2la,
                             const gmx_domdec_zones_t*          zones,
                             const std::vector<gmx_molblock_t>& molb,
                             const bool                         checkDistanceMultiBody,
                             const ivec                         rcheck,
                             const bool                         checkDistanceTwoBody,
                             const real                         cutoffSquared,
                             const t_pbc*                       pbc_null,
                             ArrayRef<const RVec>               coordinates,
                             const t_iparams*                   ip_in,
                             InteractionDefinitions*            idef,
                             int                                izone,
                             const gmx::Range<int>&             atomRange)
{
    int mb                  = 0;
    int mt                  = 0;
    int mol                 = 0;
    int atomIndexInMolecule = 0;

    const auto ddBondedChecking = rt->impl_->options.ddBondedChecking;

    int numBondedInteractions = 0;

    for (int i : atomRange)
    {
        /* Get the global atom number */
        const int globalAtomIndex = globalAtomIndices[i];
        global_atomnr_to_moltype_ind(rt->impl_->mbi, globalAtomIndex, &mb, &mt, &mol, &atomIndexInMolecule);
        /* Check all intramolecular interactions assigned to this atom */
        gmx::ArrayRef<const int>     index = rt->impl_->ril_mt[mt].index;
        gmx::ArrayRef<const t_iatom> rtil  = rt->impl_->ril_mt[mt].il;

        numBondedInteractions += assignInteractionsForAtom<InteractionConnectivity::Intramolecular>(
                i,
                globalAtomIndex,
                atomIndexInMolecule,
                index,
                rtil,
                index[atomIndexInMolecule],
                index[atomIndexInMolecule + 1],
                ga2la,
                zones,
                checkDistanceMultiBody,
                rcheck,
                checkDistanceTwoBody,
                cutoffSquared,
                pbc_null,
                coordinates,
                idef,
                izone,
                ddBondedChecking);

        // Assign position restraints, when present, for the home zone
        if (izone == 0 && rt->impl_->havePositionRestraints)
        {
            numBondedInteractions +=
                    assignPositionRestraintsForAtom(i,
                                                    mol,
                                                    atomIndexInMolecule,
                                                    rt->impl_->ril_mt[mt].numAtomsInMolecule,
                                                    rtil,
                                                    index[atomIndexInMolecule],
                                                    index[atomIndexInMolecule + 1],
                                                    molb[mb],
                                                    ip_in,
                                                    idef);
        }

        if (rt->impl_->bIntermolecularInteractions)
        {
            /* Check all intermolecular interactions assigned to this atom */
            index = rt->impl_->ril_intermol.index;
            rtil  = rt->impl_->ril_intermol.il;

            numBondedInteractions += assignInteractionsForAtom<InteractionConnectivity::Intermolecular>(
                    i,
                    -1, // not used
                    -1, // not used
                    index,
                    rtil,
                    index[globalAtomIndex],
                    index[globalAtomIndex + 1],
                    ga2la,
                    zones,
                    checkDistanceMultiBody,
                    rcheck,
                    checkDistanceTwoBody,
                    cutoffSquared,
                    pbc_null,
                    coordinates,
                    idef,
                    izone,
                    ddBondedChecking);
        }
    }

    return numBondedInteractions;
}

/*! \brief Set the exclusion data for i-zone \p iz */
static void make_exclusions_zone(ArrayRef<const int>               globalAtomIndices,
                                 const gmx_ga2la_t&                ga2la,
                                 gmx_domdec_zones_t*               zones,
                                 ArrayRef<const MolblockIndices>   molblockIndices,
                                 const std::vector<gmx_moltype_t>& moltype,
                                 const int*                        cginfo,
                                 ListOfLists<int>*                 lexcls,
                                 int                               iz,
                                 int                               at_start,
                                 int                               at_end,
                                 const gmx::ArrayRef<const int>    intermolecularExclusionGroup)
{
    const auto& jAtomRange = zones->iZones[iz].jAtomRange;

    const gmx::index oldNumLists = lexcls->ssize();

    std::vector<int> exclusionsForAtom;
    for (int at = at_start; at < at_end; at++)
    {
        exclusionsForAtom.clear();

        if (GET_CGINFO_EXCL_INTER(cginfo[at]))
        {
            int mb    = 0;
            int mt    = 0;
            int mol   = 0;
            int a_mol = 0;

            /* Copy the exclusions from the global top */
            int a_gl = globalAtomIndices[at];
            global_atomnr_to_moltype_ind(molblockIndices, a_gl, &mb, &mt, &mol, &a_mol);
            const auto excls = moltype[mt].excls[a_mol];
            for (const int aj_mol : excls)
            {
                if (const auto* jEntry = ga2la.find(a_gl + aj_mol - a_mol))
                {
                    /* This check is not necessary, but it can reduce
                     * the number of exclusions in the list, which in turn
                     * can speed up the pair list construction a bit.
                     */
                    if (jAtomRange.isInRange(jEntry->la))
                    {
                        exclusionsForAtom.push_back(jEntry->la);
                    }
                }
            }
        }

        bool isExcludedAtom = !intermolecularExclusionGroup.empty()
                              && std::find(intermolecularExclusionGroup.begin(),
                                           intermolecularExclusionGroup.end(),
                                           globalAtomIndices[at])
                                         != intermolecularExclusionGroup.end();

        if (isExcludedAtom)
        {
            for (int qmAtomGlobalIndex : intermolecularExclusionGroup)
            {
                if (const auto* entry = ga2la.find(qmAtomGlobalIndex))
                {
                    exclusionsForAtom.push_back(entry->la);
                }
            }
        }

        /* Append the exclusions for this atom to the topology */
        lexcls->pushBack(exclusionsForAtom);
    }

    GMX_RELEASE_ASSERT(
            lexcls->ssize() - oldNumLists == at_end - at_start,
            "The number of exclusion list should match the number of atoms in the range");
}

/*! \brief Generate and store all required local bonded interactions in \p idef and local exclusions in \p lexcls */
static void make_local_bondeds_excls(gmx_domdec_t*           dd,
                                     gmx_domdec_zones_t*     zones,
                                     const gmx_mtop_t&       mtop,
                                     const int*              cginfo,
                                     const bool              checkDistanceMultiBody,
                                     const ivec              rcheck,
                                     const gmx_bool          checkDistanceTwoBody,
                                     const real              cutoff,
                                     const t_pbc*            pbc_null,
                                     ArrayRef<const RVec>    coordinates,
                                     InteractionDefinitions* idef,
                                     ListOfLists<int>*       lexcls)
{
    int nzone_bondeds = 0;

    if (dd->reverse_top->impl_->bInterAtomicInteractions)
    {
        nzone_bondeds = zones->n;
    }
    else
    {
        /* Only single charge group (or atom) molecules, so interactions don't
         * cross zone boundaries and we only need to assign in the home zone.
         */
        nzone_bondeds = 1;
    }

    /* We only use exclusions from i-zones to i- and j-zones */
    const int numIZonesForExclusions = (dd->haveExclusions ? zones->iZones.size() : 0);

    gmx_reverse_top_t* rt = dd->reverse_top.get();

    const real cutoffSquared = gmx::square(cutoff);

    /* Clear the counts */
    idef->clear();
    dd->reverse_top->impl_->numBondedInteractions = 0;

    lexcls->clear();

    for (int izone = 0; izone < nzone_bondeds; izone++)
    {
        const int cg0 = zones->cg_range[izone];
        const int cg1 = zones->cg_range[izone + 1];

        const int numThreads = rt->impl_->th_work.size();
#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int thread = 0; thread < numThreads; thread++)
        {
            try
            {
                InteractionDefinitions* idef_t = nullptr;

                int cg0t = cg0 + ((cg1 - cg0) * thread) / numThreads;
                int cg1t = cg0 + ((cg1 - cg0) * (thread + 1)) / numThreads;

                if (thread == 0)
                {
                    idef_t = idef;
                }
                else
                {
                    idef_t = &rt->impl_->th_work[thread].idef;
                    idef_t->clear();
                }

                rt->impl_->th_work[thread].numBondedInteractions =
                        make_bondeds_zone(rt,
                                          dd->globalAtomIndices,
                                          *dd->ga2la,
                                          zones,
                                          mtop.molblock,
                                          checkDistanceMultiBody,
                                          rcheck,
                                          checkDistanceTwoBody,
                                          cutoffSquared,
                                          pbc_null,
                                          coordinates,
                                          idef->iparams.data(),
                                          idef_t,
                                          izone,
                                          gmx::Range<int>(cg0t, cg1t));

                if (izone < numIZonesForExclusions)
                {
                    ListOfLists<int>* excl_t = nullptr;
                    if (thread == 0)
                    {
                        // Thread 0 stores exclusions directly in the final storage
                        excl_t = lexcls;
                    }
                    else
                    {
                        // Threads > 0 store in temporary storage, starting at list index 0
                        excl_t = &rt->impl_->th_work[thread].excl;
                        excl_t->clear();
                    }

                    /* No charge groups and no distance check required */
                    make_exclusions_zone(dd->globalAtomIndices,
                                         *dd->ga2la,
                                         zones,
                                         rt->impl_->mbi,
                                         mtop.moltype,
                                         cginfo,
                                         excl_t,
                                         izone,
                                         cg0t,
                                         cg1t,
                                         mtop.intermolecularExclusionGroup);
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }

        if (rt->impl_->th_work.size() > 1)
        {
            combine_idef(idef, rt->impl_->th_work);
        }

        for (const thread_work_t& th_work : rt->impl_->th_work)
        {
            dd->reverse_top->impl_->numBondedInteractions += th_work.numBondedInteractions;
        }

        if (izone < numIZonesForExclusions)
        {
            for (std::size_t th = 1; th < rt->impl_->th_work.size(); th++)
            {
                lexcls->appendListOfLists(rt->impl_->th_work[th].excl);
            }
        }
    }

    // Note that it's possible for this to still be true from the last
    // time it was set, e.g. if repartitioning was triggered before
    // global communication that would have acted on the true
    // value. This could happen for example when replica exchange took
    // place soon after a partition.
    dd->reverse_top->impl_->shouldCheckNumberOfBondedInteractions = true;
    // Clear the old global value, which is now invalid
    dd->reverse_top->impl_->numBondedInteractionsOverAllDomains.reset();

    if (debug)
    {
        fprintf(debug, "We have %d exclusions\n", lexcls->numElements());
    }
}

bool shouldCheckNumberOfBondedInteractions(const gmx_domdec_t& dd)
{
    return dd.reverse_top->impl_->shouldCheckNumberOfBondedInteractions;
}

int numBondedInteractions(const gmx_domdec_t& dd)
{
    return dd.reverse_top->impl_->numBondedInteractions;
}

void setNumberOfBondedInteractionsOverAllDomains(const gmx_domdec_t& dd, int newValue)
{
    GMX_RELEASE_ASSERT(!dd.reverse_top->impl_->numBondedInteractionsOverAllDomains.has_value(),
                       "Cannot set number of bonded interactions because it is already set");
    dd.reverse_top->impl_->numBondedInteractionsOverAllDomains.emplace(newValue);
}

void checkNumberOfBondedInteractions(const gmx::MDLogger&           mdlog,
                                     t_commrec*                     cr,
                                     const gmx_mtop_t&              top_global,
                                     const gmx_localtop_t*          top_local,
                                     gmx::ArrayRef<const gmx::RVec> x,
                                     const matrix                   box)
{
    GMX_RELEASE_ASSERT(
            DOMAINDECOMP(cr),
            "No need to check number of bonded interactions when not using domain decomposition");
    if (cr->dd->reverse_top->impl_->shouldCheckNumberOfBondedInteractions)
    {
        GMX_RELEASE_ASSERT(cr->dd->reverse_top->impl_->numBondedInteractionsOverAllDomains.has_value(),
                           "The check for the total number of bonded interactions requires the "
                           "value to have been reduced across all domains");
        if (cr->dd->reverse_top->impl_->numBondedInteractionsOverAllDomains.value()
            != cr->dd->reverse_top->impl_->expectedNumGlobalBondedInteractions)
        {
            dd_print_missing_interactions(
                    mdlog,
                    cr,
                    cr->dd->reverse_top->impl_->numBondedInteractionsOverAllDomains.value(),
                    top_global,
                    top_local,
                    x,
                    box); // Does not return
        }
        // Now that the value is set and the check complete, future
        // global communication should not compute the value until
        // after the next partitioning.
        cr->dd->reverse_top->impl_->shouldCheckNumberOfBondedInteractions = false;
    }
}

void dd_make_local_top(gmx_domdec_t*        dd,
                       gmx_domdec_zones_t*  zones,
                       int                  npbcdim,
                       matrix               box,
                       rvec                 cellsize_min,
                       const ivec           npulse,
                       t_forcerec*          fr,
                       ArrayRef<const RVec> coordinates,
                       const gmx_mtop_t&    mtop,
                       gmx_localtop_t*      ltop)
{
    real  rc = -1;
    ivec  rcheck;
    t_pbc pbc, *pbc_null = nullptr;

    if (debug)
    {
        fprintf(debug, "Making local topology\n");
    }

    bool checkDistanceMultiBody = false;
    bool checkDistanceTwoBody   = false;

    if (dd->reverse_top->impl_->bInterAtomicInteractions)
    {
        /* We need to check to which cell bondeds should be assigned */
        rc = dd_cutoff_twobody(dd);
        if (debug)
        {
            fprintf(debug, "Two-body bonded cut-off distance is %g\n", rc);
        }

        /* Should we check distances when assigning bonded interactions? */
        for (int d = 0; d < DIM; d++)
        {
            rcheck[d] = FALSE;
            /* Only need to check for dimensions where the part of the box
             * that is not communicated is smaller than the cut-off.
             */
            if (d < npbcdim && dd->numCells[d] > 1
                && (dd->numCells[d] - npulse[d]) * cellsize_min[d] < 2 * rc)
            {
                if (dd->numCells[d] == 2)
                {
                    rcheck[d]              = TRUE;
                    checkDistanceMultiBody = TRUE;
                }
                /* Check for interactions between two atoms,
                 * where we can allow interactions up to the cut-off,
                 * instead of up to the smallest cell dimension.
                 */
                checkDistanceTwoBody = TRUE;
            }
            if (debug)
            {
                fprintf(debug,
                        "dim %d cellmin %f bonded rcheck[%d] = %d, checkDistanceTwoBody = %s\n",
                        d,
                        cellsize_min[d],
                        d,
                        rcheck[d],
                        gmx::boolToString(checkDistanceTwoBody));
            }
        }
        if (checkDistanceMultiBody || checkDistanceTwoBody)
        {
            if (fr->bMolPBC)
            {
                pbc_null = set_pbc_dd(&pbc, fr->pbcType, dd->numCells, TRUE, box);
            }
            else
            {
                pbc_null = nullptr;
            }
        }
    }

    make_local_bondeds_excls(dd,
                             zones,
                             mtop,
                             fr->cginfo.data(),
                             checkDistanceMultiBody,
                             rcheck,
                             checkDistanceTwoBody,
                             rc,
                             pbc_null,
                             coordinates,
                             &ltop->idef,
                             &ltop->excls);

    /* The ilist is not sorted yet,
     * we can only do this when we have the charge arrays.
     */
    ltop->idef.ilsort = ilsortUNKNOWN;
}

void dd_sort_local_top(const gmx_domdec_t& dd, const t_mdatoms* mdatoms, gmx_localtop_t* ltop)
{
    if (dd.reverse_top->impl_->ilsort == ilsortNO_FE)
    {
        ltop->idef.ilsort = ilsortNO_FE;
    }
    else
    {
        gmx_sort_ilist_fe(&ltop->idef, mdatoms->chargeA, mdatoms->chargeB);
    }
}

void dd_init_local_state(const gmx_domdec_t& dd, const t_state* state_global, t_state* state_local)
{
    int buf[NITEM_DD_INIT_LOCAL_STATE];

    if (DDMASTER(dd))
    {
        buf[0] = state_global->flags;
        buf[1] = state_global->ngtc;
        buf[2] = state_global->nnhpres;
        buf[3] = state_global->nhchainlength;
        buf[4] = state_global->dfhist ? state_global->dfhist->nlambda : 0;
    }
    dd_bcast(&dd, NITEM_DD_INIT_LOCAL_STATE * sizeof(int), buf);

    init_gtc_state(state_local, buf[1], buf[2], buf[3]);
    init_dfhist_state(state_local, buf[4]);
    state_local->flags = buf[0];
}

/*! \brief Check if a link is stored in \p link between charge groups \p cg_gl and \p cg_gl_j and if not so, store a link */
static void check_link(t_blocka* link, int cg_gl, int cg_gl_j)
{
    bool bFound = false;
    for (int k = link->index[cg_gl]; k < link->index[cg_gl + 1]; k++)
    {
        GMX_RELEASE_ASSERT(link->a, "Inconsistent NULL pointer while making charge-group links");
        if (link->a[k] == cg_gl_j)
        {
            bFound = TRUE;
        }
    }
    if (!bFound)
    {
        GMX_RELEASE_ASSERT(link->a || link->index[cg_gl + 1] + 1 > link->nalloc_a,
                           "Inconsistent allocation of link");
        /* Add this charge group link */
        if (link->index[cg_gl + 1] + 1 > link->nalloc_a)
        {
            link->nalloc_a = over_alloc_large(link->index[cg_gl + 1] + 1);
            srenew(link->a, link->nalloc_a);
        }
        link->a[link->index[cg_gl + 1]] = cg_gl_j;
        link->index[cg_gl + 1]++;
    }
}

t_blocka* makeBondedLinks(const gmx_mtop_t& mtop, gmx::ArrayRef<cginfo_mb_t> cginfo_mb)
{
    t_blocka* link = nullptr;

    /* For each atom make a list of other atoms in the system
     * that a linked to it via bonded interactions
     * which are also stored in reverse_top.
     */

    reverse_ilist_t ril_intermol;
    if (mtop.bIntermolecularInteractions)
    {
        t_atoms atoms;

        atoms.nr   = mtop.natoms;
        atoms.atom = nullptr;

        GMX_RELEASE_ASSERT(mtop.intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        ReverseTopOptions rtOptions(DDBondedChecking::ExcludeZeroLimit);
        make_reverse_ilist(
                *mtop.intermolecular_ilist, &atoms, rtOptions, AtomLinkRule::AllAtomsInBondeds, &ril_intermol);
    }

    snew(link, 1);
    snew(link->index, mtop.natoms + 1);
    link->nalloc_a = 0;
    link->a        = nullptr;

    link->index[0] = 0;
    int cg_offset  = 0;
    int ncgi       = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        if (molb.nmol == 0)
        {
            continue;
        }
        const gmx_moltype_t& molt = mtop.moltype[molb.type];
        /* Make a reverse ilist in which the interactions are linked
         * to all atoms, not only the first atom as in gmx_reverse_top.
         * The constraints are discarded here.
         */
        ReverseTopOptions rtOptions(DDBondedChecking::ExcludeZeroLimit);
        reverse_ilist_t   ril;
        make_reverse_ilist(molt.ilist, &molt.atoms, rtOptions, AtomLinkRule::AllAtomsInBondeds, &ril);

        cginfo_mb_t* cgi_mb = &cginfo_mb[mb];

        int mol = 0;
        for (mol = 0; mol < (mtop.bIntermolecularInteractions ? molb.nmol : 1); mol++)
        {
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                int cg_gl              = cg_offset + a;
                link->index[cg_gl + 1] = link->index[cg_gl];
                int i                  = ril.index[a];
                while (i < ril.index[a + 1])
                {
                    int ftype = ril.il[i++];
                    int nral  = NRAL(ftype);
                    /* Skip the ifunc index */
                    i++;
                    for (int j = 0; j < nral; j++)
                    {
                        int aj = ril.il[i + j];
                        if (aj != a)
                        {
                            check_link(link, cg_gl, cg_offset + aj);
                        }
                    }
                    i += nral_rt(ftype);
                }

                if (mtop.bIntermolecularInteractions)
                {
                    int i = ril_intermol.index[cg_gl];
                    while (i < ril_intermol.index[cg_gl + 1])
                    {
                        int ftype = ril_intermol.il[i++];
                        int nral  = NRAL(ftype);
                        /* Skip the ifunc index */
                        i++;
                        for (int j = 0; j < nral; j++)
                        {
                            /* Here we assume we have no charge groups;
                             * this has been checked above.
                             */
                            int aj = ril_intermol.il[i + j];
                            check_link(link, cg_gl, aj);
                        }
                        i += nral_rt(ftype);
                    }
                }
                if (link->index[cg_gl + 1] - link->index[cg_gl] > 0)
                {
                    SET_CGINFO_BOND_INTER(cgi_mb->cginfo[a]);
                    ncgi++;
                }
            }

            cg_offset += molt.atoms.nr;
        }
        int nlink_mol = link->index[cg_offset] - link->index[cg_offset - molt.atoms.nr];

        if (debug)
        {
            fprintf(debug,
                    "molecule type '%s' %d atoms has %d atom links through bonded interac.\n",
                    *molt.name,
                    molt.atoms.nr,
                    nlink_mol);
        }

        if (molb.nmol > mol)
        {
            /* Copy the data for the rest of the molecules in this block */
            link->nalloc_a += (molb.nmol - mol) * nlink_mol;
            srenew(link->a, link->nalloc_a);
            for (; mol < molb.nmol; mol++)
            {
                for (int a = 0; a < molt.atoms.nr; a++)
                {
                    int cg_gl              = cg_offset + a;
                    link->index[cg_gl + 1] = link->index[cg_gl + 1 - molt.atoms.nr] + nlink_mol;
                    for (int j = link->index[cg_gl]; j < link->index[cg_gl + 1]; j++)
                    {
                        link->a[j] = link->a[j - nlink_mol] + molt.atoms.nr;
                    }
                    if (link->index[cg_gl + 1] - link->index[cg_gl] > 0
                        && cg_gl - cgi_mb->cg_start < cgi_mb->cg_mod)
                    {
                        SET_CGINFO_BOND_INTER(cgi_mb->cginfo[cg_gl - cgi_mb->cg_start]);
                        ncgi++;
                    }
                }
                cg_offset += molt.atoms.nr;
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "Of the %d atoms %d are linked via bonded interactions\n", mtop.natoms, ncgi);
    }

    return link;
}

typedef struct
{
    real r2;
    int  ftype;
    int  a1;
    int  a2;
} bonded_distance_t;

/*! \brief Compare distance^2 \p r2 against the distance in \p bd and if larger store it along with \p ftype and atom indices \p a1 and \p a2 */
static void update_max_bonded_distance(real r2, int ftype, int a1, int a2, bonded_distance_t* bd)
{
    if (r2 > bd->r2)
    {
        bd->r2    = r2;
        bd->ftype = ftype;
        bd->a1    = a1;
        bd->a2    = a2;
    }
}

/*! \brief Set the distance, function type and atom indices for the longest distance between atoms of molecule type \p molt for two-body and multi-body bonded interactions */
static void bonded_cg_distance_mol(const gmx_moltype_t*   molt,
                                   const DDBondedChecking ddBondedChecking,
                                   gmx_bool               bExcl,
                                   ArrayRef<const RVec>   x,
                                   bonded_distance_t*     bd_2b,
                                   bonded_distance_t*     bd_mb)
{
    const ReverseTopOptions rtOptions(ddBondedChecking);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rtOptions))
        {
            const auto& il   = molt->ilist[ftype];
            int         nral = NRAL(ftype);
            if (nral > 1)
            {
                for (int i = 0; i < il.size(); i += 1 + nral)
                {
                    for (int ai = 0; ai < nral; ai++)
                    {
                        int atomI = il.iatoms[i + 1 + ai];
                        for (int aj = ai + 1; aj < nral; aj++)
                        {
                            int atomJ = il.iatoms[i + 1 + aj];
                            if (atomI != atomJ)
                            {
                                real rij2 = distance2(x[atomI], x[atomJ]);

                                update_max_bonded_distance(
                                        rij2, ftype, atomI, atomJ, (nral == 2) ? bd_2b : bd_mb);
                            }
                        }
                    }
                }
            }
        }
    }
    if (bExcl)
    {
        const auto& excls = molt->excls;
        for (gmx::index ai = 0; ai < excls.ssize(); ai++)
        {
            for (const int aj : excls[ai])
            {
                if (ai != aj)
                {
                    real rij2 = distance2(x[ai], x[aj]);

                    /* There is no function type for exclusions, use -1 */
                    update_max_bonded_distance(rij2, -1, ai, aj, bd_2b);
                }
            }
        }
    }
}

/*! \brief Set the distance, function type and atom indices for the longest atom distance involved in intermolecular interactions for two-body and multi-body bonded interactions */
static void bonded_distance_intermol(const InteractionLists& ilists_intermol,
                                     const DDBondedChecking  ddBondedChecking,
                                     ArrayRef<const RVec>    x,
                                     PbcType                 pbcType,
                                     const matrix            box,
                                     bonded_distance_t*      bd_2b,
                                     bonded_distance_t*      bd_mb)
{
    t_pbc pbc;

    set_pbc(&pbc, pbcType, box);

    const ReverseTopOptions rtOptions(ddBondedChecking);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rtOptions))
        {
            const auto& il   = ilists_intermol[ftype];
            int         nral = NRAL(ftype);

            /* No nral>1 check here, since intermol interactions always
             * have nral>=2 (and the code is also correct for nral=1).
             */
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                for (int ai = 0; ai < nral; ai++)
                {
                    int atom_i = il.iatoms[i + 1 + ai];

                    for (int aj = ai + 1; aj < nral; aj++)
                    {
                        rvec dx;

                        int atom_j = il.iatoms[i + 1 + aj];

                        pbc_dx(&pbc, x[atom_i], x[atom_j], dx);

                        const real rij2 = norm2(dx);

                        update_max_bonded_distance(rij2, ftype, atom_i, atom_j, (nral == 2) ? bd_2b : bd_mb);
                    }
                }
            }
        }
    }
}

//! Returns whether \p molt has at least one virtual site
static bool moltypeHasVsite(const gmx_moltype_t& molt)
{
    bool hasVsite = false;
    for (int i = 0; i < F_NRE; i++)
    {
        if ((interaction_function[i].flags & IF_VSITE) && !molt.ilist[i].empty())
        {
            hasVsite = true;
        }
    }

    return hasVsite;
}

//! Returns coordinates not broken over PBC for a molecule
static void getWholeMoleculeCoordinates(const gmx_moltype_t*  molt,
                                        const gmx_ffparams_t* ffparams,
                                        PbcType               pbcType,
                                        t_graph*              graph,
                                        const matrix          box,
                                        ArrayRef<const RVec>  x,
                                        ArrayRef<RVec>        xs)
{
    if (pbcType != PbcType::No)
    {
        mk_mshift(nullptr, graph, pbcType, box, as_rvec_array(x.data()));

        shift_x(graph, box, as_rvec_array(x.data()), as_rvec_array(xs.data()));
        /* By doing an extra mk_mshift the molecules that are broken
         * because they were e.g. imported from another software
         * will be made whole again. Such are the healing powers
         * of GROMACS.
         */
        mk_mshift(nullptr, graph, pbcType, box, as_rvec_array(xs.data()));
    }
    else
    {
        /* We copy the coordinates so the original coordinates remain
         * unchanged, just to be 100% sure that we do not affect
         * binary reproducibility of simulations.
         */
        for (int i = 0; i < molt->atoms.nr; i++)
        {
            copy_rvec(x[i], xs[i]);
        }
    }

    if (moltypeHasVsite(*molt))
    {
        gmx::constructVirtualSites(xs, ffparams->iparams, molt->ilist);
    }
}

void dd_bonded_cg_distance(const gmx::MDLogger&   mdlog,
                           const gmx_mtop_t&      mtop,
                           const t_inputrec&      inputrec,
                           ArrayRef<const RVec>   x,
                           const matrix           box,
                           const DDBondedChecking ddBondedChecking,
                           real*                  r_2b,
                           real*                  r_mb)
{
    bonded_distance_t bd_2b = { 0, -1, -1, -1 };
    bonded_distance_t bd_mb = { 0, -1, -1, -1 };

    bool bExclRequired = inputrecExclForces(&inputrec);

    *r_2b         = 0;
    *r_mb         = 0;
    int at_offset = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];
        if (molt.atoms.nr == 1 || molb.nmol == 0)
        {
            at_offset += molb.nmol * molt.atoms.nr;
        }
        else
        {
            t_graph graph;
            if (inputrec.pbcType != PbcType::No)
            {
                graph = mk_graph_moltype(molt);
            }

            std::vector<RVec> xs(molt.atoms.nr);
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                getWholeMoleculeCoordinates(&molt,
                                            &mtop.ffparams,
                                            inputrec.pbcType,
                                            &graph,
                                            box,
                                            x.subArray(at_offset, molt.atoms.nr),
                                            xs);

                bonded_distance_t bd_mol_2b = { 0, -1, -1, -1 };
                bonded_distance_t bd_mol_mb = { 0, -1, -1, -1 };

                bonded_cg_distance_mol(&molt, ddBondedChecking, bExclRequired, xs, &bd_mol_2b, &bd_mol_mb);

                /* Process the mol data adding the atom index offset */
                update_max_bonded_distance(bd_mol_2b.r2,
                                           bd_mol_2b.ftype,
                                           at_offset + bd_mol_2b.a1,
                                           at_offset + bd_mol_2b.a2,
                                           &bd_2b);
                update_max_bonded_distance(bd_mol_mb.r2,
                                           bd_mol_mb.ftype,
                                           at_offset + bd_mol_mb.a1,
                                           at_offset + bd_mol_mb.a2,
                                           &bd_mb);

                at_offset += molt.atoms.nr;
            }
        }
    }

    if (mtop.bIntermolecularInteractions)
    {
        GMX_RELEASE_ASSERT(mtop.intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        bonded_distance_intermol(
                *mtop.intermolecular_ilist, ddBondedChecking, x, inputrec.pbcType, box, &bd_2b, &bd_mb);
    }

    *r_2b = sqrt(bd_2b.r2);
    *r_mb = sqrt(bd_mb.r2);

    if (*r_2b > 0 || *r_mb > 0)
    {
        GMX_LOG(mdlog.info).appendText("Initial maximum distances in bonded interactions:");
        if (*r_2b > 0)
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "    two-body bonded interactions: %5.3f nm, %s, atoms %d %d",
                            *r_2b,
                            (bd_2b.ftype >= 0) ? interaction_function[bd_2b.ftype].longname : "Exclusion",
                            bd_2b.a1 + 1,
                            bd_2b.a2 + 1);
        }
        if (*r_mb > 0)
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "  multi-body bonded interactions: %5.3f nm, %s, atoms %d %d",
                            *r_mb,
                            interaction_function[bd_mb.ftype].longname,
                            bd_mb.a1 + 1,
                            bd_mb.a2 + 1);
        }
    }
}
