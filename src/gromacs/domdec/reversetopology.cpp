/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 *
 * \brief This file defines functions used in making the
 * reverse topology.
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_domdec
 */

#include "gmxpre.h"

#include "gromacs/domdec/reversetopology.h"

#include <cstdio>

#include <algorithm>
#include <memory>
#include <vector>

#include "gromacs/domdec/domdec_constraints.h"
#include "gromacs/domdec/domdec_internal.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/domdec/domdec_vsite.h"
#include "gromacs/domdec/options.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

using gmx::ArrayRef;
using gmx::DDBondedChecking;
using gmx::ListOfLists;
using gmx::RVec;

/*! \brief Struct for the reverse topology: links bonded interactions to atoms */
struct gmx_reverse_top_t::Impl
{
    //! Constructs a reverse topology from \p mtop
    Impl(const gmx_mtop_t& mtop, bool useFreeEnergy, const ReverseTopOptions& reverseTopOptions);

    //! @cond Doxygen_Suppress
    //! Options for the setup of this reverse topology
    const ReverseTopOptions options;
    //! Are there interaction of type F_POSRES and/or F_FBPOSRES
    bool hasPositionRestraints;
    //! \brief Are there bondeds/exclusions between atoms?
    bool bInterAtomicInteractions = false;
    //! \brief Reverse ilist for all moltypes
    std::vector<reverse_ilist_t> ril_mt;
    //! \brief The size of ril_mt[?].index summed over all entries
    int ril_mt_tot_size = 0;
    //! \brief Whether listed-force interaction lists should be sorted for free energy
    bool doListedForcesSorting;
    //! \brief molblock to global atom index for quick lookup of molblocks on atom index
    std::vector<MolblockIndices> mbi;

    //! \brief Do we have intermolecular interactions?
    bool bIntermolecularInteractions = false;
    //! \brief Intermolecular reverse ilist
    reverse_ilist_t ril_intermol;

    /* Work data structures for multi-threading */
    //! \brief Thread work array for local topology generation
    std::vector<thread_work_t> th_work;
    //! @endcond
};


int nral_rt(int ftype)
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

bool dd_check_ftype(const int ftype, const ReverseTopOptions& rtOptions)
{
    return ((((interaction_function[ftype].flags & IF_BOND) != 0U)
             && ((interaction_function[ftype].flags & IF_VSITE) == 0U)
             && ((rtOptions.ddBondedChecking_ == DDBondedChecking::All)
                 || ((interaction_function[ftype].flags & IF_LIMZERO) == 0U)))
            || (rtOptions.includeConstraints_ && (ftype == F_CONSTR || ftype == F_CONSTRNC))
            || (rtOptions.includeSettles_ && ftype == F_SETTLE));
}

MolecularTopologyAtomIndices globalAtomIndexToMoltypeIndices(const gmx::ArrayRef<const MolblockIndices> molblockIndices,
                                                             const int globalAtomIndex)
{
    // Find the molblock the atom belongs to using bisection
    int start = 0;
    int end   = molblockIndices.size(); /* exclusive */
    int mid   = 0;

    while (true)
    {
        mid = (start + end) / 2;
        if (globalAtomIndex >= molblockIndices[mid].a_end)
        {
            start = mid + 1;
        }
        else if (globalAtomIndex < molblockIndices[mid].a_start)
        {
            end = mid;
        }
        else
        {
            break;
        }
    }

    const MolblockIndices& mbi = molblockIndices[mid];

    MolecularTopologyAtomIndices mtai;

    mtai.blockIndex    = mid;
    mtai.moleculeType  = mbi.type;
    mtai.moleculeIndex = (globalAtomIndex - mbi.a_start) / mbi.natoms_mol;
    mtai.atomIndex     = (globalAtomIndex - mbi.a_start) - mtai.moleculeIndex * mbi.natoms_mol;

    return mtai;
}

/*! \brief Returns the maximum number of exclusions per atom */
static int getMaxNumExclusionsPerAtom(const ListOfLists<int>& excls)
{
    int maxNumExcls = 0;
    for (gmx::Index at = 0; at < excls.ssize(); at++)
    {
        const auto list     = excls[at];
        const int  numExcls = list.ssize();

        GMX_RELEASE_ASSERT(numExcls != 1 || list[0] == at,
                           "With 1 exclusion we expect a self-exclusion");

        maxNumExcls = std::max(maxNumExcls, numExcls);
    }

    return maxNumExcls;
}

/*! \brief Run the reverse ilist generation and store it in r_il when \p bAssign = TRUE */
static void low_make_reverse_ilist(const InteractionLists&  il_mt,
                                   const t_atom*            atom,
                                   int*                     count,
                                   const ReverseTopOptions& rtOptions,
                                   gmx::ArrayRef<const int> r_index,
                                   gmx::ArrayRef<int>       r_il,
                                   const AtomLinkRule       atomLinkRule,
                                   const bool               assignReverseIlist)
{
    const bool includeConstraints = rtOptions.includeConstraints_;
    const bool includeSettles     = rtOptions.includeSettles_;

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
                    count[a] += 2 + nral_rt(ftype);
                }
            }
        }
    }
}

void make_reverse_ilist(const InteractionLists&  ilist,
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
    low_make_reverse_ilist(
            ilist, atoms->atom, count.data(), rtOptions, ril_mt->index, ril_mt->il, atomLinkRule, TRUE);

    ril_mt->numAtomsInMolecule = atoms->nr;
}

gmx_reverse_top_t::gmx_reverse_top_t(const gmx_mtop_t&        mtop,
                                     bool                     useFreeEnergy,
                                     const ReverseTopOptions& reverseTopOptions) :
    impl_(std::make_unique<Impl>(mtop, useFreeEnergy, reverseTopOptions))
{
}

gmx_reverse_top_t::~gmx_reverse_top_t() {}

const ReverseTopOptions& gmx_reverse_top_t::options() const
{
    return impl_->options;
}

const reverse_ilist_t& gmx_reverse_top_t::interactionListForMoleculeType(int moleculeType) const
{
    return impl_->ril_mt[moleculeType];
}

ArrayRef<const MolblockIndices> gmx_reverse_top_t::molblockIndices() const
{
    return impl_->mbi;
}

bool gmx_reverse_top_t::hasIntermolecularInteractions() const
{
    return impl_->bIntermolecularInteractions;
}

const reverse_ilist_t& gmx_reverse_top_t::interactionListForIntermolecularInteractions() const
{
    return impl_->ril_intermol;
}

bool gmx_reverse_top_t::hasInterAtomicInteractions() const
{
    return impl_->bInterAtomicInteractions;
}

bool gmx_reverse_top_t::hasPositionRestraints() const
{
    return impl_->hasPositionRestraints;
}

ArrayRef<thread_work_t> gmx_reverse_top_t::threadWorkObjects() const
{
    return impl_->th_work;
}

bool gmx_reverse_top_t::doListedForcesSorting() const
{
    return impl_->doListedForcesSorting;
}

/*! \brief Generate the reverse topology */
gmx_reverse_top_t::Impl::Impl(const gmx_mtop_t&        mtop,
                              const bool               useFreeEnergy,
                              const ReverseTopOptions& reverseTopOptions) :
    options(reverseTopOptions),
    hasPositionRestraints(gmx_mtop_ftype_count(mtop, F_POSRES) + gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0),
    bInterAtomicInteractions(mtop.bIntermolecularInteractions)
{
    bInterAtomicInteractions = mtop.bIntermolecularInteractions;
    ril_mt.resize(mtop.moltype.size());
    ril_mt_tot_size = 0;
    for (size_t mt = 0; mt < mtop.moltype.size(); mt++)
    {
        const gmx_moltype_t& molt = mtop.moltype[mt];
        if (molt.atoms.nr > 1)
        {
            bInterAtomicInteractions = true;
        }

        /* Make the atom to interaction list for this molecule type */
        make_reverse_ilist(molt.ilist, &molt.atoms, options, AtomLinkRule::FirstAtom, &ril_mt[mt]);

        ril_mt_tot_size += ril_mt[mt].index[molt.atoms.nr];
    }
    if (debug)
    {
        fprintf(debug, "The total size of the atom to interaction index is %d integers\n", ril_mt_tot_size);
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

        make_reverse_ilist(
                *mtop.intermolecular_ilist, &atoms_global, options, AtomLinkRule::FirstAtom, &ril_intermol);
    }

    doListedForcesSorting = useFreeEnergy && gmx_mtop_bondeds_free_energy(&mtop);

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
