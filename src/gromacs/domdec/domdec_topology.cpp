/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2006 - 2014, The GROMACS development team.
 * Copyright (c) 2015,2016,2017,2018,2019,2020, by the GROMACS development team, led by
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
#include <string>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_network.h"
#include "gromacs/domdec/ga2la.h"
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

using gmx::ListOfLists;

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

    InteractionDefinitions    idef;       /**< Partial local topology */
    std::unique_ptr<VsitePbc> vsitePbc;   /**< vsite PBC structure */
    int                       nbonded;    /**< The number of bondeds in this struct */
    ListOfLists<int>          excl;       /**< List of exclusions */
    int                       excl_count; /**< The total exclusion count for \p excl */
};

/*! \brief Struct for the reverse topology: links bonded interactions to atomsx */
struct gmx_reverse_top_t
{
    //! @cond Doxygen_Suppress
    //! \brief Are there constraints in this revserse top?
    bool bConstr = false;
    //! \brief Are there settles in this revserse top?
    bool bSettle = false;
    //! \brief All bonded interactions have to be assigned?
    bool bBCheck = false;
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

    /* Work data structures for multi-threading */
    //! \brief Thread work array for local topology generation
    std::vector<thread_work_t> th_work;
    //! @endcond
};

/*! \brief Returns the number of atom entries for il in gmx_reverse_top_t */
static int nral_rt(int ftype)
{
    int nral;

    nral = NRAL(ftype);
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
static gmx_bool dd_check_ftype(int ftype, gmx_bool bBCheck, gmx_bool bConstr, gmx_bool bSettle)
{
    return ((((interaction_function[ftype].flags & IF_BOND) != 0U)
             && ((interaction_function[ftype].flags & IF_VSITE) == 0U)
             && (bBCheck || ((interaction_function[ftype].flags & IF_LIMZERO) == 0U)))
            || (bConstr && (ftype == F_CONSTR || ftype == F_CONSTRNC)) || (bSettle && ftype == F_SETTLE));
}

/*! \brief Help print error output when interactions are missing
 *
 * \note This function needs to be called on all ranks (contains a global summation)
 */
static std::string print_missing_interactions_mb(t_commrec*                    cr,
                                                 const gmx_reverse_top_t*      rt,
                                                 const char*                   moltypename,
                                                 const reverse_ilist_t*        ril,
                                                 int                           a_start,
                                                 int                           a_end,
                                                 int                           nat_mol,
                                                 int                           nmol,
                                                 const InteractionDefinitions* idef)
{
    int* assigned;
    int  nril_mol = ril->index[nat_mol];
    snew(assigned, nmol * nril_mol);
    gmx::StringOutputStream stream;
    gmx::TextWriter         log(&stream);

    gmx::ArrayRef<const int> gatindex = cr->dd->globalAtomIndices;
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, rt->bBCheck, rt->bConstr, rt->bSettle))
        {
            int                    nral = NRAL(ftype);
            const InteractionList* il   = &idef->il[ftype];
            const int*             ia   = il->iatoms.data();
            for (int i = 0; i < il->size(); i += 1 + nral)
            {
                int a0 = gatindex[ia[1]];
                /* Check if this interaction is in
                 * the currently checked molblock.
                 */
                if (a0 >= a_start && a0 < a_end)
                {
                    int  mol    = (a0 - a_start) / nat_mol;
                    int  a0_mol = (a0 - a_start) - mol * nat_mol;
                    int  j_mol  = ril->index[a0_mol];
                    bool found  = false;
                    while (j_mol < ril->index[a0_mol + 1] && !found)
                    {
                        int j       = mol * nril_mol + j_mol;
                        int ftype_j = ril->il[j_mol];
                        /* Here we need to check if this interaction has
                         * not already been assigned, since we could have
                         * multiply defined interactions.
                         */
                        if (ftype == ftype_j && ia[0] == ril->il[j_mol + 1] && assigned[j] == 0)
                        {
                            /* Check the atoms */
                            found = true;
                            for (int a = 0; a < nral; a++)
                            {
                                if (gatindex[ia[1 + a]] != a_start + mol * nat_mol + ril->il[j_mol + 2 + a])
                                {
                                    found = false;
                                }
                            }
                            if (found)
                            {
                                assigned[j] = 1;
                            }
                        }
                        j_mol += 2 + nral_rt(ftype_j);
                    }
                    if (!found)
                    {
                        gmx_incons("Some interactions seem to be assigned multiple times");
                    }
                }
                ia += 1 + nral;
            }
        }
    }

    gmx_sumi(nmol * nril_mol, assigned, cr);

    int nprint = 10;
    int i      = 0;
    for (int mol = 0; mol < nmol; mol++)
    {
        int j_mol = 0;
        while (j_mol < nril_mol)
        {
            int ftype = ril->il[j_mol];
            int nral  = NRAL(ftype);
            int j     = mol * nril_mol + j_mol;
            if (assigned[j] == 0 && !(interaction_function[ftype].flags & IF_VSITE))
            {
                if (DDMASTER(cr->dd))
                {
                    if (i == 0)
                    {
                        log.writeLineFormatted("Molecule type '%s'", moltypename);
                        log.writeLineFormatted(
                                "the first %d missing interactions, except for exclusions:", nprint);
                    }
                    log.writeStringFormatted("%20s atoms", interaction_function[ftype].longname);
                    int a;
                    for (a = 0; a < nral; a++)
                    {
                        log.writeStringFormatted("%5d", ril->il[j_mol + 2 + a] + 1);
                    }
                    while (a < 4)
                    {
                        log.writeString("     ");
                        a++;
                    }
                    log.writeString(" global");
                    for (a = 0; a < nral; a++)
                    {
                        log.writeStringFormatted(
                                "%6d", a_start + mol * nat_mol + ril->il[j_mol + 2 + a] + 1);
                    }
                    log.ensureLineBreak();
                }
                i++;
                if (i >= nprint)
                {
                    break;
                }
            }
            j_mol += 2 + nral_rt(ftype);
        }
    }

    sfree(assigned);
    return stream.toString();
}

/*! \brief Help print error output when interactions are missing */
static void print_missing_interactions_atoms(const gmx::MDLogger&          mdlog,
                                             t_commrec*                    cr,
                                             const gmx_mtop_t*             mtop,
                                             const InteractionDefinitions* idef)
{
    const gmx_reverse_top_t* rt = cr->dd->reverse_top;

    /* Print the atoms in the missing interactions per molblock */
    int a_end = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& moltype = mtop->moltype[molb.type];
        int                  a_start = a_end;
        a_end                        = a_start + molb.nmol * moltype.atoms.nr;

        auto warning = print_missing_interactions_mb(cr, rt, *(moltype.name), &rt->ril_mt[molb.type],
                                                     a_start, a_end, moltype.atoms.nr, molb.nmol, idef);

        GMX_LOG(mdlog.warning).appendText(warning);
    }
}

void dd_print_missing_interactions(const gmx::MDLogger&           mdlog,
                                   t_commrec*                     cr,
                                   int                            local_count,
                                   const gmx_mtop_t*              top_global,
                                   const gmx_localtop_t*          top_local,
                                   gmx::ArrayRef<const gmx::RVec> x,
                                   const matrix                   box)
{
    int           ndiff_tot, cl[F_NRE], n, ndiff, rest_global, rest_local;
    int           ftype, nral;
    gmx_domdec_t* dd;

    dd = cr->dd;

    GMX_LOG(mdlog.warning)
            .appendText(
                    "Not all bonded interactions have been properly assigned to the domain "
                    "decomposition cells");

    ndiff_tot = local_count - dd->nbonded_global;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        nral      = NRAL(ftype);
        cl[ftype] = top_local->idef.il[ftype].size() / (1 + nral);
    }

    gmx_sumi(F_NRE, cl, cr);

    if (DDMASTER(dd))
    {
        GMX_LOG(mdlog.warning).appendText("A list of missing interactions:");
        rest_global = dd->nbonded_global;
        rest_local  = local_count;
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            /* In the reverse and local top all constraints are merged
             * into F_CONSTR. So in the if statement we skip F_CONSTRNC
             * and add these constraints when doing F_CONSTR.
             */
            if (((interaction_function[ftype].flags & IF_BOND)
                 && (dd->reverse_top->bBCheck || !(interaction_function[ftype].flags & IF_LIMZERO)))
                || (dd->reverse_top->bConstr && ftype == F_CONSTR)
                || (dd->reverse_top->bSettle && ftype == F_SETTLE))
            {
                n = gmx_mtop_ftype_count(top_global, ftype);
                if (ftype == F_CONSTR)
                {
                    n += gmx_mtop_ftype_count(top_global, F_CONSTRNC);
                }
                ndiff = cl[ftype] - n;
                if (ndiff != 0)
                {
                    GMX_LOG(mdlog.warning)
                            .appendTextFormatted("%20s of %6d missing %6d",
                                                 interaction_function[ftype].longname, n, -ndiff);
                }
                rest_global -= n;
                rest_local -= cl[ftype];
            }
        }

        ndiff = rest_local - rest_global;
        if (ndiff != 0)
        {
            GMX_LOG(mdlog.warning).appendTextFormatted("%20s of %6d missing %6d", "exclusions", rest_global, -ndiff);
        }
    }

    print_missing_interactions_atoms(mdlog, cr, top_global, &top_local->idef);
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
                -ndiff_tot, cr->dd->nbonded_global, dd_cutoff_multibody(dd), dd_cutoff_twobody(dd));
    }
    gmx_fatal_collective(FARGS, cr->mpi_comm_mygroup, MASTER(cr), "%s", errorMessage.c_str());
}

/*! \brief Return global topology molecule information for global atom index \p i_gl */
static void global_atomnr_to_moltype_ind(const gmx_reverse_top_t* rt, int i_gl, int* mb, int* mt, int* mol, int* i_mol)
{
    const MolblockIndices* mbi   = rt->mbi.data();
    int                    start = 0;
    int                    end   = rt->mbi.size(); /* exclusive */
    int                    mid;

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

/*! \brief Run the reverse ilist generation and store it in r_il when \p bAssign = TRUE */
static int low_make_reverse_ilist(const InteractionLists&  il_mt,
                                  const t_atom*            atom,
                                  int*                     count,
                                  gmx_bool                 bConstr,
                                  gmx_bool                 bSettle,
                                  gmx_bool                 bBCheck,
                                  gmx::ArrayRef<const int> r_index,
                                  gmx::ArrayRef<int>       r_il,
                                  gmx_bool                 bLinkToAllAtoms,
                                  gmx_bool                 bAssign)
{
    int ftype, j, nlink, link;
    int a;
    int nint;

    nint = 0;
    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if ((interaction_function[ftype].flags & (IF_BOND | IF_VSITE))
            || (bConstr && (ftype == F_CONSTR || ftype == F_CONSTRNC)) || (bSettle && ftype == F_SETTLE))
        {
            const bool  bVSite = ((interaction_function[ftype].flags & IF_VSITE) != 0U);
            const int   nral   = NRAL(ftype);
            const auto& il     = il_mt[ftype];
            for (int i = 0; i < il.size(); i += 1 + nral)
            {
                const int* ia = il.iatoms.data() + i;
                if (bLinkToAllAtoms)
                {
                    if (bVSite)
                    {
                        /* We don't need the virtual sites for the cg-links */
                        nlink = 0;
                    }
                    else
                    {
                        nlink = nral;
                    }
                }
                else
                {
                    /* Couple to the first atom in the interaction */
                    nlink = 1;
                }
                for (link = 0; link < nlink; link++)
                {
                    a = ia[1 + link];
                    if (bAssign)
                    {
                        GMX_ASSERT(!r_il.empty(), "with bAssign not allowed to be empty");
                        GMX_ASSERT(!r_index.empty(), "with bAssign not allowed to be empty");
                        r_il[r_index[a] + count[a]]     = (ftype == F_CONSTRNC ? F_CONSTR : ftype);
                        r_il[r_index[a] + count[a] + 1] = ia[0];
                        for (j = 1; j < 1 + nral; j++)
                        {
                            /* Store the molecular atom number */
                            r_il[r_index[a] + count[a] + 1 + j] = ia[j];
                        }
                    }
                    if (interaction_function[ftype].flags & IF_VSITE)
                    {
                        if (bAssign)
                        {
                            /* Add an entry to iatoms for storing
                             * which of the constructing atoms are
                             * vsites again.
                             */
                            r_il[r_index[a] + count[a] + 2 + nral] = 0;
                            for (j = 2; j < 1 + nral; j++)
                            {
                                if (atom[ia[j]].ptype == eptVSite)
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
                        if (bBCheck || !(interaction_function[ftype].flags & IF_LIMZERO))
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
static int make_reverse_ilist(const InteractionLists& ilist,
                              const t_atoms*          atoms,
                              gmx_bool                bConstr,
                              gmx_bool                bSettle,
                              gmx_bool                bBCheck,
                              gmx_bool                bLinkToAllAtoms,
                              reverse_ilist_t*        ril_mt)
{
    int nat_mt, *count, i, nint_mt;

    /* Count the interactions */
    nat_mt = atoms->nr;
    snew(count, nat_mt);
    low_make_reverse_ilist(ilist, atoms->atom, count, bConstr, bSettle, bBCheck, {}, {},
                           bLinkToAllAtoms, FALSE);

    ril_mt->index.push_back(0);
    for (i = 0; i < nat_mt; i++)
    {
        ril_mt->index.push_back(ril_mt->index[i] + count[i]);
        count[i] = 0;
    }
    ril_mt->il.resize(ril_mt->index[nat_mt]);

    /* Store the interactions */
    nint_mt = low_make_reverse_ilist(ilist, atoms->atom, count, bConstr, bSettle, bBCheck,
                                     ril_mt->index, ril_mt->il, bLinkToAllAtoms, TRUE);

    sfree(count);

    ril_mt->numAtomsInMolecule = atoms->nr;

    return nint_mt;
}

/*! \brief Generate the reverse topology */
static gmx_reverse_top_t make_reverse_top(const gmx_mtop_t* mtop,
                                          gmx_bool          bFE,
                                          gmx_bool          bConstr,
                                          gmx_bool          bSettle,
                                          gmx_bool          bBCheck,
                                          int*              nint)
{
    gmx_reverse_top_t rt;

    /* Should we include constraints (for SHAKE) in rt? */
    rt.bConstr = bConstr;
    rt.bSettle = bSettle;
    rt.bBCheck = bBCheck;

    rt.bInterAtomicInteractions = mtop->bIntermolecularInteractions;
    rt.ril_mt.resize(mtop->moltype.size());
    rt.ril_mt_tot_size = 0;
    std::vector<int> nint_mt;
    for (size_t mt = 0; mt < mtop->moltype.size(); mt++)
    {
        const gmx_moltype_t& molt = mtop->moltype[mt];
        if (molt.atoms.nr > 1)
        {
            rt.bInterAtomicInteractions = true;
        }

        /* Make the atom to interaction list for this molecule type */
        int numberOfInteractions = make_reverse_ilist(
                molt.ilist, &molt.atoms, rt.bConstr, rt.bSettle, rt.bBCheck, FALSE, &rt.ril_mt[mt]);
        nint_mt.push_back(numberOfInteractions);

        rt.ril_mt_tot_size += rt.ril_mt[mt].index[molt.atoms.nr];
    }
    if (debug)
    {
        fprintf(debug, "The total size of the atom to interaction index is %d integers\n",
                rt.ril_mt_tot_size);
    }

    *nint = 0;
    for (const gmx_molblock_t& molblock : mtop->molblock)
    {
        *nint += molblock.nmol * nint_mt[molblock.type];
    }

    /* Make an intermolecular reverse top, if necessary */
    rt.bIntermolecularInteractions = mtop->bIntermolecularInteractions;
    if (rt.bIntermolecularInteractions)
    {
        t_atoms atoms_global;

        atoms_global.nr   = mtop->natoms;
        atoms_global.atom = nullptr; /* Only used with virtual sites */

        GMX_RELEASE_ASSERT(mtop->intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        *nint += make_reverse_ilist(*mtop->intermolecular_ilist, &atoms_global, rt.bConstr,
                                    rt.bSettle, rt.bBCheck, FALSE, &rt.ril_intermol);
    }

    if (bFE && gmx_mtop_bondeds_free_energy(mtop))
    {
        rt.ilsort = ilsortFE_UNSORTED;
    }
    else
    {
        rt.ilsort = ilsortNO_FE;
    }

    /* Make a molblock index for fast searching */
    int i = 0;
    for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
    {
        const gmx_molblock_t& molb           = mtop->molblock[mb];
        const int             numAtomsPerMol = mtop->moltype[molb.type].atoms.nr;
        MolblockIndices       mbi;
        mbi.a_start = i;
        i += molb.nmol * numAtomsPerMol;
        mbi.a_end      = i;
        mbi.natoms_mol = numAtomsPerMol;
        mbi.type       = molb.type;
        rt.mbi.push_back(mbi);
    }

    for (int th = 0; th < gmx_omp_nthreads_get(emntDomdec); th++)
    {
        rt.th_work.emplace_back(mtop->ffparams);
    }

    return rt;
}

void dd_make_reverse_top(FILE*              fplog,
                         gmx_domdec_t*      dd,
                         const gmx_mtop_t*  mtop,
                         const gmx_vsite_t* vsite,
                         const t_inputrec*  ir,
                         gmx_bool           bBCheck)
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

    dd->reverse_top = new gmx_reverse_top_t;
    *dd->reverse_top =
            make_reverse_top(mtop, ir->efep != efepNO, !dd->comm->systemInfo.haveSplitConstraints,
                             !dd->comm->systemInfo.haveSplitSettles, bBCheck, &dd->nbonded_global);

    dd->haveExclusions = false;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const int maxNumExclusionsPerAtom = getMaxNumExclusionsPerAtom(mtop->moltype[molb.type].excls);
        // We checked above that max 1 exclusion means only self exclusions
        if (maxNumExclusionsPerAtom > 1)
        {
            dd->haveExclusions = true;
        }
    }

    if (vsite && vsite->numInterUpdategroupVsites > 0)
    {
        if (fplog)
        {
            fprintf(fplog,
                    "There are %d inter update-group virtual sites,\n"
                    "will an extra communication step for selected coordinates and forces\n",
                    vsite->numInterUpdategroupVsites);
        }
        init_domdec_vsites(dd, vsite->numInterUpdategroupVsites);
    }

    if (dd->comm->systemInfo.haveSplitConstraints || dd->comm->systemInfo.haveSplitSettles)
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
                       const gmx_molblock_t*   molb,
                       t_iatom*                iatoms,
                       const t_iparams*        ip_in,
                       InteractionDefinitions* idef)
{
    /* This position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    const int n = idef->il[F_POSRES].size() / 2;

    /* Get the position restraint coordinates from the molblock */
    const int a_molb = mol * numAtomsInMolecule + a_mol;
    GMX_ASSERT(a_molb < ssize(molb->posres_xA),
               "We need a sufficient number of position restraint coordinates");
    /* Copy the force constants */
    t_iparams ip        = ip_in[iatoms[0]];
    ip.posres.pos0A[XX] = molb->posres_xA[a_molb][XX];
    ip.posres.pos0A[YY] = molb->posres_xA[a_molb][YY];
    ip.posres.pos0A[ZZ] = molb->posres_xA[a_molb][ZZ];
    if (!molb->posres_xB.empty())
    {
        ip.posres.pos0B[XX] = molb->posres_xB[a_molb][XX];
        ip.posres.pos0B[YY] = molb->posres_xB[a_molb][YY];
        ip.posres.pos0B[ZZ] = molb->posres_xB[a_molb][ZZ];
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
                         const gmx_molblock_t*   molb,
                         t_iatom*                iatoms,
                         const t_iparams*        ip_in,
                         InteractionDefinitions* idef)
{
    /* This flat-bottom position restraint has not been added yet,
     * so it's index is the current number of position restraints.
     */
    const int n = idef->il[F_FBPOSRES].size() / 2;

    /* Get the position restraint coordinats from the molblock */
    const int a_molb = mol * numAtomsInMolecule + a_mol;
    GMX_ASSERT(a_molb < ssize(molb->posres_xA),
               "We need a sufficient number of position restraint coordinates");
    /* Copy the force constants */
    t_iparams ip = ip_in[iatoms[0]];
    /* Take reference positions from A position of normal posres */
    ip.fbposres.pos0[XX] = molb->posres_xA[a_molb][XX];
    ip.fbposres.pos0[YY] = molb->posres_xA[a_molb][YY];
    ip.fbposres.pos0[ZZ] = molb->posres_xA[a_molb][ZZ];

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
    int     k;
    t_iatom tiatoms[1 + MAXATOMLIST];
    int     j, ftype_r, nral_r;

    /* Add this interaction to the local topology */
    add_ifunc_for_vsites(tiatoms, ga2la, nral, bHomeA, a, a_gl, a_mol, iatoms, &idef->il[ftype]);

    if (iatoms[1 + nral])
    {
        /* Check for recursion */
        for (k = 2; k < 1 + nral; k++)
        {
            if ((iatoms[1 + nral] & (2 << k)) && (tiatoms[k] < 0))
            {
                /* This construction atoms is a vsite and not a home atom */
                if (gmx_debug_at)
                {
                    fprintf(debug, "Constructing atom %d of vsite atom %d is a vsite and non-home\n",
                            iatoms[k] + 1, a_mol + 1);
                }
                /* Find the vsite construction */

                /* Check all interactions assigned to this atom */
                j = index[iatoms[k]];
                while (j < index[iatoms[k] + 1])
                {
                    ftype_r = rtil[j++];
                    nral_r  = NRAL(ftype_r);
                    if (interaction_function[ftype_r].flags & IF_VSITE)
                    {
                        /* Add this vsite (recursion) */
                        add_vsite(ga2la, index, rtil, ftype_r, nral_r, FALSE, -1,
                                  a_gl + iatoms[k] - iatoms[1], iatoms[k], rtil.data() + j, idef);
                    }
                    j += 1 + nral_rt(ftype_r);
                }
            }
        }
    }
}

/*! \brief Returns the squared distance between atoms \p i and \p j */
static real dd_dist2(t_pbc* pbc_null, const rvec* x, const int i, int j)
{
    rvec dx;

    if (pbc_null)
    {
        pbc_dx_aiuc(pbc_null, x[i], x[j], dx);
    }
    else
    {
        rvec_sub(x[i], x[j], dx);
    }

    return norm2(dx);
}

/*! \brief Append t_idef structures 1 to nsrc in src to *dest */
static void combine_idef(InteractionDefinitions* dest, gmx::ArrayRef<const thread_work_t> src)
{
    int ftype;

    for (ftype = 0; ftype < F_NRE; ftype++)
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

/*! \brief Check and when available assign bonded interactions for local atom i
 */
static inline void check_assign_interactions_atom(int                       i,
                                                  int                       i_gl,
                                                  int                       mol,
                                                  int                       i_mol,
                                                  int                       numAtomsInMolecule,
                                                  gmx::ArrayRef<const int>  index,
                                                  gmx::ArrayRef<const int>  rtil,
                                                  gmx_bool                  bInterMolInteractions,
                                                  int                       ind_start,
                                                  int                       ind_end,
                                                  const gmx_domdec_t*       dd,
                                                  const gmx_domdec_zones_t* zones,
                                                  const gmx_molblock_t*     molb,
                                                  gmx_bool                  bRCheckMB,
                                                  const ivec                rcheck,
                                                  gmx_bool                  bRCheck2B,
                                                  real                      rc2,
                                                  t_pbc*                    pbc_null,
                                                  rvec*                     cg_cm,
                                                  const t_iparams*          ip_in,
                                                  InteractionDefinitions*   idef,
                                                  int                       iz,
                                                  gmx_bool                  bBCheck,
                                                  int*                      nbonded_local)
{
    gmx::ArrayRef<const DDPairInteractionRanges> iZones = zones->iZones;

    int j = ind_start;
    while (j < ind_end)
    {
        t_iatom tiatoms[1 + MAXATOMLIST];

        const int ftype  = rtil[j++];
        auto      iatoms = gmx::constArrayRefFromArray(rtil.data() + j, rtil.size() - j);
        const int nral   = NRAL(ftype);
        if (interaction_function[ftype].flags & IF_VSITE)
        {
            assert(!bInterMolInteractions);
            /* The vsite construction goes where the vsite itself is */
            if (iz == 0)
            {
                add_vsite(*dd->ga2la, index, rtil, ftype, nral, TRUE, i, i_gl, i_mol, iatoms.data(), idef);
            }
            j += 1 + nral + 2;
        }
        else
        {
            gmx_bool bUse;

            /* Copy the type */
            tiatoms[0] = iatoms[0];

            if (nral == 1)
            {
                assert(!bInterMolInteractions);
                /* Assign single-body interactions to the home zone */
                if (iz == 0)
                {
                    bUse       = TRUE;
                    tiatoms[1] = i;
                    if (ftype == F_POSRES)
                    {
                        add_posres(mol, i_mol, numAtomsInMolecule, molb, tiatoms, ip_in, idef);
                    }
                    else if (ftype == F_FBPOSRES)
                    {
                        add_fbposres(mol, i_mol, numAtomsInMolecule, molb, tiatoms, ip_in, idef);
                    }
                }
                else
                {
                    bUse = FALSE;
                }
            }
            else if (nral == 2)
            {
                /* This is a two-body interaction, we can assign
                 * analogous to the non-bonded assignments.
                 */
                int k_gl;

                if (!bInterMolInteractions)
                {
                    /* Get the global index using the offset in the molecule */
                    k_gl = i_gl + iatoms[2] - i_mol;
                }
                else
                {
                    k_gl = iatoms[2];
                }
                if (const auto* entry = dd->ga2la->find(k_gl))
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

                        tiatoms[1] = i;
                        tiatoms[2] = entry->la;
                        /* If necessary check the cgcm distance */
                        if (bRCheck2B && dd_dist2(pbc_null, cg_cm, tiatoms[1], tiatoms[2]) >= rc2)
                        {
                            bUse = FALSE;
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
                 * the local node if we have all the atoms involved
                 * (local or communicated) and the minimum zone shift
                 * in each dimension is zero, for dimensions
                 * with 2 DD cells an extra check may be necessary.
                 */
                ivec k_zero, k_plus;
                int  k;

                bUse = TRUE;
                clear_ivec(k_zero);
                clear_ivec(k_plus);
                for (k = 1; k <= nral && bUse; k++)
                {
                    int k_gl;
                    if (!bInterMolInteractions)
                    {
                        /* Get the global index using the offset in the molecule */
                        k_gl = i_gl + iatoms[k] - i_mol;
                    }
                    else
                    {
                        k_gl = iatoms[k];
                    }
                    const auto* entry = dd->ga2la->find(k_gl);
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
                        int d;

                        tiatoms[k] = entry->la;
                        for (d = 0; d < DIM; d++)
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
                if (bRCheckMB)
                {
                    int d;

                    for (d = 0; (d < DIM && bUse); d++)
                    {
                        /* Check if the cg_cm distance falls within
                         * the cut-off to avoid possible multiple
                         * assignments of bonded interactions.
                         */
                        if (rcheck[d] && k_plus[d]
                            && dd_dist2(pbc_null, cg_cm, tiatoms[k_zero[d]], tiatoms[k_plus[d]]) >= rc2)
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
                if (bBCheck || !(interaction_function[ftype].flags & IF_LIMZERO))
                {
                    (*nbonded_local)++;
                }
            }
            j += 1 + nral;
        }
    }
}

/*! \brief This function looks up and assigns bonded interactions for zone iz.
 *
 * With thread parallelizing each thread acts on a different atom range:
 * at_start to at_end.
 */
static int make_bondeds_zone(gmx_domdec_t*                      dd,
                             const gmx_domdec_zones_t*          zones,
                             const std::vector<gmx_molblock_t>& molb,
                             gmx_bool                           bRCheckMB,
                             ivec                               rcheck,
                             gmx_bool                           bRCheck2B,
                             real                               rc2,
                             t_pbc*                             pbc_null,
                             rvec*                              cg_cm,
                             const t_iparams*                   ip_in,
                             InteractionDefinitions*            idef,
                             int                                izone,
                             const gmx::Range<int>&             atomRange)
{
    int                mb, mt, mol, i_mol;
    gmx_bool           bBCheck;
    gmx_reverse_top_t* rt;
    int                nbonded_local;

    rt = dd->reverse_top;

    bBCheck = rt->bBCheck;

    nbonded_local = 0;

    for (int i : atomRange)
    {
        /* Get the global atom number */
        const int i_gl = dd->globalAtomIndices[i];
        global_atomnr_to_moltype_ind(rt, i_gl, &mb, &mt, &mol, &i_mol);
        /* Check all intramolecular interactions assigned to this atom */
        gmx::ArrayRef<const int>     index = rt->ril_mt[mt].index;
        gmx::ArrayRef<const t_iatom> rtil  = rt->ril_mt[mt].il;

        check_assign_interactions_atom(i, i_gl, mol, i_mol, rt->ril_mt[mt].numAtomsInMolecule,
                                       index, rtil, FALSE, index[i_mol], index[i_mol + 1], dd,
                                       zones, &molb[mb], bRCheckMB, rcheck, bRCheck2B, rc2,
                                       pbc_null, cg_cm, ip_in, idef, izone, bBCheck, &nbonded_local);


        if (rt->bIntermolecularInteractions)
        {
            /* Check all intermolecular interactions assigned to this atom */
            index = rt->ril_intermol.index;
            rtil  = rt->ril_intermol.il;

            check_assign_interactions_atom(i, i_gl, mol, i_mol, rt->ril_mt[mt].numAtomsInMolecule,
                                           index, rtil, TRUE, index[i_gl], index[i_gl + 1], dd, zones,
                                           &molb[mb], bRCheckMB, rcheck, bRCheck2B, rc2, pbc_null,
                                           cg_cm, ip_in, idef, izone, bBCheck, &nbonded_local);
        }
    }

    return nbonded_local;
}

/*! \brief Set the exclusion data for i-zone \p iz */
static void make_exclusions_zone(gmx_domdec_t*                     dd,
                                 gmx_domdec_zones_t*               zones,
                                 const std::vector<gmx_moltype_t>& moltype,
                                 const int*                        cginfo,
                                 ListOfLists<int>*                 lexcls,
                                 int                               iz,
                                 int                               at_start,
                                 int                               at_end,
                                 const gmx::ArrayRef<const int>    intermolecularExclusionGroup)
{
    const gmx_ga2la_t& ga2la = *dd->ga2la;

    const auto& jAtomRange = zones->iZones[iz].jAtomRange;

    const gmx::index oldNumLists = lexcls->ssize();

    std::vector<int> exclusionsForAtom;
    for (int at = at_start; at < at_end; at++)
    {
        exclusionsForAtom.clear();

        if (GET_CGINFO_EXCL_INTER(cginfo[at]))
        {
            int a_gl, mb, mt, mol, a_mol;

            /* Copy the exclusions from the global top */
            a_gl = dd->globalAtomIndices[at];
            global_atomnr_to_moltype_ind(dd->reverse_top, a_gl, &mb, &mt, &mol, &a_mol);
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
                                           intermolecularExclusionGroup.end(), dd->globalAtomIndices[at])
                                         != intermolecularExclusionGroup.end();

        if (isExcludedAtom)
        {
            for (int qmAtomGlobalIndex : intermolecularExclusionGroup)
            {
                if (const auto* entry = dd->ga2la->find(qmAtomGlobalIndex))
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
static int make_local_bondeds_excls(gmx_domdec_t*           dd,
                                    gmx_domdec_zones_t*     zones,
                                    const gmx_mtop_t*       mtop,
                                    const int*              cginfo,
                                    gmx_bool                bRCheckMB,
                                    ivec                    rcheck,
                                    gmx_bool                bRCheck2B,
                                    real                    rc,
                                    t_pbc*                  pbc_null,
                                    rvec*                   cg_cm,
                                    InteractionDefinitions* idef,
                                    ListOfLists<int>*       lexcls,
                                    int*                    excl_count)
{
    int                nzone_bondeds;
    int                cg0, cg1;
    real               rc2;
    int                nbonded_local;
    gmx_reverse_top_t* rt;

    if (dd->reverse_top->bInterAtomicInteractions)
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

    rt = dd->reverse_top;

    rc2 = rc * rc;

    /* Clear the counts */
    idef->clear();
    nbonded_local = 0;

    lexcls->clear();
    *excl_count = 0;

    for (int izone = 0; izone < nzone_bondeds; izone++)
    {
        cg0 = zones->cg_range[izone];
        cg1 = zones->cg_range[izone + 1];

        const int numThreads = rt->th_work.size();
#pragma omp parallel for num_threads(numThreads) schedule(static)
        for (int thread = 0; thread < numThreads; thread++)
        {
            try
            {
                int                     cg0t, cg1t;
                InteractionDefinitions* idef_t;

                cg0t = cg0 + ((cg1 - cg0) * thread) / numThreads;
                cg1t = cg0 + ((cg1 - cg0) * (thread + 1)) / numThreads;

                if (thread == 0)
                {
                    idef_t = idef;
                }
                else
                {
                    idef_t = &rt->th_work[thread].idef;
                    idef_t->clear();
                }

                rt->th_work[thread].nbonded = make_bondeds_zone(
                        dd, zones, mtop->molblock, bRCheckMB, rcheck, bRCheck2B, rc2, pbc_null,
                        cg_cm, idef->iparams.data(), idef_t, izone, gmx::Range<int>(cg0t, cg1t));

                if (izone < numIZonesForExclusions)
                {
                    ListOfLists<int>* excl_t;
                    if (thread == 0)
                    {
                        // Thread 0 stores exclusions directly in the final storage
                        excl_t = lexcls;
                    }
                    else
                    {
                        // Threads > 0 store in temporary storage, starting at list index 0
                        excl_t = &rt->th_work[thread].excl;
                        excl_t->clear();
                    }

                    /* No charge groups and no distance check required */
                    make_exclusions_zone(dd, zones, mtop->moltype, cginfo, excl_t, izone, cg0t,
                                         cg1t, mtop->intermolecularExclusionGroup);
                }
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
        }

        if (rt->th_work.size() > 1)
        {
            combine_idef(idef, rt->th_work);
        }

        for (const thread_work_t& th_work : rt->th_work)
        {
            nbonded_local += th_work.nbonded;
        }

        if (izone < numIZonesForExclusions)
        {
            for (std::size_t th = 1; th < rt->th_work.size(); th++)
            {
                lexcls->appendListOfLists(rt->th_work[th].excl);
            }
            for (const thread_work_t& th_work : rt->th_work)
            {
                *excl_count += th_work.excl_count;
            }
        }
    }

    if (debug)
    {
        fprintf(debug, "We have %d exclusions, check count %d\n", lexcls->numElements(), *excl_count);
    }

    return nbonded_local;
}

void dd_make_local_top(gmx_domdec_t*       dd,
                       gmx_domdec_zones_t* zones,
                       int                 npbcdim,
                       matrix              box,
                       rvec                cellsize_min,
                       const ivec          npulse,
                       t_forcerec*         fr,
                       rvec*               cgcm_or_x,
                       const gmx_mtop_t&   mtop,
                       gmx_localtop_t*     ltop)
{
    gmx_bool bRCheckMB, bRCheck2B;
    real     rc = -1;
    ivec     rcheck;
    int      d, nexcl;
    t_pbc    pbc, *pbc_null = nullptr;

    if (debug)
    {
        fprintf(debug, "Making local topology\n");
    }

    bRCheckMB = FALSE;
    bRCheck2B = FALSE;

    if (dd->reverse_top->bInterAtomicInteractions)
    {
        /* We need to check to which cell bondeds should be assigned */
        rc = dd_cutoff_twobody(dd);
        if (debug)
        {
            fprintf(debug, "Two-body bonded cut-off distance is %g\n", rc);
        }

        /* Should we check cg_cm distances when assigning bonded interactions? */
        for (d = 0; d < DIM; d++)
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
                    rcheck[d] = TRUE;
                    bRCheckMB = TRUE;
                }
                /* Check for interactions between two atoms,
                 * where we can allow interactions up to the cut-off,
                 * instead of up to the smallest cell dimension.
                 */
                bRCheck2B = TRUE;
            }
            if (debug)
            {
                fprintf(debug, "dim %d cellmin %f bonded rcheck[%d] = %d, bRCheck2B = %s\n", d,
                        cellsize_min[d], d, rcheck[d], gmx::boolToString(bRCheck2B));
            }
        }
        if (bRCheckMB || bRCheck2B)
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

    dd->nbonded_local = make_local_bondeds_excls(dd, zones, &mtop, fr->cginfo.data(), bRCheckMB,
                                                 rcheck, bRCheck2B, rc, pbc_null, cgcm_or_x,
                                                 &ltop->idef, &ltop->excls, &nexcl);

    /* The ilist is not sorted yet,
     * we can only do this when we have the charge arrays.
     */
    ltop->idef.ilsort = ilsortUNKNOWN;
}

void dd_sort_local_top(gmx_domdec_t* dd, const t_mdatoms* mdatoms, gmx_localtop_t* ltop)
{
    if (dd->reverse_top->ilsort == ilsortNO_FE)
    {
        ltop->idef.ilsort = ilsortNO_FE;
    }
    else
    {
        gmx_sort_ilist_fe(&ltop->idef, mdatoms->chargeA, mdatoms->chargeB);
    }
}

void dd_init_local_state(gmx_domdec_t* dd, const t_state* state_global, t_state* state_local)
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
    dd_bcast(dd, NITEM_DD_INIT_LOCAL_STATE * sizeof(int), buf);

    init_gtc_state(state_local, buf[1], buf[2], buf[3]);
    init_dfhist_state(state_local, buf[4]);
    state_local->flags = buf[0];
}

/*! \brief Check if a link is stored in \p link between charge groups \p cg_gl and \p cg_gl_j and if not so, store a link */
static void check_link(t_blocka* link, int cg_gl, int cg_gl_j)
{
    int      k;
    gmx_bool bFound;

    bFound = FALSE;
    for (k = link->index[cg_gl]; k < link->index[cg_gl + 1]; k++)
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
    t_blocka*    link;
    cginfo_mb_t* cgi_mb;

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

        make_reverse_ilist(*mtop.intermolecular_ilist, &atoms, FALSE, FALSE, FALSE, TRUE, &ril_intermol);
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
        reverse_ilist_t ril;
        make_reverse_ilist(molt.ilist, &molt.atoms, FALSE, FALSE, FALSE, TRUE, &ril);

        cgi_mb = &cginfo_mb[mb];

        int mol;
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
            fprintf(debug, "molecule type '%s' %d atoms has %d atom links through bonded interac.\n",
                    *molt.name, molt.atoms.nr, nlink_mol);
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

/*! \brief Set the distance, function type and atom indices for the longest distance between charge-groups of molecule type \p molt for two-body and multi-body bonded interactions */
static void bonded_cg_distance_mol(const gmx_moltype_t* molt,
                                   gmx_bool             bBCheck,
                                   gmx_bool             bExcl,
                                   rvec*                cg_cm,
                                   bonded_distance_t*   bd_2b,
                                   bonded_distance_t*   bd_mb)
{
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, bBCheck, FALSE, FALSE))
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
                                real rij2 = distance2(cg_cm[atomI], cg_cm[atomJ]);

                                update_max_bonded_distance(rij2, ftype, atomI, atomJ,
                                                           (nral == 2) ? bd_2b : bd_mb);
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
                    real rij2 = distance2(cg_cm[ai], cg_cm[aj]);

                    /* There is no function type for exclusions, use -1 */
                    update_max_bonded_distance(rij2, -1, ai, aj, bd_2b);
                }
            }
        }
    }
}

/*! \brief Set the distance, function type and atom indices for the longest atom distance involved in intermolecular interactions for two-body and multi-body bonded interactions */
static void bonded_distance_intermol(const InteractionLists& ilists_intermol,
                                     gmx_bool                bBCheck,
                                     const rvec*             x,
                                     PbcType                 pbcType,
                                     const matrix            box,
                                     bonded_distance_t*      bd_2b,
                                     bonded_distance_t*      bd_mb)
{
    t_pbc pbc;

    set_pbc(&pbc, pbcType, box);

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (dd_check_ftype(ftype, bBCheck, FALSE, FALSE))
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
                        real rij2;

                        int atom_j = il.iatoms[i + 1 + aj];

                        pbc_dx(&pbc, x[atom_i], x[atom_j], dx);

                        rij2 = norm2(dx);

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
                                        const rvec*           x,
                                        rvec*                 xs)
{
    int n, i;

    if (pbcType != PbcType::No)
    {
        mk_mshift(nullptr, graph, pbcType, box, x);

        shift_x(graph, box, x, xs);
        /* By doing an extra mk_mshift the molecules that are broken
         * because they were e.g. imported from another software
         * will be made whole again. Such are the healing powers
         * of GROMACS.
         */
        mk_mshift(nullptr, graph, pbcType, box, xs);
    }
    else
    {
        /* We copy the coordinates so the original coordinates remain
         * unchanged, just to be 100% sure that we do not affect
         * binary reproducibility of simulations.
         */
        n = molt->atoms.nr;
        for (i = 0; i < n; i++)
        {
            copy_rvec(x[i], xs[i]);
        }
    }

    if (moltypeHasVsite(*molt))
    {
        construct_vsites(nullptr, xs, 0.0, nullptr, ffparams->iparams, molt->ilist, PbcType::No,
                         TRUE, nullptr, nullptr);
    }
}

void dd_bonded_cg_distance(const gmx::MDLogger& mdlog,
                           const gmx_mtop_t*    mtop,
                           const t_inputrec*    ir,
                           const rvec*          x,
                           const matrix         box,
                           gmx_bool             bBCheck,
                           real*                r_2b,
                           real*                r_mb)
{
    gmx_bool          bExclRequired;
    int               at_offset;
    rvec*             xs;
    bonded_distance_t bd_2b = { 0, -1, -1, -1 };
    bonded_distance_t bd_mb = { 0, -1, -1, -1 };

    bExclRequired = inputrecExclForces(ir);

    *r_2b     = 0;
    *r_mb     = 0;
    at_offset = 0;
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const gmx_moltype_t& molt = mtop->moltype[molb.type];
        if (molt.atoms.nr == 1 || molb.nmol == 0)
        {
            at_offset += molb.nmol * molt.atoms.nr;
        }
        else
        {
            t_graph graph;
            if (ir->pbcType != PbcType::No)
            {
                graph = mk_graph_moltype(molt);
            }

            snew(xs, molt.atoms.nr);
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                getWholeMoleculeCoordinates(&molt, &mtop->ffparams, ir->pbcType, &graph, box,
                                            x + at_offset, xs);

                bonded_distance_t bd_mol_2b = { 0, -1, -1, -1 };
                bonded_distance_t bd_mol_mb = { 0, -1, -1, -1 };

                bonded_cg_distance_mol(&molt, bBCheck, bExclRequired, xs, &bd_mol_2b, &bd_mol_mb);

                /* Process the mol data adding the atom index offset */
                update_max_bonded_distance(bd_mol_2b.r2, bd_mol_2b.ftype, at_offset + bd_mol_2b.a1,
                                           at_offset + bd_mol_2b.a2, &bd_2b);
                update_max_bonded_distance(bd_mol_mb.r2, bd_mol_mb.ftype, at_offset + bd_mol_mb.a1,
                                           at_offset + bd_mol_mb.a2, &bd_mb);

                at_offset += molt.atoms.nr;
            }
            sfree(xs);
        }
    }

    if (mtop->bIntermolecularInteractions)
    {
        GMX_RELEASE_ASSERT(mtop->intermolecular_ilist,
                           "We should have an ilist when intermolecular interactions are on");

        bonded_distance_intermol(*mtop->intermolecular_ilist, bBCheck, x, ir->pbcType, box, &bd_2b, &bd_mb);
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
                            "    two-body bonded interactions: %5.3f nm, %s, atoms %d %d", *r_2b,
                            (bd_2b.ftype >= 0) ? interaction_function[bd_2b.ftype].longname : "Exclusion",
                            bd_2b.a1 + 1, bd_2b.a2 + 1);
        }
        if (*r_mb > 0)
        {
            GMX_LOG(mdlog.info)
                    .appendTextFormatted(
                            "  multi-body bonded interactions: %5.3f nm, %s, atoms %d %d", *r_mb,
                            interaction_function[bd_mb.ftype].longname, bd_mb.a1 + 1, bd_mb.a2 + 1);
        }
    }
}
