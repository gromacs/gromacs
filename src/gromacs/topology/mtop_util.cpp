/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2008- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/topology/mtop_util.h"

#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/exclusionblocks.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/topology/topsort.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/range.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

void gmx_mtop_count_atomtypes(const gmx_mtop_t& mtop, int state, int typecount[])
{
    for (int i = 0; i < mtop.ffparams.atnr; ++i)
    {
        typecount[i] = 0;
    }
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const t_atoms& atoms = mtop.moltype[molb.type].atoms;
        for (int i = 0; i < atoms.nr; ++i)
        {
            const int tpi = (state == 0) ? atoms.atom[i].type : atoms.atom[i].typeB;
            typecount[tpi] += molb.nmol;
        }
    }
}

int gmx_mtop_num_molecules(const gmx_mtop_t& mtop)
{
    int numMolecules = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        numMolecules += molb.nmol;
    }
    return numMolecules;
}

int gmx_mtop_nres(const gmx_mtop_t& mtop)
{
    int nres = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        nres += molb.nmol * mtop.moltype[molb.type].atoms.nres;
    }
    return nres;
}

int gmx_mtop_ftype_count(const gmx_mtop_t& mtop, int ftype)
{
    int n = 0;

    for (const IListProxy il : IListRange(mtop))
    {
        n += il.nmol() * il.list()[ftype].size() / (1 + NRAL(ftype));
    }
    return n;
}

int gmx_mtop_interaction_count(const gmx_mtop_t& mtop, const int unsigned if_flags)
{
    int n = 0;

    for (const IListProxy il : IListRange(mtop))
    {
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if ((interaction_function[ftype].flags & if_flags) == if_flags)
            {
                n += il.nmol() * il.list()[ftype].size() / (1 + NRAL(ftype));
            }
        }
    }
    return n;
}

gmx::EnumerationArray<ParticleType, int> gmx_mtop_particletype_count(const gmx_mtop_t& mtop)
{
    gmx::EnumerationArray<ParticleType, int> count = { { 0 } };

    for (const auto& molblock : mtop.molblock)
    {
        const t_atoms& atoms = mtop.moltype[molblock.type].atoms;
        for (int a = 0; a < atoms.nr; a++)
        {
            count[atoms.atom[a].ptype] += molblock.nmol;
        }
    }

    return count;
}

static void atomcat(t_atoms* dest, const t_atoms* src, int copies, int maxres_renum, int* maxresnr)
{
    int i = 0, j = 0, l = 0, size = 0;
    int srcnr  = src->nr;
    int destnr = dest->nr;

    if (dest->nr == 0)
    {
        dest->haveMass    = src->haveMass;
        dest->haveType    = src->haveType;
        dest->haveCharge  = src->haveCharge;
        dest->haveBState  = src->haveBState;
        dest->havePdbInfo = src->havePdbInfo;
    }
    else
    {
        dest->haveMass    = dest->haveMass && src->haveMass;
        dest->haveType    = dest->haveType && src->haveType;
        dest->haveCharge  = dest->haveCharge && src->haveCharge;
        dest->haveBState  = dest->haveBState && src->haveBState;
        dest->havePdbInfo = dest->havePdbInfo && src->havePdbInfo;
    }

    if (srcnr)
    {
        size = destnr + copies * srcnr;
        srenew(dest->atom, size);
        srenew(dest->atomname, size);
        if (dest->haveType)
        {
            srenew(dest->atomtype, size);
            if (dest->haveBState)
            {
                srenew(dest->atomtypeB, size);
            }
        }
        if (dest->havePdbInfo)
        {
            srenew(dest->pdbinfo, size);
        }
    }
    if (src->nres)
    {
        size = dest->nres + copies * src->nres;
        srenew(dest->resinfo, size);
    }

    /* residue information */
    for (l = dest->nres, j = 0; (j < copies); j++, l += src->nres)
    {
        memcpy(reinterpret_cast<char*>(&(dest->resinfo[l])),
               reinterpret_cast<char*>(&(src->resinfo[0])),
               static_cast<size_t>(src->nres * sizeof(src->resinfo[0])));
    }

    for (l = destnr, j = 0; (j < copies); j++, l += srcnr)
    {
        memcpy(reinterpret_cast<char*>(&(dest->atom[l])),
               reinterpret_cast<char*>(&(src->atom[0])),
               static_cast<size_t>(srcnr * sizeof(src->atom[0])));
        memcpy(reinterpret_cast<char*>(&(dest->atomname[l])),
               reinterpret_cast<char*>(&(src->atomname[0])),
               static_cast<size_t>(srcnr * sizeof(src->atomname[0])));
        if (dest->haveType)
        {
            memcpy(reinterpret_cast<char*>(&(dest->atomtype[l])),
                   reinterpret_cast<char*>(&(src->atomtype[0])),
                   static_cast<size_t>(srcnr * sizeof(src->atomtype[0])));
            if (dest->haveBState)
            {
                memcpy(reinterpret_cast<char*>(&(dest->atomtypeB[l])),
                       reinterpret_cast<char*>(&(src->atomtypeB[0])),
                       static_cast<size_t>(srcnr * sizeof(src->atomtypeB[0])));
            }
        }
        if (dest->havePdbInfo)
        {
            memcpy(reinterpret_cast<char*>(&(dest->pdbinfo[l])),
                   reinterpret_cast<char*>(&(src->pdbinfo[0])),
                   static_cast<size_t>(srcnr * sizeof(src->pdbinfo[0])));
        }
    }

    /* Increment residue indices */
    for (l = destnr, j = 0; (j < copies); j++)
    {
        for (i = 0; (i < srcnr); i++, l++)
        {
            dest->atom[l].resind = dest->nres + j * src->nres + src->atom[i].resind;
        }
    }

    if (src->nres <= maxres_renum)
    {
        /* Single residue molecule, continue counting residues */
        for (j = 0; (j < copies); j++)
        {
            for (l = 0; l < src->nres; l++)
            {
                (*maxresnr)++;
                dest->resinfo[dest->nres + j * src->nres + l].nr = *maxresnr;
            }
        }
    }

    dest->nres += copies * src->nres;
    dest->nr += copies * src->nr;
}

t_atoms gmx_mtop_global_atoms(const gmx_mtop_t& mtop)
{
    t_atoms atoms;

    init_t_atoms(&atoms, 0, FALSE);

    int maxresnr = mtop.maxResNumberNotRenumbered();
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        atomcat(&atoms,
                &mtop.moltype[molb.type].atoms,
                molb.nmol,
                mtop.maxResiduesPerMoleculeToTriggerRenumber(),
                &maxresnr);
    }

    return atoms;
}

/*
 * The cat routines below are old code from src/kernel/topcat.c
 */

static void ilistcat(int ftype, InteractionList* dest, const InteractionList& src, int copies, int dnum, int snum)
{
    const int nral = NRAL(ftype);

    size_t destIndex = dest->iatoms.size();
    dest->iatoms.resize(dest->iatoms.size() + copies * src.size());

    for (int c = 0; c < copies; c++)
    {
        for (int i = 0; i < src.size();)
        {
            dest->iatoms[destIndex++] = src.iatoms[i++];
            for (int a = 0; a < nral; a++)
            {
                dest->iatoms[destIndex++] = dnum + src.iatoms[i++];
            }
        }
        dnum += snum;
    }
}

static void ilistcat(int ftype, t_ilist* dest, const InteractionList& src, int copies, int dnum, int snum)
{
    const int nral = NRAL(ftype);

    dest->nalloc = dest->nr + copies * src.size();
    srenew(dest->iatoms, dest->nalloc);

    for (int c = 0; c < copies; c++)
    {
        for (int i = 0; i < src.size();)
        {
            dest->iatoms[dest->nr++] = src.iatoms[i++];
            for (int a = 0; a < nral; a++)
            {
                dest->iatoms[dest->nr++] = dnum + src.iatoms[i++];
            }
        }
        dnum += snum;
    }
}

static const t_iparams& getIparams(const InteractionDefinitions& idef, const int index)
{
    return idef.iparams[index];
}

static const t_iparams& getIparams(const t_idef& idef, const int index)
{
    return idef.iparams[index];
}

static void resizeIParams(std::vector<t_iparams>* iparams, const int newSize)
{
    iparams->resize(newSize);
}

static void resizeIParams(t_iparams** iparams, const int newSize)
{
    srenew(*iparams, newSize);
}

template<typename IdefType>
static void set_posres_params(IdefType* idef, const gmx_molblock_t* molb, int i0, int a_offset)
{
    auto* il = &idef->il[F_POSRES];
    int   i1 = il->size() / 2;
    resizeIParams(&idef->iparams_posres, i1);
    for (int i = i0; i < i1; i++)
    {
        t_iparams* ip = &idef->iparams_posres[i];
        /* Copy the force constants */
        *ip        = getIparams(*idef, il->iatoms[i * 2]);
        int a_molb = il->iatoms[i * 2 + 1] - a_offset;
        if (molb->posres_xA.empty())
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        ip->posres.pos0A[XX] = molb->posres_xA[a_molb][XX];
        ip->posres.pos0A[YY] = molb->posres_xA[a_molb][YY];
        ip->posres.pos0A[ZZ] = molb->posres_xA[a_molb][ZZ];
        if (!molb->posres_xB.empty())
        {
            ip->posres.pos0B[XX] = molb->posres_xB[a_molb][XX];
            ip->posres.pos0B[YY] = molb->posres_xB[a_molb][YY];
            ip->posres.pos0B[ZZ] = molb->posres_xB[a_molb][ZZ];
        }
        else
        {
            ip->posres.pos0B[XX] = ip->posres.pos0A[XX];
            ip->posres.pos0B[YY] = ip->posres.pos0A[YY];
            ip->posres.pos0B[ZZ] = ip->posres.pos0A[ZZ];
        }
        /* Set the parameter index for idef->iparams_posre */
        il->iatoms[i * 2] = i;
    }
}

template<typename IdefType>
static void set_fbposres_params(IdefType* idef, const gmx_molblock_t* molb, int i0, int a_offset)
{
    auto* il = &idef->il[F_FBPOSRES];
    int   i1 = il->size() / 2;
    resizeIParams(&idef->iparams_fbposres, i1);
    for (int i = i0; i < i1; i++)
    {
        t_iparams* ip = &idef->iparams_fbposres[i];
        /* Copy the force constants */
        *ip        = getIparams(*idef, il->iatoms[i * 2]);
        int a_molb = il->iatoms[i * 2 + 1] - a_offset;
        if (molb->posres_xA.empty())
        {
            gmx_incons("Position restraint coordinates are missing");
        }
        /* Take flat-bottom posres reference from normal position restraints */
        ip->fbposres.pos0[XX] = molb->posres_xA[a_molb][XX];
        ip->fbposres.pos0[YY] = molb->posres_xA[a_molb][YY];
        ip->fbposres.pos0[ZZ] = molb->posres_xA[a_molb][ZZ];
        /* Note: no B-type for flat-bottom posres */

        /* Set the parameter index for idef->iparams_posre */
        il->iatoms[i * 2] = i;
    }
}

/*! \brief Copy parameters to idef structure from mtop.
 *
 * Makes a deep copy of the force field parameters data structure from a gmx_mtop_t.
 * Used to initialize legacy topology types.
 *
 * \param[in] mtop Reference to input mtop.
 * \param[in] idef Pointer to idef to populate.
 */
static void copyFFParametersFromMtop(const gmx_mtop_t& mtop, t_idef* idef)
{
    const gmx_ffparams_t* ffp = &mtop.ffparams;

    idef->ntypes = ffp->numTypes();
    idef->atnr   = ffp->atnr;
    /* we can no longer copy the pointers to the mtop members,
     * because they will become invalid as soon as mtop gets free'd.
     * We also need to make sure to only operate on valid data!
     */

    if (!ffp->functype.empty())
    {
        snew(idef->functype, ffp->functype.size());
        std::copy(ffp->functype.data(), ffp->functype.data() + ffp->functype.size(), idef->functype);
    }
    else
    {
        idef->functype = nullptr;
    }
    if (!ffp->iparams.empty())
    {
        snew(idef->iparams, ffp->iparams.size());
        std::copy(ffp->iparams.data(), ffp->iparams.data() + ffp->iparams.size(), idef->iparams);
    }
    else
    {
        idef->iparams = nullptr;
    }
    idef->iparams_posres   = nullptr;
    idef->iparams_fbposres = nullptr;
    idef->fudgeQQ          = ffp->fudgeQQ;
    idef->ilsort           = ilsortUNKNOWN;
}

/*! \brief Copy idef structure from mtop.
 *
 * Makes a deep copy of an idef data structure from a gmx_mtop_t.
 * Used to initialize legacy topology types.
 *
 * \param[in] mtop Reference to input mtop.
 * \param[in] idef Pointer to idef to populate.
 * \param[in] mergeConstr Decide if constraints will be merged.
 */
template<typename IdefType>
static void copyIListsFromMtop(const gmx_mtop_t& mtop, IdefType* idef, bool mergeConstr)
{
    int natoms = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        int srcnr  = molt.atoms.nr;
        int destnr = natoms;

        int nposre_old   = idef->il[F_POSRES].size();
        int nfbposre_old = idef->il[F_FBPOSRES].size();
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (mergeConstr && ftype == F_CONSTR && !molt.ilist[F_CONSTRNC].empty())
            {
                /* Merge all constrains into one ilist.
                 * This simplifies the constraint code.
                 */
                for (int mol = 0; mol < molb.nmol; mol++)
                {
                    ilistcat(ftype, &idef->il[F_CONSTR], molt.ilist[F_CONSTR], 1, destnr + mol * srcnr, srcnr);
                    ilistcat(ftype, &idef->il[F_CONSTR], molt.ilist[F_CONSTRNC], 1, destnr + mol * srcnr, srcnr);
                }
            }
            else if (!(mergeConstr && ftype == F_CONSTRNC))
            {
                ilistcat(ftype, &idef->il[ftype], molt.ilist[ftype], molb.nmol, destnr, srcnr);
            }
        }
        if (idef->il[F_POSRES].size() > nposre_old)
        {
            /* Executing this line line stops gmxdump -sys working
             * correctly. I'm not aware there's an elegant fix. */
            set_posres_params(idef, &molb, nposre_old / 2, natoms);
        }
        if (idef->il[F_FBPOSRES].size() > nfbposre_old)
        {
            set_fbposres_params(idef, &molb, nfbposre_old / 2, natoms);
        }

        natoms += molb.nmol * srcnr;
    }

    if (mtop.bIntermolecularInteractions)
    {
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            ilistcat(ftype, &idef->il[ftype], (*mtop.intermolecular_ilist)[ftype], 1, 0, mtop.natoms);
        }
    }

    // We have not (yet) sorted free-energy interactions to the end of the ilists
    idef->ilsort = ilsortNO_FE;
}

/*! \brief Generate a single list of lists of exclusions for the whole system
 *
 * \param[in] mtop  Reference to input mtop.
 */
static gmx::ListOfLists<int> globalExclusionLists(const gmx_mtop_t& mtop)
{
    gmx::ListOfLists<int> excls;

    int atomIndex = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        const gmx_moltype_t& molt = mtop.moltype[molb.type];

        for (int mol = 0; mol < molb.nmol; mol++)
        {
            excls.appendListOfLists(molt.excls, atomIndex);

            atomIndex += molt.atoms.nr;
        }
    }

    return excls;
}

/*! \brief Updates inter-molecular exclusion lists
 *
 * This function updates inter-molecular exclusions to exclude all
 * non-bonded interactions between a given list of atoms
 *
 * \param[inout]    excls   existing exclusions in local topology
 * \param[in]       ids     list of global IDs of atoms
 */
static void addMimicExclusions(gmx::ListOfLists<int>* excls, const gmx::ArrayRef<const int> ids)
{
    t_blocka inter_excl{};
    init_blocka(&inter_excl);
    size_t n_q = ids.size();

    inter_excl.nr  = excls->ssize();
    inter_excl.nra = n_q * n_q;

    size_t total_nra = n_q * n_q;

    snew(inter_excl.index, excls->ssize() + 1);
    snew(inter_excl.a, total_nra);

    for (int i = 0; i < inter_excl.nr; ++i)
    {
        inter_excl.index[i] = 0;
    }

    /* Here we loop over the list of QM atom ids
     *  and create exclusions between all of them resulting in
     *  n_q * n_q sized exclusion list
     */
    int prev_index = 0;
    for (int k = 0; k < inter_excl.nr; ++k)
    {
        inter_excl.index[k] = prev_index;
        for (long i = 0; i < ids.ssize(); i++)
        {
            if (k != ids[i])
            {
                continue;
            }
            size_t index             = n_q * i;
            inter_excl.index[ids[i]] = index;
            prev_index               = index + n_q;
            for (size_t j = 0; j < n_q; ++j)
            {
                inter_excl.a[n_q * i + j] = ids[j];
            }
        }
    }
    inter_excl.index[ids[n_q - 1] + 1] = n_q * n_q;

    inter_excl.index[inter_excl.nr] = n_q * n_q;

    std::vector<gmx::ExclusionBlock> qmexcl2(excls->size());
    gmx::blockaToExclusionBlocks(&inter_excl, qmexcl2);

    // Merge the created exclusion list with the existing one
    gmx::mergeExclusions(excls, qmexcl2);
}

static void sortFreeEnergyInteractionsAtEnd(const gmx_mtop_t& mtop, InteractionDefinitions* idef)
{
    std::vector<int32_t> atomInfo(mtop.natoms, 0);


    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        const gmx_moltype_t&  molt = mtop.moltype[molb.type];
        if (!molt.ilist[F_LJ14].iatoms.empty())
        {
            int offset = mtop.moleculeBlockIndices[mb].globalAtomStart;
            for (int mol = 0; mol < molb.nmol; mol++)
            {
                for (int a = 0; a < molt.atoms.nr; a++)
                {
                    if (atomHasPerturbedCharge(molt.atoms.atom[a]))
                    {
                        atomInfo[offset + a] |= gmx::sc_atomInfo_HasPerturbedCharge;
                    }
                }
                offset += molt.atoms.nr;
            }
        }
    }
    gmx_sort_ilist_fe(idef, atomInfo);
}

static void gen_local_top(const gmx_mtop_t& mtop,
                          bool              freeEnergyInteractionsAtEnd,
                          bool              bMergeConstr,
                          gmx_localtop_t*   top)
{
    copyIListsFromMtop(mtop, &top->idef, bMergeConstr);
    if (freeEnergyInteractionsAtEnd)
    {
        sortFreeEnergyInteractionsAtEnd(mtop, &top->idef);
    }
    top->excls = globalExclusionLists(mtop);
    if (!mtop.intermolecularExclusionGroup.empty())
    {
        addMimicExclusions(&top->excls, mtop.intermolecularExclusionGroup);
    }
}

void gmx_mtop_generate_local_top(const gmx_mtop_t& mtop, gmx_localtop_t* top, bool freeEnergyInteractionsAtEnd)
{
    gen_local_top(mtop, freeEnergyInteractionsAtEnd, true, top);
}

/*! \brief Fills an array with molecule begin/end atom indices
 *
 * \param[in]  mtop   The global topology
 * \param[out] index  Array of size nr. of molecules + 1 to be filled with molecule begin/end indices
 */
static void fillMoleculeIndices(const gmx_mtop_t& mtop, gmx::ArrayRef<int> index)
{
    int globalAtomIndex   = 0;
    int globalMolIndex    = 0;
    index[globalMolIndex] = globalAtomIndex;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int numAtomsPerMolecule = mtop.moltype[molb.type].atoms.nr;
        for (int mol = 0; mol < molb.nmol; mol++)
        {
            globalAtomIndex += numAtomsPerMolecule;
            globalMolIndex += 1;
            index[globalMolIndex] = globalAtomIndex;
        }
    }
}

gmx::RangePartitioning gmx_mtop_molecules(const gmx_mtop_t& mtop)
{
    gmx::RangePartitioning mols;

    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int numAtomsPerMolecule = mtop.moltype[molb.type].atoms.nr;
        for (int mol = 0; mol < molb.nmol; mol++)
        {
            mols.appendBlock(numAtomsPerMolecule);
        }
    }

    return mols;
}

std::vector<gmx::Range<int>> atomRangeOfEachResidue(const gmx_moltype_t& moltype)
{
    std::vector<gmx::Range<int>> atomRanges;
    int                          currentResidueNumber = moltype.atoms.atom[0].resind;
    int                          startAtom            = 0;
    // Go through all atoms in a molecule to store first and last atoms in each residue.
    for (int i = 0; i < moltype.atoms.nr; i++)
    {
        int residueOfThisAtom = moltype.atoms.atom[i].resind;
        if (residueOfThisAtom != currentResidueNumber)
        {
            // This atom belongs to the next residue, so record the range for the previous residue,
            // remembering that end points to one place past the last atom.
            int endAtom = i;
            atomRanges.emplace_back(startAtom, endAtom);
            // Prepare for the current residue
            startAtom            = endAtom;
            currentResidueNumber = residueOfThisAtom;
        }
    }
    // special treatment for last residue in this molecule.
    atomRanges.emplace_back(startAtom, moltype.atoms.nr);

    return atomRanges;
}

/*! \brief Creates and returns a deprecated t_block struct with molecule indices
 *
 * \param[in] mtop  The global topology
 */
static t_block gmx_mtop_molecules_t_block(const gmx_mtop_t& mtop)
{
    t_block mols;

    mols.nr           = gmx_mtop_num_molecules(mtop);
    mols.nalloc_index = mols.nr + 1;
    snew(mols.index, mols.nalloc_index);

    fillMoleculeIndices(mtop, gmx::arrayRefFromArray(mols.index, mols.nr + 1));

    return mols;
}

static void gen_t_topology(const gmx_mtop_t& mtop, bool bMergeConstr, t_topology* top)
{
    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        top->idef.il[ftype].nr     = 0;
        top->idef.il[ftype].nalloc = 0;
        top->idef.il[ftype].iatoms = nullptr;
    }
    copyFFParametersFromMtop(mtop, &top->idef);
    copyIListsFromMtop(mtop, &top->idef, bMergeConstr);

    top->name                        = mtop.name;
    top->atoms                       = gmx_mtop_global_atoms(mtop);
    top->mols                        = gmx_mtop_molecules_t_block(mtop);
    top->bIntermolecularInteractions = mtop.bIntermolecularInteractions;
    top->symtab                      = mtop.symtab;
}

t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t* mtop, bool freeMTop)
{
    t_topology top;

    gen_t_topology(*mtop, false, &top);

    if (freeMTop)
    {
        // Clear pointers and counts, such that the pointers copied to top
        // keep pointing to valid data after destroying mtop.
        mtop->symtab.symbuf = nullptr;
        mtop->symtab.nr     = 0;
    }
    return top;
}

std::vector<int> get_atom_index(const gmx_mtop_t& mtop)
{

    std::vector<int> atom_index;
    for (const AtomProxy atomP : AtomRange(mtop))
    {
        const t_atom& local = atomP.atom();
        int           index = atomP.globalAtomNumber();
        if (local.ptype == ParticleType::Atom)
        {
            atom_index.push_back(index);
        }
    }
    return atom_index;
}

void convertAtomsToMtop(t_symtab* symtab, char** name, t_atoms* atoms, gmx_mtop_t* mtop)
{
    mtop->symtab = *symtab;

    mtop->name = name;

    mtop->moltype.clear();
    mtop->moltype.resize(1);
    mtop->moltype.back().atoms = *atoms;

    mtop->molblock.resize(1);
    mtop->molblock[0].type = 0;
    mtop->molblock[0].nmol = 1;

    mtop->bIntermolecularInteractions = FALSE;

    mtop->natoms = atoms->nr;

    mtop->haveMoleculeIndices = false;

    mtop->finalize();
}

bool haveFepPerturbedNBInteractions(const gmx_mtop_t& mtop)
{
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        for (int a = 0; a < molt.atoms.nr; a++)
        {
            if (PERTURBED(molt.atoms.atom[a]))
            {
                return true;
            }
        }
    }

    return false;
}

bool haveFepPerturbedMasses(const gmx_mtop_t& mtop)
{
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        for (int a = 0; a < molt.atoms.nr; a++)
        {
            const t_atom& atom = molt.atoms.atom[a];
            if (atom.m != atom.mB)
            {
                return true;
            }
        }
    }

    return false;
}

bool haveFepPerturbedMassesInSettles(const gmx_mtop_t& mtop)
{
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        if (!molt.ilist[F_SETTLE].empty())
        {
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                const t_atom& atom = molt.atoms.atom[a];
                if (atom.m != atom.mB)
                {
                    return true;
                }
            }
        }
    }
    return false;
}

bool havePerturbedConstraints(const gmx_mtop_t& mtop)
{
    // This code assumes that all perturbed constraints parameters are actually used
    const auto& ffparams = mtop.ffparams;

    for (gmx::Index i = 0; i < gmx::ssize(ffparams.functype); i++)
    {
        if (ffparams.functype[i] == F_CONSTR || ffparams.functype[i] == F_CONSTRNC)
        {
            const auto& iparams = ffparams.iparams[i];
            if (iparams.constr.dA != iparams.constr.dB)
            {
                return true;
            }
        }
    }

    return false;
}
