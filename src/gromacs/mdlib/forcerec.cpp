/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#include "forcerec.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <array>
#include <bitset>
#include <filesystem>
#include <memory>

#include "gromacs/commandline/filenm.h"
#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/ewald/ewald.h"
#include "gromacs/ewald/pme_pp_comm_gpu.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/gpu_utils/gpu_utils.h"
#include "gromacs/hardware/hw_info.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/listed_forces/listed_forces_gpu.h"
#include "gromacs/listed_forces/pairs.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/dispersioncorrection.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/md_support.h"
#include "gromacs/mdlib/wall.h"
#include "gromacs/mdlib/wholemoleculetransform.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/fcdata.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "gpuforcereduction.h"
#include "mdgraph_gpu.h"

ForceHelperBuffers::ForceHelperBuffers(bool haveDirectVirialContributions) :
    haveDirectVirialContributions_(haveDirectVirialContributions)
{
    shiftForces_.resize(gmx::c_numShiftVectors);
}

void ForceHelperBuffers::resize(int numAtoms)
{
    if (haveDirectVirialContributions_)
    {
        forceBufferForDirectVirialContributions_.resize(numAtoms);
    }
}

std::vector<real> makeNonBondedParameterLists(const int                      numAtomTypes,
                                              gmx::ArrayRef<const t_iparams> iparams,
                                              bool                           useBuckinghamPotential)
{
    std::vector<real> nbfp;

    if (useBuckinghamPotential)
    {
        nbfp.resize(3 * numAtomTypes * numAtomTypes);
        int k = 0;
        for (int i = 0; (i < numAtomTypes); i++)
        {
            for (int j = 0; (j < numAtomTypes); j++, k++)
            {
                BHAMA(nbfp, numAtomTypes, i, j) = iparams[k].bham.a;
                BHAMB(nbfp, numAtomTypes, i, j) = iparams[k].bham.b;
                /* nbfp now includes the 6.0 derivative prefactor */
                BHAMC(nbfp, numAtomTypes, i, j) = iparams[k].bham.c * 6.0;
            }
        }
    }
    else
    {
        nbfp.resize(2 * numAtomTypes * numAtomTypes);
        int k = 0;
        for (int i = 0; (i < numAtomTypes); i++)
        {
            for (int j = 0; (j < numAtomTypes); j++, k++)
            {
                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                C6(nbfp, numAtomTypes, i, j)  = iparams[k].lj.c6 * 6.0;
                C12(nbfp, numAtomTypes, i, j) = iparams[k].lj.c12 * 12.0;
            }
        }
    }

    return nbfp;
}

std::vector<real> makeLJPmeC6GridCorrectionParameters(const int                      numAtomTypes,
                                                      gmx::ArrayRef<const t_iparams> iparams,
                                                      LongRangeVdW ljpme_combination_rule)
{
    /* For LJ-PME simulations, we correct the energies with the reciprocal space
     * inside of the cut-off. To do this the non-bonded kernels needs to have
     * access to the C6-values used on the reciprocal grid in pme.c
     */

    std::vector<real> grid(2 * numAtomTypes * numAtomTypes, 0.0);
    int               k = 0;
    for (int i = 0; (i < numAtomTypes); i++)
    {
        for (int j = 0; (j < numAtomTypes); j++, k++)
        {
            real c6i  = iparams[i * (numAtomTypes + 1)].lj.c6;
            real c12i = iparams[i * (numAtomTypes + 1)].lj.c12;
            real c6j  = iparams[j * (numAtomTypes + 1)].lj.c6;
            real c12j = iparams[j * (numAtomTypes + 1)].lj.c12;
            real c6   = std::sqrt(c6i * c6j);
            if (ljpme_combination_rule == LongRangeVdW::LB && !gmx_numzero(c6) && !gmx_numzero(c12i)
                && !gmx_numzero(c12j))
            {
                real sigmai = gmx::sixthroot(c12i / c6i);
                real sigmaj = gmx::sixthroot(c12j / c6j);
                real epsi   = c6i * c6i / c12i;
                real epsj   = c6j * c6j / c12j;
                c6          = std::sqrt(epsi * epsj) * gmx::power6(0.5 * (sigmai + sigmaj));
            }
            /* Store the elements at the same relative positions as C6 in nbfp in order
             * to simplify access in the kernels
             */
            grid[2 * (numAtomTypes * i + j)] = c6 * 6.0;
        }
    }
    return grid;
}

//! What kind of constraint affects an atom
enum class ConstraintTypeForAtom : int
{
    None,       //!< No constraint active
    Constraint, //!< F_CONSTR or F_CONSTRNC active
    Settle,     //! F_SETTLE active
};

static std::vector<gmx::AtomInfoWithinMoleculeBlock>
makeAtomInfoForEachMoleculeBlock(const gmx_mtop_t& mtop, const t_forcerec* fr)
{
    std::vector<bool> atomUsesVdw(fr->ntype, false);
    for (int ai = 0; ai < fr->ntype; ai++)
    {
        for (int j = 0; j < fr->ntype; j++)
        {
            atomUsesVdw[ai] = atomUsesVdw[ai] || fr->haveBuckingham || C6(fr->nbfp, fr->ntype, ai, j) != 0
                              || C12(fr->nbfp, fr->ntype, ai, j) != 0;
        }
    }

    std::vector<SimulationAtomGroupType> groupTypesToCheck;
    if (mtop.groups.numberOfGroupNumbers(SimulationAtomGroupType::EnergyOutput) > 0)
    {
        groupTypesToCheck.push_back(SimulationAtomGroupType::EnergyOutput);
    }
    if (mtop.groups.numberOfGroupNumbers(SimulationAtomGroupType::QuantumMechanics) > 0)
    {
        groupTypesToCheck.push_back(SimulationAtomGroupType::QuantumMechanics);
    }

    std::vector<gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock;
    int                                           indexOfFirstAtomInMoleculeBlock = 0;
    for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
    {
        const gmx_molblock_t& molb = mtop.molblock[mb];
        const gmx_moltype_t&  molt = mtop.moltype[molb.type];
        const auto&           excl = molt.excls;

        /* Check if all molecules in this block have identical
         * atominfo. (That's true unless some kind of group- or
         * distance-based algorithm is involved, e.g. QM/MM groups
         * affecting multiple molecules within a block differently.)
         * If so, we only need an array of the size of one molecule.
         * Otherwise we make an array of #mol times #atoms per
         * molecule.
         */
        bool allMoleculesWithinBlockAreIdentical = true;
        for (const auto groupType : groupTypesToCheck)
        {
            for (int m = 0; m < molb.nmol; m++)
            {
                const int numAtomsInAllMolecules = m * molt.atoms.nr;
                for (int a = 0; a < molt.atoms.nr; a++)
                {
                    if (mtop.groups.groupNumbers[groupType][indexOfFirstAtomInMoleculeBlock + numAtomsInAllMolecules + a]
                        != mtop.groups.groupNumbers[groupType][indexOfFirstAtomInMoleculeBlock + a])
                    {
                        allMoleculesWithinBlockAreIdentical = false;
                    }
                }
            }
        }

        gmx::AtomInfoWithinMoleculeBlock atomInfoOfMoleculeBlock;
        atomInfoOfMoleculeBlock.indexOfFirstAtomInMoleculeBlock = indexOfFirstAtomInMoleculeBlock;
        atomInfoOfMoleculeBlock.indexOfLastAtomInMoleculeBlock =
                indexOfFirstAtomInMoleculeBlock + molb.nmol * molt.atoms.nr;
        int atomInfoSize = (allMoleculesWithinBlockAreIdentical ? 1 : molb.nmol) * molt.atoms.nr;
        atomInfoOfMoleculeBlock.atomInfo.resize(atomInfoSize);

        /* Set constraints flags for constrained atoms */
        std::vector<ConstraintTypeForAtom> constraintTypeOfAtom(molt.atoms.nr, ConstraintTypeForAtom::None);
        for (int ftype = 0; ftype < F_NRE; ftype++)
        {
            if (interaction_function[ftype].flags & IF_CONSTRAINT)
            {
                const int nral = NRAL(ftype);
                for (int ia = 0; ia < molt.ilist[ftype].size(); ia += 1 + nral)
                {
                    for (int a = 0; a < nral; a++)
                    {
                        constraintTypeOfAtom[molt.ilist[ftype].iatoms[ia + 1 + a]] =
                                (ftype == F_SETTLE ? ConstraintTypeForAtom::Settle
                                                   : ConstraintTypeForAtom::Constraint);
                    }
                }
            }
        }

        for (int m = 0; m < (allMoleculesWithinBlockAreIdentical ? 1 : molb.nmol); m++)
        {
            const int moleculeOffsetInBlock = m * molt.atoms.nr;
            for (int a = 0; a < molt.atoms.nr; a++)
            {
                const t_atom& atom = molt.atoms.atom[a];
                int32_t& atomInfo  = atomInfoOfMoleculeBlock.atomInfo[moleculeOffsetInBlock + a];

                /* Store the energy group in atomInfo */
                int gid  = getGroupType(mtop.groups,
                                       SimulationAtomGroupType::EnergyOutput,
                                       indexOfFirstAtomInMoleculeBlock + moleculeOffsetInBlock + a);
                atomInfo = (atomInfo & ~gmx::sc_atomInfo_EnergyGroupIdMask) | gid;

                bool bHaveVDW = (atomUsesVdw[atom.type] || atomUsesVdw[atom.typeB]);
                bool bHaveQ   = (atom.q != 0 || atom.qB != 0);

                bool haveExclusions = false;
                /* Loop over all the exclusions of atom ai */
                for (const int j : excl[a])
                {
                    if (j != a)
                    {
                        haveExclusions = true;
                        break;
                    }
                }

                switch (constraintTypeOfAtom[a])
                {
                    case ConstraintTypeForAtom::Constraint:
                        atomInfo |= gmx::sc_atomInfo_Constraint;
                        break;
                    case ConstraintTypeForAtom::Settle: atomInfo |= gmx::sc_atomInfo_Settle; break;
                    default: break;
                }
                if (haveExclusions)
                {
                    atomInfo |= gmx::sc_atomInfo_Exclusion;
                }
                if (bHaveVDW)
                {
                    atomInfo |= gmx::sc_atomInfo_HasVdw;
                }
                if (bHaveQ)
                {
                    atomInfo |= gmx::sc_atomInfo_HasCharge;
                }
                if (fr->efep != FreeEnergyPerturbationType::No)
                {
                    if (PERTURBED(atom))
                    {
                        atomInfo |= gmx::sc_atomInfo_FreeEnergyPerturbation;
                    }
                    if (atomHasPerturbedCharge(atom))
                    {
                        atomInfo |= gmx::sc_atomInfo_HasPerturbedCharge;
                    }
                }
            }
        }

        atomInfoForEachMoleculeBlock.push_back(atomInfoOfMoleculeBlock);

        indexOfFirstAtomInMoleculeBlock += molb.nmol * molt.atoms.nr;
    }

    return atomInfoForEachMoleculeBlock;
}

static std::vector<int32_t> expandAtomInfo(const int nmb,
                                           gmx::ArrayRef<const gmx::AtomInfoWithinMoleculeBlock> atomInfoForEachMoleculeBlock)
{
    const int numAtoms = atomInfoForEachMoleculeBlock[nmb - 1].indexOfLastAtomInMoleculeBlock;

    std::vector<int32_t> atomInfo(numAtoms);

    int mb = 0;
    for (int a = 0; a < numAtoms; a++)
    {
        while (a >= atomInfoForEachMoleculeBlock[mb].indexOfLastAtomInMoleculeBlock)
        {
            mb++;
        }
        atomInfo[a] = atomInfoForEachMoleculeBlock[mb]
                              .atomInfo[(a - atomInfoForEachMoleculeBlock[mb].indexOfFirstAtomInMoleculeBlock)
                                        % atomInfoForEachMoleculeBlock[mb].atomInfo.size()];
    }

    return atomInfo;
}

/* Sets the sum of charges (squared) and C6 in the system in fr.
 * Returns whether the system has a net charge.
 */
static bool set_chargesum(FILE* log, t_forcerec* fr, const gmx_mtop_t& mtop)
{
    /*This now calculates sum for q and c6*/
    double qsum, q2sum, q, c6sum, c6;

    qsum  = 0;
    q2sum = 0;
    c6sum = 0;
    for (const gmx_molblock_t& molb : mtop.molblock)
    {
        int            nmol  = molb.nmol;
        const t_atoms* atoms = &mtop.moltype[molb.type].atoms;
        for (int i = 0; i < atoms->nr; i++)
        {
            q = atoms->atom[i].q;
            qsum += nmol * q;
            q2sum += nmol * q * q;
            c6 = mtop.ffparams.iparams[atoms->atom[i].type * (mtop.ffparams.atnr + 1)].lj.c6;
            c6sum += nmol * c6;
        }
    }
    fr->qsum[0]  = qsum;
    fr->q2sum[0] = q2sum;
    fr->c6sum[0] = c6sum;

    if (fr->efep != FreeEnergyPerturbationType::No)
    {
        qsum  = 0;
        q2sum = 0;
        c6sum = 0;
        for (const gmx_molblock_t& molb : mtop.molblock)
        {
            int            nmol  = molb.nmol;
            const t_atoms* atoms = &mtop.moltype[molb.type].atoms;
            for (int i = 0; i < atoms->nr; i++)
            {
                q = atoms->atom[i].qB;
                qsum += nmol * q;
                q2sum += nmol * q * q;
                c6 = mtop.ffparams.iparams[atoms->atom[i].typeB * (mtop.ffparams.atnr + 1)].lj.c6;
                c6sum += nmol * c6;
            }
            fr->qsum[1]  = qsum;
            fr->q2sum[1] = q2sum;
            fr->c6sum[1] = c6sum;
        }
    }
    else
    {
        fr->qsum[1]  = fr->qsum[0];
        fr->q2sum[1] = fr->q2sum[0];
        fr->c6sum[1] = fr->c6sum[0];
    }
    if (log)
    {
        if (fr->efep == FreeEnergyPerturbationType::No)
        {
            fprintf(log, "System total charge: %.3f\n", fr->qsum[0]);
        }
        else
        {
            fprintf(log, "System total charge, top. A: %.3f top. B: %.3f\n", fr->qsum[0], fr->qsum[1]);
        }
    }

    /* A cut-off of 1e-4 is used to catch rounding errors due to ascii input */
    return (std::abs(fr->qsum[0]) > 1e-4 || std::abs(fr->qsum[1]) > 1e-4);
}

/*!\brief If there's bonded interactions of type \c ftype1 or \c
 * ftype2 present in the topology, build an array of the number of
 * interactions present for each bonded interaction index found in the
 * topology.
 *
 * \c ftype1 or \c ftype2 may be set to -1 to disable seeking for a
 * valid type with that parameter.
 *
 * \c count will be reallocated as necessary to fit the largest bonded
 * interaction index found, and its current size will be returned in
 * \c ncount. It will contain zero for every bonded interaction index
 * for which no interactions are present in the topology.
 */
static void count_tables(int ftype1, int ftype2, const gmx_mtop_t& mtop, int* ncount, int** count)
{
    int ftype, i, j, tabnr;

    // Loop over all moleculetypes
    for (const gmx_moltype_t& molt : mtop.moltype)
    {
        // Loop over all interaction types
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            // If the current interaction type is one of the types whose tables we're trying to count...
            if (ftype == ftype1 || ftype == ftype2)
            {
                const InteractionList& il     = molt.ilist[ftype];
                const int              stride = 1 + NRAL(ftype);
                // ... and there are actually some interactions for this type
                for (i = 0; i < il.size(); i += stride)
                {
                    // Find out which table index the user wanted
                    tabnr = mtop.ffparams.iparams[il.iatoms[i]].tab.table;
                    if (tabnr < 0)
                    {
                        gmx_fatal(FARGS, "A bonded table number is smaller than 0: %d\n", tabnr);
                    }
                    // Make room for this index in the data structure
                    if (tabnr >= *ncount)
                    {
                        srenew(*count, tabnr + 1);
                        for (j = *ncount; j < tabnr + 1; j++)
                        {
                            (*count)[j] = 0;
                        }
                        *ncount = tabnr + 1;
                    }
                    // Record that this table index is used and must have a valid file
                    (*count)[tabnr]++;
                }
            }
        }
    }
}

/*!\brief If there's bonded interactions of flavour \c tabext and type
 * \c ftype1 or \c ftype2 present in the topology, seek them in the
 * list of filenames passed to mdrun, and make bonded tables from
 * those files.
 *
 * \c ftype1 or \c ftype2 may be set to -1 to disable seeking for a
 * valid type with that parameter.
 *
 * A fatal error occurs if no matching filename is found.
 */
static std::vector<bondedtable_t> make_bonded_tables(FILE*                            fplog,
                                                     int                              ftype1,
                                                     int                              ftype2,
                                                     const gmx_mtop_t&                mtop,
                                                     gmx::ArrayRef<const std::string> tabbfnm,
                                                     const char*                      tabext)
{
    std::vector<bondedtable_t> tab;

    int  ncount = 0;
    int* count  = nullptr;
    count_tables(ftype1, ftype2, mtop, &ncount, &count);

    // Are there any relevant tabulated bond interactions?
    if (ncount > 0)
    {
        tab.resize(ncount);
        for (int i = 0; i < ncount; i++)
        {
            // Do any interactions exist that requires this table?
            if (count[i] > 0)
            {
                // This pattern enforces the current requirement that
                // table filenames end in a characteristic sequence
                // before the file type extension, and avoids table 13
                // being recognized and used for table 1.
                std::string patternToFind = gmx::formatString("_%s%d.%s", tabext, i, ftp2ext(efXVG));
                bool        madeTable     = false;
                for (gmx::Index j = 0; j < tabbfnm.ssize() && !madeTable; ++j)
                {
                    if (gmx::endsWith(tabbfnm[j], patternToFind))
                    {
                        // Finally read the table from the file found
                        tab[i]    = make_bonded_table(fplog, tabbfnm[j].c_str(), NRAL(ftype1) - 2);
                        madeTable = true;
                    }
                }
                if (!madeTable)
                {
                    bool isPlural = (ftype2 != -1);
                    gmx_fatal(FARGS,
                              "Tabulated interaction of type '%s%s%s' with index %d cannot be used "
                              "because no table file whose name matched '%s' was passed via the "
                              "gmx mdrun -tableb command-line option.",
                              interaction_function[ftype1].longname,
                              isPlural ? "' or '" : "",
                              isPlural ? interaction_function[ftype2].longname : "",
                              i,
                              patternToFind.c_str());
                }
            }
        }
        sfree(count);
    }

    return tab;
}

void forcerec_set_ranges(t_forcerec* fr, int natoms_force, int natoms_force_constr, int natoms_f_novirsum)
{
    fr->natoms_force        = natoms_force;
    fr->natoms_force_constr = natoms_force_constr;

    for (auto& forceHelperBuffers : fr->forceHelperBuffers)
    {
        forceHelperBuffers.resize(natoms_f_novirsum);
    }
}

/* Generate Coulomb and/or Van der Waals Ewald long-range correction tables
 *
 * Tables are generated for one or both, depending on if the pointers
 * are non-null. The spacing for both table sets is the same and obeys
 * both accuracy requirements, when relevant.
 */
static void init_ewald_f_table(const interaction_const_t& ic,
                               const real                 rlist,
                               const real                 tabext,
                               EwaldCorrectionTables*     coulombTables,
                               EwaldCorrectionTables*     vdwTables)
{
    const bool useCoulombTable = (usingPmeOrEwald(ic.eeltype) && coulombTables != nullptr);
    const bool useVdwTable     = (usingLJPme(ic.vdwtype) && vdwTables != nullptr);

    /* Get the Ewald table spacing based on Coulomb and/or LJ
     * Ewald coefficients and rtol.
     */
    const real tableScale = ewald_spline3_table_scale(ic, useCoulombTable, useVdwTable);

    const bool havePerturbedNonbondeds = (ic.softCoreParameters != nullptr);

    real tableLen = ic.rcoulomb;
    if ((useCoulombTable || useVdwTable) && havePerturbedNonbondeds && rlist + tabext > 0.0)
    {
        /* TODO: Ideally this should also check if couple-intramol == no, but that isn't
         * stored in ir. Grompp puts that info into an opts structure that doesn't make it into the tpr.
         * The alternative is to look through all the exclusions and check if they come from
         * couple-intramol == no. Meanwhile, always having larger tables should only affect
         * memory consumption, not speed (barring cache issues).
         */
        tableLen = rlist + tabext;
    }
    const int tableSize = static_cast<int>(tableLen * tableScale) + 2;

    if (useCoulombTable)
    {
        *coulombTables =
                generateEwaldCorrectionTables(tableSize, tableScale, ic.ewaldcoeff_q, v_q_ewald_lr);
    }

    if (useVdwTable)
    {
        *vdwTables = generateEwaldCorrectionTables(tableSize, tableScale, ic.ewaldcoeff_lj, v_lj_ewald_lr);
    }
}

void init_interaction_const_tables(FILE* fp, interaction_const_t* ic, const real rlist, const real tableExtensionLength)
{
    if (usingPmeOrEwald(ic->eeltype) || usingLJPme(ic->vdwtype))
    {
        init_ewald_f_table(
                *ic, rlist, tableExtensionLength, ic->coulombEwaldTables.get(), ic->vdwEwaldTables.get());
        if (fp != nullptr)
        {
            if (usingPmeOrEwald(ic->eeltype))
            {
                fprintf(fp,
                        "Initialized non-bonded Coulomb Ewald tables, spacing: %.2e size: %zu\n\n",
                        1 / ic->coulombEwaldTables->scale,
                        ic->coulombEwaldTables->tableF.size());
            }
        }
    }
}

real cutoff_inf(real cutoff)
{
    if (cutoff == 0)
    {
        cutoff = GMX_CUTOFF_INF;
    }

    return cutoff;
}

void init_forcerec(FILE*                            fplog,
                   const gmx::MDLogger&             mdlog,
                   const gmx::SimulationWorkload&   simulationWork,
                   t_forcerec*                      forcerec,
                   const t_inputrec&                inputrec,
                   const gmx_mtop_t&                mtop,
                   const t_commrec*                 commrec,
                   matrix                           box,
                   const char*                      tabfn,
                   const char*                      tabpfn,
                   gmx::ArrayRef<const std::string> tabbfnm,
                   real                             print_force)
{
    /* The CMake default turns SIMD kernels on, but it might be turned off further down... */
    forcerec->use_simd_kernels = GMX_USE_SIMD_KERNELS;

    if (check_box(inputrec.pbcType, box))
    {
        gmx_fatal(FARGS, "%s", check_box(inputrec.pbcType, box));
    }

    /* Test particle insertion ? */
    if (EI_TPI(inputrec.eI))
    {
        /* Set to the size of the molecule to be inserted (the last one) */
        gmx::RangePartitioning molecules = gmx_mtop_molecules(mtop);
        forcerec->n_tpi                  = molecules.block(molecules.numBlocks() - 1).size();
    }
    else
    {
        forcerec->n_tpi = 0;
    }

    if (inputrec.coulombtype == CoulombInteractionType::RFNecUnsupported
        || inputrec.coulombtype == CoulombInteractionType::GRFNotused)
    {
        gmx_fatal(FARGS, "%s electrostatics is no longer supported", enumValueToString(inputrec.coulombtype));
    }

    if (inputrec.bAdress)
    {
        gmx_fatal(FARGS, "AdResS simulations are no longer supported");
    }
    if (inputrec.useTwinRange)
    {
        gmx_fatal(FARGS, "Twin-range simulations are no longer supported");
    }
    /* Copy the user determined parameters */
    forcerec->userint1  = inputrec.userint1;
    forcerec->userint2  = inputrec.userint2;
    forcerec->userint3  = inputrec.userint3;
    forcerec->userint4  = inputrec.userint4;
    forcerec->userreal1 = inputrec.userreal1;
    forcerec->userreal2 = inputrec.userreal2;
    forcerec->userreal3 = inputrec.userreal3;
    forcerec->userreal4 = inputrec.userreal4;

    /* Shell stuff */
    forcerec->fc_stepsize = inputrec.fc_stepsize;

    /* Free energy */
    forcerec->efep = inputrec.efep;

    if ((getenv("GMX_DISABLE_SIMD_KERNELS") != nullptr) || (getenv("GMX_NOOPTIMIZEDKERNELS") != nullptr))
    {
        forcerec->use_simd_kernels = FALSE;
        if (fplog != nullptr)
        {
            fprintf(fplog,
                    "\nFound environment variable GMX_DISABLE_SIMD_KERNELS.\n"
                    "Disabling the usage of any SIMD-specific non-bonded & bonded kernel routines\n"
                    "(e.g. SSE2/SSE4.1/AVX).\n\n");
        }
    }

    forcerec->haveBuckingham = (mtop.ffparams.functype[0] == F_BHAM);

    /* Neighbour searching stuff */
    forcerec->pbcType = inputrec.pbcType;

    /* Determine if we will do PBC for distances in bonded interactions */
    if (forcerec->pbcType == PbcType::No)
    {
        forcerec->bMolPBC = FALSE;
    }
    else
    {
        forcerec->bMolPBC =
                (!haveDDAtomOrdering(*commrec) || dd_bonded_molpbc(*commrec->dd, forcerec->pbcType));

        // Check and set up PBC for Ewald surface corrections or orientation restraints
        const bool useEwaldSurfaceCorrection =
                (usingPmeOrEwald(inputrec.coulombtype) && inputrec.epsilon_surface != 0);
        const bool haveOrientationRestraints = (gmx_mtop_ftype_count(mtop, F_ORIRES) > 0);
        const bool moleculesAreAlwaysWhole =
                (haveDDAtomOrdering(*commrec) && dd_moleculesAreAlwaysWhole(*commrec->dd));
        // WholeMoleculeTransform is only supported with a single PP rank
        if (!moleculesAreAlwaysWhole && !havePPDomainDecomposition(commrec)
            && (useEwaldSurfaceCorrection || haveOrientationRestraints))
        {
            if (havePPDomainDecomposition(commrec))
            {
                gmx_fatal(FARGS,
                          "You requested Ewald surface correction or orientation restraints, "
                          "but molecules are broken "
                          "over periodic boundary conditions by the domain decomposition. "
                          "Run without domain decomposition instead.");
            }

            forcerec->wholeMoleculeTransform = std::make_unique<gmx::WholeMoleculeTransform>(
                    mtop, inputrec.pbcType, haveDDAtomOrdering(*commrec));
        }

        forcerec->bMolPBC =
                !haveDDAtomOrdering(*commrec) || dd_bonded_molpbc(*commrec->dd, forcerec->pbcType);

        if (useEwaldSurfaceCorrection)
        {
            GMX_RELEASE_ASSERT(moleculesAreAlwaysWhole || forcerec->wholeMoleculeTransform,
                               "Molecules can not be broken by PBC with epsilon_surface > 0");
        }
    }

    forcerec->rc_scaling = inputrec.pressureCouplingOptions.refcoord_scaling;
    copy_rvec(inputrec.posres_com, forcerec->posres_com);
    copy_rvec(inputrec.posres_comB, forcerec->posres_comB);

    forcerec->haveBoxDeformation = ir_haveBoxDeformation(inputrec);

    forcerec->rlist                  = cutoff_inf(inputrec.rlist);
    forcerec->ljpme_combination_rule = inputrec.ljpme_combination_rule;

    /* This now calculates sum for q and c6*/
    bool systemHasNetCharge = set_chargesum(fplog, forcerec, mtop);

    /* Make data structure used by kernels */
    forcerec->ic = std::make_unique<interaction_const_t>(
            init_interaction_const(fplog, inputrec, mtop, systemHasNetCharge));
    init_interaction_const_tables(fplog, forcerec->ic.get(), forcerec->rlist, inputrec.tabext);

    const interaction_const_t* interactionConst = forcerec->ic.get();

    /* Electrostatics: Translate from interaction-setting-in-mdp-file to kernel interaction format */
    switch (interactionConst->eeltype)
    {
        case CoulombInteractionType::Cut:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::Coulomb;
            break;

        case CoulombInteractionType::RF:
        case CoulombInteractionType::RFZero:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::ReactionField;
            break;

        case CoulombInteractionType::Switch:
        case CoulombInteractionType::Shift:
        case CoulombInteractionType::User:
        case CoulombInteractionType::PmeSwitch:
        case CoulombInteractionType::PmeUser:
        case CoulombInteractionType::PmeUserSwitch:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::CubicSplineTable;
            break;

        case CoulombInteractionType::Pme:
        case CoulombInteractionType::P3mAD:
        case CoulombInteractionType::Ewald:
            forcerec->nbkernel_elec_interaction = NbkernelElecType::Ewald;
            break;

        default:
            gmx_fatal(FARGS,
                      "Unsupported electrostatic interaction: %s",
                      enumValueToString(interactionConst->eeltype));
    }
    forcerec->nbkernel_elec_modifier = interactionConst->coulomb_modifier;

    /* Vdw: Translate from mdp settings to kernel format */
    switch (interactionConst->vdwtype)
    {
        case VanDerWaalsType::Cut:
            if (forcerec->haveBuckingham)
            {
                forcerec->nbkernel_vdw_interaction = NbkernelVdwType::Buckingham;
            }
            else
            {
                forcerec->nbkernel_vdw_interaction = NbkernelVdwType::LennardJones;
            }
            break;
        case VanDerWaalsType::Pme:
            forcerec->nbkernel_vdw_interaction = NbkernelVdwType::LJEwald;
            break;

        case VanDerWaalsType::Switch:
        case VanDerWaalsType::Shift:
        case VanDerWaalsType::User:
            forcerec->nbkernel_vdw_interaction = NbkernelVdwType::CubicSplineTable;
            break;

        default:
            gmx_fatal(FARGS, "Unsupported vdw interaction: %s", enumValueToString(interactionConst->vdwtype));
    }
    forcerec->nbkernel_vdw_modifier = interactionConst->vdw_modifier;

    if (!gmx_within_tol(interactionConst->reppow, 12.0, 10 * GMX_DOUBLE_EPS))
    {
        gmx_fatal(FARGS, "Only LJ repulsion power 12 is supported");
    }
    /* Older tpr files can contain Coulomb user tables with the Verlet cutoff-scheme,
     * while mdrun does not (and never did) support this.
     */
    if (usingUserTableElectrostatics(forcerec->ic->eeltype))
    {
        gmx_fatal(FARGS,
                  "Electrostatics type %s is currently not supported",
                  enumValueToString(inputrec.coulombtype));
    }

    /* 1-4 interaction electrostatics */
    forcerec->fudgeQQ = mtop.ffparams.fudgeQQ;

    if (simulationWork.useMts)
    {
        GMX_RELEASE_ASSERT(gmx::checkMtsRequirements(inputrec).empty(),
                           "All MTS requirements should be met here");
    }

    const bool haveDirectVirialContributionsFast =
            forcerec->forceProviders->hasForceProvider() || gmx_mtop_ftype_count(mtop, F_POSRES) > 0
            || gmx_mtop_ftype_count(mtop, F_FBPOSRES) > 0 || inputrec.nwall > 0 || inputrec.bPull
            || inputrec.bRot || inputrec.bIMD;
    const bool haveDirectVirialContributionsSlow = usingFullElectrostatics(interactionConst->eeltype)
                                                   || usingLJPme(interactionConst->vdwtype);
    for (int i = 0; i < (simulationWork.useMts ? 2 : 1); i++)
    {
        bool haveDirectVirialContributions =
                (((!simulationWork.useMts || i == 0) && haveDirectVirialContributionsFast)
                 || ((!simulationWork.useMts || i == 1) && haveDirectVirialContributionsSlow));
        forcerec->forceHelperBuffers.emplace_back(haveDirectVirialContributions);
    }

    if (forcerec->shift_vec.empty())
    {
        forcerec->shift_vec.resize(gmx::c_numShiftVectors);
    }

    GMX_ASSERT(forcerec->nbfp.empty(), "The nonbonded force parameters should not be set up yet.");
    forcerec->ntype = mtop.ffparams.atnr;
    forcerec->nbfp  = makeNonBondedParameterLists(
            mtop.ffparams.atnr, mtop.ffparams.iparams, forcerec->haveBuckingham);
    if (usingLJPme(interactionConst->vdwtype))
    {
        forcerec->ljpme_c6grid = makeLJPmeC6GridCorrectionParameters(
                mtop.ffparams.atnr, mtop.ffparams.iparams, forcerec->ljpme_combination_rule);
    }

    /* Copy the energy group exclusions */
    forcerec->egp_flags = inputrec.opts.egp_flags;

    /* Van der Waals stuff */
    if ((interactionConst->vdwtype != VanDerWaalsType::Cut)
        && (interactionConst->vdwtype != VanDerWaalsType::User) && !forcerec->haveBuckingham)
    {
        if (interactionConst->rvdw_switch >= interactionConst->rvdw)
        {
            gmx_fatal(FARGS,
                      "rvdw_switch (%f) must be < rvdw (%f)",
                      interactionConst->rvdw_switch,
                      interactionConst->rvdw);
        }
        if (fplog)
        {
            fprintf(fplog,
                    "Using %s Lennard-Jones, switch between %g and %g nm\n",
                    (interactionConst->eeltype == CoulombInteractionType::Switch) ? "switched" : "shifted",
                    interactionConst->rvdw_switch,
                    interactionConst->rvdw);
        }
    }

    if (forcerec->haveBuckingham && usingLJPme(interactionConst->vdwtype))
    {
        gmx_fatal(FARGS, "LJ PME not supported with Buckingham");
    }

    if (forcerec->haveBuckingham
        && (interactionConst->vdwtype == VanDerWaalsType::Shift
            || interactionConst->vdwtype == VanDerWaalsType::Switch))
    {
        gmx_fatal(FARGS, "Switch/shift interaction not supported with Buckingham");
    }

    if (forcerec->haveBuckingham)
    {
        gmx_fatal(FARGS, "The Verlet cutoff-scheme does not (yet) support Buckingham");
    }

    if (inputrec.implicit_solvent)
    {
        gmx_fatal(FARGS, "Implict solvation is no longer supported.");
    }


    /* This code automatically gives table length tabext without cut-off's,
     * in that case grompp should already have checked that we do not need
     * normal tables and we only generate tables for 1-4 interactions.
     */
    real rtab = inputrec.rlist + inputrec.tabext;

    /* We want to use unmodified tables for 1-4 coulombic
     * interactions, so we must in general have an extra set of
     * tables. */
    if (gmx_mtop_ftype_count(mtop, F_LJ14) > 0 || gmx_mtop_ftype_count(mtop, F_LJC14_Q) > 0
        || gmx_mtop_ftype_count(mtop, F_LJC_PAIRS_NB) > 0)
    {
        forcerec->pairsTable = make_tables(fplog, interactionConst, tabpfn, rtab, GMX_MAKETABLES_14ONLY);
    }

    /* Wall stuff */
    forcerec->nwall = inputrec.nwall;
    if (inputrec.nwall && inputrec.wall_type == WallType::Table)
    {
        make_wall_tables(fplog, inputrec, tabfn, &mtop.groups, forcerec);
    }

    forcerec->fcdata = std::make_unique<t_fcdata>();

    if (!tabbfnm.empty())
    {
        t_fcdata& fcdata = *forcerec->fcdata;
        // Need to catch std::bad_alloc
        // TODO Don't need to catch this here, when merging with main branch
        try
        {
            // TODO move these tables into a separate struct and store reference in ListedForces
            fcdata.bondtab = make_bonded_tables(fplog, F_TABBONDS, F_TABBONDSNC, mtop, tabbfnm, "b");
            fcdata.angletab = make_bonded_tables(fplog, F_TABANGLES, -1, mtop, tabbfnm, "a");
            fcdata.dihtab   = make_bonded_tables(fplog, F_TABDIHS, -1, mtop, tabbfnm, "d");
        }
        GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    }
    else
    {
        if (debug)
        {
            fprintf(debug,
                    "No fcdata or table file name passed, can not read table, can not do bonded "
                    "interactions\n");
        }
    }

    /* Initialize the thread working data for bonded interactions */
    if (simulationWork.useMts)
    {
        // Add one ListedForces object for each MTS level
        bool isFirstLevel = true;
        for (const auto& mtsLevel : inputrec.mtsLevels)
        {
            ListedForces::InteractionSelection interactionSelection;
            const auto&                        forceGroups = mtsLevel.forceGroups;
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Pair)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Pairs));
            }
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Dihedral)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Dihedrals));
            }
            if (forceGroups[static_cast<int>(gmx::MtsForceGroups::Angle)])
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Angles));
            }
            if (isFirstLevel)
            {
                interactionSelection.set(static_cast<int>(ListedForces::InteractionGroup::Rest));
                isFirstLevel = false;
            }
            forcerec->listedForces.emplace_back(
                    mtop.ffparams,
                    mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                    gmx_omp_nthreads_get(ModuleMultiThread::Bonded),
                    interactionSelection,
                    fplog);
        }
    }
    else
    {
        // Add one ListedForces object with all listed interactions
        forcerec->listedForces.emplace_back(
                mtop.ffparams,
                mtop.groups.groups[SimulationAtomGroupType::EnergyOutput].size(),
                gmx_omp_nthreads_get(ModuleMultiThread::Bonded),
                ListedForces::interactionSelectionAll(),
                fplog);
    }

    // QM/MM initialization if requested
    if (inputrec.bQMMM)
    {
        gmx_incons("QM/MM was requested, but is no longer available in GROMACS");
    }

    /* Set all the static charge group info */
    forcerec->atomInfoForEachMoleculeBlock = makeAtomInfoForEachMoleculeBlock(mtop, forcerec);
    if (!haveDDAtomOrdering(*commrec))
    {
        forcerec->atomInfo = expandAtomInfo(mtop.molblock.size(), forcerec->atomInfoForEachMoleculeBlock);
    }

    if (!haveDDAtomOrdering(*commrec))
    {
        forcerec_set_ranges(forcerec, mtop.natoms, mtop.natoms, mtop.natoms);
    }

    forcerec->print_force = print_force;

    if (inputrec.eDispCorr != DispersionCorrectionType::No)
    {
        forcerec->dispersionCorrection = std::make_unique<DispersionCorrection>(
                mtop, inputrec, forcerec->haveBuckingham, *forcerec->ic, tabfn);
        forcerec->dispersionCorrection->print(mdlog);
    }

    if (fplog != nullptr)
    {
        /* Here we switch from using mdlog, which prints the newline before
         * the paragraph, to our old fprintf logging, which prints the newline
         * after the paragraph, so we should add a newline here.
         */
        fprintf(fplog, "\n");
    }
}

t_forcerec::t_forcerec() = default;

t_forcerec::~t_forcerec() = default;
