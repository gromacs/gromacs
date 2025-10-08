/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Implements embedded system topology preprocessing functions
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_topology
 */

#include "gromacs/topology/embedded_system_preprocessing.h"

#include "gromacs/fileio/warninp.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

namespace gmx
{

static bool isEmbeddedAtom(Index globalAtomIndex, const std::set<int>& embeddedIndices)
{
    return (embeddedIndices.find(globalAtomIndex) != embeddedIndices.end());
}

std::vector<bool> splitEmbeddedBlocks(gmx_mtop_t* mtop, const std::set<int>& embeddedIndices)
{
    // Global counter of atoms
    Index             iAt = 0;
    std::vector<bool> isEmbeddedBlock;

    /* Counter of molecules point to the specific moltype
     * i.e molblock 0 has 2 molecules have moltype 0 and molblock 2 has 1 additional molecule of type 0
     * then will be numMoleculesOfType[0] = 2 + 1 = 3
     * That counter is needed to decide whether we should create a new moltype for the molecule contatining embedded atoms
     * or we could modify existing moltype if there is only one molecule of that type
     */
    std::vector<int> numMoleculesOfType(mtop->moltype.size());
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        numMoleculesOfType[mtop->molblock[molBlockIndex].type] += mtop->molblock[molBlockIndex].nmol;
    }

    // Loop over all blocks in topology
    // molBlockIndex - current index of block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // Initialize block as non-embedded first
        isEmbeddedBlock.push_back(false);

        // Pointer to current block
        gmx_molblock_t* molBlock = &mtop->molblock[molBlockIndex];

        // Number of atoms in all molecules of current block
        const int numAtomsInMolecule = mtop->moltype[molBlock->type].atoms.nr;

        // Loop over all molecules in molBlock
        // mol - index of current molecule in molBlock
        for (int mol = 0; mol < molBlock->nmol; mol++)
        {

            // search for embedded atoms in current molecule
            bool isEmbedded = false;
            for (int i = 0; i < numAtomsInMolecule; i++)
            {
                if (isEmbeddedAtom(iAt, embeddedIndices))
                {
                    isEmbedded = true;
                }
                iAt++;
            }

            // Apparently current molecule (molBlockIndex, mol) contains embedded atoms
            // We should split it from the current block and create new blocks
            // For that molecule and all molecules after it new block will be created
            if (isEmbedded)
            {
                // if this block contains only 1 molecule, then splitting not needed
                if (molBlock->nmol > 1)
                {
                    // If current molecule is not the first in the block,
                    // then we need to first create block for it and all molecule before it
                    if (mol > 0)
                    {
                        // Split the molblock at this molecule
                        auto pos = mtop->molblock.begin() + molBlockIndex + 1;
                        mtop->molblock.insert(pos, mtop->molblock[molBlockIndex]);
                        mtop->molblock[molBlockIndex].nmol = mol;
                        mtop->molblock[molBlockIndex + 1].nmol -= mol;
                        isEmbeddedBlock[molBlockIndex] = false;
                        isEmbeddedBlock.push_back(true);
                        molBlockIndex++;
                        molBlock = &mtop->molblock[molBlockIndex];
                    }

                    // If current molecule is not the only one in new block,
                    // Then we split new block after that molecule
                    if (molBlock->nmol > 1)
                    {
                        auto pos = mtop->molblock.begin() + molBlockIndex + 1;
                        mtop->molblock.insert(pos, mtop->molblock[molBlockIndex]);
                        molBlock                           = &mtop->molblock[molBlockIndex];
                        mtop->molblock[molBlockIndex].nmol = 1;
                        mtop->molblock[molBlockIndex + 1].nmol -= 1;
                        isEmbeddedBlock[molBlockIndex] = true;
                    }
                }
                else
                {
                    isEmbeddedBlock[molBlockIndex] = true;
                }

                // Create a copy of a moltype for a molecule
                // containing embedded atoms and append it in the end of the moltype vector
                // if there is only 1 molecule pointing to that type then skip that step
                if (numMoleculesOfType[molBlock->type] > 1)
                {
                    // Here comes a huge piece of "not so good" code, because of deleted operator= from gmx_moltype_t
                    std::vector<gmx_moltype_t> temp(mtop->moltype.size());
                    for (size_t i = 0; i < mtop->moltype.size(); ++i)
                    {
                        copy_moltype(&mtop->moltype[i], &temp[i]);
                    }
                    mtop->moltype.resize(temp.size() + 1);
                    for (size_t i = 0; i < temp.size(); ++i)
                    {
                        copy_moltype(&temp[i], &mtop->moltype[i]);
                    }
                    copy_moltype(&mtop->moltype[molBlock->type], &mtop->moltype.back());

                    // Set the molecule type for the embedded molblock
                    molBlock->type = mtop->moltype.size() - 1;
                }
            }
        }
    }
    // Call finalize() to rebuild Block Indicies or else atoms lookup will fail
    mtop->finalize();
    return isEmbeddedBlock;
}

std::vector<real> removeEmbeddedClassicalCharges(gmx_mtop_t*              mtop,
                                                 const std::set<int>&     embeddedIndices,
                                                 const std::vector<bool>& isEmbeddedBlock,
                                                 real                     refQ,
                                                 const MDLogger&          logger,
                                                 WarningHandler*          wi)
{
    // Loop over all atoms and remove charge if they are embedded atoms.
    // Sum-up total removed charge and remaning charge on MM atoms
    // Build atomCharges_ vector
    int               molBlockIndex               = 0;
    real              totalClassicalChargeRemoved = 0.0;
    real              remainingMMCharge           = 0.0;
    std::vector<real> atomCharges;
    for (int i = 0; i < mtop->natoms; i++)
    {
        int indexInMolecule;
        mtopGetMolblockIndex(*mtop, i, &molBlockIndex, nullptr, &indexInMolecule);
        t_atom* atom = &mtop->moltype[mtop->molblock[molBlockIndex].type].atoms.atom[indexInMolecule];
        if (isEmbeddedAtom(i, embeddedIndices))
        {
            totalClassicalChargeRemoved += atom->q;
            atom->q  = 0.0;
            atom->qB = 0.0;
        }
        else
        {
            remainingMMCharge += atom->q;
        }

        atomCharges.push_back(atom->q);
    }

    // also remove charges on virtual sites
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    int numVirtualSitesModified = 0;
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If this is VSite interaction and ilist is not empty
                if (IS_VSITE(ftype) && !molType->ilist[ftype].empty())
                {
                    // Number of elements in the iatoms array for the current ftype
                    const int numInteractionElements = NRAL(ftype) + 1;

                    // Loop over all interactions of ftype
                    for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                    {
                        // Calculate number of embedded atoms in the interaction
                        int numEmbedded = 0;
                        // Here k starts from 2 because first atom in the interaction is an actuall vsite index
                        for (int k = 2; k <= NRAL(ftype); k++)
                        {
                            if (isEmbeddedAtom(molType->ilist[ftype].iatoms[j + k] + start, embeddedIndices))
                            {
                                numEmbedded++;
                            }
                        }

                        // If all atoms froming that virtual site are embedded atoms
                        // then remove classical charge from that virtual site
                        if (numEmbedded == (NRAL(ftype) - 1))
                        {
                            numVirtualSitesModified++;
                            totalClassicalChargeRemoved +=
                                    molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].q;
                            molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].q  = 0.0;
                            molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].qB = 0.0;
                        }
                    }
                }
            }
        }
    }

    // avoid negative zero charge
    real totCharge = remainingMMCharge + totalClassicalChargeRemoved;
    if (std::abs(totCharge) < 1E-5)
    {
        totCharge = 0.0;
    }
    // log the modifications
    GMX_LOG(logger.info)
            .appendTextFormatted(
                    "Total charge of the classical system (before modifications): %.5f\n", totCharge);
    GMX_LOG(logger.info).appendTextFormatted("Classical charge removed from embedded atoms: %.5f\n", totalClassicalChargeRemoved);
    if (numVirtualSitesModified > 0)
    {
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "Note: There are %d virtual sites found, which are built from embedded "
                        "atoms only. "
                        "Classical charges on them have been removed as well.\n",
                        numVirtualSitesModified);
    }

    // we should warn the user if there is inconsistence between removed classical charges
    GMX_ASSERT(wi, "WarningHandler not set.");
    if (std::abs(totalClassicalChargeRemoved - refQ) > 1E-5)
    {
        std::string msg = formatString(
                "Total charge of your embedded system differs from classical system! "
                "Consider manually spreading %.5lf charge over MM atoms near to the embedded "
                "region\n",
                totalClassicalChargeRemoved - refQ);
        wi->addWarning(msg);
    }

    return atomCharges;
}

void addEmbeddedNBExclusions(gmx_mtop_t* mtop, const std::set<int>& embeddedIndices, const MDLogger& logger)
{
    // Add all embedded atoms to the mtop->intermolecularExclusionGroup
    mtop->intermolecularExclusionGroup.reserve(mtop->intermolecularExclusionGroup.size()
                                               + embeddedIndices.size());
    int numExclusionsMade = 0;
    for (auto i : embeddedIndices)
    {
        mtop->intermolecularExclusionGroup.push_back(i);
        numExclusionsMade++;
    }
    GMX_LOG(logger.info).appendTextFormatted("Number of exclusions made: %d\n", numExclusionsMade);
}

std::vector<int> buildEmbeddedAtomNumbers(const gmx_mtop_t& mtop)
{
    // Save to atomNumbers_ atom numbers of all atoms
    std::vector<int> atomNumbers;
    AtomIterator     atoms(mtop);
    while (atoms->globalAtomNumber() < mtop.natoms)
    {
        // Check if we have valid atomnumbers
        if (atoms->atom().atomnumber < 0)
        {
            gmx_fatal(FARGS,
                      "Atoms %d does not have atomic number needed for embedded subsystem. "
                      "Check atomtypes section in your topology or forcefield.",
                      atoms->globalAtomNumber());
        }

        atomNumbers.push_back(atoms->atom().atomnumber);
        atoms++;
    }

    return atomNumbers;
}

void modifyEmbeddedTwoCenterInteractions(gmx_mtop_t*              mtop,
                                         const std::set<int>&     embeddedIndices,
                                         const std::vector<bool>& isEmbeddedBlock,
                                         const MDLogger&          logger)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    int numBondsRemoved   = 0;
    int numConnBondsAdded = 0;
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index if the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If not bonded interaction or InteractionFunction::ConnectBonds, or some form of Restraints,
                // or not pair interaction, or no interactions of that type: then go the next type
                if (!(interaction_function[ftype].flags & IF_BOND) || ftype == InteractionFunction::ConnectBonds
                    || ftype == InteractionFunction::RestraintBonds || ftype == InteractionFunction::HarmonicPotential
                    || ftype == InteractionFunction::DistanceRestraints
                    || ftype == InteractionFunction::OrientationRestraints
                    || ftype == InteractionFunction::AngleZAxisRestraints || NRAL(ftype) != 2
                    || molType->ilist[ftype].empty())
                {
                    continue;
                }

                // Number of elements in the iatoms array for the current ftype
                const int numInteractionElements = NRAL(ftype) + 1;

                // Buffer for preserving interactions that are retained
                AtomGroupIndices iatomsBuf;

                // Loop over all interactions of ftype
                for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                {

                    // If both atoms are embedded and it is IF_CHEMBOND then convert it to F_CONNBONDS
                    if (isEmbeddedAtom(molType->ilist[ftype].iatoms[j + 1] + start, embeddedIndices)
                        && isEmbeddedAtom(molType->ilist[ftype].iatoms[j + 2] + start, embeddedIndices))
                    {

                        // Add chemical bond to the InteractionFunction::ConnectBonds (bond type 5)
                        if (IS_CHEMBOND(ftype))
                        {
                            // Bond type is not used in InteractionFunction::ConnectBonds, so for generated bonds we set it to -1
                            const int connBondsType = -1;

                            // Add new InteractionFunction::ConnectBonds between atoms
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(connBondsType);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 2]);

                            numConnBondsAdded++;
                        }

                        // Since the current bond is not copied into iatomsBuf it will be removed from the topology
                        numBondsRemoved++;
                    }
                    else
                    {
                        // If one of atoms is not embedded then copy interaction into iatomsBuf
                        for (int k = 0; k < numInteractionElements; k++)
                        {
                            iatomsBuf.push_back(molType->ilist[ftype].iatoms[k + j]);
                        }
                    }
                }

                // Swap iatomsBuf and molType->ilist[ftype].iatoms
                molType->ilist[ftype].iatoms.swap(iatomsBuf);
            }
        }
    }
    if (numBondsRemoved > 0)
    {
        GMX_LOG(logger.info).appendTextFormatted("Number of bonds removed: %d\n", numBondsRemoved);
    }
    if (numConnBondsAdded > 0)
    {
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "Number of InteractionFunction::ConnectBonds (type 5 bonds) added: %d\n",
                        numConnBondsAdded);
    }
}

std::vector<LinkFrontier> buildLinkFrontier(gmx_mtop_t*              mtop,
                                            const std::set<int>&     embeddedIndices,
                                            const std::vector<bool>& isEmbeddedBlock,
                                            const MDLogger&          logger)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    std::vector<LinkFrontier> linkFrontier;
    int                       numLinkBonds = 0;
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - gloabl index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If not chemical bond interaction or not pair interaction
                // or no interactions of that type: then skip current ftype
                if (!(interaction_function[ftype].flags & IF_CHEMBOND) || NRAL(ftype) != 2
                    || molType->ilist[ftype].empty())
                {
                    continue;
                }

                // Number of elements in the iatoms array for the current ftype
                const int numInteractionElements = NRAL(ftype) + 1;

                // Loop over all interactions of ftype
                for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                {
                    // Global indexes of atoms involved into the interaction
                    int a1 = molType->ilist[ftype].iatoms[j + 1] + start;
                    int a2 = molType->ilist[ftype].iatoms[j + 2] + start;

                    // Update Link Frontier List if one of the atoms embedded and one MM
                    if (isEmbeddedAtom(a1, embeddedIndices) && !isEmbeddedAtom(a2, embeddedIndices))
                    {
                        linkFrontier.push_back({ a1, a2 });
                        numLinkBonds++;
                    }
                    if (isEmbeddedAtom(a2, embeddedIndices) && !isEmbeddedAtom(a1, embeddedIndices))
                    {
                        linkFrontier.push_back({ a2, a1 });
                        numLinkBonds++;
                    }
                }
            }
        }
    }
    if (numLinkBonds > 0)
    {
        GMX_LOG(logger.info).appendTextFormatted("Number of link bonds added: %d\n", numLinkBonds);
    }

    return linkFrontier;
}

void modifyEmbeddedThreeCenterInteractions(gmx_mtop_t*              mtop,
                                           const std::set<int>&     embeddedIndices,
                                           const std::vector<bool>& isEmbeddedBlock,
                                           const MDLogger&          logger)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    int numConnBondsAdded = 0;
    int numAnglesRemoved  = 0;
    int numSettleRemoved  = 0;
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If not bonded interaction or Restraints
                // or not three-particle interaction or no interactions of that type
                // and not Settle: then go the next type
                if ((!(interaction_function[ftype].flags & IF_BOND) || ftype == InteractionFunction::RestrictedBendingPotential
                     || NRAL(ftype) != 3 || molType->ilist[ftype].empty())
                    && (ftype != InteractionFunction::SETTLE))
                {
                    continue;
                }

                // Number of elements in the iatoms array for the current ftype
                const int numInteractionElements = NRAL(ftype) + 1;

                // Buffer for preserving interactions that are retained
                AtomGroupIndices iatomsBuf;

                // Loop over all interactions of ftype
                for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                {
                    // Calculate number of embedded atoms in the interaction
                    int numEmbedded = 0;
                    for (int k = 1; k <= NRAL(ftype); k++)
                    {
                        if (isEmbeddedAtom(molType->ilist[ftype].iatoms[j + k] + start, embeddedIndices))
                        {
                            numEmbedded++;
                        }
                    }

                    // If at least 2 atoms are embedded then remove interaction
                    if (numEmbedded >= 2)
                    {
                        // If this is SETTLE then replace it with two InteractionFunction::ConnectBonds
                        if (ftype == InteractionFunction::SETTLE)
                        {
                            // Bond type is not used in InteractionFunction::ConnectBonds, so for generated bonds we set it to -1
                            const int connBondsType = -1;

                            // Add InteractionFunction::ConnectBonds between atoms 1 and 2 first
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(connBondsType);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 2]);
                            // Then between atoms 1 and 3
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(connBondsType);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[InteractionFunction::ConnectBonds].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 3]);

                            numConnBondsAdded += 2;
                            numSettleRemoved++;
                        }
                        else
                        {
                            // If it is normal angle then it should be just removed
                            numAnglesRemoved++;
                        }
                    }
                    else
                    {
                        // If several (>1) atoms is not embedded then preserve interaction in the iatomsBuf
                        for (int k = 0; k < numInteractionElements; k++)
                        {
                            iatomsBuf.push_back(molType->ilist[ftype].iatoms[k + j]);
                        }
                    }
                }

                // Swap iatomsBuf and molType->ilist[ftype].iatoms
                molType->ilist[ftype].iatoms.swap(iatomsBuf);
            }
        }
    }
    if (numAnglesRemoved > 0)
    {
        GMX_LOG(logger.info).appendTextFormatted("Number of angles removed: %d\n", numAnglesRemoved);
    }
    if (numSettleRemoved > 0)
    {
        GMX_LOG(logger.info)
                .appendTextFormatted(
                        "Number of settles removed: %d (replaced by %d "
                        "InteractionFunction::ConnectBonds) \n",
                        numSettleRemoved,
                        numConnBondsAdded);
    }
}


void modifyEmbeddedFourCenterInteractions(gmx_mtop_t*              mtop,
                                          const std::set<int>&     embeddedIndices,
                                          const std::vector<bool>& isEmbeddedBlock,
                                          const MDLogger&          logger)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    int numDihedralsRemoved = 0;
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If not bonded interaction or Restraints
                // or not four-particle interaction or no interactions of that type: then go the next type
                if (!(interaction_function[ftype].flags & IF_BOND)
                    || ftype == InteractionFunction::RestrictedTorsionPotential || NRAL(ftype) != 4
                    || molType->ilist[ftype].empty())
                {
                    continue;
                }

                // Number of elements in the iatoms array for the current ftype
                const int numInteractionElements = NRAL(ftype) + 1;

                // Buffer for preserving interactions that are retained
                AtomGroupIndices iatomsBuf;

                // Loop over all interactions of ftype
                for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                {
                    // Calculate number of embedded atoms in the interaction
                    int numEmbedded = 0;
                    for (int k = 1; k <= NRAL(ftype); k++)
                    {
                        if (isEmbeddedAtom(molType->ilist[ftype].iatoms[j + k] + start, embeddedIndices))
                        {
                            numEmbedded++;
                        }
                    }

                    // If at least 3 atoms are embedded then remove interaction
                    if (numEmbedded >= 3)
                    {
                        numDihedralsRemoved++;
                    }
                    else
                    {
                        // If several (>1) atoms is not embedded then preserve interaction in the iatomsBuf
                        for (int k = 0; k < numInteractionElements; k++)
                        {
                            iatomsBuf.push_back(molType->ilist[ftype].iatoms[k + j]);
                        }
                    }
                }

                // Swap iatomsBuf and molType->ilist[ftype].iatoms
                molType->ilist[ftype].iatoms.swap(iatomsBuf);
            }
        }
    }
    if (numDihedralsRemoved > 0)
    {
        GMX_LOG(logger.info).appendTextFormatted("Number of dihedrals removed: %d\n", numDihedralsRemoved);
    }
}

void checkConstrainedBonds(gmx_mtop_t*              mtop,
                           const std::set<int>&     embeddedIndices,
                           const std::vector<bool>& isEmbeddedBlock,
                           WarningHandler*          wi)
{
    int numConstrainedBonds = 0;
    // Loop over all blocks in topology
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains embedded atoms
        if (isEmbeddedBlock[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - gloabl index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (const auto ftype : gmx::EnumerationWrapper<InteractionFunction>{})
            {
                // If not chemical bond interaction or not pair interaction
                // or no interactions of that type: then skip current ftype
                if (!(interaction_function[ftype].flags & IF_CHEMBOND) || NRAL(ftype) != 2
                    || molType->ilist[ftype].empty())
                {
                    continue;
                }

                // Number of elements in the iatoms array for the current ftype
                const int numInteractionElements = NRAL(ftype) + 1;

                // Loop over all interactions of ftype
                for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                {
                    // Global indexes of atoms involved into the interaction
                    int a1 = molType->ilist[ftype].iatoms[j + 1] + start;
                    int a2 = molType->ilist[ftype].iatoms[j + 2] + start;

                    // Check if it is constrained bond within embedded subsystem
                    if (isEmbeddedAtom(a2, embeddedIndices) && isEmbeddedAtom(a1, embeddedIndices)
                        && (interaction_function[ftype].flags & IF_CONSTRAINT))
                    {
                        numConstrainedBonds++;
                    }
                }
            }
        }
    }
    GMX_ASSERT(wi, "WarningHandler not set.");
    if (numConstrainedBonds > 2)
    {
        std::string msg =
                "Your embedded subsystem has a lot of constrained bonds. They probably have been "
                "generated automatically. That could produce artifacts in the simulation. "
                "Consider constraints = none in the mdp file.\n";
        wi->addWarning(msg);
    }
}

} // namespace gmx
