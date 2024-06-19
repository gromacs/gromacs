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
 * \brief
 * Implements QMMMTopologyPreprocessor
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "qmmmtopologypreprocessor.h"

#include <cstddef>

#include <filesystem>

#include "gromacs/applied_forces/qmmm/qmmmtypes.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_atomloops.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"

namespace gmx
{

QMMMTopologyPreprocessor::QMMMTopologyPreprocessor(ArrayRef<const Index> qmIndices) :
    qmIndices_(qmIndices.begin(), qmIndices.end())
{
}

void QMMMTopologyPreprocessor::preprocess(gmx_mtop_t* mtop)
{
    // 1) Split QM-containing molecules from other molecules in blocks
    splitQMblocks(mtop);

    // 2) Nullify charges on all virtual sites consisting of QM only atoms
    modifyQMMMVirtualSites(mtop);

    // 3) Nullify charges on all QM atoms
    removeQMClassicalCharges(mtop);

    // 4) Exclude LJ interactions between QM atoms
    addQMLJExclusions(mtop);

    // 5) Build atomNumbers vector with atomic numbers of all atoms
    buildQMMMAtomNumbers(mtop);

    // 6) Make F_CONNBOND between atoms within QM region
    modifyQMMMTwoCenterInteractions(mtop);

    // 7) Remove angles and settles containing 2 or more QM atoms
    modifyQMMMThreeCenterInteractions(mtop);

    // 8) Remove dihedrals containing 3 or more QM atoms
    modifyQMMMFourCenterInteractions(mtop);

    // 9) Build vector containing pairs of bonded QM - MM atoms (Link frontier)
    buildQMMMLink(mtop);

    // finalize topology
    mtop->finalize();
}

const QMMMTopologyInfo& QMMMTopologyPreprocessor::topInfo() const
{
    return topInfo_;
}

ArrayRef<const int> QMMMTopologyPreprocessor::atomNumbers() const
{
    return atomNumbers_;
}

ArrayRef<const real> QMMMTopologyPreprocessor::atomCharges() const
{
    return atomCharges_;
}

ArrayRef<const LinkFrontier> QMMMTopologyPreprocessor::linkFrontier() const
{
    return linkFrontier_;
}

bool QMMMTopologyPreprocessor::isQMAtom(Index globalAtomIndex)
{
    return (qmIndices_.find(globalAtomIndex) != qmIndices_.end());
}

void QMMMTopologyPreprocessor::splitQMblocks(gmx_mtop_t* mtop)
{

    // Global counter of atoms
    Index iAt = 0;

    /* Counter of molecules point to the specific moltype
     * i.e molblock 0 has 2 molecules have moltype 0 and molblock 2 has 1 additional molecule of type 0
     * then will be numMoleculesOfType[0] = 2 + 1 = 3
     * That counter is needed to decide whether we should create a new moltype for the molecule contatining QM atoms
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
        // Initialize block as non-QM first
        bQMBlock_.push_back(false);

        // Pointer to current block
        gmx_molblock_t* molBlock = &mtop->molblock[molBlockIndex];

        // Number of atoms in all molecules of current block
        const int numAtomsInMolecule = mtop->moltype[molBlock->type].atoms.nr;

        // Loop over all molecules in molBlock
        // mol - index of current molecule in molBlock
        for (int mol = 0; mol < molBlock->nmol; mol++)
        {

            // search for QM atoms in current molecule
            bool bQMMM = false;
            for (int i = 0; i < numAtomsInMolecule; i++)
            {
                if (isQMAtom(iAt))
                {
                    bQMMM = true;
                }
                iAt++;
            }

            // Apparently current molecule (molBlockIndex, mol) contains QM atoms
            // We should split it from the current block and create new blocks
            // For that molecule and all molecules after it new block will be created
            if (bQMMM)
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
                        bQMBlock_[molBlockIndex] = false;
                        bQMBlock_.push_back(true);
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
                        bQMBlock_[molBlockIndex] = true;
                    }
                }
                else
                {
                    bQMBlock_[molBlockIndex] = true;
                }

                // Create a copy of a moltype for a molecule
                // containing QM atoms and append it in the end of the moltype vector
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

                    // Set the molecule type for the QMMM molblock
                    molBlock->type = mtop->moltype.size() - 1;
                }
            }
        }
    }

    // Call finalize() to rebuild Block Indicies or else atoms lookup will fail
    mtop->finalize();
}

void QMMMTopologyPreprocessor::removeQMClassicalCharges(gmx_mtop_t* mtop)
{
    // Loop over all atoms and remove charge if they are QM atoms.
    // Sum-up total removed charge and remaning charge on MM atoms
    // Build atomCharges_ vector
    int molBlockIndex = 0;
    for (int i = 0; i < mtop->natoms; i++)
    {
        int indexInMolecule;
        mtopGetMolblockIndex(*mtop, i, &molBlockIndex, nullptr, &indexInMolecule);
        t_atom* atom = &mtop->moltype[mtop->molblock[molBlockIndex].type].atoms.atom[indexInMolecule];
        if (isQMAtom(i))
        {
            topInfo_.totalClassicalChargeOfQMAtoms += atom->q;
            atom->q  = 0.0;
            atom->qB = 0.0;
        }
        else
        {
            topInfo_.remainingMMCharge += atom->q;
        }

        atomCharges_.push_back(atom->q);
    }
}

void QMMMTopologyPreprocessor::addQMLJExclusions(gmx_mtop_t* mtop)
{
    // Add all QM atoms to the mtop->intermolecularExclusionGroup
    mtop->intermolecularExclusionGroup.reserve(mtop->intermolecularExclusionGroup.size()
                                               + qmIndices_.size());
    for (auto i : qmIndices_)
    {
        mtop->intermolecularExclusionGroup.push_back(i);
        topInfo_.numExclusionsMade++;
    }
}

void QMMMTopologyPreprocessor::buildQMMMAtomNumbers(gmx_mtop_t* mtop)
{
    // Save to atomNumbers_ atom numbers of all atoms
    AtomIterator atoms(*mtop);
    while (atoms->globalAtomNumber() < mtop->natoms)
    {
        // Check if we have valid atomnumbers
        if (atoms->atom().atomnumber < 0)
        {
            gmx_fatal(FARGS,
                      "Atoms %d does not have atomic number needed for QMMM. Check atomtypes "
                      "section in your topology or forcefield.",
                      atoms->globalAtomNumber());
        }

        atomNumbers_.push_back(atoms->atom().atomnumber);
        atoms++;
    }

    // Save in topInfo_ number of QM and MM atoms
    topInfo_.numQMAtoms += gmx::ssize(qmIndices_);
    topInfo_.numMMAtoms += mtop->natoms - gmx::ssize(qmIndices_);
}

void QMMMTopologyPreprocessor::modifyQMMMTwoCenterInteractions(gmx_mtop_t* mtop)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains QM atoms
        if (bQMBlock_[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index if the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (int ftype = 0; ftype < F_NRE; ftype++)
            {
                // If not bonded interaction or F_CONNBONDS, or some form of Restraints,
                // or not pair interaction, or no interactions of that type: then go the next type
                if (!(interaction_function[ftype].flags & IF_BOND) || ftype == F_CONNBONDS
                    || ftype == F_RESTRBONDS || ftype == F_HARMONIC || ftype == F_DISRES || ftype == F_ORIRES
                    || ftype == F_ANGRESZ || NRAL(ftype) != 2 || molType->ilist[ftype].empty())
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

                    // If both atoms are QM and it is IF_CHEMBOND then convert it to F_CONNBONDS
                    if (isQMAtom(molType->ilist[ftype].iatoms[j + 1] + start)
                        && isQMAtom(molType->ilist[ftype].iatoms[j + 2] + start))
                    {

                        // Add chemical bond to the F_CONNBONDS (bond type 5)
                        if (IS_CHEMBOND(ftype))
                        {
                            // Bond type is not used in F_CONNBONDS, so for generated bonds we set it to -1
                            const int connBondsType = -1;

                            // Add new CONNBOND between atoms
                            molType->ilist[F_CONNBONDS].iatoms.push_back(connBondsType);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 2]);

                            topInfo_.numConnBondsAdded++;
                        }

                        // Since the current bond is not copied into iatomsBuf it will be removed from the topology
                        topInfo_.numBondsRemoved++;
                    }
                    else
                    {
                        // If one of atoms is not QM then copy interaction into iatomsBuf
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
}

void QMMMTopologyPreprocessor::buildQMMMLink(gmx_mtop_t* mtop)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains QM atoms
        if (bQMBlock_[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - gloabl index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (int ftype = 0; ftype < F_NRE; ftype++)
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

                    // Update Link Frontier List if one of the atoms QM and one MM
                    if (isQMAtom(a1) && !isQMAtom(a2))
                    {
                        linkFrontier_.push_back({ a1, a2 });
                        topInfo_.numLinkBonds++;
                    }
                    if (isQMAtom(a2) && !isQMAtom(a1))
                    {
                        linkFrontier_.push_back({ a2, a1 });
                        topInfo_.numLinkBonds++;
                    }

                    // Check if it is constrained bond within QM subsystem
                    if (isQMAtom(a2) && isQMAtom(a1) && (interaction_function[ftype].flags & IF_CONSTRAINT))
                    {
                        topInfo_.numConstrainedBondsInQMSubsystem++;
                    }
                }
            }
        }
    }
}

void QMMMTopologyPreprocessor::modifyQMMMThreeCenterInteractions(gmx_mtop_t* mtop)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains QM atoms
        if (bQMBlock_[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (int ftype = 0; ftype < F_NRE; ftype++)
            {
                // If not bonded interaction or Restraints
                // or not three-particle interaction or no interactions of that type
                // and not Settle: then go the next type
                if ((!(interaction_function[ftype].flags & IF_BOND) || ftype == F_RESTRANGLES
                     || NRAL(ftype) != 3 || molType->ilist[ftype].empty())
                    && (ftype != F_SETTLE))
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
                    // Calculate number of qm atoms in the interaction
                    int numQm = 0;
                    for (int k = 1; k <= NRAL(ftype); k++)
                    {
                        if (isQMAtom(molType->ilist[ftype].iatoms[j + k] + start))
                        {
                            numQm++;
                        }
                    }

                    // If at least 2 atoms are QM then remove interaction
                    if (numQm >= 2)
                    {
                        // If this is SETTLE then replace it with two F_CONNBONDS
                        if (ftype == F_SETTLE)
                        {
                            // Bond type is not used in F_CONNBONDS, so for generated bonds we set it to -1
                            const int connBondsType = -1;

                            // Add CONNBOND between atoms 1 and 2 first
                            molType->ilist[F_CONNBONDS].iatoms.push_back(connBondsType);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 2]);
                            // Then between atoms 1 and 3
                            molType->ilist[F_CONNBONDS].iatoms.push_back(connBondsType);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 1]);
                            molType->ilist[F_CONNBONDS].iatoms.push_back(
                                    molType->ilist[ftype].iatoms[j + 3]);

                            topInfo_.numConnBondsAdded += 2;
                            topInfo_.numSettleRemoved++;
                        }
                        else
                        {
                            // If it is normal angle then it should be just removed
                            topInfo_.numAnglesRemoved++;
                        }
                    }
                    else
                    {
                        // If several (>1) atoms is not QM then preserve interaction in the iatomsBuf
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
}


void QMMMTopologyPreprocessor::modifyQMMMFourCenterInteractions(gmx_mtop_t* mtop)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains QM atoms
        if (bQMBlock_[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (int ftype = 0; ftype < F_NRE; ftype++)
            {
                // If not bonded interaction or Restraints
                // or not four-particle interaction or no interactions of that type: then go the next type
                if (!(interaction_function[ftype].flags & IF_BOND) || ftype == F_RESTRDIHS
                    || NRAL(ftype) != 4 || molType->ilist[ftype].empty())
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
                    // Calculate number of qm atoms in the interaction
                    int numQm = 0;
                    for (int k = 1; k <= NRAL(ftype); k++)
                    {
                        if (isQMAtom(molType->ilist[ftype].iatoms[j + k] + start))
                        {
                            numQm++;
                        }
                    }

                    // If at least 3 atoms are QM then remove interaction
                    if (numQm >= 3)
                    {
                        topInfo_.numDihedralsRemoved++;
                    }
                    else
                    {
                        // If several (>1) atoms is not QM then preserve interaction in the iatomsBuf
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
}

void QMMMTopologyPreprocessor::modifyQMMMVirtualSites(gmx_mtop_t* mtop)
{
    // Loop over all blocks in topology
    // molBlockIndex - index of current block in mtop
    for (size_t molBlockIndex = 0; molBlockIndex < mtop->molblock.size(); molBlockIndex++)
    {
        // check if current block contains QM atoms
        if (bQMBlock_[molBlockIndex])
        {
            // molType - strucutre with current block type
            gmx_moltype_t* molType = &mtop->moltype[mtop->molblock[molBlockIndex].type];

            // start - global index of the first atom in the current block
            int start = mtop->moleculeBlockIndices[molBlockIndex].globalAtomStart;

            // loop over all interaction types
            for (int ftype = 0; ftype < F_NRE; ftype++)
            {
                // If this is VSite interaction and ilist is not empty
                if (IS_VSITE(ftype) && !molType->ilist[ftype].empty())
                {
                    // Number of elements in the iatoms array for the current ftype
                    const int numInteractionElements = NRAL(ftype) + 1;

                    // Loop over all interactions of ftype
                    for (int j = 0; j < molType->ilist[ftype].size(); j += numInteractionElements)
                    {
                        // Calculate number of qm atoms in the interaction
                        int numQm = 0;
                        // Here k starts from 2 because first atom in the interaction is an actuall vsite index
                        for (int k = 2; k <= NRAL(ftype); k++)
                        {
                            if (isQMAtom(molType->ilist[ftype].iatoms[j + k] + start))
                            {
                                numQm++;
                            }
                        }

                        // If all atoms froming that virtual site are QM atoms
                        // then remove classical charge from that virtual site
                        if (numQm == (NRAL(ftype) - 1))
                        {
                            topInfo_.numVirtualSitesModified++;
                            topInfo_.totalClassicalChargeOfQMAtoms +=
                                    molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].q;
                            molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].q  = 0.0;
                            molType->atoms.atom[molType->ilist[ftype].iatoms[j + 1]].qB = 0.0;
                        }
                    }
                }
            }
        }
    }
}

} // namespace gmx
