/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
 * Helper functions for topology generation.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 */

#include "gmxpre.h"

#include "topologyhelpers.h"

#include <vector>

#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"

namespace gmx
{
namespace test
{

void addNWaterMolecules(gmx_mtop_t* mtop, int numWaters)
{
    // Target distance between oxygen and hydrogens
    constexpr real dOH = 0.09572;
    // Target distance between hydrogens
    constexpr real dHH = 0.15139;

    gmx_moltype_t waterMoleculeType;

    waterMoleculeType.atoms.nr   = NRAL(InteractionFunction::SETTLE);
    std::vector<int>& iatoms     = waterMoleculeType.ilist[InteractionFunction::SETTLE].iatoms;
    const int         settleType = 0;
    iatoms.push_back(settleType);
    iatoms.push_back(0);
    iatoms.push_back(1);
    iatoms.push_back(2);

    init_t_atoms(&waterMoleculeType.atoms, NRAL(InteractionFunction::SETTLE), false);
    for (int i = 0; i < NRAL(InteractionFunction::SETTLE); ++i)
    {
        waterMoleculeType.atoms.atom[i].m = (i % 3 == 0) ? 16 : 1;
    }
    mtop->moltype.push_back(waterMoleculeType);
    // gmx_moltype_t does not have proper move-construct semantics
    waterMoleculeType.atoms.atom     = nullptr;
    waterMoleculeType.atoms.atomname = nullptr;
    waterMoleculeType.atoms.resinfo  = nullptr;

    mtop->molblock.emplace_back();
    mtop->molblock.back().type = mtop->moltype.size() - 1;
    mtop->molblock.back().nmol = numWaters;
    mtop->natoms += waterMoleculeType.atoms.nr * mtop->molblock.back().nmol;

    // Set up the SETTLE parameters.
    t_iparams iparams;
    iparams.settle.doh = dOH;
    iparams.settle.dhh = dHH;
    mtop->ffparams.iparams.push_back(iparams);
}

void addEthane(gmx_mtop_t* mtop, const bool bondsAreFlexible)
{
    t_iparams iparams;
    if (bondsAreFlexible)
    {
        // Add an harmonic bond type
        iparams.harmonic = { 1.0, 1000.0, 1.0, 1000.0 };
    }
    else
    {
        // Add a constraint type
        iparams.constr = { 0.5, 0.5 };
    }
    mtop->ffparams.iparams.push_back(iparams);
    const int twoAtomInteractionType = 0;
    // Add an angle type
    iparams.harmonic = { 120, 1000, 120, 1000 };
    mtop->ffparams.iparams.push_back(iparams);
    const int angleInteractionType = 1;
    // Add an LJ-14 type
    iparams.lj14 = { 1, 100, 1, 100 };
    mtop->ffparams.iparams.push_back(iparams);
    const int oneFourInteractionType = 2;

    const int numAtoms = 8;
    const int carbon0  = 0;
    const int carbon1  = 4;
    // Atoms 1,2,3,5,6,7 are hydrogens.

    gmx_moltype_t ethaneMoleculeType;

    // Describe the interactions
    // clang-format off
    const InteractionFunction twoAtomInteraction = (bondsAreFlexible ? InteractionFunction::Bonds : InteractionFunction::Constraints);
    ethaneMoleculeType.ilist[twoAtomInteraction].iatoms = { twoAtomInteractionType, carbon0, 1,
                                                  twoAtomInteractionType, carbon0, 2,
                                                  twoAtomInteractionType, carbon0, 3,
                                                  twoAtomInteractionType, carbon1, 5,
                                                  twoAtomInteractionType, carbon1, 6,
                                                  twoAtomInteractionType, carbon1, 7 };
    ethaneMoleculeType.ilist[InteractionFunction::Angles].iatoms = { angleInteractionType, 1, carbon0, 2,
                                                  angleInteractionType, 1, carbon0, 3,
                                                  angleInteractionType, 2, carbon0, 3,
                                                  angleInteractionType, 5, carbon1, 6,
                                                  angleInteractionType, 5, carbon1, 7,
                                                  angleInteractionType, 6, carbon1, 7 };
    ethaneMoleculeType.ilist[InteractionFunction::LennardJones14].iatoms = { oneFourInteractionType, 1, 5,
                                                oneFourInteractionType, 1, 6,
                                                oneFourInteractionType, 1, 7,
                                                oneFourInteractionType, 2, 5,
                                                oneFourInteractionType, 2, 6,
                                                oneFourInteractionType, 2, 7,
                                                oneFourInteractionType, 3, 5,
                                                oneFourInteractionType, 3, 6,
                                                oneFourInteractionType, 3, 7 };
    // clang-format on

    // Set up the atoms, including arbitrary uniform charges
    ethaneMoleculeType.atoms.atomname = nullptr;
    ethaneMoleculeType.atoms.nres     = 0;
    ethaneMoleculeType.atoms.nr       = numAtoms;
    snew(ethaneMoleculeType.atoms.atom, numAtoms);
    for (int i = 0; i != numAtoms; ++i)
    {
        ethaneMoleculeType.atoms.atom[i].q  = 1.0;
        ethaneMoleculeType.atoms.atom[i].qB = 1.0;
    }
    mtop->moltype.push_back(ethaneMoleculeType);
    ethaneMoleculeType.atoms.atom    = nullptr;
    ethaneMoleculeType.atoms.resinfo = nullptr;

    mtop->molblock.emplace_back();
    mtop->molblock.back().type = mtop->moltype.size() - 1;
    mtop->molblock.back().nmol = 1;
    mtop->natoms += ethaneMoleculeType.atoms.nr * mtop->molblock.back().nmol;
}

void addIntermolecularInteractionBonds(gmx_mtop_t* mtop)
{
    mtop->bIntermolecularInteractions = true;
    mtop->intermolecular_ilist        = std::make_unique<InteractionLists>();
    std::vector<int>& iatoms   = (*mtop->intermolecular_ilist)[InteractionFunction::Bonds].iatoms;
    const int         bondType = 0;
    iatoms.push_back(bondType);
    iatoms.push_back(0);
    iatoms.push_back(1);
    iatoms.push_back(bondType);
    iatoms.push_back(2);
    iatoms.push_back(3);
    iatoms.push_back(bondType);
    iatoms.push_back(4);
    iatoms.push_back(5);
}

} // namespace test
} // namespace gmx
