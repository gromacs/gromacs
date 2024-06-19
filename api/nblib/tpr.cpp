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
/*! \inpublicapi \file
 * \brief
 * Implements nblib tpr reading
 *
 * \author Victor Holanda <victor.holanda@cscs.ch>
 * \author Joe Jordan <ejjordan@kth.se>
 * \author Prashanth Kanduri <kanduri@cscs.ch>
 * \author Sebastian Keller <keller@cscs.ch>
 */

#include "nblib/tpr.h"

#include <filesystem>

#include "listed_forces/conversionscommon.h"

#include "gromacs/fileio/tpxio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/listed_forces/listed_forces.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/forcerec.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/mdtypes/simulation_workload.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/logger.h"

#include "nblib/box.h"
#include "nblib/exception.h"
#include "nblib/listed_forces/bondtypes.h"

namespace nblib
{

TprReader::TprReader(std::string filename)
{
    t_inputrec              inputRecord;
    t_state                 globalState;
    gmx_mtop_t              molecularTopology;
    gmx::SimulationWorkload simulationWorkload;

    // If the file does not exist, this function will throw
    PartialDeserializedTprFile partialDeserializedTpr =
            read_tpx_state(filename, &inputRecord, &globalState, &molecularTopology);

    // init forcerec
    t_forcerec          forceRecord;
    t_commrec           commrec{};
    gmx::ForceProviders forceProviders;
    forceRecord.forceProviders = &forceProviders;
    init_forcerec(nullptr,
                  gmx::MDLogger(),
                  simulationWorkload,
                  &forceRecord,
                  inputRecord,
                  molecularTopology,
                  &commrec,
                  globalState.box,
                  nullptr,
                  nullptr,
                  gmx::ArrayRef<const std::string>{},
                  -1);

    nonbondedParameters_ = makeNonBondedParameterLists(molecularTopology.ffparams.atnr,
                                                       molecularTopology.ffparams.iparams,
                                                       forceRecord.haveBuckingham);

    gmx_localtop_t localtop(molecularTopology.ffparams);
    gmx_mtop_generate_local_top(molecularTopology, &localtop, false);
    exclusionListElements_ = std::vector<int>(localtop.excls.elementsView().begin(),
                                              localtop.excls.elementsView().end());
    exclusionListRanges_   = std::vector<int>(localtop.excls.listRangesView().begin(),
                                            localtop.excls.listRangesView().end());

    int                           ntopatoms = molecularTopology.natoms;
    std::unique_ptr<gmx::MDAtoms> mdAtoms =
            gmx::makeMDAtoms(nullptr, molecularTopology, inputRecord, false);
    atoms2md(molecularTopology, inputRecord, -1, {}, ntopatoms, mdAtoms.get());
    const double initMassLambda =
            (inputRecord.efep == FreeEnergyPerturbationType::No
                     ? 0.0
                     : inputRecord.fepvals->initialLambda(FreeEnergyPerturbationCouplingType::Mass));
    update_mdatoms(mdAtoms->mdatoms(), initMassLambda);
    auto numParticles = mdAtoms->mdatoms()->nr;
    charges_.resize(numParticles);
    particleTypeIdOfAllParticles_.resize(numParticles);
    inverseMasses_.resize(numParticles);
    for (int i = 0; i < numParticles; i++)
    {
        charges_[i]                      = mdAtoms->mdatoms()->chargeA[i];
        particleTypeIdOfAllParticles_[i] = mdAtoms->mdatoms()->typeA[i];
        inverseMasses_[i]                = mdAtoms->mdatoms()->invmass[i];
    }
    particleInteractionFlags_ = forceRecord.atomInfo;

    if (TRICLINIC(globalState.box))
    {
        throw InputException("Only rectangular unit-cells are supported here");
    }
    boxX_ = globalState.box[XX][XX];
    boxY_ = globalState.box[YY][YY];
    boxZ_ = globalState.box[ZZ][ZZ];
    coordinates_.assign(globalState.x.begin(), globalState.x.end());
    velocities_.assign(globalState.v.begin(), globalState.v.end());

    // Copy listed interactions data
    int listedInteractionCount = gmx_mtop_interaction_count(molecularTopology, IF_BOND)
                                 + gmx_mtop_interaction_count(molecularTopology, IF_PAIR)
                                 + gmx_mtop_interaction_count(molecularTopology, IF_DIHEDRAL);
    if (listedInteractionCount != 0)
    {
        InteractionDefinitions interactionDefinitions = localtop.idef;
        listedInteractionData_ = convertToNblibInteractions(interactionDefinitions);
    }
}

Box TprReader::getBox() const
{
    return Box(boxX_, boxY_, boxZ_);
}

} // namespace nblib
