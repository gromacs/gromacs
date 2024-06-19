/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Implements the force provider for colvars
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#include "colvarsforceprovider.h"

#include <cstddef>
#include <cstdint>

#include <array>
#include <string>

#include "external/colvars/colvars_memstream.h"

#include "gromacs/applied_forces/colvars/colvarproxygromacs.h"
#include "gromacs/compat/pointers.h"
#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/broadcaststructs.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/smalloc.h"

enum class PbcType : int;


namespace gmx
{
class MDLogger;

/********************************************************************
 * ColvarsForceProviderState
 */

const std::string ColvarsForceProviderState::sc_nColvarsAtomsName_ = "nColvarsAtoms";

const std::string ColvarsForceProviderState::sc_xOldWholeName_ = "xOldWhole";

const std::string ColvarsForceProviderState::sc_colvarStateFileName_ = "colvarStateFile";

const std::string ColvarsForceProviderState::sc_colvarStateFileSizeName_ = "colvarStateFileSize";

void ColvarsForceProviderState::writeState(KeyValueTreeObjectBuilder kvtBuilder,
                                           const std::string&        identifier) const
{
    writeKvtCheckpointValue(nColvarsAtoms_, sc_nColvarsAtomsName_, identifier, kvtBuilder);

    // Write colvars atoms coords
    auto DoubleArrayAdder = kvtBuilder.addUniformArray<double>(sc_xOldWholeName_);
    for (int i = 0; i < nColvarsAtoms_; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            DoubleArrayAdder.addValue(static_cast<double>(xOldWhole_[i][j]));
        }
    }


    writeKvtCheckpointValue(
            static_cast<int64_t>(colvarStateFile_.size()), sc_colvarStateFileSizeName_, identifier, kvtBuilder);

    // Write unformatted Colvars state file, one character at a time
    auto charArrayAdder = kvtBuilder.addUniformArray<unsigned char>(sc_colvarStateFileName_);
    for (const unsigned char& c : colvarStateFile_)
    {
        charArrayAdder.addValue(c);
    }
}

void ColvarsForceProviderState::readState(const KeyValueTreeObject& kvtData, const std::string& identifier)
{

    stateRead_ = true;

    readKvtCheckpointValue(
            compat::make_not_null(&nColvarsAtoms_), sc_nColvarsAtomsName_, identifier, kvtData);


    // Read colvars atoms coords
    auto kvtDoubleArray = kvtData[sc_xOldWholeName_].asArray().values();


    // Make sure the coordinates saved are consistent with the dimensions
    if (kvtDoubleArray.size() % DIM != 0)
    {
        GMX_THROW(InconsistentInputError(
                "Coordinates saved in the checkpoint file are in the wrong format."));
    }

    snew(xOldWhole_, nColvarsAtoms_);
    for (size_t i = 0; i < kvtDoubleArray.size() / DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            xOldWhole_[i][j] = static_cast<real>(kvtDoubleArray[i * DIM + j].cast<double>());
        }
    }

    int64_t colvarStateFileSize_ = 0L;
    readKvtCheckpointValue(
            compat::make_not_null(&colvarStateFileSize_), sc_colvarStateFileSizeName_, identifier, kvtData);

    // Read Colvars state file; use explicit loop because kvt types don't support std::copy
    auto charArray = kvtData[sc_colvarStateFileName_].asArray().values();
    colvarStateFile_.resize(colvarStateFileSize_);
    auto it = colvarStateFile_.begin();
    for (const auto& c : charArray)
    {
        *it = c.cast<unsigned char>();
        it++;
    }
}


/********************************************************************
 * ColvarsForceProvider
 */

ColvarsForceProvider::ColvarsForceProvider(const std::string& colvarsConfigString,
                                           t_atoms            atoms,
                                           PbcType            pbcType,
                                           const MDLogger*    logger,
                                           const std::map<std::string, std::string>& inputStrings,
                                           real                             ensembleTemperature,
                                           int                              seed,
                                           LocalAtomSetManager*             localAtomSetManager,
                                           const t_commrec*                 cr,
                                           double                           simulationTimeStep,
                                           const std::vector<RVec>&         colvarsCoords,
                                           const std::string&               outputPrefix,
                                           const ColvarsForceProviderState& state) :
    ColvarProxyGromacs(colvarsConfigString, atoms, pbcType, logger, MAIN(cr), inputStrings, ensembleTemperature, seed),
    stateToCheckpoint_(state)
{


    // From colvarproxy_system class
    // Total forces on each atom is not available in GROMACS
    total_force_requested = false;

    // Neighbor Search boolean activated during initialization
    gmxBNS = true;

    // Get GROMACS timestep (picosecond to femtosecond)
    set_integration_timestep(simulationTimeStep * 1000.0);

    // From colvaproxy_io
    output_prefix_str = outputPrefix;


    if (doParsing_)
    {
        colvars->setup_output();
    }


    // MPI initialisation

    // Initialise attributs for the MPI communication
    if (MAIN(cr))
    {
        // Retrieve the number of colvar atoms
        nColvarsAtoms = atoms_ids.size();
    }

    if (PAR(cr))
    {
        // Let the other nodes know the number of colvar atoms and their ids to construct a gmx::LocalAtomSet
        block_bc(cr->mpi_comm_mygroup, nColvarsAtoms);
        atoms_ids.resize(nColvarsAtoms);
        nblock_bc(cr->mpi_comm_mygroup, nColvarsAtoms, atoms_ids.data());

        // Initialise atoms_new_colvar_forces on non-MAIN nodes
        if (!MAIN(cr))
        {
            atoms_new_colvar_forces.resize(nColvarsAtoms);
        }
    }

    // Cast int into Index of the indices for the localAtomSetManager->add() function
    std::vector<Index> indexAtoms(atoms_ids.begin(), atoms_ids.end());
    colvarsAtoms = std::make_unique<LocalAtomSet>(localAtomSetManager->add(indexAtoms));


    snew(xColvars, nColvarsAtoms);
    snew(xColvarsShifts, nColvarsAtoms);
    snew(xColvarsEshifts, nColvarsAtoms);
    snew(fColvars, nColvarsAtoms);
    snew(xColvarsOldWhole, nColvarsAtoms);


    // Check state status (did we read a cpt file?)
    if (MAIN(cr))
    {
        if (stateToCheckpoint_.stateRead_)
        {
            if (stateToCheckpoint_.nColvarsAtoms_ != nColvarsAtoms)
            {
                cvm::error(
                        "Number of colvars atoms in the .cpt file differs from the one in .tpr "
                        "file");
            }

            // Copy back the last whole positions from the .cpt file
            for (int i = 0; i < nColvarsAtoms; i++)
            {
                copy_rvec(stateToCheckpoint_.xOldWhole_[i], xColvarsOldWhole[i]);
            }

            int errorCode = colvarproxy::setup();
            // Read input state file
            errorCode |= colvars->set_input_state_buffer(stateToCheckpoint_.colvarStateFile_);
            errorCode |= colvars->setup_input();

            if (errorCode != COLVARS_OK)
            {
                error("Error when initializing Colvars module.");
            }
        }
        else
        {
            // Initialize state variables
            stateToCheckpoint_.nColvarsAtoms_ = nColvarsAtoms;
            snew(stateToCheckpoint_.xOldWhole_, nColvarsAtoms);

            // Use input coords for the last whole positions.
            for (int i = 0; i < nColvarsAtoms; i++)
            {
                copy_rvec(colvarsCoords[i], xColvarsOldWhole[i]);
            }
        }
    }


    // // Communicate initial coordinates to all processes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, nColvarsAtoms, xColvarsOldWhole);
    }


    if (MAIN(cr) && cvm::debug())
    {
        cvm::log("atoms_ids = " + cvm::to_str(atoms_ids) + "\n");
        cvm::log("atoms_refcount = " + cvm::to_str(atoms_refcount) + "\n");
        cvm::log("positions = " + cvm::to_str(atoms_positions) + "\n");
        cvm::log("total_forces = " + cvm::to_str(atoms_total_forces) + "\n");
        cvm::log("atoms_new_colvar_forces = " + cvm::to_str(atoms_new_colvar_forces) + "\n");
        cvm::log(cvm::line_marker);
        log("Done initializing the colvars proxy object.\n");
    }

    if (MAIN(cr))
    {
        cvm::log(cvm::line_marker);
        cvm::log("End colvars Initialization.\n\n");
    }
}

ColvarsForceProvider::~ColvarsForceProvider()
{
    if (doParsing_)
    {
        post_run();
        sfree(stateToCheckpoint_.xOldWhole_);
    }
    sfree(xColvars);
    sfree(xColvarsShifts);
    sfree(xColvarsEshifts);
    sfree(fColvars);
    sfree(xColvarsOldWhole);
}

void ColvarsForceProvider::calculateForces(const ForceProviderInput& forceProviderInput,
                                           ForceProviderOutput*      forceProviderOutput)
{

    // Construct t_pbc struct
    set_pbc(&gmxPbc_, pbcType_, forceProviderInput.box_);

    const t_commrec* cr = &(forceProviderInput.cr_);
    // Local atom coords
    const gmx::ArrayRef<const gmx::RVec> x = forceProviderInput.x_;
    // Local atom coords (coerced into into old gmx type)
    const rvec* xPointer = &(x.data()->as_vec());
    const auto& box      = forceProviderInput.box_;

    colvars->it = forceProviderInput.step_;


    // Eventually there needs to be an interface to update local data upon neighbor search
    // We could check if by chance all atoms are in one node, and skip communication
    communicate_group_positions(cr,
                                xColvars,
                                xColvarsShifts,
                                xColvarsEshifts,
                                gmxBNS,
                                xPointer,
                                colvarsAtoms->numAtomsGlobal(),
                                colvarsAtoms->numAtomsLocal(),
                                colvarsAtoms->localIndex().data(),
                                colvarsAtoms->collectiveIndex().data(),
                                xColvarsOldWhole,
                                box);


    // Communicate_group_positions takes care of removing shifts (unwrapping)
    // in single node jobs, communicate_group_positions() is efficient and adds no overhead

    if (MAIN(cr))
    {
        // On non-MAIN nodes, jump directly to applying the forces

        // Zero the forces on the atoms, so that they can be accumulated by the colvars.
        for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++)
        {
            atoms_new_colvar_forces[i].x         = atoms_new_colvar_forces[i].y =
                    atoms_new_colvar_forces[i].z = 0.0;
        }

        // Copy the global Colvars atoms coordinates gathered in xColvars to atom_positions array
        // for later used in the calc() function.
        for (size_t i = 0; i < atoms_ids.size(); i++)
        {
            atoms_positions[i] = cvm::rvector(xColvars[i][0], xColvars[i][1], xColvars[i][2]);
        }

        biasEnergy = 0.0;
        // Call the collective variable module to fill atoms_new_colvar_forces
        if (colvars->calc() != COLVARS_OK)
        {
            cvm::error("Error calling colvars->calc()\n");
        }

        // Copy the forces to C array for broadcasting
        for (int i = 0; i < nColvarsAtoms; i++)
        {
            fColvars[i][0] = atoms_new_colvar_forces[i].x;
            fColvars[i][1] = atoms_new_colvar_forces[i].y;
            fColvars[i][2] = atoms_new_colvar_forces[i].z;
        }

        forceProviderOutput->enerd_.term[F_COM_PULL] += biasEnergy;

        // Copy last whole positions into State struct.
        for (int i = 0; i < nColvarsAtoms; i++)
        {
            copy_rvec(xColvarsOldWhole[i], stateToCheckpoint_.xOldWhole_[i]);
        }
    } // MAIN node


    // Broadcast the forces to all the nodes
    if (PAR(cr))
    {
        nblock_bc(cr->mpi_comm_mygroup, nColvarsAtoms, fColvars);
    }


    const gmx::ArrayRef<gmx::RVec>& fOut = forceProviderOutput->forceWithVirial_.force_;
    matrix                          localColvarsVirial     = { { 0 } };
    const auto&                     localColvarsIndex      = colvarsAtoms->localIndex();
    const auto&                     collectiveColvarsIndex = colvarsAtoms->collectiveIndex();
    // Loop through local atoms to aply the colvars forces
    for (gmx::Index l = 0; l < localColvarsIndex.ssize(); l++)
    {
        /* Get the right index of the local colvars atoms */
        int iLocal = localColvarsIndex[l];
        /* Index of this local atom in the collective colvars atom arrays */
        int iColvars = collectiveColvarsIndex[l];
        /* Add */
        rvec_inc(fOut[iLocal], fColvars[iColvars]);
        addVirialTerm(localColvarsVirial, fColvars[iColvars], xColvars[iColvars]);
    }

    forceProviderOutput->forceWithVirial_.addVirialContribution(localColvarsVirial);

    // Re-set the flag for proper update
    gmxBNS = false;
}

void ColvarsForceProvider::addVirialTerm(matrix vir, const rvec& f, const gmx::RVec& x)
{
    for (int j = 0; j < DIM; j++)
    {
        for (int m = 0; m < DIM; m++)
        {
            vir[j][m] -= 0.5 * f[j] * x[m];
        }
    }
}


void ColvarsForceProvider::add_energy(cvm::real energy)
{
    biasEnergy += energy;
}


void ColvarsForceProvider::writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting,
                                               const std::string&           moduleName)
{
    colvars->write_state_buffer(stateToCheckpoint_.colvarStateFile_);
    stateToCheckpoint_.writeState(checkpointWriting.builder_, moduleName);
}

void ColvarsForceProvider::processAtomsRedistributedSignal(const MDModulesAtomsRedistributedSignal& /*signal*/)
{
    // So far, just update the Neighbor Search boolean for the communicate_group_positions() in calculateForces()
    gmxBNS = true;
}


} // namespace gmx
