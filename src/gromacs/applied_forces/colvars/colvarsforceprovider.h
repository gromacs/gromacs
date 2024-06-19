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
 * Declares the force provider for colvars
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */

#ifndef GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
#define GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H


#include <cstdint>

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/domdec/localatomset.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/real.h"

#include "colvarproxygromacs.h"

enum class PbcType : int;
struct t_commrec;


namespace gmx
{
class KeyValueTreeObject;
class KeyValueTreeObjectBuilder;
class LocalAtomSetManager;
class MDLogger;
struct MDModulesWriteCheckpointData;


/*! \internal
 * \brief Parameters defining the internal colvars force provider state.
 */
struct ColvarsForceProviderState
{

    /*! \brief Indicate if a colvars state was read.
     */
    bool stateRead_ = false;

    /*! \brief The number of colvars atoms.
     */
    std::int64_t nColvarsAtoms_ = 0;

    /*! \brief String naming variable holding the number of colvars atoms.
     * \note Changing this name will break backwards compability for checkpoint file writing.
     */
    static const std::string sc_nColvarsAtomsName_;

    //! Last known whole positions of the colvars atoms
    //! \todo Change the type to a standard one to avoid memory leak.
    rvec* xOldWhole_ = nullptr;

    /*! \brief String naming variable holding the last known whole positions of the colvars atoms
     * \note Changing this name will break backwards compability for checkpoint file writing.
     */
    static const std::string sc_xOldWholeName_;

    /*! \brief Content of the unformatted Colvars state file.
     */
    std::vector<unsigned char> colvarStateFile_;

    /*! \brief String naming variable holding the content of the unformatted Colvars state file.
     * \note Changing this name will break backwards compability for checkpoint file writing.
     */
    static const std::string sc_colvarStateFileName_;

    /*! \brief String naming variable holding the size of the unformatted Colvars state file.
     * \note Changing this name will break backwards compability for checkpoint file writing.
     */
    static const std::string sc_colvarStateFileSizeName_;

    /*! \brief Write internal colvars data into a key value tree.
     * The entries to the kvt are identified with identifier, so that a variable
     * is indentified with the key "identifier-variablename"
     *
     * \param[in] kvtBuilder enables writing to the Key-Value-Tree
     *                              the state is written to
     *
     * \param[in] identifier denotes the module that is checkpointing the data
     */
    void writeState(KeyValueTreeObjectBuilder kvtBuilder, const std::string& identifier) const;

    /*! \brief Read the internal parameters from the checkpoint file on master
     * \param[in] kvtData holding the checkpoint information
     * \param[in] identifier identifies the data in a key-value-tree
     */
    void readState(const KeyValueTreeObject& kvtData, const std::string& identifier);
};


/*! \internal \brief
 * Implements IForceProvider for colvars.
 * Override the ColvarProxyGromacs generic class for the communication.
 */
class ColvarsForceProvider final : public ColvarProxyGromacs, public IForceProvider
{

private:
    //! The total bias energy on all colvars atoms.
    double biasEnergy;

    //! Is this a neighbor-search step?
    bool gmxBNS;


    // Node-local bookkepping data
    //! The colvars atom indices
    std::unique_ptr<gmx::LocalAtomSet> colvarsAtoms;
    //! Total number of Colvars atoms
    int nColvarsAtoms = 0;
    //! Unwrapped positions for all Colvars atoms, communicated to all nodes.
    rvec* xColvars = nullptr;
    //! Shifts for all Colvars atoms, to make molecule(s) whole.
    ivec* xColvarsShifts = nullptr;
    //! Extra shifts since last DD step.
    ivec* xColvarsEshifts = nullptr;
    //! Old positions for all Colvars atoms on master.
    rvec* xColvarsOldWhole = nullptr;
    //! Bias forces on all Colvars atoms
    rvec* fColvars = nullptr;

    //! Struct holding the information stored in the checkpoint file
    ColvarsForceProviderState stateToCheckpoint_;


public:
    friend class cvm::atom;

    /*! \brief Construct ColvarsForceProvider from its parameters
     *
     * \param[in] colvarsConfigString Content of the colvars input file.
     * \param[in] atoms Atoms topology
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] logger GROMACS logger instance
     * \param[in] inputStrings Input files stored as string in the TPR's KVT
     * \param[in] ensembleTemperature the constant ensemble temperature
     * \param[in] seed The colvars seed for random number generator
     * \param[in] localAtomSetManager Atom Manager to retrieve Colvars index atoms
     * \param[in] cr Communication Record
     * \param[in] simulationTimeStep The simulation time step
     * \param[in] colvarsCoords The colvars atoms coordinates retrived from the TPR's KVT
     * \param[in] outputPrefix The prefix for output colvars files
     * \param[in] state The state of colvars force provider to be written in the checkpoint
     */
    ColvarsForceProvider(const std::string&                        colvarsConfigString,
                         t_atoms                                   atoms,
                         PbcType                                   pbcType,
                         const MDLogger*                           logger,
                         const std::map<std::string, std::string>& inputStrings,
                         real                                      ensembleTemperature,
                         int                                       seed,
                         LocalAtomSetManager*                      localAtomSetManager,
                         const t_commrec*                          cr,
                         double                                    simulationTimeStep,
                         const std::vector<RVec>&                  colvarsCoords,
                         const std::string&                        outputPrefix,
                         const ColvarsForceProviderState&          state);

    ~ColvarsForceProvider() override;

    /*! \brief Calculate colvars forces
     * \param[in] forceProviderInput input for force provider
     * \param[out] forceProviderOutput output for force provider
     */
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    //! Compute virial tensor for position r and force f, and add to matrix vir
    static void addVirialTerm(matrix vir, const rvec& f, const gmx::RVec& x);

    /*! \brief Write internal colvars data to checkpoint file.
     * \param[in] checkpointWriting enables writing to the Key-Value-Tree
     *                              that is used for storing the checkpoint
     *                              information
     * \param[in] moduleName names the module that is checkpointing this force-provider
     *
     * \note The provided state to checkpoint has to change if checkpointing
     *       is moved before the force provider call in the MD-loop.
     */
    void writeCheckpointData(MDModulesWriteCheckpointData checkpointWriting, const std::string& moduleName);

    /*! \brief Process atomsRedistributedSignal notification during mdrun.
     * \param[in] atomsRedistributedSignal signal recieved
     */
    void processAtomsRedistributedSignal(const MDModulesAtomsRedistributedSignal& atomsRedistributedSignal);


    //! From colvarproxy

    /*! \brief add energy to the total count of bias energy biasEnergy
     * \param[in] energy the value of energy to add
     *
     */
    void add_energy(cvm::real energy) override;
};

} // namespace gmx

#endif // GMX_APPLIED_FORCES_COLVARSFORCEPROVIDER_H
