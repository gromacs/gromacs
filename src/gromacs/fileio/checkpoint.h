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

#ifndef GMX_FILEIO_CHECKPOINT_H
#define GMX_FILEIO_CHECKPOINT_H

#include <cstdint>
#include <cstdio>

#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxmpi.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

class energyhistory_t;
struct gmx_file_position_t;
struct ObservablesHistory;
struct t_commrec;
struct t_fileio;
struct t_inputrec;
class t_state;
struct t_trxframe;
enum class IntegrationAlgorithm : int;
enum class SwapType : int;
enum class LambdaWeightCalculation : int;
enum class CheckPointVersion : int;

namespace gmx
{

struct MDModulesNotifiers;
class KeyValueTreeObject;
class ReadCheckpointDataHolder;
class WriteCheckpointDataHolder;

//! The maximum number of atoms that can be stored in a checkpoint file
static constexpr int64_t sc_checkpointMaxAtomCount = std::numeric_limits<unsigned int>::max() / 3;

/*! \brief Read to a key-value-tree value used for checkpointing.
 *
 * \tparam ValueType
 *
 * \param[in] value the value to be checkpointed
 * \param[in] name name of the value to be checkpointed
 * \param[in] identifier uniquely identifies the module that is checkpointing
 *                       typically the module name
 * \param[in] kvt the key value tree to read from
 *
 * \throws InternalError if kvt does not contain requested value.
 * \note Triggers assertion if value type is not correct.
 */
template<typename ValueType>
void readKvtCheckpointValue(compat::not_null<ValueType*> value,
                            const std::string&           name,
                            const std::string&           identifier,
                            const KeyValueTreeObject&    kvt);
//! \copydoc readKvtCheckpointValue
extern template void readKvtCheckpointValue(compat::not_null<std::int64_t*> value,
                                            const std::string&              name,
                                            const std::string&              identifier,
                                            const KeyValueTreeObject&       kvt);
//! \copydoc readKvtCheckpointValue
extern template void readKvtCheckpointValue(compat::not_null<real*>   value,
                                            const std::string&        name,
                                            const std::string&        identifier,
                                            const KeyValueTreeObject& kvt);

/*! \brief Write to a key-value-tree used for checkpointing.
 *
 * \tparam ValueType
 *
 * \param[in] value name of the value to be checkpointed
 * \param[in] name the value to be checkpointed
 * \param[in] identifier uniquely identifies the module that is checkpointing
 *                       typically the module name
 * \param[in] kvtBuilder the key-value-tree builder used to store the checkpoint values
 */
template<typename ValueType>
void writeKvtCheckpointValue(const ValueType&          value,
                             const std::string&        name,
                             const std::string&        identifier,
                             KeyValueTreeObjectBuilder kvtBuilder);
//! \copydoc writeKvtCheckpointValue
extern template void writeKvtCheckpointValue(const std::int64_t&       value,
                                             const std::string&        name,
                                             const std::string&        identifier,
                                             KeyValueTreeObjectBuilder kvtBuilder);
//! \copydoc writeKvtCheckpointValue
extern template void writeKvtCheckpointValue(const real&               value,
                                             const std::string&        name,
                                             const std::string&        identifier,
                                             KeyValueTreeObjectBuilder kvtBuilder);

/*! \libinternal
 * \brief Provides the MDModules with the checkpointed data on the main rank.
 */
struct MDModulesCheckpointReadingDataOnMain
{
    //! The data of the MDModules that is stored in the checkpoint file
    const KeyValueTreeObject& checkpointedData_;
};

/*! \libinternal
 * \brief Provides the MDModules with the communication record to broadcast.
 */
struct MDModulesCheckpointReadingBroadcast
{
    //! The communicator
    MPI_Comm communicator_;
    //! Whether the run is executed in parallel
    bool isParallelRun_;
};

/*! \libinternal \brief Writing the MDModules data to a checkpoint file.
 */
struct MDModulesWriteCheckpointData
{
    //! Builder for the Key-Value-Tree to store the MDModule checkpoint data
    KeyValueTreeObjectBuilder builder_;
};

} // namespace gmx

/* the name of the environment variable to disable fsync failure checks with */
#define GMX_IGNORE_FSYNC_FAILURE_ENV "GMX_IGNORE_FSYNC_FAILURE"

// TODO Replace this mechanism with std::array<char, 1024> or similar.
#define CPTSTRLEN 1024

/*! \brief Enum of values that describe the contents of a cpt file
 * whose format matches a version number
 *
 * The enum helps the code be more self-documenting and ensure merges
 * do not silently resolve when two patches make the same bump. When
 * adding new functionality, add a new element just above Count
 * in this enumeration, and write code that does the right thing
 * according to the value of file_version.
 */
enum class CheckPointVersion : int
{
    //! Unknown initial version
    UnknownVersion0,
    //! Unknown version 1
    UnknownVersion1,
    //! Store magic number to validate file.
    AddMagicNumber,
    //! Store which sim this is.
    SafeSimulationPart,
    //! Store energy data.
    EkinDataAndFlags,
    //! Store current step.
    SafeSteps,
    //! Cutoff before version 4.5.
    Version45,
    //! Unknown version 7.
    UnknownVersion7,
    //! Store checksum.
    FileChecksumAndSize,
    //! Unknown version 9.
    UnknownVersion9,
    //! Safe NH chains for thermostat.
    NoseHooverThermostat,
    //! Safe NH chains for barostat.
    NoseHooverBarostat,
    //! Safe host info.
    HostInformation,
    //! Activate double build.
    DoublePrecisionBuild,
    //! Lambda stuff.
    LambdaStateAndHistory,
    //! ED stuff.
    EssentialDynamics,
    //! Swap state stuff.
    SwapState,
    //! AWH history flags added.
    AwhHistoryFlags,
    //! Remove functionality that makes mdrun builds non-reproducible.
    RemoveBuildMachineInformation,
    //! Allow using COM of previous step as pull group PBC reference.
    ComPrevStepAsPullGroupReference,
    //! Added possibility to output average pull force and position.
    PullAverage,
    //! Added checkpointing for MDModules.
    MDModules,
    //! Added checkpointing for modular simulator.
    ModularSimulator,
    //! Added local (per walker) weight contribution to each point in AWH.
    AwhLocalWeightSum,
    //! The total number of checkpoint versions.
    Count,
    //! Current version
    CurrentVersion = Count - 1
};


/*!
 * \brief
 * Header explaining the context of a checkpoint file.
 *
 * TODO Expand this into being a container of all data for
 * serialization of a checkpoint, which can be stored by the caller
 * (e.g. so that mdrun doesn't have to open the checkpoint twice).
 * This will separate issues of allocation from those of
 * serialization, help separate comparison from reading, and have
 * better defined transformation functions to/from trajectory frame
 * data structures.
 *
 * Several fields were once written to checkpoint file headers, but
 * have been removed. So that old files can continue to be read,
 * the names of such fields contain the string "_UNUSED" so that it
 * is clear they should not be used.
 */
struct CheckpointHeaderContents
{
    //! Version of checkpoint file read from disk.
    CheckPointVersion file_version;
    //! Version string.
    char version[CPTSTRLEN];
    //! Deprecated string for time.
    char btime_UNUSED[CPTSTRLEN];
    //! Deprecated string for user.
    char buser_UNUSED[CPTSTRLEN];
    //! Deprecated string for host.
    char bhost_UNUSED[CPTSTRLEN];
    //! Value for precision.
    int double_prec;
    //! Program string.
    char fprog[CPTSTRLEN];
    //! Time string.
    char ftime[CPTSTRLEN];
    //! Which integrator is in use.
    IntegrationAlgorithm eIntegrator;
    //! Which part of the simulation this is.
    int simulation_part;
    //! Which step the checkpoint is at.
    int64_t step;
    //! Current simulation time.
    double t;
    //! Number of nodes used for simulation,
    int nnodes;
    //! Domain decomposition settings?
    ivec dd_nc;
    //! Number of separate PME ranks.
    int npme;
    //! Number of atoms.
    int natoms;
    //! Number of temperature coupling groups.
    int ngtc;
    //! Number of Nose-Hoover pressure coupling chains.
    int nnhpres;
    //! Length of Nose-Hoover chains.
    int nhchainlength;
    //! Current FEP lambda state.
    int nlambda;
    //! Current state flags.
    int flags_state;
    //! Flags for kinetic energy.
    int flags_eks;
    //! Flags for energy history.
    int flags_enh;
    //! Flags for pull history.
    int flagsPullHistory;
    //! Flags for mystery history.
    int flags_dfh;
    //! Flags for AWH history.
    int flags_awhh;
    //! Essential dynamics states.
    int nED;
    //! Enum for coordinate swapping.
    SwapType eSwapCoords;
    //! Whether the checkpoint was written by modular simulator.
    bool isModularSimulatorCheckpoint = false;
};

/*! \brief Low-level checkpoint writing function */
void write_checkpoint_data(t_fileio*                         fp,
                           CheckpointHeaderContents          headerContents,
                           gmx_bool                          bExpanded,
                           LambdaWeightCalculation           elamstats,
                           t_state*                          state,
                           ObservablesHistory*               observablesHistory,
                           const gmx::MDModulesNotifiers&    mdModulesNotifiers,
                           std::vector<gmx_file_position_t>* outputfiles,
                           gmx::WriteCheckpointDataHolder*   modularSimulatorCheckpointData);

/* Loads a checkpoint from fn for run continuation.
 * Generates a fatal error on system size mismatch.
 * The main node reads the file
 * and communicates all the modified number of steps,
 * but not the state itself.
 * With reproducibilityRequested warns about version, build, #ranks differences.
 */
void load_checkpoint(const std::filesystem::path&   fn,
                     t_fileio*                      logfio,
                     const t_commrec*               cr,
                     const ivec                     dd_nc,
                     t_inputrec*                    ir,
                     t_state*                       state,
                     ObservablesHistory*            observablesHistory,
                     gmx_bool                       reproducibilityRequested,
                     const gmx::MDModulesNotifiers& mdModulesNotifiers,
                     gmx::ReadCheckpointDataHolder* modularSimulatorCheckpointData,
                     bool                           useModularSimulator);

/* Read everything that can be stored in t_trxframe from a checkpoint file */
void read_checkpoint_trxframe(struct t_fileio* fp, t_trxframe* fr);

/* Print the complete contents of checkpoint file fn to out */
void list_checkpoint(const std::filesystem::path& fn, FILE* out);

/*!\brief Read simulation step and part from a checkpoint file
 *
 * Used by tune_pme to handle tuning with a checkpoint file as part of the input.
 *
 * \param[in]  filename         Name of checkpoint file
 * \param[out] simulation_part  The part of the simulation that wrote the checkpoint
 * \param[out] step             The final step number of the simulation that wrote the checkpoint
 *
 * The output variables will both contain 0 if filename is NULL, the file
 * does not exist, or is not readable. */
void read_checkpoint_part_and_step(const std::filesystem::path& filename, int* simulation_part, int64_t* step);

/*!\brief Return header information from an open checkpoint file.
 *
 * Used by mdrun to handle restarts
 *
 * \param[in]  fp               Handle to open checkpoint file
 * \param[out] outputfiles      Container of output file names from the previous run. */
CheckpointHeaderContents
read_checkpoint_simulation_part_and_filenames(t_fileio* fp, std::vector<gmx_file_position_t>* outputfiles);

#endif
