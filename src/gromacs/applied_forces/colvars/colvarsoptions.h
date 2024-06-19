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
 * Declares options for Colvars. This class handles parameters set during
 * pre-processing time.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_COLVARSOPTIONS_H
#define GMX_APPLIED_FORCES_COLVARSOPTIONS_H

#include <map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"

enum class PbcType : int;
struct gmx_mtop_t;


namespace gmx
{

class KeyValueTreeObject;
class KeyValueTreeBuilder;
struct CoordinatesAndBoxPreprocessed;
struct MdRunInputFilename;
class IKeyValueTreeTransformRules;
class IOptionsContainerWithSections;
class KeyValueTreeObjectBuilder;


//! Tag with name of the Colvars MDModule
static const std::string c_colvarsModuleName = "colvars";


/*! \internal
 * \brief Input data storage for colvars
 */
class ColvarsOptions final : public IMdpOptionProvider
{
public:
    //! From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for colvars to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KV-tree.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connect option name and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    /*! \brief Store the topology of the system.
     * \param[in,out] mtop topology object
     */
    void processTopology(gmx_mtop_t* mtop);

    /*! \brief Process coordinates, PbcType and Box in order to validate the colvars input.
     * \param[in] coord structure with coordinates and box dimensions
     */
    void processCoordinates(const CoordinatesAndBoxPreprocessed& coord);

    //! Set the MDLogger instance
    void setLogger(const MDLogger& logger);

    /*! \brief Process EdrOutputFilename notification during mdrun.
     * Used to set the prefix of Colvars output files based on the .edr filename
     * \param[in] filename name of the *.edr file that mdrun will produce
     */
    void processEdrFilename(const EdrOutputFilename& filename);

    /*! \brief Store the ensemble temperature of the system if available.
     * \param[in] temp temperature object
     */
    void processTemperature(const EnsembleTemperature& temp);

    //! Report if this colvars module is active
    bool isActive() const;

    //! Return the file name of the colvars config
    const std::string& colvarsFileName() const;

    //! Return the content of the colvars config file
    const std::string& colvarsConfigContent() const;

    //! Return the colvars atoms coordinates
    const std::vector<RVec>& colvarsAtomCoords() const;

    //! Return the prefix for output colvars files
    const std::string& colvarsOutputPrefix() const;

    //! Return the ensemble temperature
    const real& colvarsEnsTemp() const;

    //! Return the map of all others colvars input files
    const std::map<std::string, std::string>& colvarsInputFiles() const;

    //! Return the colvars seed
    int colvarsSeed() const;

    /*! \brief Function to set internal paramaters outside the way done
     * through the MDModule notifiers and callbacks.
     * Use exclusively in the test framework.
     *
     * \param[in] colvarsfile Name of the colvars input file.
     * \param[in] topology Atoms topology
     * \param[in] coords Coordinates of each atom in the system
     * \param[in] pbcType Periodic boundary conditions
     * \param[in] boxValues Matrix with full box of the system
     * \param[in] temperature the constant ensemble temperature
     */
    void setParameters(const std::string&   colvarsfile,
                       const t_atoms&       topology,
                       ArrayRef<const RVec> coords,
                       PbcType              pbcType,
                       const matrix         boxValues,
                       real                 temperature);


private:
    //! Indicate if colvars module is active
    bool active_ = false;

    /*! \brief Following Tags denotes names of parameters from .mdp file
     * \note Changing this strings will break .tpr backwards compability
     */
    //! \{
    const std::string c_activeTag_          = "active";
    const std::string c_colvarsFileNameTag_ = "configfile";
    const std::string c_colvarsSeedTag_     = "seed";
    //! \}


    /*! \brief This tags for parameters which will be generated during grompp
     * and stored into *.tpr file via KVT
     */
    //! \{
    const std::string c_inputStreamsTag_   = "inputStreams";
    const std::string c_configStringTag_   = "configString";
    const std::string c_startingCoordsTag_ = "startingCoords";
    const std::string c_ensTempTag_        = "ensTemp";

    //! \}

    //! Colvars config filename
    std::string colvarsFileName_;


    //! Colvars seed for Langevin integrator
    int colvarsSeed_ = -1;

    //! Content of the colvars config file
    std::string colvarsConfigString_;
    //! Topology of the system
    t_atoms gmxAtoms_;
    //! Coordinates
    ArrayRef<const RVec> x_;
    //! PBC Type
    PbcType pbc_;
    //! Box
    matrix box_;
    //! Vector with colvars atoms coordinates
    std::vector<RVec> colvarsAtomCoords_;
    //! Inputs files saved as strings inside KVT
    std::map<std::string, std::string> inputFiles_;

    real ensembleTemperature_;


    //! Logger instance
    const MDLogger* logger_ = nullptr;

    /*! \brief String containing the prefix for output colvars files
     * default value empty, means will be deduced from *.tpr name during mdrun
     */
    std::string outputPrefix_;
};

} // namespace gmx

#endif
