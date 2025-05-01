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
 * Declares options for QM/MM
 * QMMMOptions class responsible for all parameters set up during pre-processing
 * also modificatios of topology would be done here
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMOPTIONS_H
#define GMX_APPLIED_FORCES_QMMMOPTIONS_H

#include "gromacs/mdtypes/imdpoptionprovider.h"

#include "qmmmtypes.h"

struct gmx_mtop_t;
class WarningHandler;

namespace gmx
{

class IndexGroupsAndNames;
class KeyValueTreeObject;
class KeyValueTreeBuilder;
class MDLogger;
struct MdRunInputFilename;
struct CoordinatesAndBoxPreprocessed;
struct QMInputFileName;

/*! \internal
 * \brief Input data storage for QM/MM
 */
class QMMMOptions final : public IMdpOptionProvider
{
public:
    //! Implementation of IMdpOptionProvider method
    void initMdpTransform(IKeyValueTreeTransformRules* rules) override;

    /*! \brief
     * Build mdp parameters for QMMM to be output after pre-processing.
     * \param[in, out] builder the builder for the mdp options output KVT.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    /*! \brief
     * Connects options names and data.
     */
    void initMdpOptions(IOptionsContainerWithSections* options) override;

    //! Report if this set of MDP options is active (i.e. QMMM MdModule is active)
    bool active() const;

    //! Get parameters_ instance
    const QMMMParameters& parameters();

    /*! \brief Evaluate and store atom indices.
     * During pre-processing, use the group string from the options to
     * evaluate the indices of both QM atoms and MM atoms, also stores them
     * as vectors into the parameters_
     * \param[in] indexGroupsAndNames object containing data about index groups and names
     */
    void setQMMMGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames);

    /*! \brief Process external QM input file in case it is provided with -qmi option of grompp.
     * Produces parameters_.qmInput in case parameters_.qmMethod_ = INPUT
     * \param[in] qmExternalInputFileName structure with information about external QM input
     */
    void setQMExternalInputFile(const QMInputFileName& qmExternalInputFileName);

    /*! \brief Process coordinates, PbcType and Box in order to produce CP2K sample input.
     * Produces qmPdb_ in all cases. Produces parameters_.qmInput_, parameters_.qmTrans_
     * and parameters_.qmBox_ in case parameters_.qmMethod_ != INPUT.
     * \param[in] coord structure with coordinates and box dimensions
     */
    void processCoordinates(const CoordinatesAndBoxPreprocessed& coord);

    /*! \brief Modifies topology in case of active QMMM module using QMMMTopologyPreprocessor
     * \param[in,out] mtop topology to modify for QMMM
     */
    void modifyQMMMTopology(gmx_mtop_t* mtop);

    //! Store the paramers that are not mdp options in the tpr file
    void writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder);

    //! Set the internal parameters that are stored in the tpr file
    void readInternalParametersFromKvt(const KeyValueTreeObject& tree);

    /*! \brief Process MdRunInputFilename notification during mdrun.
     * In case parameters_.qmFileNameBase_ is empty sets it to tpr name with _cp2k suffix
     * \param[in] tprFilename name of the *.tpr file that mdrun simulates
     */
    void processTprFilename(const MdRunInputFilename& tprFilename);

    //! Set the MDLogger instance
    void setLogger(const MDLogger& logger);

    //! Get the logger instance
    const MDLogger& logger() const;

    //! Set the warninp instance
    void setWarninp(WarningHandler* wi);

private:
    //! Write message to the log
    void appendLog(const std::string& msg);

    //! Write grompp warning
    void appendWarning(const std::string& msg);

    /*! \brief Processes external CP2K input file with the name qmExternalInputFileName_
     * 1) Extracts parameters_.qmCharge_ and parameters_.qmMult_ from it
     * 2) Replaces COORD_FILE_NAME parameter in it with a placeholder
     * 3) Saves it into the parameters_.qmInput_
     * \throws FileIOError if input file could not be read
     * \throws InvalidInputError if MULTIPLICITY, CHARGE or COORD_FILE_NAME not found
     */
    void processExternalInputFile();

    /*! \brief MDLogger during preprocessing
     *
     * This is a pointer only because we need an "optional reference"
     * to a const MDLogger before the notification always provides the
     * actual reference. */
    const MDLogger* logger_ = nullptr;

    //! Instance of warning bookkeeper
    WarningHandler* wi_ = nullptr;

    //! QM index group name, Default whole System
    std::string groupString_ = "System";

    //! QMMM parameters built from mdp input
    QMMMParameters parameters_;

    //! Name of the external input file provided with -qmi option of grompp
    std::string qmExternalInputFileName_;

    //! Vector with atoms point charges
    std::vector<real> atomCharges_;
};

} // namespace gmx

#endif
