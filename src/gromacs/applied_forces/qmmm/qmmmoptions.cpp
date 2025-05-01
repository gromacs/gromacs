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
 * Implements QMMMOptions class
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "qmmmoptions.h"

#include <map>

#include "gromacs/fileio/warninp.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "qmmm.h"
#include "qmmminputgenerator.h"
#include "qmmmtopologypreprocessor.h"

namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(QMMMModuleInfo::sc_name);
}

/*! \brief Following Tags denotes names of parameters from .mdp file
 * \note Changing this strings will break .tpr backwards compatibility
 */
//! \{
const std::string c_activeTag_              = "active";
const std::string c_qmGroupTag_             = "qmgroup";
const std::string c_qmChargeTag_            = "qmcharge";
const std::string c_qmMultTag_              = "qmmultiplicity";
const std::string c_qmMethodTag_            = "qmmethod";
const std::string c_qmUserInputFileNameTag_ = "qmfilenames";
//! \}

/*! \brief This tags for parameters which will be generated during grompp
 * and stored into *.tpr file via KVT
 */
//! \{
const std::string c_atomNumbersTag_ = "atomnumbers";
const std::string c_mmGroupTag_     = "mmgroup";
const std::string c_qmLinkTag_      = "qmlink";
const std::string c_mmLinkTag_      = "mmlink";
const std::string c_qmInputTag_     = "qminput";
const std::string c_qmPdbTag_       = "qmpdb";
const std::string c_qmBoxTag_       = "qmbox";
const std::string c_qmTransTag_     = "qmtrans";
//! \}

} // namespace

void QMMMOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, QMMMModuleInfo::sc_name, c_activeTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, QMMMModuleInfo::sc_name, c_qmGroupTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, QMMMModuleInfo::sc_name, c_qmMethodTag_);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, QMMMModuleInfo::sc_name, c_qmChargeTag_);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, QMMMModuleInfo::sc_name, c_qmMultTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, QMMMModuleInfo::sc_name, c_qmUserInputFileNameTag_);
}

void QMMMOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    addMdpOutputComment(builder, QMMMModuleInfo::sc_name, "empty-line", "");
    // Active flag
    addMdpOutputComment(builder, QMMMModuleInfo::sc_name, "module", "; QM/MM with CP2K");
    addMdpOutputValue(builder, QMMMModuleInfo::sc_name, c_activeTag_, parameters_.active_);

    if (parameters_.active_)
    {
        // Index group for QM atoms, default System
        addMdpOutputComment(
                builder, QMMMModuleInfo::sc_name, c_qmGroupTag_, "; Index group with QM atoms");
        addMdpOutputValue(builder, QMMMModuleInfo::sc_name, c_qmGroupTag_, groupString_);

        // QM method (DFT functional), default PBE
        addMdpOutputComment(builder,
                            QMMMModuleInfo::sc_name,
                            c_qmMethodTag_,
                            "; DFT functional for QM calculations");
        addMdpOutputValue<std::string>(
                builder, QMMMModuleInfo::sc_name, c_qmMethodTag_, c_qmmmQMMethodNames[parameters_.qmMethod_]);

        // QM charge, default 0
        addMdpOutputComment(builder, QMMMModuleInfo::sc_name, c_qmChargeTag_, "; QM charge");
        addMdpOutputValue(builder, QMMMModuleInfo::sc_name, c_qmChargeTag_, parameters_.qmCharge_);

        // QM mutiplicity, default 1
        addMdpOutputComment(builder, QMMMModuleInfo::sc_name, c_qmMultTag_, "; QM multiplicity");
        addMdpOutputValue(builder, QMMMModuleInfo::sc_name, c_qmMultTag_, parameters_.qmMultiplicity_);

        // QM input filename, default empty (will be deduced from *.tpr name during mdrun)
        addMdpOutputComment(builder,
                            QMMMModuleInfo::sc_name,
                            c_qmUserInputFileNameTag_,
                            "; Names of CP2K files during simulation");
        addMdpOutputValue(
                builder, QMMMModuleInfo::sc_name, c_qmUserInputFileNameTag_, parameters_.qmFileNameBase_);
    }
}

void QMMMOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));

    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&parameters_.active_));
    section.addOption(StringOption(c_qmGroupTag_.c_str()).store(&groupString_));
    section.addOption(EnumOption<QMMMQMMethod>(c_qmMethodTag_.c_str())
                              .enumValue(c_qmmmQMMethodNames)
                              .store(&parameters_.qmMethod_));
    section.addOption(StringOption(c_qmUserInputFileNameTag_.c_str()).store(&parameters_.qmFileNameBase_));
    section.addOption(IntegerOption(c_qmChargeTag_.c_str()).store(&parameters_.qmCharge_));
    section.addOption(IntegerOption(c_qmMultTag_.c_str()).store(&parameters_.qmMultiplicity_));
}

bool QMMMOptions::active() const
{
    return parameters_.active_;
}

const QMMMParameters& QMMMOptions::parameters()
{
    return parameters_;
}

void QMMMOptions::setLogger(const MDLogger& logger)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    logger_ = &logger;
}

const MDLogger& QMMMOptions::logger() const
{
    GMX_RELEASE_ASSERT(logger_, "Logger not set for QMMMOptions.");
    return *logger_;
}

void QMMMOptions::setWarninp(WarningHandler* wi)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    wi_ = wi;
}

void QMMMOptions::appendLog(const std::string& msg)
{
    GMX_LOG(logger().info).asParagraph().appendText(msg);
}

void QMMMOptions::appendWarning(const std::string& msg)
{
    if (wi_)
    {
        wi_->addWarning(msg);
    }
}

void QMMMOptions::setQMMMGroupIndices(const IndexGroupsAndNames& indexGroupsAndNames)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    // Create QM index
    parameters_.qmIndices_ = indexGroupsAndNames.indices(groupString_);

    // Check that group is not empty
    if (parameters_.qmIndices_.empty())
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Group %s defining QM atoms should not be empty.", groupString_.c_str())));
    }

    // Create temporary index for the whole System
    auto systemIndices = indexGroupsAndNames.indices(std::string("System"));

    // Sort qmindices_ and sindices_
    std::sort(parameters_.qmIndices_.begin(), parameters_.qmIndices_.end());
    std::sort(systemIndices.begin(), systemIndices.end());

    // Create MM index
    parameters_.mmIndices_.reserve(systemIndices.size());

    // Position in qmindicies_
    size_t j = 0;
    // Now loop over sindices_ and write to mmindices_ only the atoms which does not belong to qmindices_
    for (size_t i = 0; i < systemIndices.size(); i++)
    {
        if (systemIndices[i] != parameters_.qmIndices_[j])
        {
            parameters_.mmIndices_.push_back(systemIndices[i]);
        }
        else
        {
            if (j < parameters_.qmIndices_.size() - 1)
            {
                j++;
            }
        }
    }
}

void QMMMOptions::processTprFilename(const MdRunInputFilename& tprFilename)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    // Provided name should not be empty
    GMX_RELEASE_ASSERT(!tprFilename.mdRunFilename_.empty(),
                       "Filename of the *.tpr simulation file is empty");

    // Exit if parameters_.qmFileNameBase_ has been provided
    if (!parameters_.qmFileNameBase_.empty())
    {
        return;
    }

    parameters_.qmFileNameBase_ =
            stripExtension(std::filesystem::path(tprFilename.mdRunFilename_).filename())
                    .concat("_cp2k")
                    .string();
}

void QMMMOptions::processExternalInputFile()
{
    // First check if we could read qmExternalInputFileName_
    TextReader fInp(qmExternalInputFileName_);

    // Then we need to build a map with all CP2K parameters found in file
    std::map<std::string, std::string> cp2kParams;

    // Key-value pair of the parameters
    std::pair<std::string, std::string> kv;

    // Loop over all lines in the file
    std::string line;
    while (fInp.readLine(&line))
    {
        // Split line into words
        auto words = splitString(line);

        // If we have 2 or more words in the line then build key-value pair
        if (words.size() >= 2)
        {
            kv.first  = words[0];
            kv.second = words[1];

            // Convert kv.first to the upper case
            kv.first = toUpperCase(kv.first);

            // Put into the map
            cp2kParams.insert(kv);
        }
    }
    fInp.close();

    // Check if @INCLUDE found in file
    if (cp2kParams.count("@INCLUDE") > 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "@INCLUDE directive is not allowed but found in external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if CHARGE found in file
    if (cp2kParams.count("CHARGE") == 0)
    {
        GMX_THROW(InconsistentInputError(
                formatString("Parameter CHARGE not found in the external CP2K input file %s",
                             qmExternalInputFileName_.c_str())));
    }

    // Check if MULTIPLICITY found in file
    if (cp2kParams.count("MULTIPLICITY") == 0)
    {
        GMX_THROW(InconsistentInputError(
                formatString("Parameter MULTIPLICITY not found in the external CP2K input file %s",
                             qmExternalInputFileName_.c_str())));
    }

    // Check if RUN_TYPE in the file present
    if (cp2kParams.count("RUN_TYPE") == 0)
    {
        GMX_THROW(InconsistentInputError(
                formatString("Parameter RUN_TYPE not found in the external CP2K input file %s",
                             qmExternalInputFileName_.c_str())));
    }

    // Check if RUN_TYPE in the file is equal to ENERGY_FORCE
    if (toUpperCase(cp2kParams["RUN_TYPE"]) != "ENERGY_FORCE")
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter RUN_TYPE should be ENERGY_FORCE in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if COORD_FILE_FORMAT in the file present
    if (cp2kParams.count("COORD_FILE_FORMAT") == 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter COORD_FILE_FORMAT not found in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if COORD_FILE_FORMAT in the file is equal to PDB
    if (toUpperCase(cp2kParams["COORD_FILE_FORMAT"]) != "PDB")
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter COORD_FILE_FORMAT must be PDB in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if CHARGE_EXTENDED in the file present
    if (cp2kParams.count("CHARGE_EXTENDED") == 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter CHARGE_EXTENDED not found in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if CHARGE_EXTENDED in the file is equal to TRUE
    if (toUpperCase(cp2kParams["CHARGE_EXTENDED"]) != "TRUE")
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter CHARGE_EXTENDED must be TRUE in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Check if COORD_FILE_NAME found in file
    if (cp2kParams.count("COORD_FILE_NAME") == 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
                "Parameter COORD_FILE_NAME not found in the external CP2K input file %s",
                qmExternalInputFileName_.c_str())));
    }

    // Set parameters_.qmCharge_ and qmMult_
    parameters_.qmCharge_       = fromStdString<int>(cp2kParams["CHARGE"]);
    parameters_.qmMultiplicity_ = fromStdString<int>(cp2kParams["MULTIPLICITY"]);

    // Read the whole input as one string and replace COORD_FILE_NAME value with %s placeholder
    std::string str      = TextReader::readFileToString(qmExternalInputFileName_);
    parameters_.qmInput_ = replaceAllWords(str, cp2kParams["COORD_FILE_NAME"], "%s");
}

void QMMMOptions::setQMExternalInputFile(const QMInputFileName& qmExternalInputFileName)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    if (parameters_.qmMethod_ != QMMMQMMethod::INPUT)
    {
        if (qmExternalInputFileName.hasQMInputFileName_)
        {
            // If parameters_.qmMethod_ != INPUT then user should not provide external input file
            GMX_THROW(InconsistentInputError(
                    "External CP2K input file has been provided with -qmi option, but "
                    + moduleName() + "-" + c_qmMethodTag_ + " is not INPUT"));
        }

        // Exit if we dont need to process external input file
        return;
    }

    // Case where user should provide external input file with -qmi option
    if (parameters_.qmMethod_ == QMMMQMMethod::INPUT && !qmExternalInputFileName.hasQMInputFileName_)
    {
        GMX_THROW(InconsistentInputError(moduleName() + "-" + c_qmMethodTag_
                                         + " = INPUT requested, but external CP2K "
                                           "input file is not provided with -qmi option"));
    }

    // If external input is provided by the user then we should process it and save into the parameters_
    qmExternalInputFileName_ = qmExternalInputFileName.qmInputFileName_;
    processExternalInputFile();
}

void QMMMOptions::processCoordinates(const CoordinatesAndBoxPreprocessed& coord)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    QMMMInputGenerator inpGen(
            parameters_, coord.pbc_, coord.box_, atomCharges_, coord.coordinates_.unpaddedConstArrayRef());

    // Generate pdb file with point charges for CP2K
    parameters_.qmPdb_ = inpGen.generateCP2KPdb();

    // In case parameters_.qmMethod_ != INPUT we should generate CP2K Input, QM box and translation
    if (parameters_.qmMethod_ != QMMMQMMethod::INPUT)
    {
        /* Check if some of the box vectors dimension lower that 1 nm.
         * For SCF stability box should be big enough.
         */
        matrix box;
        copy_mat(coord.box_, box);
        if (norm(box[0]) < 1.0 || norm(box[1]) < 1.0 || norm(box[2]) < 1.0)
        {
            GMX_THROW(InconsistentInputError(
                    "One of the box vectors is shorter than 1 nm.\n"
                    "For stable CP2K SCF convergence all simulation box vectors should be "
                    ">= 1 nm. Please consider to increase simulation box or provide custom CP2K "
                    "input using "
                    + moduleName() + "-" + c_qmMethodTag_ + " = INPUT"));
        }

        parameters_.qmInput_ = inpGen.generateCP2KInput();
        copy_mat(inpGen.qmBox(), parameters_.qmBox_);
        parameters_.qmTrans_ = inpGen.qmTrans();
    }
}

void QMMMOptions::modifyQMMMTopology(gmx_mtop_t* mtop)
{
    // Exit if QMMM module is not active
    if (!parameters_.active_)
    {
        return;
    }

    // Process topology
    QMMMTopologyPreprocessor topPrep(parameters_.qmIndices_);
    topPrep.preprocess(mtop);

    // Get atom numbers
    parameters_.atomNumbers_ = copyOf(topPrep.atomNumbers());

    // Get atom point charges
    atomCharges_ = copyOf(topPrep.atomCharges());

    // Get Link Frontier
    parameters_.link_ = copyOf(topPrep.linkFrontier());

    // Get info about modifications
    QMMMTopologyInfo topInfo = topPrep.topInfo();

    // Cast int qmCharge_ to real as further calculations use floating point
    real qmC = static_cast<real>(parameters_.qmCharge_);

    // Print message to the log about performed modifications
    std::string msg = "\nQMMM Interface with CP2K is active, topology was modified!\n";

    msg += formatString(
            "Number of QM atoms: %d\nNumber of MM atoms: %d\n", topInfo.numQMAtoms, topInfo.numMMAtoms);

    msg += formatString("Total charge of the classical system (before modifications): %.5f\n",
                        topInfo.remainingMMCharge + topInfo.totalClassicalChargeOfQMAtoms);

    msg += formatString("Classical charge removed from QM atoms: %.5f\n",
                        topInfo.totalClassicalChargeOfQMAtoms);

    if (topInfo.numVirtualSitesModified > 0)
    {
        msg += formatString(
                "Note: There are %d virtual sites found, which are built from QM atoms only. "
                "Classical charges on them have been removed as well.\n",
                topInfo.numVirtualSitesModified);
    }

    msg += formatString("Total charge of QMMM system (after modifications): %.5f\n",
                        qmC + topInfo.remainingMMCharge);

    if (topInfo.numBondsRemoved > 0)
    {
        msg += formatString("Bonds removed: %d\n", topInfo.numBondsRemoved);
    }

    if (topInfo.numAnglesRemoved > 0)
    {
        msg += formatString("Angles removed: %d\n", topInfo.numAnglesRemoved);
    }

    if (topInfo.numDihedralsRemoved > 0)
    {
        msg += formatString("Dihedrals removed: %d\n", topInfo.numDihedralsRemoved);
    }

    if (topInfo.numSettleRemoved > 0)
    {
        msg += formatString("Settles removed: %d\n", topInfo.numSettleRemoved);
    }

    if (topInfo.numConnBondsAdded > 0)
    {
        msg += formatString("F_CONNBONDS (type 5 bonds) added: %d\n", topInfo.numConnBondsAdded);
    }

    if (topInfo.numLinkBonds > 0)
    {
        msg += formatString("QM-MM broken bonds found: %d\n", topInfo.numLinkBonds);
    }

    appendLog(msg + "\n");

    /* We should warn the user if there is inconsistence between removed classical charges
     * on QM atoms and total QM charge
     */
    if (std::abs(topInfo.totalClassicalChargeOfQMAtoms - qmC) > 1E-5)
    {
        msg = formatString(
                "Total charge of your QMMM system differs from classical system! "
                "Consider manually spreading %.5lf charge over MM atoms nearby to the QM "
                "region\n",
                topInfo.totalClassicalChargeOfQMAtoms - qmC);
        appendWarning(msg);
    }

    // If there are many constrained bonds in QM system then we should also warn the user
    if (topInfo.numConstrainedBondsInQMSubsystem > 2)
    {
        msg = "Your QM subsystem has a lot of constrained bonds. They probably have been "
              "generated automatically. That could produce an artifacts in the simulation. "
              "Consider constraints = none in the mdp file.\n";
        appendWarning(msg);
    }
}

void QMMMOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{
    // Write QM atoms index
    auto GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_qmGroupTag_);
    for (const auto& indexValue : parameters_.qmIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write MM atoms index
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_mmGroupTag_);
    for (const auto& indexValue : parameters_.mmIndices_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write atoms numbers
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_atomNumbersTag_);
    for (const auto& indexValue : parameters_.atomNumbers_)
    {
        GroupIndexAdder.addValue(indexValue);
    }

    // Write link
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_qmLinkTag_);
    for (const auto& indexValue : parameters_.link_)
    {
        GroupIndexAdder.addValue(indexValue.qm);
    }
    GroupIndexAdder = treeBuilder.addUniformArray<std::int64_t>(moduleName() + "-" + c_mmLinkTag_);
    for (const auto& indexValue : parameters_.link_)
    {
        GroupIndexAdder.addValue(indexValue.mm);
    }

    // Write CP2K input file content
    treeBuilder.addValue<std::string>(moduleName() + "-" + c_qmInputTag_, parameters_.qmInput_);

    // Write CP2K pdb file content
    treeBuilder.addValue<std::string>(moduleName() + "-" + c_qmPdbTag_, parameters_.qmPdb_);

    // Write QM box matrix
    auto DoubleArrayAdder = treeBuilder.addUniformArray<double>(moduleName() + "-" + c_qmBoxTag_);
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            DoubleArrayAdder.addValue(static_cast<double>(parameters_.qmBox_[i][j]));
        }
    }

    // Write QM Translation vector
    DoubleArrayAdder = treeBuilder.addUniformArray<double>(moduleName() + "-" + c_qmTransTag_);
    for (int i = 0; i < DIM; i++)
    {
        DoubleArrayAdder.addValue(static_cast<double>(parameters_.qmTrans_[i]));
    }
}

void QMMMOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{
    // Check if active
    if (!parameters_.active_)
    {
        return;
    }

    // Try to read QM atoms index
    if (!tree.keyExists(moduleName() + "-" + c_qmGroupTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find QM atoms index vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    auto kvtIndexArray = tree[moduleName() + "-" + c_qmGroupTag_].asArray().values();
    parameters_.qmIndices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(parameters_.qmIndices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    // Try to read MM atoms index
    if (!tree.keyExists(moduleName() + "-" + c_mmGroupTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find MM atoms index vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_mmGroupTag_].asArray().values();
    parameters_.mmIndices_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(parameters_.mmIndices_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    // Try to read atoms numbers
    if (!tree.keyExists(moduleName() + "-" + c_atomNumbersTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find Atom Numbers vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_atomNumbersTag_].asArray().values();
    parameters_.atomNumbers_.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(parameters_.atomNumbers_),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    // Try to read Link Frontier (two separate vectors and then combine)
    std::vector<Index> qmLink;
    std::vector<Index> mmLink;

    if (!tree.keyExists(moduleName() + "-" + c_qmLinkTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find QM Link Frontier vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_qmLinkTag_].asArray().values();
    qmLink.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(qmLink),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    if (!tree.keyExists(moduleName() + "-" + c_mmLinkTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find MM Link Frontier vector required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    kvtIndexArray = tree[moduleName() + "-" + c_mmLinkTag_].asArray().values();
    mmLink.resize(kvtIndexArray.size());
    std::transform(std::begin(kvtIndexArray),
                   std::end(kvtIndexArray),
                   std::begin(mmLink),
                   [](const KeyValueTreeValue& val) { return val.cast<std::int64_t>(); });

    parameters_.link_.resize(qmLink.size());
    for (size_t i = 0; i < qmLink.size(); i++)
    {
        parameters_.link_[i].qm = qmLink[i];
        parameters_.link_[i].mm = mmLink[i];
    }

    // Try to read CP2K input and pdb strings from *.tpr
    if (!tree.keyExists(moduleName() + "-" + c_qmInputTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find CP2K input string required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    parameters_.qmInput_ = tree[moduleName() + "-" + c_qmInputTag_].cast<std::string>();

    if (!tree.keyExists(moduleName() + "-" + c_qmPdbTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find CP2K pdb string required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    parameters_.qmPdb_ = tree[moduleName() + "-" + c_qmPdbTag_].cast<std::string>();

    // Try to read QM box
    if (!tree.keyExists(moduleName() + "-" + c_qmBoxTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find QM box matrix required for QM/MM simulation.\nThis could be "
                "caused by incompatible or corrupted tpr input file."));
    }
    auto kvtDoubleArray = tree[moduleName() + "-" + c_qmBoxTag_].asArray().values();
    for (int i = 0; i < DIM; i++)
    {
        for (int j = 0; j < DIM; j++)
        {
            parameters_.qmBox_[i][j] = static_cast<real>(kvtDoubleArray[i * 3 + j].cast<double>());
        }
    }

    // Try to read QM translation vector
    if (!tree.keyExists(moduleName() + "-" + c_qmTransTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find QM subsystem centering information for QM/MM simulation.\nThis could "
                "be caused by incompatible or corrupted tpr input file."));
    }
    kvtDoubleArray = tree[moduleName() + "-" + c_qmTransTag_].asArray().values();
    for (int i = 0; i < DIM; i++)
    {
        parameters_.qmTrans_[i] = static_cast<real>(kvtDoubleArray[i].cast<double>());
    }
}

} // namespace gmx
