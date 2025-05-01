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
 * Implements options for Colvars.
 * Keep a minimal ColvarsOptions implementation in the case Colvars is not activated
 * to ensure compatbility (with gmx tools for example).
 */
#include "gmxpre.h"

#include "colvarsoptions.h"

#include "config.h"

#include <cstddef>

#include <filesystem>
#include <fstream>
#include <optional>

#include "gromacs/math/arrayrefwithpadding.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/imdpoptionprovider_helpers.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/textreader.h"

#include "colvarsMDModule.h"

#if GMX_HAVE_COLVARS
#    include "colvarspreprocessor.h"
#endif

enum class PbcType : int;
struct gmx_mtop_t;


namespace gmx
{

namespace
{

//! Helper function to make a std::string containing the module name
std::string moduleName()
{
    return std::string(ColvarsModuleInfo::sc_name);
}

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

} // namespace


void ColvarsOptions::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    const auto& stringIdentityTransform = [](std::string s) { return s; };
    addMdpTransformFromString<bool>(rules, &fromStdString<bool>, ColvarsModuleInfo::sc_name, c_activeTag_);
    addMdpTransformFromString<std::string>(
            rules, stringIdentityTransform, ColvarsModuleInfo::sc_name, c_colvarsFileNameTag_);
    addMdpTransformFromString<int>(rules, &fromStdString<int>, ColvarsModuleInfo::sc_name, c_colvarsSeedTag_);
}


void ColvarsOptions::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    // new empty line before writing colvars mdp values
    addMdpOutputComment(builder, ColvarsModuleInfo::sc_name, "empty-line", "");
    addMdpOutputComment(builder, ColvarsModuleInfo::sc_name, "module", "; Colvars bias");
    addMdpOutputValue(builder, ColvarsModuleInfo::sc_name, c_activeTag_, active_);

    if (active_)
    {
        addMdpOutputComment(
                builder, ColvarsModuleInfo::sc_name, c_colvarsFileNameTag_, "; colvars config file");
        addMdpOutputValue<std::string>(
                builder, ColvarsModuleInfo::sc_name, c_colvarsFileNameTag_, colvarsFileName_);

        addMdpOutputComment(
                builder, ColvarsModuleInfo::sc_name, c_colvarsSeedTag_, "; Colvars seed");
        addMdpOutputValue<int>(builder, ColvarsModuleInfo::sc_name, c_colvarsSeedTag_, colvarsSeed_);
    }
}


void ColvarsOptions::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection(moduleName().c_str()));
    section.addOption(BooleanOption(c_activeTag_.c_str()).store(&active_));
    section.addOption(StringOption(c_colvarsFileNameTag_.c_str()).store(&colvarsFileName_));
    section.addOption(IntegerOption(c_colvarsSeedTag_.c_str()).store(&colvarsSeed_));
}


bool ColvarsOptions::isActive() const
{
    return active_;
}

//! Get the logger instance
const MDLogger& ColvarsOptions::logger() const
{
    GMX_RELEASE_ASSERT(logger_, "Logger not set for ColvarsOptions.");
    return *logger_;
}

#if GMX_HAVE_COLVARS


void ColvarsOptions::writeInternalParametersToKvt(KeyValueTreeObjectBuilder treeBuilder)
{

    // Copy the content of the colvars input file into a string for latter save in KVT
    if (!colvarsFileName_.empty())
    {
        colvarsConfigString_ = TextReader::readFileToString(colvarsFileName_);
    }

    // Write colvars input file as a string
    treeBuilder.addValue<std::string>(moduleName() + "-" + c_configStringTag_, colvarsConfigString_);


    ColvarsPreProcessor colvarsPreProcess(
            colvarsConfigString_, gmxAtoms_, pbc_, logger(), ensembleTemperature_, colvarsSeed_, box_, x_);
    //! Vector with colvars atoms coordinates
    colvarsAtomCoords_ = colvarsPreProcess.getColvarsCoords();

    // Save other colvars input files into the KVT
    if (!colvarsPreProcess.inputStreamsToKVT(treeBuilder, moduleName() + "-" + c_inputStreamsTag_))
    {
        GMX_THROW(InternalError("Cannot save colvars input files into the tpr."));
    }

    // Write colvars atoms coords
    auto DoubleArrayAdder = treeBuilder.addUniformArray<double>(moduleName() + "-" + c_startingCoordsTag_);
    for (const auto& indexValue : colvarsAtomCoords_)
    {
        for (int j = 0; j < DIM; j++)
        {
            DoubleArrayAdder.addValue(static_cast<double>(indexValue[j]));
        }
    }

    // Write ensemble temperature
    treeBuilder.addValue<real>(moduleName() + "-" + c_ensTempTag_, ensembleTemperature_);

    // Write seed
    treeBuilder.addValue<int>(moduleName() + "-" + c_colvarsSeedTag_, colvarsSeed_);
}


void ColvarsOptions::readInternalParametersFromKvt(const KeyValueTreeObject& tree)
{

    if (!active_)
    {
        return;
    }


    // Retrieve the content of all inputfiles listed in the KVT as "colvars-inputStreams-filename"
    const std::string inputStreamsKeySubstring = moduleName() + "-" + c_inputStreamsTag_;
    for (const auto& a : tree.properties())
    {
        std::size_t pos = a.key().find(inputStreamsKeySubstring);
        if (pos != std::string::npos)
        {
            const std::string filename = a.key().substr(pos + inputStreamsKeySubstring.size() + 1);
            inputFiles_[filename]      = tree[a.key()].cast<std::string>();
        }
    }

    if (!tree.keyExists(moduleName() + "-" + c_configStringTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find colvars-configString required for colvars simulation."));
    }
    colvarsConfigString_ = tree[moduleName() + "-" + c_configStringTag_].cast<std::string>();


    if (!tree.keyExists(moduleName() + "-" + c_startingCoordsTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find colvars-startingCoords required for colvars simulation."));
    }

    auto kvtDoubleArray = tree[moduleName() + "-" + c_startingCoordsTag_].asArray().values();


    // Make sure the coordinates saved are consistent with the dimensions
    if (kvtDoubleArray.size() % DIM != 0)
    {
        GMX_THROW(InconsistentInputError(
                "Coordinates saved in colvars-startingCoords are in the wrong format."));
    }

    for (size_t i = 0; i < kvtDoubleArray.size() / DIM; i++)
    {
        RVec x;
        for (int j = 0; j < DIM; j++)
        {
            x[j] = static_cast<real>(kvtDoubleArray[i * DIM + j].cast<double>());
        }
        colvarsAtomCoords_.push_back(x);
    }

    if (!tree.keyExists(moduleName() + "-" + c_ensTempTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find ensemble temperature required for colvars simulation."));
    }
    ensembleTemperature_ = tree[moduleName() + "-" + c_ensTempTag_].cast<real>();


    if (!tree.keyExists(moduleName() + "-" + c_colvarsSeedTag_))
    {
        GMX_THROW(InconsistentInputError(
                "Cannot find colvars seed required for colvars simulation."));
    }
    colvarsSeed_ = tree[moduleName() + "-" + c_colvarsSeedTag_].cast<int>();
}


void ColvarsOptions::processTopology(gmx_mtop_t* mtop)
{
    gmxAtoms_ = gmx_mtop_global_atoms(*mtop);
}


void ColvarsOptions::processCoordinates(const CoordinatesAndBoxPreprocessed& coord)
{

    x_   = coord.coordinates_.unpaddedConstArrayRef();
    pbc_ = coord.pbc_;
    copy_mat(coord.box_, box_);
}

void ColvarsOptions::setLogger(const MDLogger& logger)
{
    logger_ = &logger;
}

void ColvarsOptions::processEdrFilename(const EdrOutputFilename& filename)
{
    // Do nothing if Colvars is not active
    if (!active_)
    {
        return;
    }

    // Provided name should not be empty
    GMX_RELEASE_ASSERT(!filename.edrOutputFilename_.empty(), "Empty name for the *.edr output file");

    outputPrefix_ = stripExtension(std::filesystem::path(filename.edrOutputFilename_)).string();
}


void ColvarsOptions::processTemperature(const EnsembleTemperature& temp)
{
    if (temp.constantEnsembleTemperature_)
    {
        ensembleTemperature_ = temp.constantEnsembleTemperature_.value();
    }
    else
    {
        ensembleTemperature_ = -1;
    }
}

const std::string& ColvarsOptions::colvarsFileName() const
{
    return colvarsFileName_;
}


const std::string& ColvarsOptions::colvarsConfigContent() const
{
    return colvarsConfigString_;
}

const std::vector<RVec>& ColvarsOptions::colvarsAtomCoords() const
{
    return colvarsAtomCoords_;
}

const std::string& ColvarsOptions::colvarsOutputPrefix() const
{
    return outputPrefix_;
}

const real& ColvarsOptions::colvarsEnsTemp() const
{
    return ensembleTemperature_;
}

const std::map<std::string, std::string>& ColvarsOptions::colvarsInputFiles() const
{
    return inputFiles_;
}

int ColvarsOptions::colvarsSeed() const
{
    return colvarsSeed_;
}

void ColvarsOptions::setParameters(const std::string&   colvarsfile,
                                   const t_atoms&       topology,
                                   ArrayRef<const RVec> coords,
                                   PbcType              pbcType,
                                   const matrix         boxValues,
                                   real                 temperature)
{
    colvarsFileName_ = colvarsfile;
    gmxAtoms_        = topology;
    x_               = coords;
    pbc_             = pbcType;
    copy_mat(boxValues, box_);
    ensembleTemperature_ = temperature;
}


#endif // GMX_HAVE_COLVARS

} // namespace gmx
