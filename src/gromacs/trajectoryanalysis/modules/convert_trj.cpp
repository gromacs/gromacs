/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
 * Implements gmx::analysismodules::ConvertTrj.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "convert_trj.h"

#include <algorithm>
#include <memory>
#include <string>

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/requirements.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"

struct t_pbc;
struct t_trxframe;

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * ConvertTrj
 */

class ConvertTrj : public TrajectoryAnalysisModule
{
public:
    ConvertTrj();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    TrajectoryFrameWriterPointer    output_;
    Selection                       sel_;
    std::string                     name_;
    OutputRequirementOptionDirector requirementsBuilder_;
};

ConvertTrj::ConvertTrj() {}


void ConvertTrj::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] converts trajectory files between different formats.",
        "The module supports writing all GROMACS supported file formats from",
        "the supported input formats.",
        "",
        "Included is also a selection of possible options to modify individual",
        "trajectory frames, including options to produce slimmer",
        "output files. It is also possible to replace the particle information stored",
        "in the input trajectory with those from a structure file",
        "",
        "The module can also generate subsets of trajectories based on user supplied",
        "selections.",
    };

    options->addOption(SelectionOption("select").store(&sel_).onlyAtoms().description(
            "Selection of particles to write to the file"));

    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Trajectory)
                               .outputFile()
                               .store(&name_)
                               .defaultBasename("trajout")
                               .required()
                               .description("Output trajectory"));

    requirementsBuilder_.initOptions(options);

    settings->setHelpText(desc);
}

void ConvertTrj::optionsFinished(TrajectoryAnalysisSettings* settings)
{
    int frameFlags = TRX_NEED_X;

    frameFlags |= TRX_READ_V;
    frameFlags |= TRX_READ_F;

    settings->setFrameFlags(frameFlags);
}


void ConvertTrj::initAnalysis(const TrajectoryAnalysisSettings& /*settings*/, const TopologyInformation& top)
{
    output_ = createTrajectoryFrameWriter(top.mtop(),
                                          sel_,
                                          name_,
                                          top.hasTopology() ? top.copyAtoms() : nullptr,
                                          requirementsBuilder_.process());
}

void ConvertTrj::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* /* pbc */, TrajectoryAnalysisModuleData* /*pdata*/)
{
    output_->prepareAndWriteFrame(frnr, fr);
}

void ConvertTrj::finishAnalysis(int /*nframes*/) {}


void ConvertTrj::writeOutput() {}

} // namespace

const char ConvertTrjInfo::name[]             = "convert-trj";
const char ConvertTrjInfo::shortDescription[] = "Converts between different trajectory types";

TrajectoryAnalysisModulePointer ConvertTrjInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ConvertTrj);
}

} // namespace analysismodules

} // namespace gmx
