/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::ConvertTrj.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "convert-trj.h"

#include <algorithm>

#include "gromacs/coordinateio/builder.h"
#include "gromacs/coordinateio/flags.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

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

        void initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings) override;
        void optionsFinished(TrajectoryAnalysisSettings *settings) override;
        void initAnalysis(const TrajectoryAnalysisSettings  &settings,
                          const TopologyInformation         &top) override;
        void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata) override;

        void finishAnalysis(int nframes) override;
        void writeOutput() override;

    private:
        OutputManagerPointer  output_;
        Selection             sel_;
        std::string           name_;
        OutputRequirements    requirements_;
};

ConvertTrj::ConvertTrj()
{
}


void
ConvertTrj::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] converts trajectory files between different formats.",
        "The module supports writing all GROMACS supported file formats from",
        "the supported input formats.",
        "",
        "Included is also a selection of possible options to change meta-information",
        "in the underlying trajectory data to obtain slimmer output files, or to",
        "direct requirements for which outputs are necessary.",
        "",
        "The module also can generate subsets of trajectories based on user supplied",
        "selections.",
    };

    options->addOption(SelectionOption("select")
                           .store(&sel_)
                           .onlyAtoms()
                           .description("Selection of atoms to write to the file"));

    options->addOption(FileNameOption("o")
                           .filetype(eftTrajectory)
                           .outputFile()
                           .store(&name_).defaultBasename("trajout")
                           .required()
                           .description("Output trajectory after conversion"));

    initOutputRequirementFileOptions(options, &requirements_);

    settings->setHelpText(desc);
    // set correct flags to indicate we need a proper topology for the analysis
    // settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
ConvertTrj::optionsFinished(TrajectoryAnalysisSettings * settings)
{
    checkRequirementsOptions(&requirements_);
    int frameFlags = TRX_NEED_X;

    frameFlags |= TRX_READ_V;
    frameFlags |= TRX_READ_F;

    settings->setFrameFlags(frameFlags);

}


void
ConvertTrj::initAnalysis(const TrajectoryAnalysisSettings    & /*settings*/,
                         const TopologyInformation          &top)
{
    output_ = createOutputManager(top.hasTopology() ? top.mtop() : nullptr,
                                  sel_,
                                  name_,
                                  top.copyAtoms(),
                                  requirements_);
}

void
ConvertTrj::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /* pbc */,
                         TrajectoryAnalysisModuleData * /*pdata*/)
{
    // modify frame to write out correct number of coords
    // and actually write out
    output_->prepareFrame(frnr, fr);
}

void
ConvertTrj::finishAnalysis(int /*nframes*/)
{
}



void
ConvertTrj::writeOutput()
{
}

}       // namespace

const char ConvertTrjInfo::name[]             = "convert-trj";
const char ConvertTrjInfo::shortDescription[] =
    "Converts between different trajectory types";

TrajectoryAnalysisModulePointer ConvertTrjInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ConvertTrj);
}

} // namespace analysismodules

} // namespace gmx
