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
 * Implements gmx::analysismodules::Convert.
 *
 * \author 
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "convert.h"

#include <algorithm>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/analysisdata/modules/write.h"
#include "gromacs/analysisdata/modules/framehandler.h"
#include "gromacs/analysisdata/modules/filehandler.h"
#include "gromacs/analysisdata/modules/settings.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * Convert
 */

class Convert : public TrajectoryAnalysisModule
{
    public:
        Convert();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        Selection                       sel_;
        std::string                         fnOutput_;
	AnalysisData				write_;
        TrajectoryDataWriteModulePointer       writeModule_;

        TrajectoryDataWriteSettings              writeSettings_;

};

Convert::Convert() : sel_(nullptr)
{
    registerAnalysisDataset(&write_, "Trajectory");
}


void
Convert::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] converts trajectory files between different formats."
    };

    settings->setHelpText(desc);

    options->addOption(SelectionOption("Output").store(&sel_).dynamicMask()
                            .required()
                           .description("Selection to write out"));

    options->addOption(FileNameOption("o").filetype(eftTrajectory).outputFile()
                            .store(&fnOutput_).defaultBasename("trajout")
                            .required()
                            .description("Output trajectory after conversion"));
    // set correct flags to indicate we need a proper topology for the analysis
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
Convert::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    // set up the correct file names and conversion thingies for new output
}


void
Convert::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation         &top)
{
    // add settings from write module
    writeSettings_ = settings.writeSettings();
    TrajectoryDataWriteModulePointer writeModule(
            new TrajectoryDataWriteModule(writeSettings_));
    writeModule_ = writeModule;
    // initialize the write with the minimal info needed in one function instead of five
    writeModule_->setExternal(&sel_, fnOutput_, top.mtop());

    write_.addModule(writeModule_);

    // one data set only, the frames to be written out
    write_.setDataSetCount(1);
    // one column for the data set, the actual output frames after manipulation
    write_.setColumnCount(0,1);
}

void
Convert::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /* pbc */,
        TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh = pdata->dataHandle(write_);
    dh.startCoordFrame(frnr, fr);
    // work on frame time data set
    dh.selectDataSet(0);
    // set data point with frame time for each frame
    dh.setCoordPoint(0, fr, true);
    // writeModule_->writeFrame(&fr);
    dh.finishFrame();
    // can have the write method in the writetrj class now?
    
    
}

void
Convert::finishAnalysis(int /*nframes*/)
{
}



void
Convert::writeOutput()
{
}

}       // namespace

const char ConvertInfo::name[]             = "convert";
const char ConvertInfo::shortDescription[] =
    "Converts between different trajectory types";

TrajectoryAnalysisModulePointer ConvertInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Convert);
}

} // namespace analysismodules

} // namespace gmx
