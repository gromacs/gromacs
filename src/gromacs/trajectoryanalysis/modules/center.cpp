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
 * Implements gmx::analysismodules::Center.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "center.h"

#include <algorithm>

#include "gromacs/coordinatedata.h"
#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/coordinateio/flags.h"
#include "gromacs/coordinateio/outputmanager.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * Center
 */

class Center : public TrajectoryAnalysisModule
{
    public:
        Center();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings  &settings,
                                  const TopologyInformation         &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        OutputManagerPointer                  output_;
        Selection                             sel_;
        Selection  selCenter_;
        Selection *selCenterPointer_ = nullptr;
        bool       haveSelCenter_    = false;
        std::string                                        name_;
        CoordinateFileWriteFlags                           flags_;
        CenteringType    centerFlag_ = CenteringType::ecNR;
        SetCenterPointer setCenter_;
};

Center::Center()
{
}


void
Center::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] converts trajectory files between different formats."
    };

    options->addOption(SelectionOption("select")
                           .store(&sel_)
                           .dynamicMask()
                           .required()
                           .description("Selection of atoms to write to the file"));

    options->addOption(FileNameOption("o").filetype(eftTrajectory).outputFile()
                           .store(&name_)
                           .defaultBasename("trajout")
                           .required()
                           .description("Output trajectory after conversion"));
    options->addOption(SelectionOption("center")
                           .store(&selCenter_)
                           .dynamicMask()
                           .storeIsSet(&haveSelCenter_)
                           .description("Selection of atoms used for centering system, default is to use all atoms"));
    options->addOption(EnumOption<CenteringType>("type")
                           .enumValue(cCenterTypeEnum)
                           .store(&centerFlag_)
                           .description("How the system should be centered"));
    flags_.initFileOptions(options);

    settings->setHelpText(desc);
    // set correct flags to indicate we need a proper topology for the analysis
    // settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
Center::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    flags_.checkOptions();
    if (!haveSelCenter_)
    {
        selCenterPointer_ = &sel_;
    }
    else
    {
        selCenterPointer_ = &selCenter_;
    }
}


void
Center::initAnalysis(const TrajectoryAnalysisSettings    & /*settings*/,
                     const TopologyInformation          &top)
{
    output_    = compat::make_unique<OutputManager>(name_, &sel_, top.mtop());
    setCenter_ = compat::make_unique<SetCenter>(selCenterPointer_, centerFlag_);
    flags_.registerModules(output_);
}

void
Center::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /* pbc */,
                     TrajectoryAnalysisModuleData * /*pdata*/)
{
    // modify frame to write out correct number of coords
    // and actually write out
    setCenter_->setFrame(fr);
    setCenter_->convertFrame(fr);
    output_->prepareFrame(frnr, setCenter_->getFrame());
}

void
Center::finishAnalysis(int /*nframes*/)
{
}



void
Center::writeOutput()
{
}

}       // namespace

const char CenterInfo::name[]             = "center";
const char CenterInfo::shortDescription[] =
    "Centers atom selection in box during trajectory";

TrajectoryAnalysisModulePointer CenterInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Center);
}

} // namespace analysismodules

} // namespace gmx
