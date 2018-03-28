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
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

#include "filehandler.h"

namespace gmx
{

struct t_writeFileBool;
struct t_writeFileDoubles;

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
        Filehandler                          file_;

};

Convert::Convert()
{
}


void
Convert::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] converts trajectory files between different formats."
    };

    file_.initFileOptions(options);
    settings->setHelpText(desc);
    // set correct flags to indicate we need a proper topology for the analysis
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void
Convert::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    file_.checkOptions();
}


void
Convert::initAnalysis(const TrajectoryAnalysisSettings   & /*settings*/,
                      const TopologyInformation         &top)
{
    file_.setMtop(top.mtop());
    file_.initOutput();
}

void
Convert::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /* pbc */,
                      TrajectoryAnalysisModuleData * /*pdata*/)
{
    file_.setBox(fr.box);
    file_.modifyFrame(&fr);
    file_.writeFrame();
}

void
Convert::finishAnalysis(int /*nframes*/)
{
    file_.closeFile();
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
