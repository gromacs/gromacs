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
 * Implements gmx::analysismodules::ChangePBC.
 *
 * \author
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "changepbc.h"

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
#include "gromacs/topology/mtop_lookup.h"
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
#include "pbchandler.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * ChangePBC
 */

class ChangePBC : public TrajectoryAnalysisModule
{
    public:
        ChangePBC();

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
        Pbchandler                           pbc_;
};

ChangePBC::ChangePBC()
{
}


void
ChangePBC::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] allows the conversion between different representations of periodic boundary conditions.",
        "It supports all of the coordinate file formats that are native to GROMACS. File format types are",
        "detected according to the file extension [REF].xtc[ref], [REF].trr[ref], [REF].pdb[ref], [REF].gro[ref]",
        "[REF].g96[ref] and [REF].tng[ref].",
        "The new periodic boundary conditions can be defined according to the following options:",
        "",
        " * [TT]mol[tt] puts the center of mass of molecules in the box,",
        "   and requires a run input file to be supplied with [TT]-s[tt].",
        " * [TT]res[tt] puts the center of mass of residues in the box.",
        " * [TT]atom[tt] puts all the atoms in the box.",
        " * [TT]nojump[tt] checks if atoms jump across the box and then puts",
        "   them back. This has the effect that all molecules",
        "   will remain whole (provided they were whole in the initial",
        "   conformation). [BB]Note[bb] that this ensures a continuous trajectory but",
        "   molecules may diffuse out of the box. The starting configuration",
        "   for this procedure is taken from the structure file, if one is",
        "   supplied, otherwise it is the first frame.",
        " * [TT]cluster[tt] clusters all the atoms in the selected index",
        "   such that they are all closest to the center of mass of the cluster,",
        "   which is iteratively updated. [BB]Note[bb] that this will only give meaningful",
        "   results if you in fact have a cluster. Luckily that can be checked",
        "   afterwards using a trajectory viewer. Note also that if your molecules",
        "   are broken this will not work either.",
        " * [TT]whole[tt] only makes broken molecules whole.",
        "",
        "Option [TT]-ur[tt] sets the unit cell representation for options",
        "[TT]mol[tt], [TT]res[tt] and [TT]atom[tt] of [TT]-pbc[tt].",
        "All three options give different results for triclinic boxes and",
        "identical results for rectangular boxes.",
        "[TT]rect[tt] is the ordinary brick shape.",
        "[TT]tric[tt] is the triclinic unit cell.",
        "[TT]compact[tt] puts all atoms at the closest distance from the center",
        "of the box. This can be useful for visualizing e.g. truncated octahedra",
        "or rhombic dodecahedra. The center for options [TT]tric[tt] and [TT]compact[tt]",
        "is [TT]tric[tt] (see below), unless the option [TT]-boxcenter[tt]",
        "is set differently.[PAR]",

        "Option [TT]-center[tt] centers the system in the box. The user can",
        "select the group which is used to determine the geometrical center.",
        "Option [TT]-boxcenter[tt] sets the location of the center of the box",
        "for options [TT]-pbc[tt] and [TT]-center[tt]. The center options are:",
        "[TT]tric[tt]: half of the sum of the box vectors,",
        "[TT]rect[tt]: half of the box diagonal,",
        "[TT]zero[tt]: zero.",
        "Use option [TT]-pbc mol[tt] in addition to [TT]-center[tt] when you",
        "want all molecules in the box after the centering.[PAR]",

        "It is not always possible to use combinations of [TT]-pbc[tt],",
        "[TT]-fit[tt], [TT]-ur[tt] and [TT]-center[tt] to do exactly what",
        "you want in one call to [THISMODULE]. Consider using multiple",
        "calls, and check out the GROMACS website for suggestions.[PAR]",
    };

    file_.initFileOptions(options);
    pbc_.initPBCOptions(options);
    settings->setHelpText(desc);
    // set correct flags to indicate we need a proper topology for the analysis
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop |
                      TrajectoryAnalysisSettings::efUseTopX);

}

void
ChangePBC::optionsFinished(TrajectoryAnalysisSettings *settings)
{
    file_.checkOptions();
    pbc_.checkOptions(file_.getSel(), settings);

}


void
ChangePBC::initAnalysis(const TrajectoryAnalysisSettings   & /*settings*/,
                        const TopologyInformation         &top)
{
    file_.setMtop(top.mtop());
    file_.initOutput();
    pbc_.setReferenceCoordinates(top);
}

void
ChangePBC::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /* pbc */,
                        TrajectoryAnalysisModuleData * /*pdata*/)
{
    // first, set correct data in file handler for later output
    file_.modifyFrame(&fr);
    if (file_.getUseNewBox())
    {
        matrix             newMatrix;
        clear_mat(newMatrix);
        std::vector<float> newBox = file_.getNewBoxDimensions();
        for (size_t i = 0; i < DIM; i++)
        {
            newMatrix[i][i] = newBox[i];
        }
        file_.setBox(newMatrix);
        pbc_.setBox(newMatrix);
    }
    else
    {
        file_.setBox(fr.box);
        pbc_.setBox(fr.box);
    }
    // hand of the real work to the pbchandler class
    pbc_.doPBC(&fr, file_.getFrameForModification(), file_.getMtop());
    file_.writeFrame();
}

void
ChangePBC::finishAnalysis(int /*nframes*/)
{
    file_.closeFile();
}


void
ChangePBC::writeOutput()
{
}

}       // namespace

const char ChangePBCInfo::name[]             = "changePBC";
const char ChangePBCInfo::shortDescription[] =
    "Allows changing of the PBC representation between different types for coordinate frames";

TrajectoryAnalysisModulePointer ChangePBCInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ChangePBC);
}

} // namespace analysismodules

} // namespace gmx
