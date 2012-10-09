/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
 * Implements gmx::analysismodules::ProteinFilamentGeometry
 *
 * \author Alexey Shvetsov <alexxy@omrb.pnpi.spb.ru>
 * \ingroup module_trajectoryanalysis
 */
#include "pfg.h"

#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/princ.h"
#include "gromacs/legacyheaders/do_fit.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

const char ProteinFilamentGeometry::name[]             = "pfg";
const char ProteinFilamentGeometry::shortDescription[] =
    "Calculates protein filament geometry";

ProteinFilamentGeometry::ProteinFilamentGeometry()
    : TrajectoryAnalysisModule(name, shortDescription)
{
    pitch_.setColumnCount(1);
    registerAnalysisDataset(&pitch_, "pitch");
    rotation_.setColumnCount(1);
    registerAnalysisDataset(&rotation_, "rotation");
    nmon_.setColumnCount(1);
    registerAnalysisDataset(&nmon_, "nmon");
}


ProteinFilamentGeometry::~ProteinFilamentGeometry()
{
}

void
ProteinFilamentGeometry::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "This tool compute protein filament properties as a quantity of time such as:",
        "[PAR]",
        "helical pitch",
        "[PAR]",
        "number of monomers per turn",
        "[PAR]",
        "helical rotation",
        "[PAR]",
        "As input you need at lest one filament period and any two neighbour monomers"
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("pitch").filetype(eftPlot).outputFile()
                           .store(&fnPitch_).defaultBasename("pitch")
                           .description("Helical pitch as a function of time"));
    options->addOption(FileNameOption("rot").filetype(eftPlot).outputFile()
                           .store(&fnRotation_).defaultBasename("rotation")
                           .description("Helical rotation per monomer as a function of time"));
    options->addOption(FileNameOption("nmon").filetype(eftPlot).outputFile()
                           .store(&fnNMonomers_).defaultBasename("nmonomers")
                           .description("Number of monomers per turn"));

    options->addOption(SelectionOption("sel").required().valueCount(2)
                           .store(sel_).description("Select two monomers"));

}

void
ProteinFilamentGeometry::checkSelections(const Selection &sel1,
                                         const Selection &sel2) const
{
    if (sel1.atomCount() != sel2.atomCount() )
    {
        GMX_THROW(InconsistentInputError(
                          "Both monomers should have same number of analyzed atoms"));
    }
}


void
ProteinFilamentGeometry::initAnalysis(const TrajectoryAnalysisSettings &settings,
                                      const TopologyInformation        &top)
{
    checkSelections(sel_[0], sel_[1]);

    if (!fnPitch_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnPitch_);
        plotm->setTitle("Helical Pitch");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Helical Pitch (nm)");
        pitch_.addModule(plotm);
    }

    if (!fnRotation_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnRotation_);
        plotm->setTitle("Helical rotation per monomer");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Rotation (degrees)");
        rotation_.addModule(plotm);
    }

    if (!fnNMonomers_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnNMonomers_);
        plotm->setTitle("Helical number of monomers per turn");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Number of monomers");
        nmon_.addModule(plotm);
    }
}

void ProteinFilamentGeometry::calcCenterOfMass(const Selection &sel, rvec xcm, gmx_bool bMass, gmx_bool bQ)
{
    rvec x;
    real cm = 0.0;

    clear_rvec(xcm);

    for (int i = 0; i < sel.posCount(); ++i)
    {
        clear_rvec(x);
        if (bMass)
        {
            cm += sel.position(i).mass();
            svmul(sel.position(i).mass(), sel.position(i).x(), x);
        }
        else if (bQ)
        {
            cm += sel.position(i).charge();
            svmul(sel.position(i).charge(), sel.position(i).x(), x);
        }
        else
        {
            cm += 1.0;
            copy_rvec(sel.position(i).x(), x);
        }

        rvec_add(xcm, x, xcm);
    }
    svmul(1.0/cm, xcm, xcm);
}

void ProteinFilamentGeometry::resetCenterOfMass(const Selection &sel, rvec xreset[],  gmx_bool bMass, gmx_bool bQ)
{
    rvec xcm;
    calcCenterOfMass(sel, xcm, bMass, bQ);

    for (int i = 0; i < sel.posCount(); ++i)
    {
        rvec_sub(sel.position(i).x(), xcm, xreset[i]);
    }
}

void
ProteinFilamentGeometry::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                      TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle       dhpitch = pdata->dataHandle(pitch_);
    AnalysisDataHandle       dhnmon  = pdata->dataHandle(nmon_);
    AnalysisDataHandle       dhrot   = pdata->dataHandle(rotation_);
    const Selection         &sel1    = pdata->parallelSelection(sel_[0]);
    const Selection         &sel2    = pdata->parallelSelection(sel_[1]);

    rvec                     xcm1;
    rvec                     xcm2;
    rvec                     dx;
    rvec                     e;

    real                     h, npturn, pitch, theta;

    matrix                   R;

    const int                NDIM = 3;

    checkSelections(sel1, sel2);

    dhpitch.startFrame(frnr, fr.time);
    dhnmon.startFrame(frnr, fr.time);
    dhrot.startFrame(frnr, fr.time);

    /*
     * calculate Center of mass for monomers
     */
    calcCenterOfMass(sel1, xcm1, TRUE, FALSE);
    calcCenterOfMass(sel2, xcm2, TRUE, FALSE);

    rvec_sub(xcm2, xcm1, dx);
    /*
     * reset CoM to coordinate origin
     */
    rvec *xmon1 = new rvec[sel1.atomCount()];
    rvec *xmon2 = new rvec[sel2.atomCount()];
    resetCenterOfMass(sel1, xmon1, TRUE, FALSE);
    resetCenterOfMass(sel2, xmon2, TRUE, FALSE);
    /*
     * Make and w_rls array for calc_fit_R
     */
    real  *w_rls = new real[sel1.atomCount()];
    for (int i = 0; i < sel1.posCount(); ++i)
    {
        w_rls[i] = sel1.position(i).mass();
    }

    calc_fit_R(NDIM, sel1.atomCount(), w_rls, xmon1, xmon2, R);
    delete[] xmon1;
    delete[] xmon2;
    delete[] w_rls;
    /*
     * now get theta rotation angle
     */
    theta = acos(0.5*(R[XX][XX] + R[YY][YY] + R[ZZ][ZZ] - 1.0));
    /* now we can get rotation axis if angle wasnt nPi */
    e[XX] = (R[ZZ][YY] - R[YY][ZZ])/(2.0*sin(theta));
    e[YY] = (R[XX][ZZ] - R[ZZ][XX])/(2.0*sin(theta));
    e[ZZ] = (R[YY][XX] - R[XX][YY])/(2.0*sin(theta));
    /*
     * Now we want to compute h per monomer
     */
    h      = iprod(dx, e);
    npturn = M_2PI/theta;
    pitch  = fabs(h)*npturn;

    /*
     * Output
     */
    dhpitch.setPoint(0, pitch);
    dhnmon.setPoint(0, npturn);
    dhrot.setPoint(0, theta*RAD2DEG);

    dhpitch.finishFrame();
    dhnmon.finishFrame();
    dhrot.finishFrame();
}


void
ProteinFilamentGeometry::finishAnalysis(int /*nframes*/)
{
}


void
ProteinFilamentGeometry::writeOutput()
{
}

} // namespace modules

} // namespace gmxana
