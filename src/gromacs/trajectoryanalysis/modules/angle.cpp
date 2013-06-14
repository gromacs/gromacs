/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::Angle.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "angle.h"

#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"

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

const char Angle::name[]             = "gangle";
const char Angle::shortDescription[] =
    "Calculate angles";

Angle::Angle()
    : TrajectoryAnalysisModule(name, shortDescription),
      sel1info_(NULL), sel2info_(NULL), binWidth_(1.0), natoms1_(0), natoms2_(0)
{
    averageModule_.reset(new AnalysisDataFrameAverageModule());
    angles_.addModule(averageModule_);
    histogramModule_.reset(new AnalysisDataSimpleHistogramModule());
    angles_.addModule(histogramModule_);

    registerAnalysisDataset(&angles_, "angle");
    registerBasicDataset(averageModule_.get(), "average");
    registerBasicDataset(&histogramModule_->averager(), "histogram");
}


Angle::~Angle()
{
    for (size_t g = 0; g < vt0_.size(); ++g)
    {
        delete [] vt0_[g];
    }
}


void
Angle::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "g_angle computes different types of angles between vectors.",
        "It supports both vectors defined by two positions and normals of",
        "planes defined by three positions.",
        "The z axis or the local normal of a sphere can also be used as",
        "one of the vectors.",
        "There are also convenience options 'angle' and 'dihedral' for",
        "calculating bond angles and dihedrals defined by three/four",
        "positions.[PAR]",
        "The type of the angle is specified with [TT]-g1[tt] and [TT]-g2[tt].",
        "If [TT]-g1[tt] is [TT]angle[tt] or [TT]dihedral[tt], [TT]-g2[tt]",
        "should not be specified.",
        "In this case, [TT]-group1[tt] should specify one or more selections,",
        "and each should contain triplets or quartets of positions that define",
        "the angles to be calculated.[PAR]",
        "If [TT]-g1[tt] is [TT]vector[tt] or [TT]plane[tt], [TT]-group1[tt]",
        "should specify selections that contain either pairs ([TT]vector[tt])",
        "or triplets ([TT]plane[tt]) of positions. For vectors, the positions",
        "set the endpoints of the vector, and for planes, the three positions",
        "are used to calculate the normal of the plane. In both cases,",
        "[TT]-g2[tt] specifies the other vector to use (see below).[PAR]",
        "With [TT]-g2 vector[tt] or [TT]-g2 plane[tt], [TT]-group2[tt] should",
        "specify another set of vectors. [TT]-group1[tt] and [TT]-group2[tt]",
        "should specify the same number of selections, and for each selection",
        "in [TT]-group1[tt], the corresponding selection in [TT]-group2[tt]",
        "should specify the same number of vectors.[PAR]",
        "With [TT]-g2 sphnorm[tt], each selection in [TT]-group2[tt] should",
        "specify a single position that is the center of the sphere.",
        "The second vector is calculated as the vector from the center to the",
        "midpoint of the positions specified by [TT]-group1[tt].[PAR]",
        "With [TT]-g2 z[tt], [TT]-group2[tt] is not necessary, and angles",
        "between the first vectors and the positive Z axis are calculated.[PAR]",
        "With [TT]-g2 t0[tt], [TT]-group2[tt] is not necessary, and angles",
        "are calculated from the vectors as they are in the first frame.[PAR]",
        "There are three options for output:",
        "[TT]-oav[tt] writes an xvgr file with the time and the average angle",
        "for each frame.",
        "[TT]-oall[tt] writes all the individual angles.",
        "[TT]-oh[tt] writes a histogram of the angles. The bin width can be",
        "set with [TT]-binw[tt].",
        "For [TT]-oav[tt] and [TT]-oh[tt], separate average/histogram is",
        "computed for each selection in [TT]-group1[tt]."
        /* TODO: Consider if the dump option is necessary and how to best
         * implement it.
           "[TT]-od[tt] can be used to dump all the individual angles,",
           "each on a separate line. This format is better suited for",
           "further processing, e.g., if angles from multiple runs are needed."
         */
    };
    static const char *const cGroup1TypeEnum[] =
    { "angle", "dihedral", "vector", "plane" };
    static const char *const cGroup2TypeEnum[] =
    { "none", "vector", "plane", "t0", "z", "sphnorm" };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("oav").filetype(eftPlot).outputFile()
                           .store(&fnAverage_).defaultBasename("angaver")
                           .description("Average angles as a function of time"));
    options->addOption(FileNameOption("oall").filetype(eftPlot).outputFile()
                           .store(&fnAll_).defaultBasename("angles")
                           .description("All angles as a function of time"));
    options->addOption(FileNameOption("oh").filetype(eftPlot).outputFile()
                           .store(&fnHistogram_).defaultBasename("anghist")
                           .description("Histogram of the angles"));

    options->addOption(StringOption("g1").enumValue(cGroup1TypeEnum)
                           .defaultEnumIndex(0).store(&g1type_)
                           .description("Type of analysis/first vector group"));
    options->addOption(StringOption("g2").enumValue(cGroup2TypeEnum)
                           .defaultEnumIndex(0).store(&g2type_)
                           .description("Type of second vector group"));
    options->addOption(DoubleOption("binw").store(&binWidth_)
                           .description("Binwidth for -oh in degrees"));

    // TODO: Consider what is the best way to support dynamic selections.
    // Again, most of the code already supports it, but it needs to be
    // considered how should -oall work, and additional checks should be added.
    sel1info_ = options->addOption(SelectionOption("group1")
                                       .required().onlyStatic().storeVector(&sel1_)
                                       .multiValue()
                                       .description("First analysis/vector selection"));
    sel2info_ = options->addOption(SelectionOption("group2")
                                       .onlyStatic().storeVector(&sel2_)
                                       .multiValue()
                                       .description("Second analysis/vector selection"));
}


void
Angle::optionsFinished(Options *options, TrajectoryAnalysisSettings *settings)
{
    bool bSingle = (g1type_[0] == 'a' || g1type_[0] == 'd');

    if (bSingle && g2type_[0] != 'n')
    {
        GMX_THROW(InconsistentInputError("Cannot use a second group (-g2) with "
                                         "-g1 angle or dihedral"));
    }
    if (bSingle && options->isSet("group2"))
    {
        GMX_THROW(InconsistentInputError("Cannot provide a second selection "
                                         "(-group2) with -g1 angle or dihedral"));
    }
    if (!bSingle && g2type_[0] == 'n')
    {
        GMX_THROW(InconsistentInputError("Should specify a second group (-g2) "
                                         "if the first group is not an angle or a dihedral"));
    }

    // Set up the number of positions per angle.
    switch (g1type_[0])
    {
        case 'a': natoms1_ = 3; break;
        case 'd': natoms1_ = 4; break;
        case 'v': natoms1_ = 2; break;
        case 'p': natoms1_ = 3; break;
        default:
            GMX_THROW(InternalError("invalid -g1 value"));
    }
    switch (g2type_[0])
    {
        case 'n': natoms2_ = 0; break;
        case 'v': natoms2_ = 2; break;
        case 'p': natoms2_ = 3; break;
        case 't': natoms2_ = 0; break;
        case 'z': natoms2_ = 0; break;
        case 's': natoms2_ = 1; break;
        default:
            GMX_THROW(InternalError("invalid -g2 value"));
    }
    if (natoms2_ == 0 && options->isSet("group2"))
    {
        GMX_THROW(InconsistentInputError("Cannot provide a second selection (-group2) with -g2 t0 or z"));
    }
}


void
Angle::checkSelections(const SelectionList &sel1,
                       const SelectionList &sel2) const
{
    if (natoms2_ > 0 && sel1.size() != sel2.size())
    {
        GMX_THROW(InconsistentInputError(
                          "-group1 and -group2 should specify the same number of selections"));
    }

    for (size_t g = 0; g < sel1.size(); ++g)
    {
        int na1 = sel1[g].posCount();
        int na2 = (natoms2_ > 0) ? sel2[g].posCount() : 0;
        if (natoms1_ > 1 && na1 % natoms1_ != 0)
        {
            GMX_THROW(InconsistentInputError(formatString(
                                                     "Number of positions in selection %d in the first group not divisible by %d",
                                                     static_cast<int>(g + 1), natoms1_)));
        }
        if (natoms2_ > 1 && na2 % natoms2_ != 0)
        {
            GMX_THROW(InconsistentInputError(formatString(
                                                     "Number of positions in selection %d in the second group not divisible by %d",
                                                     static_cast<int>(g + 1), natoms2_)));
        }
        if (natoms1_ > 0 && natoms2_ > 1 && na1 / natoms1_ != na2 / natoms2_)
        {
            GMX_THROW(InconsistentInputError(
                              "Number of vectors defined by the two groups are not the same"));
        }
        if (g2type_[0] == 's' && sel2[g].posCount() != 1)
        {
            GMX_THROW(InconsistentInputError(
                              "The second group should contain a single position with -g2 sphnorm"));
        }
    }
}


void
Angle::initAnalysis(const TrajectoryAnalysisSettings &settings,
                    const TopologyInformation        &top)
{
    checkSelections(sel1_, sel2_);

    // checkSelections() ensures that both selection lists have the same size.
    angles_.setDataSetCount(sel1_.size());
    for (size_t i = 0; i < sel1_.size(); ++i)
    {
        angles_.setColumnCount(i, sel1_[i].posCount() / natoms1_);
    }
    double histogramMin = (g1type_ == "dihedral" ? -180.0 : 0);
    histogramModule_->init(histogramFromRange(histogramMin, 180.0)
                               .binWidth(binWidth_).includeAll());

    if (g2type_ == "t0")
    {
        vt0_.resize(sel1_.size());
        for (size_t g = 0; g < sel1_.size(); ++g)
        {
            vt0_[g] = new rvec[sel1_[g].posCount() / natoms1_];
        }
    }

    if (!fnAverage_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAverage_);
        plotm->setTitle("Average angle");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Angle (degrees)");
        // TODO: Add legends
        averageModule_->addModule(plotm);
    }

    if (!fnAll_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAll_);
        plotm->setTitle("Angle");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Angle (degrees)");
        // TODO: Add legends? (there can be a massive amount of columns)
        angles_.addModule(plotm);
    }

    if (!fnHistogram_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnHistogram_);
        plotm->setTitle("Angle histogram");
        plotm->setXLabel("Angle (degrees)");
        plotm->setYLabel("Probability");
        // TODO: Add legends
        histogramModule_->averager().addModule(plotm);
    }
}


//! Helper method to process selections into an array of coordinates.
static void
copy_pos(const SelectionList &sel, int natoms, int g, int first, rvec x[])
{
    for (int k = 0; k < natoms; ++k)
    {
        copy_rvec(sel[g].position(first + k).x(), x[k]);
    }
}


//! Helper method to calculate a vector from two or three positions.
static void
calc_vec(int natoms, rvec x[], t_pbc *pbc, rvec xout, rvec cout)
{
    switch (natoms)
    {
        case 2:
            if (pbc)
            {
                pbc_dx(pbc, x[1], x[0], xout);
            }
            else
            {
                rvec_sub(x[1], x[0], xout);
            }
            svmul(0.5, xout, cout);
            rvec_add(x[0], cout, cout);
            break;
        case 3:
        {
            rvec v1, v2;
            if (pbc)
            {
                pbc_dx(pbc, x[1], x[0], v1);
                pbc_dx(pbc, x[2], x[0], v2);
            }
            else
            {
                rvec_sub(x[1], x[0], v1);
                rvec_sub(x[2], x[0], v2);
            }
            cprod(v1, v2, xout);
            rvec_add(x[0], x[1], cout);
            rvec_add(cout, x[2], cout);
            svmul(1.0/3.0, cout, cout);
            break;
        }
        default:
            GMX_RELEASE_ASSERT(false, "Incorrectly initialized number of atoms");
    }
}


void
Angle::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                    TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle       dh   = pdata->dataHandle(angles_);
    const SelectionList     &sel1 = pdata->parallelSelections(sel1_);
    const SelectionList     &sel2 = pdata->parallelSelections(sel2_);

    checkSelections(sel1, sel2);

    dh.startFrame(frnr, fr.time);

    for (size_t g = 0; g < sel1_.size(); ++g)
    {
        rvec  v1, v2;
        rvec  c1, c2;
        switch (g2type_[0])
        {
            case 'z':
                clear_rvec(v2);
                v2[ZZ] = 1.0;
                clear_rvec(c2);
                break;
            case 's':
                copy_rvec(sel2_[g].position(0).x(), c2);
                break;
        }
        dh.selectDataSet(g);
        for (int i = 0, j = 0, n = 0;
             i < sel1[g].posCount();
             i += natoms1_, j += natoms2_, ++n)
        {
            rvec x[4];
            real angle;
            copy_pos(sel1, natoms1_, g, i, x);
            switch (g1type_[0])
            {
                case 'a':
                    if (pbc)
                    {
                        pbc_dx(pbc, x[0], x[1], v1);
                        pbc_dx(pbc, x[2], x[1], v2);
                    }
                    else
                    {
                        rvec_sub(x[0], x[1], v1);
                        rvec_sub(x[2], x[1], v2);
                    }
                    angle = gmx_angle(v1, v2);
                    break;
                case 'd':
                {
                    rvec dx[3];
                    if (pbc)
                    {
                        pbc_dx(pbc, x[0], x[1], dx[0]);
                        pbc_dx(pbc, x[2], x[1], dx[1]);
                        pbc_dx(pbc, x[2], x[3], dx[2]);
                    }
                    else
                    {
                        rvec_sub(x[0], x[1], dx[0]);
                        rvec_sub(x[2], x[1], dx[1]);
                        rvec_sub(x[2], x[3], dx[2]);
                    }
                    cprod(dx[0], dx[1], v1);
                    cprod(dx[1], dx[2], v2);
                    angle = gmx_angle(v1, v2);
                    real ipr = iprod(dx[0], v2);
                    if (ipr < 0)
                    {
                        angle = -angle;
                    }
                    break;
                }
                case 'v':
                case 'p':
                    calc_vec(natoms1_, x, pbc, v1, c1);
                    switch (g2type_[0])
                    {
                        case 'v':
                        case 'p':
                            copy_pos(sel2, natoms2_, 0, j, x);
                            calc_vec(natoms2_, x, pbc, v2, c2);
                            break;
                        case 't':
                            // FIXME: This is not parallelizable.
                            if (frnr == 0)
                            {
                                copy_rvec(v1, vt0_[g][n]);
                            }
                            copy_rvec(vt0_[g][n], v2);
                            break;
                        case 'z':
                            c1[XX] = c1[YY] = 0.0;
                            break;
                        case 's':
                            if (pbc)
                            {
                                pbc_dx(pbc, c1, c2, v2);
                            }
                            else
                            {
                                rvec_sub(c1, c2, v2);
                            }
                            break;
                        default:
                            GMX_THROW(InternalError("invalid -g2 value"));
                    }
                    angle = gmx_angle(v1, v2);
                    break;
                default:
                    GMX_THROW(InternalError("invalid -g1 value"));
            }
            /* TODO: Should we also calculate distances like g_sgangle?
             * Could be better to leave that for a separate tool.
               real dist = 0.0;
               if (bDumpDist_)
               {
                if (pbc)
                {
                    rvec dx;
                    pbc_dx(pbc, c2, c1, dx);
                    dist = norm(dx);
                }
                else
                {
                    dist = sqrt(distance2(c1, c2));
                }
               }
             */
            dh.setPoint(n, angle * RAD2DEG);
        }
    }
    dh.finishFrame();
}


void
Angle::finishAnalysis(int /*nframes*/)
{
    AbstractAverageHistogram &averageHistogram = histogramModule_->averager();
    averageHistogram.normalizeProbability();
    averageHistogram.done();
}


void
Angle::writeOutput()
{
}

} // namespace modules

} // namespace gmxana
