/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Angle.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "angle.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "pbc.h"
#include "vec.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/selection/selectionoptioninfo.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

const char Angle::name[] = "angle";
const char Angle::shortDescription[] =
    "Calculate angles";

Angle::Angle()
    : options_(name, shortDescription),
      sel1info_(NULL), sel2info_(NULL),
      bSplit1_(false), bSplit2_(false), bMulti_(false), bAll_(false),
      bDumpDist_(false), natoms1_(0), natoms2_(0), vt0_(NULL)
{
    registerAnalysisDataset(&data_, "angle");
}


Angle::~Angle()
{
    delete[] vt0_;
}


Options &
Angle::initOptions(TrajectoryAnalysisSettings *settings)
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
        "In this case, one selection is required, and it should contain",
        "triplets or quartets of positions that define the angles to be",
        "calculated.",
        "If [TT]-g1[tt] is not [TT]angle[tt] or [TT]dihedral[tt], [TT]-g2[tt]",
        "should not be [TT]none[tt], and the two options define two vectors",
        "for the calculation. For vectors ([TT]vector[tt]), a selection with",
        "pairs of positions is required, and for planes ([TT]plane[tt]),",
        "triplets of positions are required.",
        "If both vectors are specified by positions, the number of vectors",
        "should be the same in both selections.",
        "[TT]-g2 sphnorm[tt] requires a reference selection that defines",
        "the center of the sphere.",
        "[TT]-g2 z[tt] does not require any selection.[PAR]",
        "With [TT]-split1[tt], the positions for [TT]-g1[tt] are specified",
        "using N separate selections with M positions each, instead of the",
        "default M*N positions in one selection.",
        "[TT]-split2[tt] does the same for [TT]-g2[tt].[PAR]",
        "There are two options for output:",
        "[TT]-o[tt] writes an xvgr file with the time and the average angle",
        "for each frame.",
        "With [TT]-all[tt], also the individual angles are written (only",
        "supported for static selections).",
        "[TT]-od[tt] can be used to dump all the individual angles,",
        "each on a separate line. This format is better suited for",
        "further processing, e.g., if angles from multiple runs are needed."
    };
    static const char *const cGroup1TypeEnum[] =
        { "angle", "dihedral", "vector", "plane", NULL };
    static const char *const cGroup2TypeEnum[] =
        { "none", "vector", "plane", "t0", "z", "sphnorm", NULL };

    options_.setDescription(concatenateStrings(desc));

    options_.addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnAngle_).defaultValueIfSet("angle"));
    options_.addOption(FileNameOption("od").filetype(eftPlot).outputFile()
                           .store(&fnDump_).defaultValueIfSet("angdump"));

    options_.addOption(StringOption("g1").enumValue(cGroup1TypeEnum)
        .defaultEnumIndex(0).store(&g1type_)
        .description("Type of analysis/first vector group"));
    options_.addOption(StringOption("g2").enumValue(cGroup2TypeEnum)
        .defaultEnumIndex(0).store(&g2type_)
        .description("Type of second vector group"));
    options_.addOption(BooleanOption("split1").store(&bSplit1_)
        .description("Each position of first group in separate selection"));
    options_.addOption(BooleanOption("split2").store(&bSplit2_)
        .description("Each position of second group in separate selection"));
    options_.addOption(BooleanOption("multi").store(&bMulti_)
        .description("Analyze multiple sets of angles/dihedrals"));
    options_.addOption(BooleanOption("all").store(&bAll_)
        .description("Print individual angles together with the average"));
    options_.addOption(BooleanOption("dumpd").store(&bDumpDist_)
        .description("Write also distances with -od"));

    options_.addOption(SelectionOption("group1").multiValue().required()
        .dynamicOnlyWhole().storeVector(&sel1_).getAdjuster(&sel1info_)
        .description("First analysis/vector selection"));
    options_.addOption(SelectionOption("group2").multiValue()
        .dynamicOnlyWhole().storeVector(&sel2_).getAdjuster(&sel2info_)
        .description("Second analysis/vector selection"));

    return options_;
}


void
Angle::initOptionsDone(TrajectoryAnalysisSettings *settings)
{
    // Validity checks.
    bool bSingle = (g1type_[0] == 'a' || g1type_[0] == 'd');

    if (bSingle && g2type_[0] != 'n')
    {
        GMX_THROW(InconsistentInputError("Cannot use a second group (-g2) with "
                                         "-g1 angle or dihedral"));
    }
    if (bSingle && options_.isSet("group2"))
    {
        GMX_THROW(InconsistentInputError("Cannot provide a second selection "
                                         "(-group2) with -g1 angle or dihedral"));
    }
    if (!bSingle && g2type_[0] == 'n')
    {
        GMX_THROW(InconsistentInputError("Should specify a second group (-g2) "
                                         "if the first group is not an angle or a dihedral"));
    }
    if (bSingle && bDumpDist_)
    {
        GMX_THROW(InconsistentInputError("Cannot calculate distances with -g1 angle or dihedral"));
        // bDumpDist_ = false;
    }
    if (bMulti_ && !bSingle)
    {
        GMX_THROW(InconsistentInputError("-mult can only be combined with -g1 angle or dihedral"));
    }
    if (bMulti_ && bSplit1_)
    {
        GMX_THROW(InconsistentInputError("-mult can not be combined with -split1"));
    }
    if (bMulti_ && bAll_)
    {
        GMX_THROW(InconsistentInputError("-mult and -all are mutually exclusive options"));
    }

    if (bAll_)
    {
        sel1info_->setOnlyStatic(true);
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
    if (natoms2_ == 0 && options_.isSet("group2"))
    {
        GMX_THROW(InconsistentInputError("Cannot provide a second selection (-group2) with -g2 t0 or z"));
    }

    if (!bMulti_)
    {
        sel1info_->setValueCount(bSplit1_ ? natoms1_ : 1);
    }
    if (natoms2_ > 0)
    {
        sel2info_->setValueCount(bSplit2_ ? natoms2_ : 1);
    }
}


void
Angle::checkSelections(const SelectionList &sel1,
                       const SelectionList &sel2) const
{
    if (bMulti_)
    {
        for (size_t g = 0; g < sel1.size(); ++g)
        {
            if (sel1[g].posCount() % natoms1_ != 0)
            {
                GMX_THROW(InconsistentInputError(formatString(
                    "Number of positions in selection %d not divisible by %d",
                    static_cast<int>(g + 1), natoms1_)));
            }
        }
        return;
    }

    int na1 = sel1[0].posCount();
    int na2 = (natoms2_ > 0) ? sel2[0].posCount() : 0;

    if (!bSplit1_ && natoms1_ > 1 && na1 % natoms1_ != 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
            "Number of positions in the first group not divisible by %d",
            natoms1_)));
    }
    if (!bSplit2_ && natoms2_ > 1 && na2 % natoms2_ != 0)
    {
        GMX_THROW(InconsistentInputError(formatString(
            "Number of positions in the second group not divisible by %d",
            natoms2_)));
    }

    if (bSplit1_)
    {
        for (int g = 1; g < natoms1_; ++g)
        {
            if (sel1[g].posCount() != na1)
            {
                GMX_THROW(InconsistentInputError(
                          "All selections in the first group should contain "
                          "the same number of positions"));
            }
        }
    }
    else
    {
        na1 /= natoms1_;
    }
    if (natoms2_ > 1)
    {
        if (bSplit2_)
        {
            for (int g = 1; g < natoms2_; ++g)
            {
                if (sel2[g].posCount() != na2)
                {
                    GMX_THROW(InconsistentInputError(
                              "All selections in the second group should contain "
                              "the same number of positions"));
                }
            }
        }
        else
        {
            na2 /= natoms2_;
        }
    }
    if (natoms1_ > 0 && natoms2_ > 1 && na1 != na2)
    {
        GMX_THROW(InconsistentInputError(
                  "Number of vectors defined by the two groups are not the same"));
    }
    if (g2type_[0] == 's' && sel2[0].posCount() != 1)
    {
        GMX_THROW(InconsistentInputError(
                  "The second group should contain a single position with -g2 sphnorm"));
    }
}


void
Angle::initAnalysis(const TrajectoryAnalysisSettings &settings,
                    const TopologyInformation &top)
{
    checkSelections(sel1_, sel2_);

    if (bMulti_)
    {
        data_.setColumnCount(sel1_.size());
    }
    else if (bAll_)
    {
        int na = sel1_[0].posCount();
        if (!bSplit1_)
        {
            na /= natoms1_;
        }
        data_.setColumnCount(na + 1);
    }
    else
    {
        data_.setColumnCount(1);
    }

    if (g2type_ == "t0")
    {
        int na = sel1_[0].posCount();
        if (!bSplit1_)
        {
            na /= natoms1_;
        }
        vt0_ = new rvec[na];
    }

    AnalysisDataPlotModulePointer plotm(
        new AnalysisDataPlotModule(settings.plotSettings()));
    plotm->setFileName(fnAngle_);
    plotm->setTitle("Angle");
    plotm->setXAxisIsTime();
    plotm->setYLabel("Angle (degrees)");
    data_.addModule(plotm);
}


//! Helper method to process selections into an array of coordinates.
static void
copy_pos(const SelectionList &sel, bool bSplit, int natoms,
         int firstg, int first, rvec x[])
{
    if (bSplit)
    {
        for (int k = 0; k < natoms; ++k)
        {
            copy_rvec(sel[firstg + k].position(first).x(), x[k]);
        }
    }
    else
    {
        for (int k = 0; k < natoms; ++k)
        {
            copy_rvec(sel[firstg].position(first + k).x(), x[k]);
        }
    }
}


//! Helper method to calculate a vector from two or three positions..
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
        case 3: {
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
    AnalysisDataHandle       dh = pdata->dataHandle(data_);
    const SelectionList     &sel1 = pdata->parallelSelections(sel1_);
    const SelectionList     &sel2 = pdata->parallelSelections(sel2_);

    checkSelections(sel1, sel2);

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
            copy_rvec(sel2_[0].position(0).x(), c2);
            break;
    }

    dh.startFrame(frnr, fr.time);

    int incr1 = bSplit1_ ? 1 : natoms1_;
    int incr2 = bSplit2_ ? 1 : natoms2_;
    int ngrps = bMulti_ ? sel1_.size() : 1;

    for (int g = 0; g < ngrps; ++g)
    {
        real ave = 0.0;
        int n = 0;
        int i, j;
        for (i = j = 0; i < sel1[g].posCount(); i += incr1)
        {
            rvec x[4];
            real angle;
            copy_pos(sel1, bSplit1_, natoms1_, g, i, x);
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
                case 'd': {
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
                            copy_pos(sel2, bSplit2_, natoms2_, 0, j, x);
                            calc_vec(natoms2_, x, pbc, v2, c2);
                            j += incr2;
                            break;
                        case 't':
                            // FIXME: This is not parallelizable.
                            if (frnr == 0)
                            {
                                copy_rvec(v1, vt0_[n]);
                            }
                            copy_rvec(vt0_[n], v2);
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
            angle *= RAD2DEG;
            /* TODO: add support for -od and -dumpd 
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
            if (bAll_)
            {
                dh.setPoint(n + 1, angle);
            }
            ave += angle;
            ++n;
        }
        if (n > 0)
        {
            ave /= n;
        }
        dh.setPoint(g, ave);
    }
    dh.finishFrame();
}


void
Angle::finishAnalysis(int /*nframes*/)
{
}


void
Angle::writeOutput()
{
}

} // namespace modules

} // namespace gmxana
