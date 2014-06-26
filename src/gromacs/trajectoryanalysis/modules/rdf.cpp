/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::Rdf.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com> (C++ conversion)
 * \ingroup module_trajectoryanalysis
 */
#include "rdf.h"

#include <cmath>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! \addtogroup module_trajectoryanalysis
//! \{

/********************************************************************
 * Actual analysis module
 */

/*! \brief
 * Implements `gmx rdf` trajectory analysis module.
 */
class Rdf : public TrajectoryAnalysisModule
{
    public:
        Rdf();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void optionsFinished(Options                    *options,
                                     TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void initAfterFirstFrame(const TrajectoryAnalysisSettings &settings,
                                         const t_trxframe                 &fr);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        std::string                               fnRdf_;
        std::string                               fnCumulative_;
        std::string                               surface_;
        AnalysisDataPlotSettings                  plotSettings_;

        Selection                                 refSel_;
        SelectionList                             sel_;

        AnalysisData                              pairDist_;
        AnalysisDataSimpleHistogramModulePointer  rdf_;
        AnalysisNeighborhood                      nb_;

        double                                    binwidth_;
        double                                    cutoff_;
        bool                                      bNormalize_;
        bool                                      bXY_;

        real                                      cut2_;
        real                                      rmax2_;

        // TODO: Use frame-local data for accumulating these.
        double                                    inverseVolumeSum_;
        std::vector<gmx_int64_t>                  pairCountSum_;
        gmx_int64_t                               refCountSum_;

        std::vector<real>                         surfaceDist2_;

        // Copy and assign disallowed by base.
};

Rdf::Rdf()
    : TrajectoryAnalysisModule(RdfInfo::name, RdfInfo::shortDescription),
      rdf_(new AnalysisDataSimpleHistogramModule()),
      binwidth_(0.002), cutoff_(0.0),
      bNormalize_(true), bXY_(false),
      cut2_(0.0), rmax2_(0.0),
      inverseVolumeSum_(0.0), refCountSum_(0)
{
    pairDist_.setMultipoint(true);
    registerAnalysisDataset(&pairDist_, "pairdist");
    pairDist_.addModule(rdf_);
}

void
Rdf::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates radial distribution functions in different ways.",
        "The normal method is around a (set of) particle(s), the other methods",
        "are around the center of mass of a set of particles ([TT]-com[tt])",
        "or to the closest particle in a set ([TT]-surf[tt]).",
        "With all methods, the RDF can also be calculated around axes parallel",
        "to the [IT]z[it]-axis with option [TT]-xy[tt].",
        "With option [TT]-surf[tt] normalization can not be used.[PAR]",
        "For [TT]-surf[tt], the selection provided to [TT]-refsel[tt] must select",
        "atoms, i.e., centers of mass are not supported.",
        //"If a run input file is supplied ([TT]-s[tt]) and [TT]-rdf[tt] is set",
        //"to [TT]atom[tt], exclusions defined",
        //"in that file are taken into account when calculating the RDF.",
        //"The option [TT]-cut[tt] is meant as an alternative way to avoid",
        //"intramolecular peaks in the RDF plot.",
        //"It is however better to supply a run input file with a higher number of",
        //"exclusions. For e.g. benzene a topology, setting nrexcl to 5",
        //"would eliminate all intramolecular contributions to the RDF.",
        "Note that all atoms in the selected groups are used, also the ones",
        "that don't have Lennard-Jones interactions.[PAR]",
        "Option [TT]-cn[tt] produces the cumulative number RDF,",
        "i.e. the average number of particles within a distance r.[PAR]"
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile().required()
                           .store(&fnRdf_).defaultBasename("rdf")
                           .description("Computed RDFs"));
    options->addOption(FileNameOption("cn").filetype(eftPlot).outputFile()
                           .store(&fnCumulative_).defaultBasename("rdf_cn")
                           .description("Cumulative RDFs"));

    options->addOption(DoubleOption("bin").store(&binwidth_)
                           .description("Bin width (nm)"));
    options->addOption(BooleanOption("norm").store(&bNormalize_)
                           .description("Normalize for volume and density"));
    options->addOption(BooleanOption("xy").store(&bXY_)
                           .description("Use only the x andy components of the distance"));
    options->addOption(DoubleOption("cut").store(&cutoff_)
                           .description("Shortest distance (nm) to be considered"));

    const char *const cSurfaceEnum[] = { "no", "mol", "res" };
    options->addOption(StringOption("surf").enumValue(cSurfaceEnum)
                           .defaultEnumIndex(0).store(&surface_)
                           .description("RDF with respect to the surface of the reference"));

    options->addOption(SelectionOption("ref").store(&refSel_).required()
                           .description("Reference selection for RDF computation"));
    options->addOption(SelectionOption("sel").storeVector(&sel_)
                           .required().multiValue()
                           .description("Selections to compute RDFs for from the reference"));
}

void
Rdf::optionsFinished(Options *options, TrajectoryAnalysisSettings *settings)
{
    if (surface_ != "none")
    {
        settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

        if (options->isSet("norm") && bNormalize_)
        {
            GMX_THROW(InconsistentInputError("-surf cannot be combined with -norm"));
        }
        bNormalize_ = false;
    }
    if (cutoff_ < 0.0)
    {
        cutoff_ = 0.0;
    }
}

void
Rdf::initAnalysis(const TrajectoryAnalysisSettings &settings,
                  const TopologyInformation        &top)
{
    pairDist_.setDataSetCount(sel_.size());
    for (size_t i = 0; i < sel_.size(); ++i)
    {
        pairDist_.setColumnCount(i, 1);
    }
    pairCountSum_.resize(sel_.size());
    plotSettings_ = settings.plotSettings();
    nb_.setXYMode(bXY_);

    const bool bSurface = (surface_ != "none");
    if (bSurface)
    {
        if (!refSel_.hasOnlyAtoms())
        {
            GMX_THROW(InconsistentInputError("-surf only works with -refsel that consists of atoms"));
        }
        // TODO: Make a reusable function for this (indexutil.cpp has something
        // like this already).
        const bool         bMol       = (surface_ == "mol");
        // TODO: Check that molecule information is present.
        const t_topology  *topology   = top.topology();
        int                groupCount = -1;
        int                prevId     = -1;
        int                currId     = -1;
        ConstArrayRef<int> atoms      = refSel_.atomIndices();
        for (int i = 0; i < refSel_.posCount(); ++i)
        {
            const int atomId = atoms[i];
            if (bMol)
            {
                while (atomId >= topology->mols.index[currId + 1])
                {
                    ++currId;
                }
            }
            else
            {
                currId = topology->atoms.atom[atomId].resind;
            }
            if (currId != prevId)
            {
                ++groupCount;
                prevId = currId;
            }
            refSel_.setOriginalId(i, groupCount);
        }
        ++groupCount;
        surfaceDist2_.resize(groupCount);
    }

    /* TODO: Check that these are covered by the framework
       if (fnTPS)
       {
        if (natoms > top->atoms.nr)
        {
            gmx_fatal(FARGS, "Trajectory (%d atoms) does not match topology (%d atoms)",
                      natoms, top->atoms.nr);
        }
       }
       for (i = 0; i < ng+1; i++)
       {
        for (j = 0; j < isize[i]; j++)
        {
            if (index[i][j] >= natoms)
            {
                gmx_fatal(FARGS, "Atom index (%d) in index group %s (%d atoms) larger "
                          "than number of atoms in trajectory (%d atoms)",
                          index[i][j], grpname[i], isize[i], natoms);
            }
        }
       }
     */

    // TODO: Implement exclusions
}

void
Rdf::initAfterFirstFrame(const TrajectoryAnalysisSettings &settings,
                         const t_trxframe                 &fr)
{
    // TODO: Add an -rmax option.
    matrix box;
    copy_mat(fr.box, box);
    if (settings.hasPBC())
    {
        if (bXY_)
        {
            box[ZZ][ZZ] = 2*std::max(box[XX][XX], box[YY][YY]);
        }
        rmax2_ = 0.99*0.99*max_cutoff2(bXY_ ? epbcXY : epbcXYZ, box);
    }
    else
    {
        if (bXY_)
        {
            clear_rvec(box[ZZ]);
        }
        rmax2_ = sqr(3*std::max(box[XX][XX], std::max(box[YY][YY], box[ZZ][ZZ])));
    }
    cut2_  = sqr(cutoff_);
    // We use the double amount of bins, so we can correctly
    // write the rdf and rdf_cn output at i*binwidth values.
    rdf_->init(histogramFromRange(0.0, sqrt(rmax2_)).binWidth(binwidth_ / 2.0));
}

void
Rdf::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                  TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle   dh       = pdata->dataHandle(pairDist_);
    const Selection     &refSel   = pdata->parallelSelection(refSel_);
    const SelectionList &sel      = pdata->parallelSelections(sel_);
    const bool           bSurface = !surfaceDist2_.empty();

    matrix               boxForVolume;
    copy_mat(fr.box, boxForVolume);
    if (bXY_)
    {
        // Set z-size to 1 so we get the surface are iso the volume
        clear_rvec(boxForVolume[ZZ]);
        boxForVolume[ZZ][ZZ] = 1;
    }
    const real invvol = 1.0 / det(boxForVolume);
    inverseVolumeSum_ += invvol;

    dh.startFrame(frnr, fr.time);
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refSel);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        dh.selectDataSet(g);

        if (bSurface)
        {
            for (int i = 0; i < sel[g].posCount(); ++i)
            {
                std::fill(surfaceDist2_.begin(), surfaceDist2_.end(),
                          std::numeric_limits<real>::max());
                AnalysisNeighborhoodPairSearch pairSearch =
                    nbsearch.startPairSearch(sel[g].position(i));
                AnalysisNeighborhoodPair       pair;
                while (pairSearch.findNextPair(&pair))
                {
                    const real r2    = pair.distance2();
                    const int  refId = refSel.position(pair.refIndex()).mappedId();
                    if (r2 < surfaceDist2_[refId])
                    {
                        surfaceDist2_[refId] = r2;
                    }
                }
                for (size_t i = 0; i < surfaceDist2_.size(); ++i)
                {
                    const real r2 = surfaceDist2_[i];
                    if (r2 > cut2_ && r2 <= rmax2_)
                    {
                        dh.setPoint(0, std::sqrt(r2));
                        dh.finishPointSet();
                    }
                }
            }
        }
        else
        {
            AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(sel[g]);
            AnalysisNeighborhoodPair       pair;
            while (pairSearch.findNextPair(&pair))
            {
                const real r2 = pair.distance2();
                if (r2 > cut2_ && r2 <= rmax2_)
                {
                    // TODO: Consider whether the histogramming could be done with
                    // less overhead (after first measuring the overhead).
                    // TODO: Single-precision accumulation may not be sufficient in
                    // all cases.
                    dh.setPoint(0, std::sqrt(r2));
                    dh.finishPointSet();
                }
            }
            pairCountSum_[g] += refSel.posCount() * sel[g].posCount();
        }
    }
    refCountSum_ += refSel.posCount();
    dh.finishFrame();
}

void
Rdf::finishAnalysis(int nframes)
{
    /* Average volume */
    const real invvol = inverseVolumeSum_ / nframes;
    refCountSum_ /= nframes;

    // TODO: Consider how these could be exposed to the testing framework
    // through the dataset registration mechanism.
    AverageHistogramPointer finalRdf
        = rdf_->averager().resampleDoubleBinWidth(true);
    const int               nbin = finalRdf->settings().binCount();

    /* Calculate volume of sphere segments or length of circle segments */
    if (bNormalize_)
    {
        std::vector<real> invBinVolume;
        invBinVolume.resize(nbin);
        real              prevSphereVolume = 0.0;
        for (int i = 0; i < nbin; ++i)
        {
            const real r = (i + 0.5)*binwidth_;
            real       sphereVolume;
            if (bXY_)
            {
                sphereVolume = M_PI*r*r;
            }
            else
            {
                sphereVolume = (4.0/3.0)*M_PI*r*r*r;
            }
            const real binVolume = sphereVolume - prevSphereVolume;
            invBinVolume[i]  = 1.0 / binVolume;
            prevSphereVolume = sphereVolume;
        }
        finalRdf->scaleAllByVector(&invBinVolume[0]);
    }

    for (size_t g = 0; g < sel_.size(); ++g)
    {
        real normfac;
        if (bNormalize_)
        {
            normfac = 1.0 / (invvol * pairCountSum_[g] / nframes);
        }
        else
        {
            normfac = 1.0 / (binwidth_ * refCountSum_);
        }
        finalRdf->scaleSingle(g, normfac);
    }

    // TODO: Consider if some of this should be done in writeOutput().
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(plotSettings_));
        plotm->setFileName(fnRdf_);
        plotm->setTitle("Radial distribution");
        // TODO: Add an overload.
        plotm->setSubtitle(formatString("reference %s", refSel_.name()).c_str());
        plotm->setXLabel("r (nm)");
        for (size_t i = 0; i < sel_.size(); ++i)
        {
            plotm->appendLegend(sel_[i].name());
        }
        finalRdf->addModule(plotm);
    }
    finalRdf->done();

    if (!fnCumulative_.empty())
    {
        AverageHistogramPointer cumulativeRdf
            = rdf_->averager().resampleDoubleBinWidth(false);
        cumulativeRdf->scaleAll(1.0 / refCountSum_);
        cumulativeRdf->makeCumulative();

        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(plotSettings_));
        plotm->setFileName(fnRdf_);
        plotm->setTitle("Cumulative Number RDF");
        // TODO: Add an overload.
        plotm->setSubtitle(formatString("reference %s", refSel_.name()).c_str());
        plotm->setXLabel("r (nm)");
        plotm->setYLabel("number");
        for (size_t i = 0; i < sel_.size(); ++i)
        {
            plotm->appendLegend(sel_[i].name());
        }
        cumulativeRdf->addModule(plotm);

        cumulativeRdf->done();
    }
}

void
Rdf::writeOutput()
{
}

//! \}

}       // namespace

const char RdfInfo::name[]             = "rdf_new";
const char RdfInfo::shortDescription[] =
    "Calculate radial distribution functions";

TrajectoryAnalysisModulePointer RdfInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Rdf);
}

} // namespace analysismodules

} // namespace gmx
