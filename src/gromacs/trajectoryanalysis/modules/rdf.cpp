/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * Implements gmx::analysismodules::Rdf.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com> (C++ conversion)
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "rdf.h"

#include <cmath>
#include <cstddef>

#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class AnalysisDataParallelOptions;
class SelectionCollection;

namespace analysismodules
{

namespace
{

//! \addtogroup module_trajectoryanalysis
//! \{

/********************************************************************
 * Actual analysis module
 */

//! Normalization for the computed distribution.
enum class Normalization : int
{
    Rdf,
    NumberDensity,
    None,
    Count
};
//! String values corresponding to Normalization.
const EnumerationArray<Normalization, const char*> c_normalizationNames = {
    { "rdf", "number_density", "none" }
};
//! Whether to compute RDF wrt. surface of the reference group.
enum class SurfaceType : int
{
    None,
    Molecule,
    Residue,
    Count
};
//! String values corresponding to SurfaceType.
const EnumerationArray<SurfaceType, const char*> c_surfaceTypeNames = { { "no", "mol", "res" } };

/*! \brief
 * Implements `gmx rdf` trajectory analysis module.
 */
class Rdf : public TrajectoryAnalysisModule
{
public:
    Rdf();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void initAfterFirstFrame(const TrajectoryAnalysisSettings& settings, const t_trxframe& fr) override;

    TrajectoryAnalysisModuleDataPointer startFrames(const AnalysisDataParallelOptions& opt,
                                                    const SelectionCollection& selections) override;
    void                                analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    std::string              fnRdf_;
    std::string              fnCumulative_;
    SurfaceType              surface_;
    AnalysisDataPlotSettings plotSettings_;

    /*! \brief
     * Reference selection to compute RDFs around.
     *
     * With -surf, Selection::originalIds() and Selection::mappedIds()
     * store the index of the surface group to which that position belongs.
     * The RDF is computed by finding the nearest position from each
     * surface group for each position, and then binning those distances.
     */
    Selection refSel_;
    /*! \brief
     * Selections to compute RDFs for.
     */
    SelectionList sel_;

    /*! \brief
     * Raw pairwise distance data from which the RDF is computed.
     *
     * There is a data set for each selection in `sel_`, with a single
     * column.  Each point set will contain a single pairwise distance
     * that contributes to the RDF.
     */
    AnalysisData pairDist_;
    /*! \brief
     * Normalization factors for each frame.
     *
     * The first column contains the number of positions in `refSel_` for
     * that frame (with surface RDF, the number of groups).  There are
     * `sel_.size()` more columns, each containing the number density of
     * positions for one selection.
     */
    AnalysisData normFactors_;
    /*! \brief
     * Histogram module that computes the actual RDF from `pairDist_`.
     *
     * The per-frame histograms are raw pair counts in each bin;
     * the averager is normalized by the average number of reference
     * positions (average of the first column of `normFactors_`).
     */
    AnalysisDataSimpleHistogramModulePointer pairCounts_;
    /*! \brief
     * Average normalization factors.
     */
    AnalysisDataAverageModulePointer normAve_;
    //! Neighborhood search with `refSel_` as the reference positions.
    AnalysisNeighborhood nb_;
    //! Topology exclusions used by neighborhood searching.
    const gmx_localtop_t* localTop_;

    // User input options.
    double        binwidth_;
    double        cutoff_;
    double        rmax_;
    Normalization normalization_;
    bool          bNormalizationSet_;
    bool          bXY_;
    bool          bExclusions_;

    // Pre-computed values for faster access during analysis.
    real cut2_;
    real rmax2_;
    int  surfaceGroupCount_;

    // Copy and assign disallowed by base.
};

Rdf::Rdf() :
    surface_(SurfaceType::None),
    pairCounts_(new AnalysisDataSimpleHistogramModule()),
    normAve_(new AnalysisDataAverageModule()),
    localTop_(nullptr),
    binwidth_(0.002),
    cutoff_(0.0),
    rmax_(0.0),
    normalization_(Normalization::Rdf),
    bNormalizationSet_(false),
    bXY_(false),
    bExclusions_(false),
    cut2_(0.0),
    rmax2_(0.0),
    surfaceGroupCount_(0)
{
    pairDist_.setMultipoint(true);
    pairDist_.addModule(pairCounts_);
    registerAnalysisDataset(&pairDist_, "pairdist");
    registerBasicDataset(pairCounts_.get(), "paircount");

    normFactors_.addModule(normAve_);
    registerAnalysisDataset(&normFactors_, "norm");
}

void Rdf::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    const char* const desc[] = {
        "[THISMODULE] calculates radial distribution functions from one",
        "reference set of position (set with [TT]-ref[tt]) to one or more",
        "sets of positions (set with [TT]-sel[tt]).  To compute the RDF with",
        "respect to the closest position in a set in [TT]-ref[tt] instead, use",
        "[TT]-surf[tt]: if set, then [TT]-ref[tt] is partitioned into sets",
        "based on the value of [TT]-surf[tt], and the closest position in each",
        "set is used. To compute the RDF around axes parallel to the",
        "[IT]z[it]-axis, i.e., only in the [IT]x[it]-[IT]y[it] plane, use",
        "[TT]-xy[tt].",
        "",
        "To set the bin width and maximum distance to use in the RDF, use",
        "[TT]-bin[tt] and [TT]-rmax[tt], respectively. The latter can be",
        "used to limit the computational cost if the RDF is not of interest",
        "up to the default (half of the box size with PBC, three times the",
        "box size without PBC).",
        "",
        "To use exclusions from the topology ([TT]-s[tt]), set [TT]-excl[tt]",
        "and ensure that both [TT]-ref[tt] and [TT]-sel[tt] only select atoms.",
        "A rougher alternative to exclude intra-molecular peaks is to set",
        "[TT]-cut[tt] to a non-zero value to clear the RDF at small",
        "distances.",
        "",
        "The RDFs are normalized by 1) average number of positions in",
        "[TT]-ref[tt] (the number of groups with [TT]-surf[tt]), 2) volume",
        "of the bin, and 3) average particle density of [TT]-sel[tt] positions",
        "for that selection. To change the normalization, use [TT]-norm[tt]:",
        "",
        "* [TT]rdf[tt]: Use all factors for normalization.",
        "  This produces a normal RDF.",
        "* [TT]number_density[tt]: Use the first two factors.",
        "  This produces a number density as a function of distance.",
        "* [TT]none[tt]: Use only the first factor.",
        "  In this case, the RDF is only scaled with the bin width to make",
        "  the integral of the curve represent the number of pairs within a",
        "  range.",
        "",
        "Note that exclusions do not affect the normalization: even if",
        "[TT]-excl[tt] is set, or [TT]-ref[tt] and",
        "[TT]-sel[tt] contain the same selection, the normalization factor",
        "is still N*M, not N*(M-excluded).",
        "",
        "For [TT]-surf[tt], the selection provided to [TT]-ref[tt] must",
        "select atoms, i.e., centers of mass are not supported. Further,",
        "[TT]-nonorm[tt] is implied, as the bins have irregular shapes and",
        "the volume of a bin is not easily computable.",
        "",
        "Option [TT]-cn[tt] produces the cumulative number RDF,",
        "i.e. the average number of particles within a distance r."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .required()
                               .store(&fnRdf_)
                               .defaultBasename("rdf")
                               .description("Computed RDFs"));
    options->addOption(FileNameOption("cn")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnCumulative_)
                               .defaultBasename("rdf_cn")
                               .description("Cumulative RDFs"));

    options->addOption(DoubleOption("bin").store(&binwidth_).description("Bin width (nm)"));
    options->addOption(EnumOption<Normalization>("norm")
                               .enumValue(c_normalizationNames)
                               .store(&normalization_)
                               .storeIsSet(&bNormalizationSet_)
                               .description("Normalization"));
    options->addOption(BooleanOption("xy").store(&bXY_).description(
            "Use only the x and y components of the distance"));
    options->addOption(
            BooleanOption("excl").store(&bExclusions_).description("Use exclusions from topology"));
    options->addOption(DoubleOption("cut").store(&cutoff_).description(
            "Shortest distance (nm) to be considered"));
    options->addOption(
            DoubleOption("rmax").store(&rmax_).description("Largest distance (nm) to calculate"));

    options->addOption(EnumOption<SurfaceType>("surf")
                               .enumValue(c_surfaceTypeNames)
                               .store(&surface_)
                               .description("RDF with respect to the surface of the reference"));

    options->addOption(SelectionOption("ref").store(&refSel_).required().description(
            "Reference selection for RDF computation"));
    options->addOption(SelectionOption("sel").storeVector(&sel_).required().multiValue().description(
            "Selections to compute RDFs for from the reference"));
}

void Rdf::optionsFinished(TrajectoryAnalysisSettings* settings)
{
    if (surface_ != SurfaceType::None)
    {
        settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

        if (bNormalizationSet_ && normalization_ != Normalization::None)
        {
            GMX_THROW(InconsistentInputError("-surf cannot be combined with -norm"));
        }
        normalization_ = Normalization::None;
        if (bExclusions_)
        {
            GMX_THROW(InconsistentInputError("-surf cannot be combined with -excl"));
        }
    }
    if (bExclusions_)
    {
        settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    }
    if (cutoff_ < 0.0)
    {
        cutoff_ = 0.0;
    }
}

void Rdf::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    pairDist_.setDataSetCount(sel_.size());
    for (size_t i = 0; i < sel_.size(); ++i)
    {
        pairDist_.setColumnCount(i, 1);
    }
    plotSettings_ = settings.plotSettings();
    nb_.setXYMode(bXY_);

    normFactors_.setColumnCount(0, sel_.size() + 1);

    const bool bSurface = (surface_ != SurfaceType::None);
    if (bSurface)
    {
        if (!refSel_.hasOnlyAtoms())
        {
            GMX_THROW(InconsistentInputError("-surf only works with -ref that consists of atoms"));
        }
        const e_index_t type = (surface_ == SurfaceType::Molecule ? INDEX_MOL : INDEX_RES);
        surfaceGroupCount_   = refSel_.initOriginalIdsToGroup(top.mtop(), type);
    }

    if (bExclusions_)
    {
        if (!refSel_.hasOnlyAtoms() || !refSel_.hasSortedAtomIndices())
        {
            GMX_THROW(
                    InconsistentInputError("-excl only works with a -ref selection that consist of "
                                           "atoms in ascending (sorted) order"));
        }
        for (size_t i = 0; i < sel_.size(); ++i)
        {
            if (!sel_[i].hasOnlyAtoms())
            {
                GMX_THROW(InconsistentInputError(
                        "-excl only works with selections that consist of atoms"));
            }
        }
        localTop_ = top.expandedTopology();
        if (localTop_->excls.empty())
        {
            GMX_THROW(InconsistentInputError(
                    "-excl is set, but the file provided to -s does not define exclusions"));
        }
        nb_.setTopologyExclusions(&localTop_->excls);
    }
}

void Rdf::initAfterFirstFrame(const TrajectoryAnalysisSettings& settings, const t_trxframe& fr)
{
    // If -rmax is not provided, determine one from the box for the first frame.
    if (rmax_ <= 0.0)
    {
        matrix box;
        copy_mat(fr.box, box);
        if (settings.hasPBC())
        {
            if (bXY_)
            {
                box[ZZ][ZZ] = 2 * std::max(box[XX][XX], box[YY][YY]);
            }
            rmax_ = std::sqrt(0.99 * 0.99 * max_cutoff2(bXY_ ? PbcType::XY : PbcType::Xyz, box));
        }
        else
        {
            if (bXY_)
            {
                clear_rvec(box[ZZ]);
            }
            rmax_ = 3 * std::max(box[XX][XX], std::max(box[YY][YY], box[ZZ][ZZ]));
        }
    }
    cut2_  = gmx::square(cutoff_);
    rmax2_ = gmx::square(rmax_);
    nb_.setCutoff(rmax_);
    // We use the double amount of bins, so we can correctly
    // write the rdf and rdf_cn output at i*binwidth values.
    pairCounts_->init(histogramFromRange(0.0, rmax_).binWidth(binwidth_ / 2.0));
}

/*! \brief
 * Temporary memory for use within a single-frame calculation.
 */
class RdfModuleData : public TrajectoryAnalysisModuleData
{
public:
    /*! \brief
     * Reserves memory for the frame-local data.
     *
     * `surfaceGroupCount` will be zero if -surf is not specified.
     */
    RdfModuleData(TrajectoryAnalysisModule*          module,
                  const AnalysisDataParallelOptions& opt,
                  const SelectionCollection&         selections,
                  int                                surfaceGroupCount) :
        TrajectoryAnalysisModuleData(module, opt, selections)
    {
        surfaceDist2_.resize(surfaceGroupCount);
    }

    void finish() override { finishDataHandles(); }

    /*! \brief
     * Minimum distance to each surface group.
     *
     * One entry for each group (residue/molecule, per -surf) in the
     * reference selection.
     * This is needed to support neighborhood searching, which may not
     * return the reference positions in order: for each position, we need
     * to search through all the reference positions and update this array
     * to find the minimum distance to each surface group, and then compute
     * the RDF from these numbers.
     */
    std::vector<real> surfaceDist2_;
};

TrajectoryAnalysisModuleDataPointer Rdf::startFrames(const AnalysisDataParallelOptions& opt,
                                                     const SelectionCollection&         selections)
{
    return TrajectoryAnalysisModuleDataPointer(new RdfModuleData(this, opt, selections, surfaceGroupCount_));
}

void Rdf::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   dh        = pdata->dataHandle(pairDist_);
    AnalysisDataHandle   nh        = pdata->dataHandle(normFactors_);
    const Selection&     refSel    = TrajectoryAnalysisModuleData::parallelSelection(refSel_);
    const SelectionList& sel       = TrajectoryAnalysisModuleData::parallelSelections(sel_);
    RdfModuleData&       frameData = *static_cast<RdfModuleData*>(pdata);
    const bool           bSurface  = !frameData.surfaceDist2_.empty();

    matrix boxForVolume;
    copy_mat(fr.box, boxForVolume);
    if (bXY_)
    {
        // Set z-size to 1 so we get the surface are iso the volume
        clear_rvec(boxForVolume[ZZ]);
        boxForVolume[ZZ][ZZ] = 1;
    }
    const real inverseVolume = 1.0 / det(boxForVolume);

    nh.startFrame(frnr, fr.time);
    // Compute the normalization factor for the number of reference positions.
    if (bSurface)
    {
        if (refSel.isDynamic())
        {
            // Count the number of distinct groups.
            // This assumes that each group is continuous, which is currently
            // the case.
            int count  = 0;
            int prevId = -1;
            for (int i = 0; i < refSel.posCount(); ++i)
            {
                const int id = refSel.position(i).mappedId();
                if (id != prevId)
                {
                    ++count;
                    prevId = id;
                }
            }
            nh.setPoint(0, count);
        }
        else
        {
            nh.setPoint(0, surfaceGroupCount_);
        }
    }
    else
    {
        nh.setPoint(0, refSel.posCount());
    }

    dh.startFrame(frnr, fr.time);
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, refSel);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        dh.selectDataSet(g);

        if (bSurface)
        {
            // Special loop for surface calculation, where a separate neighbor
            // search is done for each position in the selection, and the
            // nearest position from each surface group is tracked.
            std::vector<real>& surfaceDist2 = frameData.surfaceDist2_;
            for (int i = 0; i < sel[g].posCount(); ++i)
            {
                std::fill(surfaceDist2.begin(), surfaceDist2.end(), std::numeric_limits<real>::max());
                AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(sel[g].position(i));
                AnalysisNeighborhoodPair       pair;
                while (pairSearch.findNextPair(&pair))
                {
                    const real r2    = pair.distance2();
                    const int  refId = refSel.position(pair.refIndex()).mappedId();
                    if (r2 < surfaceDist2[refId])
                    {
                        surfaceDist2[refId] = r2;
                    }
                }
                // Accumulate the RDF from the distances to the surface.
                for (size_t i = 0; i < surfaceDist2.size(); ++i)
                {
                    const real r2 = surfaceDist2[i];
                    // Here, we need to check for rmax, since the value might
                    // be above the cutoff if no points were close to some
                    // surface positions.
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
            // Standard neighborhood search over all pairs within the cutoff
            // for the -surf no case.
            AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(sel[g]);
            AnalysisNeighborhoodPair       pair;
            while (pairSearch.findNextPair(&pair))
            {
                const real r2 = pair.distance2();
                if (r2 > cut2_)
                {
                    // TODO: Consider whether the histogramming could be done with
                    // less overhead (after first measuring the overhead).
                    dh.setPoint(0, std::sqrt(r2));
                    dh.finishPointSet();
                }
            }
        }
        // Normalization factor for the number density (only used without
        // -surf, but does not hurt to populate otherwise).
        nh.setPoint(g + 1, sel[g].posCount() * inverseVolume);
    }
    dh.finishFrame();
    nh.finishFrame();
}

void Rdf::finishAnalysis(int /*nframes*/)
{
    // Normalize the averager with the number of reference positions,
    // from where the normalization propagates to all the output.
    const real refPosCount = normAve_->average(0, 0);
    pairCounts_->averager().scaleAll(1.0 / refPosCount);
    pairCounts_->averager().done();

    // TODO: Consider how these could be exposed to the testing framework
    // through the dataset registration mechanism.
    AverageHistogramPointer finalRdf = pairCounts_->averager().resampleDoubleBinWidth(true);

    if (normalization_ != Normalization::None)
    {
        // Normalize by the volume of the bins (volume of sphere segments or
        // length of circle segments).
        std::vector<real> invBinVolume;
        const int         nbin = finalRdf->settings().binCount();
        invBinVolume.resize(nbin);
        real prevSphereVolume = 0.0;
        for (int i = 0; i < nbin; ++i)
        {
            const real r            = (i + 0.5) * binwidth_;
            const real sphereVolume = (bXY_) ? M_PI * r * r : (4.0 / 3.0) * M_PI * r * r * r;
            const real binVolume    = sphereVolume - prevSphereVolume;
            invBinVolume[i]         = 1.0 / binVolume;
            prevSphereVolume        = sphereVolume;
        }
        finalRdf->scaleAllByVector(invBinVolume.data());

        if (normalization_ == Normalization::Rdf)
        {
            // Normalize by particle density.
            for (size_t g = 0; g < sel_.size(); ++g)
            {
                finalRdf->scaleSingle(g, 1.0 / normAve_->average(0, g + 1));
            }
        }
    }
    else
    {
        // With no normalization, just scale with bin width to make the
        // integral of the curve (instead of raw bin sum) represent the pair
        // count.
        finalRdf->scaleAll(1.0 / binwidth_);
    }
    finalRdf->done();

    // TODO: Consider if some of this should be done in writeOutput().
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(plotSettings_));
        plotm->setFileName(fnRdf_);
        plotm->setTitle("Radial distribution");
        plotm->setSubtitle(formatString("reference %s", refSel_.name()));
        plotm->setXLabel("r (nm)");
        plotm->setYLabel("g(r)");
        plotm->setXFormat(11, 6);
        plotm->setYFormat(11, 6);
        for (size_t i = 0; i < sel_.size(); ++i)
        {
            plotm->appendLegend(sel_[i].name());
        }
        finalRdf->addModule(plotm);
    }

    if (!fnCumulative_.empty())
    {
        AverageHistogramPointer cumulativeRdf = pairCounts_->averager().resampleDoubleBinWidth(false);
        cumulativeRdf->makeCumulative();
        cumulativeRdf->done();

        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(plotSettings_));
        plotm->setFileName(fnCumulative_);
        plotm->setTitle("Cumulative Number RDF");
        plotm->setSubtitle(formatString("reference %s", refSel_.name()));
        plotm->setXLabel("r (nm)");
        plotm->setYLabel("number");
        for (size_t i = 0; i < sel_.size(); ++i)
        {
            plotm->appendLegend(sel_[i].name());
        }
        cumulativeRdf->addModule(plotm);
    }
}

void Rdf::writeOutput() {}

//! \}

} // namespace

const char RdfInfo::name[]             = "rdf";
const char RdfInfo::shortDescription[] = "Calculate radial distribution functions";

TrajectoryAnalysisModulePointer RdfInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Rdf);
}

} // namespace analysismodules

} // namespace gmx
