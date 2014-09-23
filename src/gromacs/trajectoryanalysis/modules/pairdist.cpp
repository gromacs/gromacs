/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::PairDistance.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "pairdist.h"

#include <cmath>

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/trx.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
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

//! Enum value to store the selected value for `-type`.
enum DistanceType
{
    eDistanceType_Min,
    eDistanceType_Max
};

//! Enum value to store the selected value for `-refgrouping`/`-selgrouping`.
enum GroupType
{
    eGroupType_All,
    eGroupType_Residue,
    eGroupType_Molecule,
    eGroupType_None
};

//! Strings corresponding to DistanceType.
const char *const           c_distanceTypes[] = { "min", "max" };
//! Strings corresponding to GroupType.
const char *const           c_groupTypes[]    = { "all", "res", "mol", "none" };

/*! \brief
 * Implements `gmx pairdist` trajectory analysis module.
 */
class PairDistance : public TrajectoryAnalysisModule
{
    public:
        PairDistance();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual TrajectoryAnalysisModuleDataPointer startFrames(
            const AnalysisDataParallelOptions &opt,
            const SelectionCollection         &selections);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        /*! \brief
         * Computed distances as a function of time.
         *
         * There is one data set for each selection in `sel_`.
         * Within each data set, there is one column for each distance to be
         * computed, as explained in the `-h` text.
         */
        AnalysisData            distances_;

        /*! \brief
         * Reference selection to compute distances to.
         *
         * mappedId() identifies the group (of type `refGroupType_`) into which
         * each position belogs.
         */
        Selection               refSel_;
        /*! \brief
         * Selections to compute distances from.
         *
         * mappedId() identifies the group (of type `selGroupType_`) into which
         * each position belogs.
         */
        SelectionList           sel_;

        std::string             fnDist_;

        double                  cutoff_;
        int                     distanceType_;
        int                     refGroupType_;
        int                     selGroupType_;

        //! Number of groups in `refSel_`.
        int                     refGroupCount_;
        //! Maximum number of pairs of groups for one selection.
        int                     maxGroupCount_;
        //! Initial squared distance for distance accumulation.
        real                    initialDist2_;
        //! Cutoff squared for use in the actual calculation.
        real                    cutoff2_;

        //! Neighborhood search object for the pair search.
        AnalysisNeighborhood    nb_;

        // Copy and assign disallowed by base.
};

PairDistance::PairDistance()
    : TrajectoryAnalysisModule(PairDistanceInfo::name, PairDistanceInfo::shortDescription),
      cutoff_(0.0), distanceType_(eDistanceType_Min),
      refGroupType_(eGroupType_All), selGroupType_(eGroupType_All),
      refGroupCount_(0), maxGroupCount_(0), initialDist2_(0.0), cutoff2_(0.0)
{
    registerAnalysisDataset(&distances_, "dist");
}


void
PairDistance::initOptions(Options *options, TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates pairwise distances between one reference",
        "selection (given with [TT]-ref[tt]) and one or more other selections",
        "(given with [TT]-sel[tt]).  It can calculate either the minimum",
        "distance (the default), or the maximum distance (with",
        "[TT]-type max[tt]).  Distances to each selection provided with",
        "[TT]-sel[tt] are computed independently.[PAR]",
        "By default, the global minimum/maximum distance is computed.",
        "To compute more distances (e.g., minimum distances to each residue",
        "in [TT]-ref[tt]), use [TT]-refgrouping[tt] and/or [TT]-selgrouping[tt]",
        "to specify how the positions within each selection should be",
        "grouped.[PAR]",
        "Computed distances are written to the file specified with [TT]-o[tt].",
        "If there are N groups in [TT]-ref[tt] and M groups in the first",
        "selection in [TT]-sel[tt], then the output contains N*M columns",
        "for the first selection. The columns contain distances like this:",
        "r1-s1, r2-s1, ..., r1-s2, r2-s2, ..., where rn is the n'th group",
        "in [TT]-ref[tt] and sn is the n'th group in the other selection.",
        "The distances for the second selection comes as separate columns",
        "after the first selection, and so on.  If some selections are",
        "dynamic, only the selected positions are used in the computation",
        "but the same number of columns is always written out.  If there",
        "are no positions contributing to some group pair, then the cutoff",
        "value is written (see below).[PAR]",
        "[TT]-cutoff[tt] sets a cutoff for the computed distances.",
        "If the result would contain a distance over the cutoff, the cutoff",
        "value is written to the output file instead. By default, no cutoff",
        "is used, but if you are not interested in values beyond a cutoff,",
        "or if you know that the minimum distance is smaller than a cutoff,",
        "you should set this option to allow the tool to use grid-based",
        "searching and be significantly faster.[PAR]",
        "If you want to compute distances between fixed pairs,",
        "[gmx-distance] may be a more suitable tool."
    };

    options->setDescription(desc);

    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile().required()
                           .store(&fnDist_).defaultBasename("dist")
                           .description("Distances as function of time"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                           .description("Maximum distance to consider"));
    options->addOption(StringOption("type").storeEnumIndex(&distanceType_)
                           .defaultEnumIndex(0).enumValue(c_distanceTypes)
                           .description("Type of distances to calculate"));
    options->addOption(StringOption("refgrouping").storeEnumIndex(&refGroupType_)
                           .defaultEnumIndex(0).enumValue(c_groupTypes)
                           .description("Grouping of -ref positions to compute the min/max over"));
    options->addOption(StringOption("selgrouping").storeEnumIndex(&selGroupType_)
                           .defaultEnumIndex(0).enumValue(c_groupTypes)
                           .description("Grouping of -sel positions to compute the min/max over"));

    options->addOption(SelectionOption("ref").store(&refSel_).required()
                           .description("Reference positions to calculate distances from"));
    options->addOption(SelectionOption("sel").storeVector(&sel_).required().multiValue()
                           .description("Positions to calculate distances for"));
}

//! Helper function to initialize the grouping for a selection.
int initSelectionGroups(Selection *sel, t_topology *top, int type)
{
    e_index_t indexType = INDEX_UNKNOWN;
    switch (type)
    {
        case eGroupType_All:      indexType = INDEX_ALL; break;
        case eGroupType_Residue:  indexType = INDEX_RES; break;
        case eGroupType_Molecule: indexType = INDEX_MOL; break;
        case eGroupType_None:     indexType = INDEX_ATOM; break;
    }
    return sel->initOriginalIdsToGroup(top, indexType);
}


void
PairDistance::initAnalysis(const TrajectoryAnalysisSettings &settings,
                           const TopologyInformation        &top)
{
    refGroupCount_ = initSelectionGroups(&refSel_, top.topology(), refGroupType_);

    maxGroupCount_ = 0;
    distances_.setDataSetCount(sel_.size());
    for (size_t i = 0; i < sel_.size(); ++i)
    {
        const int selGroupCount
            = initSelectionGroups(&sel_[i], top.topology(), selGroupType_);
        const int columnCount = refGroupCount_ * selGroupCount;
        maxGroupCount_ = std::max(maxGroupCount_, columnCount);
        distances_.setColumnCount(i, columnCount);
    }

    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        if (distanceType_ == eDistanceType_Max)
        {
            plotm->setTitle("Maximum distance");
        }
        else
        {
            plotm->setTitle("Minimum distance");
        }
        // TODO: Figure out and add a descriptive subtitle and/or a longer
        // title and/or better legends based on the grouping and the reference
        // selection.
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plotm->appendLegend(sel_[g].name());
        }
        distances_.addModule(plotm);
    }

    nb_.setCutoff(cutoff_);
    if (cutoff_ <= 0.0)
    {
        cutoff_       = 0.0;
        initialDist2_ = std::numeric_limits<real>::max();
    }
    else
    {
        initialDist2_ = cutoff_ * cutoff_;
    }
    if (distanceType_ == eDistanceType_Max)
    {
        initialDist2_ = 0.0;
    }
    cutoff2_ = cutoff_ * cutoff_;
}

/*! \brief
 * Temporary memory for use within a single-frame calculation.
 */
class PairDistanceModuleData : public TrajectoryAnalysisModuleData
{
    public:
        /*! \brief
         * Reserves memory for the frame-local data.
         */
        PairDistanceModuleData(TrajectoryAnalysisModule          *module,
                               const AnalysisDataParallelOptions &opt,
                               const SelectionCollection         &selections,
                               int                                refGroupCount,
                               const Selection                   &refSel,
                               int                                maxGroupCount)
            : TrajectoryAnalysisModuleData(module, opt, selections)
        {
            distArray_.resize(maxGroupCount);
            countArray_.resize(maxGroupCount);
            refCountArray_.resize(refGroupCount);
            if (!refSel.isDynamic())
            {
                initRefCountArray(refSel);
            }
        }

        virtual void finish() { finishDataHandles(); }

        /*! \brief
         * Computes the number of positions in each group in \p refSel
         * and stores them into `refCountArray_`.
         */
        void initRefCountArray(const Selection &refSel)
        {
            std::fill(refCountArray_.begin(), refCountArray_.end(), 0);
            int refPos = 0;
            while (refPos < refSel.posCount())
            {
                const int refIndex = refSel.position(refPos).mappedId();
                const int startPos = refPos;
                ++refPos;
                while (refPos < refSel.posCount()
                       && refSel.position(refPos).mappedId() == refIndex)
                {
                    ++refPos;
                }
                refCountArray_[refIndex] = refPos - startPos;
            }
        }

        /*! \brief
         * Squared distance between each group
         *
         * One entry for each group pair for the current selection.
         * Enough memory is allocated to fit the largest calculation selection.
         * This is needed to support neighborhood searching, which may not
         * return the pairs in order: for each group pair, we need to search
         * through all the position pairs and update this array to find the
         * minimum/maximum distance between them.
         */
        std::vector<real> distArray_;
        /*! \brief
         * Number of pairs within the cutoff that have contributed to the value
         * in `distArray_`.
         *
         * This is needed to identify whether there were any pairs inside the
         * cutoff and whether there were additional pairs outside the cutoff
         * that were not covered by the neihborhood search.
         */
        std::vector<int>  countArray_;
        /*! \brief
         * Number of positions within each reference group.
         *
         * This is used to more efficiently compute the total number of pairs
         * (for comparison with `countArray_`), as otherwise these numbers
         * would need to be recomputed for each selection.
         */
        std::vector<int>  refCountArray_;
};

TrajectoryAnalysisModuleDataPointer PairDistance::startFrames(
        const AnalysisDataParallelOptions &opt,
        const SelectionCollection         &selections)
{
    return TrajectoryAnalysisModuleDataPointer(
            new PairDistanceModuleData(this, opt, selections, refGroupCount_,
                                       refSel_, maxGroupCount_));
}

void
PairDistance::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                           TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle         dh            = pdata->dataHandle(distances_);
    const Selection           &refSel        = pdata->parallelSelection(refSel_);
    const SelectionList       &sel           = pdata->parallelSelections(sel_);
    PairDistanceModuleData    &frameData     = *static_cast<PairDistanceModuleData *>(pdata);
    std::vector<real>         &distArray     = frameData.distArray_;
    std::vector<int>          &countArray    = frameData.countArray_;

    if (cutoff_ > 0.0 && refSel.isDynamic())
    {
        // Count the number of reference positions in each group, so that
        // this does not need to be computed again for each selection.
        // This is needed only if it is possible that the neighborhood search
        // does not cover all the pairs, hence the cutoff > 0.0 check.
        // If refSel is static, then the array contents are static as well,
        // and it has been initialized in the constructor of the data object.
        frameData.initRefCountArray(refSel);
    }
    const std::vector<int>    &refCountArray = frameData.refCountArray_;

    AnalysisNeighborhoodSearch nbsearch  = nb_.initSearch(pbc, refSel);
    dh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        const int columnCount = distances_.columnCount(g);
        std::fill(distArray.begin(), distArray.begin() + columnCount, initialDist2_);
        std::fill(countArray.begin(), countArray.begin() + columnCount, 0);

        // Accumulate the number of position pairs within the cutoff and the
        // min/max distance for each group pair.
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(sel[g]);
        AnalysisNeighborhoodPair       pair;
        while (pairSearch.findNextPair(&pair))
        {
            const SelectionPosition &refPos   = refSel.position(pair.refIndex());
            const SelectionPosition &selPos   = sel[g].position(pair.testIndex());
            const int                refIndex = refPos.mappedId();
            const int                selIndex = selPos.mappedId();
            const int                index    = selIndex * refGroupCount_ + refIndex;
            const real               r2       = pair.distance2();
            if (distanceType_ == eDistanceType_Min)
            {
                if (distArray[index] > r2)
                {
                    distArray[index] = r2;
                }
            }
            else
            {
                if (distArray[index] < r2)
                {
                    distArray[index] = r2;
                }
            }
            ++countArray[index];
        }

        // If it is possible that positions outside the cutoff (or lack of
        // them) affects the result, then we need to check whether there were
        // any.  This is necessary for two cases:
        //  - With max distances, if there are pairs outside the cutoff, then
        //    the computed distance should be equal to the cutoff instead of
        //    the largest distance that was found above.
        //  - With either distance type, if all pairs are outside the cutoff,
        //    then countArray must be updated so that the presence flag
        //    in the output data reflects the dynamic selection status, not
        //    whether something was inside the cutoff or not.
        if (cutoff_ > 0.0)
        {
            int selPos = 0;
            // Loop over groups in this selection (at start, selPos is always
            // the first position in the next group).
            while (selPos < sel[g].posCount())
            {
                // Count the number of positions in this group.
                const int selIndex = sel[g].position(selPos).mappedId();
                const int startPos = selPos;
                ++selPos;
                while (selPos < sel[g].posCount()
                       && sel[g].position(selPos).mappedId() == selIndex)
                {
                    ++selPos;
                }
                const int count = selPos - startPos;
                // Check all group pairs that contain this group.
                for (int i = 0; i < refGroupCount_; ++i)
                {
                    const int index      = selIndex * refGroupCount_ + i;
                    const int totalCount = refCountArray[i] * count;
                    // If there were positions outside the cutoff,
                    // update the distance if necessary and the count.
                    if (countArray[index] < totalCount)
                    {
                        if (distanceType_ == eDistanceType_Max)
                        {
                            distArray[index] = cutoff2_;
                        }
                        countArray[index] = totalCount;
                    }
                }
            }
        }

        // Write the computed distances to the output data.
        dh.selectDataSet(g);
        for (int i = 0; i < columnCount; ++i)
        {
            if (countArray[i] > 0)
            {
                dh.setPoint(i, std::sqrt(distArray[i]));
            }
            else
            {
                // If there are no contributing positions, write out the cutoff
                // value.
                dh.setPoint(i, cutoff_, false);
            }
        }
    }
    dh.finishFrame();
}

void
PairDistance::finishAnalysis(int /*nframes*/)
{
}

void
PairDistance::writeOutput()
{
}

//! \}

}       // namespace

const char PairDistanceInfo::name[]             = "pairdist";
const char PairDistanceInfo::shortDescription[] =
    "Calculate pairwise distances between groups of positions";

TrajectoryAnalysisModulePointer PairDistanceInfo::create()
{
    return TrajectoryAnalysisModulePointer(new PairDistance);
}

} // namespace analysismodules

} // namespace gmx
