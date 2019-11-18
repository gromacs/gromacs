/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::ExtractCluster.
 *
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "extract_cluster.h"

#include <algorithm>

#include "gromacs/coordinateio/coordinatefile.h"
#include "gromacs/coordinateio/requirements.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/index.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*
 * ExtractCluster
 */

class ExtractCluster : public TrajectoryAnalysisModule
{
public:
    ExtractCluster();

    ~ExtractCluster() override;

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int frameNumber, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    //! Storage of objects that handle output files.
    std::vector<TrajectoryFrameWriterPointer> writers_;
    //! Selection used for output.
    Selection sel_;
    //! Name for output file.
    std::string outputNamePrefix_;
    //! Name for index file.
    std::string indexFileName_;
    //! Storage of requirements for creating output files.
    OutputRequirementOptionDirector requirementsBuilder_;
    //! Stores the index information for the clusters. TODO refactor this!
    t_cluster_ndx* clusterIndex_ = nullptr;
};

ExtractCluster::ExtractCluster() {}

ExtractCluster::~ExtractCluster()
{
    if (clusterIndex_ != nullptr)
    {
        if (clusterIndex_->grpname != nullptr)
        {
            for (int i = 0; i < clusterIndex_->clust->nr; i++)
            {
                sfree(clusterIndex_->grpname[i]);
            }
            sfree(clusterIndex_->grpname);
        }
        if (clusterIndex_->clust != nullptr)
        {
            sfree(clusterIndex_->inv_clust);
            done_blocka(clusterIndex_->clust);
            sfree(clusterIndex_->clust);
        }
        sfree(clusterIndex_);
    }
}


void ExtractCluster::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] can be used to extract trajectory frames that correspond to clusters ",
        "obtained from running gmx cluster with the -clndx option.",
        "The module supports writing all GROMACS supported trajectory file formats.",
        "",
        "Included is also a selection of possible options to change additional information.",
        "",
        "It is possible to write only a selection of atoms to the output trajectory",
        "files for each cluster.",
    };

    options->addOption(FileNameOption("clusters")
                               .filetype(eftIndex)
                               .inputFile()
                               .required()
                               .store(&indexFileName_)
                               .defaultBasename("cluster")
                               .description("Name of index file containing frame indices for each "
                                            "cluster, obtained from gmx cluster -clndx."));

    options->addOption(SelectionOption("select").store(&sel_).onlyAtoms().description(
            "Selection of atoms to write to the file"));

    options->addOption(FileNameOption("o")
                               .filetype(eftTrajectory)
                               .outputFile()
                               .store(&outputNamePrefix_)
                               .defaultBasename("trajout")
                               .required()
                               .description("Prefix for the name of the trajectory file written "
                                            "for each cluster."));

    requirementsBuilder_.initOptions(options);

    settings->setHelpText(desc);
}

void ExtractCluster::optionsFinished(TrajectoryAnalysisSettings* settings)
{
    int frameFlags = TRX_NEED_X;

    frameFlags |= TRX_READ_V | TRX_READ_F;

    settings->setFrameFlags(frameFlags);

    clusterIndex_ = cluster_index(nullptr, indexFileName_.c_str());
}


void ExtractCluster::initAnalysis(const TrajectoryAnalysisSettings& /*settings*/,
                                  const TopologyInformation& top)
{
    int numberOfClusters = clusterIndex_->clust->nr;
    for (int i = 0; i < numberOfClusters; i++)
    {
        std::string outputName = Path::concatenateBeforeExtension(
                outputNamePrefix_, formatString("_%s", clusterIndex_->grpname[i]));
        writers_.emplace_back(createTrajectoryFrameWriter(top.mtop(), sel_, outputName,
                                                          top.hasTopology() ? top.copyAtoms() : nullptr,
                                                          requirementsBuilder_.process()));
    }
}

void ExtractCluster::analyzeFrame(int               frameNumber,
                                  const t_trxframe& frame,
                                  t_pbc* /* pbc */,
                                  TrajectoryAnalysisModuleData* /*pdata*/)
{
    // modify frame to write out correct number of coords
    // and actually write out
    int clusterToWriteTo = clusterIndex_->inv_clust[frameNumber];
    // Check for valid entry in cluster list, otherwise skip frame.
    if (clusterToWriteTo != -1 && clusterToWriteTo < clusterIndex_->clust->nr)
    {
        writers_[clusterToWriteTo]->prepareAndWriteFrame(frameNumber, frame);
    }
    else
    {
        printf("Frame %d was not found in any cluster!", frameNumber);
    }
}

void ExtractCluster::finishAnalysis(int /*nframes*/) {}


void ExtractCluster::writeOutput() {}

} // namespace

const char ExtractClusterInfo::name[] = "extract-cluster";
const char ExtractClusterInfo::shortDescription[] =
        "Allows extracting frames corresponding to clusters from trajectory";

TrajectoryAnalysisModulePointer ExtractClusterInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ExtractCluster);
}

} // namespace analysismodules

} // namespace gmx
