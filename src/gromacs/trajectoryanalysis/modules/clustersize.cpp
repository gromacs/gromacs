/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::ClusterSize.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "clustersize.h"

#include <string>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/threefry.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*! \brief
 * Class used to compute free volume in a simulations box.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 * Does not implement any new functionality.
 *
 * \ingroup module_trajectoryanalysis
 */
class ClusterSize : public TrajectoryAnalysisModule
{
    public:
        ClusterSize();
        virtual ~ClusterSize() {};

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        /*! \brief Convenience function for adding a plotfile
         *
         * \param[in] fn       The filename
         * \param[in] title    Title of the xvg file
         * \param[in] settings Used forplot settings
         */
        void addPlotModule(const std::string                &fn,
                           const char                       *title,
                           const TrajectoryAnalysisSettings &settings);
        std::string                       fnMaxClust_;
        std::string                       fnNClust_;
        std::string                       fnAvClust_;
        std::string                       fnHistoClust_;
        std::string                       fnTemperature_;
        std::string                       fnMaxClustNdx_;
        // Selection to use either atom or mol
        Selection                         sel_;
        AnalysisData                      data_;
        AnalysisDataAverageModulePointer  adata_;
        bool                              mol_;
        double                            cutoff_;
        RVec                              rlo_;
        RVec                              rhi_;
        gmx::DefaultRandomEngine          rng_;
        int                               nlevels_, nskip_;
        AnalysisNeighborhood              nb_;

        // Copy and assign disallowed by base.
};

// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.
ClusterSize::ClusterSize()
    : adata_(new AnalysisDataAverageModule()), rlo_({1, 0, 0}
                                                    ), rhi_({0, 0, 1})
{
    // We compute many numbers per frame
    data_.setColumnCount(0, 1);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "clustersize");
    cutoff_       = 0;
    nlevels_      = 10;
    nskip_        = 0;
    mol_          = false;
}

void
ClusterSize::initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes the size distributions of molecular/atomic clusters in",
        "the gas phase. The output is given in the form of an [REF].xpm[ref] file.",
        "The total number of clusters is written to an [REF].xvg[ref] file.[PAR]",
        "When the [TT]-mol[tt] option is given clusters will be made out of",
        "molecules rather than atoms, which allows clustering of large molecules.",
        "In this case an index file would still contain atom numbers",
        "or your calculation will die with a SEGV.[PAR]",
        "When velocities are present in your trajectory, the temperature of",
        "the largest cluster will be printed in a separate [REF].xvg[ref] file assuming",
        "that the particles are free to move. If you are using constraints,",
        "please correct the temperature. For instance water simulated with SHAKE",
        "or SETTLE will yield a temperature that is 1.5 times too low. You can",
        "compensate for this with the [TT]-ndf[tt] option. Remember to take the removal",
        "of center of mass motion into account.[PAR]",
        "The [TT]-mc[tt] option will produce an index file containing the",
        "atom numbers of the largest cluster."
    };

    settings->setHelpText(desc);

    // Add option for optional output file
    options->addOption(FileNameOption("mc").filetype(eftPlot).outputFile()
                           .store(&fnMaxClust_).defaultBasename("maxclust")
                           .description("Size of largest cluster"));
    options->addOption(FileNameOption("nc").filetype(eftPlot).outputFile()
                           .store(&fnNClust_).defaultBasename("nclust")
                           .description("Number of clusters"));
    options->addOption(FileNameOption("ac").filetype(eftPlot).outputFile()
                           .store(&fnAvClust_).defaultBasename("avclust")
                           .description("Average cluster size"));
    options->addOption(FileNameOption("hc").filetype(eftPlot).outputFile()
                           .store(&fnHistoClust_).defaultBasename("histo-clust")
                           .description("Cluster size histogram"));
    options->addOption(FileNameOption("temp").filetype(eftPlot).outputFile()
                           .store(&fnTemperature_).defaultBasename("temp")
                           .description("Temperature of largest cluster"));
    options->addOption(FileNameOption("mcn").filetype(eftIndex).outputFile()
                           .store(&fnMaxClustNdx_).defaultBasename("maxclust")
                           .description("Atoms in largest cluster"));
    options->addOption(SelectionOption("selection").store(&sel_)
                           .required()
                           .description("Cluster particle selection"));
    // Add option for selecting a subset of atoms
    options->addOption(BooleanOption("mol")
                           .store(&mol_)
                           .description("Cluster molecules instead of atoms"));

    // Add option for the probe radius and initialize it
    options->addOption(DoubleOption("cut").store(&cutoff_)
                           .description("Largest distance (nm) to be considered in a cluster."));
    //options->addOption(BoolOption("pbc").store(&pbc_)
    //                       .description("Use periodic boundary conditions."));
    options->addOption(IntegerOption("nskip").store(&nskip_)
                           .description("Number of frames to skip between writing."));
    options->addOption(IntegerOption("nlevels").store(&nlevels_)
                           .description("Number of levels of grey in [REF].xpm[ref] output."));
    //    options->addOption(RvecOption("rgblo").store(rlo_)
    //                  .description("RGB values for the color of the lowest occupied cluster size"));
    //options->addOption(RvecOption("rgbhi").store(rhi_)
    //                 .description("RGB values for the color of the highest occupied cluster size"));

    // Control input settings
    //    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
    //                 TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);
}

void ClusterSize::addPlotModule(const std::string                &fn,
                                const char                       *title,
                                const TrajectoryAnalysisSettings &settings)
{
    if (!fn.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule());
        plotm->setSettings(settings.plotSettings());
        plotm->setFileName(fn);
        plotm->setTitle(title);
        plotm->setXAxisIsTime();
        if (mol_)
        {
            plotm->setYLabel("Molecules");
        }
        else
        {
            plotm->setYLabel("Atoms");
        }
    }
}

int initSelectionGroups(Selection *sel, const gmx_mtop_t *top, bool mol)
{
    if (mol)
    {
        return sel->initOriginalIdsToGroup(top, INDEX_MOL);
    }
    else
    {
        return sel->initOriginalIdsToGroup(top, INDEX_ATOM);
    }
}

void
ClusterSize::initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top)
{
    (void) initSelectionGroups(&sel_, top.mtop(), mol_);

    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);

    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data et.
    addPlotModule(fnMaxClust_, "Size of Largest Cluster", settings);
    addPlotModule(fnNClust_, "Number of Clusters", settings);
    addPlotModule(fnAvClust_, "Average Cluster Size", settings);

    // Increase cutoff by proberadius to make sure we do not miss
    // anything

    // verbosity flag?
    printf("cutoff       = %g nm\n", cutoff_);

    if (mol_)
    {
        // /

    }
    // Initiate the neighborsearching code
    nb_.setCutoff(cutoff_);
}

void
ClusterSize::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh   = pdata->dataHandle(data_);
    const Selection    &sel  = pdata->parallelSelection(sel_);

    GMX_RELEASE_ASSERT(nullptr != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Use neighborsearching tools!
    AnalysisNeighborhoodSearch     nbsearch = nb_.initSearch(pbc, sel);
    AnalysisNeighborhoodPair       pair;
    AnalysisNeighborhoodPositions  positions(fr.x, fr.natoms);
    AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(positions);
    std::vector<int>               cluster_index, cluster_size;
    cluster_index.resize(sel.atomCount(), -1);
    cluster_size.resize(sel.atomCount(), 0);
    int max_clust_index = 0;
    while (pairSearch.findNextPair(&pair))
    {
        int pi = pair.refIndex();
        int pj = pair.testIndex();
        if ((cluster_index[pi] == -1) && (cluster_index[pj] == -1))
        {
            cluster_index[pi] = cluster_index[pj] = max_clust_index++;
        }
        else if (cluster_index[pi] == -1)
        {
            cluster_index[pi] = cluster_index[pj];
        }
        else if (cluster_index[pj] == -1)
        {
            cluster_index[pj] = cluster_index[pi];
        }
        else
        {
            /* Merge clusters: check for all atoms whether they are in
             * cluster cluster_index[pj] and if so, put them in cluster_index[pi]
             */
            for (int k = 0; (k < sel.atomCount()); k++)
            {
                if (cluster_index[k] == cluster_index[pj])
                {
                    if (cluster_size[cluster_index[pj]] <= 0)
                    {
                        gmx_fatal(FARGS, "negative cluster size %d for element %d",
                                  cluster_size[cluster_index[pj]], cluster_index[pj]);
                    }
                    cluster_size[cluster_index[pj]]--;
                    cluster_index[k] = cluster_index[pi];
                    cluster_size[cluster_index[pi]]++;
                }
            }
        }
    }
    for (size_t k = 0; k < cluster_index.size(); k++)
    {
        printf("cluster %d size %d\n", static_cast<int>(k), cluster_size[k]);
    }
    // Magic
    dh.finishFrame();
}


void
ClusterSize::finishAnalysis(int /* nframes */)
{
    please_cite(stdout, "Bondi1964a");
    please_cite(stdout, "Lourenco2013a");
}

void
ClusterSize::writeOutput()
{
    // Final results come from statistics module in analysis framework
}

}       // namespace

const char ClusterSizeInfo::name[]             = "clustersize";
const char ClusterSizeInfo::shortDescription[] =
    "Calculate cluster size distributions";

TrajectoryAnalysisModulePointer ClusterSizeInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ClusterSize);
}

} // namespace analysismodules

} // namespace gmx
