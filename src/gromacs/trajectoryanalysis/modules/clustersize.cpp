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
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

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
        //! Filename for output file containing size of largest cluster
        std::string                       fnMaxClust_;
        //! Filename for output file containing number of clusters
        std::string                       fnNClust_;
        //! Filename for output file containing average size of clusters
        std::string                       fnAvClust_;
        //! Filename for output file containing histogram of cluster sizes
        std::string                       fnHistoClust_;
        //! Filename for output file containing temperature of largest cluster
        std::string                       fnTemperature_;
        //! Filename for output file containing atom numbers in largest cluster
        std::string                       fnMaxClustNdx_;
        //! Atom to molecule number translation
        std::vector<int>                  molIndex_;
        //! Selection to use either atom or mol
        Selection                         sel_;
        AnalysisData                      maxclust_;
        AnalysisData                      nclust_;
        AnalysisData                      avclust_;
        AnalysisDataAverageModulePointer  adata_;
        //! Count molecules rather than atoms
        bool                              mol_;
        //! Number of molecules
        size_t                            nmol_;
        double                            cutoff_;
        RVec                              rlo_;
        RVec                              rhi_;
        int                               nlevels_;
        int                               nskip_;
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
    cutoff_       = 0.35;
    nlevels_      = 10;
    nskip_        = 0;
    mol_          = false;
    nmol_         = 0;
    // Tell the analysis framework that these components exist
    registerAnalysisDataset(&maxclust_, "maxclust");
    registerAnalysisDataset(&nclust_, "nclust");
    registerAnalysisDataset(&avclust_, "avclust");
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

    // Add options for optional output files
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

/*! \brief Convenience function for adding a plotfile
 *
 * \param[out] data     The data handle to add stuff to
 * \param[in]  fn       The filename
 * \param[in]  title    Title of the xvg file
 * \param[in]  settings Used forplot settings
 * \param[in]  mol      Put molecules or atoms on the Y axis
 */
void addPlotModule(AnalysisData                     *data,
                   const std::string                &fn,
                   const char                       *title,
                   const TrajectoryAnalysisSettings &settings,
                   bool                              mol)
{
    if (!fn.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule());
        plotm->setSettings(settings.plotSettings());
        plotm->setFileName(fn);
        plotm->setTitle(title);
        plotm->setXAxisIsTime();
        if (mol)
        {
            plotm->setYLabel("Molecules");
        }
        else
        {
            plotm->setYLabel("Atoms");
        }
        data->setColumnCount(0, 1);
        data->addModule(plotm);
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

    // Add modules for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data set.
    addPlotModule(&maxclust_, fnMaxClust_, "Size of Largest Cluster", settings, mol_);
    addPlotModule(&nclust_, fnNClust_, "Number of Clusters", settings, mol_);
    addPlotModule(&avclust_, fnAvClust_, "Average Cluster Size", settings, mol_);

    // verbosity flag?
    printf("cutoff       = %g nm\n", cutoff_);

    if (mol_)
    {
        // Fill molIndex_ vector with molecules number for each atom
        const gmx_mtop_t *mtop = top.mtop();
        nmol_                  = 0;
        for (int i = 0; i < mtop->nmolblock; i++)
        {
            for (int j = 0; j < mtop->molblock[i].nmol; j++)
            {
                for (int k = 0; k < mtop->molblock[i].natoms_mol; k++)
                {
                    molIndex_.push_back(nmol_);
                }
                nmol_++;
            }

        }
        printf("There are %d molecules in your system\n", static_cast<int>(nmol_));
    }
    // Initiate the neighborsearching code
    nb_.setCutoff(cutoff_);
}

void
ClusterSize::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata)
{
    if (nskip_ > 0 && frnr % nskip_ != 0)
    {
        return;
    }
    GMX_RELEASE_ASSERT(nullptr != pbc,
                       "You have no periodic boundary conditions");

    AnalysisDataHandle dhmc, dhnc, dhac;
    if (!fnMaxClust_.empty())
    {
        dhmc = pdata->dataHandle(maxclust_);
        dhmc.startFrame(frnr, fr.time);
    }
    if (!fnNClust_.empty())
    {
        dhnc = pdata->dataHandle(nclust_);
        dhnc.startFrame(frnr, fr.time);
    }
    if (!fnAvClust_.empty())
    {
        dhac = pdata->dataHandle(avclust_);
        dhac.startFrame(frnr, fr.time);
    }
    const Selection    &sel  = pdata->parallelSelection(sel_);

    // Use neighborsearching tools!
    AnalysisNeighborhoodSearch     nbsearch = nb_.initSearch(pbc, sel);
    AnalysisNeighborhoodPair       pair;
    AnalysisNeighborhoodPositions  positions(fr.x, fr.natoms);
    AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(positions);
    std::vector<int>               cluster_index, cluster_size;
    cluster_index.resize(sel.atomCount(), -1);
    cluster_size.resize(sel.atomCount(), 0);
    int max_clust_index = 0;
    if (mol_)
    {
        // Put each whole molecule in a cluster of it's own
        int k = 0;
        for (auto &mi : molIndex_)
        {
            cluster_index[k] = mi;
            cluster_size[mi]++;
            k++;
        }
    }
    while (pairSearch.findNextPair(&pair))
    {
        int pi = pair.refIndex();
        int pj = pair.testIndex();
        if (pj <= pi)
        {
            continue;
        }
        printf("Found pair %d %d cluster_index before %d %d", pi, pj, cluster_index[pi], cluster_index[pj]);
        if ((cluster_index[pi] == -1) && (cluster_index[pj] == -1))
        {
            cluster_index[pi]                = cluster_index[pj] = max_clust_index++;
            cluster_size[cluster_index[pi]] += 2;
        }
        else if (cluster_index[pi] == -1)
        {
            cluster_index[pi]                = cluster_index[pj];
            cluster_size[cluster_index[pi]] += 1;
        }
        else if (cluster_index[pj] == -1)
        {
            cluster_index[pj]                = cluster_index[pi];
            cluster_size[cluster_index[pi]] += 1;
        }
        else if (cluster_index[pi] != cluster_index[pj])
        {
            /* Merge clusters: check for all atoms whether they are in
             * cluster cluster_index[pj] and if so,
             * put them in cluster_index[pi]
             */
            for (int k = 0; (k < sel.atomCount()); k++)
            {
                if (cluster_index[k] == cluster_index[pj])
                {
                    GMX_RELEASE_ASSERT(cluster_size[cluster_index[pj]] > 0,
                                       gmx::formatString("negative cluster size %d for element %d", cluster_size[cluster_index[pj]], cluster_index[pj]).c_str());

                    cluster_size[cluster_index[pj]]--;
                    cluster_index[k] = cluster_index[pi];
                    cluster_size[cluster_index[pi]]++;
                }
            }
        }
        printf(" after %d %d\n", cluster_index[pi], cluster_index[pj]);
    }
    int nclust     = 0;
    int maxclust   = 0;
    int maxclustID = -1;
    int k          = 0;
    for (auto &cs : cluster_size)
    {
        if (cs > 0)
        {
            nclust++;
            if (cs > maxclust)
            {
                maxclust   = cs;
                maxclustID = k;
            }
        }
        k++;
    }
    if (mol_)
    {
        std::vector<bool> inMaxClust;
        inMaxClust.resize(nmol_, false);
        int               k = 0;
        for (auto &ci : cluster_index)
        {
            if (ci == maxclustID)
            {
                inMaxClust[molIndex_[k]] = true;
            }
            k++;
        }
        maxclust = std::count_if(inMaxClust.begin(),
                                 inMaxClust.end(),
                                 [](bool b) { return b; });
    }
    if (!fnMaxClust_.empty())
    {
        dhmc.setPoint(0, maxclust);
        dhmc.finishFrame();
    }
    if (!fnNClust_.empty())
    {
        dhnc.setPoint(0, nclust);
        dhnc.finishFrame();
    }
    if (!fnAvClust_.empty())
    {
        real avclust = sel.atomCount()/(1.0*nclust);
        dhac.setPoint(0, avclust);
        dhac.finishFrame();
    }
}


void
ClusterSize::finishAnalysis(int /* nframes */)
{
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
