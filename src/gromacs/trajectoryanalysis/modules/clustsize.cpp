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
 * Implements gmx::analysismodules::ClustSize.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "clustsize.h"

#include <algorithm>
#include <string>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/matio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
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
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

#include "unionfind.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/*! \brief
 * Class used to compute cluster sizes in the gas phase.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 *
 * \ingroup module_trajectoryanalysis
 */
class ClustSize : public TrajectoryAnalysisModule
{
    public:
        ClustSize();

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
        //! Filename for output file containing cluster size distribution
        std::string                       fnCSize_;
        //! Filename for output file containing weighed cluster size distribution
        std::string                       fnCSizeW_;
        //! Filename for output file containing atom numbers in largest cluster
        std::string                       fnMaxClustNdx_;
        //! Number of "particles"
        int                               nparticles_;
        //! Selection to use either atom or mol
        Selection                         sel_;
        //! Data storage for maxclust.xvg
        AnalysisData                      maxClust_;
        //! Data storage for nclust.xvg
        AnalysisData                      nClust_;
        //! Data storage for avclust.xvg
        AnalysisData                      avClust_;
        //! Data storage for histoclust.xvg
        AnalysisData                      clustSizeDist_;
        /*! \brief
         * Histogram module that computes the actual histogram
         * from `clustSizeDist_`.
         *
         * The per-frame histograms are cluster counts for each
         * cluster size.
         */
        AnalysisDataSimpleHistogramModulePointer pairCounts_;
        // Plot settings
        AnalysisDataPlotSettings                 plotSettings_;
        //! Count molecules rather than atoms
        bool                                     mol_;
        //! Molecule index for each atom
        std::vector<int>                         molIndex_;
        //! Cut-off distance
        double                                   cutoff_;
        //! For xpm files the starting color (zero == white)
        RVec                                     rmid_;
        //! For xpm files the ending color
        RVec                                     rhi_;
        //! Number of levels in the xpm color map for the -o option
        int                                      nlevels_;
        //! Number of levels in the xpm color map for the -ow option
        int                                      nlevelsW_;
        //! Neighborsearching data structure
        AnalysisNeighborhood                     nb_;

        // Copy and assign disallowed by base.
};

// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.
ClustSize::ClustSize()
    : pairCounts_(new AnalysisDataSimpleHistogramModule())
{
    cutoff_       = 0.35;
    mol_          = false;
    nlevels_      = 20;
    nlevelsW_     = 20;
    rmid_         = { 1, 1, 0 };
    rhi_          = { 0, 0, 1 };
    nparticles_   = 0;
    // Tell the analysis framework that these components exist
    registerAnalysisDataset(&maxClust_, "maxclust");
    registerAnalysisDataset(&nClust_, "nclust");
    registerAnalysisDataset(&avClust_, "avclust");
    registerAnalysisDataset(&clustSizeDist_, "clustSizeDist");
    registerBasicDataset(pairCounts_.get(), "pairCounts");
}

void
ClustSize::initOptions(IOptionsContainer          *options,
                       TrajectoryAnalysisSettings *settings)
{
    const char *const desc[] = {
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
                           .description("Max cluster size"));
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
    options->addOption(FileNameOption("o").filetype(eftXpm).outputFile()
                           .store(&fnCSize_).defaultBasename("csize")
                           .description("Cluster size distribution"));
    options->addOption(FileNameOption("ow").filetype(eftXpm).outputFile()
                           .store(&fnCSizeW_).defaultBasename("csizew")
                           .description("Weighted cluster size distribution"));
    options->addOption(FileNameOption("mcn").filetype(eftIndex).outputFile()
                           .store(&fnMaxClustNdx_).defaultBasename("maxclust")
                           .description("Atoms in largest cluster"));
    options->addOption(SelectionOption("sel").store(&sel_)
                           .required()
                           .description("Cluster particle selection"));
    // Add option for selecting a subset of atoms
    options->addOption(BooleanOption("mol")
                           .store(&mol_)
                           .description("Cluster molecules instead of atoms"));

    // Add option for the probe radius and initialize it
    options->addOption(DoubleOption("cut").store(&cutoff_)
                           .description("Largest distance (nm) to be considered in a cluster."));
    options->addOption(IntegerOption("nlevels").store(&nlevels_)
                           .description("Number of levels of grey in csize.xpm output (-o flag) "));
    options->addOption(IntegerOption("nlevelsW").store(&nlevelsW_)
                           .description("Number of levels of grey in csizew.xpm output (-ow flag)"));
    //    options->addOption(RvecOption("rgblo").store(rmid_)
    //                  .description("RGB values for the color of the lowest occupied cluster size"));
    //options->addOption(RvecOption("rgbhi").store(rhi_)
    //                 .description("RGB values for the color of the highest occupied cluster size"));

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop);
    settings->setPBC(true);
}

/*! \brief Convenience function for adding a plotfile
 *
 * \param[out] data     The data handle to add stuff to
 * \param[in]  fn       The filename
 * \param[in]  title    Title of the xvg file
 * \param[in]  settings Used forplot settings
 * \param[in]  ylabel   Label for Y axis
 */
void addPlotModule(AnalysisData                     *data,
                   const std::string                &fn,
                   const char                       *title,
                   const TrajectoryAnalysisSettings &settings,
                   const char                       *ylabel)
{
    if (!fn.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule());
        plotm->setSettings(settings.plotSettings());
        plotm->setFileName(fn);
        plotm->setTitle(title);
        plotm->setXAxisIsTime();
        plotm->setYLabel(ylabel);
        data->setColumnCount(0, 1);
        data->addModule(plotm);
    }
}

void
ClustSize::initAnalysis(const TrajectoryAnalysisSettings &settings,
                        const TopologyInformation        &top)
{
    // Initiate selection and compute number of particles
    nparticles_   = sel_.initOriginalIdsToGroup(top.mtop(),
                                                mol_ ? INDEX_MOL : INDEX_ATOM);
    // Store plot settings
    plotSettings_ = settings.plotSettings();

    // Set the total number of data points in the histogram to the
    // number of atoms in the selection, alternative number of molecules
    clustSizeDist_.setMultipoint(false);
    clustSizeDist_.setDataSetCount(1);
    clustSizeDist_.setColumnCount(0, nparticles_);
    clustSizeDist_.requestStorage(-1);
    clustSizeDist_.addModule(pairCounts_);
    // Add modules for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data set.
    addPlotModule(&maxClust_, fnMaxClust_, "Max cluster size",
                  settings, "#molecules");
    addPlotModule(&nClust_, fnNClust_, "Number of clusters", settings, "N");
    addPlotModule(&avClust_, fnAvClust_, "Average cluster size",
                  settings, "#molecules");

    // Initiate the neighborsearching code
    nb_.setCutoff(cutoff_);

    // Initiate the histoqgram
    pairCounts_->init(histogramFromRange(0.0, nparticles_).binWidth(1.0).integerBins(true));

    // Initiate molecule index
    t_topology *tp = top.topology();
    for (int i = 0; i < tp->mols.nr; i++)
    {
        for (int j = tp->mols.index[i]; j < tp->mols.index[i+1]; j++)
        {
            molIndex_.push_back(i);
        }
    }
}

void
ClustSize::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                        TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle dhmc, dhnc, dhac, dhhc;

    dhmc = pdata->dataHandle(maxClust_);
    if (dhmc.isValid())
    {
        dhmc.startFrame(frnr, fr.time);
    }
    dhnc = pdata->dataHandle(nClust_);
    if (dhnc.isValid())
    {
        dhnc.startFrame(frnr, fr.time);
    }
    dhac = pdata->dataHandle(avClust_);
    if (dhac.isValid())
    {
        dhac.startFrame(frnr, fr.time);
    }
    dhhc = pdata->dataHandle(clustSizeDist_);
    if (dhhc.isValid())
    {
        dhhc.startFrame(frnr, fr.time);
    }
    const Selection               &sel        = pdata->parallelSelection(sel_);
    AnalysisNeighborhoodSearch     nbsearch   = nb_.initSearch(pbc, sel);
    AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startSelfPairSearch();

    // Class to manage sets of indices
    MappedUnionFinder              finder;
    finder.initWithGroupIndices(sel.mappedIds());

    // Initiate the cluster index to each particle having its own cluster
    AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        int mi = pair.refIndex();
        int mj = pair.testIndex();

        if (mol_)
        {
            mi = molIndex_[mi];
            mj = molIndex_[mj];
        }
        /* Merge clusters */
        finder.mergeGroups(mi, mj);
    }

    std::vector<int> cluster_size = finder.allSizes();
    if (dhhc.isValid())
    {
        int k = 0;
        for (auto &cs : cluster_size)
        {
            dhhc.setPoint(k, cs);
            k++;
        }
    }

    if (dhmc.isValid())
    {
        if (cluster_size.size() > 0)
        {
            dhmc.setPoint(0, *std::max_element(cluster_size.begin(),
                                               cluster_size.end()));
        }
        else
        {
            dhmc.setPoint(0, 0);
        }
        dhmc.finishFrame();
    }
    if (dhnc.isValid())
    {
        dhnc.setPoint(0, cluster_size.size());
        dhnc.finishFrame();
    }
    if (dhhc.isValid())
    {
        dhhc.finishFrame();
    }
    if (dhac.isValid())
    {
        int inCluster = 0;
        int nCluster  = 0;
        for (auto &cs : cluster_size)
        {
            if (cs > 1)
            {
                inCluster += cs;
                nCluster++;
            }
        }
        if (nCluster > 0)
        {
            real avclust = 1.0*inCluster/nCluster;
            dhac.setPoint(0, avclust);
        }
        dhac.finishFrame();
    }
}

void
ClustSize::finishAnalysis(int /* nframes */)
{
    if (!fnCSize_.empty() || !fnCSizeW_.empty())
    {
        int                             max_clust_size = 0;
        int                             max_num_clust  = 0;
        std::vector<std::vector<real> > csreal;
        std::vector<real>               xaxis;

        for (int i = 0; i < clustSizeDist_.frameCount(); i++)
        {
            auto              frame = clustSizeDist_.getDataFrame(i);
            std::vector<real> cs;
            int               k;
            xaxis.push_back(frame.x());
            for (k = 0; k < clustSizeDist_.columnCount(); k++)
            {
                auto value = frame.y(k);
                if (value == 0)
                {
                    break;
                }
                int ival = static_cast<int>(std::round(value));
                max_clust_size = std::max(max_clust_size, ival);
                if (cs.size() <= value)
                {
                    cs.resize(ival, 0);
                }
                int index = ival - 1;
                cs[index] += 1;
                if (cs[index] > max_num_clust)
                {
                    max_num_clust = std::round(cs[index]);
                }
            }
            csreal.push_back(cs);
        }
        if (max_clust_size == 0)
        {
            printf("All data is zero, not generating %s\n",
                   fnCSize_.c_str());
            return;
        }

        std::vector<real> yaxis;
        yaxis.resize(max_clust_size, 0);
        for (size_t i = 0; i < yaxis.size(); i++)
        {
            yaxis[i] = i+1;
        }
        std::vector<real *> cs_xpm;
        for (size_t k = 0; k < csreal.size(); k++)
        {
            csreal[k].resize(max_clust_size+1, 0);
            cs_xpm.push_back(csreal[k].data());
        }

        t_rgb rlo, rhi, rmid;
        // TODO: implement user input for colours
        rlo.r  = 1;
        rlo.g  = 1;
        rlo.b  = 1;
        rmid.r = rmid_[0];
        rmid.g = rmid_[1];
        rmid.b = rmid_[2];
        rhi.r  = rhi_[0];
        rhi.g  = rhi_[1];
        rhi.b  = rhi_[2];
        // TODO: insert correct time unit
        std::string timebuf = formatString("Time (ps)");

        if (!fnCSize_.empty())
        {
            real cmid = 1;
            real cmax = 0.0;
            for (size_t i = 0; (i < cs_xpm.size()); i++)
            {
                for (size_t j = 0; (j < yaxis.size()); j++)
                {
                    if ((cs_xpm[i][j] > 0) && (cs_xpm[i][j] < cmid))
                    {
                        cmid = cs_xpm[i][j];
                    }
                    cmax = std::max(cs_xpm[i][j], cmax);
                }
            }

            FILE *fp  = gmx_ffopen(fnCSize_.c_str(), "w");
            write_xpm3(fp, 0, "Cluster size distribution", "# clusters",
                       timebuf, "Size",
                       xaxis.size(), yaxis.size(),
                       xaxis.data(), yaxis.data(),
                       cs_xpm.data(), 0, cmid, max_num_clust,
                       rlo, rmid, rhi, &nlevels_);
            gmx_ffclose(fp);
        }

        if (!fnCSizeW_.empty())
        {
            real cmid = 100.0;
            real cmax = 0.0;
            for (size_t i = 0; (i < cs_xpm.size()); i++)
            {
                for (size_t j = 0; (j < yaxis.size()); j++)
                {
                    cs_xpm[i][j] *= (j+1);
                    if ((cs_xpm[i][j] > 0) && (cs_xpm[i][j] < cmid))
                    {
                        cmid = cs_xpm[i][j];
                    }
                    cmax = std::max(cs_xpm[i][j], cmax);
                }
            }

            FILE *fp = gmx_ffopen(fnCSizeW_.c_str(), "w");
            write_xpm3(fp, 0, "Weighted cluster size distribution", "Fraction",
                       timebuf, "Size",
                       xaxis.size(), yaxis.size(),
                       xaxis.data(), yaxis.data(),
                       cs_xpm.data(), 0, cmid, cmax,
                       rlo, rmid, rhi, &nlevelsW_);
            gmx_ffclose(fp);
        }
    }

    if (!fnHistoClust_.empty())
    {
        AbstractAverageHistogram &finalHisto
            = pairCounts_->averager();
        finalHisto.done();

        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(plotSettings_));
        plotm->setFileName(fnHistoClust_);
        plotm->setTitle("Cluster size distribution");
        plotm->setXLabel("Cluster size");
        plotm->setYLabel("()");
        finalHisto.addModule(plotm);
    }
}

void
ClustSize::writeOutput()
{
    // Final results come from statistics module in analysis framework
}

}       // namespace

const char ClustSizeInfo::name[]             = "clustsize";
const char ClustSizeInfo::shortDescription[] =
    "Calculate cluster size distributions";

TrajectoryAnalysisModulePointer ClustSizeInfo::create()
{
    return TrajectoryAnalysisModulePointer(new ClustSize);
}

} // namespace analysismodules

} // namespace gmx
