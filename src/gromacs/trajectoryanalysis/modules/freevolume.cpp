/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::Freevolume.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "freevolume.h"

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
class FreeVolume : public TrajectoryAnalysisModule
{
    public:
        FreeVolume();
        virtual ~FreeVolume() {};

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        std::string                       fnFreevol_;
        Selection                         sel_;
        AnalysisData                      data_;
        AnalysisDataAverageModulePointer  adata_;

        int                               nmol_;
        double                            mtot_;
        double                            cutoff_;
        double                            probeRadius_;
        gmx::DefaultRandomEngine          rng_;
        int                               seed_, ninsert_;
        AnalysisNeighborhood              nb_;
        //! The van der Waals radius per atom
        std::vector<double>               vdw_radius_;

        // Copy and assign disallowed by base.
};

// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.
FreeVolume::FreeVolume()
    : adata_(new AnalysisDataAverageModule())
{
    // We only compute two numbers per frame
    data_.setColumnCount(0, 2);
    // Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "freevolume");
    nmol_        = 0;
    mtot_        = 0;
    cutoff_      = 0;
    probeRadius_ = 0;
    seed_        = 0;
    ninsert_     = 1000;
}


void
FreeVolume::initOptions(IOptionsContainer          *options,
                        TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates the free volume in a box as",
        "a function of time. The free volume is",
        "plotted as a fraction of the total volume.",
        "The program tries to insert a probe with a given radius,",
        "into the simulations box and if the distance between the",
        "probe and any atom is less than the sums of the",
        "van der Waals radii of both atoms, the position is",
        "considered to be occupied, i.e. non-free. By using a",
        "probe radius of 0, the true free volume is computed.",
        "By using a larger radius, e.g. 0.14 nm, roughly corresponding",
        "to a water molecule, the free volume for a hypothetical",
        "particle with that size will be produced.",
        "Note however, that since atoms are treated as hard-spheres",
        "these number are very approximate, and typically only",
        "relative changes are meaningful, for instance by doing a",
        "series of simulations at different temperature.[PAR]",
        "The group specified by the selection is considered to",
        "delineate non-free volume.",
        "The number of insertions per unit of volume is important",
        "to get a converged result. About 1000/nm^3 yields an overall",
        "standard deviation that is determined by the fluctuations in",
        "the trajectory rather than by the fluctuations due to the",
        "random numbers.[PAR]",
        "The results are critically dependent on the van der Waals radii;",
        "we recommend to use the values due to Bondi (1964).[PAR]",
        "The Fractional Free Volume (FFV) that some authors like to use",
        "is given by 1 - 1.3*(1-Free Volume). This value is printed on",
        "the terminal."
    };

    settings->setHelpText(desc);

    // Add option for optional output file
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnFreevol_).defaultBasename("freevolume")
                           .description("Computed free volume"));

    // Add option for selecting a subset of atoms
    options->addOption(SelectionOption("select")
                           .store(&sel_).defaultSelectionText("all")
                           .onlyAtoms()
                           .description("Atoms that are considered as part of the excluded volume"));

    // Add option for the probe radius and initialize it
    options->addOption(DoubleOption("radius").store(&probeRadius_)
                           .description("Radius of the probe to be inserted (nm, 0 yields the true free volume)"));

    // Add option for the random number seed and initialize it to
    // generate a value automatically
    options->addOption(IntegerOption("seed").store(&seed_)
                           .description("Seed for random number generator (0 means generate)."));

    // Add option to determine number of insertion trials per frame
    options->addOption(IntegerOption("ninsert").store(&ninsert_)
                           .description("Number of probe insertions per cubic nm to try for each frame in the trajectory."));

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);
}


void
FreeVolume::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation        &top)
{
    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(adata_);

    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data et.
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(fnFreevol_);
    plotm_->setTitle("Free Volume");
    plotm_->setXAxisIsTime();
    plotm_->setYLabel("Free Volume (%)");
    plotm_->appendLegend("Free Volume");
    plotm_->appendLegend("Volume");

    data_.addModule(plotm_);

    // Initiate variable
    cutoff_               = 0;
    int            nnovdw = 0;
    gmx_atomprop_t aps    = gmx_atomprop_init();
    t_atoms       *atoms  = &(top.topology()->atoms);

    // Compute total mass
    mtot_ = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        mtot_ += atoms->atom[i].m;
    }

    // Extracts number of molecules
    nmol_ = top.topology()->mols.nr;

    // Loop over atoms in the selection using an iterator
    const int           maxnovdw = 10;
    ArrayRef<const int> atomind  = sel_.atomIndices();
    for (ArrayRef<const int>::iterator ai = atomind.begin(); (ai < atomind.end()); ++ai)
    {
        // Dereference the iterator to obtain an atom number
        int  i = *ai;
        real value;

        // Lookup the Van der Waals radius of this atom
        int resnr = atoms->atom[i].resind;
        if (TRUE == gmx_atomprop_query(aps, epropVDW,
                                       *(atoms->resinfo[resnr].name),
                                       *(atoms->atomname[i]),
                                       &value))
        {
            vdw_radius_.push_back(value);
            if (value > cutoff_)
            {
                cutoff_ = value;
            }
        }
        else
        {
            nnovdw++;
            if (nnovdw < maxnovdw)
            {
                fprintf(stderr, "Could not determine VDW radius for %s-%s. Set to zero.\n",
                        *(atoms->resinfo[resnr].name),
                        *(atoms->atomname[i]));
            }
            vdw_radius_.push_back(0.0);
        }
    }
    gmx_atomprop_destroy(aps);

    // Increase cutoff by proberadius to make sure we do not miss
    // anything
    cutoff_ += probeRadius_;

    if (nnovdw >= maxnovdw)
    {
        fprintf(stderr, "Could not determine VDW radius for %d particles. These were set to zero.\n", nnovdw);
    }

    if (seed_ == 0)
    {
        seed_ = static_cast<int>(gmx::makeRandomSeed());
    }

    // Print parameters to output. Maybe should make dependent on
    // verbosity flag?
    printf("cutoff       = %g nm\n", cutoff_);
    printf("probe_radius = %g nm\n", probeRadius_);
    printf("seed         = %d\n", seed_);
    printf("ninsert      = %d probes per nm^3\n", ninsert_);

    // Initiate the random number generator
    rng_.seed(seed_);

    // Initiate the neighborsearching code
    nb_.setCutoff(cutoff_);
}

void
FreeVolume::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                         TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle                   dh   = pdata->dataHandle(data_);
    const Selection                     &sel  = pdata->parallelSelection(sel_);
    gmx::UniformRealDistribution<real>   dist;

    GMX_RELEASE_ASSERT(nullptr != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Compute volume and number of insertions to perform
    real V       = det(fr.box);
    int  Ninsert = static_cast<int>(ninsert_*V);

    // Use neighborsearching tools!
    AnalysisNeighborhoodSearch nbsearch = nb_.initSearch(pbc, sel);

    // Then loop over insertions
    int NinsTot = 0;
    for (int i = 0; (i < Ninsert); i++)
    {
        rvec rand, ins, dx;

        for (int m = 0; (m < DIM); m++)
        {
            // Generate random number between 0 and 1
            // cppcheck-suppress uninitvar
            rand[m] = dist(rng_);
        }
        // Generate random 3D position within the box
        mvmul(fr.box, rand, ins);

        // Find the first reference position within the cutoff.
        bool                           bOverlap = false;
        AnalysisNeighborhoodPair       pair;
        AnalysisNeighborhoodPairSearch pairSearch = nbsearch.startPairSearch(ins);
        while (!bOverlap && pairSearch.findNextPair(&pair))
        {
            int jp = pair.refIndex();
            // Compute distance vector to first atom in the neighborlist
            pbc_dx(pbc, ins, sel.position(jp).x(), dx);

            // See whether the distance is smaller than allowed
            bOverlap = (norm(dx) <
                        probeRadius_+vdw_radius_[sel.position(jp).refId()]);

        }

        if (!bOverlap)
        {
            // We found some free volume!
            NinsTot++;
        }
    }
    // Compute total free volume for this frame
    double frac = 0;
    if (Ninsert > 0)
    {
        frac = (100.0*NinsTot)/Ninsert;
    }
    // Add the free volume fraction to the data set in column 0
    dh.setPoint(0, frac);
    // Add the total volume to the data set in column 1
    dh.setPoint(1, V);

    // Magic
    dh.finishFrame();
}


void
FreeVolume::finishAnalysis(int /* nframes */)
{
    please_cite(stdout, "Bondi1964a");
    please_cite(stdout, "Lourenco2013a");
}

void
FreeVolume::writeOutput()
{
    // Final results come from statistics module in analysis framework
    double FVaver  = adata_->average(0, 0);
    double FVerror = adata_->standardDeviation(0, 0);
    printf("Free volume %.2f +/- %.2f %%\n", FVaver, FVerror);

    double Vaver  = adata_->average(0, 1);
    double Verror = adata_->standardDeviation(0, 1);
    printf("Total volume %.2f +/- %.2f nm^3\n", Vaver, Verror);

    printf("Number of molecules %d total mass %.2f Dalton\n", nmol_, mtot_);
    double RhoAver  = mtot_ / (Vaver * 1e-24 * AVOGADRO);
    double RhoError = gmx::square(RhoAver / Vaver)*Verror;
    printf("Average molar mass: %.2f Dalton\n", mtot_/nmol_);

    double VmAver  = Vaver/nmol_;
    double VmError = Verror/nmol_;
    printf("Density rho: %.2f +/- %.2f nm^3\n", RhoAver, RhoError);
    printf("Molecular volume Vm assuming homogeneity: %.4f +/- %.4f nm^3\n",
           VmAver, VmError);

    double VvdWaver  = (1-FVaver/100)*VmAver;
    double VvdWerror = 0;
    printf("Molecular van der Waals volume assuming homogeneity:  %.4f +/- %.4f nm^3\n",
           VvdWaver, VvdWerror);

    double FFVaver  = 1-1.3*((100-FVaver)/100);
    double FFVerror = (FVerror/FVaver)*FFVaver;
    printf("Fractional free volume %.3f +/- %.3f\n", FFVaver, FFVerror);
}

}       // namespace

const char FreeVolumeInfo::name[]             = "freevolume";
const char FreeVolumeInfo::shortDescription[] =
    "Calculate free volume";

TrajectoryAnalysisModulePointer FreeVolumeInfo::create()
{
    return TrajectoryAnalysisModulePointer(new FreeVolume);
}

} // namespace analysismodules

} // namespace gmx
