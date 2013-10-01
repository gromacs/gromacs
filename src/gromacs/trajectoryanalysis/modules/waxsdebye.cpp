/*
 * This file is part of the GROMACS molecular simulation package.
 *
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
 * Implements gmx::analysismodules::Waxsdebye.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "waxsdebye.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/topology.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/vec.h"
#include "gromacs/legacyheaders/vsite.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/random/random.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/waxsdebye/waxs_debye_force.h"

namespace gmx
{

namespace analysismodules
{

/*! \brief
 * Class used to compute waxs scatering in a simulations box.
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 * Does not implement any new functionality.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class WaxsDebye : public TrajectoryAnalysisModule
{
    public:
        //! Name of the tool
        static const char name[];

        //! One line description
        static const char shortDescription[];

        //! Constructor
        WaxsDebye();

        //! Destructor
        virtual ~WaxsDebye();

        //! Set the options and setting
        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);

        //! First routine called by the analysis frame work
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        //! Call for each frame of the trajectory
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        //! Last routine called by the analysis frame work
        virtual void finishAnalysis(int nframes);

        //! Routine to write output, that is additional over the built-in
        virtual void writeOutput();

    private:
        //! File names for input and output
        std::string                       fnSfactor_, fnSqref_, fnSqdiff_, fnSqcalc_, fnEner_, fnAlpha_;
        //! Selection - is this necessary?
        //Selection                         sel_;

        //! Number of interactions
        int                               nbonds_;

        //! Interaction definition for waxs debye
        t_iatom                          *iatoms_wd_;

        //! Topology data structure
        t_topology                       *top_;

        //! Virtual site creation structure
        gmx_vsite_t                      *vsite_;

        //! Force vector
        rvec                             *f_;

        //! Analysis data structure
        AnalysisData                      data_;

        //! Analysis data structures
        //AnalysisDataAverageModulePointer  alphadata_;
        AnalysisDataAverageModulePointer  enerdata_;

        //! WAXS data structure that does the real work
        WaxsDebyeForce                   *wdf_;

        //! Communication Record
        t_commrec                        *cr_;
        // Copy and assign disallowed by base.
};


/*! \brief Constructor.
 * Here it is important to initialize the pointer to
 * subclasses that are elements of the main class. Here we have only
 * one. The type of this depends on what kind of tool you need.
 * Here we only have simple value/time kind of data.
 */
WaxsDebye::WaxsDebye()
    : TrajectoryAnalysisModule(WaxsDebyeInfo::name, WaxsDebyeInfo::shortDescription),
      enerdata_(new AnalysisDataAverageModule())
{
    //! We only compute one number per frame
    data_.setColumnCount(0, 2);

    //! Tell the analysis framework that this component exists
    registerAnalysisDataset(&data_, "waxsdebye");

    //! Initialize private variables to NULL/0
    wdf_        = NULL;
    f_          = NULL;
    nbonds_     = 0;
    iatoms_wd_  = NULL;
    cr_         = init_commrec();
    top_        = NULL;
    vsite_      = NULL;
}

WaxsDebye::~WaxsDebye()
{
    // C++ takes care of memory in classes (hopefully)
    if (NULL != wdf_)
    {
        delete wdf_;
    }
    if (NULL != iatoms_wd_)
    {
        sfree(iatoms_wd_);
    }
    if (NULL != top_)
    {
        //done_top(top_);
    }
}


void
WaxsDebye::initOptions(Options                    *options,
                       TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "gmx waxsdebye can compute WAXS/SAXS scattering curves",
        "for each frame in a trajectory, and plot those as",
        "line graphs in an xvg file. It can also compute the deviation",
        "from a reference curve if a reference curve is given.",
        "The code is based on the same data structures that are used",
        "for WAXS/SAXS refinement in mdrun and hence the results",
        "are compatible with the refinement. The value of alpha is computed",
        "using a bisection algorithm for alpha that minimizes the energy,",
        "in effect minimizing the energy with respect to alpha, within",
        "the bounds given in the [TT]tpr[tt] file. Finally the tool can",
        "do some data management, like scaling a reference (e.g. experimental)",
        "data file to match the expected S(0) based on the Debye formula."
    };

    // Add the descriptive text (program help text) to the options
    options->setDescription(concatenateStrings(desc));

    // Add option for optional input file
    options->addOption(FileNameOption("sfac").filetype(eftXML)
                           .inputFile().required()
                           .store(&fnSfactor_).defaultBasename("sfactor")
                           .description("Structure factors"));

    // Add option for optional input file
    options->addOption(FileNameOption("waxs_ref").filetype(eftGenericData)
                           .inputFile().required()
                           .store(&fnSqref_).defaultBasename("sq")
                           .description("Reference S(q)"));

    options->addOption(FileNameOption("waxs_diff").filetype(eftGenericData).inputFile()
                           .store(&fnSqdiff_).defaultBasename("sq")
                           .description("Difference S(q)"));

    // Add option for optional output file
    options->addOption(FileNameOption("waxs_out").filetype(eftPlot)
                           .outputFile().required()
                           .store(&fnSqcalc_).defaultBasename("sqt")
                           .description("Waxs scattering as function of time"));

    // Add option for optional output file
    options->addOption(FileNameOption("was_ener").filetype(eftPlot).outputFile()
                           .store(&fnEner_).defaultBasename("ener")
                           .description("Waxs-Debye energy (see manual)"));

    // Add option for optional output file
    options->addOption(FileNameOption("waxs_alpha").filetype(eftPlot).outputFile()
                           .store(&fnAlpha_).defaultBasename("alpha")
                           .description("Excited state population (see manual)"));

    // Add option for selecting a subset of atoms
    //options->addOption(SelectionOption("select").required().valueCount(1)
    //                     .store(&sel_)
    //                     .onlyAtoms());

    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efNoUserPBC);

    // Should this be true or false?
    settings->setPBC(false);
}

static gmx_vsite_t *init_vsite_top(t_topology *top)
{
    int            nvsite, i;
    gmx_vsite_t   *vsite;

    /* check if there are vsites */
    nvsite = 0;
    for (i = 0; i < F_NRE; i++)
    {
        if (interaction_function[i].flags & IF_VSITE)
        {
            nvsite += top->idef.il[i].nr;
        }
    }

    if (nvsite == 0)
    {
        return NULL;
    }

    snew(vsite, 1);

    vsite->n_intercg_vsite = 0;

    /* If we don't have charge groups, the vsite follows its own pbc */
    vsite->nvsite_pbc_molt = 0;

    snew(vsite->vsite_pbc_loc_nalloc, F_VSITEN-F_VSITE2+1);
    snew(vsite->vsite_pbc_loc, F_VSITEN-F_VSITE2+1);

    vsite->nthreads      = 1;
    vsite->th_ind        = NULL;
    vsite->th_ind_nalloc = 0;

    return vsite;
}

void
WaxsDebye::initAnalysis(const TrajectoryAnalysisSettings &settings,
                        const TopologyInformation        &top)
{
    // Add the module that will contain the averaging and the time series
    // for our calculation
    data_.addModule(enerdata_);

    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data et.
    AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
    plotm_->setSettings(settings.plotSettings());
    plotm_->setFileName(fnEner_);
    plotm_->setTitle("WAXS Energy");
    plotm_->setXAxisIsTime();
    plotm_->setYLabel("Energy (kJ/mol)");

    data_.addModule(plotm_);

    // Add a module for plotting the data automatically at the end of
    // the calculation. With this in place you only have to add data
    // points to the data et.
    /*
       AnalysisDataPlotModulePointer plotm_(new AnalysisDataPlotModule());
       plotm_->setSettings(settings.plotSettings());
       plotm_->setFileName(fnAlpha_);
       plotm_->setTitle("WAXS Alphay");
       plotm_->setXAxisIsTime();
       plotm_->setYLabel("Alpha ()");

       data_.addModule(plotm_);
     */
    // Initiate the calculation engine. Think about parallellism!
    t_inputrec ir;

    init_inputrec(&ir);
    ir.waxs.debye_r_min = 0;
    ir.waxs.debye_r_max = 1000;
    ir.waxs.kwaxs       = 1;

    /* Starting values for the optimization of alpha and other run parameters */
    ir.userint1             = 0;
    ir.waxs.debye_alpha_min = 1;
    ir.waxs.debye_alpha_max = 1;
    ir.waxs.nstout          = 1;

    wdf_ = new WaxsDebyeForce(stdout,
                              ((fnSfactor_.length() > 0) ? fnSfactor_.c_str() : NULL),
                              ((fnSqref_.length() > 0) ? fnSqref_.c_str() : NULL),
                              (fnSqdiff_.length() > 0)  ? fnSqdiff_.c_str() : NULL,
                              ((fnSqcalc_.length() > 0) ? fnSqcalc_.c_str() : NULL),
                              NULL,
                              cr_,
                              &ir);

    top_ = top.topology();

    // Count the number of WAXS "bonds"
    nbonds_    = 0;
    iatoms_wd_ = NULL;
    int nn = top_->idef.il[F_WAXS_DEBYE].nr;
    if (nn > 0)
    {
        nbonds_ += nn;
        srenew(iatoms_wd_, nbonds_);
        memcpy(&iatoms_wd_[nbonds_ - nn],
               top_->idef.il[F_WAXS_DEBYE].iatoms,
               nn*sizeof(iatoms_wd_[0]));
    }
    vsite_ = init_vsite_top(top_);
}

void
WaxsDebye::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                        TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle       dh   = pdata->dataHandle(data_);
    //const Selection         &sel  = pdata->parallelSelection(sel_);

    //GMX_RELEASE_ASSERT(NULL != pbc, "You have no periodic boundary conditions");

    // Analysis framework magic
    dh.startFrame(frnr, fr.time);

    // Generate virtual sites if needed
    if (NULL != vsite_)
    {
        construct_vsites(vsite_, fr.x, 1, NULL, top_->idef.iparams,
                         top_->idef.il, epbcNONE, FALSE,
                         cr_, (rvec *) fr.box);
    }
    // Compute the WAXS energy
    snew(f_, fr.natoms);
    double ener = wdf_->calc(stdout,
                             nbonds_,
                             iatoms_wd_,
                             top_->idef.iparams,
                             fr.x,
                             f_,
                             pbc,
                             cr_,
                             NULL);
    sfree(f_);
    // Retrieve present value of Alpha
    double alpha = wdf_->getAlpha();

    // Add the energy to the data set in column 1
    dh.setPoint(0, ener);
    // Add alpha to the data set in column 1
    dh.setPoint(1, alpha);

    // Magic
    dh.finishFrame();
}


void
WaxsDebye::finishAnalysis(gmx_unused int nframes)
{
    please_cite(stdout, "Bjorling2014a");
}

void
WaxsDebye::writeOutput()
{
    // Final results come from statistics module in analysis framework
}

const char WaxsDebyeInfo::name[]             = "waxsdebye";
const char WaxsDebyeInfo::shortDescription[] =
    "Calculate X-ray scattering curves";

TrajectoryAnalysisModulePointer WaxsDebyeInfo::create()
{
    return TrajectoryAnalysisModulePointer(new WaxsDebye);
}

} // namespace analysismodules

} // namespace gmx
