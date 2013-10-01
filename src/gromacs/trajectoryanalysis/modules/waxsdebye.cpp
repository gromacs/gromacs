/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#include "gromacs/legacyheaders/pbc.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/smalloc.h"
#include "gromacs/legacyheaders/atomprop.h"
#include "gromacs/legacyheaders/copyrite.h"
#include "gromacs/legacyheaders/vsite.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

/*! \file\internal\brief
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


        //! Interaction definition for vsites
        t_ilist                          *ilist_vs_;

        //! Number of force field parameters
        int                               ntypes_;

        //! Force field parameters
        t_iparams                        *iparams_;

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
    ilist_vs_   = NULL;
    ntypes_     = 0;
    iparams_    = NULL;
    cr_         = init_commrec();
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
    if (NULL != ilist_vs_)
    {
        if (ilist_vs_->nr > 0)
        {
            sfree(ilist_vs_->iatoms);
        }
        sfree(ilist_vs_);
    }
    if (NULL != iparams_)
    {
        sfree(iparams_);
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

    const t_topology *topology = top.topology();
    nbonds_ = topology->idef.il[F_WAXS_DEBYE].nr;
    snew(iatoms_wd_, nbonds_);
    memcpy(iatoms_wd_, topology->idef.il[F_WAXS_DEBYE].iatoms,
           nbonds_*sizeof(iatoms_wd_[0]));
    ntypes_ = topology->idef.ntypes;

    snew(ilist_vs_, 1);
    ilist_vs_->nr              = topology->idef.il[F_VSITEN].nr;
    ilist_vs_->nr_nonperturbed = topology->idef.il[F_VSITEN].nr_nonperturbed;
    ilist_vs_->nalloc          = topology->idef.il[F_VSITEN].nalloc;
    // Check consistency
    GMX_RELEASE_ASSERT((ilist_vs_->nr > ilist_vs_->nalloc),
                       "Inconsistency in VSITEN data structure");
    if (0 < ilist_vs_->nalloc)
    {
        snew(ilist_vs_->iatoms, ilist_vs_->nalloc);
        memcpy(ilist_vs_->iatoms, topology->idef.il[F_VSITEN].iatoms,
               ilist_vs_->nalloc*sizeof(ilist_vs_->iatoms[0]));
    }
    gmx_mtop_t *mtop = NULL; //top.mtop();
    vsite_ = init_vsite(mtop, cr_, TRUE);

    ntypes_ = topology->idef.ntypes;
    snew(iparams_, ntypes_);
    memcpy(iparams_, topology->idef.iparams,
           ntypes_*sizeof(iparams_[0]));
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
    if (0 < ilist_vs_->nr)
    {
        matrix box;
        for (int i = 0; (i < DIM); i++)
        {
            for (int j = 0; (j < DIM); j++)
            {
                box[i][j] = fr.box[i][j];
            }
        }
        construct_vsites(vsite_,
                         fr.x,
                         0,
                         NULL,
                         iparams_,
                         ilist_vs_,
                         pbc->ePBC,
                         FALSE,
                         (t_graph *)NULL,
                         cr_,
                         box);
    }
    // Compute the WAXS energy
    double ener = wdf_->calc(stdout,
                             nbonds_,
                             iatoms_wd_,
                             iparams_,
                             fr.x,
                             f_,
                             pbc,
                             cr_,
                             NULL);

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
WaxsDebye::finishAnalysis(int nframes)
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
