/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::FMA.
 *
 * \author Jan Henning Peters <JanHPeters@gmx.net>
 * \ingroup module_trajectoryanalysis
 */

#include "fma.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include "gromacs/utility/smalloc.h"
#include <gromacs/trajectoryanalysis.h>
#include "gromacs/linearalgebra/pls.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/legacyheaders/copyrite.h"

namespace gmx
{

namespace analysismodules
{

namespace
{


/*! \brief
 * Implementation of partial least squares (PLS) based functional mode analysis (FMA).
 */
class Fma : public TrajectoryAnalysisModule
{
    public:
        Fma();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;

        std::string            fnFMAvec_;    /*!File Name for saving FMA result vectors */
        std::string            fnFMAvecEw_;  /*!File Name for saving enesemble-weighted FMA result vectors */
        std::string            fnFMAmode_;
        std::string            fnFMAmodeEw_;
        std::string            fnFMAref_;
        std::string            fnModel_;
        std::string            fnY_;
        bool                   mwfit_;
        FMatrix                y_;
        RVector                yTrain_;
        RVector                yValidate_;
        real                   avgy_;
        int                    dimY_;
        int                    nSamples_;
        int                    nFrames_;
        int                    nTFrames_;
        int                    nVFrames_;
        int                    nTrainParam_;
        int                    nTrain_;
        int                    nValidate_;
        int                    nPLSvec_;
        int                    total_n_atoms;
        int					   nModeFrames_;
        FMatrix                fframes_;
        RVector                times_;
        std::vector<int>       frameOrder_;
        RVector                avg_;
        FMatrix                w_;
        RVector                q_;
        Selection              fitsel_;
        Selection              anasel_;
        rvec                  *refstruct_;
        t_atoms                refatoms_;
        RVector                fitmass;
        std::vector<atom_id>   fitindices;
        std::vector<atom_id>   anaindices;

        real vnorm(real *x, int n);
        real vmean(real *x, int n, int step);
        real corrcoeff(RVector x, RVector y, int n, int xstep, int ystep);


};

Fma::Fma()
    : TrajectoryAnalysisModule("fma", "PLS-based FMA")
{
    /* Initialization here is rather pointless, as it happens in
     * initAnalysis - but it's required by Jenkins
     */
    mwfit_        = false;
    avgy_         = 0.0;
    dimY_         = 0;
    nSamples_     = 0;
    nFrames_      = 0;
    nTFrames_     = 0;
    nVFrames_     = 0;
    nTrainParam_  = 0;
    nTrain_       = 0;
    nValidate_    = 0;
    nPLSvec_      = 0;
    total_n_atoms = 0;
    nModeFrames_  = 0;
    refstruct_    = NULL;
    init_t_atoms(&refatoms_,0,FALSE);
}

void
Fma::initOptions(Options                    *options,
                 TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Implementation of partial least squares (PLS) based functional mode",
        "analysis (FMA). The tool requires a trajectory and a text file",
        "containing a single scalar value for each structure. It will then use",
        "PLS to identify the collective mode that best models the scalar value.",
        "\n",
        "The text file (-y) must contain two or three columns:",
        "First, the timestamp of the structure (which will not be used).",
        "Second, the scalar value to be modeled.",
        "Third (optional), a flag how to handle the structure-value pair ",
        "(0 - ignore, 1 - use for training, 2 - use for validation.\n"
        "If a value is given for -ntrain, a third column will be ignored,",
        "if it is not given and no third column exists, all structure-value",
        "pairs will be used for training."
    };

    options->setDescription(desc);

    options->addOption(FileNameOption("y")
                           .filetype(eftPlot).inputFile().required()
                           .store(&fnY_).defaultBasename("y").defaultExtension(".xvg")
                           .description("Target values for the model"));
    options->addOption(FileNameOption("v")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvec_).defaultBasename("FMAvec")
                           .defaultExtension(".trr").description("FMA result vectors"));
    options->addOption(FileNameOption("ew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvecEw_).defaultBasename("FMAvec_ew")
                           .defaultExtension(".trr")
                           .description("Ensemble weighted FMA result vectors"));
    options->addOption(FileNameOption("om")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAmode_).defaultBasename("FMAmode")
                           .defaultExtension(".pdb").description("FMA mode"));
    options->addOption(FileNameOption("omew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAmodeEw_).defaultBasename("FMAmode_ew")
                           .defaultExtension(".pdb").description("Ensemble weighted FMA mode"));
    options->addOption(FileNameOption("ref")
                           .filetype(eftPDB).outputFile().required()
                           .store(&fnFMAref_).defaultBasename("FMAref")
                           .description("Centered reference structure used for fitting."));
    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile().required()
                           .store(&fnModel_).defaultBasename("FMAmodel")
                           .description("Projection of structures from the input trajectory onto the model."));

    options->addOption(IntegerOption("dim").store(&nPLSvec_).required()
                           .defaultValue(10)
                           .description("Number of PLS-vectors to use."));
    options->addOption(IntegerOption("nframes").store(&nModeFrames_).required()
                           .defaultValue(20)
                           .description("Number of frames for FMAmode output."));
    options->addOption(IntegerOption("ntrain").store(&nTrainParam_)
                           .defaultValue(-1).required()
                           .description("Use the first ntrain structures for training, all others for validation."));

    options->addOption(BooleanOption("mwfit").store(&mwfit_).required()
                           .defaultValue(FALSE).description("Use mass-weighted fit."));

    options->addOption(SelectionOption("fit")
                           .store(&fitsel_).required()
                           .description("Group for least squares fitting before analysis"));
    options->addOption(SelectionOption("analysis")
                           .store(&anasel_).required().onlyStatic()
                           .description("Group to use in analysis. Only static selections allowed."));


    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}


void
Fma::initAnalysis(const TrajectoryAnalysisSettings   & /*settings*/,
                  const TopologyInformation         &top)
{
    /********************************************************************
     * reading y data                                                   *
     *                                                                  *
     * The y data of the model is saved in a .xvg text file, which is   *
     * supplied with the -y parameter. Comment lines are to be ignored. *
     * While reading, the number of data points is counted (it should   *
     * correspond with the number of frames in the trajectory) and the  *
     * dimensionality of the data is determined as the smallest number  *
     * of columns (ignoring the first columns, which represents the     *
     * time) in the file. Warnings are printed out if the number of     *
     * columns is not consistent over the length of the file and also   *
     * (at the moment) if the dimensionality his higher than one as the *
     * current implementation can only handle one.                      *
     ********************************************************************/

    int      i, j;
    double **ydata;

    nTrain_    = 0;
    nValidate_ =  0;

    fitindices.resize(fitsel_.atomCount());

    for (i = 0; i < fitsel_.atomCount(); i++)
    {
        fitindices[i] = fitsel_.atomIndices()[i];
    }

    anaindices.resize(anasel_.atomCount());
    init_t_atoms(&refatoms_,anasel_.atomCount(),FALSE);

    for (i = 0; i < anasel_.atomCount(); i++)
    {
        anaindices[i] = anasel_.atomIndices()[i];

        //refatoms_.atom[i] = top.topology()->atoms.atom[anaindices[i]];
        //refatoms_.atomname[i] = top.topology()->atoms.atomname[anaindices[i]];
        //refatoms_.atomtype[i] = top.topology()->atoms.atomtype[anaindices[i]];
        //refatoms_.atomtypeB[i] = top.topology()->atoms.atomtypeB[anaindices[i]];
    }

    refatoms_ = top.topology()->atoms;

    total_n_atoms = top.topology()->atoms.nr;

    fitmass.resize(total_n_atoms);
    for (i = 0; i < fitsel_.atomCount(); i++)
    {
        if (mwfit_)
        {
            fitmass[fitindices[i]] = refatoms_.atom[fitindices[i]].m;
        }
        else
        {
            fitmass[fitindices[i]] = 1.0;
        }
    }

    top.getTopologyConf(&refstruct_, NULL);
    int ifit = fitsel_.atomCount();
    reset_x(ifit, &fitindices.front(),
            total_n_atoms, NULL,
            refstruct_, &fitmass.front());

    nSamples_ = read_xvg(fnY_.c_str(), &ydata, &dimY_);

    y_.resize(dimY_,nSamples_);

    dimY_--;

    if ((nSamples_ == 0) || (dimY_ == 0))
    {
        GMX_THROW(FileIOError("Error. File " + fnY_ + " did not contain any data."));
    }

    for (i = 0; i < dimY_+1; i++)
    {
        for (j = 0; j < nSamples_; j++)
        {
            y_(i,j) = ydata[i][j];
        }
    }

    if ((dimY_ > 1) && (nTrainParam_ == -1))
    {
        printf("Found at least three columns and no explicit value given\n");
        printf("for -ntrain. Using third column to label structure-value\n");
        printf("pairs for ignore/train/validate.\n\n");
    }
    else if ((dimY_ > 1) && (nTrainParam_ != -1))
    {
        gmx_warning("at least three columns found in %s but explicit\n"
                    " -ntrain given. Ignoring third column and using first\n"
                    "%d structure-value pairs for training.\n\n",
                    fnY_.c_str(), nTrainParam_);
    }


    if (nTrainParam_ != -1)
    {
        nTrain_    = std::min(nTrainParam_, nSamples_);
        nValidate_ = nSamples_ - nTrain_;

        yTrain_.resize(nTrain_);
        yValidate_.resize(nValidate_);

        for (i = 0; i < nTrain_; i++)
        {
            yTrain_[i] = y_(1,i);
        }
        for (i = 0; i < nValidate_; i++)
        {
            yValidate_[i] = y_(1,nTrain_+i);
        }
    }
    else if (dimY_ == 1)
    {
        nTrain_    = nSamples_;
        nValidate_ = 0;

        yTrain_.resize(nTrain_);
        for (i = 0; i < nTrain_; i++)
        {
            yTrain_[i] = y_(1,i);
        }
    }
    else
    {
        for (i = 0; i < nSamples_; i++)
        {
            if (y_(2,i) == 1)
            {
                nTrain_++;
            }
            if (y_(2,i) == 2)
            {
                nValidate_++;
            }
        }

        yTrain_.resize(nTrain_);
        yValidate_.resize(nValidate_);

        j = 0;
        for (i = 0; i < nSamples_; i++)
        {
            if (y_(2,i) == 1)
            {
                yTrain_[j] = y_(1,i);
                j++;
            }
        }

        j = 0;
        for (i = 0; i < nSamples_; i++)
        {
            if (y_(2,i) == 2)
            {
                yValidate_[j] = y_(1,i);
                j++;
            }
        }

    }

    printf("Read %d Lines of %d column(s).\n", nSamples_, dimY_);
    printf("Using %d for training and %d for validation.\n", nTrain_, nValidate_);

    if (nTrain_ == 0)
    {
        GMX_THROW(gmx::FileIOError("Error. File " + fnY_ + " did not contain any training-data."));
    }

    for (i = 0; i < nTrain_; i++)
    {
        avgy_ += yTrain_[i];
    }

    avgy_ = avgy_  / (real) nTrain_;

    printf("\nThe average value of y (in the training set) is %.3f .\n\n",avgy_);

    for (i = 0; i < nTrain_; i++)
    {
        yTrain_[i] -= avgy_;
    }
    for (i = 0; i < nValidate_; i++)
    {
        yValidate_[i] -= avgy_;
    }


    /**********************************************************************
     * initialize memory for the trajectory                               *
     *                                                                    *
     * from the input file, we know how many frames there are to be       *
     * expected, so now we can reserve the appropriate ammount of memory, *
     * which is a                                                         *
     * (3 * nAtoms) by (nFrames)                                          *
     * array. If this prediction does not hold (because the trajectory is *
     * either too long or too short), something went wrong and the tool   *
     * exits with an error.                                               *
     **********************************************************************/

    fframes_.resize(nTrain_+nValidate_, DIM * anasel_.atomCount());

    avg_.resize(DIM * anasel_.atomCount());

    /*
    if (frames_ == NULL)
    {
        GMX_THROW(gmx::InternalError("Failed to initialize memory. Maybe use a shorter trajectory ?"));
    }*/
    nFrames_  = 0;
    nTFrames_ = 0;
    nVFrames_ = 0;

    frameOrder_.resize(nTrain_ + nValidate_);
    times_.resize(nTrain_ + nValidate_);

}


void
Fma::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /*pbc*/,
                  TrajectoryAnalysisModuleData * /*pdata*/)
{
    /**********************************************************************
     * reading in the trajectory                                          *
     *                                                                    *
     * For PLS-FMA, the the coordinates of the whole trajectory have to   *
     * be saved in memory (at least for the atoms that are used in the    *
     * analysis. This is accomplished here. The real analysis is then     *
     * performed once the whole trajectory has been read                  *
     **********************************************************************/

    int   i, j;
    int   frameNumber;
    int   keepFrame     = 1;
    int   useForAverage = 1;

    if ((nTrainParam_ != -1) || (dimY_ < 2))
    {
        if (nTFrames_ < nTrain_)
        {
            frameNumber = nTFrames_;
            nTFrames_++;
        }
        else
        {
            frameNumber   = nTrain_ + nVFrames_;
            useForAverage = 0;
            nVFrames_++;
        }
    }
    else
    {
        if (y_(2,nFrames_) == 1)
        {
            frameNumber = nTFrames_;
            nTFrames_++;
        }
        else if (y_(2,nFrames_) == 2)
        {
            frameNumber   = nTrain_ + nVFrames_;
            useForAverage = 0;
            nVFrames_++;
        }
        else
        {
            keepFrame = 0;
        }
    }

    if (keepFrame)
    {
        if ((nTFrames_ + nVFrames_ - 1) == (nTrain_ + nValidate_))
        {
            GMX_THROW(gmx::InconsistentInputError("File " + fnY_ + " does not contain data for all frames of the trajectory."));
        }

        /* For analysis in pls_denham, the training structures have to be saved in the beginning of
           the coordinate matrix, while the validation structures (which are not used in pls_denham)
           have to be in the end. To restore the original order when calculating model values (and
           correlation coefficients), the position of the frame at <originalPosition> it is stored in
           frameOrder_[originalPosition].*/
        frameOrder_[nTFrames_+nVFrames_-1] = frameNumber;

        /* center and fit the structure (using reference group) */
        reset_x(fitsel_.atomCount(), &fitindices.front(),
                total_n_atoms, NULL,
                fr.x, &fitmass.front());

        do_fit(total_n_atoms, &fitmass.front(), refstruct_, fr.x);

        for (i = 0; i < anasel_.atomCount(); i++)
        {
            for (j = 0; j < DIM; j++)
            {

                /* already saving the coordinates in FORTRAN compatible COLUMN MAJOR order
                   each ROW describes one structure, while each COLUMN describes one coordinate */
                fframes_(frameNumber,(DIM*i) + j) = fr.x[anaindices[i]][j];

                times_[frameNumber] = fr.time;
                if (useForAverage)
                {
                    avg_[(DIM*i) + j] += fframes_(frameNumber, (DIM*i) + j)/real(nTrain_);
                }
            }
        }
    }
    nFrames_++;
}


void
Fma::finishAnalysis(int /*nframes*/)
{
    int i, j, k;

    nFrames_  = nTFrames_ + nVFrames_;
    nSamples_ = nTrain_ + nValidate_;

    if (nFrames_ < nSamples_)
    {
        GMX_THROW(gmx::InconsistentInputError("Trajectory containts fewer frames than there are datapoints in File " + fnY_));
    }


    /* center data */
    for (k = 0; k < (nSamples_); k++)
    {
        for (i = 0; i < anasel_.atomCount(); i++)
        {
            for (j = 0; j < DIM; j++)
            {
                {
                    fframes_(k, (DIM*i) + j) -= avg_[(DIM*i) + j];
                }
            }
        }
    }

    printf("Running PLS on %d training structures (of %d total structures)\n",nTrain_,nSamples_);

    w_.resize(anasel_.atomCount() * DIM,nPLSvec_);
    q_.resize(nPLSvec_);

    pls_denham(&fframes_, &yTrain_, nTrain_, anasel_.atomCount() * DIM, nSamples_, nPLSvec_, &w_, &q_);

    please_cite(stdout, "HubdeGroot2009");
    please_cite(stdout, "KrivobokovaEtAl2012");
}


void
Fma::writeOutput()
{
    int          i, j, k;
    int          nonZeroPLSvectorComponents = 0;
    real 		 minProjection = GMX_REAL_MAX;
    real 		 maxProjection = GMX_REAL_MIN;
    real         v2 = 0.0;
    real         vEw2 = 0.0;

    output_env_t oenv;
    output_env_init_default(&oenv);

    /* the lambda value saved in the vector .trr-file determines if g_anaeig uses mass weighting in
       the fit or not.*/
    real lambda;
    lambda = mwfit_ ? 1.0 : 0.0;

    RVector PLSvec(anasel_.atomCount()*DIM);

    RVector mYt(nTrain_);
    RVector mYv(nValidate_);

    std::vector<std::vector<real>> modelY(2,std::vector<real> (nSamples_));

    t_trxframe  fr;
    clear_trxframe(&fr, TRUE);
    fr.bAtoms = TRUE;
    fr.atoms = &refatoms_;

    rvec *structure;
    snew(structure, anasel_.atomCount());

    rvec *box;
    snew(box, 3);
    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < DIM; j++)
        {
            box[i][j] = 0;
            fr.box[i][j] = 0;
        }
    }

    t_fileio *outputfile;
    t_fileio *outputfileEw;
    t_trxstatus *modefile;
    t_trxstatus *modefileEw;

    outputfile   = open_trn(fnFMAvec_.c_str(), "w");
    outputfileEw = open_trn(fnFMAvecEw_.c_str(), "w");
    modefile   = open_trx(fnFMAmode_.c_str(), "w");
    modefileEw = open_trx(fnFMAmodeEw_.c_str(), "w");


    for (i = 0; i < anasel_.atomCount(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            for (k = 0; k < nPLSvec_; k++)
            {
                PLSvec[i*DIM +j] += q_[k] * w_((i*DIM + j),k);
            }
        }
    }

    if (!fnFMAref_.empty())
    {
        write_sto_conf(fnFMAref_.c_str(), "Reference structure used in FMA",
                       &refatoms_, refstruct_, NULL, -1, box);
}

    for (i = 0; i < anasel_.atomCount(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            structure[i][j] = refstruct_[anasel_.atomIndices()[i]][j];
        }
    }

    for (i = 0; i < anasel_.atomCount(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            if ( (PLSvec[i*DIM + j] > GMX_REAL_EPS) || (PLSvec[i*DIM + j] < - GMX_REAL_EPS) )
            {
            	nonZeroPLSvectorComponents++;
            }
        }
    }

    for (i = 0; i < anasel_.atomCount(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
        	structure[i][j] = avg_[i*DIM + j];
            if ( (PLSvec[i*DIM + j] > GMX_REAL_EPS) || (PLSvec[i*DIM + j] < - GMX_REAL_EPS) )
            {
                /* To account for a possible non-zero mean of the training values of y, the average
                   structure written to the vector file is shifted for all coordinate corresponding
                   to non-zero PLS vector components*/
            	structure[i][j] = structure[i][j] - avgy_/(PLSvec[i*DIM +j])/nonZeroPLSvectorComponents;
            }

        }
    }

    fwrite_trn(outputfile, 0, 0.0, lambda,
               box, anasel_.atomCount(), structure, NULL, NULL);
    fwrite_trn(outputfileEw, 0, 0.0, lambda,
               box, anasel_.atomCount(), structure, NULL, NULL);


    for (i = 0; i < anasel_.atomCount(); i++)
    {
        for (j = 0; j < DIM; j++)
        {
            structure[i][j] = PLSvec[i*DIM +j];
            v2 = v2 + (PLSvec[i*DIM +j] * PLSvec[i*DIM +j]);
        }
    }

    fwrite_trn(outputfile, 1, 1.0, lambda,
               box, anasel_.atomCount(), structure, NULL, NULL);

    for (k = 0; k < nPLSvec_; k++)
    {
        for (i = 0; i < anasel_.atomCount(); i++)
        {
            for (j = 0; j < DIM; j++)
            {
                structure[i][j] = q_[k] * w_((i*DIM + j),k);
                vEw2 = vEw2 + structure[i][j] * structure[i][j];
            }
        }
        fwrite_trn(outputfileEw, k+1, q_[k], lambda,
                   box, anasel_.atomCount(), structure, NULL, NULL);
    }


    close_trn(outputfile);
    close_trn(outputfileEw);

    for (i = 0; i < nSamples_; i++)
    {
        modelY[0][i] = times_[frameOrder_[i]];
        modelY[1][i] = avgy_;
        for (j = 0; j < anasel_.atomCount()*DIM; j++)
        {
            modelY[1][i] += PLSvec[j] * (fframes_(frameOrder_[i],j));
        }
        if (frameOrder_[i] < nTrain_)
        {
            mYt[frameOrder_[i]] = modelY[1][i];
        }
        else
        {
            mYv[frameOrder_[i]-nTrain_] = modelY[1][i];
        }

        if (modelY[1][i] > maxProjection)
        {
        	maxProjection = modelY[1][i];
        }

        if (modelY[1][i] < minProjection)
        {
        	minProjection = modelY[1][i];
        }

    }

    printf("minimal projection = %8.3f\n", minProjection);
    printf("maximal projection = %8.3f\n", maxProjection);


    /*
     * Create a trajectory to visualize the FMA mode
     */

    fr.natoms = total_n_atoms;
    fr.lambda = lambda;


    for (k=0; k < nModeFrames_; k++)
    {
        for (i = 0; i < anasel_.atomCount(); i++)
        {
            for (j = 0; j < DIM; j++)
            {
                refstruct_[anasel_.atomIndices()[i]][j] = avg_[i*DIM + j] + (minProjection + (maxProjection-minProjection)/nModeFrames_*k - avgy_) / v2 * PLSvec[i*DIM +j];
            }
        }

        fr.time = k+1;
        fr.step = k+1;
        fr.bX = TRUE;
        fr.x = refstruct_;

        write_trxframe_indexed(modefile,&fr, anasel_.atomCount(),
        		&anaindices.front(), NULL);
    }

    close_trx(modefile);

    for (k=0; k < nModeFrames_; k++)
    {
        for (i = 0; i < anasel_.atomCount(); i++)
        {
            for (j = 0; j < DIM; j++)
            {
                refstruct_[anasel_.atomIndices()[i]][j] = avg_[i*DIM + j] + (minProjection + (maxProjection-minProjection)/nModeFrames_*k - avgy_) / vEw2 * q_[0] * w_((i*DIM + j),0);;
            }
        }

        fr.time = k+1;
        fr.step = k+1;
        fr.bX = TRUE;
        fr.x = refstruct_;

        write_trxframe_indexed(modefileEw,&fr, anasel_.atomCount(),
        		&anaindices.front(), NULL);
    }

    close_trx(modefileEw);


    real cct = corrcoeff(yTrain_, mYt, nTrain_, 1, 1);
    printf("r(training)   = %8.3f\n",cct);
    if (nValidate_ > 0)
    {
    	real ccv = corrcoeff(yValidate_, mYv, nValidate_, 1, 1);
    	printf("r(validation) = %8.3f\n", ccv);
    }

    FILE *modelout;
    modelout = xvgropen(fnModel_.c_str(), "Projection of structure onto model.", "X", "Y", oenv);
    for (i = 0; (i < nSamples_); i++)
    {
    	fprintf(modelout, " %12.5e %12.5e\n", modelY[0][i], modelY[1][i]);
    }
    xvgrclose(modelout);
}

/*
 * Vector norm calculation
 *
 * Calculation of the norm of a vector.
 *
 * \param x (n)-vector to calculate the norm of
 * \param n length of n
 */
real Fma::vnorm(real *x, int n)
{
    int  i;
    real sn = 0.0;

    for (i = 0; i < n; i++)
    {
        sn += x[i]*x[i];
    }

    return sqrt(sn);
}

/*
 * Mean value calculation
 *
 * Calculation of the mean value of a vector.
 *
 * \param x (n*step)-vector to calculate the mean of
 * \param n number of elements to use for averaging
 * \param step use only every step-th element of x
 */
real Fma::vmean(real *x, int n, int step)
{
    int  i;
    real m = 0.0;

    for (i = 0; i < n; i++)
    {
        m += x[i*step];
    }
    return m / (real) n;
}

/*
 * Calculation of the correlation coefficient between two vectors
 *
 * Simple calculation of the Pearson product-moment correlation coefficient between two
 * (n)-vectory x and y
 *
 * \param x first (n*step)-vector
 * \param y second (n*step)-vector
 * \param number of values to use in calculation
 * \param xstep use only every step-th element of x
 * \param ystep use only every step-th element of y
 */
real Fma::corrcoeff(RVector x, RVector y, int n, int xstep, int ystep)
{
    int  i;
    real cxy = 0;
    real cxx = 0;
    real cyy = 0;
    real mx, my;

    mx = vmean(&x.front(), n, xstep);
    my = vmean(&y.front(), n, ystep);

    for (i = 0; i < n; i++)
    {
        cxy += (x[i*xstep]-mx) * (y[i*ystep]-my);
        cxx += (x[i*xstep]-mx) * (x[i*xstep]-mx);
        cyy += (y[i*ystep]-my) * (y[i*ystep]-my);
    }

    return cxy / (sqrt(cxx) * sqrt(cyy));
}

}

const char FMAInfo::name[]             = "fma";
const char FMAInfo::shortDescription[] =
    "Functional Mode Analysis";

TrajectoryAnalysisModulePointer FMAInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Fma);
}

}

}
