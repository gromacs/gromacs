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
#include "smalloc.h"
#include <gromacs/trajectoryanalysis.h>
#include "gromacs/linearalgebra/gmx_blas.h"
#include "gromacs/linearalgebra/gmx_lapack.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trnio.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/legacyheaders/xvgr.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/do_fit.h"
#include "gromacs/legacyheaders/gmx_fatal.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

using namespace gmx;

/*! \brief
 * Implementation of partial least squares (PLS) based functional mode analysis (FMA).
 */
class FMA : public TrajectoryAnalysisModule
{
    public:
        FMA();

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

        std::string        fnFMAvec_;      /*!File Name for saving FMA result vectors */
        std::string        fnFMAvecEw_;    /*!File Name for saving enesemble-weighted FMA result vectors */
        std::string        fnFMAref_;
        std::string        fnModel_;
        std::string        fnY_;
        bool               mwfit_;
        double             **y_;
        double             *yTrain_;
        double             *yValidate_;
        double             avgy_;
        int                dimY_;
        int                nSamples_;
        int                nFrames_;
        int                nTFrames_;
        int                nVFrames_;
        int                nTrainParam_;
        int                nTrain_;
        int                nValidate_;
        int                nPLSvec_;
        int                total_n_atoms;
        double             *frames_;
        double             *times_;
        int                *frameOrder_;
        double             *avg_;
        double             *w_;
        double             *q_;
        Selection          fitsel_;
        Selection          anasel_;
        rvec               *refstruct_;
        t_atoms            refatoms_;
        real               *fitmass;
        real               *anamass;
        atom_id            *fitindices;
        atom_id            *anaindices;

        double vnorm(double *x, int n);
        double vmean(double *x, int n, int step);
        double corrcoeff(double *x, double *y, int n, int xstep, int ystep);
        void pls_denham(double *x, double *y, int n, int k, int ln, int a, double *w, double *q);


};


/*
 * Vector norm calculation
 *
 * Calculation of the norm of a vector.
 *
 * \param x (n)-vector to calculate the norm of
 * \param n length of n
 */
double FMA::vnorm(double *x, int n)
{
    int i;
    double sn=0.0;

    for (i=0; i<n; i++) {
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
double FMA::vmean(double *x, int n, int step)
{
    int i;
    double m=0.0;

    for (i=0; i<n; i++) {
        m += x[i*step];
    }
    return m / (double) n;
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
double FMA::corrcoeff(double *x, double *y, int n, int xstep, int ystep)
{
    int i;
    double cxy = 0;
    double cxx = 0;
    double cyy = 0;
    double mx,my;

    mx = vmean(x,n,xstep);
    my = vmean(y,n,ystep);

    for (i=0; i<n; i++)
    {
        cxy += (x[i*xstep]-mx) * (y[i*ystep]-my);
        cxx += (x[i*xstep]-mx) * (x[i*xstep]-mx);
        cyy += (y[i*ystep]-my) * (y[i*ystep]-my);
    }

    return cxy / (sqrt(cxx) * sqrt(cyy));
}

/*
 * Partial least squares (PLS) regression
 *
 * Implemented as described in:
 *
 * Denham, M. C. "Implementing Partial Least Squares."
 * Statistics and Computing 5, no. 3 (September 1995): 191–202. doi:10.1007/BF00142661.
 *
 * The algoritm aims to find a vector t to minimize
 *
 * abs(y-X*t)
 *
 * where the vector t=w*q is a linear combination of vectors.
 */
void FMA::pls_denham(double *x, double *y, int n, int k, int ln, int a, double *w, double *q)
{

    /*! \param x  contains the centered (ln,k)-matrix X.
        \param y  contains the centered (n)-vector y.
        \param n  the number of "relevant" rows of the matrix X.
        \param lk the "real" number of rows in the matrix x, or the "leading dimension"
                     in FORTRAN terms. If all structure-value pairs are used for training, this number
                  is equal to n, otherwise bigger.
        \param k  the number of columns (to be used) of the matrix x.
        \param a  < MIN(n-1,k) - the number of PLS factors to include in the regression
                  of y on the matrix X.
        \param w  (output) a (k,a)-matrix containing the coefficient vectors stored by column of
                  the a PLS regression factors obtained from the matrix X on .
        \param q  a (a)-vector containing the least squares regression coefficients of
                  the ordinary least squares regression of y on the PLS factor matrix t.

       Arrays entered into this function should be in FORTRAN-compatible
       column major order.
    */

    double *t;
    snew(t,n*a);
    double *qrt;
    snew(qrt,n*a);
    double *b;
    snew(b,k);
    double *dum;
    snew(dum,n);
    double *rsd;
    snew(rsd,n);

    int i;
    // As FORTRAN functions in C require call by reference, all constants have
    // to be defined explicitly
    char T = 'T';
    char N = 'N';
    double c_done = 1.0;
    double c_dnegone = -1.0;
    double c_dzero = 0.0;
    int c_ione = 1;

    int idum;
    int info = 0;
    int lwork = -1;

    // Note: as BLAS and LAPACK functions are less than self-explaining, this function
    // contains comments that explain what the calls do. The algorithm itself is described
    // in the paper cited above.
    // w_1 = x' * y
    F77_FUNC(dgemv, DGEMV) (    &T, &n, &k, &c_done, x, &ln, y, &c_ione, &c_dzero, w, &c_ione);

    // t_1 = x * w_1
    F77_FUNC(dgemv, DGEMV) (&N,&n,&k,&c_done,x,&ln,w,&c_ione,&c_dzero,t,&c_ione);

    F77_FUNC(dcopy, DCOPY) (&n,t,&c_ione,qrt,&c_ione);
    F77_FUNC(dcopy, DCOPY) (&n,y,&c_ione,rsd,&c_ione);

    // perform regression
    // min(rsd) abs(y - qrt * rsd)
    F77_FUNC(dgels, DGELS) (    &N,      &n, &c_ione, &c_ione, qrt,  &n,  rsd,  &n,  dum, &lwork, &info );

    q[0] = rsd[0];
    F77_FUNC(dcopy, DCOPY) (&n,y,&c_ione,rsd,&c_ione);

    // calculate residuals rsd := y - t_1*q_1
    // using DGEMV          (y := alpha*A*x + beta*y)
    F77_FUNC(dgemv, DGEMV) (    &N,      &n, &c_ione, &c_dnegone,   t,      &n,     q, &c_ione, &c_done, rsd, &c_ione);


    for (i=1; i<a; i++)
    {
        // w_i := X' * rsd
        F77_FUNC(dgemv, DGEMV) (&T, &n, &k, &c_done, x, &ln, rsd, &c_ione, &c_dzero, &w[i*k], &c_ione);
        // t_i := X * w_i
        F77_FUNC(dgemv, DGEMV) (&N, &n, &k, &c_done, x, &ln, &w[i*k], &c_ione, &c_dzero, &t[n*i], &c_ione);

        idum = n*a;
        F77_FUNC(dcopy, DCOPY) (&idum, t, &c_ione, qrt, &c_ione);
        F77_FUNC(dcopy, DCOPY) (&n, y, &c_ione, rsd, &c_ione);

        idum = i + 1;
        // run regression
        F77_FUNC(dgels, DGELS) ( &N, &n, &idum, &c_ione, qrt,  &n, rsd,  &n,  dum, &lwork, &info );

        F77_FUNC(dcopy, DCOPY) (&a,     rsd,&c_ione,     q,&c_ione);

        F77_FUNC(dcopy, DCOPY) (&n,y,&c_ione,rsd,&c_ione);

        // calculate residuals rsd := y - t_1*q_1
        // using DGEMV          (y := alpha*A*x + beta*y)
        idum = i + 1;
        F77_FUNC(dgemv, DGEMV) (    &N,      &n, &idum, &c_dnegone,   t,      &n,     q, &c_ione, &c_done, rsd, &c_ione);

    }

    F77_FUNC(dgemv, DGEMV) (  &N, &k, &a, &c_done,      w, &k, q, &c_ione, &c_dzero,   b, &c_ione);
}

FMA::FMA()
    : TrajectoryAnalysisModule("g_fma", "PLS-based FMA")
{
}


void
FMA::initOptions(Options                    *options,
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
                           .store(&fnY_).defaultBasename("y")
                           .description("Target values for the model"));
    options->addOption(FileNameOption("v")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvec_).defaultBasename("FMAvec.trr")
                           .description("FMA result vectors"));
    options->addOption(FileNameOption("ew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvecEw_).defaultBasename("FMAvec_ew.trr")
                           .description("Ensemble weighted FMA result vectors"));
    options->addOption(FileNameOption("ref")
                           .filetype(eftPDB).outputFile()
                           .store(&fnFMAref_).defaultBasename("FMAref.pdb")
                           .description("Centered reference structure used for fitting."));
    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile().required()
                           .store(&fnModel_).defaultBasename("FMAmodel.xvg")
                           .description("Projection of structures from the input trajectory onto the model."));

    options->addOption(IntegerOption("dim").store(&nPLSvec_).required()
                           .defaultValue(10)
                           .description("Number of PLS-vectors to use."));
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
FMA::initAnalysis(const TrajectoryAnalysisSettings &settings,
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

    int i,j;

    nTrain_ = 0;
    nValidate_ =  0;

    snew(fitindices,fitsel_.atomCount());

    for (i=0; i<fitsel_.atomCount(); i++) {
        fitindices[i] = fitsel_.atomIndices()[i];
    }

    snew(anaindices,anasel_.atomCount());

    for (i=0; i<anasel_.atomCount(); i++) {
        anaindices[i] = anasel_.atomIndices()[i];
    }

    refatoms_ = top.topology()->atoms;

    total_n_atoms = refatoms_.nr;

    snew(fitmass,total_n_atoms);
    for (i=0; i<fitsel_.atomCount(); i++) {
        if (mwfit_) {
            fitmass[fitindices[i]] = refatoms_.atom[fitindices[i]].m;
        }
        else {
            fitmass[fitindices[i]] = 1.0;
        }
    }

    snew(anamass,total_n_atoms);
    for (i=0; i<anasel_.atomCount(); i++) {
        anamass[anaindices[i]] = 1.0;
    }

    top.getTopologyConf(&refstruct_,NULL);
    int ifit = fitsel_.atomCount();
    reset_x(ifit, fitindices,
            total_n_atoms, NULL,
                refstruct_, fitmass);

    nSamples_ = read_xvg(fnY_.c_str(),&y_,&dimY_);
    dimY_--;

    if ((nSamples_ == 0) || (dimY_ == 0)) {
        throw gmx::FileIOError("Error. File " + fnY_ + " did not contain any data.");
    }

    if ((dimY_ > 1) && (nTrainParam_ == -1)) {
        std::cout << "Found at least three columns and no explicit value given\n";
        std::cout << "for -ntrain. Using third column to label structure-value\n";
        std::cout << "pairs for ignore/train/validate.\n\n";
    }
    else if ((dimY_ > 1) && (nTrainParam_ != -1)) {
        gmx_warning("at least three columns found in %s but explicit\n"
                " -ntrain given. Ignoring third column and using first\n"
                "%d structure-value pairs for training.\n\n"
                , fnY_.c_str(), nTrainParam_);
    }


    if (nTrainParam_ != -1) {
        nTrain_ = min(nTrainParam_,nSamples_);
        nValidate_ = nSamples_ - nTrain_;

        snew(yTrain_,nTrain_);
        snew(yValidate_,nValidate_);

        for (i=0; i<nTrain_;i++) {
            yTrain_[i] = y_[1][i];
        }
        for (i=0; i<nValidate_;i++) {
            yValidate_[i] = y_[1][nTrain_+i];
        }
    }
    else if (dimY_ == 1) {
        nTrain_ = nSamples_;
        nValidate_ = 0;

        snew(yTrain_,nTrain_);
        for (i=0; i<nTrain_;i++) {
            yTrain_[i] = y_[1][i];
        }
    }
    else {
        for (i=0; i<nSamples_; i++) {
            if (y_[2][i] == 1) nTrain_++;
            if (y_[2][i] == 2) nValidate_++;
        }

        snew(yTrain_,nTrain_);
        snew(yValidate_,nValidate_);

        j=0;
        for (i=0; i<nSamples_; i++) {
            if (y_[2][i] == 1) {
                yTrain_[j] = y_[1][i];
                j++;
            }
        }

        j=0;
        for (i=0; i<nSamples_; i++) {
            if (y_[2][i] == 2) {
                yValidate_[j] = y_[1][i];
                j++;
            }
        }

    }

    std::cout << "Read "  << nSamples_ << " Lines of "  << dimY_ << " column(s).\n";
    std::cout << "Using " << nTrain_ << " for training and  " << nValidate_ << " for validation.\n";

    if ((nTrain_ == 0)) {
        throw gmx::FileIOError("Error. File " + fnY_ + " did not contain any training-data.");
    }

    for (i=0; i < nTrain_; i++) {
           avgy_ += yTrain_[i];
    }

      avgy_ = avgy_  / (double) nTrain_;

    std::cout << "\nThe average value of y (in the training set) is " << avgy_ << " .\n\n";

    for (i=0; i < nTrain_; i++) {
           yTrain_[i] -= avgy_;
    }
    for (i=0; i < nValidate_; i++) {
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

    snew(frames_, (nTrain_ + nValidate_) * DIM * anasel_.atomCount());

    snew(avg_,DIM * anasel_.atomCount());

    if (frames_ == NULL) {
    throw gmx::InternalError("Failed to initialize memory. Maybe use a shorter trajectory ?");
    }
    nFrames_ = 0;
    nTFrames_ = 0;
    nVFrames_ = 0;

    snew(frameOrder_, (nTrain_ + nValidate_));
    snew(times_, (nTrain_ + nValidate_));

}


void
FMA::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    /**********************************************************************
     * reading in the trajectory                                          *
     *                                                                    *
     * For PLS-FMA, the the coordinates of the whole trajectory have to   *
     * be saved in memory (at least for the atoms that are used in the    *
     * analysis. This is accomplished here. The real analysis is then     *
     * performed once the whole trajectory has been read                  *
     **********************************************************************/

    int i,j;
    int frameNumber;
    int keepFrame = 1;
    int useForAverage = 1;

    double *avgx;
    snew(avgx,DIM);

    if ((nTrainParam_ != -1) || (dimY_ < 2)) {
        if (nTFrames_ < nTrain_) {
            frameNumber = nTFrames_;
            nTFrames_++;
        }
        else {
            frameNumber = nTrain_ + nVFrames_;
            useForAverage = 0;
            nVFrames_++;
        }
    }
    else {
        if (y_[2][nFrames_] == 1) {
            frameNumber = nTFrames_;
            nTFrames_++;
        }
        else if (y_[2][nFrames_] == 2) {
            frameNumber = nTrain_ + nVFrames_;
            useForAverage = 0;
            nVFrames_++;
        }
        else {
            keepFrame = 0;
        }
    }

    if (keepFrame) {
        if ((nTFrames_ + nVFrames_ - 1) == (nTrain_ + nValidate_)) {
            throw gmx::InconsistentInputError("File " + fnY_ + " does not contain data for all frames of the trajectory.");
        }

        // For analysis in pls_denham, the training structures have to be saved in the beginning of
        // the coordinate matrix, while the validation structures (which are not used in pls_denham)
        // have to be in the end. To restore the original order when calculating model values (and
        // correlation coefficients), the position of the frame at <originalPosition> it is stored in
        // frameOrder_[originalPosition].
        frameOrder_[nTFrames_+nVFrames_-1] = frameNumber;

        // center and fit the structure (using reference group)
        reset_x(fitsel_.atomCount(), fitindices,
            total_n_atoms, NULL,
                    fr.x, fitmass);

        do_fit(total_n_atoms, fitmass, refstruct_, fr.x);

        for (i=0; i < anasel_.atomCount(); i++) {
            for (j=0; j < DIM; j++) {

                // already saving the coordinates in FORTRAN compatible COLUMN MAJOR order
                // each ROW describes one structure, while each COLUMN describes one coordinate
                frames_[frameNumber  + ((DIM*i) + j) * (nTrain_ + nValidate_)] = fr.x[anaindices[i]][j];

                times_[frameNumber] = fr.time;
                if (useForAverage) {
                    avg_[(DIM*i) + j] += frames_[frameNumber  + ((DIM*i) + j) * (nTrain_ + nValidate_)]/double(nTrain_);
                }
            }
           }
       }
    nFrames_++;
}


void
FMA::finishAnalysis(int /*nframes*/)
{
    // counting variables:
    int i,j,k;

    nFrames_ = nTFrames_ + nVFrames_;
    nSamples_ = nTrain_ + nValidate_;

    if (nFrames_ < nSamples_)
    {
        throw gmx::InconsistentInputError("Trajectory containts fewer frames than there are datapoints in File " + fnY_);
    }


    // center data
    for (k=0; k<(nSamples_); k++)
        {
        for (i=0; i <anasel_.atomCount(); i++)
            for (j=0; j < DIM; j++) {
                {
                    frames_[k + ((DIM*i) + j) * nSamples_] -= avg_[(DIM*i) + j];
                }
        }
    }

    std::cout << "Running PLS on " << nTrain_ << " training structures (of " << nSamples_ << " total structures)\n";

    // create vector to hold results of PLS
    snew(w_,anasel_.atomCount() * DIM*nPLSvec_);
    snew(q_,nPLSvec_);

    // call PLS function
    pls_denham(frames_,yTrain_, nTrain_, anasel_.atomCount() * DIM, nSamples_, nPLSvec_, w_, q_);


}


void
FMA::writeOutput()
{
    // counting variables:
    int i,j,k;

    output_env_t oenv;
    output_env_init_default(&oenv);

    // the lambda value saved in the vector .trr-file determines if g_anaeig uses mass weighting in
    // the fit or not.
    real lambda;
    lambda = mwfit_ ? 1.0 : 0.0;

    double *PLSvec;
    snew(PLSvec,anasel_.atomCount()*DIM);

    double *mYt;
    snew(mYt,nTrain_);
    double *mYv;
    snew(mYv,nValidate_);

    real **modelY;
    snew(modelY,2);
    for (i=0; i<2; i++) {
        snew(modelY[i],nSamples_);
    }
    double *avgx;
    snew(avgx,DIM);

    // working structure
    rvec *structure;
    snew(structure,anasel_.atomCount());

    //dummy box
    rvec *box;
    snew(box,3);
    for (i=0; i<3; i++) {
        for (j=0; j<DIM; j++) {
            box[i][j] = 0;
        }
    }

    // the "eigenvalues"
    real *values;
    snew(values,(nPLSvec_ + 1));

    // the output files
    t_fileio *outputfile;
    t_fileio *outputfileEw;
    outputfile = open_trn(fnFMAvec_.c_str(), "w");
    outputfileEw = open_trn(fnFMAvecEw_.c_str(), "w");


    // calculate the PLS vector from the PLS result
    for (i=0; i<anasel_.atomCount(); i++) {
        for (j=0; j<DIM; j++) {
            for (k=0; k<nPLSvec_; k++) {
                PLSvec[i*DIM +j] += q_[k] * w_[(i*DIM + j) + (k * anasel_.atomCount() * DIM)];
            }
        }
    }

    for (i=0; i<DIM; i++) {
        avgx[i] = vmean(&PLSvec[i], anasel_.atomCount(),DIM);
    }

    if (!fnFMAref_.empty())
    {
        write_sto_conf(fnFMAref_.c_str(), "Reference structure used in FMA",
                        &refatoms_,
                        refstruct_, NULL, -1, box);
    }

    // write the reference structure to the vector-file
    for (i=0; i<anasel_.atomCount(); i++) {
        for (j=0; j<DIM; j++) {
            structure[i][j] = refstruct_[anasel_.atomIndices()[i]][j];
        }
    }

    // To account for a possible non-zero mean of the training values of y, the average
    // structure written to the vector file is shifted
    for (i=0; i<anasel_.atomCount(); i++) {
        for (j=0; j<DIM; j++) {
            structure[i][j] = avg_[i*DIM + j] - avgy_/(PLSvec[i*DIM +j])/anasel_.atomCount()/DIM;
        }
    }

    fwrite_trn(outputfile, 0, 0.0, lambda,
                    box, anasel_.atomCount(), structure, NULL, NULL);
    fwrite_trn(outputfileEw, 0, 0.0, lambda,
                    box, anasel_.atomCount(), structure, NULL, NULL);


    // Write out PLS result vector (linear combination) to one file
    for (i=0; i<anasel_.atomCount(); i++) {
        for (j=0; j<DIM; j++) {
            structure[i][j] = PLSvec[i*DIM +j];// - avgx[j];
        }
    }

    fwrite_trn(outputfile, 1, 1.0, lambda,
                    box, anasel_.atomCount(), structure, NULL, NULL);
    // Write the component vectors to the ensemble weighted file
    for (k=0; k<nPLSvec_; k++) {
        for (i=0; i<anasel_.atomCount(); i++) {
            for (j=0; j<DIM; j++) {
                structure[i][j] = q_[k] * w_[(i*DIM + j) + (k * anasel_.atomCount() * DIM)];
            }
        }
        fwrite_trn(outputfileEw, k+1, q_[k], lambda,
                        box, anasel_.atomCount(), structure, NULL, NULL);
    }


    close_trn(outputfile);
    close_trn(outputfileEw);

    //calculate the model-values

    for (i=0; i<nSamples_; i++)
    {
        modelY[0][i] = times_[frameOrder_[i]];
        modelY[1][i] = avgy_;
        for(j=0; j<anasel_.atomCount()*DIM; j++)
        {
            modelY[1][i] += PLSvec[j] * (frames_[frameOrder_[i]  + j * nSamples_]);
        }
        if (frameOrder_[i]<nTrain_) {
            mYt[frameOrder_[i]] = modelY[1][i];
        }
        else
        {
            mYv[frameOrder_[i]-nTrain_] = modelY[1][i];
        }

    }

    std::cout << "r(training)   = " << corrcoeff(yTrain_, mYt, nTrain_, 1, 1) <<  ".\n";
    std::cout << "r(validation) = " << corrcoeff(yValidate_, mYv, nValidate_, 1, 1) <<  ".\n";

    write_xvg(fnModel_.c_str(), "Projection of structure onto model.", nSamples_, 2, modelY,
                   NULL, oenv);
}

}       // namespace


const char FMAInfo::name[]             = "fma";
const char FMAInfo::shortDescription[] =
    "Functional Mode Analysis";

TrajectoryAnalysisModulePointer FMAInfo::create()
{
    return TrajectoryAnalysisModulePointer(new FMA);
}

} // namespace analysismodules

} // namespace gmx
