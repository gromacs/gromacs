/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2017,2018,2019, by the GROMACS development team, led by
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
 * \author Berenger Bramas <berenger.bramas@mpcdf.mpg.de>
 * \author R. Thomas Ullmann <thomas.ullmann@mpibpc.mpg.de>
 * \author Luka Stanisic <luka.stanisic@mpcdf.mpg.de>
 *
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "fma.h"

#include <math.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/trajectoryanalysis.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trrio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/linearalgebra/partial_least_squares/partial_least_squares.h"
#include "gromacs/math/do_fit.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"

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

        void initOptions(IOptionsContainer                  *options,
                         TrajectoryAnalysisSettings         *settings) final;
        void initAnalysis(const TrajectoryAnalysisSettings  &settings,
                          const TopologyInformation         &top) final;

        void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata) final;

        void finishAnalysis(int nframes) final;
        void writeOutput() final;

        //! the allocator type for allocating real
        using AllocatorClass   = gmx::SimdSetup<real>::allocator_type;
        //! the matrix datatype to be used in pls_denham
        using MatrixClass      = gmx::IrregArray2D<real, AllocatorClass>;
        //! the vector datatype to be used in pls_denham
        using VectorClass = std::vector<real, AllocatorClass>;
        //! the datatype to use to represent indices
        using index_type = int64_t;
        //! the datatype to use to represent sizes
        using size_type = size_t;
        //! type of vector for storing RVecs, can be treated as rvec array
        using RVecVector = gmx::PaddedVector<gmx::RVec>;
    private:
        class ModuleData;

        //! filename for saving FMA result vectors
        std::string                             fnFMAvec_;
        //! filename for saving the ensemble-weighted FMA result vector
        std::string                             fnFMAvecEw_;
        //! filename for saving the FMA result vector in anaeig compatible format
        std::string                             fnFMAvecAnaeig_;
        //! filename for saving the ensemble-weighted FMA result vector in anaeig compatible format
        std::string                             fnFMAvecEwAnaeig_;
        //! filename for saving the FMA result vector in terms of the separate PLS vectors of order 1 to k
        std::string                             fnFMAmode_;
        //! filename for saving the ensemble-weighted FMA result vector (first PLS component only)
        std::string                             fnFMAmodeEw_;
        //! the reference structure used for least-squares superposition of the trajectory frames before FMA
        std::string                             fnFMAref_;
        //! filename for saving the projection of the input structures on the FMA model (in the ideal case equal to the input observable)
        std::string                             fnModel_;
        //! filename for saving the projection of the input structures on the ensemble-weighted FMA model (first PLS component only)
        std::string                             fnModelEw_;
        //! filename for reading the input observable values
        std::string                             fnY_;
        //! if true, perform mass-weighted least-squares fit of the trajectory frames onto the reference structure
        bool                                    mwfit_  = false;
        //! if true, optimize the number of PLS components to obtain the maximum correlation of the predicted and observed observable values for the validation set
        bool                                    optDim_ = false;
        //! if true, shift the average structure in the output vector files such that it corresponds to the minimum observable value in the training set
        bool                                    shiftAveStructToMinimumObs_ = false;
        //! the input observable values and their assignement to training or validation set
        MatrixClass                             y_;
        //! the observable values of the training set
        VectorClass                             yTrain_;
        //! the observable values of the validation set
        VectorClass                             yValidate_;
        //! the average value of the input observable in the training set
        real                                    avgy_        = 0.0;
        //! the minimum value of the input observable in the training set
        real                                    miny_        = 0.0;
        //! the number of input observables (currently restricted to one)
        int                                     dimY_        = 0;
        //! the number of observable samples
        int64_t                                 nSamples_    = 0;
        //! the number of trajectory frames
        int64_t                                 nFrames_     = 0;
        //! number of frames in the training set
        int64_t                                 nTFrames_    = 0;
        //! number of frames in the validation set
        int64_t                                 nVFrames_    = 0;
        //! if set to a value > 0 triggers using the first nTrainParam_ frames for training and the rest for validation
        int64_t                                 nTrainParam_ = 0;
        //! the number of frames in the training set
        int64_t                                 nTrain_      = 0;
        //! the number of frames in the validation set
        int64_t                                 nValidate_   = 0;
        //! the number of PLS vectors to use in the FMA
        int                                     nPLSvec_     = 0;
        //! the total number of atoms in the input topology (reference structure)
        size_t                                  totalNAtoms_ = 0;
        //! the number of frames to be used in generating the pseudo-trajectory for visualizing the FMA mode
        int                                     nModeFrames_ = 0;
        //! the trajectory of input structures centered around the average structure
        MatrixClass                             fframes_;
        //! the time values of the trajectory frames
        VectorClass                             times_;
        //! indices of the input frames in the original input order
        std::vector<index_type>                 frameOrder_;
        //! the average input structure in the training set
        VectorClass                             avg_;
        //! coefficients for the calculation of the k partial least squares regressors from a linear combination of the input coordinates
        MatrixClass                             w_;
        //! the weights of the PLS factors/regressors
        VectorClass                             q_;
        //! atom selection used for least-squares superposition of the trajectory structures to the reference structure
        Selection                               fitsel_;
        //! atom selection used in the FMA
        Selection                               anasel_;
        //! coordinates of the reference structure
        RVecVector                              refstruct_;
        //! the atoms of the reference structure
        AtomsDataPtr                            refatoms_ = nullptr;
        //! the masses of the atoms in refatoms_
        VectorClass                             fitmass_;
        //! indices of the atoms used for least-squares superposition
        std::vector<int>                        fitindices_;
        //! indices of the atoms considered in the FMA
        std::vector<int>                        anaindices_;

        /*! \brief compute the arithmetic mean (average) of a set of values

            \param[in]   x       data set x
            \param[in]   n       use x values with indices < n
            \param[in]   step    use every step^th value of x

            \returns     the arithmetic  mean of the values in x
         */
        real vmean(const VectorClass &x, size_t n, size_t step) const;

        /*! \brief compute the Pearson correlation coefficient of the values in x and y

            \param[in]   x       data set x
            \param[in]   y       data set y
            \param[in]   n       use x and y values with indices < n
            \param[in]   xstep   use every xstep^th value of x
            \param[in]   ystep   use every ystep^th value of y

            \returns     the Pearson correlation coefficient between x and y
         */
        real corrcoeff(const VectorClass &x, const VectorClass &y, size_t n, size_t xstep, size_t ystep) const;

        /*! \brief helper function for writing output for k PLS components

            This function avoids code duplication for writing data of the
            full FMA mode (nPLSVec >= 1) and the ensemble-weighted FMA mode (nPLSVec = 1).

            \param[in]   nPLSVec          the number of PLS components to use for the linear model
            \param[in]   fnModel          filename for the observable values predicted by the model
            \param[in]   fnVec            filename for the vector output in the original format
            \param[in]   fnVecAnaeig      filename for the vector output in anaeig-compatible format
            \param[in]   fnMode           filename for the psdeudo-trajectory for visualizing the FMA mode
            \param[in]   writeRefStruct   write the reference structure to fnFMAref_ if true
         */
        void writeOutput(size_t             nPLSVec,
                         const std::string &fnModel,
                         const std::string &fnVec,
                         const std::string &fnVecAnaeig,
                         const std::string &fnMode,
                         bool               writeRefStruct);
};

Fma::Fma()
{
}

void
Fma::initOptions(IOptionsContainer          *options,
                 TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] performs partial least squares (PLS) based functional",
        "mode analysis (FMA). The tool requires a trajectory and a text file",
        "containing a single scalar observable value for each structure as",
        "input. It will then use PLS to identify the collective mode that is",
        "most highly correlated to the observable. This mode can be used to",
        "predict or model the observable value for any structure. The input",
        "is divided into training set and validation set. The training set",
        "will be used to learn the collective mode. The validation set is",
        "excluded from the training so that it can be used for cross-validation",
        "of the predictive value of the obtained functional mode. This",
        "predictive value is quantified by computing the Pearson-correlation",
        "coefficient between the actual input observable values and the",
        "observable value predicted by projecting the structures of the",
        "validation set on the functional mode.[PAR]",

        "The text file ([TT]-y[tt]) must contain two or three columns:",
        "First, the timestamp of the structure (which will not be used).",
        "Second, the scalar value to be modeled.",
        "Third (optional), a flag how to handle the structure-value pair",
        "(0 - ignore, 1 - use for training, 2 - use for validation.\n",
        "If a value is given for [TT]-ntrain[tt], a third column will be",
        "ignored if it is not given and no third column exists, all",
        "structure-value pairs will be used for training.[PAR]",

        "The number of PLS components to use in training the model is set by",
        "[TT]-dim[tt], or automatically optimized when using [TT]-optDim[tt].",
        "Atom selections to be used for least-squares superposition prior to FMA",
        "and for the FMA are specified with [TT]-fit[tt] and [TT]-analysis[tt],",
        "respectively.[PAR]",

        "The obtained output vectors and the separate PLS vectors can be",
        "saved with [TT]-v[tt], [TT]-ew[tt], [TT]-vec[tt] and [TT]-ewvec[tt].",
        "The observable values predicted by the model for the input trajectory",
        "frames can be written with [TT]-o[tt].[PAR]",

        "Pseudo-trajectories for visualization of the FMA modes that linearly",
        "interpolate between hypothetical average structures that correspond",
        "to the maximum and minimum observable values can be written with",
        "[TT]-om[tt] and [TT]-omew[tt].[PAR]",

        "Use [TT]-shiftAve[tt] to trigger shifting the average structure written",
        "with [TT]-v[tt], [TT]-ew[tt], [TT]-vec[tt], [TT]-ewvec[tt] along the",
        "respective output vector such that it corresponds to the minimum",
        "observable value in the training set. A projection value of zero for",
        "structures on these vectors obtained with [gmx-anaeig] or the",
        "essential dynamics module onto the FMA vectors will then correspond",
        "to this minimum instead of the average observable value in the",
        "training set.[PAR]",
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("y")
                           .filetype(eftPlot).inputFile().required()
                           .store(&fnY_).defaultBasename("y")
                           .description("Target values for the model."));
    options->addOption(FileNameOption("v")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvec_).defaultBasename("FMAvec")
                           .defaultType(efTRR).description("FMA result vectors (average structure followed by the separate PLS vectors)"));
    options->addOption(FileNameOption("ew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvecEw_).defaultBasename("FMAvec_ew")
                           .defaultType(efTRR)
                           .description("Ensemble weighted FMA result vector (average structure followed by the 1st PLS vector)"));
    options->addOption(FileNameOption("vec")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvecAnaeig_).defaultBasename("FMAvec_anaeig")
                           .defaultType(efTRR).description("FMA result vector, anaeig/make_edi compatible"));
    options->addOption(FileNameOption("vecew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAvecEwAnaeig_).defaultBasename("FMAvec_ew_anaeig")
                           .defaultType(efTRR)
                           .description("Ensemble weighted FMA result vector, anaeig/make_edi compatible"));
    options->addOption(FileNameOption("om")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAmode_).defaultBasename("FMAmode")
                           .defaultType(efPDB).description("FMA mode"));
    options->addOption(FileNameOption("omew")
                           .filetype(eftTrajectory).outputFile().required()
                           .store(&fnFMAmodeEw_).defaultBasename("FMAmode_ew")
                           .defaultType(efPDB).description("Ensemble weighted FMA mode"));
    options->addOption(FileNameOption("ref")
                           .filetype(eftPDB).outputFile().required()
                           .store(&fnFMAref_).defaultBasename("FMAref")
                           .description("Centered reference structure used for fitting."));
    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile().required()
                           .store(&fnModel_).defaultBasename("FMAmodel")
                           .description("Projection of structures from the input trajectory onto the model."));
    options->addOption(FileNameOption("oew")
                           .filetype(eftPlot).outputFile().required()
                           .store(&fnModelEw_).defaultBasename("FMAmodel_ew")
                           .description("Projection of structures from the input trajectory onto the model."));
    options->addOption(IntegerOption("dim").store(&nPLSvec_).required()
                           .defaultValue(10)
                           .description("Number of PLS-vectors to use."));
    options->addOption(IntegerOption("nframes").store(&nModeFrames_).required()
                           .defaultValue(20)
                           .description("Number of frames for FMA mode output."));
    options->addOption(Int64Option("ntrain").store(&nTrainParam_)
                           .defaultValue(-1).required()
                           .description("Use the first ntrain structures for training, all others for validation."));
    options->addOption(BooleanOption("mwfit").store(&mwfit_).required()
                           .defaultValue(false).description("Use mass-weighted fit."));
    options->addOption(BooleanOption("optDim").store(&optDim_).required()
                           .defaultValue(false).description("Optimize number of used PLS-vectors."));
    options->addOption(BooleanOption("shiftAve").store(&shiftAveStructToMinimumObs_).required()
                           .defaultValue(true).description("Shift the average structure, written to the output vector files, along the respective vector such that it corresponds to the minimum observable value."));
    options->addOption(SelectionOption("fit")
                           .store(&fitsel_).required()
                           .description("Group for least squares fitting before analysis."));
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
     * The y data of the model is saved in an .xvg text file, which is  *
     * supplied with the -y parameter. Comment lines are to be ignored. *
     * While reading, the number of data points is counted (it should   *
     * correspond with the number of frames in the trajectory) and the  *
     * dimensionality of the data is determined as the smallest number  *
     * of columns (ignoring the first columns, which represents the     *
     * time) in the file. Warnings are printed out if the number of     *
     * columns is not consistent over the length of the file and also   *
     * (at the moment) if the dimensionality is higher than one, as the *
     * current implementation can only handle one.                      *
     ********************************************************************/

    double **ydata = nullptr;

    nTrain_    = 0;
    nValidate_ =  0;

    const size_t nAtomsFitSel = fitsel_.atomCount();
    fitindices_.resize(nAtomsFitSel);

    for (size_t i = 0; i < nAtomsFitSel; i++)
    {
        fitindices_[i] = fitsel_.atomIndices()[i];
    }

    const size_t nAtomsAnaSel = fitsel_.atomCount();
    anaindices_.resize(nAtomsAnaSel);
    for (size_t i = 0; i < nAtomsAnaSel; i++)
    {
        anaindices_[i] = anasel_.atomIndices()[i];

        //refatoms_.atom[i] = top.topology()->atoms.atom[anaindices_[i]];
        //refatoms_.atomname[i] = top.topology()->atoms.atomname[anaindices_[i]];
        //refatoms_.atomtype[i] = top.topology()->atoms.atomtype[anaindices_[i]];
        //refatoms_.atomtypeB[i] = top.topology()->atoms.atomtypeB[anaindices_[i]];
    }

    refatoms_ = top.copyAtoms();

    totalNAtoms_ = refatoms_->nr;

    fitmass_.resize(totalNAtoms_);
    for (size_t i = 0; i < nAtomsFitSel; i++)
    {
        if (mwfit_)
        {
            fitmass_[fitindices_[i]] = refatoms_->atom[fitindices_[i]].m;
        }
        else
        {
            fitmass_[fitindices_[i]] = 1.0;
        }
    }

    refstruct_.resizeWithPadding(top.x().size());
    std::copy_n(top.x().begin(), top.x().size(), refstruct_.begin());
    reset_x(nAtomsFitSel, &fitindices_.front(),
            totalNAtoms_, nullptr,
            refstruct_.rvec_array(), fitmass_.data());

    nSamples_ = read_xvg(fnY_.c_str(), &ydata, &dimY_);

    y_.initArray(dimY_, nSamples_);

    dimY_--;

    if ((nSamples_ == 0) || (dimY_ == 0))
    {
        GMX_THROW(FileIOError("Error. File " + fnY_ + " did not contain any data."));
    }

    for (index_type i = 0; i < dimY_+1; i++)
    {
        for (index_type j = 0; j < static_cast<index_type>(nSamples_); j++)
        {
            y_(i, j) = ydata[i][j];
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
                    "%li structure-value pairs for training.\n\n",
                    fnY_.c_str(), nTrainParam_);
    }


    if (nTrainParam_ != -1)
    {
        nTrain_    = std::min(nTrainParam_, nSamples_);
        nValidate_ = nSamples_ - nTrain_;

        yTrain_.resize(nTrain_);
        yValidate_.resize(nValidate_);

        for (index_type i = 0; i < nTrain_; i++)
        {
            yTrain_[i] = y_(1, i);
        }
        for (index_type i = 0; i < nValidate_; i++)
        {
            yValidate_[i] = y_(1, nTrain_+i);
        }
    }
    else if (dimY_ == 1)
    {
        nTrain_    = nSamples_;
        nValidate_ = 0;

        yTrain_.resize(nTrain_);
        for (index_type i = 0; i < nTrain_; i++)
        {
            yTrain_[i] = y_(1, i);
        }
    }
    else
    {
        for (index_type i = 0; i < static_cast<index_type>(nSamples_); i++)
        {
            if (y_(2, i) == 1)
            {
                nTrain_++;
            }
            if (y_(2, i) == 2)
            {
                nValidate_++;
            }
        }

        yTrain_.resize(nTrain_);
        yValidate_.resize(nValidate_);

        size_t j = 0;
        for (index_type i = 0; i < static_cast<index_type>(nSamples_); i++)
        {
            if (y_(2, i) == 1)
            {
                yTrain_[j] = y_(1, i);
                j++;
            }
        }

        size_t k = 0;
        for (index_type i = 0; i < static_cast<index_type>(nSamples_); i++)
        {
            if (y_(2, i) == 2)
            {
                yValidate_[k] = y_(1, i);
                k++;
            }
        }

    }

    if (optDim_ && (nValidate_ < 2))
    {
        GMX_THROW(gmx::InconsistentInputError("PLS component optimization needs at least 2 (better more) datapoints assigned to validation."));
    }

    printf("Read %li lines of %i column(s).\n", nSamples_, dimY_);
    printf("Using %li for training and %li for validation.\n", nTrain_, nValidate_);

    if (nTrain_ == 0)
    {
        GMX_THROW(gmx::FileIOError("Error. File " + fnY_ + " did not contain any training-data."));
    }

    avgy_ = 0;
    miny_ = *(std::min_element(std::begin(yTrain_), std::end(yTrain_)));
    for (index_type i = 0; i < nTrain_; i++)
    {
        avgy_ += yTrain_[i];
    }

    avgy_ = avgy_  / static_cast<real>(nTrain_);

    printf("\nThe average value of y (in the training set) is %.3f .\n\n", avgy_);

    for (index_type i = 0; i < nTrain_; i++)
    {
        yTrain_[i] -= avgy_;
    }
    for (index_type i = 0; i < nValidate_; i++)
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

    fframes_.initArray((nTrain_+nValidate_), (DIM * anasel_.atomCount()));

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

    if (ydata != nullptr)
    {
        for (index_type i = 0; i < dimY_+1; ++i)
        {
            sfree(ydata[i]);
        }
        sfree(ydata);
    }
}


void
Fma::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /*pbc*/,
                  TrajectoryAnalysisModuleData * /*pdata*/)
{
    /**********************************************************************
     * reading in the trajectory                                          *
     *                                                                    *
     * For PLS-FMA, the coordinates of the whole trajectory have to be    *
     * saved in memory (at least for the atoms that are used in the       *
     * analysis). This is accomplished here, and the real analysis is     *
     * then performed once the whole trajectory has been read.            *
     **********************************************************************/

    size_type frameNumber   = 0;
    bool      keepFrame     = true;
    bool      useForAverage = true;

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
            useForAverage = false;
            nVFrames_++;
        }
    }
    else
    {
        if (y_(2, nFrames_) == 1)
        {
            frameNumber = nTFrames_;
            nTFrames_++;
        }
        else if (y_(2, nFrames_) == 2)
        {
            frameNumber   = nTrain_ + nVFrames_;
            useForAverage = false;
            nVFrames_++;
        }
        else
        {
            keepFrame = false;
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
        reset_x(fitsel_.atomCount(), &fitindices_.front(),
                totalNAtoms_, nullptr,
                fr.x, fitmass_.data());

        do_fit(totalNAtoms_, fitmass_.data(), refstruct_.rvec_array(), fr.x);

        for (index_type i = 0; i < anasel_.atomCount(); i++)
        {
            for (size_t j = 0; j < DIM; j++)
            {
                /* already saving the coordinates in FORTRAN compatible COLUMN MAJOR order
                   each ROW describes one structure, while each COLUMN describes one coordinate */
                fframes_(frameNumber, (DIM*i) + j) = fr.x[anaindices_[i]][j];

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
    nFrames_  = nTFrames_ + nVFrames_;
    nSamples_ = nTrain_ + nValidate_;

    if (nFrames_ < nSamples_)
    {
        GMX_THROW(gmx::InconsistentInputError("Trajectory containts fewer frames than there are datapoints in File " + fnY_));
    }


    /* center data */
    for (index_type k = 0; k < nSamples_; k++)
    {
        for (index_type i = 0; i < anasel_.atomCount(); i++)
        {
            for (index_type j = 0; j < DIM; j++)
            {
                fframes_(k, (DIM*i) + j) -= avg_[(DIM*i) + j];
            }
        }
    }

    printf("Running PLS on %li training structures (of %li total structures)\n", nTrain_, nSamples_);

    real bestCC     = -1.0;
    int  optimalDim = 0;
    PartialLeastSquares<MatrixClass, VectorClass> pls;

    if (optDim_)
    {
        printf("\nDetermining optimal number of PLS components.\n");
        printf("  n r(training) r(validation)\n");

        for (index_type l = 1; l <= nPLSvec_; l++)
        {
            Fma::MatrixClass w_l((anasel_.atomCount() * DIM), l);
            q_.resize(l);

            VectorClass      PLSvec(anasel_.atomCount()*DIM);
            VectorClass      mYt(nTrain_);
            VectorClass      mYv(nValidate_);

            Fma::MatrixClass X(nTrain_, (DIM * anasel_.atomCount()));
            for (size_type idx_row = 0; idx_row < X.length1(); ++idx_row)
            {
                for (size_type idx_col = 0; idx_col < X.length2(); ++idx_col)
                {
                    X(idx_row, idx_col) = fframes_(idx_row, idx_col);
                }
            }

            // perform PLS regression where the numbers of analyzed variables,
            // structure-observable pair samples, and regressors is deduced from
            // the dimensions of the input arrays
            pls.pls_denham(X, yTrain_, &w_l, &q_);
            for (index_type i = 0; i < anasel_.atomCount(); i++)
            {
                for (index_type j = 0; j < DIM; j++)
                {
                    for (index_type k = 0; k < l; k++)
                    {
                        PLSvec[i*DIM +j] += q_[k] * w_l((i*DIM + j), k);
                    }
                }
            }

            for (index_type i = 0; i < nSamples_; i++)
            {
                if (frameOrder_[i] < nTrain_)
                {
                    mYt[frameOrder_[i]] = avgy_;
                    for (index_type j = 0; j < anasel_.atomCount()*DIM; j++)
                    {
                        mYt[frameOrder_[i]] += PLSvec[j] * (fframes_(frameOrder_[i], j));
                    }
                }
                else
                {
                    mYv[frameOrder_[i]-nTrain_] = avgy_;
                    for (index_type j = 0; j < anasel_.atomCount()*DIM; j++)
                    {
                        mYv[frameOrder_[i]-nTrain_] += PLSvec[j] * (fframes_(frameOrder_[i], j));
                    }
                }

            }

            real cct = corrcoeff(yTrain_, mYt, nTrain_, 1, 1);
            real ccv = corrcoeff(yValidate_, mYv, nValidate_, 1, 1);
            printf("%li %9.3f %13.3f\n", l, cct, ccv);
            if (ccv > bestCC)
            {
                bestCC     = ccv;
                optimalDim = l;
            }
        }
        printf("\nUsing %d PLS components for best correlation on validation set (r=%8.3f)\n", optimalDim, bestCC);
        nPLSvec_ = optimalDim;
    }

    w_.initArray((anasel_.atomCount() * DIM), nPLSvec_);
    q_.resize(nPLSvec_);

    Fma::MatrixClass X(nTrain_, (DIM * anasel_.atomCount()));
    for (size_type idx_row = 0; idx_row < X.length1(); ++idx_row)
    {
        for (size_type idx_col = 0; idx_col < X.length2(); ++idx_col)
        {
            X(idx_row, idx_col) = fframes_(idx_row, idx_col);
        }
    }

    // the number of input coordinates, samples and regressors is auto-determined
    // from the array/vector sizes
    pls.pls_denham(X, yTrain_, &w_, &q_);

    please_cite(stdout, "HubdeGroot2009");
    please_cite(stdout, "KrivobokovaEtAl2012");
}

void
Fma::writeOutput()
{
    // write output for the ensemble-weighted FMA mode (use 1st PLS component only)
    writeOutput(1,        fnModelEw_, fnFMAvecEw_, fnFMAvecEwAnaeig_, fnFMAmodeEw_, false);
    // write output for the full FMA mode (use the preset or optimized nPLSvec_ PLS components)
    writeOutput(nPLSvec_, fnModel_,   fnFMAvec_, fnFMAvecAnaeig_, fnFMAmode_, true);
}

void Fma::writeOutput(size_t             nPLSvec,
                      const std::string &fnModel,
                      const std::string &fnVec,
                      const std::string &fnVecAnaeig,
                      const std::string &fnMode,
                      bool               writeRefStruct)
{
    real minProjection   = std::numeric_limits<real>::max();
    real maxProjection   = std::numeric_limits<real>::min();
    real v2              = 0.0;

    // the lambda value saved in the vector .trr-files determines whether
    // g_anaeig uses mass weighting in the structure superposition or not
    const real lambda = mwfit_ ? 1.0 : 0.0;

    // the PLS vector using k PLS components
    VectorClass                     PLSvec(anasel_.atomCount()*DIM, 0);

    // the observable values predicted by the PLS model for the training set
    VectorClass                     mYt(nTrain_);
    // the observable values predicted by the PLS model for the validation set
    VectorClass                     mYv(nValidate_);

    // element 0 contains the times
    // element 1 the observable values predicted by the PLS model for the full data set
    std::vector<VectorClass>        modelY(2, VectorClass(nSamples_));

    // make a copy of the full structure which we can manipulate while creating the
    // pseudo-trajectory for visualization of the FMA mode
    RVecVector fullStruct(totalNAtoms_, {0, 0, 0});
    for (size_type i = 0; i < totalNAtoms_; i++)
    {
        for (size_type j = 0; j < DIM; j++)
        {
            fullStruct[i][j] = refstruct_[i][j];
        }
    }
    // structure of the atom subset analyzed, used as temporary storage for output of
    // the reference structure, the average structure and the vector components
    RVecVector structure(anasel_.atomCount(), {0, 0, 0});

    t_trxframe fr;
    clear_trxframe(&fr, TRUE);
    fr.bAtoms = TRUE;
    fr.atoms  = refatoms_.get();

    matrix box;
    for (index_type i = 0; i < 3; i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            box[i][j]    = 0;
            fr.box[i][j] = 0;
        }
    }

    // the additional flag avoids writing the reference structure twice
    if (!fnFMAref_.empty() && writeRefStruct)
    {
        write_sto_conf(fnFMAref_.c_str(), "Reference structure used in FMA",
                       refatoms_.get(), refstruct_.rvec_array(), nullptr, -1, box);
    }

    // output vectors describing the functional modes in anaeig, make_edi compatible standard format
    // step  -1 reference structure for superpositiion (t = -1),
    // step   0 average structure as reference for computing projections (t = 0),
    // step   1 the vector itself (t = squared vector length (eigenvalue in the PCA case))
    t_fileio    *outputfile_anaeig   = gmx_trr_open(fnVecAnaeig.c_str(), "w");
    // output vectors describing the functional modes in the generic format used previously
    t_fileio    *outputfile          = gmx_trr_open(fnVec.c_str(), "w");
    // the pseudo-trajectories for visualizing the functional modes
    t_trxstatus *modefile            = open_trx(fnMode.c_str(), "w");;

    // compute the PLS vector and its squared length
    for (index_type i = 0; i < anasel_.atomCount(); i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            const size_t tmpIndex = i*DIM + j;
            for (size_type k = 0; k < nPLSvec; k++)
            {
                PLSvec[tmpIndex] += q_[k] * w_(tmpIndex, k);
            }
            v2 += PLSvec[tmpIndex] * PLSvec[tmpIndex];
        }
    }

    const real one_over_v2 = (v2 > 0 ? 1.0 / v2 : std::numeric_limits<real>::max());

    // write the reference structure for superposition to file
    for (index_type i = 0; i < anasel_.atomCount(); i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            structure.rvec_array()[i][j] = refstruct_.rvec_array()[anasel_.atomIndices()[i]][j];
        }
    }
    gmx_trr_write_frame(outputfile_anaeig, -1, -1.0, lambda,
                        box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);


    // write the average structure to file which is really the reference structure
    // for computing the displacement of the particle coordinates that is then
    // projected onto the vector(s) via y - yave = (x - xave) * v
    for (index_type i = 0; i < anasel_.atomCount(); i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            structure.rvec_array()[i][j] = avg_[i*DIM + j];
            if (shiftAveStructToMinimumObs_)
            {
                // Shift the average structure along the PLS vector such that it
                // corresponds to the minimum observable value y_min in the training set.
                // y^0 - y^ave = (x^0 - x^ave) v, with y^0 = 0 and x^0 to be determined
                // the following Ansatz fullfils the equation:
                // x^0 = x^ave + v / |v|^2 * (y^0 - y^ave).
                structure.rvec_array()[i][j] += PLSvec[i*DIM +j] * (miny_ - avgy_) * one_over_v2;
            }
        }
    }

    gmx_trr_write_frame(outputfile, 0, 0.0, lambda,
                        box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);
    gmx_trr_write_frame(outputfile_anaeig, 0, 0.0, lambda,
                        box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);

    // now the full FMA vector
    for (index_type i = 0; i < anasel_.atomCount(); i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            structure.rvec_array()[i][j] = PLSvec[i*DIM +j];
        }
    }

    gmx_trr_write_frame(outputfile, 1, 1.0, lambda,
                        box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);

    // write the separate PLS vectors and their weights to the non-standard output file
    for (size_type k = 0; k < nPLSvec; k++)
    {
        for (index_type i = 0; i < anasel_.atomCount(); i++)
        {
            for (index_type j = 0; j < DIM; j++)
            {
                structure.rvec_array()[i][j] = q_[k] * w_((i*DIM + j), k);
            }
        }
        gmx_trr_write_frame(outputfile, k+1, q_[k], lambda,
                            box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);
    }

    // write the normalized FMA vector to the anaeig/make_edi compatible file
    const real          v          = std::sqrt(v2);
    const real          one_over_v = (v > 0 ? 1 / v : std::numeric_limits<real>::max());
    for (index_type i = 0; i < anasel_.atomCount(); i++)
    {
        for (index_type j = 0; j < DIM; j++)
        {
            structure.rvec_array()[i][j] = PLSvec[i*DIM +j] * one_over_v;
        }
    }
    gmx_trr_write_frame(outputfile_anaeig, 1, v2, lambda,
                        box, anasel_.atomCount(), structure.rvec_array(), nullptr, nullptr);

    // done writing the vector files
    gmx_trr_close(outputfile);
    gmx_trr_close(outputfile_anaeig);

    for (index_type i = 0; i < nSamples_; i++)
    {
        modelY[0][i] = times_[frameOrder_[i]];
        modelY[1][i] = avgy_;
        for (index_type j = 0; j < anasel_.atomCount()*DIM; j++)
        {
            modelY[1][i] += PLSvec[j]   * (fframes_(frameOrder_[i], j));
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

    printf("minimal projection using %zu PLS component%s = %8.3f (same unit as -y)\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), minProjection);
    printf("maximal projection using %zu PLS component%s = %8.3f (same unit as -y)\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), maxProjection);

    // put the origin of the Cartesian coordinate to the minimum or average observable value in the training set
    const real referenceObs = (shiftAveStructToMinimumObs_ ? miny_ : avgy_);
    // shift the projections to the new reference observable value and convert to length in [nm]
    // this is the same range of values that should be obtained when projecting
    // the trajectory on the output vector with anaeig
    const real minProjectionMetric = (minProjection - referenceObs) * one_over_v;
    const real maxProjectionMetric = (maxProjection - referenceObs) * one_over_v;

    printf("minimal projection using %zu PLS component%s = %8.3f nm\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), minProjectionMetric);
    printf("maximal projection using %zu PLS component%s = %8.3f nm\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), maxProjectionMetric);

    // Create a trajectory to visualize the FMA mode by linear interpolation along the FMA vector
    fr.natoms = totalNAtoms_;
    fr.lambda = lambda;

    for (index_type k = 0; k < nModeFrames_; k++)
    {
        // need this workaround because RVec operator[] has two ambigous overloads
        rvec * const fsPtr = fullStruct.rvec_array();
        for (index_type i = 0; i < anasel_.atomCount(); i++)
        {
            for (index_type j = 0; j < DIM; j++)
            {
                fsPtr[anasel_.atomIndices()[i]][j] = avg_[i*DIM + j] + (minProjection + (maxProjection-minProjection)/nModeFrames_*k - avgy_) / v2 * PLSvec[i*DIM +j];
            }
        }

        fr.time = k+1;
        fr.step = k+1;
        fr.bX   = true;
        fr.x    = fullStruct.rvec_array();

        write_trxframe_indexed(modefile, &fr, anasel_.atomCount(),
                               &anaindices_.front(), nullptr);
    }

    close_trx(modefile);

    real cct = corrcoeff(yTrain_, mYt, nTrain_, 1, 1);
    printf("r(training) using %zu PLS component%s   = %8.3f\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), cct);
    if (nValidate_ > 0)
    {
        real ccv = corrcoeff(yValidate_, mYv, nValidate_, 1, 1);
        printf("r(validation) using %zu PLS component%s = %8.3f\n", nPLSvec, (nPLSvec > 1 ? "s" : ""), ccv);
    }

    gmx_output_env_t* oenv;
    output_env_init_default(&oenv);
    FILE             *modelout;
    modelout = xvgropen(fnModel.c_str(), "Projection of structure onto model.", "X", "Y", oenv);
    for (index_type i = 0; i < nSamples_; i++)
    {
        fprintf(modelout, " %12.5e %12.5e\n", modelY[0][i], modelY[1][i]);
    }
    xvgrclose(modelout);
    output_env_done(oenv);
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
real Fma::vmean(const VectorClass &x, size_t n, size_t step) const
{
    real m = 0.0;

    for (size_t i = 0; i < n; i++)
    {
        m += x[i*step];
    }
    return m / static_cast<real>(n);
}

/*
 * Calculation of the correlation coefficient between two vectors
 *
 * Simple calculation of the Pearson product-moment correlation coefficient between two
 * (n)-vectors x and y
 *
 * \param x first (n*step)-vector
 * \param y second (n*step)-vector
 * \param number of values to use in calculation
 * \param xstep use only every step-th element of x
 * \param ystep use only every step-th element of y
 */
real Fma::corrcoeff(const VectorClass &x, const VectorClass &y, size_t n, size_t xstep, size_t ystep) const
{
    real       cxy = 0;
    real       cxx = 0;
    real       cyy = 0;

    const real mx = vmean(x, n, xstep);
    const real my = vmean(y, n, ystep);

    for (size_t i = 0; i < n; i++)
    {
        cxy += (x[i*xstep]-mx) * (y[i*ystep]-my);
        cxx += (x[i*xstep]-mx) * (x[i*xstep]-mx);
        cyy += (y[i*ystep]-my) * (y[i*ystep]-my);
    }

    return cxy / (std::sqrt(cxx) * std::sqrt(cyy));
}

} /* namespace */

const char FMAInfo::name[]             = "fma";
const char FMAInfo::shortDescription[] =
    "Functional Mode Analysis";

TrajectoryAnalysisModulePointer FMAInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Fma);
}

} /* namespace analysismodules */

} /* namespace gmx */
