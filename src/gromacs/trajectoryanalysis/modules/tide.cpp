/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * Implements gmx::analysismodules::tide.
 *
 * \author Camilo Aponte <ca.aponte.uniandes.edu.co>
 * \author Rodolfo Briones <fitobriones@gmail.com>
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "tide.h"

#include <map>
#include <string.h>

#include "gromacs/math/container/mask.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/encompassinggrid.h"
#include "gromacs/math/griddata/operations/densityspreader.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/trajectoryanalysis.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strdb.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/* enumerate physical models for spreading density */
enum DensityType {
    enXRay, enCryoEM, enVdW, enCharge
};
const char *const cDensityTypeEnum[] = {"x-ray", "cryo-EM", "vdW", "charge"};

/*!\brief Class to control command line printing.
 * Control verobsity.
 */
class CommandLineOutput
{
    public:
        void operator()(std::string output)
        {
            if (verbose_)
            {
                fprintf(stderr, "\n%s\n", output.c_str());
            }
        };
        void warning(std::string warningText)
        {
            if (verbose_)
            {
                fprintf(stderr, "\nWARNING: %s\n", warningText.c_str());
            }
        }
        /*!\brief Results are always printed, no matter the verbosity.
         */
        void result(std::string resultText)
        {
            fprintf(stdout, "\n\t%s\n", resultText.c_str());
        }
        void setVerbosity(bool verbose)
        {
            verbose_ = verbose;
        }
    private:
        bool verbose_ = false;
};


struct GaussianParameters{
    // less-than operator enables storage in std::map
    bool operator<(const GaussianParameters &other) const
    {
        return (other.amplitude < amplitude) || ((other.amplitude == amplitude) && (other.precision < precision));
    }
    float amplitude;
    float precision;
};

//! \brief TODO
float maskedNormalizedInnerProduct(const GridDataReal3D::container_type &vec1,
                                   const GridDataReal3D::container_type &vec2,
                                   GridDataReal3D::value_type            mean1,
                                   GridDataReal3D::value_type            mean2,
                                   const Mask                           &mask)
{
    double     innerProduct = 0.;
    const auto minMapSize   = std::min(vec1.size(), vec2.size());
#pragma omp parallel shared(innerProduct)
    {

#pragma omp for reduction(+:innerProduct)
        for (size_t l = 0; l < minMapSize; l++)
        {
            if (mask[l])
            {
                innerProduct += (vec1[l]-mean1) * (vec2[l]-mean2);
            }
        }
    }

    return innerProduct;
}

GridDataReal3D::value_type maskedMean(const GridDataReal3D::container_type &v, const Mask &mask)
{
    double sum = 0.;
#pragma omp parallel shared(sum)
    {
#pragma omp for reduction(+:sum)
        for (size_t l = 0; l < v.size(); ++l)
        {
            if (mask[l])
            {
                sum += v[l];
            }
        }
    }
    return sum/v.size();
}

//! \brief TODO
float correlationCoefficient(const GridDataReal3D::container_type &map1,
                             const GridDataReal3D::container_type &map2,
                             const Mask                           &mask)
{
    const auto mean1 = maskedMean(map1, mask);
    const auto mean2 = maskedMean(map2, mask);
    return maskedNormalizedInnerProduct(map1, map2, mean1, mean2, mask) / sqrt(maskedNormalizedInnerProduct(map2, map2, mean2, mean2, mask) *
                                                                               maskedNormalizedInnerProduct(map1, map1, mean1, mean1, mask));
}

/*! \brief
 * Time averaged density analysis.
 */
class TimeAveragedDensity final : public TrajectoryAnalysisModule
{
    public:
        TimeAveragedDensity();

        void initOptions(IOptionsContainer          *options,
                         TrajectoryAnalysisSettings *settings) override;
        void optionsFinished(TrajectoryAnalysisSettings *settings) override;
        void initAnalysis(const TrajectoryAnalysisSettings &settings,
                          const TopologyInformation        &top) override;

        void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                          TrajectoryAnalysisModuleData *pdata) override;

        void finishAnalysis(int nframes) override;
        void writeOutput() override;

    private:
        std::map < int, std::vector < GaussianParameters>> setGaussianCoeff(bool hasSigma);
        void setAttypeStringsFromTop(const t_topology &top);

        float                    maskCutoff_                   = 1e-5; //< rho < mask_rheshold will be ignored from the calculation
        int                      blocksize_                    = 1;    //< for time-resolved global correlation
        float                    sigma_                        = 0.3;  //< spread with in nm
        bool                     sigmaIsSet_                   = false;
        bool                     spacingIsSet_                 = false;
        float                    spacing_                      = 0.2; // grid spacing in nm
        float                    margin_                       = 0.5; // grid margin around atoms in nm
        float                    nThreadsOmp_                  = 1;   // number of openMP threads in analysis
        float                    nSigma_                       = 5.0; // width of spreading in sigma

        SelectionList            sel_;
        std::vector<real>        atomCharges_;
        std::vector<real>        vdWRadii_;

        AnalysisData             goodnessOfFitData_;

        std::string              fnGaussParameters_;      // input structure factor file
        std::string              fnMapAverage_;           // average map (3d)
        std::string              fnMapStandardDeviation_; // standard deviation (3d)
        std::string              fnReferenceMap_;
        std::string              fnMaskMap_;
        std::string              fnCorrelation_;

        bool                     verbose_          = false;
        bool                     integrateDensity_ = false;
        std::map < int, std::vector < GaussianParameters>> gaussianParametersPerAtomicNumber_; // Ai factors (Ngaussians x Natoms 1D array)
        std::vector<int>         atomicNumbersIndex_;                                          // array will store the index of the

        // === define the lattice that will store the density for each point===
        GridDataReal3D    referenceMap_;                      /* The reference=experimental map to compare to        */
        GridDataReal3D    averageMap_;                        /* The computed timed-averaged map*/
        GridDataReal3D    stddevMap_;                         /* The computed timed-averaged standard deviation map */
        GridDataReal3D    blockAverageMap_;                   /* map for block averaging*/
        Mask              mask_;                              //< Masking out density values
        CommandLineOutput output_;                            //< output to command line

        float             nFramesBlock_             = 0;      // number of frames for block averaging
        int               dataFrameNumber_          = 0;      // enumerates output data frames, increases in frames where data is output
        DensityType       densityType_              = enXRay; //< Physical model for spreading the atoms
};

//! \brief Set the Gaussian coefficients for cryo-EM scattering,
std::map < int, std::vector < GaussianParameters>> TimeAveragedDensity::setGaussianCoeff(bool hasSigma)
{
    // read gaussian parameters from file if provided, else use defaults
    FILE * parameterFile = nullptr;
    if (!fnGaussParameters_.empty())
    {
        output_("Reading Gaussian coefficients from " + fnGaussParameters_ );
        parameterFile = gmx_ffopen(fnGaussParameters_.c_str(), "r");
    }
    else
    {
        switch (densityType_)
        {
            case enXRay:
                output_("Reading default Gaussian coefficients from 'electronscattering.dat'");
                parameterFile = libopen("electronscattering.dat");
                break;
            case enCryoEM:
                output_("Reading electron microscopy Gauss weights from 'emcrosssections.dat");
                parameterFile = libopen("emcrosssections.dat");
                break;
            default:
                break;
        }
    }
    if (parameterFile == nullptr && !hasSigma && densityType_ != enVdW)
    {
        GMX_THROW(FileIOError("Cannot open file " + fnGaussParameters_ + " or library file for reading and no spreading option -sigma was given. Exiting now because cannot reasonably guess how to spread atoms density."));
    }
    std::map < int, std::vector < GaussianParameters>> gaussianParametersLookup;
    constexpr size_t maxLineLength = 1000;
    float            precision     = 1/(2.*sigma_*sigma_);
    if (parameterFile != nullptr)
    {
        char line[maxLineLength];
        while (get_a_line(parameterFile, line, maxLineLength))
        {
            int   currentAtomType;
            float amplitude;
            if (hasSigma)
            {
                if (sscanf(line, " %d %e", &currentAtomType, &amplitude) == 2)
                {
                    gaussianParametersLookup[currentAtomType].push_back({amplitude, precision});
                }
            }
            else
            {
                if (sscanf(line, " %d %e %e", &currentAtomType, &amplitude, &precision) == 3)
                {
                    gaussianParametersLookup[currentAtomType].push_back({amplitude, precision});
                }
            }
        }
        gmx_ffclose(parameterFile);
    }
    else
    {
        constexpr size_t maxElementType = 100;
        for (size_t currentAtomType = 0; currentAtomType < maxElementType; ++currentAtomType)
        {
            gaussianParametersLookup[currentAtomType].push_back({1, precision});
        }
    }
    return gaussianParametersLookup;
}

TimeAveragedDensity::TimeAveragedDensity()
{
    registerAnalysisDataset(&goodnessOfFitData_, "fitness");
}

void TimeAveragedDensity::initOptions(IOptionsContainer          *options,
                                      TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes density and its properties of a group of "
        "spread atoms from a trajectory or structure. The density is calculated as a linear combination of Gaussians",
        "rho(R)=sum_{i=1}^{i=Ngaussians} Ai*exp(-Bi*R^2)[PAR]",

        "By default, 2-Gaussian coefficients from Electron Microscopy are considered for H,C,N,O,P,S.",
        "Non-default values can be specified via the -gaussparameters input file as",
        "\"atomic-number A B\". Precomputed A and B parameters, for 2- and 4-Gaussian sums, from X-ray "
        "and EM, as well as for coarse-grained systems, are available at XXXX.",
        "From this, the following output may be generated:",
        "",
        "* [TT]-mo[tt]: time average of density",
        "* [TT]-stddevmap[tt]: time standard deviation of density",
        "",
        "A reference map ([TT]-refmap[tt]) enables the following comparison options",
        "",
        "* [TT]-corr[tt]: Time resolved correlation between computed and the reference map",
        "* [TT]-blocksize[tt] The correlation is computed in time windows",
        "A masking map can be provided with the -mask option. Grid points in the ",
        "mask, with density smaller than a user-specified",
        "threshold value (-maskcutoff), will be ignored. In the output maps, ",
        "the masked points will be labeled as NAN[PAR]",
        "    NOTE: If the density is intended to be calculated with respect to ",
        "other reference group, a previous fitting step (e.g. with trjconv) ",
        "should precede the gtide calculation.[PAR]",
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("gaussparameters").filetype(eftGenericData)
                           .inputFile().store(&fnGaussParameters_).libraryFile()
                           .defaultBasename("electronscattering")
                           .description("Gaussian coefficients"));

    options->addOption(FileNameOption("refmap").filetype(eftCCP4).inputFile()
                           .store(&fnReferenceMap_).defaultBasename("refmap")
                           .description("Reference map"));

    options->addOption(FileNameOption("mask").filetype(eftCCP4).inputFile()
                           .store(&fnMaskMap_).defaultBasename("mask")
                           .description("Map for masking above threshold"));

    options->addOption(FileNameOption("mo").filetype(eftCCP4).outputFile()
                           .store(&fnMapAverage_).defaultBasename("average")
                           .description("Average density map"));

    options->addOption(FileNameOption("stddevmap").filetype(eftCCP4).outputFile()
                           .store(&fnMapStandardDeviation_).defaultBasename("stddev")
                           .description("Standard deviation map"));

    options->addOption(FileNameOption("correlation").filetype(eftPlot).outputFile()
                           .store(&fnCorrelation_).defaultBasename("correlation")
                           .description("Correlation between computed and reference map"));

    options->addOption(SelectionOption("select").storeVector(&sel_).required()
                           .description("Selection to calculate the density"));

    options->addOption(FloatOption("maskcutoff").store(&maskCutoff_)
                           .description("Voxels smaller than maskcutoff in mask map are excluded from calculations"));

    options->addOption(BooleanOption("integrate").store(&integrateDensity_)
                           .description("Integrate density over voxel Volumes instead of Gauss transform"));

    options->addOption(FloatOption("spacing").store(&spacing_).storeIsSet(&spacingIsSet_)
                           .description("Map spacing in nm"));

    options->addOption(IntegerOption("blocksize").store(&blocksize_)
                           .description("Number of frames for block-averaging correlation coefficient"));
    options->addOption(FloatOption("sigma").store(&sigma_).storeIsSet(&sigmaIsSet_).
                           description("Spread width for atoms in nm"));
    options->addOption(FloatOption("nsigma").store(&nSigma_).
                           description("Range in sigma until spreading is truncated"));
    options->addOption(BooleanOption("v").store(&verbose_).description("Verboseness"));

    options->addOption(FloatOption("ntomp").store(&nThreadsOmp_).description("Number of OpenMP threads for analysis"));

    options->addOption(EnumOption<DensityType>("type").enumValue(cDensityTypeEnum)
                           .store(&densityType_)
                           .description("Physical model for spreading"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}


void TimeAveragedDensity::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    // Set gaussian parameters per atomic numbers from external files if
    // a physical model that weights density based on atom species is chosen
    gaussianParametersPerAtomicNumber_ = setGaussianCoeff(sigmaIsSet_);
}

void TimeAveragedDensity::initAnalysis(
        const gmx_unused TrajectoryAnalysisSettings &settings,
        const TopologyInformation                   &top)
{
    output_.setVerbosity(verbose_);

    goodnessOfFitData_.setColumnCount(0, sel_.size());

    if (!spacingIsSet_ && sigmaIsSet_)
    {
        spacing_ = sigma_ / 3.0;
    }

    if (densityType_ == enCharge && !top.hasFullTopology())
    {
        GMX_THROW(InconsistentInputError("Cannot use ' -type charge ' without full topology information (i.e. a tpr file), because this is the only way to determine the partial charge per atom."));
    }
    if (top.hasFullTopology())
    {
        ArrayRef<const t_atom> atoms({top.topology()->atoms.atom, top.topology()->atoms.atom+top.topology()->atoms.nr});
        for (const auto &atom : atoms)
        {
            atomicNumbersIndex_.push_back(atom.atomnumber);
            atomCharges_.push_back(atom.q);
        }
    }
    else
    {
        if (top.hasTopology())
        {
            output_("WARNING: Guessing atomic number from atom names");
            const auto &atoms = top.topology()->atoms;
            auto        aps   = gmx_atomprop_init();
            for (int i = 0; i < atoms.nr; ++i)
            {
                real atomicNumber;
                gmx_atomprop_query(aps, epropElement, "???", *atoms.atomname[i], &atomicNumber);
                if (atomicNumber < 0)
                {
                    output_.warning("Cannot guess atomic number for '" + std::string(*atoms.atomname[i]) + "'  - ignoring for spreading. Is this a vitual site?");
                }
                if (densityType_ == enVdW)
                {
                    real vdWRadius;
                    gmx_atomprop_query(aps, epropVDW, "???", *atoms.atomname[i], &vdWRadius);
                    if (vdWRadius > 0)
                    {
                        vdWRadii_.push_back(vdWRadius);
                    }
                    else
                    {
                        output_.warning("Cannot find vdW radius for this atom type, setting radius to zero");
                        vdWRadii_.push_back(0);
                    }
                }
                atomicNumbersIndex_.push_back(static_cast<int>(atomicNumber));
            }
            gmx_atomprop_destroy(aps);
        }
        else
        {
            output_.warning("Cannot find topology information, setting all atom weights to unity");
        }
    }
    std::unique_ptr < IGrid < DIM>> referenceGrid;
    // create grid from reference map if provided, alternatively from structure
    if (!fnReferenceMap_.empty())
    {
        output_("Reading reference map "+ fnReferenceMap_ );
        referenceMap_       = MrcFile().read(fnReferenceMap_);
        output_("Using grid from reference map as reference");
        referenceGrid       = referenceMap_.getGrid().duplicate();
    }
    else
    {
        rvec      * coordinates;
        matrix      box;
        top.getTopologyConf(&coordinates, box);
        output_("Using reference topologly to construct density grid");
        referenceGrid = encompassingGridFromCoordinates({coordinates, coordinates+top.topology()->atoms.nr}, spacing_, margin_).duplicate();
    }

    averageMap_.setGrid(referenceGrid->duplicate());
    stddevMap_.setGrid(referenceGrid->duplicate());

    /* read mask map if provided */
    if (!fnMaskMap_.empty())
    {
        mask_ = Mask(MrcFile().read(fnMaskMap_), maskCutoff_);
    }

    if (!fnCorrelation_.empty())
    {
        blockAverageMap_.setGrid(referenceGrid->duplicate());
        if (fnReferenceMap_.empty())
        {
            GMX_THROW(FileIOError("Correlation calculation requires a reference map. \n"
                                  "If you intended to calculate the density correlation between a reference structure"
                                  " and a trajectory, first generate a reference density from the reference structure ."));
        }
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnCorrelation_);
        plotm->setTitle("Time-resolved correlation between computed and reference maps");
        plotm->setSubtitle("Number of frames for block-averaging:" + std::to_string(blocksize_)+" frames" );
        plotm->setXAxisIsTime();
        plotm->setYLabel("Correlation");
        goodnessOfFitData_.addModule(plotm);
    }
}

void TimeAveragedDensity::analyzeFrame(int /*frnr*/, const t_trxframe & /*fr*/, t_pbc * /*pbc*/, TrajectoryAnalysisModuleData *pdata)
{
    const size_t       group = 0; // selection-group index (only one considered thus it is 0)
    const Selection   &sel   = pdata->parallelSelection(sel_[group]);

    // Gather all atoms with the same Gaussian Parameters, also considering multiple Gaussian Parameters per atom
    std::map < float, std::vector < RVec>> atomCoordinatesPerPrecision;
    std::map < float, std::vector < real>> atomWeightsPerPrecision;

    for (int atomIndex = 0; atomIndex < sel.atomCount(); atomIndex++)
    {
        // try to look up the gaussian parameters corresponding to the current atoms atomic number
        int atomicNumber = 0;
        // if atomic numbers are not available, the atomic number stays zero
        try
        {
            atomicNumber = atomicNumbersIndex_.at(sel.atomIndices()[atomIndex]);
        }
        catch (const std::out_of_range &oor)
        {
            output_.warning("Cannot determine atomic number of atom #" + std::to_string(sel.atomIndices()[atomIndex]) +  ". Falling back to atomic number 0.");
        }
        try
        {
            const auto &gaussianParameter = gaussianParametersPerAtomicNumber_.at(atomicNumber);
            for (const auto &gp : gaussianParameter)
            {
                atomCoordinatesPerPrecision[gp.precision].push_back(sel.position(atomIndex).x());
                real weight;
                switch (densityType_)
                {
                    case enCharge:
                        weight = atomCharges_[sel.atomIndices()[atomIndex]];
                        break;
                    case enVdW:
                        weight = 1;
                        break;
                    default:
                        weight = gp.amplitude;
                        break;
                }
                atomWeightsPerPrecision[gp.precision].push_back(weight);
            }
        }
        catch  (const std::out_of_range &oor)
        {
            output_.warning("Cannot find atomic number " + std::to_string(atomicNumber) + " in spread parameters.");
            output_.warning("Ignoring atom number " + std::to_string(sel.atomIndices()[atomIndex]) + " from selection.");
        }
    }

    DensitySpreader spreader(averageMap_.getGrid(), nThreadsOmp_, nSigma_, 1., integrateDensity_);

    // Spread all atoms with the same Gaussian parameters as one block
    for (const auto &para : atomCoordinatesPerPrecision)
    {
        const auto precision = para.first;
        if (precision > 0)
        {
            auto              weights = atomWeightsPerPrecision[precision];
            const auto        sigma   = sqrt(1/(2*precision));
            spreader.setSigma(sigma);
            // ensure integral over atom density contribution equals its weight
            for (auto &w : weights)
            {
                const auto gaussNorm = sqrt(2*M_PI*sigma);
                w *= gaussNorm * gaussNorm * gaussNorm;
            }
            spreader.spreadLocalAtoms(
                    { para.second.data(), para.second.data()+para.second.size() },
                    { weights.data(), weights.data()+weights.size() });
        }
    }

    // update average and stdev maps
    containeroperation::add(&averageMap_, spreader.getSpreadGrid());
    std::transform(std::begin(stddevMap_), std::end(stddevMap_),
                   std::begin(spreader.getSpreadGrid()), std::begin(stddevMap_),
                   [](float oldSD, float voxelValue) {
                       return oldSD + voxelValue * voxelValue;
                   });

    if (!fnCorrelation_.empty())
    {
        containeroperation::add(&blockAverageMap_, spreader.getSpreadGrid());
        nFramesBlock_++;

        if (nFramesBlock_ == blocksize_)
        {
            containeroperation::divide(&blockAverageMap_, blocksize_);
            // Add data points only for the frames where a block average is calculated
            AnalysisDataHandle dh = pdata->dataHandle(goodnessOfFitData_);
            dh.startFrame(dataFrameNumber_, dataFrameNumber_ * blocksize_);
            dh.setPoint(0, correlationCoefficient(referenceMap_, blockAverageMap_, mask_));
            dh.finishFrame();
            ++dataFrameNumber_;
            containeroperation::setZero(&blockAverageMap_);
            nFramesBlock_ = 0;
        }
    }

}

void TimeAveragedDensity::finishAnalysis(int nframes)
{
    containeroperation::divide(&averageMap_, nframes);
    // re-normalize the accumulated standard deviation map using the average map
    if (!fnMapStandardDeviation_.empty() && nframes > 1)
    {
        std::transform(std::begin(averageMap_), std::end(averageMap_),
                       std::begin(stddevMap_), std::begin(stddevMap_),
                       [nframes](float av, float sd) {
                           return sqrt(sd / (nframes - 1.0) -
                                       av * av * nframes / (nframes - 1.0));
                       });
    }
}

void TimeAveragedDensity::writeOutput()
{
    if (!fnMapAverage_.empty())
    {
        if (!fnMaskMap_.empty())
        {
            mask_.maskToValue(&averageMap_, 0.);
        }
        MrcFile().write(fnMapAverage_, averageMap_);
    }
    if (!fnMapStandardDeviation_.empty())
    {
        if (!fnMaskMap_.empty())
        {
            mask_.maskToValue(&stddevMap_, 0.);
        }
        MrcFile().write(fnMapStandardDeviation_, stddevMap_);
    }
}

}       // namespace

const char TimeAveragedDensityInfo::name[]             = "tide";
const char TimeAveragedDensityInfo::shortDescription[] =
    "Calculate TIme averaged electron DEnsities";

TrajectoryAnalysisModulePointer TimeAveragedDensityInfo::create()
{
    return TrajectoryAnalysisModulePointer(new TimeAveragedDensity);
}

} // namespace analysismodules

} // namespace gmx
