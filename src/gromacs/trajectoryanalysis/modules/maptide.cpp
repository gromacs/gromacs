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
 * Implements gmx::analysismodules::maptide.
 *
 * \author Camilo Aponte <ca.aponte.uniandes.edu.co>
 * \author Rodolfo Briones <fitobriones@gmail.com>
 * \author Carsten Kutzner <ckutzne@gwdg.de>
 * \author Christian Blau <cblau@gwdg.de>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "maptide.h"

#include <map>
#include <string.h>

#include "gromacs/math/container/mask.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/encompassinggrid.h"
#include "gromacs/math/griddata/operations/densityspreader.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/trajectoryanalysis.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
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

        //! \brief Set the Gaussian coefficients
        std::map < int, std::vector < GaussianParameters>> setGaussianCoeff(bool hasSigma);
        void setAttypeStringsFromTop(const t_topology &top);

        //! Verbosity of the tool
        bool                        verbose_          = false;
        //! Output channel to command line.
        CommandLineOutput           output_;

        //! The library file with the spreading parameters for X-ray
        const std::string eScatteringDatabase_ = "electronscattering.dat";
        //! The library file with the spreading files for cryo-EM conditions (150eV)
        const std::string cryoEmScatteringDatabase_ = "emcrosssections.dat";

        //! rho < mask_rheshold will be ignored from the calculation
        float                       maskCutoff_   = 1e-5;
        //! number of frames to average for time-resolved global correlation
        int                         blocksize_    = 1;
        //! spread width in nm
        float                       sigma_        = 0.3;
        //! tracks if -sigma flag was set
        bool                        sigmaIsSet_   = false;
        //!  grid spacing in nm
        float                       spacing_      = 0.2;
        //! track if -spacing flag was set
        bool                        spacingIsSet_ = false;
        //! grid margin around atoms in structure in nm
        float                       margin_       = 0.5;
        //! number of openMP threads in analysis
        int                         nThreadsOmp_  = 1;
        //! width of spreading in sigma
        float                       nSigma_       = 5.0;

        //! Selection of atoms used for spreading
        SelectionList                 sel_;
        //! The atom charges when using -type charge
        std::vector<real>             atomCharges_;
        //! The van der Waals radii of all atoms
        std::map<int, real>           vdWRadii_;

        //! The cross-correlation that is output to the xvg file
        AnalysisData                goodnessOfFitData_;

        //! input structure factor file name
        std::string                 fnGaussParameters_;
        //! average map output file name
        std::string                 fnMapAverage_;
        //! standard deviation output file name
        std::string                 fnMapStandardDeviation_;
        //! Reference map file name
        std::string                 fnReferenceMap_;
        //! Mask map file name
        std::string                 fnMaskMap_;
        //! File name for correlation output file
        std::string                 fnCorrelation_;

        //! Value used to fill regions not in the mask
        float maskValue_ = std::numeric_limits<float>::quiet_NaN();

        //! Indicates if Gaussians will be intergated over voxels or value at voxel center is evaluated
        bool                        integrateDensity_ = false;
        //! Amplitude and precision for of the Gaussian used for spreading
        std::map < int, std::vector < GaussianParameters>> gaussianParametersPerAtomicNumber_;
        //! Look-up table that translates between atom name and atomic number
        std::map<std::string, real> atomicNumberTable_;
        /*!\brief The atomic number of each atom in the input.
         * Determined only once when reading the topolgy.
         */
        std::vector<int>            atomicNumbersIndex_;
        //! The reference=experimental map to compare to
        GridDataReal3D              referenceMap_;
        //! The computed timed-averaged map
        GridDataReal3D              averageMap_;
        //! The computed timed-averaged standard deviation map
        GridDataReal3D              stddevMap_;
        //! The block-averaged map
        GridDataReal3D              blockAverageMap_;
        //< The mask used for masking out density values
        Mask                        mask_;

        //! Number of frames per averaging block
        float             nFramesBlock_    = 0;
        //! enumerates output data frames, increases in frames where data is output
        int               dataFrameNumber_ = 0;
        //! Physical model for spreading the atoms
        DensityType       densityType_     = enXRay;
};

TimeAveragedDensity::TimeAveragedDensity()
{
    registerAnalysisDataset(&goodnessOfFitData_, "fitness");
}

void TimeAveragedDensity::initOptions(IOptionsContainer          *options,
                                      TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] transforms structures or trajectories into"
        " three-dimensional density data."
        " The density is calculated as a linear combination of Gaussians"
        " rho(R)=sum_{i=1}^{nGaussians} Ai*exp(-Bi*R^2),"
        " where amplitude ([IT]A[it]) and precision ([IT]B[it])"
        " depend on the atomic number of the atom to be spread."
        " By default ([TT]-type x-ray[tt]), 4 Gaussian coefficients from"
        " are read for each atomic number from a database file for H,C,N,O,P,S.",
        "",
        "For [TT]-type cryo-EM[tt], electron scattering factors at 150keV"
        " for all elements are read. Chosing [TT]-type charge[tt] uses the partial"
        " charges from a force field as atom weights for scattering and a"
        " constant scattering Gaussian width which may be controlled with the"
        " [TT]-sigma[tt] flag. The [TT]-type vdW[tt] uses the van-der-Waals"
        " radii as used elsewhere in gromacs as spreading width and all atom"
        " weights set to unity.",
        "",
        "Other spreading parameters A, B, are read via the "
        " [TT]-gaussparameters[tt] input file as \"atomic-number A B\". "
        " Manual asigments of atom names to atomic numbers may be enforced by"
        " adding a line \"atomic-number name\". Atom names match if they"
        " are longer than name, e.g., NH2 will match N. Non-matched atoms are"
        " silently ignored unless the [TT]-v[tt] flag is set. Precision is"
        " always interpreted as given in nm.",
        "",
        "If Ai or Bi is not implicated via [TT]-type[tt] or"
        " [TT]-gaussparameters[tt], [TT]-sigma[tt] is used for the spreading",
        " width and one for atom weights."
        "",
        "[TT]-integrate[tt] replaces Gaussian spreading with integration over"
        " Gaussians. This is beneficial for coarse grids where"
        " the spreading width is smaller than the grid spacing.",
        "",
        "The [TT]-refmap[tt] is used to evaluate the cross-correlation and"
        " determines the three dimensional grid for spreading the atoms."
        " If no reference map is given, the grid is determined from the"
        " input structure extend. [TT]-spacing[tt] may then be used to"
        " adjust grid spacing. With a fixed spread width [TT]-sigma[tt],"
        " grid spacing defaults to one third of [TT]-sigma[tt], ensuring"
        " a small numercial difference between integrated Gaussians and"
        " Gaussian spreading. [TT]-margin[tt] sets the margin around the"
        " structure when determining the grid (in nm).",
        "",
        "Atom spreading is truncated after [TT]-nsigma[tt] multiples of"
        " [TT]-sigma[tt] or the corresponding Precision B if no [TT]-sigma[tt]"
        " is given. Setting [TT]-nsigma[tt] to a lower value than five, will"
        " significantly increase performance at the cost of truncation errors."
        " [TT]-nsigma[tt] larger than five will yield almost no benefit due to"
        " limited floating point precision.",
        "",
        "Performance may be increased [TT]-ntomp[tt] by using more than one"
        " openmp thread for spreading the density, note though that scaling"
        " is very sub-linear.",
        "",
        "From a trajectory or structure alone, the following output "
        " may be generated:",
        "",
        "* [TT]-mo[tt]: time average of density.",
        "* [TT]-stddevmap[tt]: time standard deviation of density.",
        "",
        "A reference map ([TT]-refmap[tt]) enables additional comparison"
        " options",
        "",
        "* [TT]-corr[tt]: Time resolved correlation between computed and the"
        " reference map.",
        "* [TT]-blocksize[tt] Block averaged density correlation with windows"
        " of blocksize. Block size one corresponds to not averaging.",
        "",
        "A masking map is provided with [TT]-mask[tt]. Grid points in the "
        " mask, with density smaller than a user-specified"
        " threshold value ([TT]-maskcutoff[tt]), will be ignored. In the output"
        " maps, masked points will be set to nan or any other, user specified"
        " value ([TT]-maskvalue[tt]).",
        "",
        "NOTE: [THISMODULE] will not perform any translational or rotational "
        " fitting to the reference density or along the trajectory. "
        " A previous fitting step (e.g. with trjconv) should precede the calculation.",
        "",
        "NOTE: [THISMODULE] tries hard to infer the atomic number for all atoms"
        " in the structure file, beware however that it does not know the"
        " chemical context nor user intent and thus might fail with exotic"
        " naming schemes or virtual sites. Enable [TT]-v[tt] to gain more"
        " information about which atom names could not be inferred."
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
                           .description("Voxels smaller than maskcutoff in mask"
                                        " map are excluded from calculations"));

    options->addOption(FloatOption("maskvalue").store(&maskValue_).
                           description("Value that indicates out-of-mask regions."));

    options->addOption(BooleanOption("integrate").store(&integrateDensity_)
                           .description("Integrate Gaussians over voxel volumes"
                                        " instead of Gauss transform"));

    options->addOption(FloatOption("spacing").store(&spacing_).storeIsSet(&spacingIsSet_)
                           .description("Map spacing in nm"));

    options->addOption(FloatOption("margin").store(&margin_)
                           .description("Margin of the grid boundary from the structure if no refernce grid is given."));

    options->addOption(IntegerOption("blocksize").store(&blocksize_)
                           .description("Number of frames for block-averaging"
                                        " correlation coefficient"));

    options->addOption(FloatOption("sigma").store(&sigma_).storeIsSet(&sigmaIsSet_).
                           description("Spread width for atoms in nm"));

    options->addOption(FloatOption("nsigma").store(&nSigma_).
                           description("Range in sigma until spreading is truncated"));

    options->addOption(BooleanOption("v").store(&verbose_).description("Verboseness"));

    options->addOption(IntegerOption("ntomp").store(&nThreadsOmp_)
                           .description("Number of OpenMP threads for analysis"));

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
        GMX_THROW(InconsistentInputError("Cannot use ' -type charge ' without "
                                         "full topology information (i.e. a tpr file), because this is the "
                                         "only way to determine the partial charge per atom."));
    }
    if (top.hasFullTopology() && atomicNumberTable_.empty())
    {
        ArrayRef<const t_atom> atoms({top.topology()->atoms.atom,
                                      top.topology()->atoms.atom+top.topology()->atoms.nr});
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
            output_.warning("Guessing atomic number from atom names");
            const auto &atoms = top.topology()->atoms;
            auto        aps   = gmx_atomprop_init();
            for (int i = 0; i < atoms.nr; ++i)
            {
                real atomicNumber;
                gmx_atomprop_query(aps, epropElement, "???", *atoms.atomname[i],
                                   &atomicNumber);
                if (atomicNumber < 0 || !atomicNumberTable_.empty())
                {
                    // try to match the atom name to the longest prefix in the atomicNumberTable
                    std::string atomName = *atoms.atomname[i];
                    while (atomicNumberTable_.find(atomName) ==
                           atomicNumberTable_.end() && atomName.size() != 0)
                    {
                        atomName.pop_back();
                    }
                    if (atomName.size() == 0)
                    {
                        output_.warning("Cannot guess atomic number for '" +
                                        std::string(*atoms.atomname[i]) +
                                        "' - ignoring for spreading. "
                                        "Is this a virtual site?");
                    }
                    else
                    {
                        atomicNumber = atomicNumberTable_.find(atomName)->second;
                    }
                }
                if (densityType_ == enVdW)
                {
                    real vdWRadius;
                    gmx_atomprop_query(aps, epropVDW, "???",
                                       *atoms.atomname[i], &vdWRadius);
                    if (vdWRadius > 0)
                    {
                        vdWRadii_[atomicNumber] = vdWRadius;
                    }
                    else
                    {
                        output_.warning("Cannot find vdW radius for this atom type, setting radius to zero");
                        vdWRadii_[atomicNumber] = 0;
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
            GMX_THROW(FileIOError("Correlation calculation requires a reference map.\n"
                                  " If you intended to calculate the density correlation"
                                  " between a reference structure and a trajectory, first"
                                  " generate a reference density from the reference structure."));
        }
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnCorrelation_);
        plotm->setTitle("Time-resolved correlation between computed and reference maps");
        plotm->setSubtitle("Number of frames for block-averaging:" + std::to_string(blocksize_)+" frames" );
        plotm->setXLabel("Frame Number");
        plotm->setYLabel("Correlation");
        goodnessOfFitData_.addModule(plotm);
    }
}

void TimeAveragedDensity::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc * /*pbc*/, TrajectoryAnalysisModuleData *pdata)
{
    const size_t       group = 0; // selection-group index (only one considered thus it is 0)
    const Selection   &sel   = pdata->parallelSelection(sel_[group]);

    // Gather all atoms with the same Gaussian Parameters, also considering multiple Gaussian Parameters per atom
    std::map < float, std::vector < RVec>> atomCoordinatesPerPrecision;
    std::map < float, std::vector < real>> atomWeightsPerPrecision;

    // For the density spreading we always use all atoms in the selection,
    // regdless of selection type, because it is unclear how spreading parameters for different atoms should be combined
    for (auto selectedAtomIndex : sel.atomIndices())
    {
        // try to look up the gaussian parameters corresponding to the current atoms atomic number
        int atomicNumber = 0;
        // if atomic numbers are not available, the atomic number stays zero
        if (0 <= selectedAtomIndex && selectedAtomIndex < static_cast<int>(atomicNumbersIndex_.size()))
        {
            atomicNumber = atomicNumbersIndex_[selectedAtomIndex];
        }
        else
        {
            output_.warning("Cannot determine atomic number of atom #" + std::to_string(selectedAtomIndex) +  ". Falling back to atomic number 0.");
        }

        const auto &gaussianParameter = gaussianParametersPerAtomicNumber_.find(atomicNumber);
        if (gaussianParameter !=  gaussianParametersPerAtomicNumber_.end())
        {
            for (const auto &gp : gaussianParameter->second)
            {
                // use fr.x instead of sel.coordinates,
                // because for some selection types (e.g., center of mass),
                // sel.coordinates() are not the atom positions of sel.atomIndices()
                atomCoordinatesPerPrecision[gp.precision].push_back(fr.x[selectedAtomIndex]);
                real weight;
                switch (densityType_)
                {
                    case enCharge:
                        weight = atomCharges_[selectedAtomIndex];
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
        else
        {
            output_.warning("Cannot find atomic number " + std::to_string(atomicNumber) + " in spread parameters.");
            output_.warning("Ignoring atom number " + std::to_string(selectedAtomIndex) + " from selection.");
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
            mask_.maskToValue(&averageMap_, maskValue_);
        }
        MrcFile().write(fnMapAverage_, averageMap_);
    }

    if (!fnMapStandardDeviation_.empty())
    {
        if (!fnMaskMap_.empty())
        {
            mask_.maskToValue(&stddevMap_, maskValue_);
        }
        MrcFile().write(fnMapStandardDeviation_, stddevMap_);
    }
}

std::map < int, std::vector < GaussianParameters>>
TimeAveragedDensity::setGaussianCoeff(bool hasSigma)
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
                output_("Reading default Gaussian coefficients from '"
                        + eScatteringDatabase_ + "' .");
                parameterFile = libopen(eScatteringDatabase_.c_str());
                break;
            case enCryoEM:
                output_("Reading electron microscopy Gauss weights from '"
                        + cryoEmScatteringDatabase_ + "'");
                parameterFile = libopen(cryoEmScatteringDatabase_.c_str());
                break;
            default:
                break;
        }
    }
    if (parameterFile == nullptr && !hasSigma && densityType_ != enVdW)
    {
        GMX_THROW(FileIOError("Cannot open file " + fnGaussParameters_ +
                              " or library file for reading and no spreading"
                              " option -sigma was given. Exiting now because"
                              " cannot reasonably guess how to spread atom density."));
    }

    std::map < int, std::vector < GaussianParameters>> gaussianParametersLookup;

    if (parameterFile != nullptr)
    {
        constexpr size_t maxLineLength = 1000;
        char             line[maxLineLength];
        while (get_a_line(parameterFile, line, maxLineLength))
        {
            float precision;
            int   iAtomType;
            float amplitude;
            if (sscanf(line, " %d %e %e", &iAtomType, &amplitude, &precision) != 3)
            {
                if (sscanf(line, " %d %e", &iAtomType, &amplitude) != 2)
                {
                    char atomName[maxLineLength];
                    if (sscanf(line, "%d %s", &iAtomType, atomName) == 2)
                    {
                        atomicNumberTable_[atomName] = iAtomType;
                    }
                }
                else
                {
                    precision     = 1/(2.*sigma_*sigma_);
                    gaussianParametersLookup[iAtomType].push_back({amplitude, precision});
                }
            }
            else
            {
                gaussianParametersLookup[iAtomType].push_back({amplitude, precision});
            }
        }
        gmx_ffclose(parameterFile);
    }
    else
    {
        constexpr int    maxElementType = 100;
        real             weight         = 1;
        float            precision      = 1/(2.*sigma_*sigma_);
        for (int iAtomType = 0; iAtomType < maxElementType; ++iAtomType)
        {
            if (densityType_ == DensityType::enVdW)
            {
                const auto currentVdWR = vdWRadii_.find(iAtomType);
                if (currentVdWR != vdWRadii_.end())
                {
                    precision = currentVdWR->second;
                }
                else
                {
                    precision = 1.;
                    weight    = 0.;
                }
            }
            gaussianParametersLookup[iAtomType].push_back(
                    {weight, precision});
        }
    }
    return gaussianParametersLookup;
}


}       // namespace

const char TimeAveragedDensityInfo::name[]             = "maptide";
const char TimeAveragedDensityInfo::shortDescription[] =
    "Calculate TIme averaged electron DEnsities";

TrajectoryAnalysisModulePointer TimeAveragedDensityInfo::create()
{
    return TrajectoryAnalysisModulePointer(new TimeAveragedDensity);
}

} // namespace analysismodules

} // namespace gmx
