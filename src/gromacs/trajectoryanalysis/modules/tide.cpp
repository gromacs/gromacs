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

#include "gromacs/applied-forces/densityfitting/densityspreader.h"
#include "gromacs/math/container/mask.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/griddata/encompassinggrid.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/fileio/griddataio.h"
#include "gromacs/gmxana/gmx_ana.h"
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
    // less-than operator enables use with std::map
    bool operator<(const GaussianParameters &other) const
    {
        return (other.amplitude < amplitude) || ((other.amplitude == amplitude) && (other.precision < precision));
    }
    float amplitude;
    float precision;
};

//! \brief TODO
float maskedInnerProduct(const GridDataReal3D::container_type &vec1,
                         const GridDataReal3D::container_type &vec2,
                         const Mask                           &mask)
{
    double     result     = 0.0;
    const auto minMapSize = std::min(vec1.size(), vec2.size());
    for (size_t l = 0; l < minMapSize; l++)
    {
        if (mask[l])
        {
            result += vec1[l] * vec2[l];
        }
    }
    return result;
}

//! \brief TODO
float correlationCoefficient(const GridDataReal3D::container_type &map1,
                             const GridDataReal3D::container_type &map2,
                             const Mask                           &mask)
{
    return maskedInnerProduct(map1, map2, mask) / sqrt(maskedInnerProduct(map2, map2, mask) *
                                                       maskedInnerProduct(map1, map1, mask));
}

//! \brief Calculate the local correlation
void local_correlation_coeff(const GridDataReal3D &map1, const GridDataReal3D &map2, const Mask &mask,
                             GridDataReal3D *localcorrelationmap, int Nodd)
{
    int Nhalf = static_cast<int>(Nodd / 2); // center point of the sliding grid (e.g. if Nodd=5
    // ; then Nhalf=2, Note that it starts counting from 0)

    // Initialize the output localcorrelationmap map, setting the borders and masked points,
    // with non-assigned NAN values
    std::fill(std::begin(*localcorrelationmap), std::end(*localcorrelationmap), nanf(""));
    const auto &mapExtend = map1.getGrid().lattice().extend();
    const auto &lattice   = map1.getGrid().lattice();
    // loop over the map1 reduced by Nhalf at both ends
    for (int i = Nhalf; i < mapExtend[XX] - Nhalf; i++)
    {
        for (int j = Nhalf; j < mapExtend[YY] - Nhalf; j++)
        {
            for (int k = Nhalf; k < mapExtend[ZZ] - Nhalf; k++)
            {
                ColumnMajorLattice<DIM>::MultiIndex currentIndex({{i, j, k}});

                // reset correlation variables
                float sum12 = 0.0;
                float sum11 = 0.0;
                float sum22 = 0.0;
                // slide over Nodd points, centered at Nhalf to accumulate the density
                for (int i_sli = i - Nhalf; i_sli <= i + Nhalf; i_sli++)
                {
                    for (int j_sli = j - Nhalf; j_sli <= j + Nhalf; j_sli++)
                    {
                        for (int k_sli = k - Nhalf; k_sli <= k + Nhalf; k_sli++)
                        {
                            ColumnMajorLattice<DIM>::MultiIndex localSurroundingIndex({{i_sli, j_sli, k_sli}});
                            if (mask[lattice.lineariseVectorIndex(localSurroundingIndex)])
                            {
                                const auto v1 = *(map1.iteratorAtMultiIndex(localSurroundingIndex));
                                const auto v2 = *(map2.iteratorAtMultiIndex(localSurroundingIndex));
                                sum12 += v1 * v2;
                                sum11 += v1 * v1;
                                sum22 += v2 * v2;
                            }
                        }
                    }
                }

                *(localcorrelationmap->iteratorAtMultiIndex(currentIndex)) = sum12 / sqrt(sum22 * sum11);
            }
        }
    }
}

/*! \brief
 * Template class to serve as a basis for user analysis tools.
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
        void writeMaskedMap(const std::string &outputFileName, const GridDataReal3D &map, const Mask &mask);
        void setAttypeStringsFromTop(const t_topology &top);

        float                    maskCutoff_   = 1e-5; //< rho < mask_rheshold will be ignored from the calculation
        float                    timeWindow_   = 1.0;  //< for time-resolved global correlation
        float                    sigma_        = 0.3;  //< spread with in nm
        bool                     sigmaIsSet_   = false;
        bool                     spacingIsSet_ = false;
        int                      nSlide_       = 5;   // Size of the sliding grid for local-correlation calculation
        float                    spacing_      = 0.2; // grid spacing in nm
        float                    margin_       = 0.5; // grid margin around atoms in nm
        float                    nThreadsOmp_  = 1;   // number of openMP threads in analysis
        float                    nSigma_       = 5.0; // width of spreading in sigma

        SelectionList            sel_;

        AnalysisData             goodnessOfFitData_;

        std::string              fnGaussParameters_;      // input structure factor file
        std::string              fnMapAverage_;           // average map (3d)
        std::string              fnMapStandardDeviation_; // standard deviation (3d)
        std::string              fnReferenceMap_;
        std::string              fnDifferenceMap_;
        std::string              fnMaskMap_;
        std::string              fnLogRatioMap_;
        std::string              fnLocalCorrelationMap_;
        std::string              fnCorrelation_;

        bool                     verbose_          = false;
        bool                     normalizeMaps_    = true;                                     // Normalization
        bool                     shiftLogMap_      = true;                                     // remove mean value from mlog map
        bool                     integrateDensity_ = false;
        std::map < int, std::vector < GaussianParameters>> gaussianParametersPerAtomicNumber_; // Ai factors (Ngaussians x Natoms 1D array)
        std::vector<int>         atomicNumbersIndex_;                                          // array will store the index of the

        // === define the lattice that will store the density for each point===
        GridDataReal3D    referenceMap_;     /* The reference=experimental map to compare to        */
        GridDataReal3D    averageMap_;       /* The computed timed-averaged map*/
        GridDataReal3D    stddevMap_;        /* The computed timed-averaged standard deviation map */
        GridDataReal3D    blockAverageMap_;  /* map for block averaging*/
        Mask              mask_;             //< Masking out density values
        CommandLineOutput output_;           //< Controls output to command line

        float             nFramesBlock_ = 0; // Nframes in the timeWindow_ for block averaging
};

//! \brief Set the Gaussian coefficients for cryo-EM scattering,
std::map < int, std::vector < GaussianParameters>> TimeAveragedDensity::setGaussianCoeff(bool hasSigma)
{
    // read gaussian parameters from file if provided, else use defaults
    FILE * parameterFile;
    if (!fnGaussParameters_.empty())
    {
        output_("Reading Gaussian coefficients from " + fnGaussParameters_ );
        parameterFile = gmx_ffopen(fnGaussParameters_.c_str(), "r");
    }
    else
    {
        if (hasSigma)
        {
            output_("Reading electron microscopy Gauss weights from 'emcrosssections.dat");
            parameterFile = libopen("emcrosssections.dat");
        }
        else
        {
            output_("Reading default Gaussian coefficients from 'electronscattering.dat'");
            parameterFile = libopen("electronscattering.dat");
        }
    }
    if (parameterFile == nullptr)
    {
        FileIOError("Cannot open file " + fnGaussParameters_ + " or library file for reading.");
    }
    std::map < int, std::vector < GaussianParameters>> gaussianParametersLookup;
    char line[1000];
    // all non-header lines
    while (get_a_line(parameterFile, line, 1000))
    {
        int   currentAtomType;
        float amplitude;
        float precision;
        if (sscanf(line, " %d %e %e", &currentAtomType, &amplitude, &precision) == 3)
        {
            if (hasSigma)
            {
                gaussianParametersLookup[currentAtomType].push_back({amplitude, static_cast<float>(precision *(1/(2.*sigma_*sigma_)))});
            }
            else
            {
                gaussianParametersLookup[currentAtomType].push_back({amplitude, precision});
            }
        }
    }
    gmx_ffclose(parameterFile);
    return gaussianParametersLookup;
}

/*! \brief Outputs the density map.
 */
void TimeAveragedDensity::writeMaskedMap(const std::string &outputFileName, const GridDataReal3D &map,
                                         const Mask &mask)
{
    if (outputFileName.empty())
    {
        InconsistentInputError("No filename given for writing masked file.");
    }

    auto mapCopy = map; // copy map before masking
    mask.maskToNan(&mapCopy);
    MrcFile().write(outputFileName, mapCopy);
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
        "These two maps can be normalized to sigma units and shifted to zero "
        "mean for easier visualization ([TT]-normalize[tt])",
        "",
        "A reference map ([TT]-refmap[tt]) enables the following comparison options",
        "",
        "* [TT]-corr[tt]: Time resolved global correlation between the computed and the reference map",
        "* [TT]-timeWindow_[tt] The correlation is computed in time windows",
        "* [TT]-diff[tt]: Difference map:  average - reference ",
        "* [TT]-logmap[tt]: log (average /reference) or log (average) without reference",
        "* [TT]-shiftlogmap[tt]: sets the mean value of logmap to zero",
        "* [TT]-localcorrelationmap[tt]: Local correlation between the computed and the reference map",
        "At all grid points, the correlation between average and reference within a cubic grid of (nSlide_,nSlide_,nSlide_) is stored[PAR]",
        "",
        "A masking map can be provided with the -mask option. Grid points in the ",
        "mask, with density smaller than a user-specified",
        "threshold value (-maskcutoff), will be ignored. In the output maps, ",
        "the masked points will be labeled as NAN[PAR]",
        "    NOTE: If the density is intended to be calculated with respect to ",
        "other reference group, a previous fitting step (e.g. with trjconv) ",
        "should precede the gtide calculation.[PAR]",

        "===Examples===:\n", "(1) Time-average density of the Protein in a box "
        "gmx tide -f traj.xtc -s conf.gro -select \'group \"Protein\"\' -map_av\n",
        "(2)Reading non-default Gaussian parameters  gmx tide -f traj.xtc -s conf.gro -Gauss_coef Gauss_coef.dat -map_av",
        "(3)Time-average density of the Protein, without normalization",
        "including zeros in the spatial mean, printing as well the difference to"
        "a reference map",
        "(lattice parameters taken from the reference map):\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -nonorm "
        "-select \'group \"Protein\"\' -map_av -diff\n",
        "(4) Same as example 3, but considering a mask, to ignore grid points "
        "which have a density smaller than 1e-5 in the mask:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor -nonorm "
        " -maskcutoff 1e-5 -select \'group \"Protein\"\' -map_av "
        "-diff\n",
        "(5) Same as example 4, but computing -log[rho(map_av)/rho(map_ref)], "
        "shifting the mlog map such that the background spatial mean is zero:",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -mask mask.xplor "
        "-maskcutoff 1e-5 -shiftlogmap -select \'group \"Protein\"\' -mlog \n",
        "(6) Time-resolved global correlation of the computed map to the "
        "reference map, computing the density in time blocks of 1ps:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -timeWindow_ 1.0 -select "
        "\'group \"Protein\"\' -corr \n",
        "(7) Local correlation of the computed map to the reference map, sliding "
        "a sub-lattice of size (5,5,5) bins through the main lattice:\n",
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -nslide 5 -select \'group "
        "\"Protein\"\' -localcorrelationmap \n",
        "(8) Evaluate the density, difference and -log to reference map  "
        "gmx tide -f traj.xtc -s conf.gro -mi mi.xplor -select \'group "
        "\"Protein\"\' -diff -mo -logmap\n",
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

    options->addOption(FileNameOption("difference").filetype(eftCCP4).outputFile()
                           .store(&fnDifferenceMap_).defaultBasename("difference")
                           .description("Difference density map"));

    options->addOption(FileNameOption("logmap").filetype(eftCCP4).outputFile()
                           .store(&fnLogRatioMap_).defaultBasename("logmap")
                           .description("Map ln(reference/computed)"));

    options->addOption(FileNameOption("correlation").filetype(eftPlot).outputFile()
                           .store(&fnCorrelation_).defaultBasename("correlation")
                           .description("Correlation between computed and reference map"));

    options->addOption(FileNameOption("localcorrelationmap").filetype(eftCCP4).outputFile()
                           .store(&fnLocalCorrelationMap_).defaultBasename("localcorrelation")
                           .description("Local correlation between computed and reference map"));

    options->addOption(SelectionOption("select").storeVector(&sel_).required()
                           .description("Selection to calculate the density"));

    options->addOption(FloatOption("maskcutoff").store(&maskCutoff_)
                           .description("Voxels smaller than maskcutoff in mask map are excluded from calculations"));

    options->addOption(BooleanOption("normalize").store(&normalizeMaps_)
                           .description( "Normalize the output maps to unit standard deviation and zero mean"));

    options->addOption(BooleanOption("integrate").store(&integrateDensity_)
                           .description("Integrate density over voxel Volumes instead of Gauss transform"));

    options->addOption(BooleanOption("shiftlogmap").store(&shiftLogMap_)
                           .description("remove mean value from mlog map"));

    options->addOption(FloatOption("spacing").store(&spacing_).storeIsSet(&spacingIsSet_)
                           .description("Map spacing in nm"));

    options->addOption(FloatOption("timewindow").store(&timeWindow_).timeValue()
                           .description("Time window for block-averaging correlation coefficient"));
    options->addOption(FloatOption("sigma").store(&sigma_).storeIsSet(&sigmaIsSet_).
                           description("Spread width for atoms in nm"));
    options->addOption(FloatOption("nsigma").store(&nSigma_).
                           description("Range in sigma until spreading is truncated"));
    options->addOption(BooleanOption("v").store(&verbose_).description("Be verbose or not."));

    options->addOption(FloatOption("ntomp").store(&nThreadsOmp_).description("Number of OpenMP threads used for analysis"));

    options->addOption(IntegerOption("nslide").store(&nSlide_)
                           .description("Size of the sliding grid for local correlation calculation (Must be an odd number)"));
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
    settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}


void TimeAveragedDensity::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    // make sure we exit when parameter file reading fails before selections are made.
    // to fail as early as possible
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
    if (top.hasFullTopology())
    {
        ArrayRef<const t_atom> atoms({top.topology()->atoms.atom, top.topology()->atoms.atom+top.topology()->atoms.nr});
        for (const auto &atom : atoms)
        {
            atomicNumbersIndex_.push_back(atom.atomnumber);
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
        referenceMap_       = MrcFile().read(fnReferenceMap_);
        referenceGrid       = referenceMap_.getGrid().duplicate();
    }
    else
    {
        rvec      * coordinates;
        matrix      box;
        top.getTopologyConf(&coordinates, box);
        referenceGrid = encompassingGridFromCoordinates({coordinates, coordinates+top.topology()->atoms.nr}, spacing_, margin_).duplicate();
    }

    averageMap_.setGrid(referenceGrid->duplicate());
    stddevMap_.setGrid(referenceGrid->duplicate());
    blockAverageMap_.setGrid(referenceGrid->duplicate());

    /* read mask map if provided */
    if (!fnMaskMap_.empty())
    {
        mask_ = Mask(MrcFile().read(fnMaskMap_), maskCutoff_);
    }

    if (!fnCorrelation_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnCorrelation_);
        plotm->setTitle("Time-resolved correlation between "
                        "computed and reference maps");
        plotm->setSubtitle("Time-window for block-averaging:" + std::to_string(timeWindow_)+" ps" );
        plotm->setXAxisIsTime();
        plotm->setYLabel("Correlation");
        goodnessOfFitData_.addModule(plotm);
    }
}

void TimeAveragedDensity::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc * /*pbc*/, TrajectoryAnalysisModuleData *pdata)
{


    const size_t       group = 0; // selection-group index (only one considered thus it is 0)
    const Selection   &sel   = pdata->parallelSelection(sel_[group]);

    // Gather all atoms with the same Gaussian Parameters, also considering multiple Gaussian Parameters per atom
    std::map < GaussianParameters, std::vector < RVec>> atomIndexPerParameter;
    for (int atomIndex = 0; atomIndex < sel.atomCount(); atomIndex++)
    {
        // try to look up the gaussian parameters corresponding to the current atoms atomic number
        int atomicNumber = 0;
        // if atomic numbers are not available, the atomic number stays zero
        if (!atomicNumbersIndex_.empty())
        {
            try
            {
                atomicNumber = atomicNumbersIndex_.at(sel.atomIndices()[atomIndex]);
            }
            catch (const std::out_of_range &oor)
            {
                output_.warning("Cannot determine atomic number of atom #" + std::to_string(sel.atomIndices()[atomIndex]) +  ". Falling back to atomic number 0.");
            }
        }
        try
        {
            const auto &gaussianParameter = gaussianParametersPerAtomicNumber_.at(atomicNumber);
            for (const auto &gp : gaussianParameter)
            {
                atomIndexPerParameter[gp].push_back(sel.position(atomIndex).x());
            }
        }
        catch  (const std::out_of_range &oor)
        {
            output_.warning("Cannot find atomic number " + std::to_string(atomicNumber) + " in spread parameters.");
            output_.warning("Ignoring atom number " + std::to_string(sel.atomIndices()[atomIndex]) + " from selection.");
        }
    }

    DensitySpreader spreader(averageMap_.getGrid(), nThreadsOmp_, nSigma_, 1., integrateDensity_);

    for (const auto &para : atomIndexPerParameter)
    {
        if (para.first.amplitude > 0 && para.first.precision > 0)
        {
            const auto        sigma = sqrt(1/(2*para.first.precision));
            spreader.setSigma(sigma);
            std::vector<real> weights(para.second.size(), para.first.amplitude * sqrt(8*M_PI*M_PI*M_PI*sigma*sigma*sigma));
            ArrayRef<real>    weightsRef {
                weights.data(), weights.data()+weights.size()
            };
            PaddedArrayRef<const RVec> xRef {
                para.second.data(), para.second.data()+para.second.size()
            };
            spreader.spreadLocalAtoms(xRef, weightsRef);
        }
    }

    // update average and stdev maps
    containeroperation::add(&averageMap_, spreader.getSpreadGrid());
    containeroperation::add(&blockAverageMap_, spreader.getSpreadGrid());
    nFramesBlock_++;
    std::transform(std::begin(stddevMap_), std::end(stddevMap_), std::begin(spreader.getSpreadGrid()), std::begin(stddevMap_), [](float oldSD, float voxelValue) { return oldSD + voxelValue * voxelValue; });

    if (!fnCorrelation_.empty())
    {
        // cast to double to ensure condition evaluates correctly for large fr.time and small timeWindow_
        if (static_cast<double>(fr.time) / timeWindow_ == floor(static_cast<double>(fr.time) / timeWindow_) && fr.time > 0)
        {
            containeroperation::divide(&blockAverageMap_, nFramesBlock_);

            // Add data points only for the frames where actually a block average is calculated
            AnalysisDataHandle dh = pdata->dataHandle(goodnessOfFitData_);
            dh.startFrame(frnr, fr.time);
            dh.setPoint(0, correlationCoefficient(referenceMap_, blockAverageMap_, mask_));
            dh.finishFrame();

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
    if (!fnDifferenceMap_.empty())
    {
        auto map_diff = containermeasure::subtract(referenceMap_, averageMap_);
        if (!fnDifferenceMap_.empty())
        {
            writeMaskedMap(fnDifferenceMap_, map_diff, mask_);
        }
    }

    if (!fnLogRatioMap_.empty())
    {
        // with reference map, -ln( averageMap_ / referenceMap_ ), else -ln(averageMap_)
        auto map_mlog = fnReferenceMap_.empty() ? containermeasure::negLog(averageMap_) : containermeasure::logRatio(referenceMap_, averageMap_);
        if (shiftLogMap_)
        {
            containeroperation::addScalar(&map_mlog, -containermeasure::mean(map_mlog));
        }
        writeMaskedMap(fnLogRatioMap_, map_mlog, mask_);
    }

    if (!fnReferenceMap_.empty())
    {
        output_.result("Correlation coefficient between computed and reference maps: " + std::to_string(correlationCoefficient(referenceMap_, averageMap_, mask_)));
        if (!fnLocalCorrelationMap_.empty())
        {
            auto map_localcorrelationmap = averageMap_;
            local_correlation_coeff(referenceMap_, averageMap_, mask_, &map_localcorrelationmap,
                                    nSlide_);
            writeMaskedMap(fnLocalCorrelationMap_, map_localcorrelationmap, mask_);
        }
    }

    if (!fnMapAverage_.empty())
    {
        if (normalizeMaps_)
        {
            containeroperation::addScalar(&averageMap_, -containermeasure::mean(averageMap_));
            containeroperation::divide(&averageMap_, containermeasure::rms(averageMap_));
        }
        writeMaskedMap(fnMapAverage_, averageMap_, mask_);
    }

    if (!fnMapStandardDeviation_.empty())
    {
        if (normalizeMaps_)
        {
            containeroperation::addScalar(&stddevMap_, -containermeasure::mean(stddevMap_));
            containeroperation::divide(&stddevMap_, containermeasure::rms(stddevMap_));
        }
        writeMaskedMap(fnMapStandardDeviation_, stddevMap_, mask_);
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
