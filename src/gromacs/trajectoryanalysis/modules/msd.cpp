/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Defines the trajectory analysis module for mean squared displacement calculations.
 *
 * \author Kevin Boyd <kevin44boyd@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "msd.h"

#include <numeric>
#include <optional>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/analysisdata/paralleloptions.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace analysismodules
{

namespace
{

//! Convert nm^2/ps to 10e-5 cm^2/s
constexpr double c_diffusionConversionFactor = 1000.0;
//! Three dimensional diffusion coefficient multiplication constant.
constexpr double c_3DdiffusionDimensionFactor = 6.0;
//! Two dimensional diffusion coefficient multiplication constant.
constexpr double c_2DdiffusionDimensionFactor = 4.0;
//! One dimensional diffusion coefficient multiplication constant.
constexpr double c_1DdiffusionDimensionFactor = 2.0;


/*! \brief Mean Squared Displacement data accumulator
 *
 * This class is used to accumulate individual MSD data points
 * and emit tau-averaged results once data is finished collecting. Displacements at each observed
 * time difference (tau) are recorded from the trajectory. Because it is not known in advance which
 * time differences will be observed from the trajectory, this data structure is built adaptively.
 * New columns corresponding to observed time differences are added as needed, and additional
 * observations at formerly observed time differences are added to those columns. Separate time lags
 * will likely have differing total data points.
 *
 * Data columns per tau are accessed via operator[], which always guarantees
 * a column is initialized and returns an MsdColumProxy to the column that can push data.
 */
class MsdData
{
public:
    //! Proxy to a MsdData tau column vector. Supports only push_back.
    class MsdColumnProxy
    {
    public:
        MsdColumnProxy(std::vector<double>* column) : column_(column) {}

        void push_back(double value) { column_->push_back(value); }

    private:
        std::vector<double>* column_;
    };
    //! Returns a proxy to the column for the given tau index. Guarantees that the column is initialized.
    MsdColumnProxy operator[](size_t index)
    {
        if (msds_.size() <= index)
        {
            msds_.resize(index + 1);
        }
        return MsdColumnProxy(&msds_[index]);
    }
    /*! \brief Compute per-tau MSDs averaged over all added points.
     *
     * The resulting vector is size(max tau index). Any indices
     * that have no data points have MSD set to 0.
     *
     * \return Average MSD per tau
     */
    [[nodiscard]] std::vector<real> averageMsds() const;

private:
    //! Results - first indexed by tau, then data points
    std::vector<std::vector<double>> msds_;
};


std::vector<real> MsdData::averageMsds() const
{
    std::vector<real> msdSums;
    msdSums.reserve(msds_.size());
    for (gmx::ArrayRef<const double> msdValues : msds_)
    {
        if (msdValues.empty())
        {
            msdSums.push_back(0.0);
            continue;
        }
        msdSums.push_back(std::accumulate(msdValues.begin(), msdValues.end(), 0.0, std::plus<>())
                          / msdValues.size());
    }
    return msdSums;
}

/*! \brief Calculates 1,2, or 3D distance for two vectors.
 *
 * \todo Remove NOLINTs once clang-tidy is updated to v11, it should be able to handle constexpr.
 *
 * \tparam x If true, calculate x dimension of displacement
 * \tparam y If true, calculate y dimension of displacement
 * \tparam z If true, calculate z dimension of displacement
 * \param[in] c1 First point
 * \param[in] c2 Second point
 * \return Euclidian distance for the given dimension.
 */
template<bool x, bool y, bool z>
inline double calcSingleSquaredDistance(const RVec c1, const RVec c2)
{
    static_assert(x || y || z, "zero-dimensional MSD selected");
    const DVec firstCoords  = c1.toDVec();
    const DVec secondCoords = c2.toDVec();
    double     result       = 0;
    if constexpr (x)
    {
        result += (firstCoords[XX] - secondCoords[XX]) * (firstCoords[XX] - secondCoords[XX]);
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (y)
    {
        result += (firstCoords[YY] - secondCoords[YY]) * (firstCoords[YY] - secondCoords[YY]);
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    if constexpr (z)
    {
        result += (firstCoords[ZZ] - secondCoords[ZZ]) * (firstCoords[ZZ] - secondCoords[ZZ]);
    }
    // NOLINTNEXTLINE(readability-misleading-indentation) remove when clang-tidy-13 is required
    return result;
}

/*! \brief Calculate average displacement between sets of points
 *
 * Each displacement c1[i] - c2[i] is calculated and the distances
 * are averaged.
 *
 * \tparam x If true, calculate x dimension of displacement
 * \tparam y If true, calculate y dimension of displacement
 * \tparam z If true, calculate z dimension of displacement
 * \param[in] c1 First vector
 * \param[in] c2 Second vector
 * \return Per-particle averaged distance
 */
template<bool x, bool y, bool z>
double calcAverageDisplacement(ArrayRef<const RVec> c1, ArrayRef<const RVec> c2)
{
    double result = 0;
    for (size_t i = 0; i < c1.size(); i++)
    {
        result += calcSingleSquaredDistance<x, y, z>(c1[i], c2[i]);
    }
    return result / c1.size();
}


//! Describes 1D MSDs, in the given dimension.
enum class SingleDimDiffType : int
{
    X = 0,
    Y,
    Z,
    Unused,
    Count,
};

//! Describes 2D MSDs, in the plane normal to the given dimension.
enum class TwoDimDiffType : int
{
    NormalToX = 0,
    NormalToY,
    NormalToZ,
    Unused,
    Count,
};

/*! \brief Removes jumps across periodic boundaries for currentFrame, based on the positions in
 * previousFrame. Updates currentCoords in place.
 */
void removePbcJumps(ArrayRef<RVec> currentCoords, ArrayRef<const RVec> previousCoords, t_pbc* pbc)
{
    // There are two types of "pbc removal" in gmx msd. The first happens in the trajectoryanalysis
    // framework, which makes molecules whole across periodic boundaries and is done
    // automatically where the inputs support it. This lambda performs the second PBC correction, where
    // any "jump" across periodic boundaries BETWEEN FRAMES is put back. The order of these
    // operations is important - since the first transformation may only apply to part of a
    // molecule (e.g., one half in/out of the box is put on one side of the box), the
    // subsequent step needs to be applied to the molecule COM rather than individual atoms, or
    // we'd have a clash where the per-mol PBC removal moves an atom that gets put back into
    // it's original position by the second transformation. Therefore, this second transformation
    // is applied *after* per molecule coordinates have been consolidated into COMs.
    auto pbcRemover = [pbc](RVec in, RVec prev) {
        rvec dx;
        pbc_dx(pbc, in, prev, dx);
        return prev + dx;
    };
    std::transform(
            currentCoords.begin(), currentCoords.end(), previousCoords.begin(), currentCoords.begin(), pbcRemover);
}

//! Holds data needed for MSD calculations for a single molecule, if requested.
struct MoleculeData
{
    //! Number of atoms in the molecule.
    int atomCount = 0;
    //! Total mass.
    double mass = 0;
    //! MSD accumulator and calculator for the molecule
    MsdData msdData;
    //! Calculated diffusion coefficient
    real diffusionCoefficient = 0;
};

/*! \brief Handles coordinate operations for MSD calculations.
 *
 * Can be used to hold coordinates for individual atoms as well as molecules COMs. Handles PBC
 * jump removal between consecutive frames.
 *
 * Only previous_ contains valid data between calls to buildCoordinates(), although both vectors
 * are maintained at the correct size for the number of coordinates needed.
 */
class MsdCoordinateManager
{
public:
    MsdCoordinateManager(const int                    numAtoms,
                         ArrayRef<const MoleculeData> molecules,
                         ArrayRef<const int>          moleculeIndexMapping) :
        current_(numAtoms), previous_(numAtoms), molecules_(molecules), moleculeIndexMapping_(moleculeIndexMapping)
    {
    }
    /*! \brief Prepares coordinates for the current frame.
     *
     * Reads in selection data, and returns an ArrayRef of the particle positions or molecule
     * centers of mass (if the molecules input is not empty). Removes jumps across periodic
     * boundaries based on the previous frame coordinates, except for the first frame built with
     * builtCoordinates(), which has no previous frame as a reference.
     *
     * \param[in] sel                   The selection object which holds coordinates
     * \param[in] pbc                   Information about periodicity.
     * \returns                         The current frames coordinates in proper format.
     */
    ArrayRef<const RVec> buildCoordinates(const Selection& sel, t_pbc* pbc);

private:
    //! The current coordinates.
    std::vector<RVec> current_;
    //! The previous frame's coordinates.
    std::vector<RVec> previous_;
    //! Molecule data.
    ArrayRef<const MoleculeData> molecules_;
    //! Mapping of atom indices to molecle indices;
    ArrayRef<const int> moleculeIndexMapping_;
    //! Tracks if a previous frame exists to compare with for PBC handling.
    bool containsPreviousFrame_ = false;
};

ArrayRef<const RVec> MsdCoordinateManager::buildCoordinates(const Selection& sel, t_pbc* pbc)
{

    if (molecules_.empty())
    {
        std::copy(sel.coordinates().begin(), sel.coordinates().end(), current_.begin());
    }
    else
    {
        // Prepare for molecule COM calculation, then sum up all positions per molecule.

        std::fill(current_.begin(), current_.end(), RVec(0, 0, 0));
        gmx::ArrayRef<const real> masses = sel.masses();
        for (int i = 0; i < sel.posCount(); i++)
        {
            const int moleculeIndex = moleculeIndexMapping_[i];
            // accumulate ri * mi, and do division at the end to minimize number of divisions.
            current_[moleculeIndex] += RVec(sel.position(i).x()) * masses[i];
        }
        // Divide accumulated mass * positions to get COM.
        std::transform(current_.begin(),
                       current_.end(),
                       molecules_.begin(),
                       current_.begin(),
                       [](const RVec& position, const MoleculeData& molecule) -> RVec {
                           return position / molecule.mass;
                       });
    }

    if (containsPreviousFrame_)
    {
        removePbcJumps(current_, previous_, pbc);
    }
    else
    {
        containsPreviousFrame_ = true;
    }

    // Previous is no longer needed, swap with current and return "current" coordinates which
    // now reside in previous.
    current_.swap(previous_);
    return previous_;
}


//! Holds per-group coordinates, analysis, and results.
struct MsdGroupData
{
    explicit MsdGroupData(const Selection&             inputSel,
                          ArrayRef<const MoleculeData> molecules,
                          ArrayRef<const int>          moleculeAtomMapping) :
        sel(inputSel),
        coordinateManager_(molecules.empty() ? sel.posCount() : molecules.size(), molecules, moleculeAtomMapping)
    {
    }

    //! Selection associated with this group.
    const Selection& sel;

    //! Stored coordinates, indexed by frame then atom number.
    std::vector<std::vector<RVec>> frames;

    //! MSD result accumulator
    MsdData msds;
    //! Coordinate handler for the group.
    MsdCoordinateManager coordinateManager_;
    //! Collector for processed MSD averages per tau
    std::vector<real> msdSums;
    //! Fitted diffusion coefficient
    real diffusionCoefficient = 0.0;
    //! Uncertainty of diffusion coefficient
    double sigma = 0.0;
};

} // namespace

/*! \brief Implements the gmx msd module
 *
 * \todo Implement -tensor for full MSD tensor calculation
 * \todo Implement -rmcomm for total-frame COM removal
 * \todo Implement -pdb for molecule B factors
 */
class Msd : public TrajectoryAnalysisModule
{
public:
    Msd();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAfterFirstFrame(const TrajectoryAnalysisSettings& settings, const t_trxframe& fr) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int                           frameNumber,
                      const t_trxframe&             frame,
                      t_pbc*                        pbc,
                      TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    //! Selections for MSD output
    SelectionList selections_;

    //! MSD type information, for -type {x,y,z}
    SingleDimDiffType singleDimType_ = SingleDimDiffType::Unused;
    //! MSD type information, for -lateral {x,y,z}
    TwoDimDiffType twoDimType_ = TwoDimDiffType::Unused;
    //! Diffusion coefficient conversion factor
    double diffusionCoefficientDimensionFactor_ = c_3DdiffusionDimensionFactor;
    //! Method used to calculate MSD - changes based on dimensonality.
    std::function<double(ArrayRef<const RVec>, ArrayRef<const RVec>)> calcMsd_ =
            calcAverageDisplacement<true, true, true>;

    //! Picoseconds between restarts
    double trestart_ = 10.0;
    //! Initial time
    double t0_ = 0;
    //! Inter-frame delta-t
    std::optional<double> dt_ = std::nullopt;
    //! Inter-frame maximum time delta.
    double maxTau_ = std::numeric_limits<double>::max();
    //! Tracks the index of the first valid frame. Starts at 0, increases as old frames are deleted
    size_t firstValidFrame_ = 0;

    //! First tau value to fit from for diffusion coefficient, defaults to 0.1 * max tau
    real beginFit_ = -1.0;
    //! Final tau value to fit to for diffusion coefficient, defaults to 0.9 * max tau
    real endFit_ = -1.0;

    //! All selection group-specific data stored here.
    std::vector<MsdGroupData> groupData_;

    //! Time of all frames.
    std::vector<double> times_;
    //! Taus for output - won't know the size until the end.
    std::vector<double> taus_;
    //! Tau indices for fitting.
    size_t beginFitIndex_ = 0;
    size_t endFitIndex_   = 0;

    // MSD per-molecule stuff
    //! Are we doing molecule COM-based MSDs?
    bool molSelected_ = false;
    //! Per molecule topology information and MSD accumulators.
    std::vector<MoleculeData> molecules_;
    //! Atom index -> mol index map
    std::vector<int> moleculeIndexMappings_;

    // Output stuff
    AnalysisData msdPlotData_;
    AnalysisData msdMoleculePlotData_;

    AnalysisDataPlotSettings plotSettings_;
    //! Per-tau MSDs for each selected group
    std::string output_;
    //! Per molecule diffusion coefficients if -mol is selected.
    std::string moleculeOutput_;
};


Msd::Msd() = default;

void Msd::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] computes the mean square displacement (MSD) of atoms from",
        "a set of initial positions. This provides an easy way to compute",
        "the diffusion constant using the Einstein relation.",
        "The time between the reference points for the MSD calculation",
        "is set with [TT]-trestart[tt].",
        "The diffusion constant is calculated by least squares fitting a",
        "straight line (D*t + c) through the MSD(t) from [TT]-beginfit[tt] to",
        "[TT]-endfit[tt] (note that t is time from the reference positions,",
        "not simulation time). An error estimate given, which is the difference",
        "of the diffusion coefficients obtained from fits over the two halves",
        "of the fit interval.[PAR]",
        "There are three, mutually exclusive, options to determine different",
        "types of mean square displacement: [TT]-type[tt], [TT]-lateral[tt]",
        "and [TT]-ten[tt]. Option [TT]-ten[tt] writes the full MSD tensor for",
        "each group, the order in the output is: trace xx yy zz yx zx zy.[PAR]",
        "If [TT]-mol[tt] is set, [THISMODULE] plots the MSD for individual molecules",
        "(including making molecules whole across periodic boundaries): ",
        "for each individual molecule a diffusion constant is computed for ",
        "its center of mass. The chosen index group will be split into ",
        "molecules. With -mol, only one index group can be selected.[PAR]",
        "The diffusion coefficient is determined by linear regression of the MSD.",
        "When [TT]-beginfit[tt] is -1, fitting starts at 10%",
        "and when [TT]-endfit[tt] is -1, fitting goes to 90%.",
        "Using this option one also gets an accurate error estimate",
        "based on the statistics between individual molecules.",
        "Note that this diffusion coefficient and error estimate are only",
        "accurate when the MSD is completely linear between",
        "[TT]-beginfit[tt] and [TT]-endfit[tt].[PAR]",
        "By default, [THISMODULE] compares all trajectory frames against every frame stored at",
        "[TT]-trestart[TT] intervals, so the number of frames stored scales linearly with the",
        "number of frames processed. This can lead to long analysis times and out-of-memory errors",
        "for long/large trajectories, and often the data at higher time deltas lacks sufficient ",
        "sampling, often manifesting as a wobbly line on the MSD plot after a straighter region at",
        "lower time deltas. The [TT]-maxtau[TT] option can be used to cap the maximum time delta",
        "for frame comparison, which may improve performance and can be used to avoid",
        "out-of-memory issues.[PAR]"
    };
    settings->setHelpText(desc);

    // Selections
    options->addOption(SelectionOption("sel")
                               .storeVector(&selections_)
                               .required()
                               .onlyStatic()
                               .multiValue()
                               .description("Selections to compute MSDs for from the reference"));

    // Select MSD type - defaults to 3D if neither option is selected.
    EnumerationArray<SingleDimDiffType, const char*> enumTypeNames    = { "x", "y", "z", "unused" };
    EnumerationArray<TwoDimDiffType, const char*>    enumLateralNames = { "x", "y", "z", "unused" };
    options->addOption(EnumOption<SingleDimDiffType>("type")
                               .enumValue(enumTypeNames)
                               .store(&singleDimType_)
                               .defaultValue(SingleDimDiffType::Unused));
    options->addOption(EnumOption<TwoDimDiffType>("lateral")
                               .enumValue(enumLateralNames)
                               .store(&twoDimType_)
                               .defaultValue(TwoDimDiffType::Unused));

    options->addOption(DoubleOption("trestart")
                               .description("Time between restarting points in trajectory (ps)")
                               .defaultValue(10.0)
                               .store(&trestart_));
    options->addOption(
            DoubleOption("maxtau")
                    .description("Maximum time delta between frames to calculate MSDs for (ps)")
                    .store(&maxTau_));
    options->addOption(
            RealOption("beginfit").description("Time point at which to start fitting.").store(&beginFit_));
    options->addOption(RealOption("endfit").description("End time for fitting.").store(&endFit_));

    // Output options
    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&output_)
                               .defaultBasename("msdout")
                               .description("MSD output"));
    options->addOption(
            FileNameOption("mol")
                    .filetype(OptionFileType::Plot)
                    .outputFile()
                    .store(&moleculeOutput_)
                    .storeIsSet(&molSelected_)
                    .defaultBasename("diff_mol")
                    .description("Report diffusion coefficients for each molecule in selection"));
}

void Msd::optionsFinished(TrajectoryAnalysisSettings gmx_unused* settings)
{
    if (singleDimType_ != SingleDimDiffType::Unused && twoDimType_ != TwoDimDiffType::Unused)
    {
        std::string errorMessage =
                "Options -type and -lateral are mutually exclusive. Choose one or neither (for 3D "
                "MSDs).";
        GMX_THROW(InconsistentInputError(errorMessage));
    }
    if (selections_.size() > 1 && molSelected_)
    {
        std::string errorMessage =
                "Cannot have multiple groups selected with -sel when using -mol.";
        GMX_THROW(InconsistentInputError(errorMessage));
    }
}


void Msd::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    plotSettings_ = settings.plotSettings();

    // Enumeration helpers for dispatching the right MSD calculation type.
    const EnumerationArray<SingleDimDiffType, decltype(&calcAverageDisplacement<true, true, true>)>
            oneDimensionalMsdFunctions = { &calcAverageDisplacement<true, false, false>,
                                           &calcAverageDisplacement<false, true, false>,
                                           &calcAverageDisplacement<false, false, true> };
    const EnumerationArray<TwoDimDiffType, decltype(&calcAverageDisplacement<true, true, true>)>
            twoDimensionalMsdFunctions = { calcAverageDisplacement<false, true, true>,
                                           calcAverageDisplacement<true, false, true>,
                                           calcAverageDisplacement<true, true, false> };

    // Parse dimensionality and assign the MSD calculating function.
    // Note if we don't hit either of these cases, we're computing 3D MSDs.
    if (singleDimType_ != SingleDimDiffType::Unused)
    {
        calcMsd_                             = oneDimensionalMsdFunctions[singleDimType_];
        diffusionCoefficientDimensionFactor_ = c_1DdiffusionDimensionFactor;
    }
    else if (twoDimType_ != TwoDimDiffType::Unused)
    {
        calcMsd_                             = twoDimensionalMsdFunctions[twoDimType_];
        diffusionCoefficientDimensionFactor_ = c_2DdiffusionDimensionFactor;
    }

    // TODO validate that we have mol info and not atom only - and masses, and topology.
    if (molSelected_)
    {
        Selection& sel  = selections_[0];
        const int  nMol = sel.initOriginalIdsToGroup(top.mtop(), INDEX_MOL);

        gmx::ArrayRef<const int> mappedIds = selections_[0].mappedIds();
        moleculeIndexMappings_.resize(selections_[0].posCount());
        std::copy(mappedIds.begin(), mappedIds.end(), moleculeIndexMappings_.begin());

        // Precalculate each molecules mass for speeding up COM calculations.
        ArrayRef<const real> masses = sel.masses();

        molecules_.resize(nMol);
        for (int i = 0; i < sel.posCount(); i++)
        {
            molecules_[mappedIds[i]].atomCount++;
            molecules_[mappedIds[i]].mass += masses[i];
        }
    }

    // Accumulated frames and results
    for (const Selection& sel : selections_)
    {
        groupData_.emplace_back(sel, molecules_, moleculeIndexMappings_);
    }
}

//! Helper function to round \c frame.time to nearest integer and check that it is lossless
static real roundedFrameTime(const t_trxframe& frame)
{
    constexpr real relTolerance = 0.1;
    const real     roundedTime  = std::round(frame.time);
    if (std::fabs(roundedTime - frame.time) > relTolerance * std::fabs(roundedTime + frame.time))
    {
        GMX_THROW(gmx::ToleranceError(gmx::formatString(
                "Frame %" PRId64
                " has non-integral time %f. 'gmx msd' uses time discretization internally and "
                "cannot work if the time (usually measured in ps) is not integral. You can use "
                "'gmx convert-trj -dt 1' to subsample your trajectory before the analysis.",
                frame.step,
                frame.time)));
    }
    return roundedTime;
}

void Msd::initAfterFirstFrame(const TrajectoryAnalysisSettings gmx_unused& settings, const t_trxframe& fr)
{
    t0_ = roundedFrameTime(fr);
}

void Msd::analyzeFrame(int gmx_unused                frameNumber,
                       const t_trxframe&             frame,
                       t_pbc*                        pbc,
                       TrajectoryAnalysisModuleData* pdata)
{
    const real time = roundedFrameTime(frame);
    // Need to populate dt on frame 2;
    if (!dt_.has_value() && !times_.empty())
    {
        dt_ = time - times_[0];
        // Place conditions so they are only checked once
        if (*dt_ > trestart_)
        {
            std::string errorMessage = "-dt cannot be larger than -trestart (default 10 ps).";
            GMX_THROW(InconsistentInputError(errorMessage));
        }
        if (!bRmod(trestart_, 0, *dt_))
        {
            std::string errorMessage =
                    "-trestart (default 10 ps) must be divisible by -dt for useful results.";
            GMX_THROW(InconsistentInputError(errorMessage));
        }
        if (*dt_ == trestart_)
        {
            //\n included below to avoid conflict with other output
            fprintf(stderr,
                    "\nWARNING: -dt and -trestart are equal. Statistics for each tau data point "
                    "will not be independent.\n");
        }
    }

    // Each frame gets an entry in times, but frameTimes only updates if we're at a restart.
    times_.push_back(time);

    // Each frame will get a tau between it and frame 0, and all other frame combos should be
    // covered by this.
    if (const double tau = time - times_[0]; tau <= maxTau_)
    {
        taus_.push_back(time - times_[0]);
    }

    for (MsdGroupData& msdData : groupData_)
    {
        //NOLINTNEXTLINE(readability-static-accessed-through-instance)
        const Selection& sel = pdata->parallelSelection(msdData.sel);

        ArrayRef<const RVec> coords = msdData.coordinateManager_.buildCoordinates(sel, pbc);

        // For each preceding frame, calculate tau and do comparison.
        for (size_t i = firstValidFrame_; i < msdData.frames.size(); i++)
        {

            const double tau = time - (t0_ + trestart_ * i);
            if (tau > maxTau_)
            {
                // The (now empty) entry is no longer needed, so over time the outer vector will
                // grow with extraneous empty elements persisting, but the alternative would require
                // some more complicated remapping of tau to frame index.
                msdData.frames[i].clear();
                firstValidFrame_ = i + 1;
                continue;
            }
            int64_t tauIndex = gmx::roundToInt64(tau / *dt_);
            msdData.msds[tauIndex].push_back(calcMsd_(coords, msdData.frames[i]));

            for (size_t molInd = 0; molInd < molecules_.size(); molInd++)
            {
                molecules_[molInd].msdData[tauIndex].push_back(
                        calcMsd_(arrayRefFromArray(&coords[molInd], 1),
                                 arrayRefFromArray(&msdData.frames[i][molInd], 1)));
            }
        }


        // We only store the frame for the future if it's a restart per -trestart.
        if (bRmod(time, t0_, trestart_))
        {
            msdData.frames.emplace_back(coords.begin(), coords.end());
        }
    }
}

//! Calculate the tau index for fitting. If userFitTau < 0, uses the default fraction of max tau.
static size_t calculateFitIndex(const int    userFitTau,
                                const double defaultTauFraction,
                                const int    numTaus,
                                const double dt)
{
    if (userFitTau < 0)
    {
        return gmx::roundToInt((numTaus - 1) * defaultTauFraction);
    }
    return std::min<size_t>(numTaus - 1, gmx::roundToInt(static_cast<double>(userFitTau) / dt));
}


void Msd::finishAnalysis(int gmx_unused nframes)
{
    static constexpr double c_defaultStartFitIndexFraction = 0.1;
    static constexpr double c_defaultEndFitIndexFraction   = 0.9;
    beginFitIndex_ = calculateFitIndex(beginFit_, c_defaultStartFitIndexFraction, taus_.size(), *dt_);
    endFitIndex_   = calculateFitIndex(endFit_, c_defaultEndFitIndexFraction, taus_.size(), *dt_);
    const int numTausForFit = 1 + endFitIndex_ - beginFitIndex_;

    // These aren't used, except for correlationCoefficient, which is used to estimate error if
    // enough points are available.
    real b = 0.0, correlationCoefficient = 0.0, chiSquared = 0.0;

    for (MsdGroupData& msdData : groupData_)
    {
        msdData.msdSums = msdData.msds.averageMsds();

        if (numTausForFit >= 4)
        {
            const int halfNumTaus         = numTausForFit / 2;
            const int secondaryStartIndex = beginFitIndex_ + halfNumTaus;
            // Split the fit in 2, and compare the results of each fit;
            real a = 0.0, a2 = 0.0;
            lsq_y_ax_b_xdouble(halfNumTaus,
                               &taus_[beginFitIndex_],
                               &msdData.msdSums[beginFitIndex_],
                               &a,
                               &b,
                               &correlationCoefficient,
                               &chiSquared);
            lsq_y_ax_b_xdouble(halfNumTaus,
                               &taus_[secondaryStartIndex],
                               &msdData.msdSums[secondaryStartIndex],
                               &a2,
                               &b,
                               &correlationCoefficient,
                               &chiSquared);
            msdData.sigma = std::abs(a - a2);
        }
        lsq_y_ax_b_xdouble(numTausForFit,
                           &taus_[beginFitIndex_],
                           &msdData.msdSums[beginFitIndex_],
                           &msdData.diffusionCoefficient,
                           &b,
                           &correlationCoefficient,
                           &chiSquared);
        msdData.diffusionCoefficient *= c_diffusionConversionFactor / diffusionCoefficientDimensionFactor_;
        msdData.sigma *= c_diffusionConversionFactor / diffusionCoefficientDimensionFactor_;
    }

    for (MoleculeData& molecule : molecules_)
    {
        std::vector<real> msds = molecule.msdData.averageMsds();
        lsq_y_ax_b_xdouble(numTausForFit,
                           &taus_[beginFitIndex_],
                           &msds[beginFitIndex_],
                           &molecule.diffusionCoefficient,
                           &b,
                           &correlationCoefficient,
                           &chiSquared);
        molecule.diffusionCoefficient *= c_diffusionConversionFactor / diffusionCoefficientDimensionFactor_;
    }
}

void Msd::writeOutput()
{
    // AnalysisData currently doesn't support changing column counts after analysis has started.
    // We can't determine the number of tau values until the trajectory is fully read, so analysis
    // data construction and plotting are done here.
    AnalysisDataPlotModulePointer msdPlotModule(new AnalysisDataPlotModule(plotSettings_));
    msdPlotModule->setFileName(output_);
    msdPlotModule->setTitle("Mean Squared Displacement");
    msdPlotModule->setXLabel("tau (ps)");
    msdPlotModule->setYLabel(R"(MSD (nm\\S2\\N))");
    msdPlotModule->setYFormat(10, 6, 'g');
    for (const auto& group : groupData_)
    {
        const real D = group.diffusionCoefficient;
        if (D > 0.01 && D < 1e4)
        {
            msdPlotModule->appendLegend(formatString(
                    "D[%10s] = %.4f (+/- %.4f) (1e-5 cm^2/s)", group.sel.name(), D, group.sigma));
        }
        else
        {
            msdPlotModule->appendLegend(formatString(
                    "D[%10s] = %.4g (+/- %.4f) (1e-5 cm^2/s)", group.sel.name(), D, group.sigma));
        }
    }
    msdPlotData_.addModule(msdPlotModule);
    msdPlotData_.setDataSetCount(groupData_.size());
    for (size_t i = 0; i < groupData_.size(); i++)
    {
        msdPlotData_.setColumnCount(i, 1);
    }
    AnalysisDataHandle dh = msdPlotData_.startData({});
    for (size_t tauIndex = 0; tauIndex < taus_.size(); tauIndex++)
    {
        dh.startFrame(tauIndex, taus_[tauIndex]);
        for (size_t dataSetIndex = 0; dataSetIndex < groupData_.size(); dataSetIndex++)
        {
            dh.selectDataSet(dataSetIndex);
            dh.setPoint(0, groupData_[dataSetIndex].msdSums[tauIndex]);
        }
        dh.finishFrame();
    }
    dh.finishData();

    if (molSelected_)
    {
        AnalysisDataPlotModulePointer molPlotModule(new AnalysisDataPlotModule(plotSettings_));
        molPlotModule->setFileName(moleculeOutput_);
        molPlotModule->setTitle("Mean Squared Displacement / Molecule");
        molPlotModule->setXLabel("Molecule");
        molPlotModule->setYLabel("D(1e-5 cm^2/s)");
        molPlotModule->setYFormat(10, 0, 'g');
        msdMoleculePlotData_.addModule(molPlotModule);
        msdMoleculePlotData_.setDataSetCount(1);
        msdMoleculePlotData_.setColumnCount(0, 1);
        AnalysisDataHandle molDh = msdMoleculePlotData_.startData({});
        for (size_t moleculeIndex = 0; moleculeIndex < molecules_.size(); moleculeIndex++)
        {
            molDh.startFrame(moleculeIndex, moleculeIndex);
            molDh.setPoint(0, molecules_[moleculeIndex].diffusionCoefficient);
            molDh.finishFrame();
        }
        molDh.finishData();
    }
}


const char                      MsdInfo::name[]             = "msd";
const char                      MsdInfo::shortDescription[] = "Compute mean squared displacements";
TrajectoryAnalysisModulePointer MsdInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Msd);
}


} // namespace analysismodules
} // namespace gmx
