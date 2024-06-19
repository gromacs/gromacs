/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
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
 * Implements gmx::analysismodules::Angle.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "angle.h"

#include <cstddef>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class TopologyInformation;

namespace analysismodules
{

namespace
{

/********************************************************************
 * Helper classes
 */

/*! \brief
 * Helper to encapsulate logic for looping over input selections.
 *
 * This class provides two-dimensional iteration:
 *  - Over _angle groups_, corresponding to an input selection.  If the input
 *    selection list contains a single selection, that selection gets used
 *    for all angle groups.
 *  - Within a group, over _values_, each consisting of a fixed number of
 *    selection positions.  If there is only a single value within a selection,
 *    that value is returned over and over again.
 * This transparently provides the semantics of using a single selection/vector
 * to compute angles against multiple selections/vectors as described in the
 * tool help text.
 *
 * This class isn't perferctly self-contained and requires the caller to know
 * some of the internals to use it properly, but it serves its purpose for this
 * single analysis tool by simplifying the loops.
 * Some methods have also been tailored to allow the caller to use it a bit
 * more easily.
 *
 * \ingroup module_trajectoryanalysis
 */
class AnglePositionIterator
{
public:
    /*! \brief
     * Creates an iterator to loop over input selection positions.
     *
     * \param[in] selections       List of selections.
     * \param[in] posCountPerValue Number of selection positions that
     *     constitute a single value for the iteration.
     *
     * If \p selections is empty, and/or \p posCountPerValue is zero, the
     * iterator can still be advanced and hasValue()/hasSingleValue()
     * called, but values cannot be accessed.
     */
    AnglePositionIterator(const SelectionList& selections, int posCountPerValue) :
        selections_(selections), posCountPerValue_(posCountPerValue), currentSelection_(0), nextPosition_(0)
    {
    }

    //! Advances the iterator to the next group of angles.
    void nextGroup()
    {
        if (selections_.size() > 1)
        {
            ++currentSelection_;
        }
        nextPosition_ = 0;
    }
    //! Advances the iterator to the next angle in the current group.
    void nextValue()
    {
        if (!hasSingleValue())
        {
            nextPosition_ += posCountPerValue_;
        }
    }

    /*! \brief
     * Returns whether this iterator represents any values.
     *
     * If the return value is `false`, only nextGroup(), nextValue() and
     * hasSingleValue() are allowed to be called.
     */
    bool hasValue() const { return !selections_.empty(); }
    /*! \brief
     * Returns whether the current selection only contains a single value.
     *
     * Returns `false` if hasValue() returns false, which allows cutting
     * some corners in consistency checks.
     */
    bool hasSingleValue() const
    {
        return hasValue() && currentSelection().posCount() == posCountPerValue_;
    }
    //! Returns whether the current selection is dynamic.
    bool isDynamic() const { return currentSelection().isDynamic(); }
    /*! \brief
     * Returns whether positions in the current value are either all
     * selected or all unselected.
     */
    bool allValuesConsistentlySelected() const
    {
        if (posCountPerValue_ <= 1)
        {
            return true;
        }
        const bool bSelected = currentPosition(0).selected();
        for (int i = 1; i < posCountPerValue_; ++i)
        {
            if (currentPosition(i).selected() != bSelected)
            {
                return false;
            }
        }
        return true;
    }
    /*! \brief
     * Returns whether positions in the current value are selected.
     *
     * Only works reliably if allValuesConsistentlySelected() returns
     * `true`.
     */
    bool currentValuesSelected() const
    {
        return selections_.empty() || currentPosition(0).selected();
    }

    //! Returns the currently active selection.
    const Selection& currentSelection() const
    {
        GMX_ASSERT(currentSelection_ < gmx::ssize(selections_), "Accessing an invalid selection");
        return selections_[currentSelection_];
    }
    //! Returns the `i`th position for the current value.
    SelectionPosition currentPosition(int i) const
    {
        return currentSelection().position(nextPosition_ + i);
    }
    /*! \brief
     * Extracts all coordinates corresponding to the current value.
     *
     * \param[out] x  Array to which the positions are extracted.
     *
     * \p x should contain at minimum the number of positions per value
     * passed to the constructor.
     */
    void getCurrentPositions(rvec x[]) const
    {
        GMX_ASSERT(posCountPerValue_ > 0, "Accessing positions for an invalid angle type");
        GMX_ASSERT(nextPosition_ + posCountPerValue_ <= currentSelection().posCount(),
                   "Accessing an invalid position");
        for (int i = 0; i < posCountPerValue_; ++i)
        {
            copy_rvec(currentPosition(i).x(), x[i]);
        }
    }

private:
    const SelectionList& selections_;
    const int            posCountPerValue_;
    int                  currentSelection_;
    int                  nextPosition_;

    GMX_DISALLOW_COPY_AND_ASSIGN(AnglePositionIterator);
};

/********************************************************************
 * Actual analysis module
 */

//! How to interpret the selections in -group1.
enum class Group1Type : int
{
    Angle,
    Dihedral,
    Vector,
    Plane,
    Count
};
//! How to interpret the selections in -group2.
enum class Group2Type : int
{
    None,
    Vector,
    Plane,
    TimeZero,
    Z,
    SphereNormal,
    Count
};
//! String values corresponding to Group1Type.
const EnumerationArray<Group1Type, const char*> c_group1TypeEnumNames = {
    { "angle", "dihedral", "vector", "plane" }
};
//! String values corresponding to Group2Type.
const EnumerationArray<Group2Type, const char*> c_group2TypeEnumNames = {
    { "none", "vector", "plane", "t0", "z", "sphnorm" }
};

class Angle : public TrajectoryAnalysisModule
{
public:
    Angle();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;

    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    void initFromSelections(const SelectionList& sel1, const SelectionList& sel2);
    void checkSelections(const SelectionList& sel1, const SelectionList& sel2) const;

    SelectionList        sel1_;
    SelectionList        sel2_;
    SelectionOptionInfo* sel1info_;
    SelectionOptionInfo* sel2info_;
    std::string          fnAverage_;
    std::string          fnAll_;
    std::string          fnHistogram_;

    Group1Type g1type_;
    Group2Type g2type_;
    double     binWidth_;

    AnalysisData                             angles_;
    AnalysisDataFrameAverageModulePointer    averageModule_;
    AnalysisDataSimpleHistogramModulePointer histogramModule_;

    std::vector<int>               angleCount_;
    int                            natoms1_;
    int                            natoms2_;
    std::vector<std::vector<RVec>> vt0_;

    // Copy and assign disallowed by base.
};

Angle::Angle() :
    sel1info_(nullptr),
    sel2info_(nullptr),
    g1type_(Group1Type::Angle),
    g2type_(Group2Type::None),
    binWidth_(1.0),
    natoms1_(0),
    natoms2_(0)
{
    averageModule_ = std::make_unique<AnalysisDataFrameAverageModule>();
    angles_.addModule(averageModule_);
    histogramModule_ = std::make_unique<AnalysisDataSimpleHistogramModule>();
    angles_.addModule(histogramModule_);

    registerAnalysisDataset(&angles_, "angle");
    registerBasicDataset(averageModule_.get(), "average");
    registerBasicDataset(&histogramModule_->averager(), "histogram");
}


void Angle::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] computes different types of angles between vectors.",
        "It supports both vectors defined by two positions and normals of",
        "planes defined by three positions.",
        "The z axis or the local normal of a sphere can also be used as",
        "one of the vectors.",
        "There are also convenience options 'angle' and 'dihedral' for",
        "calculating bond angles and dihedrals defined by three/four",
        "positions.[PAR]",
        "The type of the angle is specified with [TT]-g1[tt] and [TT]-g2[tt].",
        "If [TT]-g1[tt] is [TT]angle[tt] or [TT]dihedral[tt], [TT]-g2[tt]",
        "should not be specified.",
        "In this case, [TT]-group1[tt] should specify one or more selections,",
        "and each should contain triplets or quartets of positions that define",
        "the angles to be calculated.[PAR]",
        "If [TT]-g1[tt] is [TT]vector[tt] or [TT]plane[tt], [TT]-group1[tt]",
        "should specify selections that contain either pairs ([TT]vector[tt])",
        "or triplets ([TT]plane[tt]) of positions. For vectors, the positions",
        "set the endpoints of the vector, and for planes, the three positions",
        "are used to calculate the normal of the plane. In both cases,",
        "[TT]-g2[tt] specifies the other vector to use (see below).[PAR]",
        "With [TT]-g2 vector[tt] or [TT]-g2 plane[tt], [TT]-group2[tt] should",
        "specify another set of vectors. [TT]-group1[tt] and [TT]-group2[tt]",
        "should specify the same number of selections. It is also allowed to",
        "only have a single selection for one of the options, in which case",
        "the same selection is used with each selection in the other group.",
        "Similarly, for each selection in [TT]-group1[tt], the corresponding",
        "selection in [TT]-group2[tt] should specify the same number of",
        "vectors or a single vector. In the latter case, the angle is",
        "calculated between that single vector and each vector from the other",
        "selection.[PAR]",
        "With [TT]-g2 sphnorm[tt], each selection in [TT]-group2[tt] should",
        "specify a single position that is the center of the sphere.",
        "The second vector is calculated as the vector from the center to the",
        "midpoint of the positions specified by [TT]-group1[tt].[PAR]",
        "With [TT]-g2 z[tt], [TT]-group2[tt] is not necessary, and angles",
        "between the first vectors and the positive Z axis are calculated.[PAR]",
        "With [TT]-g2 t0[tt], [TT]-group2[tt] is not necessary, and angles",
        "are calculated from the vectors as they are in the first frame.[PAR]",
        "There are three options for output:",
        "[TT]-oav[tt] writes an xvg file with the time and the average angle",
        "for each frame.",
        "[TT]-oall[tt] writes all the individual angles.",
        "[TT]-oh[tt] writes a histogram of the angles. The bin width can be",
        "set with [TT]-binw[tt].",
        "For [TT]-oav[tt] and [TT]-oh[tt], separate average/histogram is",
        "computed for each selection in [TT]-group1[tt]."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("oav")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAverage_)
                               .defaultBasename("angaver")
                               .description("Average angles as a function of time"));
    options->addOption(FileNameOption("oall")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAll_)
                               .defaultBasename("angles")
                               .description("All angles as a function of time"));
    options->addOption(FileNameOption("oh")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnHistogram_)
                               .defaultBasename("anghist")
                               .description("Histogram of the angles"));

    options->addOption(EnumOption<Group1Type>("g1")
                               .enumValue(c_group1TypeEnumNames)
                               .store(&g1type_)
                               .description("Type of analysis/first vector group"));
    options->addOption(
            EnumOption<Group2Type>("g2").enumValue(c_group2TypeEnumNames).store(&g2type_).description("Type of second vector group"));
    options->addOption(
            DoubleOption("binw").store(&binWidth_).description("Binwidth for -oh in degrees"));

    sel1info_ = options->addOption(
            SelectionOption("group1").required().dynamicMask().storeVector(&sel1_).multiValue().description(
                    "First analysis/vector selection"));
    sel2info_ = options->addOption(
            SelectionOption("group2").dynamicMask().storeVector(&sel2_).multiValue().description(
                    "Second analysis/vector selection"));
}


void Angle::optionsFinished(TrajectoryAnalysisSettings* /* settings */)
{
    const bool bSingle = (g1type_ == Group1Type::Angle || g1type_ == Group1Type::Dihedral);

    if (bSingle && g2type_ != Group2Type::None)
    {
        GMX_THROW(
                InconsistentInputError("Cannot use a second group (-g2) with "
                                       "-g1 angle or dihedral"));
    }
    if (bSingle && sel2info_->isSet())
    {
        GMX_THROW(
                InconsistentInputError("Cannot provide a second selection "
                                       "(-group2) with -g1 angle or dihedral"));
    }
    if (!bSingle && g2type_ == Group2Type::None)
    {
        GMX_THROW(
                InconsistentInputError("Should specify a second group (-g2) "
                                       "if the first group is not an angle or a dihedral"));
    }

    // Set up the number of positions per angle.
    switch (g1type_)
    {
        case Group1Type::Angle: natoms1_ = 3; break;
        case Group1Type::Dihedral: natoms1_ = 4; break;
        case Group1Type::Vector: natoms1_ = 2; break;
        case Group1Type::Plane: natoms1_ = 3; break;
        default: GMX_THROW(InternalError("invalid -g1 value"));
    }
    switch (g2type_)
    {
        case Group2Type::None: natoms2_ = 0; break;
        case Group2Type::Vector: natoms2_ = 2; break;
        case Group2Type::Plane: natoms2_ = 3; break;
        case Group2Type::TimeZero: // Intended to fall through
        case Group2Type::Z: natoms2_ = 0; break;
        case Group2Type::SphereNormal: natoms2_ = 1; break;
        default: GMX_THROW(InternalError("invalid -g2 value"));
    }
    if (natoms2_ == 0 && sel2info_->isSet())
    {
        GMX_THROW(InconsistentInputError(
                "Cannot provide a second selection (-group2) with -g2 t0 or z"));
    }
    // TODO: If bSingle is not set, the second selection option should be
    // required.
}


void Angle::initFromSelections(const SelectionList& sel1, const SelectionList& sel2)
{
    const int  angleGroups         = std::max(sel1.size(), sel2.size());
    const bool bHasSecondSelection = natoms2_ > 0;

    if (bHasSecondSelection && sel1.size() != sel2.size() && std::min(sel1.size(), sel2.size()) != 1)
    {
        GMX_THROW(InconsistentInputError(
                "-group1 and -group2 should specify the same number of selections"));
    }

    AnglePositionIterator iter1(sel1, natoms1_);
    AnglePositionIterator iter2(sel2, natoms2_);
    for (int g = 0; g < angleGroups; ++g, iter1.nextGroup(), iter2.nextGroup())
    {
        const int posCount1 = iter1.currentSelection().posCount();
        if (natoms1_ > 1 && posCount1 % natoms1_ != 0)
        {
            GMX_THROW(InconsistentInputError(formatString(
                    "Number of positions in selection %d in the first group not divisible by %d",
                    static_cast<int>(g + 1),
                    natoms1_)));
        }
        const int angleCount1 = posCount1 / natoms1_;
        int       angleCount  = angleCount1;

        if (bHasSecondSelection)
        {
            const int posCount2 = iter2.currentSelection().posCount();
            if (natoms2_ > 1 && posCount2 % natoms2_ != 0)
            {
                GMX_THROW(InconsistentInputError(
                        formatString("Number of positions in selection %d in the second group not "
                                     "divisible by %d",
                                     static_cast<int>(g + 1),
                                     natoms2_)));
            }
            if (g2type_ == Group2Type::SphereNormal && posCount2 != 1)
            {
                GMX_THROW(InconsistentInputError(
                        "The second group should contain a single position with -g2 sphnorm"));
            }

            const int angleCount2 = posCount2 / natoms2_;
            angleCount            = std::max(angleCount1, angleCount2);
            if (angleCount1 != angleCount2 && std::min(angleCount1, angleCount2) != 1)
            {
                GMX_THROW(InconsistentInputError(
                        "Number of vectors defined by the two groups are not the same"));
            }
        }
        angleCount_.push_back(angleCount);
    }
}


void Angle::checkSelections(const SelectionList& sel1, const SelectionList& sel2) const
{
    AnglePositionIterator iter1(sel1, natoms1_);
    AnglePositionIterator iter2(sel2, natoms2_);
    for (size_t g = 0; g < angleCount_.size(); ++g, iter1.nextGroup(), iter2.nextGroup())
    {
        if (iter1.isDynamic() || (iter2.hasValue() && iter2.isDynamic()))
        {
            for (int n = 0; n < angleCount_[g]; ++n, iter1.nextValue(), iter2.nextValue())
            {
                bool bOk = true;
                if (!iter1.allValuesConsistentlySelected())
                {
                    bOk = false;
                }
                if (!iter2.allValuesConsistentlySelected())
                {
                    bOk = false;
                }
                if (angleCount_[g] > 1)
                {
                    if (iter1.hasSingleValue() && !iter1.currentValuesSelected())
                    {
                        bOk = false;
                    }
                    if (iter2.hasSingleValue() && !iter2.currentValuesSelected())
                    {
                        bOk = false;
                    }
                }
                if (iter2.hasValue()
                    && (angleCount_[g] == 1 || (!iter1.hasSingleValue() && !iter2.hasSingleValue()))
                    && iter1.currentValuesSelected() != iter2.currentValuesSelected())
                {
                    bOk = false;
                }
                if (!bOk)
                {
                    std::string message = formatString(
                            "Dynamic selection %d does not select "
                            "a consistent set of angles over the frames",
                            static_cast<int>(g + 1));
                    GMX_THROW(InconsistentInputError(message));
                }
            }
        }
    }
}


void Angle::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& /* top */)
{
    initFromSelections(sel1_, sel2_);

    // checkSelections() ensures that both selection lists have the same size.
    angles_.setDataSetCount(angleCount_.size());
    for (size_t i = 0; i < angleCount_.size(); ++i)
    {
        angles_.setColumnCount(i, angleCount_[i]);
    }
    double histogramMin = (g1type_ == Group1Type::Dihedral ? -180.0 : 0);
    histogramModule_->init(histogramFromRange(histogramMin, 180.0).binWidth(binWidth_).includeAll());

    if (g2type_ == Group2Type::TimeZero)
    {
        vt0_.resize(sel1_.size());
        for (size_t g = 0; g < sel1_.size(); ++g)
        {
            vt0_[g].resize(sel1_[g].posCount() / natoms1_);
        }
    }

    if (!fnAverage_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAverage_);
        plotm->setTitle("Average angle");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Angle (degrees)");
        // TODO: Consider adding information about the second selection,
        // and/or a subtitle describing what kind of angle this is.
        for (size_t g = 0; g < sel1_.size(); ++g)
        {
            plotm->appendLegend(sel1_[g].name());
        }
        averageModule_->addModule(plotm);
    }

    if (!fnAll_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAll_);
        plotm->setTitle("Angle");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Angle (degrees)");
        // TODO: Add legends? (there can be a massive amount of columns)
        angles_.addModule(plotm);
    }

    if (!fnHistogram_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnHistogram_);
        plotm->setTitle("Angle histogram");
        plotm->setXLabel("Angle (degrees)");
        plotm->setYLabel("Probability");
        // TODO: Consider adding information about the second selection,
        // and/or a subtitle describing what kind of angle this is.
        for (size_t g = 0; g < sel1_.size(); ++g)
        {
            plotm->appendLegend(sel1_[g].name());
        }
        histogramModule_->averager().addModule(plotm);
    }
}


//! Helper method to calculate a vector from two or three positions.
void calc_vec(int natoms, rvec x[], t_pbc* pbc, rvec xout, rvec cout)
{
    switch (natoms)
    {
        case 2:
            if (pbc)
            {
                pbc_dx(pbc, x[1], x[0], xout);
            }
            else
            {
                rvec_sub(x[1], x[0], xout);
            }
            svmul(0.5, xout, cout);
            rvec_add(x[0], cout, cout);
            break;
        case 3:
        {
            rvec v1, v2;
            if (pbc)
            {
                pbc_dx(pbc, x[1], x[0], v1);
                pbc_dx(pbc, x[2], x[0], v2);
            }
            else
            {
                rvec_sub(x[1], x[0], v1);
                rvec_sub(x[2], x[0], v2);
            }
            cprod(v1, v2, xout);
            rvec_add(x[0], x[1], cout);
            rvec_add(cout, x[2], cout);
            svmul(1.0 / 3.0, cout, cout);
            break;
        }
        default: GMX_RELEASE_ASSERT(false, "Incorrectly initialized number of atoms");
    }
}


void Angle::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   dh   = pdata->dataHandle(angles_);
    const SelectionList& sel1 = TrajectoryAnalysisModuleData::parallelSelections(sel1_);
    const SelectionList& sel2 = TrajectoryAnalysisModuleData::parallelSelections(sel2_);

    checkSelections(sel1, sel2);

    dh.startFrame(frnr, fr.time);

    AnglePositionIterator iter1(sel1, natoms1_);
    AnglePositionIterator iter2(sel2, natoms2_);
    for (size_t g = 0; g < angleCount_.size(); ++g, iter1.nextGroup(), iter2.nextGroup())
    {
        rvec v1, v2;
        rvec c1, c2;

        // v2 & c2 are conditionally set in the switch statement below, and conditionally
        // used in a different switch statement later. Apparently the clang static analyzer
        // thinks there are cases where they can be used uninitialized (which I cannot find),
        // but to avoid trouble if we ever change just one of the switch statements it
        // makes sense to clear them outside the first switch.

        clear_rvec(v2);
        clear_rvec(c2);

        switch (g2type_)
        {
            case Group2Type::Z: v2[ZZ] = 1.0; break;
            case Group2Type::SphereNormal: copy_rvec(sel2_[g].position(0).x(), c2); break;
            default:
                // do nothing
                break;
        }

        dh.selectDataSet(g);
        for (int n = 0; n < angleCount_[g]; ++n, iter1.nextValue(), iter2.nextValue())
        {
            rvec x[4];
            // x[] will be assigned below based on the number of atoms used to initialize iter1,
            // which in turn should correspond perfectly to g1type_ (which determines how many we
            // read), but unsurprisingly the static analyzer chokes a bit on that.
            clear_rvecs(4, x);

            real angle = 0;
            // checkSelections() ensures that this reflects all the involved
            // positions.
            const bool bPresent = iter1.currentValuesSelected() && iter2.currentValuesSelected();
            iter1.getCurrentPositions(x);
            switch (g1type_)
            {
                case Group1Type::Angle:
                    if (pbc)
                    {
                        pbc_dx(pbc, x[0], x[1], v1);
                        pbc_dx(pbc, x[2], x[1], v2);
                    }
                    else
                    {
                        rvec_sub(x[0], x[1], v1);
                        rvec_sub(x[2], x[1], v2);
                    }
                    angle = gmx_angle(v1, v2);
                    break;
                case Group1Type::Dihedral:
                {
                    rvec dx[3];
                    if (pbc)
                    {
                        pbc_dx(pbc, x[0], x[1], dx[0]);
                        pbc_dx(pbc, x[2], x[1], dx[1]);
                        pbc_dx(pbc, x[2], x[3], dx[2]);
                    }
                    else
                    {
                        rvec_sub(x[0], x[1], dx[0]);
                        rvec_sub(x[2], x[1], dx[1]);
                        rvec_sub(x[2], x[3], dx[2]);
                    }
                    cprod(dx[0], dx[1], v1);
                    cprod(dx[1], dx[2], v2);
                    angle    = gmx_angle(v1, v2);
                    real ipr = iprod(dx[0], v2);
                    if (ipr < 0)
                    {
                        angle = -angle;
                    }
                    break;
                }
                case Group1Type::Vector:
                case Group1Type::Plane:
                    calc_vec(natoms1_, x, pbc, v1, c1);
                    switch (g2type_)
                    {
                        case Group2Type::Vector:
                        case Group2Type::Plane:
                            iter2.getCurrentPositions(x);
                            calc_vec(natoms2_, x, pbc, v2, c2);
                            break;
                        case Group2Type::TimeZero:
                            // FIXME: This is not parallelizable.
                            if (frnr == 0)
                            {
                                copy_rvec(v1, vt0_[g][n]);
                            }
                            copy_rvec(vt0_[g][n], v2);
                            break;
                        case Group2Type::Z: c1[XX] = c1[YY] = 0.0; break;
                        case Group2Type::SphereNormal:
                            if (pbc)
                            {
                                pbc_dx(pbc, c1, c2, v2);
                            }
                            else
                            {
                                rvec_sub(c1, c2, v2);
                            }
                            break;
                        default: GMX_THROW(InternalError("invalid -g2 value"));
                    }
                    angle = gmx_angle(v1, v2);
                    break;
                default: GMX_THROW(InternalError("invalid -g1 value"));
            }
            dh.setPoint(n, angle * gmx::c_rad2Deg, bPresent);
        }
    }
    dh.finishFrame();
}


void Angle::finishAnalysis(int /*nframes*/)
{
    AbstractAverageHistogram& averageHistogram = histogramModule_->averager();
    averageHistogram.normalizeProbability();
    averageHistogram.done();
}


void Angle::writeOutput() {}

} // namespace

const char AngleInfo::name[]             = "gangle";
const char AngleInfo::shortDescription[] = "Calculate angles";

TrajectoryAnalysisModulePointer AngleInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Angle);
}

} // namespace analysismodules

} // namespace gmx
