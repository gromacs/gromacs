/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
 * Implements gmx::analysismodules::Select.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "select.h"

#include <cstdio>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/lifetime.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionenums.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct t_pbc;

namespace gmx
{
class AbstractAnalysisData;

namespace analysismodules
{

namespace
{

/*! \internal \brief
 * Data module for writing index files.
 *
 * \ingroup module_analysisdata
 */
class IndexFileWriterModule : public AnalysisDataModuleSerial
{
public:
    IndexFileWriterModule();
    ~IndexFileWriterModule() override;

    //! Sets the file name to write the index file to.
    void setFileName(const std::string& fnm);
    /*! \brief
     * Adds information about a group to be printed.
     *
     * Must be called for each group present in the input data.
     */
    void addGroup(const std::string& name, bool bDynamic);

    int flags() const override;

    void dataStarted(AbstractAnalysisData* data) override;
    void frameStarted(const AnalysisDataFrameHeader& header) override;
    void pointsAdded(const AnalysisDataPointSetRef& points) override;
    void frameFinished(const AnalysisDataFrameHeader& header) override;
    void dataFinished() override;

private:
    void closeFile();

    struct GroupInfo
    {
        GroupInfo(const std::string& name, bool bDynamic) : name_(name), bDynamic_(bDynamic) {}

        std::string name_;
        bool        bDynamic_;
    };

    std::string            fnm_;
    std::vector<GroupInfo> groups_;
    FILE*                  fp_;
    int                    currentGroup_;
    int                    currentSize_;
    bool                   bAnyWritten_;
};

/********************************************************************
 * IndexFileWriterModule
 */

IndexFileWriterModule::IndexFileWriterModule() :
    fp_(nullptr), currentGroup_(-1), currentSize_(0), bAnyWritten_(false)
{
}


IndexFileWriterModule::~IndexFileWriterModule()
{
    closeFile();
}


void IndexFileWriterModule::closeFile()
{
    if (fp_ != nullptr)
    {
        gmx_fio_fclose(fp_);
        fp_ = nullptr;
    }
}


void IndexFileWriterModule::setFileName(const std::string& fnm)
{
    fnm_ = fnm;
}


void IndexFileWriterModule::addGroup(const std::string& name, bool bDynamic)
{
    std::string newName(name);
    std::replace(newName.begin(), newName.end(), ' ', '_');
    groups_.emplace_back(newName, bDynamic);
}


int IndexFileWriterModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void IndexFileWriterModule::dataStarted(AbstractAnalysisData* /*data*/)
{
    if (!fnm_.empty())
    {
        fp_ = gmx_fio_fopen(fnm_.c_str(), "w");
    }
}


void IndexFileWriterModule::frameStarted(const AnalysisDataFrameHeader& /*header*/)
{
    bAnyWritten_  = false;
    currentGroup_ = -1;
}


void IndexFileWriterModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (fp_ == nullptr)
    {
        return;
    }
    bool bFirstFrame = (points.frameIndex() == 0);
    if (points.firstColumn() == 0)
    {
        ++currentGroup_;
        GMX_RELEASE_ASSERT(currentGroup_ < gmx::ssize(groups_), "Too few groups initialized");
        if (bFirstFrame || groups_[currentGroup_].bDynamic_)
        {
            if (!bFirstFrame || currentGroup_ > 0)
            {
                std::fprintf(fp_, "\n\n");
            }
            std::string name = groups_[currentGroup_].name_;
            if (groups_[currentGroup_].bDynamic_)
            {
                name += formatString("_f%d_t%.3f", points.frameIndex(), points.x());
            }
            std::fprintf(fp_, "[ %s ]", name.c_str());
            bAnyWritten_ = true;
            currentSize_ = 0;
        }
    }
    else
    {
        if (bFirstFrame || groups_[currentGroup_].bDynamic_)
        {
            if (currentSize_ % 15 == 0)
            {
                std::fprintf(fp_, "\n");
            }
            std::fprintf(fp_, "%4d ", static_cast<int>(points.y(0)));
            ++currentSize_;
        }
    }
}


void IndexFileWriterModule::frameFinished(const AnalysisDataFrameHeader& /*header*/) {}


void IndexFileWriterModule::dataFinished()
{
    if (fp_ != nullptr)
    {
        std::fprintf(fp_, "\n");
    }
    closeFile();
}

/********************************************************************
 * Select
 */

//! How to identify residues in output files.
enum class ResidueNumbering : int
{
    ByNumber,
    ByIndex,
    Count
};
//! Which atoms to write out to PDB files.
enum class PdbAtomsSelection : int
{
    All,
    MaxSelection,
    Selected,
    Count
};
//! String values corresponding to ResidueNumbering.
const EnumerationArray<ResidueNumbering, const char*> c_residueNumberingTypeNames = { { "number",
                                                                                        "index" } };
//! String values corresponding to PdbAtomsSelection.
const EnumerationArray<PdbAtomsSelection, const char*> c_pdbAtomsTypeNames = {
    { "all", "maxsel", "selected" }
};

class Select : public TrajectoryAnalysisModule
{
public:
    Select();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;

    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    SelectionList sel_;

    std::string       fnSize_;
    std::string       fnFrac_;
    std::string       fnIndex_;
    std::string       fnNdx_;
    std::string       fnMask_;
    std::string       fnOccupancy_;
    std::string       fnPDB_;
    std::string       fnLifetime_;
    bool              bTotNorm_;
    bool              bFracNorm_;
    bool              bResInd_;
    bool              bCumulativeLifetimes_;
    ResidueNumbering  resNumberType_;
    PdbAtomsSelection pdbAtoms_;

    //! The input topology.
    const TopologyInformation*        top_;
    std::vector<int>                  totsize_;
    AnalysisData                      sdata_;
    AnalysisData                      cdata_;
    AnalysisData                      idata_;
    AnalysisData                      mdata_;
    AnalysisDataAverageModulePointer  occupancyModule_;
    AnalysisDataLifetimeModulePointer lifetimeModule_;
};

Select::Select() :
    bTotNorm_(false),
    bFracNorm_(false),
    bResInd_(false),
    bCumulativeLifetimes_(true),
    resNumberType_(ResidueNumbering::ByNumber),
    pdbAtoms_(PdbAtomsSelection::All),
    top_(nullptr),
    occupancyModule_(new AnalysisDataAverageModule()),
    lifetimeModule_(new AnalysisDataLifetimeModule())
{
    mdata_.addModule(occupancyModule_);
    mdata_.addModule(lifetimeModule_);

    registerAnalysisDataset(&sdata_, "size");
    registerAnalysisDataset(&cdata_, "cfrac");
    idata_.setColumnCount(0, 2);
    idata_.setMultipoint(true);
    registerAnalysisDataset(&idata_, "index");
    registerAnalysisDataset(&mdata_, "mask");
    occupancyModule_->setXAxis(1.0, 1.0);
    registerBasicDataset(occupancyModule_.get(), "occupancy");
    registerBasicDataset(lifetimeModule_.get(), "lifetime");
}


void Select::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] writes out basic data about dynamic selections.",
        "It can be used for some simple analyses, or the output can",
        "be combined with output from other programs and/or external",
        "analysis programs to calculate more complex things.",
        "For detailed help on the selection syntax, please use",
        "[TT]gmx help selections[tt].",
        "",
        "Any combination of the output options is possible, but note",
        "that [TT]-om[tt] only operates on the first selection.",
        "Also note that if you provide no output options, no output is",
        "produced.",
        "",
        "With [TT]-os[tt], calculates the number of positions in each",
        "selection for each frame. With [TT]-norm[tt], the output is",
        "between 0 and 1 and describes the fraction from the maximum",
        "number of positions (e.g., for selection 'resname RA and x < 5'",
        "the maximum number of positions is the number of atoms in",
        "RA residues). With [TT]-cfnorm[tt], the output is divided",
        "by the fraction covered by the selection.",
        "[TT]-norm[tt] and [TT]-cfnorm[tt] can be specified independently",
        "of one another.",
        "",
        "With [TT]-oc[tt], the fraction covered by each selection is",
        "written out as a function of time.",
        "",
        "With [TT]-oi[tt], the selected atoms/residues/molecules are",
        "written out as a function of time. In the output, the first",
        "column contains the frame time, the second contains the number",
        "of positions, followed by the atom/residue/molecule numbers.",
        "If more than one selection is specified, the size of the second",
        "group immediately follows the last number of the first group",
        "and so on.",
        "",
        "With [TT]-on[tt], the selected atoms are written as a index file",
        "compatible with [TT]make_ndx[tt] and the analyzing tools. Each selection",
        "is written as a selection group and for dynamic selections a",
        "group is written for each frame.",
        "",
        "For residue numbers, the output of [TT]-oi[tt] can be controlled",
        "with [TT]-resnr[tt]: [TT]number[tt] (default) prints the residue",
        "numbers as they appear in the input file, while [TT]index[tt] prints",
        "unique numbers assigned to the residues in the order they appear",
        "in the input file, starting with 1. The former is more intuitive,",
        "but if the input contains multiple residues with the same number,",
        "the output can be less useful.",
        "",
        "With [TT]-om[tt], a mask is printed for the first selection",
        "as a function of time. Each line in the output corresponds to",
        "one frame, and contains either 0/1 for each atom/residue/molecule",
        "possibly selected. 1 stands for the atom/residue/molecule being",
        "selected for the current frame, 0 for not selected.",
        "",
        "With [TT]-of[tt], the occupancy fraction of each position (i.e.,",
        "the fraction of frames where the position is selected) is",
        "printed.",
        "",
        "With [TT]-ofpdb[tt], a PDB file is written out where the occupancy",
        "column is filled with the occupancy fraction of each atom in the",
        "selection. The coordinates in the PDB file will be those from the",
        "input topology. [TT]-pdbatoms[tt] can be used to control which atoms",
        "appear in the output PDB file: with [TT]all[tt] all atoms are",
        "present, with [TT]maxsel[tt] all atoms possibly selected by the",
        "selection are present, and with [TT]selected[tt] only atoms that are",
        "selected at least in one frame are present.",
        "",
        "With [TT]-olt[tt], a histogram is produced that shows the number of",
        "selected positions as a function of the time the position was",
        "continuously selected. [TT]-cumlt[tt] can be used to control whether",
        "subintervals of longer intervals are included in the histogram.[PAR]",
        "[TT]-om[tt], [TT]-of[tt], and [TT]-olt[tt] only make sense with",
        "dynamic selections.",
        "",
        "To plot coordinates for selections, use [gmx-trajectory]."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("os")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnSize_)
                               .defaultBasename("size")
                               .description("Number of positions in each selection"));
    options->addOption(FileNameOption("oc")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnFrac_)
                               .defaultBasename("cfrac")
                               .description("Covered fraction for each selection"));
    options->addOption(FileNameOption("oi")
                               .filetype(OptionFileType::GenericData)
                               .outputFile()
                               .store(&fnIndex_)
                               .defaultBasename("index")
                               .description("Indices selected by each selection"));
    options->addOption(FileNameOption("on")
                               .filetype(OptionFileType::AtomIndex)
                               .outputFile()
                               .store(&fnNdx_)
                               .defaultBasename("index")
                               .description("Index file from the selection"));
    options->addOption(FileNameOption("om")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnMask_)
                               .defaultBasename("mask")
                               .description("Mask for selected positions"));
    options->addOption(FileNameOption("of")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnOccupancy_)
                               .defaultBasename("occupancy")
                               .description("Occupied fraction for selected positions"));
    options->addOption(
            FileNameOption("ofpdb")
                    .filetype(OptionFileType::PDB)
                    .outputFile()
                    .store(&fnPDB_)
                    .defaultBasename("occupancy")
                    .description("PDB file with occupied fraction for selected positions"));
    options->addOption(FileNameOption("olt")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnLifetime_)
                               .defaultBasename("lifetime")
                               .description("Lifetime histogram"));

    options->addOption(SelectionOption("select").storeVector(&sel_).required().multiValue().description(
            "Selections to analyze"));

    options->addOption(
            BooleanOption("norm").store(&bTotNorm_).description("Normalize by total number of positions with -os"));
    options->addOption(
            BooleanOption("cfnorm").store(&bFracNorm_).description("Normalize by covered fraction with -os"));
    options->addOption(EnumOption<ResidueNumbering>("resnr")
                               .store(&resNumberType_)
                               .enumValue(c_residueNumberingTypeNames)
                               .description("Residue number output type with -oi and -on"));
    options->addOption(EnumOption<PdbAtomsSelection>("pdbatoms")
                               .store(&pdbAtoms_)
                               .enumValue(c_pdbAtomsTypeNames)
                               .description("Atoms to write with -ofpdb"));
    options->addOption(BooleanOption("cumlt")
                               .store(&bCumulativeLifetimes_)
                               .description("Cumulate subintervals of longer intervals in -olt"));
}

void Select::optionsFinished(TrajectoryAnalysisSettings* settings)
{
    if (!fnPDB_.empty())
    {
        settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
        settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
    }
}

void Select::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    bResInd_ = (resNumberType_ == ResidueNumbering::ByIndex);

    for (auto& i : sel_)
    {
        i.initCoveredFraction(CFRAC_SOLIDANGLE);
    }

    // TODO: For large systems, a float may not have enough precision
    sdata_.setColumnCount(0, sel_.size());
    totsize_.reserve(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        totsize_.push_back(sel_[g].posCount());
    }
    if (!fnSize_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnSize_);
        plot->setTitle("Selection size");
        plot->setXAxisIsTime();
        plot->setYLabel("Number");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plot->appendLegend(sel_[g].name());
        }
        sdata_.addModule(plot);
    }

    cdata_.setColumnCount(0, sel_.size());
    if (!fnFrac_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnFrac_);
        plot->setTitle("Covered fraction");
        plot->setXAxisIsTime();
        plot->setYLabel("Fraction");
        plot->setYFormat(6, 4);
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plot->appendLegend(sel_[g].name());
        }
        cdata_.addModule(plot);
    }

    // TODO: For large systems, a float may not have enough precision
    if (!fnIndex_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnIndex_);
        plot->setPlainOutput(true);
        plot->setYFormat(4, 0);
        idata_.addModule(plot);
    }
    if (!fnNdx_.empty())
    {
        std::shared_ptr<IndexFileWriterModule> writer(new IndexFileWriterModule());
        writer->setFileName(fnNdx_);
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            writer->addGroup(sel_[g].name(), sel_[g].isDynamic());
        }
        idata_.addModule(writer);
    }

    mdata_.setDataSetCount(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        mdata_.setColumnCount(g, sel_[g].posCount());
    }
    lifetimeModule_->setCumulative(bCumulativeLifetimes_);
    if (!fnMask_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnMask_);
        plot->setTitle("Selection mask");
        plot->setXAxisIsTime();
        plot->setYLabel("Occupancy");
        plot->setYFormat(1, 0);
        // TODO: Add legend? (there can be massive amount of columns)
        mdata_.addModule(plot);
    }
    if (!fnOccupancy_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnOccupancy_);
        plot->setTitle("Fraction of time selection matches");
        plot->setXLabel("Selected position");
        plot->setYLabel("Occupied fraction");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plot->appendLegend(sel_[g].name());
        }
        occupancyModule_->addModule(plot);
    }
    if (!fnLifetime_.empty())
    {
        AnalysisDataPlotModulePointer plot(new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnLifetime_);
        plot->setTitle("Lifetime histogram");
        plot->setXAxisIsTime();
        plot->setYLabel("Number of occurrences");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plot->appendLegend(sel_[g].name());
        }
        lifetimeModule_->addModule(plot);
    }

    top_ = &top;
}


void Select::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* /* pbc */, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   sdh = pdata->dataHandle(sdata_);
    AnalysisDataHandle   cdh = pdata->dataHandle(cdata_);
    AnalysisDataHandle   idh = pdata->dataHandle(idata_);
    AnalysisDataHandle   mdh = pdata->dataHandle(mdata_);
    const SelectionList& sel = TrajectoryAnalysisModuleData::parallelSelections(sel_);

    sdh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        real normfac = bFracNorm_ ? 1.0 / sel[g].coveredFraction() : 1.0;
        if (bTotNorm_)
        {
            normfac /= totsize_[g];
        }
        sdh.setPoint(g, sel[g].posCount() * normfac);
    }
    sdh.finishFrame();

    cdh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        cdh.setPoint(g, sel[g].coveredFraction());
    }
    cdh.finishFrame();

    idh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        idh.setPoint(0, sel[g].posCount());
        idh.finishPointSet();
        for (int i = 0; i < sel[g].posCount(); ++i)
        {
            const SelectionPosition& p = sel[g].position(i);
            if (sel[g].type() == INDEX_RES && !bResInd_)
            {
                idh.setPoint(1, top_->atoms()->resinfo[p.mappedId()].nr);
            }
            else
            {
                idh.setPoint(1, p.mappedId() + 1);
            }
            idh.finishPointSet();
        }
    }
    idh.finishFrame();

    mdh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        mdh.selectDataSet(g);
        for (int i = 0; i < totsize_[g]; ++i)
        {
            mdh.setPoint(i, 0);
        }
        for (int i = 0; i < sel[g].posCount(); ++i)
        {
            mdh.setPoint(sel[g].position(i).refId(), 1);
        }
    }
    mdh.finishFrame();
}


void Select::finishAnalysis(int /*nframes*/) {}


void Select::writeOutput()
{
    if (!fnPDB_.empty())
    {
        GMX_RELEASE_ASSERT(top_->hasTopology(),
                           "Topology should have been loaded or an error given earlier");
        auto atoms = top_->copyAtoms();
        if (!atoms->havePdbInfo)
        {
            snew(atoms->pdbinfo, atoms->nr);
            atoms->havePdbInfo = TRUE;
        }
        for (int i = 0; i < atoms->nr; ++i)
        {
            atoms->pdbinfo[i].occup = 0.0;
        }
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            for (int i = 0; i < sel_[g].posCount(); ++i)
            {
                ArrayRef<const int>                 atomIndices = sel_[g].position(i).atomIndices();
                ArrayRef<const int>::const_iterator ai;
                for (ai = atomIndices.begin(); ai != atomIndices.end(); ++ai)
                {
                    atoms->pdbinfo[*ai].occup += occupancyModule_->average(g, i);
                }
            }
        }

        std::vector<RVec> x = copyOf(top_->x());
        t_trxframe        fr;
        clear_trxframe(&fr, TRUE);
        fr.bAtoms = TRUE;
        fr.atoms  = atoms.get();
        fr.bX     = TRUE;
        fr.bBox   = TRUE;
        fr.x      = as_rvec_array(x.data());
        top_->getBox(fr.box);

        switch (pdbAtoms_)
        {
            case PdbAtomsSelection::All:
            {
                t_trxstatus* status = open_trx(fnPDB_.c_str(), "w");
                write_trxframe(status, &fr, nullptr);
                close_trx(status);
                break;
            }
            case PdbAtomsSelection::MaxSelection:
            {
                std::set<int> atomIndicesSet;
                for (size_t g = 0; g < sel_.size(); ++g)
                {
                    ArrayRef<const int> atomIndices = sel_[g].atomIndices();
                    atomIndicesSet.insert(atomIndices.begin(), atomIndices.end());
                }
                std::vector<int> allAtomIndices(atomIndicesSet.begin(), atomIndicesSet.end());
                t_trxstatus*     status = open_trx(fnPDB_.c_str(), "w");
                write_trxframe_indexed(status, &fr, allAtomIndices.size(), allAtomIndices.data(), nullptr);
                close_trx(status);
                break;
            }
            case PdbAtomsSelection::Selected:
            {
                std::vector<int> indices;
                for (int i = 0; i < atoms->nr; ++i)
                {
                    if (atoms->pdbinfo[i].occup > 0.0)
                    {
                        indices.push_back(i);
                    }
                }
                t_trxstatus* status = open_trx(fnPDB_.c_str(), "w");
                write_trxframe_indexed(status, &fr, indices.size(), indices.data(), nullptr);
                close_trx(status);
                break;
            }
            case PdbAtomsSelection::Count:
                GMX_RELEASE_ASSERT(false,
                                   "Mismatch between -pdbatoms enum values and implementation");
        }
    }
}

} // namespace

const char SelectInfo::name[]             = "select";
const char SelectInfo::shortDescription[] = "Print general information about selections";

TrajectoryAnalysisModulePointer SelectInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Select);
}

} // namespace analysismodules

} // namespace gmx
