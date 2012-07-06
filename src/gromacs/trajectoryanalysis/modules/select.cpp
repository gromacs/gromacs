/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Select.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_trajectoryanalysis
 */
#include "select.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

#include "gmxfio.h"

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
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

namespace
{

/*! \internal \brief
 * Data module for writing index files.
 *
 * \ingroup module_analysisdata
 */
class IndexFileWriterModule : public AnalysisDataModuleInterface
{
    public:
        IndexFileWriterModule();
        virtual ~IndexFileWriterModule();

        //! Sets the file name to write the index file to.
        void setFileName(const std::string &fnm);
        void addGroup(const std::string &name, bool bDynamic);

        virtual int flags() const;

        virtual void dataStarted(AbstractAnalysisData *data);
        virtual void frameStarted(const AnalysisDataFrameHeader &header);
        virtual void pointsAdded(const AnalysisDataPointSetRef &points);
        virtual void frameFinished(const AnalysisDataFrameHeader &header);
        virtual void dataFinished();

    private:
        void closeFile();

        struct GroupInfo
        {
            GroupInfo(const std::string &name, bool bDynamic)
                : name(name), bDynamic(bDynamic)
            { }

            std::string         name;
            bool                bDynamic;
        };

        std::string             fnm_;
        std::vector<GroupInfo>  groups_;
        FILE                   *fp_;
        int                     currentGroup_;
        int                     currentSize_;
        bool                    bAnyWritten_;
};

/********************************************************************
 * IndexFileWriterModule
 */

IndexFileWriterModule::IndexFileWriterModule()
    : fp_(NULL), currentGroup_(-1), currentSize_(0), bAnyWritten_(false)
{
}


IndexFileWriterModule::~IndexFileWriterModule()
{
    closeFile();
}


void IndexFileWriterModule::closeFile()
{
    if (fp_ != NULL)
    {
        gmx_fio_fclose(fp_);
        fp_ = NULL;
    }
}


void IndexFileWriterModule::setFileName(const std::string &fnm)
{
    fnm_ = fnm;
}


void IndexFileWriterModule::addGroup(const std::string &name, bool bDynamic)
{
    std::string newName(name);
    std::replace(newName.begin(), newName.end(), ' ', '_');
    groups_.push_back(GroupInfo(newName, bDynamic));
}


int IndexFileWriterModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void IndexFileWriterModule::dataStarted(AbstractAnalysisData * /*data*/)
{
    if (!fnm_.empty())
    {
        fp_ = gmx_fio_fopen(fnm_.c_str(), "w");
    }
}


void IndexFileWriterModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
    bAnyWritten_ = false;
    currentGroup_ = -1;
}


void
IndexFileWriterModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (fp_ == NULL)
    {
        return;
    }
    bool bFirstFrame = (points.frameIndex() == 0);
    if (points.firstColumn() == 0)
    {
        ++currentGroup_;
        if (bFirstFrame || groups_[currentGroup_].bDynamic)
        {
            if (!bFirstFrame || currentGroup_ > 0)
            {
                std::fprintf(fp_, "\n\n");
            }
            std::string name = groups_[currentGroup_].name;
            if (groups_[currentGroup_].bDynamic)
            {
                name += formatString("f_%d_t%.3f", points.frameIndex(), points.x());
            }
            std::fprintf(fp_, "[ %s ]", name.c_str());
            bAnyWritten_ = true;
            currentSize_ = 0;
        }
    }
    else
    {
        if (bFirstFrame || groups_[currentGroup_].bDynamic)
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


void IndexFileWriterModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
}


void IndexFileWriterModule::dataFinished()
{
    if (fp_ != NULL)
    {
        std::fprintf(fp_, "\n");
    }
    closeFile();
}

} // namespace


/********************************************************************
 * Select
 */

const char Select::name[] = "select";
const char Select::shortDescription[] =
    "Print general information about selections";

Select::Select()
    : options_(name, shortDescription),
      bDump_(false), bTotNorm_(false), bFracNorm_(false), bResInd_(false),
      top_(NULL)
{
    registerAnalysisDataset(&sdata_, "size");
    registerAnalysisDataset(&cdata_, "cfrac");
    idata_.setColumnCount(2);
    idata_.setMultipoint(true);
    registerAnalysisDataset(&idata_, "index");
    registerAnalysisDataset(&mdata_, "mask");
}


Select::~Select()
{
}


Options &
Select::initOptions(TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[TT]g_select[tt] writes out basic data about dynamic selections.",
        "It can be used for some simple analyses, or the output can",
        "be combined with output from other programs and/or external",
        "analysis programs to calculate more complex things.",
        "Any combination of the output options is possible, but note",
        "that [TT]-om[tt] only operates on the first selection.",
        "Also note that if you provide no output options, no output is",
        "produced.[PAR]",
        "With [TT]-os[tt], calculates the number of positions in each",
        "selection for each frame. With [TT]-norm[tt], the output is",
        "between 0 and 1 and describes the fraction from the maximum",
        "number of positions (e.g., for selection 'resname RA and x < 5'",
        "the maximum number of positions is the number of atoms in",
        "RA residues). With [TT]-cfnorm[tt], the output is divided",
        "by the fraction covered by the selection.",
        "[TT]-norm[tt] and [TT]-cfnorm[tt] can be specified independently",
        "of one another.[PAR]",
        "With [TT]-oc[tt], the fraction covered by each selection is",
        "written out as a function of time.[PAR]",
        "With [TT]-oi[tt], the selected atoms/residues/molecules are",
        "written out as a function of time. In the output, the first",
        "column contains the frame time, the second contains the number",
        "of positions, followed by the atom/residue/molecule numbers.",
        "If more than one selection is specified, the size of the second",
        "group immediately follows the last number of the first group",
        "and so on. With [TT]-dump[tt], the frame time and the number",
        "of positions is omitted from the output. In this case, only one",
        "selection can be given.[PAR]",
        "With [TT]-on[tt], the selected atoms are written as a index file",
        "compatible with [TT]make_ndx[tt] and the analyzing tools. Each selection",
        "is written as a selection group and for dynamic selections a",
        "group is written for each frame.[PAR]",
        "For residue numbers, the output of [TT]-oi[tt] can be controlled",
        "with [TT]-resnr[tt]: [TT]number[tt] (default) prints the residue",
        "numbers as they appear in the input file, while [TT]index[tt] prints",
        "unique numbers assigned to the residues in the order they appear",
        "in the input file, starting with 1. The former is more intuitive,",
        "but if the input contains multiple residues with the same number,",
        "the output can be less useful.[PAR]",
        "With [TT]-om[tt], a mask is printed for the first selection",
        "as a function of time. Each line in the output corresponds to",
        "one frame, and contains either 0/1 for each atom/residue/molecule",
        "possibly selected. 1 stands for the atom/residue/molecule being",
        "selected for the current frame, 0 for not selected.",
        "With [TT]-dump[tt], the frame time is omitted from the output."
    };

    options_.setDescription(concatenateStrings(desc));

    options_.addOption(FileNameOption("os").filetype(eftPlot).outputFile()
                           .store(&fnSize_).defaultValueIfSet("size"));
    options_.addOption(FileNameOption("oc").filetype(eftPlot).outputFile()
                           .store(&fnFrac_).defaultValueIfSet("frac"));
    options_.addOption(FileNameOption("oi").filetype(eftGenericData).outputFile()
                           .store(&fnIndex_).defaultValueIfSet("index"));
    options_.addOption(FileNameOption("on").filetype(eftIndex).outputFile()
                           .store(&fnNdx_).defaultValueIfSet("index"));
    options_.addOption(FileNameOption("om").filetype(eftPlot).outputFile()
                           .store(&fnMask_).defaultValueIfSet("mask"));

    options_.addOption(SelectionOption("select").storeVector(&sel_)
        .required().multiValue()
        .description("Selections to analyze"));

    options_.addOption(BooleanOption("dump").store(&bDump_)
        .description("Do not print the frame time (-om, -oi) or the index size (-oi)"));
    options_.addOption(BooleanOption("norm").store(&bTotNorm_)
        .description("Normalize by total number of positions with -os"));
    options_.addOption(BooleanOption("cfnorm").store(&bFracNorm_)
        .description("Normalize by covered fraction with -os"));
    const char *const cResNumberEnum[] = { "number", "index", NULL };
    options_.addOption(StringOption("resnr").store(&resNumberType_)
        .enumValue(cResNumberEnum).defaultEnumIndex(0)
        .description("Residue number output type"));

    return options_;
}


void
Select::initAnalysis(const TrajectoryAnalysisSettings &settings,
                     const TopologyInformation &top)
{
    if (!fnIndex_.empty() && bDump_ && sel_.size() > 1U)
    {
        GMX_THROW(InconsistentInputError("With -oi and -dump, there can be only one selection"));
    }
    bResInd_ = (resNumberType_ == "index");

    for (SelectionList::iterator i = sel_.begin(); i != sel_.end(); ++i)
    {
        i->initCoveredFraction(CFRAC_SOLIDANGLE);
    }

    // TODO: For large systems, a float may not have enough precision
    sdata_.setColumnCount(sel_.size());
    totsize_.reserve(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        totsize_.push_back(sel_[g].posCount());
    }
    if (!fnSize_.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnSize_);
        plot->setTitle("Selection size");
        plot->setXAxisIsTime();
        plot->setYLabel("Number");
        sdata_.addModule(plot);
    }

    cdata_.setColumnCount(sel_.size());
    if (!fnFrac_.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnFrac_);
        plot->setTitle("Covered fraction");
        plot->setXAxisIsTime();
        plot->setYLabel("Fraction");
        plot->setYFormat(6, 4);
        cdata_.addModule(plot);
    }

    // TODO: For large systems, a float may not have enough precision
    if (!fnIndex_.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(fnIndex_);
        plot->setPlainOutput(true);
        plot->setYFormat(4, 0);
        if (bDump_)
        {
            plot->setOmitX(bDump_);
            idata_.addColumnModule(1, 1, plot);
        }
        else
        {
            idata_.addModule(plot);
        }
    }
    if (!fnNdx_.empty())
    {
        boost::shared_ptr<IndexFileWriterModule> writer(new IndexFileWriterModule());
        writer->setFileName(fnNdx_);
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            writer->addGroup(sel_[g].name(), sel_[g].isDynamic());
        }
        idata_.addModule(writer);
    }

    mdata_.setColumnCount(sel_[0].posCount());
    if (!fnMask_.empty())
    {
        if (sel_.size() > 1U)
        {
            fprintf(stderr, "WARNING: the mask (-om) will only be written for the first group\n");
        }
        if (!sel_[0].isDynamic())
        {
            fprintf(stderr, "WARNING: will not write the mask (-om) for a static selection\n");
        }
        else
        {
            AnalysisDataPlotModulePointer plot(
                new AnalysisDataPlotModule(settings.plotSettings()));
            plot->setFileName(fnMask_);
            plot->setPlainOutput(bDump_);
            plot->setOmitX(bDump_);
            plot->setTitle("Selection mask");
            plot->setXAxisIsTime();
            plot->setYLabel("Occupancy");
            plot->setYFormat(1, 0);
            mdata_.addModule(plot);
        }
    }

    top_ = top.topology();
}


void
Select::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                     TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle sdh = pdata->dataHandle(sdata_);
    AnalysisDataHandle cdh = pdata->dataHandle(cdata_);
    AnalysisDataHandle idh = pdata->dataHandle(idata_);
    AnalysisDataHandle mdh = pdata->dataHandle(mdata_);
    const SelectionList &sel = pdata->parallelSelections(sel_);

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
            const SelectionPosition &p = sel[g].position(i);
            if (sel[g].type() == INDEX_RES && !bResInd_)
            {
                idh.setPoint(1, top_->atoms.resinfo[p.mappedId()].nr);
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
    for (int i = 0; i < totsize_[0]; ++i)
    {
        mdh.setPoint(i, 0);
    }
    for (int i = 0; i < sel[0].posCount(); ++i)
    {
        mdh.setPoint(sel[0].position(i).refId(), 1);
    }
    mdh.finishFrame();
}


void
Select::finishAnalysis(int /*nframes*/)
{
}


void
Select::writeOutput()
{
}

} // namespace analysismodules

} // namespace gmx
