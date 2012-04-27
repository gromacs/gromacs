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
#include "gromacs/utility/format.h"

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

        std::string             _fnm;
        std::vector<GroupInfo>  _groups;
        FILE                   *_fp;
        int                     _currentGroup;
        int                     _currentSize;
        bool                    _bAnyWritten;
};

/********************************************************************
 * IndexFileWriterModule
 */

IndexFileWriterModule::IndexFileWriterModule()
    : _fp(NULL), _currentGroup(-1), _currentSize(0), _bAnyWritten(false)
{
}


IndexFileWriterModule::~IndexFileWriterModule()
{
    closeFile();
}


void IndexFileWriterModule::closeFile()
{
    if (_fp != NULL)
    {
        gmx_fio_fclose(_fp);
        _fp = NULL;
    }
}


void IndexFileWriterModule::setFileName(const std::string &fnm)
{
    _fnm = fnm;
}


void IndexFileWriterModule::addGroup(const std::string &name, bool bDynamic)
{
    std::string newName(name);
    std::replace(newName.begin(), newName.end(), ' ', '_');
    _groups.push_back(GroupInfo(newName, bDynamic));
}


int IndexFileWriterModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void IndexFileWriterModule::dataStarted(AbstractAnalysisData * /*data*/)
{
    if (!_fnm.empty())
    {
        _fp = gmx_fio_fopen(_fnm.c_str(), "w");
    }
}


void IndexFileWriterModule::frameStarted(const AnalysisDataFrameHeader & /*header*/)
{
    _bAnyWritten = false;
    _currentGroup = -1;
}


void
IndexFileWriterModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (_fp == NULL)
    {
        return;
    }
    bool bFirstFrame = (points.frameIndex() == 0);
    if (points.firstColumn() == 0)
    {
        ++_currentGroup;
        if (bFirstFrame || _groups[_currentGroup].bDynamic)
        {
            if (!bFirstFrame || _currentGroup > 0)
            {
                std::fprintf(_fp, "\n\n");
            }
            std::string name = _groups[_currentGroup].name;
            if (_groups[_currentGroup].bDynamic)
            {
                name += formatString("_f%d_t%.3f", points.frameIndex(), points.x());
            }
            std::fprintf(_fp, "[ %s ]", name.c_str());
            _bAnyWritten = true;
            _currentSize = 0;
        }
    }
    else
    {
        if (bFirstFrame || _groups[_currentGroup].bDynamic)
        {
            if (_currentSize % 15 == 0)
            {
                std::fprintf(_fp, "\n");
            }
            std::fprintf(_fp, "%4d ", static_cast<int>(points.y(0)));
            ++_currentSize;
        }
    }
}


void IndexFileWriterModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
}


void IndexFileWriterModule::dataFinished()
{
    if (_fp != NULL)
    {
        std::fprintf(_fp, "\n");
    }
    closeFile();
}

} // namespace


/********************************************************************
 * Select
 */

Select::Select()
    : _options("select", "Selection information"),
      _bDump(false), _bTotNorm(false), _bFracNorm(false), _bResInd(false),
      _top(NULL)
{
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
        "that [TT]-om[tt] only operates on the first selection.[PAR]",
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
        "With [TT]-dump[tt], the frame time is omitted from the output.",
        NULL
    };

    _options.setDescription(desc);

    _options.addOption(FileNameOption("os").filetype(eftPlot).outputFile()
                           .store(&_fnSize).defaultValueIfSet("size"));
    _options.addOption(FileNameOption("oc").filetype(eftPlot).outputFile()
                           .store(&_fnFrac).defaultValueIfSet("frac"));
    _options.addOption(FileNameOption("oi").filetype(eftGenericData).outputFile()
                           .store(&_fnIndex).defaultValueIfSet("index"));
    _options.addOption(FileNameOption("on").filetype(eftIndex).outputFile()
                           .store(&_fnNdx).defaultValueIfSet("index"));
    _options.addOption(FileNameOption("om").filetype(eftPlot).outputFile()
                           .store(&_fnMask).defaultValueIfSet("mask"));

    _options.addOption(SelectionOption("select").required().multiValue()
                           .storeVector(&_sel));

    _options.addOption(BooleanOption("dump").store(&_bDump)
        .description("Do not print the frame time (-om, -oi) or the index size (-oi)"));
    _options.addOption(BooleanOption("norm").store(&_bTotNorm)
        .description("Normalize by total number of positions with -os"));
    _options.addOption(BooleanOption("cfnorm").store(&_bFracNorm)
        .description("Normalize by covered fraction with -os"));
    const char *const cResNumberEnum[] = { "number", "index", NULL };
    _options.addOption(StringOption("resnr").store(&_resNumberType)
        .enumValue(cResNumberEnum).defaultEnumIndex(0)
        .description("Residue number output type"));

    return _options;
}


void
Select::initAnalysis(const TrajectoryAnalysisSettings &settings,
                     const TopologyInformation &top)
{
    if (!_fnIndex.empty() && _bDump && _sel.size() > 1U)
    {
        GMX_THROW(InconsistentInputError("With -oi and -dump, there can be only one selection"));
    }
    _bResInd = (_resNumberType == "index");

    for (SelectionList::iterator i = _sel.begin(); i != _sel.end(); ++i)
    {
        i->initCoveredFraction(CFRAC_SOLIDANGLE);
    }

    // TODO: For large systems, a float may not have enough precision
    _sdata.setColumnCount(_sel.size());
    registerAnalysisDataset(&_sdata, "size");
    _totsize.reserve(_sel.size());
    for (size_t g = 0; g < _sel.size(); ++g)
    {
        _totsize.push_back(_sel[g].posCount());
    }
    if (!_fnSize.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(_fnSize);
        plot->setTitle("Selection size");
        plot->setXAxisIsTime();
        plot->setYLabel("Number");
        _sdata.addModule(plot);
    }

    _cdata.setColumnCount(_sel.size());
    registerAnalysisDataset(&_cdata, "cfrac");
    if (!_fnFrac.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(_fnFrac);
        plot->setTitle("Covered fraction");
        plot->setXAxisIsTime();
        plot->setYLabel("Fraction");
        plot->setYFormat(6, 4);
        _cdata.addModule(plot);
    }

    // TODO: For large systems, a float may not have enough precision
    _idata.setColumnCount(2);
    _idata.setMultipoint(true);
    registerAnalysisDataset(&_idata, "index");
    if (!_fnIndex.empty())
    {
        AnalysisDataPlotModulePointer plot(
            new AnalysisDataPlotModule(settings.plotSettings()));
        plot->setFileName(_fnIndex);
        plot->setPlainOutput(true);
        plot->setYFormat(4, 0);
        if (_bDump)
        {
            plot->setOmitX(_bDump);
            _idata.addColumnModule(1, 1, plot);
        }
        else
        {
            _idata.addModule(plot);
        }
    }
    if (!_fnNdx.empty())
    {
        boost::shared_ptr<IndexFileWriterModule> writer(new IndexFileWriterModule());
        writer->setFileName(_fnNdx);
        for (size_t g = 0; g < _sel.size(); ++g)
        {
            writer->addGroup(_sel[g].name(), _sel[g].isDynamic());
        }
        _idata.addModule(writer);
    }

    _mdata.setColumnCount(_sel[0].posCount());
    registerAnalysisDataset(&_mdata, "mask");
    if (!_fnMask.empty())
    {
        if (_sel.size() > 1U)
        {
            fprintf(stderr, "WARNING: the mask (-om) will only be written for the first group\n");
        }
        if (!_sel[0].isDynamic())
        {
            fprintf(stderr, "WARNING: will not write the mask (-om) for a static selection\n");
        }
        else
        {
            AnalysisDataPlotModulePointer plot(
                new AnalysisDataPlotModule(settings.plotSettings()));
            plot->setFileName(_fnMask);
            plot->setPlainOutput(_bDump);
            plot->setOmitX(_bDump);
            plot->setTitle("Selection mask");
            plot->setXAxisIsTime();
            plot->setYLabel("Occupancy");
            plot->setYFormat(1, 0);
            _mdata.addModule(plot);
        }
    }

    _top = top.topology();
}


void
Select::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                     TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle sdh = pdata->dataHandle(_sdata);
    AnalysisDataHandle cdh = pdata->dataHandle(_cdata);
    AnalysisDataHandle idh = pdata->dataHandle(_idata);
    AnalysisDataHandle mdh = pdata->dataHandle(_mdata);
    const SelectionList &sel = pdata->parallelSelections(_sel);

    sdh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        real normfac = _bFracNorm ? 1.0 / sel[g].coveredFraction() : 1.0;
        if (_bTotNorm)
        {
            normfac /= _totsize[g];
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
            if (sel[g].type() == INDEX_RES && !_bResInd)
            {
                idh.setPoint(1, _top->atoms.resinfo[p.mappedId()].nr);
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
    int nextRefId = sel[0].position(0).refId();
    for (int b = 0, i = 0; b < _totsize[0]; ++b)
    {
        mdh.setPoint(b, b == nextRefId ? 1 : 0);
        if (b == nextRefId)
        {
            ++i;
            nextRefId = sel[0].position(i).refId();
        }
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


TrajectoryAnalysisModulePointer
Select::create()
{
    return TrajectoryAnalysisModulePointer(new Select());
}

} // namespace analysismodules

} // namespace gmx
