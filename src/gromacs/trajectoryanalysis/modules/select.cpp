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
#include "select.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <index.h>
#include <smalloc.h>
#include <string2.h>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fatalerror/fatalerror.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"

namespace gmx
{

namespace analysismodules
{

class Select::ModuleData : public TrajectoryAnalysisModuleData
{
    public:
        ModuleData() : _mmap(NULL)
        {
        }

        virtual ~ModuleData()
        {
            if (_mmap)
            {
                gmx_ana_indexmap_deinit(_mmap);
                sfree(_mmap);
            }
        }

        virtual int finish()
        {
            return finishDataHandles();
        }

        gmx_ana_indexmap_t  *_mmap;
};


Select::Select()
    : _options("select", "Selection information"),
      _bTotNorm(false), _bFracNorm(false), _bResInd(false),
      _block(NULL), _gnames(NULL)
{
}


Select::~Select()
{
    if (_block != NULL)
    {
        done_blocka(_block);
        sfree(_block);
    }
}


Options *
Select::initOptions(TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "g_select",
        NULL
    };

    _options.setDescription(desc);

    _options.addOption(FileNameOption("os").filetype(eftPlot).writeOnly()
                           .store(&_fnSize).defaultValueIfSet("size"));
    _options.addOption(FileNameOption("oc").filetype(eftPlot).writeOnly()
                           .store(&_fnFrac).defaultValueIfSet("frac"));
    _options.addOption(FileNameOption("oi").filetype(eftPlot).writeOnly()
                           .store(&_fnIndex).defaultValueIfSet("index"));
    _options.addOption(FileNameOption("on").filetype(eftIndex).writeOnly()
                           .store(&_fnNdx).defaultValueIfSet("index"));
    _options.addOption(FileNameOption("om").filetype(eftPlot).writeOnly()
                           .store(&_fnMask).defaultValueIfSet("mask"));

    _options.addOption(SelectionOption("select").required().multiValue()
                           .storeVector(&_sel));

    _options.addOption(BooleanOption("norm").store(&_bTotNorm)
        .description("Normalize by total number of positions with -os"));
    _options.addOption(BooleanOption("cfnorm").store(&_bFracNorm)
        .description("Normalize by covered fraction with -os"));
    _options.addOption(BooleanOption("resind").store(&_bResInd)
        .description("Write unique residue index instead of residue number"));

    return &_options;
}


int
Select::initAnalysis(const TopologyInformation &top)
{
    for (std::vector<Selection *>::const_iterator i = _sel.begin(); i != _sel.end(); ++i)
    {
        (*i)->initCoveredFraction(CFRAC_SOLIDANGLE);
    }

    _sdata.setColumns(_sel.size());
    registerAnalysisDataset(&_sdata, "size");
    snew(_totsize, _sel.size());
    for (size_t g = 0; g < _sel.size(); ++g)
    {
        _totsize[g] = _bTotNorm ? _sel[g]->posCount() : 1;
    }
    if (!_fnSize.empty())
    {
        AnalysisDataPlotModule *plot = new AnalysisDataPlotModule(_options);
        plot->setFileName(_fnSize);
        plot->setTitle("Selection size");
        plot->setXLabel("Time [ps]");
        plot->setYLabel("Number");
        _sdata.addModule(plot);
    }

    _cdata.setColumns(_sel.size());
    registerAnalysisDataset(&_cdata, "cfrac");
    if (!_fnFrac.empty())
    {
        AnalysisDataPlotModule *plot = new AnalysisDataPlotModule(_options);
        plot->setFileName(_fnFrac);
        plot->setTitle("Covered fraction");
        plot->setXLabel("Time [ps]");
        plot->setYLabel("Fraction");
        plot->setYFormat(6, 4);
        _cdata.addModule(plot);
    }

    _idata.setColumns(2, true);
    registerAnalysisDataset(&_idata, "index");
    if (!_fnIndex.empty())
    {
        AnalysisDataPlotModule *plot = new AnalysisDataPlotModule(_options);
        plot->setFileName(_fnIndex);
        plot->setPlainOutput(true);
        plot->setYFormat(4, 0);
        _idata.addModule(plot);
    }

    if (!_fnNdx.empty())
    {
        _block = new_blocka();
        _gnames = NULL;
        _modnames.reserve(_sel.size());
        for (size_t g = 0; g < _sel.size(); ++g)
        {
            _modnames.push_back(_sel[g]->name());
            size_t pos;
            while ((pos = _modnames[g].find(' ')) != std::string::npos)
            {
                _modnames[g][pos] = '_';
            }

            if (!_sel[g]->isDynamic())
            {
                add_grp(_block, &_gnames, _sel[g]->posCount(),
                        _sel[g]->mapIds(), _modnames[g].c_str());
                _modnames[g].clear();
            }
        }
    }

    _mdata.setColumns(_sel[0]->posCount());
    registerAnalysisDataset(&_mdata, "mask");
    if (!_fnMask.empty())
    {
        if (_sel.size() > 1U)
        {
            fprintf(stderr, "WARNING: the mask (-om) will only be written for the first group\n");
        }
        if (!_sel[0]->isDynamic())
        {
            fprintf(stderr, "WARNING: will not write the mask (-om) for a static selection\n");
        }
        else
        {
            AnalysisDataPlotModule *plot = new AnalysisDataPlotModule(_options);
            plot->setFileName(_fnMask);
            plot->setTitle("Selection mask");
            plot->setXLabel("Time [ps]");
            plot->setYLabel("Occupancy");
            plot->setYFormat(1, 0);
            _mdata.addModule(plot);
        }
    }

    _top = top.topology();

    return 0;
}


int
Select::startFrames(AnalysisDataParallelOptions opt,
                    const SelectionCollection &selections,
                    TrajectoryAnalysisModuleData **pdatap)
{
    ModuleData *pdata = new ModuleData();

    *pdatap = pdata;
    int rc = pdata->init(this, opt, selections);
    if (rc != 0)
    {
        delete pdata;
        *pdatap = NULL;
        return rc;
    }
    snew(pdata->_mmap, 1);
    gmx_ana_indexmap_init(pdata->_mmap, pdata->parallelSelection(_sel[0])->indexGroup(),
                          _top, _sel[0]->type());
    return 0;
}


int
Select::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                     TrajectoryAnalysisModuleData *pdata)
{
    ModuleData *d = static_cast<ModuleData *>(pdata);
    AnalysisDataHandle *sdh = pdata->dataHandle("size");
    AnalysisDataHandle *cdh = pdata->dataHandle("cfrac");
    AnalysisDataHandle *idh = pdata->dataHandle("index");
    AnalysisDataHandle *mdh = pdata->dataHandle("mask");
    std::vector<Selection *> sel(pdata->parallelSelections(_sel));

    if (sdh != NULL)
    {
        sdh->startFrame(frnr, fr.time);
        for (size_t g = 0; g < sel.size(); ++g)
        {
            real normfac = _bFracNorm ? 1.0 / sel[g]->cfrac() : 1.0;
            normfac /= _totsize[g];
            sdh->addPoint(g, sel[g]->posCount() * normfac);
        }
        sdh->finishFrame();
    }

    if (cdh != NULL)
    {
        cdh->startFrame(frnr, fr.time);
        for (size_t g = 0; g < sel.size(); ++g)
        {
            cdh->addPoint(g, sel[g]->cfrac());
        }
        cdh->finishFrame();
    }

    if (idh != NULL)
    {
        idh->startFrame(frnr, fr.time);
        for (size_t g = 0; g < sel.size(); ++g)
        {
            idh->addPoint(0, sel[g]->posCount());
            for (int i = 0; i < sel[g]->posCount(); ++i)
            {
                if (sel[g]->type() == INDEX_RES && !_bResInd)
                {
                    idh->addPoint(1, _top->atoms.resinfo[sel[g]->mapId(i)].nr);
                }
                else
                {
                    idh->addPoint(1, sel[g]->mapId(i) + 1);
                }
            }
        }
        idh->finishFrame();
    }

    /** TODO: This is not thread-safe */
    if (_block != NULL)
    {
        for (size_t g = 0; g < sel.size(); ++g)
        {
            if (sel[g]->isDynamic())
            {
                char tbuf[50];
                char *buf;

                sprintf(tbuf, "_%.3f", fr.time);
                snew(buf, _modnames[g].size() + strlen(tbuf) + 1);
                strcpy(buf, _modnames[g].c_str());
                strcat(buf, tbuf);
                add_grp(_block, &_gnames, sel[g]->posCount(),
                        sel[g]->mapIds(), buf);
                sfree(buf);
            }
        }
    }

    if (mdh != NULL)
    {
        gmx_ana_indexmap_update(d->_mmap, sel[0]->indexGroup(), TRUE);
        mdh->startFrame(frnr, fr.time);
        for (int b = 0; b < d->_mmap->nr; ++b)
        {
            mdh->addPoint(b, d->_mmap->refid[b] == -1 ? 0 : 1);
        }
        mdh->finishFrame();
    }
    return 0;
}


int
Select::finishAnalysis(int /*nframes*/)
{
    return 0;
}


int
Select::writeOutput()
{
    if (!_fnNdx.empty())
    {
        write_index(_fnNdx.c_str(), _block, _gnames);
    }

    return 0;
}


TrajectoryAnalysisModule *
Select::create()
{
    return new Select();
}

} // namespace analysismodules

} // namespace gmx
