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
 * Implements classes in plot.h.
 *
 * \ingroup module_analysisdata
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#include "gromacs/analysisdata/modules/plot.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include <cstdio>
#include <cstring>

#include <boost/shared_ptr.hpp>

#include "gromacs/legacyheaders/gmxfio.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/vec.h"
#include "gromacs/legacyheaders/xvgr.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/options/options.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"
#include "gromacs/utility/gmxassert.h"

static const char *const g_plotFormats[] = {
    "none", "xmgrace", "xmgr", NULL
};

namespace gmx
{

/********************************************************************
 * AnalysisDataPlotSettings
 */

AnalysisDataPlotSettings::AnalysisDataPlotSettings()
    : selections_(NULL), timeUnit_(eTimeUnit_ps), plotFormat_(1)
{
}

void
AnalysisDataPlotSettings::setSelectionCollection(const SelectionCollection *selections)
{
    selections_ = selections;
}


void
AnalysisDataPlotSettings::addOptions(Options *options)
{
    options->addOption(StringOption("xvg").enumValue(g_plotFormats)
                           .defaultValue("xmgrace")
                           .storeEnumIndex(&plotFormat_)
                           .description("Plot formatting"));
}


/********************************************************************
 * AbstractPlotModule::Impl
 */

class AbstractPlotModule::Impl
{
    public:
        explicit Impl(const AnalysisDataPlotSettings &settings);
        ~Impl();

        void closeFile();

        AnalysisDataPlotSettings settings;
        std::string             fnm;
        FILE                   *fp;

        bool                    bPlain;
        bool                    bOmitX;
        std::string             title;
        std::string             subtitle;
        std::string             xlabel;
        std::string             ylabel;
        std::vector<std::string>  leg;
        char                    xfmt[15];
        char                    yfmt[15];
        real                    xscale;
};

AbstractPlotModule::Impl::Impl(const AnalysisDataPlotSettings &settings)
    : settings(settings), fp(NULL), bPlain(false), bOmitX(false), xscale(1.0)
{
    strcpy(xfmt, "%11.3f");
    strcpy(yfmt, " %8.3f");
}

AbstractPlotModule::Impl::~Impl()
{
    closeFile();
}


void
AbstractPlotModule::Impl::closeFile()
{
    if (fp != NULL)
    {
        if (bPlain)
        {
            gmx_fio_fclose(fp);
        }
        else
        {
            xvgrclose(fp);
        }
        fp = NULL;
    }
}


/********************************************************************
 * AbstractPlotModule
 */
/*! \cond libapi */
AbstractPlotModule::AbstractPlotModule()
    : _impl(new Impl(AnalysisDataPlotSettings()))
{
}

AbstractPlotModule::AbstractPlotModule(const AnalysisDataPlotSettings &settings)
    : _impl(new Impl(settings))
{
}
//! \endcond

AbstractPlotModule::~AbstractPlotModule()
{
}


void
AbstractPlotModule::setSettings(const AnalysisDataPlotSettings &settings)
{
    _impl->settings = settings;
}


void
AbstractPlotModule::setFileName(const std::string &fnm)
{
    _impl->fnm = fnm;
}


void
AbstractPlotModule::setPlainOutput(bool bPlain)
{
    _impl->bPlain = bPlain;
}


void
AbstractPlotModule::setOmitX(bool bOmitX)
{
    _impl->bOmitX = bOmitX;
}


void
AbstractPlotModule::setTitle(const char *title)
{
    _impl->title = title;
}


void
AbstractPlotModule::setSubtitle(const char *subtitle)
{
    _impl->subtitle = subtitle;
}


void
AbstractPlotModule::setXLabel(const char *label)
{
    _impl->xlabel = label;
}


void
AbstractPlotModule::setXAxisIsTime()
{
    TimeUnitManager manager(_impl->settings.timeUnit());
    _impl->xlabel = formatString("Time (%s)", manager.timeUnitAsString());
    _impl->xscale = manager.inverseTimeScaleFactor();
}


void
AbstractPlotModule::setYLabel(const char *label)
{
    _impl->ylabel = label;
}


void
AbstractPlotModule::setLegend(int nsets, const char * const *setname)
{
    _impl->leg.reserve(_impl->leg.size() + nsets);
    for (int i = 0; i < nsets; ++i)
    {
        appendLegend(setname[i]);
    }
}


void
AbstractPlotModule::appendLegend(const char *setname)
{
    _impl->leg.push_back(setname);
}


void
AbstractPlotModule::setXFormat(int width, int prec, char fmt)
{
    GMX_RELEASE_ASSERT(width >= 0 && prec >= 0 && width <= 99 && prec <= 99,
                       "Invalid width or precision");
    GMX_RELEASE_ASSERT(strchr("eEfFgG", fmt) != NULL,
                       "Invalid format specifier");
    sprintf(_impl->xfmt, "%%%d.%d%c", width, prec, fmt);
}


void
AbstractPlotModule::setYFormat(int width, int prec, char fmt)
{
    GMX_RELEASE_ASSERT(width >= 0 && prec >= 0 && width <= 99 && prec <= 99,
                       "Invalid width or precision");
    GMX_RELEASE_ASSERT(strchr("eEfFgG", fmt) != NULL,
                       "Invalid format specifier");
    sprintf(_impl->yfmt, " %%%d.%d%c", width, prec, fmt);
}


int
AbstractPlotModule::flags() const
{
    return efAllowMulticolumn | efAllowMultipoint;
}


void
AbstractPlotModule::dataStarted(AbstractAnalysisData *data)
{
    if (!_impl->fnm.empty())
    {
        if (_impl->bPlain)
        {
            _impl->fp = gmx_fio_fopen(_impl->fnm.c_str(), "w");
        }
        else
        {
            time_unit_t time_unit
                = static_cast<time_unit_t>(_impl->settings.timeUnit() + 1);
            xvg_format_t xvg_format
                = (_impl->settings.plotFormat() > 0
                    ? static_cast<xvg_format_t>(_impl->settings.plotFormat())
                    : exvgNONE);
            output_env_t oenv;
            output_env_init(&oenv, 0, NULL, time_unit, FALSE, xvg_format, 0, 0);
            boost::shared_ptr<output_env> oenvGuard(oenv, &output_env_done);
            _impl->fp = xvgropen(_impl->fnm.c_str(), _impl->title.c_str(),
                                 _impl->xlabel.c_str(), _impl->ylabel.c_str(),
                                 oenv);
            const SelectionCollection *selections
                = _impl->settings.selectionCollection();
            if (selections != NULL)
            {
                selections->printXvgrInfo(_impl->fp, oenv);
            }
            if (!_impl->subtitle.empty())
            {
                xvgr_subtitle(_impl->fp, _impl->subtitle.c_str(), oenv);
            }
            if (output_env_get_print_xvgr_codes(oenv)
                && !_impl->leg.empty())
            {
                std::vector<const char *> leg;
                leg.reserve(_impl->leg.size());
                for (size_t i = 0; i < _impl->leg.size(); ++i)
                {
                    leg.push_back(_impl->leg[i].c_str());
                }
                xvgr_legend(_impl->fp, leg.size(), &leg[0], oenv);
            }
        }
    }
}


void
AbstractPlotModule::frameStarted(const AnalysisDataFrameHeader &frame)
{
    if (!isFileOpen())
    {
        return;
    }
    if (!_impl->bOmitX)
    {
        std::fprintf(_impl->fp, _impl->xfmt, frame.x() * _impl->xscale);
    }
}


void
AbstractPlotModule::frameFinished(const AnalysisDataFrameHeader & /*header*/)
{
    if (!isFileOpen())
    {
        return;
    }
    std::fprintf(_impl->fp, "\n");
}


void
AbstractPlotModule::dataFinished()
{
    _impl->closeFile();
}

/*! \cond libapi */
bool
AbstractPlotModule::isFileOpen() const
{
    return _impl->fp != NULL;
}


void
AbstractPlotModule::writeValue(real value) const
{
    GMX_ASSERT(isFileOpen(), "File not opened, but write attempted");
    std::fprintf(_impl->fp, _impl->yfmt, value);
}
//! \endcond

/********************************************************************
 * DataPlotModule
 */

AnalysisDataPlotModule::AnalysisDataPlotModule()
{
}

AnalysisDataPlotModule::AnalysisDataPlotModule(
        const AnalysisDataPlotSettings &settings)
    : AbstractPlotModule(settings)
{
}


void
AnalysisDataPlotModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (!isFileOpen())
    {
        return;
    }
    for (int i = 0; i < points.columnCount(); ++i)
    {
        writeValue(points.y(i));
    }
}


/********************************************************************
 * DataVectorPlotModule
 */

AnalysisDataVectorPlotModule::AnalysisDataVectorPlotModule()
{
    for (int i = 0; i < DIM; ++i)
    {
        _bWrite[i] = true;
    }
    _bWrite[DIM] = false;
}


AnalysisDataVectorPlotModule::AnalysisDataVectorPlotModule(
        const AnalysisDataPlotSettings &settings)
    : AbstractPlotModule(settings)
{
    for (int i = 0; i < DIM; ++i)
    {
        _bWrite[i] = true;
    }
    _bWrite[DIM] = false;
}


void
AnalysisDataVectorPlotModule::setWriteX(bool bWrite)
{
    _bWrite[XX] = bWrite;
}


void
AnalysisDataVectorPlotModule::setWriteY(bool bWrite)
{
    _bWrite[YY] = bWrite;
}


void
AnalysisDataVectorPlotModule::setWriteZ(bool bWrite)
{
    _bWrite[ZZ] = bWrite;
}


void
AnalysisDataVectorPlotModule::setWriteNorm(bool bWrite)
{
    _bWrite[DIM] = bWrite;
}


void
AnalysisDataVectorPlotModule::setWriteMask(bool bWrite[DIM + 1])
{
    for (int i = 0; i < DIM + 1; ++i)
    {
        _bWrite[i] = bWrite[i];
    }
}


void
AnalysisDataVectorPlotModule::pointsAdded(const AnalysisDataPointSetRef &points)
{
    if (points.firstColumn() % DIM != 0)
    {
        GMX_THROW(APIError("Partial data points"));
    }
    if (!isFileOpen())
    {
        return;
    }
    for (int i = 0; i < points.columnCount(); i += 3)
    {
        for (int d = 0; d < DIM; ++d)
        {
            if (_bWrite[i])
            {
                writeValue(points.y(i + d));
            }
        }
        if (_bWrite[DIM])
        {
            rvec y = { points.y(i), points.y(i + 1), points.y(i + 2) };
            writeValue(norm(y));
        }
    }
}

} // namespace gmx
