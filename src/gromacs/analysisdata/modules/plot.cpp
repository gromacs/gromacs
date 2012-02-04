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

#include <gmxfio.h>
#include <statutil.h>
#include <vec.h>
#include <xvgr.h>

#include "gromacs/options/globalproperties.h"
#include "gromacs/options/options.h"
#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"
#include "gromacs/selection/selectioncollection.h"

#include "plot-impl.h"

namespace gmx
{

/********************************************************************
 * AbstractPlotModule::Impl
 */

AbstractPlotModule::Impl::Impl(const Options &options)
    : fp(NULL), bPlain(false), bOmitX(false),
      oenv(options.globalProperties().output_env()),
      sel(NULL)
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

AbstractPlotModule::AbstractPlotModule(const Options &options)
    : _impl(new Impl(options))
{
}


AbstractPlotModule::~AbstractPlotModule()
{
    delete _impl;
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
AbstractPlotModule::setXTimeLabel()
{
    _impl->xlabel = output_env_get_xvgr_tlabel(_impl->oenv);
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
            _impl->fp = xvgropen(_impl->fnm.c_str(), _impl->title.c_str(),
                                 _impl->xlabel.c_str(), _impl->ylabel.c_str(),
                                 _impl->oenv);
            if (_impl->sel != NULL)
            {
                _impl->sel->printXvgrInfo(_impl->fp, _impl->oenv);
            }
            if (!_impl->subtitle.empty())
            {
                xvgr_subtitle(_impl->fp, _impl->subtitle.c_str(), _impl->oenv);
            }
            if (output_env_get_print_xvgr_codes(_impl->oenv)
                && !_impl->leg.empty())
            {
                const char **leg;

                leg = new const char *[_impl->leg.size()];
                for (size_t i = 0; i < _impl->leg.size(); ++i)
                {
                    leg[i] = _impl->leg[i].c_str();
                }
                xvgr_legend(_impl->fp, _impl->leg.size(), leg, _impl->oenv);
                delete [] leg;
            }
        }
    }
}


void
AbstractPlotModule::frameStarted(real x, real dx)
{
    if (!isFileOpen())
    {
        return;
    }
    if (!_impl->bOmitX)
    {
        std::fprintf(_impl->fp, _impl->xfmt, x);
    }
}


void
AbstractPlotModule::frameFinished()
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


/********************************************************************
 * DataPlotModule
 */

AnalysisDataPlotModule::AnalysisDataPlotModule(const Options &options)
    : AbstractPlotModule(options)
{
}


void
AnalysisDataPlotModule::pointsAdded(real x, real dx, int firstcol, int n,
                                    const real *y, const real *dy,
                                    const bool *present)
{
    if (!isFileOpen())
    {
        return;
    }
    for (int i = 0; i < n; ++i)
    {
        writeValue(y[i]);
    }
}


/********************************************************************
 * DataVectorPlotModule
 */

AnalysisDataVectorPlotModule::AnalysisDataVectorPlotModule(const Options &options)
    : AbstractPlotModule(options)
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
AnalysisDataVectorPlotModule::pointsAdded(real x, real dx, int firstcol, int n,
                                          const real *y, const real *dy,
                                          const bool *present)
{
    if (firstcol % DIM != 0)
    {
        GMX_THROW(APIError("Partial data points"));
    }
    if (!isFileOpen())
    {
        return;
    }
    for (int i = 0; i < n; i += 3)
    {
        for (int d = 0; d < DIM; ++d)
        {
            if (_bWrite[i])
            {
                writeValue(y[i + d]);
            }
        }
        if (_bWrite[DIM])
        {
            writeValue(norm(&y[i]));
        }
    }
}

} // namespace gmx
