/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Implements classes in plot.h.
 *
 * \ingroup module_analysisdata
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "gmxpre.h"

#include "gromacs/analysisdata/modules/plot.h"

#include <cstdio>
#include <cstring>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/unique_cptr.h"

struct gmx_output_env_t;

namespace gmx
{
class AbstractAnalysisData;

/********************************************************************
 * AnalysisDataPlotSettings
 */

AnalysisDataPlotSettings::AnalysisDataPlotSettings() :
    selections_(nullptr), timeUnit_(TimeUnit::Default), plotFormat_(XvgFormat::Xmgrace)
{
}

void AnalysisDataPlotSettings::setSelectionCollection(const SelectionCollection* selections)
{
    selections_ = selections;
}

/*! \brief Names for XvgFormat
 *
 * Technically this duplicates a definition in pargs.cpp for legacy
 * support code, but as the latter will go away and the alternatives
 * are ugly, the duplication is acceptable. */
const gmx::EnumerationArray<XvgFormat, const char*> c_xvgFormatNames = {
    { "xmgrace", "xmgr", "none" }
};

void AnalysisDataPlotSettings::initOptions(IOptionsContainer* options)
{
    options->addOption(
            EnumOption<XvgFormat>("xvg").enumValue(c_xvgFormatNames).store(&plotFormat_).description("Plot formatting"));
}


/********************************************************************
 * AbstractPlotModule::Impl
 */

class AbstractPlotModule::Impl
{
public:
    explicit Impl(const AnalysisDataPlotSettings& settings);
    ~Impl();

    void closeFile();

    AnalysisDataPlotSettings settings_;
    std::string              filename_;
    FILE*                    fp_;

    bool                     bPlain_;
    bool                     bOmitX_;
    bool                     bErrorsAsSeparateColumn_;
    std::string              title_;
    std::string              subtitle_;
    std::string              xlabel_;
    std::string              ylabel_;
    std::vector<std::string> legend_;
    std::string              xformat_;
    std::string              yformat_;
    real                     xscale_;
};

AbstractPlotModule::Impl::Impl(const AnalysisDataPlotSettings& settings) :
    settings_(settings),
    fp_(nullptr),
    bPlain_(false),
    bOmitX_(false),
    bErrorsAsSeparateColumn_(false),
    xformat_("%11.3f"),
    yformat_(" %8.3f"),
    xscale_(1.0)
{
}

AbstractPlotModule::Impl::~Impl()
{
    closeFile();
}


void AbstractPlotModule::Impl::closeFile()
{
    if (fp_ != nullptr)
    {
        if (bPlain_)
        {
            gmx_fio_fclose(fp_);
        }
        else
        {
            xvgrclose(fp_);
        }
        fp_ = nullptr;
    }
}


/********************************************************************
 * AbstractPlotModule
 */
/*! \cond libapi */
AbstractPlotModule::AbstractPlotModule() : impl_(new Impl(AnalysisDataPlotSettings())) {}

AbstractPlotModule::AbstractPlotModule(const AnalysisDataPlotSettings& settings) :
    impl_(new Impl(settings))
{
}
//! \endcond

AbstractPlotModule::~AbstractPlotModule() {}


void AbstractPlotModule::setSettings(const AnalysisDataPlotSettings& settings)
{
    impl_->settings_ = settings;
}


void AbstractPlotModule::setFileName(const std::string& filename)
{
    impl_->filename_ = filename;
}


void AbstractPlotModule::setPlainOutput(bool bPlain)
{
    impl_->bPlain_ = bPlain;
}


void AbstractPlotModule::setErrorsAsSeparateColumn(bool bSeparate)
{
    impl_->bErrorsAsSeparateColumn_ = bSeparate;
}


void AbstractPlotModule::setOmitX(bool bOmitX)
{
    impl_->bOmitX_ = bOmitX;
}


void AbstractPlotModule::setTitle(const char* title)
{
    impl_->title_ = title;
}

void AbstractPlotModule::setTitle(const std::string& title)
{
    impl_->title_ = title;
}


void AbstractPlotModule::setSubtitle(const char* subtitle)
{
    impl_->subtitle_ = subtitle;
}


void AbstractPlotModule::setSubtitle(const std::string& subtitle)
{
    impl_->subtitle_ = subtitle;
}


void AbstractPlotModule::setXLabel(const char* label)
{
    impl_->xlabel_ = label;
}


void AbstractPlotModule::setXAxisIsTime()
{
    TimeUnitManager manager(impl_->settings_.timeUnit());
    impl_->xlabel_ = formatString("Time (%s)", manager.timeUnitAsString());
    impl_->xscale_ = manager.inverseTimeScaleFactor();
}


void AbstractPlotModule::setYLabel(const char* label)
{
    impl_->ylabel_ = label;
}


void AbstractPlotModule::setLegend(int nsets, const char* const* setname)
{
    impl_->legend_.reserve(impl_->legend_.size() + nsets);
    for (int i = 0; i < nsets; ++i)
    {
        appendLegend(setname[i]);
    }
}


void AbstractPlotModule::appendLegend(const char* setname)
{
    impl_->legend_.emplace_back(setname);
}


void AbstractPlotModule::appendLegend(const std::string& setname)
{
    impl_->legend_.push_back(setname);
}


void AbstractPlotModule::setXFormat(int width, int precision, char format)
{
    GMX_RELEASE_ASSERT(width >= 0 && precision >= 0 && width <= 99 && precision <= 99,
                       "Invalid width or precision");
    GMX_RELEASE_ASSERT(strchr("eEfFgG", format) != nullptr, "Invalid format specifier");
    impl_->xformat_ = formatString("%%%d.%d%c", width, precision, format);
}


void AbstractPlotModule::setYFormat(int width, int precision, char format)
{
    GMX_RELEASE_ASSERT(width >= 0 && precision >= 0 && width <= 99 && precision <= 99,
                       "Invalid width or precision");
    GMX_RELEASE_ASSERT(strchr("eEfFgG", format) != nullptr, "Invalid format specifier");
    impl_->yformat_ = formatString(" %%%d.%d%c", width, precision, format);
}


int AbstractPlotModule::flags() const
{
    return efAllowMissing | efAllowMulticolumn | efAllowMultipoint | efAllowMultipleDataSets;
}


void AbstractPlotModule::dataStarted(AbstractAnalysisData* /* data */)
{
    if (!impl_->filename_.empty())
    {
        if (impl_->bPlain_)
        {
            impl_->fp_ = gmx_fio_fopen(impl_->filename_.c_str(), "w");
        }
        else
        {
            const TimeUnit    timeUnit  = impl_->settings_.timeUnit();
            const XvgFormat   xvgFormat = impl_->settings_.plotFormat();
            gmx_output_env_t* oenv;
            output_env_init(&oenv, getProgramContext(), timeUnit, FALSE, xvgFormat, 0);
            const unique_cptr<gmx_output_env_t, output_env_done> oenvGuard(oenv);
            impl_->fp_ = xvgropen(
                    impl_->filename_.c_str(), impl_->title_.c_str(), impl_->xlabel_, impl_->ylabel_, oenv);
            const SelectionCollection* selections = impl_->settings_.selectionCollection();
            if (selections != nullptr && output_env_get_xvg_format(oenv) != XvgFormat::None)
            {
                selections->printXvgrInfo(impl_->fp_);
            }
            if (!impl_->subtitle_.empty())
            {
                xvgr_subtitle(impl_->fp_, impl_->subtitle_.c_str(), oenv);
            }
            if (output_env_get_print_xvgr_codes(oenv) && !impl_->legend_.empty())
            {
                xvgrLegend(impl_->fp_, impl_->legend_, oenv);
            }
        }
    }
}


void AbstractPlotModule::frameStarted(const AnalysisDataFrameHeader& header)
{
    if (!isFileOpen())
    {
        return;
    }
    if (!impl_->bOmitX_)
    {
        std::fprintf(impl_->fp_, impl_->xformat_.c_str(), header.x() * impl_->xscale_);
    }
}


void AbstractPlotModule::frameFinished(const AnalysisDataFrameHeader& /*header*/)
{
    if (!isFileOpen())
    {
        return;
    }
    std::fprintf(impl_->fp_, "\n");
}


void AbstractPlotModule::dataFinished()
{
    impl_->closeFile();
}

/*! \cond libapi */
bool AbstractPlotModule::isFileOpen() const
{
    return impl_->fp_ != nullptr;
}


void AbstractPlotModule::writeValue(const AnalysisDataValue& value) const
{
    GMX_ASSERT(isFileOpen(), "File not opened, but write attempted");
    const real y = value.isSet() ? value.value() : 0.0;
    std::fprintf(impl_->fp_, impl_->yformat_.c_str(), y);
    if (impl_->bErrorsAsSeparateColumn_)
    {
        const real dy = value.isSet() ? value.error() : 0.0;
        std::fprintf(impl_->fp_, impl_->yformat_.c_str(), dy);
    }
}
//! \endcond

/********************************************************************
 * DataPlotModule
 */

AnalysisDataPlotModule::AnalysisDataPlotModule() {}

AnalysisDataPlotModule::AnalysisDataPlotModule(const AnalysisDataPlotSettings& settings) :
    AbstractPlotModule(settings)
{
}


void AnalysisDataPlotModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (!isFileOpen())
    {
        return;
    }
    for (int i = 0; i < points.columnCount(); ++i)
    {
        writeValue(points.values()[i]);
    }
}


/********************************************************************
 * DataVectorPlotModule
 */

AnalysisDataVectorPlotModule::AnalysisDataVectorPlotModule() : bWrite_{ true, true, true, false } {}


AnalysisDataVectorPlotModule::AnalysisDataVectorPlotModule(const AnalysisDataPlotSettings& settings) :
    AbstractPlotModule(settings), bWrite_{ true, true, true, false }
{
}


void AnalysisDataVectorPlotModule::setWriteX(bool bWrite)
{
    bWrite_[XX] = bWrite;
}


void AnalysisDataVectorPlotModule::setWriteY(bool bWrite)
{
    bWrite_[YY] = bWrite;
}


void AnalysisDataVectorPlotModule::setWriteZ(bool bWrite)
{
    bWrite_[ZZ] = bWrite;
}


void AnalysisDataVectorPlotModule::setWriteNorm(bool bWrite)
{
    bWrite_[DIM] = bWrite;
}


void AnalysisDataVectorPlotModule::setWriteMask(const bool bWrite[DIM + 1])
{
    for (int i = 0; i < DIM + 1; ++i)
    {
        bWrite_[i] = bWrite[i];
    }
}


void AnalysisDataVectorPlotModule::pointsAdded(const AnalysisDataPointSetRef& points)
{
    if (points.firstColumn() % DIM != 0 || points.columnCount() % DIM != 0)
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
            if (bWrite_[d])
            {
                writeValue(points.values()[i + d]);
            }
        }
        if (bWrite_[DIM])
        {
            const rvec        y = { points.y(i), points.y(i + 1), points.y(i + 2) };
            AnalysisDataValue value(norm(y));
            writeValue(value);
        }
    }
}

} // namespace gmx
