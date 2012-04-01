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
 * Implements gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/abstractdata.h"

#include "gromacs/fatalerror/exceptions.h"
#include "gromacs/fatalerror/gmxassert.h"

#include "abstractdata-impl.h"
#include "dataframe.h"
#include "dataproxy.h"

namespace gmx
{

/********************************************************************
 * AbstractAnalysisData::Impl
 */

AbstractAnalysisData::Impl::Impl()
    : _bAllowMissing(true), _bDataStart(false), _bInData(false), _bInFrame(false),
      _currIndex(-1), _nframes(0)
{
}

AbstractAnalysisData::Impl::~Impl()
{
}


void
AbstractAnalysisData::Impl::presentData(AbstractAnalysisData *data,
                                        AnalysisDataModuleInterface *module)
{
    module->dataStarted(data);
    bool bCheckMissing = _bAllowMissing
        && !(module->flags() & AnalysisDataModuleInterface::efAllowMissing);
    for (int i = 0; i < data->frameCount(); ++i)
    {
        AnalysisDataFrameRef frame = data->getDataFrame(i);
        GMX_RELEASE_ASSERT(frame.isValid(), "Invalid data frame returned");
        // TODO: Check all frames before doing anything for slightly better
        // exception behavior.
        if (bCheckMissing && !frame.allPresent())
        {
            GMX_THROW(APIError("Missing data not supported by a module"));
        }
        module->frameStarted(frame.header());
        module->pointsAdded(frame.points());
        module->frameFinished(frame.header());
    }
    if (!_bInData)
    {
        module->dataFinished();
    }
}


/********************************************************************
 * AbstractAnalysisData
 */
/*! \cond libapi */
AbstractAnalysisData::AbstractAnalysisData()
    : _impl(new Impl()), _ncol(0), _bMultiPoint(false)
{
}
//! \endcond

AbstractAnalysisData::~AbstractAnalysisData()
{
}


int
AbstractAnalysisData::frameCount() const
{
    return _impl->_nframes;
}


AnalysisDataFrameRef
AbstractAnalysisData::tryGetDataFrame(int index) const
{
    if (index < 0 || index >= frameCount())
    {
        return AnalysisDataFrameRef();
    }
    return tryGetDataFrameInternal(index);
}


AnalysisDataFrameRef
AbstractAnalysisData::getDataFrame(int index) const
{
    AnalysisDataFrameRef frame = tryGetDataFrame(index);
    if (!frame.isValid())
    {
        GMX_THROW(APIError("Invalid frame accessed"));
    }
    return frame;
}


bool
AbstractAnalysisData::requestStorage(int nframes)
{
    GMX_RELEASE_ASSERT(nframes >= -1, "Invalid number of frames requested");
    if (nframes == 0)
    {
        return true;
    }
    return requestStorageInternal(nframes);
}


void
AbstractAnalysisData::addModule(AnalysisDataModulePointer module)
{
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }

    if (_impl->_bDataStart)
    {
        GMX_RELEASE_ASSERT(!_impl->_bInFrame,
                           "Cannot add data modules in mid-frame");
        _impl->presentData(this, module.get());
    }
    if (!(module->flags() & AnalysisDataModuleInterface::efAllowMissing))
    {
        _impl->_bAllowMissing = false;
    }
    _impl->_modules.push_back(move(module));
}


void
AbstractAnalysisData::addColumnModule(int col, int span,
                                      AnalysisDataModulePointer module)
{
    GMX_RELEASE_ASSERT(col >= 0 && span >= 1 && col + span <= _ncol,
                       "Invalid columns specified for a column module");
    if (_impl->_bDataStart)
    {
        GMX_THROW(NotImplementedError("Cannot add column modules after data"));
    }

    boost::shared_ptr<AnalysisDataProxy> proxy(
            new AnalysisDataProxy(col, span, this));
    proxy->addModule(module);
    addModule(proxy);
}


void
AbstractAnalysisData::applyModule(AnalysisDataModuleInterface *module)
{
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_THROW(APIError("Data module not compatible with data object properties"));
    }
    GMX_RELEASE_ASSERT(_impl->_bDataStart && !_impl->_bInData,
                       "Data module can only be applied to ready data");

    _impl->presentData(this, module);
}

/*! \cond libapi */
void
AbstractAnalysisData::setColumnCount(int ncol)
{
    GMX_RELEASE_ASSERT(ncol > 0, "Invalid data column count");
    GMX_RELEASE_ASSERT(_ncol == 0 || _impl->_modules.empty(),
                       "Data column count cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "Data column count cannot be changed after data has been added");
    _ncol = ncol;
}


void
AbstractAnalysisData::setMultipoint(bool multipoint)
{
    GMX_RELEASE_ASSERT(_impl->_modules.empty(),
                       "Data type cannot be changed after modules are added");
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "Data type cannot be changed after data has been added");
    _bMultiPoint = multipoint;
}


/*! \internal
 * This method is not const because the dataStarted() methods of the attached
 * modules can request storage of the data.
 */
void
AbstractAnalysisData::notifyDataStart()
{
    GMX_RELEASE_ASSERT(!_impl->_bDataStart,
                       "notifyDataStart() called more than once");
    GMX_RELEASE_ASSERT(_ncol > 0, "Data column count is not set");
    _impl->_bDataStart = _impl->_bInData = true;

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        if (_ncol > 1 && !((*i)->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        {
            GMX_THROW(APIError("Data module not compatible with data object properties"));
        }
        (*i)->dataStarted(this);
    }
}


void
AbstractAnalysisData::notifyFrameStart(const AnalysisDataFrameHeader &header) const
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(!_impl->_bInFrame,
               "notifyFrameStart() called while inside a frame");
    GMX_ASSERT(header.index() == _impl->_nframes,
               "Out of order frames");
    _impl->_bInFrame = true;
    _impl->_currIndex = header.index();

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->frameStarted(header);
    }
}


void
AbstractAnalysisData::notifyPointsAdd(const AnalysisDataPointSetRef &points) const
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(_impl->_bInFrame, "notifyFrameStart() not called");
    GMX_ASSERT(points.lastColumn() < columnCount(), "Invalid columns");
    GMX_ASSERT(points.frameIndex() == _impl->_currIndex,
               "Points do not correspond to current frame");
    if (!_impl->_bAllowMissing && !points.allPresent())
    {
        GMX_THROW(APIError("Missing data not supported by a module"));
    }

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->pointsAdded(points);
    }
}


void
AbstractAnalysisData::notifyFrameFinish(const AnalysisDataFrameHeader &header)
{
    GMX_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_ASSERT(_impl->_bInFrame, "notifyFrameStart() not called");
    GMX_ASSERT(header.index() == _impl->_currIndex,
               "Header does not correspond to current frame");
    _impl->_bInFrame = false;
    _impl->_currIndex = -1;

    // Increment the counter before notifications to allow frame access from
    // modules.
    ++_impl->_nframes;

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->frameFinished(header);
    }
}


void
AbstractAnalysisData::notifyDataFinish() const
{
    GMX_RELEASE_ASSERT(_impl->_bInData, "notifyDataStart() not called");
    GMX_RELEASE_ASSERT(!_impl->_bInFrame,
                       "notifyDataFinish() called while inside a frame");
    _impl->_bInData = false;

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        (*i)->dataFinished();
    }
}
//! \endcond

} // namespace gmx
