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
 * Implements gmx::AbstractAnalysisData and gmx::AbstractAnalysisDataStored.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/abstractdata.h"

#include <cassert>

// Legacy header.
#include "smalloc.h"

#include "gromacs/fatalerror/fatalerror.h"

#include "abstractdata-impl.h"
#include "dataproxy.h"

namespace gmx
{

/********************************************************************
 * AbstractAnalysisData::Impl
 */

AbstractAnalysisData::Impl::Impl()
    : _bDataStart(false), _bInData(false), _bInFrame(false),
      _bAllowMissing(true)
{
}

AbstractAnalysisData::Impl::~Impl()
{
    ModuleList::const_iterator i;
    for (i = _modules.begin(); i != _modules.end(); ++i)
    {
        delete *i;
    }
}


int
AbstractAnalysisData::Impl::presentData(AbstractAnalysisData *data,
                                        AnalysisDataModuleInterface *module)
{
    int rc = module->dataStarted(data);
    if (rc != 0)
    {
        return rc;
    }
    bool bCheckMissing = _bAllowMissing
        && !(module->flags() & AnalysisDataModuleInterface::efAllowMissing);
    int ncol = data->columnCount();
    for (int i = 0; i < data->frameCount(); ++i)
    {
        real        x, dx;
        const real *y, *dy;
        const bool *present;

        rc = data->getDataWErr(i, &x, &dx, &y, &dy, &present);
        if (rc == 0)
        {
            if (bCheckMissing && present)
            {
                for (int j = 0; j < ncol; ++j)
                {
                    if (!present[j])
                    {
                        GMX_ERROR(eeInvalidValue,
                                  "Missing data not supported by a module");
                    }
                }
            }
            rc = module->frameStarted(x, dx);
        }
        if (rc == 0)
        {
            rc = module->pointsAdded(x, dx, 0, ncol, y, dy, present);
        }
        if (rc == 0)
        {
            rc = module->frameFinished();
        }
        if (rc != 0)
        {
            return rc;
        }
    }
    if (!_bInData)
    {
        rc = module->dataFinished();
    }
    return rc;
}


/********************************************************************
 * AbstractAnalysisData
 */

AbstractAnalysisData::AbstractAnalysisData()
    : _impl(new Impl()), _ncol(0), _bMultiPoint(false)
{
}


AbstractAnalysisData::~AbstractAnalysisData()
{
    delete _impl;
}


int
AbstractAnalysisData::getData(int index, real *x, const real **y,
                              const bool **missing) const
{
    return getDataWErr(index, x, 0, y, 0, missing);
}


int
AbstractAnalysisData::getErrors(int index, real *dx, const real **dy) const
{
    return getDataWErr(index, 0, dx, 0, dy, 0);
}


int
AbstractAnalysisData::addModule(AnalysisDataModuleInterface *module)
{
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_ERROR(eeInvalidValue,
                  "Data module not compatible with data object properties");
    }

    if (_impl->_bDataStart)
    {
        if (_impl->_bInFrame)
        {
            GMX_ERROR(eeInvalidCall,
                      "Cannot add data modules in mid-frame");
        }
        int rc = _impl->presentData(this, module);
        if (rc != 0)
        {
            return rc;
        }
    }
    if (!(module->flags() & AnalysisDataModuleInterface::efAllowMissing))
    {
        _impl->_bAllowMissing = false;
    }
    _impl->_modules.push_back(module);
    return 0;
}


int
AbstractAnalysisData::addColumnModule(int col, int span,
                                      AnalysisDataModuleInterface *module)
{
    assert(col >= 0 && span >= 1 && col + span <= _ncol);
    if (_impl->_bDataStart)
    {
        GMX_ERROR(eeNotImplemented,
                  "Cannot add column modules after data");
    }

    AnalysisDataProxy *proxy = new AnalysisDataProxy(col, span, this);
    int rc = proxy->addModule(module);
    if (rc == 0)
    {
        rc = addModule(proxy);
    }
    if (rc != 0)
    {
        delete proxy;
    }
    return rc;
}


int
AbstractAnalysisData::applyModule(AnalysisDataModuleInterface *module)
{
    if ((columnCount() > 1 && !(module->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        || (isMultipoint() && !(module->flags() & AnalysisDataModuleInterface::efAllowMultipoint))
        || (!isMultipoint() && (module->flags() & AnalysisDataModuleInterface::efOnlyMultipoint)))
    {
        GMX_ERROR(eeInvalidValue,
                  "Data module not compatible with data object properties");
    }
    if (!_impl->_bDataStart || _impl->_bInData)
    {
        GMX_ERROR(eeInvalidCall,
                  "Data module can only be applied to ready data");
    }

    return _impl->presentData(this, module);
}


void
AbstractAnalysisData::setColumnCount(int ncol)
{
    assert(ncol > 0);
    assert(_ncol == 0 || _impl->_modules.empty());
    assert(!_impl->_bDataStart);
    _ncol = ncol;
}


void
AbstractAnalysisData::setMultipoint(bool multipoint)
{
    assert(_impl->_modules.empty());
    assert(!_impl->_bDataStart);
    _bMultiPoint = multipoint;
}


/*! \internal
 * This method is not const because the dataStarted() methods of the attached
 * modules can request storage of the data.
 */
int
AbstractAnalysisData::notifyDataStart()
{
    assert(!_impl->_bDataStart);
    assert(_ncol > 0);
    _impl->_bDataStart = _impl->_bInData = true;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        if (_ncol > 1 && !((*i)->flags() & AnalysisDataModuleInterface::efAllowMulticolumn))
        {
            GMX_ERROR(eeInvalidValue,
                      "Data module not compatible with data object properties");
        }
        int rc = (*i)->dataStarted(this);
        if (rc != 0)
        {
            return rc;
        }
    }
    return 0;
}


int
AbstractAnalysisData::notifyFrameStart(real x, real dx) const
{
    assert(_impl->_bInData && !_impl->_bInFrame);
    _impl->_bInFrame = true;
    _impl->_currx  = x;
    _impl->_currdx = dx;

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        int rc = (*i)->frameStarted(x, dx);
        if (rc != 0)
        {
            return rc;
        }
    }
    return 0;
}


int
AbstractAnalysisData::notifyPointsAdd(int firstcol, int n,
                                      const real *y, const real *dy,
                                      const bool *present) const
{
    assert(_impl->_bInData && _impl->_bInFrame);
    assert(firstcol >= 0 && n > 0 && firstcol + n <= _ncol);
    if (present && !_impl->_bAllowMissing)
    {
        for (int i = 0; i < n; ++i)
        {
            if (!present[i])
            {
                GMX_ERROR(eeInvalidValue,
                          "Missing data not supported by a module");
            }
        }
    }

    Impl::ModuleList::const_iterator i;
    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        int rc = (*i)->pointsAdded(_impl->_currx, _impl->_currdx, firstcol, n,
                                   y, dy, present);
        if (rc != 0)
        {
            return rc;
        }
    }
    return 0;
}


int
AbstractAnalysisData::notifyFrameFinish() const
{
    assert(_impl->_bInData && _impl->_bInFrame);
    _impl->_bInFrame = false;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        int rc = (*i)->frameFinished();
        if (rc != 0)
        {
            return rc;
        }
    }
    return 0;
}


int
AbstractAnalysisData::notifyDataFinish() const
{
    assert(_impl->_bInData && !_impl->_bInFrame);
    _impl->_bInData = false;

    Impl::ModuleList::const_iterator i;

    for (i = _impl->_modules.begin(); i != _impl->_modules.end(); ++i)
    {
        int rc = (*i)->dataFinished();
        if (rc != 0)
        {
            return rc;
        }
    }
    return 0;
}


/********************************************************************
 * AnalysisDataFrame
 */

AnalysisDataFrame::AnalysisDataFrame()
    : _y(NULL), _dy(NULL), _present(NULL)
{
}

AnalysisDataFrame::~AnalysisDataFrame()
{
    sfree(_y);
    sfree(_dy);
    sfree(_present);
}


void AnalysisDataFrame::allocate(int ncol)
{
    snew(_y, ncol);
    snew(_dy, ncol);
    snew(_present, ncol);
}


/********************************************************************
 * AbstractAnalysisDataStored::Impl
 */

AbstractAnalysisDataStored::Impl::Impl()
    : _nframes(0), _nalloc(0), _bStoreAll(false), _nextind(-1)
{
}


AbstractAnalysisDataStored::Impl::~Impl()
{
    FrameList::const_iterator i;
    for (i = _store.begin(); i != _store.end(); ++i)
    {
        delete *i;
    }
}


int
AbstractAnalysisDataStored::Impl::getStoreIndex(int index) const
{
    // Check that the requested index is available.
    if ((index < 0 && (-index > _nalloc || -index > _nframes))
        || index >= _nframes || (index >= 0 && index < _nframes - _nalloc))
    {
        return -1;
    }
    // Calculate the index into the storage array.
    if (index < 0)
    {
        index = _nextind + index;
        if (index < 0)
        {
            index += _nalloc;
        }
    }
    else if (_nframes > _nalloc)
    {
        index %= _nalloc;
    }
    return index;
}


/********************************************************************
 * AbstractAnalysisDataStored
 */

AbstractAnalysisDataStored::AbstractAnalysisDataStored()
    : _impl(new Impl())
{
}


AbstractAnalysisDataStored::~AbstractAnalysisDataStored()
{
    delete _impl;
}


int
AbstractAnalysisDataStored::frameCount() const
{
    return _impl->_nframes;
}


int
AbstractAnalysisDataStored::getDataWErr(int index, real *x, real *dx,
                                        const real **y, const real **dy,
                                        const bool **present) const
{
    index = _impl->getStoreIndex(index);
    if (index < 0)
    {
        return eedataDataNotAvailable;
    }

    // Retrieve the data.
    AnalysisDataFrame *fr = _impl->_store[index];
    if (x)
    {
        *x = fr->_x;
    }
    if (dx)
    {
        *dx = fr->_dx;
    }
    if (y)
    {
        *y = fr->_y;
    }
    if (dy)
    {
        *dy = fr->_dy;
    }
    if (present)
    {
        *present = fr->_present;
    }
    return 0;
}


int
AbstractAnalysisDataStored::requestStorage(int nframes)
{
    assert(nframes >= -1);
    if (nframes == 0)
    {
        return 0;
    }
    if (isMultipoint())
    {
        GMX_ERROR(eeNotImplemented,
                  "Storage of multipoint data not supported");
    }

    // Handle the case when everything needs to be stored.
    if (nframes == -1)
    {
        _impl->_bStoreAll = true;
        _impl->_nalloc = 1;
        return 0;
    }
    // Check whether an earier call has requested more storage.
    if (_impl->_bStoreAll || nframes < _impl->_nalloc)
    {
        return 0;
    }
    _impl->_nalloc = nframes;
    return 0;
}


void
AbstractAnalysisDataStored::setMultipoint(bool multipoint)
{
    assert(_impl->_nalloc == 0 || !multipoint);
    AbstractAnalysisData::setMultipoint(multipoint);
}


int
AbstractAnalysisDataStored::startDataStore()
{
    // We first notify any attached modules, because they also might request
    // some storage.
    int rc = notifyDataStart();
    if (rc != 0)
    {
        return rc;
    }

    int ncol = columnCount();

    // If any storage has been requested, preallocate it.
    if (_impl->_nalloc > 0)
    {
        _impl->_store.resize(_impl->_nalloc);
        for (int i = 0; i < _impl->_nalloc; ++i)
        {
            _impl->_store[i] = new AnalysisDataFrame();
            _impl->_store[i]->allocate(ncol);
        }
        _impl->_nextind = 0;
    }
    return 0;
}


int
AbstractAnalysisDataStored::startNextFrame(real x, real dx)
{
    // Start storing the frame if needed.
    if (_impl->_nalloc > 0)
    {
        if (_impl->_nextind >= _impl->_nalloc)
        {
            if (_impl->_bStoreAll)
            {
                int ncol = columnCount();

                _impl->_nalloc = _impl->_nextind + 1;
                _impl->_store.resize(_impl->_nalloc);
                for (int i = _impl->_nextind; i < _impl->_nalloc; ++i)
                {
                    _impl->_store[i] = new AnalysisDataFrame();
                    _impl->_store[i]->allocate(ncol);
                }
            }
            else
            {
                _impl->_nextind = 0;
            }
        }

        _impl->_store[_impl->_nextind]->_x  = x;
        _impl->_store[_impl->_nextind]->_dx = dx;
    }

    // Notify any modules.
    return notifyFrameStart(x, dx);
}


int
AbstractAnalysisDataStored::storeThisFrame(const real *y, const real *dy,
                                           const bool *present)
{
    int ncol = columnCount();

    // Store the values if required.
    if (_impl->_nextind >= 0)
    {
        AnalysisDataFrame *fr = _impl->_store[_impl->_nextind];

        for (int i = 0; i < ncol; ++i)
        {
            fr->_y[i] = y[i];
            if (dy)
            {
                fr->_dy[i] = dy[i];
            }
            if (present)
            {
                fr->_present[i] = present[i];
            }
        }
    }
    ++_impl->_nframes;

    // Notify modules of new data.
    int rc = notifyPointsAdd(0, ncol, y, dy, present);
    // The index needs to be incremented after the notifications to allow
    // the modules to use getData() properly.
    if (_impl->_nextind >= 0)
    {
        ++_impl->_nextind;
    }
    if (rc != 0)
    {
        return rc;
    }
    return notifyFrameFinish();
}


int
AbstractAnalysisDataStored::storeNextFrame(real x, real dx, const real *y,
                                           const real *dy, const bool *present)
{
    int rc = startNextFrame(x, dx);
    if (rc != 0)
    {
        return rc;
    }
    return storeThisFrame(y, dy, present);
}

} // namespace gmx
