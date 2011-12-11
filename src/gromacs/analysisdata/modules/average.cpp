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
 * Implements gmx::AnalysisDataAverageModule.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#include "gromacs/analysisdata/modules/average.h"

#include <cmath>

// Legacy include
#include "smalloc.h"

#include "gromacs/basicmath.h"

namespace gmx
{

AnalysisDataAverageModule::AnalysisDataAverageModule()
    : _nsamples(NULL)
{
    setColumnCount(2);
}


AnalysisDataAverageModule::~AnalysisDataAverageModule()
{
    sfree(_nsamples);
}


int
AnalysisDataAverageModule::flags() const
{
    return efAllowMultipoint | efAllowMulticolumn | efAllowMissing;
}


void
AnalysisDataAverageModule::dataStarted(AbstractAnalysisData *data)
{
    int nrows = data->columnCount();
    setRowCount(nrows);
    snew(_nsamples, nrows);
}


void
AnalysisDataAverageModule::frameStarted(real x, real dx)
{
}


void
AnalysisDataAverageModule::pointsAdded(real x, real dx, int firstcol, int n,
                                       const real *y, const real *dy,
                                       const bool *present)
{
    for (int i = 0; i < n; ++i)
    {
        if (!present || present[i])
        {
            value(firstcol + i, 0)  += y[i];
            value(firstcol + i, 1)  += y[i] * y[i];
            _nsamples[firstcol + i] += 1;
        }
    }
}


void
AnalysisDataAverageModule::frameFinished()
{
}


void
AnalysisDataAverageModule::dataFinished()
{
    for (int i = 0; i < rowCount(); ++i)
    {
        real ave = value(i, 0) / _nsamples[i];
        real std = sqrt(value(i, 1) / _nsamples[i] - ave * ave);
        setValue(i, 0, ave);
        setValue(i, 1, std);
    }
    valuesReady();
}


real
AnalysisDataAverageModule::average(int index) const
{
    return value(index, 0);
}


real
AnalysisDataAverageModule::stddev(int index) const
{
    return value(index, 1);
}

} // namespace gmx
