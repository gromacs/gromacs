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
/*! \file
 * \brief
 * Declares gmx::AbstractAnalysisArrayData and gmx::AnalysisArrayData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ARRAYDATA_H
#define GMX_ANALYSISDATA_ARRAYDATA_H

#include <cstddef>

#include "../fatalerror/gmxassert.h"

#include "abstractdata.h"

namespace gmx
{

/*! \brief
 * Abstract base class for data objects that present in-memory data.
 *
 * This class implements a subclass of AbstractAnalysisData that presents an
 * in-memory array through the AbstractAnalysisData interface.  Subclasses
 * should initialize the in-memory array through the provided protected member
 * functions.
 *
 * \ingroup module_analysisdata
 */
class AbstractAnalysisArrayData : public AbstractAnalysisData
{
    public:
        virtual ~AbstractAnalysisArrayData();

        virtual int frameCount() const;
        virtual bool getDataWErr(int index, real *x, real *dx,
                                 const real **y, const real **dy,
                                 const bool **present = 0) const;
        virtual bool requestStorage(int nframes = -1);

    protected:
        AbstractAnalysisArrayData();

        /*! \brief
         * Returns the number of rows in the data array.
         *
         * This function is identical to frameCount(), except that frameCount()
         * returns 0 before valuesReady() has been called.
         */
        int rowCount() const { return _nrows; }
        /*! \brief
         * Sets the number of columns in the data array.
         */
        void setColumnCount(int ncols);
        /*! \brief
         * Sets the number of rows in the data array.
         */
        void setRowCount(int nrows);
        //! Returns the x value of the first frame.
        real xstart() const { return _xstart; }
        //! Returns the step between frame x values.
        real xstep() const { return _xstep; }
        /*! \brief
         * Sets the values reported as x values for frames.
         */
        void setXAxis(real start, real step);
        //! Returns a reference to a given array element.
        real &value(int row, int col)
        {
            GMX_ASSERT(row >= 0 && row < _nrows, "Row index out of range");
            GMX_ASSERT(col >= 0 && col < columnCount(), "Column index out of range");
            GMX_ASSERT(_value != NULL, "Data array not allocated");
            return _value[row * columnCount() + col];
        }
        //! Returns a given array element.
        const real &value(int row, int col) const
        {
            GMX_ASSERT(row >= 0 && row < _nrows, "Row index out of range");
            GMX_ASSERT(col >= 0 && col < columnCount(), "Column index out of range");
            GMX_ASSERT(_value != NULL, "Data array not allocated");
            return _value[row * columnCount() + col];
        }
        /*! \brief
         * Sets the value of an element in the array.
         */
        void setValue(int row, int col, real val)
        {
            value(row, col) = val;
        }
        /*! \brief
         * Notifies modules of the data.
         *
         * This function should be called once the values in the array
         * have been initialized. The values should not be changed after this
         * function has been called.
         */
        void valuesReady();

        /*! \brief
         * Copies the contents into a new object.
         *
         * \param[in]     src  Object to copy data from.
         * \param[in,out] dest Empty array data object to copy data to.
         *
         * \p dest should not have previous contents.
         */
        static void copyContents(const AbstractAnalysisArrayData *src,
                                 AbstractAnalysisArrayData *dest);

    private:
        int                  _nrows;
        real                *_value;
        real                 _xstart;
        real                 _xstep;
        bool                 _bReady;

        // Copy and assign disallowed by base.
};

/*! \brief
 * Simple in-memory data array.
 *
 * This class simply exposes the protected functions of
 * AbstractAnalysisArrayData to allow construction of simple in-memory data
 * arrays for input into data modules.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisArrayData : public AbstractAnalysisArrayData
{
    public:
        AnalysisArrayData() {}

        using AbstractAnalysisArrayData::rowCount;
        using AbstractAnalysisArrayData::setColumnCount;
        using AbstractAnalysisArrayData::setRowCount;
        using AbstractAnalysisArrayData::xstart;
        using AbstractAnalysisArrayData::xstep;
        using AbstractAnalysisArrayData::setXAxis;
        using AbstractAnalysisArrayData::value;
        using AbstractAnalysisArrayData::setValue;
        using AbstractAnalysisArrayData::valuesReady;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
