/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
 * \brief
 * Declares gmx::AbstractAnalysisArrayData and gmx::AnalysisArrayData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ARRAYDATA_H
#define GMX_ANALYSISDATA_ARRAYDATA_H

#include <vector>

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief
 * Abstract base class for data objects that present in-memory data.
 *
 * This class implements a subclass of AbstractAnalysisData that presents an
 * in-memory array through the AbstractAnalysisData interface.  Subclasses
 * should initialize the in-memory array through the provided protected member
 * functions.  This class provides public accessor methods for read access to
 * the data.
 *
 * Public accessor methods in this class do not throw, but assert if data is
 * accessed before it is available.
 *
 * \todo
 * Add support for multiple data sets.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AbstractAnalysisArrayData : public AbstractAnalysisData
{
    public:
        virtual ~AbstractAnalysisArrayData();

        virtual int frameCount() const
        {
            return bReady_ ? rowCount_ : 0;
        }

        /*! \brief
         * Returns the number of rows in the data array.
         *
         * This function is identical to frameCount(), except that frameCount()
         * returns 0 before valuesReady() has been called.
         */
        int rowCount() const { return rowCount_; }
        //! Returns true if values have been allocated.
        bool isAllocated() const { return !value_.empty(); }
        //! Returns the x value of the first frame.
        real xstart() const { return xvalue_[0]; }
        //! Returns the step between frame x values.
        real xstep() const
        {
            GMX_ASSERT(bUniformX_, "Accessing x step for non-uniform data");
            return xstep_;
        }
        //! Returns the x value of a row.
        real xvalue(int row) const
        {
            GMX_ASSERT(row >= 0 && row < rowCount(), "Row index out of range");
            return xvalue_[row];
        }
        //! Returns a given array element.
        const AnalysisDataValue &value(int row, int col) const
        {
            GMX_ASSERT(row >= 0 && row < rowCount(), "Row index out of range");
            GMX_ASSERT(col >= 0 && col < columnCount(), "Column index out of range");
            GMX_ASSERT(isAllocated(), "Data array not allocated");
            return value_[row * columnCount() + col];
        }

    protected:
        /*! \brief
         * Initializes an empty array data object.
         *
         * \throws std::bad_alloc if out of memory.
         */
        AbstractAnalysisArrayData();

        /*! \brief
         * Sets the number of columns in the data array.
         *
         * \param[in] ncols  Number of columns in the data.
         *
         * Cannot be called after allocateValues().
         *
         * See AbstractAnalysisData::setColumnCount() for exception behavior.
         */
        void setColumnCount(int ncols);
        /*! \brief
         * Sets the number of rows in the data array.
         *
         * \param[in] rowCount  Number of rows in the data.
         *
         * Cannot be called after allocateValues().
         *
         * Does not throw.
         */
        void setRowCount(int rowCount);
        /*! \brief
         * Allocates memory for the values.
         *
         * \throws std::bad_alloc if memory allocation fails.
         *
         * setColumnCount() and setRowCount() must have been called.
         *
         * Strong exception safety guarantee.
         */
        void allocateValues();
        /*! \brief
         * Sets the values reported as x values for frames.
         *
         * \param[in] start  x value for the first frame.
         * \param[in] step   Step between x values of successive frames.
         *
         * Must not be called after valuesReady().
         * Any values set with setXAxisValue() are overwritten.
         *
         * Does not throw.
         */
        void setXAxis(real start, real step);
        /*! \brief
         * Sets a single value reported as x value for frames.
         *
         * \param[in] row    Row/frame for which to set the value.
         * \param[in] value  x value for the frame specified by \p row.
         *
         * Must not be called after valuesReady().
         *
         * Does not throw.
         */
        void setXAxisValue(int row, real value);
        //! Returns a reference to a given array element.
        AnalysisDataValue &value(int row, int col)
        {
            GMX_ASSERT(row >= 0 && row < rowCount(), "Row index out of range");
            GMX_ASSERT(col >= 0 && col < columnCount(), "Column index out of range");
            GMX_ASSERT(isAllocated(), "Data array not allocated");
            return value_[row * columnCount() + col];
        }
        /*! \brief
         * Notifies modules of the data.
         *
         * \throws    unspecified Any exception thrown by attached data modules
         *      in data notification methods.
         *
         * This function should be called once the values in the array
         * have been initialized.  The values should not be changed after this
         * function has been called.
         */
        void valuesReady();

        /*! \brief
         * Copies the contents into a new object.
         *
         * \param[in]     src  Object to copy data from.
         * \param[in,out] dest Empty array data object to copy data to.
         * \throws std::bad_alloc if memory allocation for \p dest fails.
         *
         * \p dest should not have previous contents.
         */
        static void copyContents(const AbstractAnalysisArrayData *src,
                                 AbstractAnalysisArrayData       *dest);

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        int                            rowCount_;
        AnalysisDataPointSetInfo       pointSetInfo_;
        std::vector<AnalysisDataValue> value_;
        std::vector<real>              xvalue_;
        real                           xstep_;
        bool                           bUniformX_;
        bool                           bReady_;

        // Copy and assign disallowed by base.
};

/*! \brief
 * Simple in-memory data array.
 *
 * This class is a simple alternative to AnalysisData for in-memory data arrays
 * that are constructed in-place.
 *
 * Public accessor methods in this class do not throw, but assert if data is
 * accessed before it is available.
 *
 * \if libapi
 * This class exposes the protected functions of AbstractAnalysisArrayData for
 * users.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisArrayData : public AbstractAnalysisArrayData
{
    public:
        /*! \brief
         * Initializes an empty array data object.
         *
         * \throws std::bad_alloc if out of memory.
         */
        AnalysisArrayData() {}

        // TODO: These statements cause Doxygen to generate confusing
        // documentation.
        using AbstractAnalysisArrayData::setColumnCount;
        using AbstractAnalysisArrayData::setRowCount;
        using AbstractAnalysisArrayData::allocateValues;
        using AbstractAnalysisArrayData::setXAxis;
        using AbstractAnalysisArrayData::setXAxisValue;
        using AbstractAnalysisArrayData::value;
        using AbstractAnalysisArrayData::valuesReady;

        // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
