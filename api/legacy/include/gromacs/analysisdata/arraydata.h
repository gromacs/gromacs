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
 * The x axis defaults to one with uniform spacing between values,
 * which caters to a typical case of e.g. a time series equally spaced
 * in time across a series of trajectory frames. However, non-uniform
 * cases are supported via setXAxisValue(), which can be useful when
 * handling other kinds of data.
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
    ~AbstractAnalysisArrayData() override;

    size_t frameCount() const override { return bReady_ ? rowCount_ : 0; }

    /*! \brief
     * Returns the number of rows in the data array.
     *
     * This function is identical to frameCount(), except that frameCount()
     * returns 0 before valuesReady() has been called.
     */
    size_t rowCount() const { return rowCount_; }
    //! Returns true if values have been allocated.
    bool isAllocated() const { return !value_.empty(); }
    //! Returns the step between frame x values when the x axis is uniform.
    real xstep() const
    {
        GMX_ASSERT(bUniformX_, "Accessing x step for non-uniform data");
        return xstep_;
    }
    //! Returns the x value of a row.
    real xvalue(size_t row) const
    {
        GMX_ASSERT(row < rowCount(), "Row index out of range");
        return xvalue_[row];
    }
    //! Returns a given array element.
    const AnalysisDataValue& value(size_t row, size_t col) const
    {
        GMX_ASSERT(row < rowCount(), "Row index out of range");
        GMX_ASSERT(col < columnCount(), "Column index out of range");
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
    void setColumnCount(size_t ncols);
    /*! \brief
     * Sets the number of rows in the data array.
     *
     * \param[in] rowCount  Number of rows in the data.
     *
     * Cannot be called after allocateValues().
     *
     * Cannot be called after setXAxisValues() made a non-uniform X
     * axis unless \c ncols equals the largest such X-axis value
     * previously set.
     *
     * Does not throw.
     */
    void setRowCount(size_t rowCount);
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
     * Afterwards, the X axis is uniform.
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
     * Afterwards, the X axis is non-uniform. Can be used to adjust
     * the values of an X axis created with setXAxis().
     *
     * The row count is never changed, and might need to be managed
     * explicitly with setRowCount() if needed.
     *
     * \param[in] row    Row/frame for which to set the value.
     * \param[in] value  x value for the frame specified by \p row.
     *
     * When the row count is already set, \c row must be in range.
     *
     * Must not be called after valuesReady().
     *
     * Does not throw.
     */
    void setXAxisValue(size_t row, real value);
    //! Returns a reference to a given array element.
    AnalysisDataValue& value(size_t row, size_t col)
    {
        GMX_ASSERT(row < rowCount(), "Row index out of range");
        GMX_ASSERT(col < columnCount(), "Column index out of range");
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
    static void copyContents(const AbstractAnalysisArrayData* src, AbstractAnalysisArrayData* dest);

private:
    AnalysisDataFrameRef tryGetDataFrameInternal(size_t index) const override;
    bool                 requestStorageInternal(size_t nframes) override;

    //! The number of rows
    size_t                   rowCount_;
    AnalysisDataPointSetInfo pointSetInfo_;
    //! The values of the columns of data
    std::vector<AnalysisDataValue> value_;
    //! The values of the X axis
    std::vector<real> xvalue_;
    //! Starting value for a uniform X axis
    real xstart_;
    //! Interval between x values for a uniform X axis
    real xstep_;
    //! Whether the X axis is uniform
    bool bUniformX_;
    //! Whether the data set is ready, ie. valuesReady() has been called
    bool bReady_;

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
    using AbstractAnalysisArrayData::allocateValues;
    using AbstractAnalysisArrayData::setColumnCount;
    using AbstractAnalysisArrayData::setRowCount;
    using AbstractAnalysisArrayData::setXAxis;
    using AbstractAnalysisArrayData::setXAxisValue;
    using AbstractAnalysisArrayData::value;
    using AbstractAnalysisArrayData::valuesReady;

    // Copy and assign disallowed by base.
};

} // namespace gmx

#endif
