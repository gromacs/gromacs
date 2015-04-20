/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Declares classes for accessing data frame information.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAFRAME_H
#define GMX_ANALYSISDATA_DATAFRAME_H

#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/flags.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/*! \brief
 * Value type for representing a single value in analysis data objects.
 *
 * Default copy constructor and assignment operator are used and work as
 * intended.
 *
 * Methods in this class do not throw.
 *
 * Non-const methods are provided for use within the library only; currently
 * it is not possible to access a non-const AnalysisDataValue through the
 * public interface.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataValue
{
    public:
        /*! \brief
         * Constructs an unset value.
         */
        AnalysisDataValue() : value_(0.0), error_(0.0) {}
        /*! \brief
         * Constructs a value object with the given value.
         *
         * The constructed object is marked as set and present.
         */
        explicit AnalysisDataValue(real value)
            : value_(value), error_(0.0)
        {
            flags_.set(efSet);
            flags_.set(efPresent);
        }

        /*! \brief
         * Direct access to the value.
         *
         * Assigning a value to this does not mark the value as set; setValue()
         * must be used for this.
         */
        real &value() { return value_; }
        /*! \brief
         * Direct access to the error estimate.
         *
         * Assigning a value to this does not mark the error estimate as set;
         * setValue() must be used for this.
         */
        real &error() { return error_; }
        //! Returns the value for this value.
        real value() const { return value_; }
        //! Returns the error estimate for this value, or zero if not set.
        real error() const { return error_; }
        /*! \brief
         * Returns whether this value has been set.
         *
         * If this method returns false, the return value of value() and
         * error() are undefined.
         */
        bool isSet() const { return flags_.test(efSet); }
        /*! \brief
         * Returns whether the error estimate for this value has been set.
         *
         * If this method returns false, but isSet() returns true, error()
         * returns zero.
         */
        bool hasError() const { return flags_.test(efErrorSet); }
        /*! \brief
         * Returns whether this value has been marked as present.
         *
         * If this method returns false, it is up to the source data to define
         * whether isSet() may return true.
         */
        bool isPresent() const { return flags_.test(efPresent); }

        //! Clears and unsets this value.
        void clear()
        {
            *this = AnalysisDataValue();
        }
        //! Sets this value.
        void setValue(real value, bool bPresent = true)
        {
            value_ = value;
            flags_.set(efSet);
            flags_.set(efPresent, bPresent);
        }
        //! Sets this value and its error estimate.
        void setValue(real value, real error, bool bPresent = true)
        {
            value_ = value;
            error_ = error;
            flags_.set(efSet);
            flags_.set(efErrorSet);
            flags_.set(efPresent, bPresent);
        }
        //! Set only error estimate for this value.
        void setError(real error)
        {
            error_ = error;
            flags_.set(efErrorSet);
        }

    private:
        //! Possible flags for \a flags_.
        enum Flag
        {
            efSet       = 1<<0, //!< Value has been set.
            efErrorSet  = 1<<1, //!< Error estimate has been set.
            efPresent   = 1<<2  //!< Value is set as present.
        };

        //! Value for this value.
        real                    value_;
        //! Error estimate for this value, zero if not set.
        real                    error_;
        //! Status flags for thise value.
        FlagsTemplate<Flag>     flags_;
};

//! Shorthand for reference to an array of data values.
typedef ConstArrayRef<AnalysisDataValue> AnalysisDataValuesRef;


/*! \brief
 * Value type for storing frame-level information for analysis data.
 *
 * Default copy constructor and assignment operator are used and work as
 * intended.
 * Typically new objects of this type are only constructed internally by the
 * library and in classes that are derived from AbstractAnalysisData.
 *
 * Methods in this class do not throw, but may contain asserts for incorrect
 * usage.
 *
 * Note that it is not possible to change the contents of an initialized
 * object, except by assigning a new object to replace it completely.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataFrameHeader
{
    public:
        /*! \brief
         * Constructs an invalid frame header.
         *
         * Return values of other methods than isValid() are unspecified for
         * the constructed object.
         */
        AnalysisDataFrameHeader();
        /*! \brief
         * Constructs a frame header from given values.
         *
         * \param[in] index  Index of the frame. Must be >= 0.
         * \param[in] x      x coordinate for the frame.
         * \param[in] dx     Error estimate for x.
         */
        AnalysisDataFrameHeader(int index, real x, real dx);

        /*! \brief
         * Returns whether the frame header corresponds to a valid frame.
         *
         * If returns false, return values of other methods are not specified.
         */
        bool isValid() const
        {
            return index_ >= 0;
        }
        /*! \brief
         * Returns zero-based index of the frame.
         *
         * The return value is >= 0 for valid frames.
         * Should not be called for invalid frames.
         */
        int index() const
        {
            GMX_ASSERT(isValid(), "Tried to access invalid frame header");
            return index_;
        }
        /*! \brief
         * Returns the x coordinate for the frame.
         *
         * Should not be called for invalid frames.
         */
        real x() const
        {
            GMX_ASSERT(isValid(), "Tried to access invalid frame header");
            return x_;
        }
        /*! \brief
         * Returns error in the x coordinate for the frame (if applicable).
         *
         * All data do not provide error estimates.
         * Typically returns zero in those cases.
         *
         * Should not be called for invalid frames.
         */
        real dx() const
        {
            GMX_ASSERT(isValid(), "Tried to access invalid frame header");
            return dx_;
        }

    private:
        int                     index_;
        real                    x_;
        real                    dx_;
};


/*! \cond libinternal */
/*! \libinternal \brief
 * Value type for internal indexing of point sets.
 *
 * This class contains the necessary data to split an array of
 * AnalysisDataValue objects into point sets.  It is always specified in the
 * context of an array of AnalysisDataValues: the point set specified by this
 * class contains valueCount() values, starting from the array index
 * valueOffset().
 * The value at location valueOffset() corresponds to column firstColumn().
 * It is not necessary for code using the analysis data framework to know of
 * this class, but it is declared in a public header to allow using it in other
 * types.
 *
 * Default copy constructor and assignment operator are used and work as
 * intended.
 * Typically new objects of this type are only constructed internally by the
 * library and in classes that are derived from AbstractAnalysisData.
 *
 * Methods in this class do not throw, but may contain asserts for incorrect
 * usage.
 *
 * Note that it is not possible to change the contents of an initialized
 * object, except by assigning a new object to replace it completely.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataPointSetInfo
{
    public:
        //! Construct point set data object with the given values.
        AnalysisDataPointSetInfo(int valueOffset, int valueCount,
                                 int dataSetIndex, int firstColumn)
            : valueOffset_(valueOffset), valueCount_(valueCount),
              dataSetIndex_(dataSetIndex), firstColumn_(firstColumn)
        {
            GMX_ASSERT(valueOffset  >= 0, "Negative value offsets are invalid");
            GMX_ASSERT(valueCount   >= 0, "Negative value counts are invalid");
            GMX_ASSERT(dataSetIndex >= 0, "Negative data set indices are invalid");
            GMX_ASSERT(firstColumn  >= 0, "Negative column indices are invalid");
        }

        //! Returns the offset of the first value in the referenced value array.
        int valueOffset() const { return valueOffset_; }
        //! Returns the number of values in this point set.
        int valueCount() const { return valueCount_; }
        //! Returns the data set index for this point set.
        int dataSetIndex() const { return dataSetIndex_; }
        //! Returns the index of the first column in this point set.
        int firstColumn() const { return firstColumn_; }

    private:
        int                     valueOffset_;
        int                     valueCount_;
        int                     dataSetIndex_;
        int                     firstColumn_;
};

//! Shorthand for reference to an array of point set data objects.
typedef ConstArrayRef<AnalysisDataPointSetInfo> AnalysisDataPointSetInfosRef;

//! \endcond


/*! \brief
 * Value type wrapper for non-mutable access to a set of data column values.
 *
 * Default copy constructor and assignment operator are used and work as
 * intended.
 * Typically new objects of this type are only constructed internally by the
 * library and in classes that are derived from AbstractAnalysisData.
 *
 * Methods in this class do not throw, but may contain asserts for incorrect
 * usage.
 *
 * The design of the interfaces is such that all objects of this type should be
 * valid, i.e., header().isValid() should always return true.
 *
 * Note that it is not possible to change the contents of an initialized
 * object, except by assigning a new object to replace it completely.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataPointSetRef
{
    public:
        /*! \brief
         * Constructs a point set reference from given values.
         *
         * \param[in] header       Header for the frame.
         * \param[in] pointSetInfo Information about the point set.
         * \param[in] values       Values for each column.
         *
         * The first element of the point set should be found from \p values
         * using the offset in \p pointSetInfo.
         */
        AnalysisDataPointSetRef(const AnalysisDataFrameHeader  &header,
                                const AnalysisDataPointSetInfo &pointSetInfo,
                                const AnalysisDataValuesRef    &values);
        /*! \brief
         * Constructs a point set reference from given values.
         *
         * \param[in] header      Header for the frame.
         * \param[in] values      Values for each column.
         *
         * The first element in \p values should correspond to the first
         * column.
         */
        AnalysisDataPointSetRef(const AnalysisDataFrameHeader        &header,
                                const std::vector<AnalysisDataValue> &values);
        /*! \brief
         * Constructs a point set reference to a subset of columns.
         *
         * \param[in] points      Point set to use as source.
         * \param[in] firstColumn First column index to include.
         * \param[in] columnCount Number of columns to include.
         *
         * Creates a point set that contains \p columnCount columns starting
         * from \p firstColumn from \p points, or a subset if all requested
         * columns are not present in \p points.  If the requested column range
         * and the range in \p points do not intersect, the result has
         * columnCount() == 0.
         *
         * \p firstColumn is relative to the whole data set, i.e., not relative
         * to points.firstColumn().
         *
         * Mainly intended for internal use.
         */
        AnalysisDataPointSetRef(const AnalysisDataPointSetRef &points,
                                int firstColumn, int columnCount);

        /*! \brief
         * Returns the frame header for the frame of this point set.
         */
        const AnalysisDataFrameHeader &header() const
        {
            return header_;
        }
        //! \copydoc AnalysisDataFrameHeader::index()
        int frameIndex() const
        {
            return header_.index();
        }
        //! \copydoc AnalysisDataFrameHeader::x()
        real x() const
        {
            return header_.x();
        }
        //! \copydoc AnalysisDataFrameHeader::dx()
        real dx() const
        {
            return header_.dx();
        }
        //! Returns zero-based index of the dataset that this set is part of.
        int dataSetIndex() const
        {
            return dataSetIndex_;
        }
        //! Returns zero-based index of the first column included in this set.
        int firstColumn() const
        {
            return firstColumn_;
        }
        //! Returns the number of columns included in this set.
        int columnCount() const
        {
            return values().size();
        }
        //! Returns zero-based index of the last column included in this set (inclusive).
        int lastColumn() const
        {
            return firstColumn_ + columnCount() - 1;
        }
        /*! \brief
         * Returns reference container for all values.
         *
         * First value in the returned container corresponds to firstColumn().
         */
        const AnalysisDataValuesRef &values() const
        {
            return values_;
        }
        /*! \brief
         * Returns data value for a column in this set.
         *
         * \param[in] i  Zero-based column index relative to firstColumn().
         *     Should be >= 0 and < columnCount().
         */
        real y(int i) const
        {
            GMX_ASSERT(i >= 0 && i < columnCount(), "Out of range data access");
            return values()[i].value();
        }
        /*! \brief
         * Returns error estimate for a column in this set if applicable.
         *
         * \param[in] i  Zero-based column index relative to firstColumn().
         *     Should be >= 0 and < columnCount().
         *
         * Currently, this method returns zero if the source data does not
         * specify errors.
         */
        real dy(int i) const
        {
            GMX_ASSERT(i >= 0 && i < columnCount(), "Out of range data access");
            return values()[i].error();
        }
        /*! \brief
         * Returns whether a column is present in this set.
         *
         * \param[in] i  Zero-based column index relative to firstColumn().
         *     Should be >= 0 and < columnCount().
         *
         * If present(i) returns false, it is depends on the source data
         * whether y(i) and/or dy(i) are defined.
         */
        bool present(int i) const
        {
            GMX_ASSERT(i >= 0 && i < columnCount(), "Out of range data access");
            return values()[i].isPresent();
        }
        /*! \brief
         * Returns true if all points in this point set are present.
         *
         * That is, if present() would return true for all points.
         */
        bool allPresent() const;

    private:
        AnalysisDataFrameHeader header_;
        int                     dataSetIndex_;
        int                     firstColumn_;
        AnalysisDataValuesRef   values_;
};


/*! \brief
 * Value type wrapper for non-mutable access to a data frame.
 *
 * Default copy constructor and assignment operator are used and work as
 * intended.
 * Typically new objects of this type are only constructed internally by the
 * library and in classes that are derived from AbstractAnalysisData.
 *
 * Methods in this class do not throw, but may contain asserts for incorrect
 * usage.
 *
 * Note that it is not possible to change the contents of an initialized
 * object, except by assigning a new object to replace it completely.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataFrameRef
{
    public:
        /*! \brief
         * Constructs an invalid frame reference.
         *
         * Return values of other methods than isValid() are unspecified for
         * the constructed object.
         */
        AnalysisDataFrameRef();
        /*! \brief
         * Constructs a frame reference from given values.
         *
         * \param[in] header      Header for the frame.
         * \param[in] values      Values for each column.
         * \param[in] pointSets   Point set data.
         */
        AnalysisDataFrameRef(const AnalysisDataFrameHeader      &header,
                             const AnalysisDataValuesRef        &values,
                             const AnalysisDataPointSetInfosRef &pointSets);
        /*! \brief
         * Constructs a frame reference from given values.
         *
         * \param[in] header      Header for the frame.
         * \param[in] values      Values for each column.
         * \param[in] pointSets   Point set data.
         */
        AnalysisDataFrameRef(const AnalysisDataFrameHeader               &header,
                             const std::vector<AnalysisDataValue>        &values,
                             const std::vector<AnalysisDataPointSetInfo> &pointSets);
        /*! \brief
         * Constructs a frame reference to a subset of columns.
         *
         * \param[in] frame       Frame to use as source.
         * \param[in] firstColumn First column index to include.
         * \param[in] columnCount Number of columns to include.
         *
         * Creates a frame reference that contains \p columnCount columns
         * starting from \p firstColumn from \p frame, or a subset if all
         * requested columns are not present in \p frame.
         *
         * Mainly intended for internal use.
         */
        AnalysisDataFrameRef(const AnalysisDataFrameRef &frame,
                             int firstColumn, int columnCount);

        /*! \brief
         * Returns whether the object refers to a valid frame.
         *
         * If returns false, return values of other methods are not specified.
         */
        bool isValid() const
        {
            return header().isValid();
        }
        //! Returns the header for this frame.
        const AnalysisDataFrameHeader &header() const
        {
            return header_;
        }
        //! \copydoc AnalysisDataFrameHeader::index()
        int frameIndex() const
        {
            return header().index();
        }
        //! \copydoc AnalysisDataFrameHeader::x()
        real x() const
        {
            return header().x();
        }
        //! \copydoc AnalysisDataFrameHeader::dx()
        real dx() const
        {
            return header().dx();
        }
        /*! \brief
         * Returns the number of point sets for this frame.
         *
         * Returns zero for an invalid frame.
         */
        int pointSetCount() const
        {
            return pointSets_.size();
        }
        /*! \brief
         * Returns point set reference for a given point set.
         *
         * Should not be called for invalid frames.
         */
        AnalysisDataPointSetRef pointSet(int index) const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            GMX_ASSERT(index >= 0 && index < pointSetCount(),
                       "Out of range data access");
            return AnalysisDataPointSetRef(header_, pointSets_[index], values_);
        }
        /*! \brief
         * Convenience method for accessing a column value in simple data.
         *
         * \copydetails AnalysisDataPointSetRef::y()
         */
        real y(int i) const
        {
            return singleColumnValue(i).value();
        }
        /*! \brief
         * Convenience method for accessing error for a column value in simple
         * data.
         *
         * \copydetails AnalysisDataPointSetRef::dy()
         */
        real dy(int i) const
        {
            return singleColumnValue(i).error();
        }
        /*! \brief
         * Convenience method for accessing present status for a column in
         * simple data.
         *
         * \copydetails AnalysisDataPointSetRef::present()
         */
        bool present(int i) const
        {
            return singleColumnValue(i).isPresent();
        }
        /*! \brief
         * Returns true if all points in this frame are present.
         */
        bool allPresent() const;

    private:
        //! Helper method for accessing single columns in simple data.
        const AnalysisDataValue &singleColumnValue(int i) const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            GMX_ASSERT(pointSets_.size() == 1U && pointSets_[0].firstColumn() == 0,
                       "Convenience method not available for multiple point sets");
            GMX_ASSERT(i >= 0 && i < static_cast<int>(values_.size()),
                       "Out of range data access");
            return values_[i];
        }

        AnalysisDataFrameHeader      header_;
        AnalysisDataValuesRef        values_;
        AnalysisDataPointSetInfosRef pointSets_;
};

} // namespace gmx

#endif
