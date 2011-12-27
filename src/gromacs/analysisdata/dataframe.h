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
 * Declares classes for accessing data frame information.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAFRAME_H
#define GMX_ANALYSISDATA_DATAFRAME_H

#include <cstddef>

#include "../legacyheaders/types/simple.h"

#include "../fatalerror/gmxassert.h"

namespace gmx
{

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
 * valid, i.e., header().isValid() should always return true.  This is
 * currently not strictly enforced in the constructors because of an
 * implementation detail of AnalysisDataFrameRef, but this is subject to
 * change.
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
         * \param[in] index       Index of the frame. Must be >= 0.
         * \param[in] x           x coordinate for the frame.
         * \param[in] dx          Error estimate for x.
         * \param[in] firstColumn Zero-based index of the first column.
         *     Must be >= 0.
         * \param[in] columnCount Number of columns to include.
         * \param[in] y           Array of values for each column.
         *     Must not be NULL if columnCount > 0.
         * \param[in] dy          Array of error estimates for corresponding y.
         *     Can be NULL, in which case errors cannot be accessed.
         * \param[in] present     Array of flags giving presence of each point.
         *     Can be NULL, in which case all values are treated as present.
         *
         * Arrays \p y, \p dy and \p dy should all have \p columnCount
         * elements.  The first elements in these arrays should correspond to
         * \p firstColumn.
         */
        AnalysisDataPointSetRef(int index, real x, real dx,
                                int firstColumn, int columnCount,
                                const real *y, const real *dy,
                                const bool *present);
        /*! \brief
         * Constructs a point set reference from given values.
         *
         * \param[in] header      Header for the frame.
         * \param[in] firstColumn Zero-based index of the first column.
         *     Must be >= 0.
         * \param[in] columnCount Number of columns to include.
         * \param[in] y           Array of values for each column.
         *     Must not be NULL if columnCount > 0.
         * \param[in] dy          Array of error estimates for corresponding y.
         *     Can be NULL, in which case errors cannot be accessed.
         * \param[in] present     Array of flags giving presence of each point.
         *     Can be NULL, in which case all values are treated as present.
         *
         * Arrays \p y, \p dy and \p dy should all have \p columnCount
         * elements.  The first elements in these arrays should correspond to
         * \p firstColumn.
         */
        AnalysisDataPointSetRef(const AnalysisDataFrameHeader &header,
                                int firstColumn, int columnCount,
                                const real *y, const real *dy,
                                const bool *present);
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
        //! Returns zero-based index of the first column included in this set.
        int firstColumn() const
        {
            return firstColumn_;
        }
        //! Returns the number of columns included in this set.
        int columnCount() const
        {
            return columnCount_;
        }
        //! Returns zero-based index of the last column included in this set (inclusive).
        int lastColumn() const
        {
            return firstColumn_ + columnCount_ - 1;
        }
        /*! \brief
         * Returns data value for a column in this set.
         *
         * \param[in] i  Zero-based column index relative to firstColumn().
         *     Should be >= 0 and < columnCount().
         */
        real y(int i) const
        {
            GMX_ASSERT(i >= 0 && i < columnCount_, "Out of range data access");
            return y_[i];
        }
        /*! \brief
         * Returns error estimate for a column in this set if applicable.
         *
         * \param[in] i  Zero-based column index relative to firstColumn().
         *     Should be >= 0 and < columnCount().
         *
         * Currently, this method either asserts or returns zero if the source
         * data does not specify errors.
         */
        real dy(int i) const
        {
            GMX_ASSERT(dy_ != NULL, "Errors not present, but accessed");
            GMX_ASSERT(i >= 0 && i < columnCount_, "Out of range data access");
            return dy_[i];
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
            GMX_ASSERT(i >= 0 && i < columnCount_, "Out of range data access");
            return present_ == NULL || present_[i];
        }
        /*! \brief
         * Returns true if all points in this point set are present.
         *
         * That is, if present() would return true for all points.
         */
        bool allPresent() const;

    private:
        AnalysisDataFrameHeader header_;
        int                     firstColumn_;
        int                     columnCount_;
        const real             *y_;
        const real             *dy_;
        const bool             *present_;
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
 * \todo
 * Support for multipoint data.
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
         * \param[in] index       Index of the frame. Must be >= 0.
         * \param[in] x           x coordinate for the frame.
         * \param[in] dx          Error estimate for x.
         * \param[in] columnCount Number of columns to include. Must be >= 0.
         * \param[in] y           Array of values for each column.
         *     Must not be NULL if columnCount > 0.
         * \param[in] dy          Array of error estimates for corresponding y.
         *     Can be NULL, in which case errors cannot be accessed.
         * \param[in] present     Array of flags giving presence of each point.
         *     Can be NULL, in which case all values are treated as present.
         *
         * Arrays \p y, \p dy and \p dy should all have \p columnCount
         * elements.
         */
        AnalysisDataFrameRef(int index, real x, real dx,
                             int columnCount,
                             const real *y, const real *dy,
                             const bool *present);
        /*! \brief
         * Constructs a frame reference from given values.
         *
         * \param[in] header      Header for the frame.
         * \param[in] columnCount Number of columns to include.
         * \param[in] y           Array of values for each column.
         *     Must not be NULL if columnCount > 0.
         * \param[in] dy          Array of error estimates for corresponding y.
         *     Can be NULL, in which case errors cannot be accessed.
         * \param[in] present     Array of flags giving presence of each point.
         *     Can be NULL, in which case all values are treated as present.
         *
         * Arrays \p y, \p dy and \p dy should all have \p columnCount
         * elements.
         */
        AnalysisDataFrameRef(const AnalysisDataFrameHeader &header,
                             int columnCount,
                             const real *y, const real *dy,
                             const bool *present);
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
            return points_.header();
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
         * Returns point set reference to the column values of this frame.
         *
         * Should not be called for invalid frames.
         */
        const AnalysisDataPointSetRef &points() const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            return points_;
        }
        /*! \brief
         * Returns number of columns in this frame.
         *
         * Returns zero for an invalid frame.
         */
        int columnCount() const
        {
            return points().columnCount();
        }
        /*! \brief
         * Convenience method for accessing a column value.
         *
         * \copydetails AnalysisDataPointSetRef::y()
         */
        real y(int i) const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            return points().y(i);
        }
        /*! \brief
         * Convenience method for accessing error for a column value.
         *
         * \copydetails AnalysisDataPointSetRef::dy()
         */
        real dy(int i) const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            return points().dy(i);
        }
        /*! \brief
         * Convenience method for accessing present status for a column.
         *
         * \copydetails AnalysisDataPointSetRef::present()
         */
        bool present(int i) const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            return points().present(i);
        }
        /*! \brief
         * Returns true if all points in this frame are present.
         */
        bool allPresent() const
        {
            GMX_ASSERT(isValid(), "Invalid data frame accessed");
            return points().allPresent();
        }

    private:
        AnalysisDataPointSetRef points_;
};

} // namespace gmx

#endif
