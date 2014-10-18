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
 * Declares gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ABSTRACTDATA_H
#define GMX_ANALYSISDATA_ABSTRACTDATA_H

#include <boost/shared_ptr.hpp>

#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class AnalysisDataModuleInterface;
class AnalysisDataModuleManager;
class AnalysisDataFrameHeader;
class AnalysisDataFrameRef;
class AnalysisDataPointSetRef;

//! Smart pointer for managing a generic analysis data module.
typedef boost::shared_ptr<AnalysisDataModuleInterface> AnalysisDataModulePointer;

/*! \brief
 * Abstract base class for all objects that provide data.
 *
 * The public interface includes methods for querying the data (isMultipoint(),
 * dataSetCount(), columnCount(), frameCount(), tryGetDataFrame(),
 * getDataFrame(), requestStorage()) and methods for using modules for
 * processing the data (addModule(), addColumnModule(), applyModule()).
 *
 * Notice that even for non-const objects, the interface does not provide any
 * means of altering the data.  It is only possible to add modules, making it
 * relatively safe to return a non-const pointer of this type pointing to an
 * internal data structure without worrying about possible modifications of the
 * data.
 *
 * \if libapi
 * This class also provides protected methods for use in derived classes.
 * The properties returned by isMultipoint(), dataSetCount(), and columnCount()
 * must be set using setMultipoint(), setDataSetCount(), and setColumnCount().
 * notify*() methods in the AnalysisDataModuleManager returned by
 * moduleManager() must be used to report when data becomes available for
 * modules to process it.
 * There are also three pure virtual methods that need to be implemented to
 * provide access to stored data: one public (frameCount()) and two protected
 * ones (requestStorageInternal() and tryGetDataFrameInternal()).
 *
 * It is up to subclasses to ensure that the virtual methods and the
 * notifications in AnalysisDataModuleManager are called in a
 * correct sequence (the methods will assert in most incorrect use cases), and
 * that the data provided through the public interface matches that passed to
 * the modules with the notify methods.
 * Helper class AnalysisDataStorage provides a default implementation for
 * storing data (calls to the pure virtual methods can simply be forwarded to
 * appropriate methods in the helper class), and takes care of correctly
 * calling the notification methods when new data is added to the storage.
 * In most cases, it should be used to implement the derived classes.
 * \endif
 *
 * Currently, it is not possible to continue using the data object if an
 * attached module throws an exception during data processing; it is only safe
 * to destroy such data object.
 *
 * \todo
 * Improve the exception-handling semantics.  In most cases, it doesn't make
 * much sense to continue data processing after one module fails, but having
 * the alternative would not hurt.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AbstractAnalysisData
{
    public:
        virtual ~AbstractAnalysisData();

        /*! \brief
         * Whether the data can have multiple points in the same column
         * in the same frame.
         *
         * \returns \c true if multiple points in the same column are
         *     allowed within a single frame.
         *
         * This kind of data can appear in many histogramming applications
         * (e.g., RDFs), where each trajectory frame has several data points
         * (possibly a different number for each frame). The current interface
         * doesn't support storing such data, but this should rarely be
         * necessary.
         *
         * The returned value does not change after modules have been notified
         * of data start.
         * \if libapi
         * Derived classes can change the type by calling setMultipoint()
         * subject to the above restriction.
         * If this is not done, the function always returns false.
         * \endif
         *
         * Does not throw.
         */
        bool isMultipoint() const;
        /*! \brief
         * Returns the number of data sets in the data object.
         *
         * \returns The number of data sets in the data.
         *
         * If the number is not yet known, returns 0.
         * The returned value does not change after modules have been notified
         * of data start, but may change multiple times before that, depending
         * on the actual data class.
         * \if libapi
         * Derived classes should set the number of columns with
         * setDataSetCount(), within the above limitations.
         * \endif
         *
         * Does not throw.
         */
        int dataSetCount() const;
        /*! \brief
         * Returns the number of columns in a data set.
         *
         * \param[in] dataSet Zero-based index of the data set to query.
         * \returns The number of columns in the data.
         *
         * If the number of columns is not yet known, returns 0.
         * The returned value does not change after modules have been notified
         * of data start, but may change multiple times before that, depending
         * on the actual data class.
         * \if libapi
         * Derived classes should set the number of columns with
         * setColumnCount(), within the above limitations.
         * \endif
         *
         * Does not throw.
         */
        int columnCount(int dataSet) const;
        /*! \brief
         * Returns the number of columns in the data.
         *
         * \returns The number of columns in the data.
         *
         * This is a convenience method for data objects with a single data set.
         * Can only be called if dataSetCount() == 1.
         *
         * Does not throw.
         *
         * \see columnCount(int)
         */
        int columnCount() const;
        /*! \brief
         * Returns the total number of frames in the data.
         *
         * \returns The total number of frames in the data.
         *
         * This function returns the number of frames that the object has
         * produced.  If requestStorage() has been successfully called,
         * tryGetDataframe() or getDataFrame() can be used to access some or
         * all of these frames.
         *
         * Does not throw.
         *
         * \if libapi
         * Derived classes should implement this to return the number of
         * frames.  The frame count should not be incremented before
         * tryGetDataFrameInternal() can return the new frame.
         * The frame count must be incremented before
         * AnalysisDataModuleManager::notifyFrameFinish() is called.
         * \endif
         */
        virtual int frameCount() const = 0;
        /*! \brief
         * Access stored data.
         *
         * \param[in] index  Zero-based frame index to access.
         * \returns   Frame reference to frame \p index, or an invalid
         *      reference if no such frame is available.
         *
         * Does not throw.  Failure to access a frame with the given index is
         * indicated through the return value.  Negative \p index is allowed,
         * and will always result in an invalid reference being returned.
         *
         * \see requestStorage()
         * \see getDataFrame()
         */
        AnalysisDataFrameRef tryGetDataFrame(int index) const;
        /*! \brief
         * Access stored data.
         *
         * \param[in] index  Zero-based frame index to access.
         * \returns   Frame reference to frame \p index.
         * \throws    APIError if the requested frame is not accessible.
         *
         * If the data is not certainly available, use tryGetDataFrame().
         *
         * \see requestStorage()
         * \see tryGetDataFrame()
         */
        AnalysisDataFrameRef getDataFrame(int index) const;
        /*! \brief
         * Request storage of frames.
         *
         * \param[in] nframes  Request storing at least \c nframes previous
         *     frames (-1 = request storing all). Must be >= -1.
         * \returns true if the request could be satisfied.
         *
         * If called multiple times, the largest request is honored.
         *
         * Does not throw.  Failure to honor the request is indicated through
         * the return value.
         *
         * \see getDataFrame()
         * \see tryGetDataFrame()
         */
        bool requestStorage(int nframes);

        /*! \brief
         * Adds a module to process the data.
         *
         * \param     module  Module to add.
         * \throws    std::bad_alloc if out of memory.
         * \throws    APIError if
         *      - \p module is not compatible with the data object
         *      - data has already been added to the data object and everything
         *        is not available through getDataFrame().
         * \throws    unspecified Any exception thrown by \p module in its
         *      notification methods (if data has been added).
         *
         * If data has already been added to the data, the new module
         * immediately processes all existing data.  APIError is thrown
         * if all data is not available through getDataFrame().
         *
         * The caller can keep a copy of the module pointer if it requires
         * later access to the module.
         *
         * If the method throws, the state of the data object is not changed.
         * The state of the data module is indeterminate.
         */
        void addModule(AnalysisDataModulePointer module);
        /*! \brief
         * Adds a module that processes only a subset of the columns.
         *
         * \param[in] col     First column.
         * \param[in] span    Number of columns.
         * \param     module  Module to add.
         *
         * Throws in the same situations as addModule().
         *
         * Currently, all data sets are filtered using the same column mask.
         *
         * \todo
         * This method doesn't currently work in all cases with multipoint
         * data or with multiple data sets.  In particular, if the added module
         * requests storage and uses getDataFrame(), it will behave
         * unpredictably (most likely asserts).
         *
         * \todo
         * Generalize this method to multiple data sets (e.g., for adding
         * modules that only process a single data set).
         *
         * \see addModule()
         */
        void addColumnModule(int col, int span, AnalysisDataModulePointer module);
        /*! \brief
         * Applies a module to process data that is ready.
         *
         * \param     module  Module to apply.
         * \throws    APIError in same situations as addModule().
         * \throws    unspecified Any exception thrown by \p module in its
         *      notification methods.
         *
         * This function works as addModule(), except that it does not keep a
         * reference to \p module within the data object after it returns.
         * Also, it can only be called after the data is ready, and only if
         * getDataFrame() gives access to all of the data.
         * It is provided for additional flexibility in postprocessing
         * in-memory data.
         *
         * \todo
         * Currently, this method may not work correctly if \p module requests
         * storage (addModule() has the same problem if called after data is
         * started).
         */
        void applyModule(AnalysisDataModuleInterface *module);

    protected:
        /*! \cond libapi */
        /*! \brief
         * Initializes a new analysis data object.
         *
         * \throws std::bad_alloc if out of memory.
         */
        AbstractAnalysisData();

        /*! \brief
         * Sets the number of data sets.
         *
         * \param[in] dataSetCount  Number of data sets (must be > 0).
         * \throws    std::bad_alloc if out of memory.
         * \throws    APIError if modules have been added that are not
         *      compatible with the new data set count.
         *
         * It not called, the data object has a single data set.  Can be called
         * only before AnalysisDataModuleManager::notifyDataStart().
         * Multiple calls are allowed before that point; the last call takes
         * effect.
         *
         * Strong exception safety.
         *
         * \see dataSetCount()
         */
        void setDataSetCount(int dataSetCount);
        /*! \brief
         * Sets the number of columns for a data set.
         *
         * \param[in] dataSet      Zero-based index of the data set.
         * \param[in] columnCount  Number of columns in \p dataSet (must be > 0).
         * \throws    APIError if modules have been added that are not
         *      compatible with the new column count.
         *
         * Must be called at least once for each data set before
         * AnalysisDataModuleManager::notifyDataStart().  Can be called only
         * before AnalysisDataModuleManager::notifyDataStart().
         * Multiple calls are allowed before that point; the last call takes
         * effect.
         *
         * Strong exception safety.
         *
         * \see columnCount()
         */
        void setColumnCount(int dataSet, int columnCount);
        /*! \brief
         * Sets whether the data has multiple points per column in a frame.
         *
         * \param[in] bMultipoint  Whether multiple points per column are
         *     possible.
         * \throws    APIError if modules have been added that are not
         *      compatible with the new setting.
         *
         * If not called, only a single point per column is allowed.  Can be
         * called only before AnalysisDataModuleManager::notifyDataStart().
         * Multiple calls are allowed before that point; the last call takes
         * effect.
         *
         * Strong exception safety.
         *
         * \see isMultipoint()
         */
        void setMultipoint(bool bMultipoint);

        /*! \brief
         * Implements access to data frames.
         *
         * \param[in] index  Zero-based frame index to access.
         * \returns   Frame reference to frame \p index, or an invalid
         *      reference if no such frame is available.
         *
         * Must not throw.  Failure to access a frame with the given index is
         * indicated through the return value.
         *
         * Code in derived classes can assume that \p index is non-negative and
         * less than frameCount().
         *
         * Derived classes can choose to return an invalid reference if
         * requestStorageInternal() has not been called at all, or if the frame
         * is too old (compared to the value given to requestStorageInternal()).
         *
         * This method is called internally by tryGetDataFrame() and
         * getDataFrame().
         *
         * \see AnalysisDataStorage
         */
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const = 0;
        /*! \brief
         * Implements storage requests.
         *
         * \param[in] nframes  Request storing at least \c nframes previous
         *     frames (-1 = request storing all). Will be either -1 or >0.
         * \returns   true if the request could be satisfied.
         *
         * Must not throw.  Failure to access a frame with the given index is
         * indicated through the return value.
         *
         * Derived classes should be prepared for any number of calls to this
         * method before notifyDataStart() is called (and during that call).
         *
         * This method is called internally by requestStorage().
         *
         * \see AnalysisDataStorage
         */
        virtual bool requestStorageInternal(int nframes) = 0;

        //! Returns the module manager to use for calling notification methods.
        AnalysisDataModuleManager       &moduleManager();
        //! Returns the module manager to use for calling notification methods.
        const AnalysisDataModuleManager &moduleManager() const;
        //! \endcond

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
