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
 * Declares gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ABSTRACTDATA_H
#define GMX_ANALYSISDATA_ABSTRACTDATA_H

#include "../legacyheaders/types/simple.h"

#include "../utility/common.h"

namespace gmx
{

class AnalysisDataModuleInterface;
class AnalysisDataFrameHeader;
class AnalysisDataFrameRef;
class AnalysisDataPointSetRef;
class AnalysisDataStorage;

/*! \brief
 * Abstract base class for all objects that provide data.
 *
 * The public interface includes functions for querying the data
 * (isMultipoint(), columnCount(), frameCount(), tryGetDataFrame(),
 * getDataFrame(), requestStorage()) and functions for using modules for
 * processing the data (addModule(), addColumnModule(), applyModule()).
 *
 * \if libapi
 * There are also protected functions for use in derived classes:
 * setting the properties returned by isMultipoint() and columnCount()
 * through setMultipoint() and setColumnCount(), and notifying attached
 * modules of added data.
 * \endif
 *
 * Notice that even for non-const objects, the interface does not provide
 * any means of altering the data. It is only possible to add modules,
 * making it relatively safe to return a non-const pointer of this type
 * pointing to an internal data structure without worrying about possible
 * modifications of the data.
 *
 * \if libapi
 * It is up to subclasses to ensure that the protected functions are called
 * in a correct sequence (the functions will assert in most incorrect use
 * cases), and that the data provided through the public interface matches
 * that passed to the modules with the notify methods.
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
         * \if libapi
         * Derived classes can change the type by calling setMultipoint().
         * If this is not done, the function always returns false.
         * \endif
         */
        bool isMultipoint() const { return _bMultiPoint; }
        /*! \brief
         * Returns the number of columns in the data.
         *
         * \returns The number of columns in the data.
         *
         * If the number of columns is yet known, returns 0.
         * \if libapi
         * Derived classes should set the number of columns with
         * setColumnCount().
         * \endif
         */
        int columnCount() const { return _ncol; }
        /*! \brief
         * Returns the total number of frames in the data.
         *
         * \returns The total number of frames in the data.
         *
         * This function returns the number of frames that the object has
         * produced. If requestStorage() has been successfully called,
         * getDataFrame() can be used to access some or all of these frames.
         */
        int frameCount() const;
        /*! \brief
         * Access stored data.
         *
         * \param[in] index  Zero-based frame index to access.
         * \retval    Frame reference to frame \p index, or an invalid
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
         * \retval    Frame reference to frame \p index.
         * \exception APIError if the requested frame is accessible.
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
         * \retval true if the request could be satisfied.
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
         * \param  module  Module to add.
         * \exception APIError if
         *      - \p module is not compatible with the data object
         *      - data has already been added to the data object and everything
         *        is not available through getDataFrame().
         *
         * If data has already been added to the module, the new module
         * immediately processes all existing data.  APIError is thrown
         * if all data is not available through getDataFrame().
         *
         * When this function is entered, the data object takes ownership of the
         * module, and automatically destructs it when the data object itself
         * is destroyed.
         *
         * \todo
         * Provide additional semantics that does not acquire ownership of the
         * data object.
         */
        void addModule(AnalysisDataModuleInterface *module);
        /*! \brief
         * Adds a module that processes only a subset of the columns.
         *
         * \param[in] col     First column.
         * \param[in] span    Number of columns.
         * \param     module  Module to add.
         *
         * \see addModule()
         */
        void addColumnModule(int col, int span, AnalysisDataModuleInterface *module);
        /*! \brief
         * Applies a module to process data that is ready.
         *
         * \param  module  Module to apply.
         *
         * This function works as addModule(), except that it does not take
         * ownership of \p module. Also, it can only be called after the data
         * is ready, and only if getDataFrame() gives access to all of the data.
         * It is provided for additional flexibility in postprocessing
         * in-memory data.
         */
        void applyModule(AnalysisDataModuleInterface *module);

    protected:
        /*! \cond libapi */
        AbstractAnalysisData();

        /*! \brief
         * Sets the number of columns.
         *
         * \param[in] ncol  Number of columns in the data (must be > 0).
         *
         * Can be called only before notifyDataStart(), otherwise asserts.
         * Multiple calls are only allowed if all of them occur before
         * addModule() has been called, otherwise asserts (a single call
         * can occur after addModule() if no calls have been made earlier).
         *
         * \see columnCount()
         */
        void setColumnCount(int ncol);
        /*! \brief
         * Sets whether the data has multiple points per column in a frame.
         *
         * \param[in] multipoint  Whether multiple points per column are
         *     possible.
         *
         * Can be called only before addModule() or notifyDataStart(),
         * otherwise asserts.
         *
         * \see isMultipoint()
         */
        void setMultipoint(bool multipoint);

        /*! \brief
         * Implements access to data frames.
         *
         * \param[in] index  Zero-based frame index to access.
         * \retval    Frame reference to frame \p index, or an invalid
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
         */
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const = 0;
        /*! \brief
         * Implements storage requests.
         *
         * \param[in] nframes  Request storing at least \c nframes previous
         *     frames (-1 = request storing all). Will be either -1 or >0.
         * \retval true if the request could be satisfied.
         *
         * Must not throw.  Failure to access a frame with the given index is
         * indicated through the return value.
         *
         * Derived classes should be prepared for any number of calls to this
         * method before notifyDataStart() is called (and during that call).
         *
         * This method is called internally by requestStorage().
         */
        virtual bool requestStorageInternal(int nframes) = 0;

        /*! \brief
         * Notifies attached modules of the start of data.
         *
         * Should be called once, after data properties have been set with
         * setColumnCount() and isMultipoint(), and before any of the
         * notification functions. The derived class should prepare for
         * requestStorage() calls from the attached modules.
         */
        void notifyDataStart();
        /*! \brief
         * Notifies attached modules of the start of a frame.
         *
         * Should be called once for each frame, before notifyPointsAdd() calls
         * for that frame.
         */
        void notifyFrameStart(const AnalysisDataFrameHeader &header) const;
        /*! \brief
         * Notifies attached modules of the addition of points to the
         * current frame.
         *
         * Can be called zero or more times for each frame.
         * The caller should ensure that any column occurs at most once in the
         * calls, unless the data is multipoint.
         * For efficiency reasons, calls to this method should be aggregated
         * whenever possible, i.e., it's better to handle multiple columns or
         * even the whole frame in a single call rather than calling the method
         * for each column separately.
         */
        void notifyPointsAdd(const AnalysisDataPointSetRef &points) const;
        /*! \brief
         * Notifies attached modules of the end of a frame.
         *
         * Should be called once for each call of notifyFrameStart(), after any
         * notifyPointsAdd() calls for the frame.
         */
        void notifyFrameFinish(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Notifies attached modules of the end of data.
         *
         * Should be called once, after all the other notification calls.
         */
        void notifyDataFinish() const;
        //! \endcond

    private:
        class Impl;

        PrivateImplPointer<Impl> _impl;
        int                     _ncol;
        bool                    _bMultiPoint;

        /*! \brief
         * Needed to provide access to notification methods.
         */
        friend class AnalysisDataStorage;
};

} // namespace gmx

#endif
