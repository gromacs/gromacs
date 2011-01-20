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
/*! \libinternal \file
 * \brief
 * Declares gmx::AbstractAnalysisData and gmx::AbstractAnalysisDataStored.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ABSTRACTDATA_H
#define GMX_ANALYSISDATA_ABSTRACTDATA_H

#include "../legacyheaders/types/simple.h"

namespace gmx
{

//! Additional error codes for functions in the analysis data module.
enum AnalysisDataErrorCode
{
    eedataDataNotAvailable = -1,
};

class AnalysisDataModuleInterface;

/*! \libinternal \brief
 * Abstract base class for all objects that provide data.
 *
 * The public interface includes functions for querying the data
 * (isMultipoint(), columnCount(), frameCount(), getDataWErr(), getData(),
 * getErrors(), requestStorage()) and functions for using modules for
 * processing the data (addModule(), addColumnModule(), applyModule()).
 *
 * There are also protected functions for use in derived classes:
 * setting the properties returned by isMultipoint() and columnCount()
 * through setMultipoint() and setColumnCount(), and notifying attached
 * modules of added data.
 *
 * Notice that even for non-const objects, the interface does not provide
 * any means of altering the data. It is only possible to add modules,
 * making it relatively safe to return a non-const pointer of this type
 * pointing to an internal data structure without worrying about possible
 * modifications of the data.
 *
 * It is up to subclasses to ensure that the protected functions are called
 * in a correct sequence (the functions will assert in most incorrect use
 * cases), and that the data provided through the public interface matches
 * that passed to the modules with the notify methods.
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
         * Derived classes can change the type by calling setMultipoint().
         * If this is not done, the function always returns false.
         */
        bool isMultipoint() const { return _bMultiPoint; }
        /*! \brief
         * Returns the number of columns in the data.
         *
         * \returns The number of columns in the data.
         *
         * Derived classes should set the number of columns with
         * setColumnCount().
         * If the number of columns is not set, returns 0.
         */
        int columnCount() const { return _ncol; }
        /*! \brief
         * Returns the total number of frames in the data.
         *
         * \returns The total number of frames in the data.
         *
         * This function returns the number of frames that the object has
         * produced. If requestStorage() has been successfully called,
         * getData() can be used to access some or all of these frames.
         */
        virtual int frameCount() const = 0;
        /*! \brief
         * Access stored data.
         *
         * \param[in]  index   Frame index to access
         *      (negative indices count backwards from the current frame).
         * \param[out] x
         * \param[out] dx
         * \param[out] y
         * \param[out] dy
         * \param[out] present Returns a pointer to an array that tells
         *      whether the corresponding column is present in that frame.
         *      If NULL, no missing information is returned.
         * \retval 0 on success.
         * \retval ::eedataDataNotAvailable if data for the requested frame is
         *      no longer available (this is regarded as a normal outcome,
         *      i.e., no error handlers are invoked).
         *
         * Derived classes can choose to return ::eedataDataNotAvailable if
         * requestStorage() has not been called at all, or if the frame is
         * too old (compared to the value given to requestStorage()).
         */
        virtual int getDataWErr(int index, real *x, real *dx,
                                const real **y, const real **dy,
                                const bool **present = 0) const = 0;
        /*! \brief
         * Convenience function for accessing stored data.
         *
         * \see getDataWErr()
         */
        int getData(int index, real *x, const real **y,
                    const bool **present = 0) const;
        /*! \brief
         * Convenience function for accessing errors for stored data.
         *
         * \see getDataWErr()
         */
        int getErrors(int index, real *dx, const real **dy) const;
        /*! \brief
         * Request storage of frames.
         *
         * \param[in] nframes  Request storing at least \c nframes previous
         *     frames (-1 = request storing all).
         * \retval 0 on success.
         * \retval eedataInvalidCall if data has already been added and
         *      cannot be stored.
         * \retval eedataNotSupported if the object does not support storage.
         *
         * If called multiple times, the largest request should be honored.
         *
         * \see getData()
         */
        virtual int requestStorage(int nframes = -1) = 0;

        /*! \brief
         * Adds a module to process the data.
         *
         * \param  module  Module to add.
         * \retval 0 on success.
         * \retval ::eeInvalidValue if \p module is not compatible with the
         *      data object.
         * \retval ::eedataDataNotAvailable if data has already been added to
         *      the data object and everything is not available through
         *      getData().
         *
         * If data has already been added to the module, the new module
         * immediately processes all existing data.  ::eedataDataNotAvailable
         * is returned if all data is not available through getData().
         *
         * If the call successful, the data object takes ownership of the
         * module, and automatically destructs it when the data object itself
         * is destroyed.
         */
        int addModule(AnalysisDataModuleInterface *module);
        /*! \brief
         * Adds a module that processes only a subset of the columns.
         *
         * \param[in] col     First column.
         * \param[in] span    Number of columns.
         * \param     module  Module to add.
         * \retval 0 on success.
         *
         * Can return any of the return codes for addModule().
         */
        int addColumnModule(int col, int span, AnalysisDataModuleInterface *module);
        /*! \brief
         * Applies a module to process data that is ready.
         *
         * \param  module  Module to apply.
         *
         * This function works as addModule(), except that it does not take
         * ownership of \p module. Also, it can only be called after the data
         * is ready, and only if getData() gives access to all of the data.
         * It is provided for additional flexibility in postprocessing
         * in-memory data.
         */
        int applyModule(AnalysisDataModuleInterface *module);

    protected:
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
         * Notifies attached modules of the start of data.
         *
         * Should be called once, after data properties have been set with
         * setColumnCount() and isMultipoint(), and before any of the
         * notification functions. The derived class should prepare for
         * requestStorage() calls from the attached modules.
         */
        int notifyDataStart();
        /*! \brief
         * Notifies attached modules of the start of a frame.
         *
         * Should be called once for each frame, before notifyPointsAdd() calls
         * for thet frame.
         */
        int notifyFrameStart(real x, real dx) const;
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
        int notifyPointsAdd(int firstcol, int n,
                            const real *y, const real *dy,
                            const bool *present) const;
        /*! \brief
         * Notifies attached modules of the end of a frame.
         *
         * Should be called once for each call of notifyFrameStart(), after any
         * notifyPointsAdd() calls for the frame.
         */
        int notifyFrameFinish() const;
        /*! \brief
         * Notifies attached modules of the end of data.
         *
         * Should be called once, after all the other notification calls.
         */
        int notifyDataFinish() const;

    private:
        class Impl;

        Impl                   *_impl;
        int                     _ncol;
        bool                    _bMultiPoint;

        // Disallow copy and assign.
        AbstractAnalysisData(const AbstractAnalysisData &);
        void operator =(const AbstractAnalysisData &);
};


/*! \libinternal \brief
 * Abstract class that implements storage of data.
 *
 * This class implements a standard way of storing data, to avoid implementing
 * storage in each derived class separately. All the pure virtual methods of
 * AbstractData are implemented, and protected methods are provided to add data
 * to the storage. These protected methods automatically call notifyDataStart(),
 * notifyFrameStart(), notifyPointsAdd() and notifyFrameFinish()
 * functions in AbstractAnalysisData, but the derived class should still
 * call notifyDataFinish().
 *
 * Additional protected functions could be implemented to allow optimization:
 * in the current interface, some data copying is unavoidable.
 * Some changes could make it possible to obtain a pointer to the
 * storage, allowing the calculated values to be stored there directly
 * instead of a temporary array.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AbstractAnalysisDataStored : public AbstractAnalysisData
{
    public:
        virtual ~AbstractAnalysisDataStored();

        virtual int frameCount() const;
        virtual int getDataWErr(int index, real *x, real *dx,
                                const real **y, const real **dy,
                                const bool **present = 0) const;
        virtual int requestStorage(int nframes = -1);

    protected:
        AbstractAnalysisDataStored();

        /*! \copydoc AbstractAnalysisData::setMultipoint()
         *
         * The overridden method also asserts if
         * storage has been requested and \p multipoint is \c true.
         */
        void setMultipoint(bool multipoint);

        //! Start storing data.
        int startDataStore();
        //! Starts storing a next frame.
        int startNextFrame(real x, real dx);
        //! Stores the whole frame in a single call after start_next_frame().
        int storeThisFrame(const real *y, const real *dy, const bool *present);
        //! Convenience function for storing a whole frame in a single call.
        int storeNextFrame(real x, real dx, const real *y, const real *dy,
                           const bool *present);

    private:
        class Impl;

        Impl                   *_impl;

        // Copy and assign disallowed by base class.
};

} // namespace gmx

#endif
