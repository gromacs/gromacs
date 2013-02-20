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
 * Declares gmx::AnalysisData and gmx::AnalysisDataHandle.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ANALYSISDATA_H
#define GMX_ANALYSISDATA_ANALYSISDATA_H

#include "abstractdata.h"

namespace gmx
{

class AnalysisDataHandle;
class AnalysisDataParallelOptions;

/*! \brief
 * Parallelizable data container for raw data.
 *
 * This is the main class used to implement parallelizable data processing in
 * analysis tools.  It is used by first creating an object and setting its
 * properties using setColumnCount() and setMultipoint(), and attaching
 * necessary modules using addModule() etc.  Then one or more
 * AnalysisDataHandle objects can be created using startData().  Each data
 * handle can then be independently used to provide data frames (each frame
 * must be provided by a single handle, but different frames can be freely
 * mixed between the handles).  When all data has been provided, the handles
 * are destroyed using finishData() (or AnalysisDataHandle::finishData()).
 * The AnalysisData object takes care of internally sorting the frames and
 * passing them to the attached modules in the order in which the modules
 * expect them.
 *
 * \todo
 * Currently, multiple handles with multipoint data are not implemented.
 *
 * \todo
 * Parallel implementation is not complete.
 *
 * \if internal
 * Special note for MPI implementation: assuming that the initialization of
 * data objects is identical in all processes, associating the data objects
 * in different MPI processes should be possible without changes in the
 * interface.
 * Alternative, more robust implementation could get a unique ID as parameter
 * to the constructor or a separate function, but would require all tools to
 * provide it.  With the current registration mechanism in
 * TrajectoryAnalysisModule, this should be straightforward.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisData : public AbstractAnalysisData
{
    public:
        /*! \brief
         * Creates an empty analysis data object.
         *
         * \throws std::bad_alloc if out of memory.
         */
        AnalysisData();
        virtual ~AnalysisData();

        /*! \brief
         * Sets the number of columns in the data.
         *
         * \param[in] ncol  Number of columns in the data (must be > 0).
         *
         * Must be called before startData(), and can be called multiple times
         * before modules are added.
         * Must not be called after startData() has been called.
         *
         * Does not currently throw, but this may change for the case that
         * modules have already been added.
         */
        void setColumnCount(int ncol);
        /*! \brief
         * Sets whether the data contains multiple points per column per frame.
         *
         * \param[in] bMultipoint  Whether the data will allow multiple points
         *      per column within a single frame.
         *
         * If this method is not called, the data is not multipoint.
         *
         * Must not be called after modules have been added or startData() has
         * been called.
         *
         * Does not currently throw, but this may change for the case that
         * modules have already been added.
         *
         * \see isMultipoint()
         */
        void setMultipoint(bool bMultipoint);

        /*! \brief
         * Create a handle for adding data.
         *
         * \param[in]  opt     Options for setting how this handle will be
         *     used.
         * \returns The created handle.
         * \throws  std::bad_alloc if out of memory.
         * \throws  APIError if any attached data module is not compatible.
         * \throws  unspecified  Any exception thrown by attached data modules
         *      in AnalysisDataModuleInterface::dataStarted().
         *
         * The caller should retain the returned handle (or a copy of it), and
         * pass it to finishData() after successfully adding all data.
         * The caller should discard the returned handle if an error occurs;
         * memory allocated for the handle will be freed when the AnalysisData
         * object is destroyed.
         *
         * The \p opt options should be the same for all calls to this method,
         * and the number of calls should match the parallelization factor
         * defined in \p opt.
         */
        AnalysisDataHandle startData(const AnalysisDataParallelOptions &opt);
        /*! \brief
         * Destroy a handle after all data has been added.
         *
         * \param[in]  handle  Handle to destroy.
         * \throws  unspecified  Any exception thrown by attached data modules
         *      in AnalysisDataModuleInterface::dataFinished().
         *
         * \p handle must have been obtained from startData() of this object.
         * The order of the calls with respect to the corresponding startData()
         * calls is not important.
         *
         * The \p handle (and any copies) are invalid after the call.
         */
        void finishData(AnalysisDataHandle handle);

    private:
        virtual AnalysisDataFrameRef tryGetDataFrameInternal(int index) const;
        virtual bool requestStorageInternal(int nframes);

        class Impl;

        PrivateImplPointer<Impl> impl_;

        friend class AnalysisDataHandle;
};

namespace internal
{
class AnalysisDataHandleImpl;
}   // namespace internal

/*! \brief
 * Handle for inserting data into AnalysisData.
 *
 * This class provides an interface for adding data frames into an AnalysisData
 * object.  After a handle is obtained from AnalysisData::startData(), new
 * frames can be added using startFrame().  Then values for that frame are set
 * using provided methods (see below), and finishFrame() is called.  After all
 * frames have been added, finishData() (or AnalysisData::finishData()) must be
 * called.
 *
 * For simple (non-multipoint) data, within a frame values can be set using
 * setPoint() and setPoints().  Setting the same column multiple times
 * overrides previously set values.  When the frame is finished, attached
 * modules are notified.
 *
 * Multipoint data works otherwise similarly, but requires finishPointSet() to
 * be called for each set of points for which the modules need to be notified.
 * Each point set starts empty (after startFrame() or finishPointSet()), and
 * values can be set using setPoint()/setPoints().  finishPointSet() must also
 * be called for the last point set just before finishFrame().
 *
 * This class works like a pointer type: copying and assignment is lightweight,
 * and all copies work interchangeably, accessing the same internal handle.
 * However, normally you should only keep one copy of a handle, i.e., treat
 * this type as movable.
 * Several handles created from the same AnalysisData object can exist
 * concurrently, but must currently operate on separate frames.
 *
 * \inpublicapi
 * \ingroup module_analysisdata
 */
class AnalysisDataHandle
{
    public:
        /*! \brief
         * Constructs an invalid data handle.
         *
         * This constructor is provided for convenience in cases where it is
         * easiest to declare an AnalysisDataHandle without immediately
         * assigning a value to it.  Any attempt to call methods without first
         * assigning a value from AnalysisData::startData() to the handle
         * causes an assert.
         *
         * Does not throw.
         */
        AnalysisDataHandle();

        /*! \brief
         * Start data for a new frame.
         *
         * \param[in] index  Zero-based index for the frame to start.
         * \param[in] x      x value for the frame.
         * \param[in] dx     Error in x for the frame if applicable.
         *
         * \throws    unspecified  Any exception thrown by attached data
         *      modules in AnalysisDataModuleInterface::frameStarted().
         *
         * Each \p index value 0, 1, ..., N (where N is the total number of
         * frames) should be started exactly once by exactly one handle of an
         * AnalysisData object.  The frames may be started out of order, but
         * currently the implementation places some limitations on how far
         * the index can be in the future (as counted from the first frame that
         * is not finished).
         */
        void startFrame(int index, real x, real dx = 0.0);
        /*! \brief
         * Set a value for a single column for the current frame.
         *
         * \param[in] column  Zero-based column index.
         * \param[in] value   Value to set for the column.
         * \param[in] bPresent Present flag to set for the column.
         *
         * If called multiple times for a column (within one point set for
         * multipoint data), old values are overwritten.
         *
         * Does not throw.
         */
        void setPoint(int column, real value, bool bPresent = true);
        /*! \brief
         * Set a value and its error estimate for a single column for the
         * current frame.
         *
         * \param[in] column  Zero-based column index.
         * \param[in] value   Value to set for the column.
         * \param[in] error   Error estimate to set for the column.
         * \param[in] bPresent Present flag to set for the column.
         *
         * If called multiple times for a column (within one point set for
         * multipoint data), old values are overwritten.
         *
         * Does not throw.
         */
        void setPoint(int column, real value, real error, bool bPresent = true);
        /*! \brief
         * Set values for consecutive columns for the current frame.
         *
         * \param[in] firstColumn  Zero-based column index.
         * \param[in] count        Number of columns to set.
         * \param[in] values       Value array of \p column items.
         *
         * Equivalent to calling setPoint(firstColumn + i, values[i]) for
         * i from 0 to count.
         *
         * Does not throw.
         */
        void setPoints(int firstColumn, int count, const real *values);
        /*! \brief
         * Finish data for the current point set.
         *
         * \throws    APIError if any attached data module is not compatible.
         * \throws    unspecified  Any exception thrown by attached data
         *      modules in AnalysisDataModuleInterface::pointsAdded().
         *
         * Must be called after each point set for multipoint data, including
         * the last (i.e., no values must be set between the last call to this
         * method and AnalysisDataStorage::finishFrame()).
         * Must not be called for non-multipoint data.
         */
        void finishPointSet();
        /*! \brief
         * Finish data for the current frame.
         *
         * \throws    APIError if any attached data module is not compatible.
         * \throws    unspecified  Any exception thrown by attached data
         *      modules in frame notification methods.
         */
        void finishFrame();
        //! Calls AnalysisData::finishData() for this handle.
        void finishData();

    private:
        /*! \brief
         * Creates a new data handle associated with \p data.
         *
         * \param  data Data to associate the handle with.
         *
         * The constructor is private because data handles should only be
         * constructed through AnalysisData::startData().
         *
         * Does not throw.
         */
        explicit AnalysisDataHandle(internal::AnalysisDataHandleImpl *impl);

        /*! \brief
         * Pointer to the internal implementation class.
         *
         * The memory for this object is managed by the AnalysisData object,
         * and AnalysisDataHandle simply provides a public interface for
         * accessing the implementation.
         */
        internal::AnalysisDataHandleImpl *impl_;

        /*! \brief
         * Needed to access the non-public implementation.
         */
        friend class AnalysisData;
};

} // namespace gmx

#endif
