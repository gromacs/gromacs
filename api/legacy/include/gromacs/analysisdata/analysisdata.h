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
/*! \defgroup module_analysisdata Parallelizable Handling of Output Data (analysisdata)
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for handling and processing output data from
 * analysis.
 *
 * <H3>Overview</H3>
 *
 * This module provides functionality to do common processing for tabular data
 * in analysis tools.  In addition to providing this common functionality, one
 * major driver for this module is to make it simple to write analysis tools
 * that process frames in parallel: the functionality in this module takes care
 * of necessary synchronization and communication such that output from the
 * frames is collected and output in the correct order.
 * See \ref page_analysisdata for an overview of the high-level functionality
 * and the terminology used.
 *
 * This module consists of two main parts.  The first is formed by the
 * gmx::AbstractAnalysisData class and classes that derive from it:
 * gmx::AnalysisData and gmx::AnalysisArrayData.  These classes are used to
 * process and store raw data as produced by the analysis tool.  They also
 * provide an interface to attach data modules that implement
 * gmx::IAnalysisDataModule.
 *
 * Modules that implement gmx::IAnalysisDataModule form the second part
 * of the module, and they provide functionality to do processing on the data.
 * These modules can also derive from gmx::AbstractAnalysisData, allowing other
 * modules to be attached to them to form a processing chain that best suits
 * the analysis tool.  Typically, such a processing chain ends in a plotting
 * module that writes the data into a file, but the final module can also
 * provide direct access to the processed data, allowing the analysis tool to
 * do custom postprocessing outside the module framework.
 *
 * <H3>Using Data Objects and Modules</H3>
 *
 * To use the functionality in this module, you typically declare one or more
 * AnalysisData objects and set its properties.  You then create some module
 * objects and set their properties (see the list of classes that implement
 * gmx::IAnalysisDataModule) and attach them to the data objects or to
 * one another using gmx::AbstractAnalysisData::addModule().  Then you add the
 * actual data values to the gmx::AnalysisData object, which automatically
 * passes it on to the modules.
 * After all data is added, you may optionally access some results directly
 * from the module objects or from the gmx::AnalysisData object itself.
 * However, in many cases it is sufficient to initially add a plotting module
 * to the processing chain, which will then automatically write the results
 * into a file.
 *
 * For simple processing needs with a small amount of data, an
 * gmx::AnalysisArrayData class is also provided, which keeps all the data in an
 * in-memory array and allows you to manipulate the data as you wish before you
 * pass the data to the attached modules.
 *
 * <H3>Data Modules</H3>
 *
 * Modules that derive from gmx::IAnalysisDataModule can operate in two
 * modes:
 *  - In _serial_ mode, the frames are presented to the module always in the
 *    order of increasing indices, even if they become ready in a different
 *    order in the attached data.
 *  - In _parallel_ mode, the frames are presented in the order that they
 *    become available in the input data, which may not be sequential.
 *    This mode allows the input data to optimize its behavior if it does not
 *    need to store and sort the frames.
 *
 * The figure below shows the sequence of callbacks that the module receives.
 * Arrows show a dependency between callbacks: the event at the start of the
 * arrow always occurs before the event at the end.  The events in the box are
 * repeated for each frame.  Dashed lines within this box show dependencies
 * between these frames:
 *  - In serial mode, all the events are called in a deterministic order, with
 *    each frame completely processed before the next starts.
 *  - In parallel mode, multiple frames can be in progress simultaneously, and
 *    the events for different frames can occur even concurrently on different
 *    threads.  However, frameFinishSerial() events will always occur in
 *    deterministic, sequential order for the frames.  Also, the number of
 *    concurrent frames is limited by the parallelization factor passed to
 *    parallelDataStarted(): only M frames after the last frame for which
 *    frameFinishSerial() has been called can be in progress
 *
 * \dot
 *     digraph datamodule_events {
 *         rankdir = LR
 *         node [ shape=box ]
 *
 *         start  [ label="dataStarted()",
 *                  URL="\ref gmx::IAnalysisDataModule::dataStarted()" ]
 *         pstart [ label="parallelDataStarted()",
 *                  URL="\ref gmx::IAnalysisDataModule::parallelDataStarted()" ]
 *         subgraph cluster_frame {
 *             label = "for each frame"
 *             framestart   [ label="frameStarted()",
 *                            URL="\ref gmx::IAnalysisDataModule::frameStarted()" ]
 *             pointsadd    [ label="pointsAdded()",
 *                            URL="\ref gmx::IAnalysisDataModule::pointsAdded()" ]
 *             framefinish  [ label="frameFinished()",
 *                            URL="\ref gmx::IAnalysisDataModule::frameFinished()" ]
 *             serialfinish [ label="frameFinishedSerial()",
 *                            URL="\ref gmx::IAnalysisDataModule::frameFinishedSerial()" ]
 *         }
 *         finish [ label="dataFinished()",
 *                  URL="\ref gmx::IAnalysisDataModule::dataFinished()" ]
 *
 *         start -> framestart
 *         pstart -> framestart
 *         framestart -> pointsadd
 *         pointsadd -> pointsadd [ label="0..*", dir=back ]
 *         pointsadd -> framefinish
 *         framefinish -> serialfinish
 *         serialfinish -> finish
 *
 *         framestart:se -> serialfinish:sw [ dir=back, style=dashed, weight=0,
 *                                            label="serial: frame n+1\nparallel: frame n+M" ]
 *         serialfinish -> serialfinish [ dir=back, style=dashed,
 *                                        label="frame n+1" ]
 *     }
 * \enddot
 *
 * If the input data supports parallel mode, it calls parallelDataStarted().
 * If the module returns `true` from this method, then it will process the
 * frames in the parallel mode.  If the module returns `false`, it will get the
 * frames in serial order.
 * If the input data does not support parallel mode, it calls dataStarted(),
 * and the module will always get the frames in order.
 *
 * The sequence of when the module methods are called with respect to when data
 * is added to the data object depends on the type of the module and the type
 * of the data.  However, generally the modules do not need to know the details
 * of how this happens, as long as they work with the above state diagram.
 *
 * For parallel processing, the gmx::AnalysisData object itself only provides
 * the infrastructure to support all of the above, including the reordering of
 * the frames for serial processing.  However, the caller is still responsible
 * of the actual thread synchronization, and must call
 * gmx::AnalysisData::finishFrameSerial() for each frame from a suitable
 * context where the serial processing for that frame can be done.  When using
 * the data objects as part of the trajectory analysis framework
 * (\ref page_trajectoryanalysisframework)
 * or energy analysis framework (\ref page_energyanalysisframework),
 * these calls are handled by the framework.
 *
 * \if libapi
 * <H3>Writing New Data and Module Objects</H3>
 *
 * New data modules can be implemented to perform custom operations that are
 * not supported by the modules provided in this module.  This is done by
 * creating a new class that implements gmx::IAnalysisDataModule.
 * If the new module computes values that can be used as input for other
 * modules, the new class should also derive from gmx::AbstractAnalysisData, and
 * preferably use gmx::AnalysisDataStorage internally to implement storage of
 * values.  See the documentation of the mentioned classes for more details on
 * how to implement custom modules.
 * When implementing a new module, it should be considered whether it can be of
 * more general use, and if so, it should be added to this module.
 *
 * It is also possible to implement new data source objects by deriving a class
 * from gmx::AbstractAnalysisData.  This should not normally be necessary, since
 * this module provides general data source objects for most typical uses.
 * If the classes in this module are not suitable for some specific use, it
 * should be considered whether a new generic class could be added (or an
 * existing extended) instead of implementing a local custom solution.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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

#include "gromacs/analysisdata/abstractdata.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class AnalysisDataHandle;
class AnalysisDataParallelOptions;

/*! \brief
 * Parallelizable data container for raw data.
 *
 * This is the main class used to implement parallelizable data processing in
 * analysis tools.  It is used by first creating an object and setting its
 * properties using setDataSetCount(), setColumnCount() and setMultipoint(),
 * and attaching necessary modules using addModule() etc.  Then one or more
 * AnalysisDataHandle objects can be created using startData().  Each data
 * handle can then be independently used to provide data frames (each frame
 * must be provided by a single handle, but different frames can be freely
 * mixed between the handles).  The finishFrameSerial() method must be called
 * in serial for each frame, after one of the handles has been used to provide
 * the data for that frame.  When all data has been provided, the handles
 * are destroyed using finishData() (or AnalysisDataHandle::finishData()).
 *
 * When used through the trajectory analysis framework, calls to startData(),
 * finishFrameSerial(), and finishData() are handled by the framework.
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
    ~AnalysisData() override;

    /*! \brief
     * Sets the number of data sets.
     *
     * \param[in] dataSetCount  Number of data sets (must be > 0).
     * \throws    std::bad_alloc if out of memory.
     * \throws    APIError if modules have been added that are not
     *      compatible with the new data set count.
     *
     * Must not be called after startData() has been called.
     * If not called, a single data set is assumed.
     * If called multiple times, the last call takes effect.
     */
    void setDataSetCount(size_t dataSetCount);
    /*! \brief
     * Sets the number of columns in a data set.
     *
     * \param[in] dataSet      Zero-based data set index.
     * \param[in] columnCount  Number of columns in the data (must be > 0).
     * \throws    APIError if modules have been added that are not
     *      compatible with the new column count.
     *
     * Must be called before startData() for each data set.
     * Must not be called after startData() has been called.
     * If called multiple times for a data set, the last call takes effect.
     */
    void setColumnCount(size_t dataSet, size_t columnCount);
    /*! \brief
     * Sets whether the data contains multiple points per column per frame.
     *
     * \param[in] bMultipoint  Whether the data will allow multiple points
     *      per column within a single frame.
     * \throws    APIError if modules have been added that are not
     *      compatible with the new setting.
     *
     * If this method is not called, the data is not multipoint.
     *
     * Must not be called after startData() has been called.
     *
     * \see isMultipoint()
     */
    void setMultipoint(bool bMultipoint);

    size_t frameCount() const override;

    /*! \brief
     * Creates a handle for adding data.
     *
     * \param[in]  opt     Options for setting how this handle will be
     *     used.
     * \returns The created handle.
     * \throws  std::bad_alloc if out of memory.
     * \throws  APIError if any attached data module is not compatible.
     * \throws  unspecified  Any exception thrown by attached data modules
     *      in IAnalysisDataModule::dataStarted().
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
    AnalysisDataHandle startData(const AnalysisDataParallelOptions& opt);
    /*! \brief
     * Performs in-order sequential processing for the next frame.
     *
     * \param[in]  frameIndex Index of the frame that has been finished.
     * \throws  unspecified  Any exception thrown by attached data modules
     *      in IAnalysisDataModule::frameFinishedSerial().
     *
     * This method should be called sequentially for each frame, after data
     * for that frame has been produced.  It is not necessary to call this
     * method if there is no parallelism, i.e., if only a single data
     * handle is created and the parallelization options provided at that
     * time do not indicate parallelism.
     */
    void finishFrameSerial(size_t frameIndex);
    /*! \brief
     * Destroys a handle after all data has been added.
     *
     * \param[in]  handle  Handle to destroy.
     * \throws  unspecified  Any exception thrown by attached data modules
     *      in IAnalysisDataModule::dataFinished().
     *
     * \p handle must have been obtained from startData() of this object.
     * The order of the calls with respect to the corresponding startData()
     * calls is not important.
     *
     * The \p handle (and any copies) are invalid after the call.
     */
    void finishData(AnalysisDataHandle handle);

private:
    AnalysisDataFrameRef tryGetDataFrameInternal(size_t index) const override;
    bool                 requestStorageInternal(size_t nframes) override;

    class Impl;

    std::unique_ptr<Impl> impl_;

    friend class AnalysisDataHandle;
};

namespace internal
{
class AnalysisDataHandleImpl;
} // namespace internal

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
 * selectDataSet(), setPoint() and setPoints().  Setting the same column in the
 * same data set multiple times overrides previously set values.
 * When the frame is finished, attached modules are notified.
 *
 * Multipoint data works otherwise similarly, but requires finishPointSet() to
 * be called for each set of points for which the modules need to be notified.
 * Each point set starts empty (after startFrame() or finishPointSet()), and
 * values can be set using setPoint()/setPoints().
 * A single point set can contain values only for a single data set, which must
 * be selected with selectDataSet() before setting any values.
 * finishPointSet() must also be called for the last point set just before
 * finishFrame().
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

    //! Returns whether this data handle is valid.
    bool isValid() const { return impl_ != nullptr; }

    /*! \brief
     * Start data for a new frame.
     *
     * \param[in] index  Zero-based index for the frame to start.
     * \param[in] x      x value for the frame.
     * \param[in] dx     Error in x for the frame if applicable.
     *
     * \throws    unspecified  Any exception thrown by attached data
     *      modules in IAnalysisDataModule::frameStarted().
     *
     * Each \p index value 0, 1, ..., N (where N is the total number of
     * frames) should be started exactly once by exactly one handle of an
     * AnalysisData object.  The frames may be started out of order, but
     * currently the implementation places some limitations on how far
     * the index can be in the future (as counted from the first frame that
     * is not finished).
     */
    void startFrame(size_t index, real x, real dx = 0.0);
    /*! \brief
     * Selects a data set for subsequent setPoint()/setPoints() calls.
     *
     * \param[in] index  Zero-based data set index.
     *
     * After startFrame(), the first data set is always selected.
     * The set value is remembered until the end of the current frame, also
     * across finishPointSet() calls.
     *
     * Does not throw.
     */
    void selectDataSet(size_t index);
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
    void setPoint(size_t column, real value, bool bPresent = true);
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
    void setPoint(size_t column, real value, real error, bool bPresent = true);
    /*! \brief
     * Set values for consecutive columns for the current frame.
     *
     * \param[in] firstColumn  Zero-based column index.
     * \param[in] count        Number of columns to set.
     * \param[in] values       Value array of \p column items.
     * \param[in] bPresent     Present flag to set for the column.
     *
     * Equivalent to calling setPoint(firstColumn + i, values[i], bPresent) for
     * i from 0 to count.
     *
     * Does not throw.
     */
    void setPoints(size_t firstColumn, size_t count, const real* values, bool bPresent = true);
    /*! \brief
     * Finish data for the current point set.
     *
     * \throws    APIError if any attached data module is not compatible.
     * \throws    unspecified  Any exception thrown by attached data
     *      modules in IAnalysisDataModule::pointsAdded().
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
     * \param  impl Data to associate the handle with.
     *
     * The constructor is private because data handles should only be
     * constructed through AnalysisData::startData().
     *
     * Does not throw.
     */
    explicit AnalysisDataHandle(internal::AnalysisDataHandleImpl* impl);

    /*! \brief
     * Pointer to the internal implementation class.
     *
     * The memory for this object is managed by the AnalysisData object,
     * and AnalysisDataHandle simply provides a public interface for
     * accessing the implementation.
     */
    internal::AnalysisDataHandleImpl* impl_;

    /*! \brief
     * Needed to access the non-public implementation.
     */
    friend class AnalysisData;
};

} // namespace gmx

#endif
