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
 * Declares gmx::TrajectoryAnalysisModule and
 * gmx::TrajectoryAnalysisModuleData.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "gromacs/selection/selection.h" // For gmx::SelectionList
#include "gromacs/utility/classhelpers.h"

struct t_pbc;
struct t_trxframe;

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;
class AnalysisDataParallelOptions;
class Options;
class SelectionCollection;
class TopologyInformation;
class TrajectoryAnalysisModule;
class TrajectoryAnalysisSettings;

/*! \brief
 * Base class for thread-local data storage during trajectory analysis.
 *
 * Thread-local storage of data handles and selections is implemented in this
 * class; TrajectoryAnalysisModule instances can access the thread-local values
 * in their TrajectoryAnalysisModule::analyzeFrame() method using dataHandle()
 * and parallelSelection().
 *
 * \see TrajectoryAnalysisModule::startFrames()
 * \see TrajectoryAnalysisModule::analyzeFrame()
 * \see TrajectoryAnalysisModule::finishFrames()
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisModuleData
{
    public:
        virtual ~TrajectoryAnalysisModuleData();

        /*! \brief
         * Performs any finishing actions after all frames have been processed.
         *
         * \throws  unspecified Implementation may throw exceptions to indicate
         *      errors.
         *
         * This function is called immediately before the destructor, after
         * TrajectoryAnalysisModule::finishFrames().
         * Derived classes should implement any final operations that need to
         * be done after successful analysis.
         * All implementations should call finishDataHandles().
         */
        virtual void finish() = 0;

        /*! \brief
         * Returns a data handle for a given dataset.
         *
         * \param[in] data  Analysis data object.
         * \returns   Data handle for \p data stored in this thread-local data.
         *
         * \p data should have previously been registered with
         * TrajectoryAnalysisModule::registerAnalysisDataset().
         * If \p data has zero columns in all data sets, the returned data
         * handle is invalid.
         *
         * Does not throw.
         */
        AnalysisDataHandle dataHandle(const AnalysisData &data);
        /*! \brief
         * Returns a selection that corresponds to the given selection.
         *
         * \param[in] selection Global selection object.
         * \returns   Selection object corresponding to this thread-local data.
         *
         * \p selection is the selection object that was obtained from
         * SelectionOption.  The return value is the corresponding selection
         * in the selection collection with which this data object was
         * constructed with.
         *
         * Does not throw.
         */
        Selection parallelSelection(const Selection &selection);
        /*! \brief
         * Returns a set of selection that corresponds to the given selections.
         *
         * \throws std::bad_alloc if out of memory.
         *
         * Works as parallelSelection(), but for a list of selections at once.
         *
         * \see parallelSelection()
         */
        SelectionList parallelSelections(const SelectionList &selections);

    protected:
        /*! \brief
         * Initializes thread-local storage for data handles and selections.
         *
         * \param[in] module     Analysis module to use for data objects.
         * \param[in] opt        Data parallelization options.
         * \param[in] selections Thread-local selection collection.
         * \throws  std::bad_alloc if out of memory.
         * \throws  unspecified Can throw any exception thrown by
         *      AnalysisData::startData().
         *
         * Calls AnalysisData::startData() on all data objects registered with
         * TrajectoryAnalysisModule::registerAnalysisDataset() in \p module.
         * The handles are accessible through dataHandle().
         */
        TrajectoryAnalysisModuleData(TrajectoryAnalysisModule          *module,
                                     const AnalysisDataParallelOptions &opt,
                                     const SelectionCollection         &selections);

        /*! \brief
         * Calls finishData() on all data handles.
         *
         * \throws  unspecified Can throw any exception thrown by
         *      AnalysisDataHandle::finishData().
         *
         * This function should be called from the implementation of finish()
         * in all subclasses.
         */
        void finishDataHandles();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

//! Smart pointer to manage a TrajectoryAnalysisModuleData object.
typedef boost::shared_ptr<TrajectoryAnalysisModuleData>
    TrajectoryAnalysisModuleDataPointer;

/*! \brief
 * Base class for trajectory analysis modules.
 *
 * Trajectory analysis methods should derive from this class and override the
 * necessary virtual methods to implement initialization (initOptions(),
 * optionsFinished(), initAnalysis(), initAfterFirstFrame()), per-frame analysis
 * (analyzeFrame()), and final processing (finishFrames(), finishAnalysis(),
 * writeOutput()).
 *
 * For parallel analysis using threads, only a single object is constructed,
 * but the methods startFrames(), analyzeFrame() and finishFrames() are called
 * in each thread.  Frame-local data should be initialized in startFrames() and
 * stored in a class derived from TrajectoryAnalysisModuleData that is passed
 * to the other methods.  The default implementation of startFrames() can be
 * used if only data handles and selections need to be thread-local.
 *
 * To get the full benefit from this class,
 * \ref module_analysisdata "analysis data objects" and
 * \ref module_selection "selections" should be used in the implementation.
 * See the corresponding modules' documentation for details of how they work.
 *
 * Typical way of using AnalysisData in derived classes is to have the
 * AnalysisData object as a member variable and register it using
 * registerAnalysisDataset().  Analysis modules are initialized in
 * initAnalysis() and the processing chain is initialized.  If any of the
 * modules is required, e.g., for post-processing in finishAnalysis(), it can
 * be stored in a member variable.  To add data to the data object in
 * analyzeFrame(), a data handle is obtained using
 * TrajectoryAnalysisModuleData::dataHandle().
 *
 * Typical way of using selections in derived classes is to have the required
 * \ref Selection objects (or ::SelectionList objects) as member variables, and
 * add the required selection options in initOptions().  These member variables
 * can be accessed in initAnalysis() to get general information about the
 * selections.  In analyzeFrame(), these selection objects should not be used
 * directly, but instead TrajectoryAnalysisModuleData::parallelSelection()
 * should be used to obtain a selection object that works correctly also for
 * parallel analysis.
 *
 * Derived classes should use exceptions to indicate errors in the virtual
 * methods.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisModule
{
    public:
        virtual ~TrajectoryAnalysisModule();

        /*! \brief
         * Initializes options understood by the module.
         *
         * \param[in,out] options  Options object to add the options to.
         * \param[in,out] settings Settings to pass to and from the module.
         *
         * This method is called first after the constructor, and it should
         * add options understood by the module to \p options.  Output values
         * from options (including selections) should be stored in member
         * variables.
         *
         * In addition to initializing the options, this method can also
         * provide information about the module's requirements using the
         * \p settings object; see TrajectoryAnalysisSettings for more details.
         *
         * If settings depend on the option values provided by the user, see
         * optionsFinished().
         */
        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings) = 0;
        /*! \brief
         * Called after all option values have been set.
         *
         * \param[in,out] options  Options object in which options are stored.
         * \param[in,out] settings Settings to pass to and from the module.
         *
         * This method is called after option values have been assigned (but
         * interactive selection input has not yet been performed).
         *
         * If the module needs to change settings that affect topology loading
         * (can be done using the \p settings object) or selection
         * initialization (can be done using SelectionOptionInfo) based on
         * option values, this method has to be overridden.
         *
         * The default implementation does nothing.
         */
        virtual void optionsFinished(Options                    *options,
                                     TrajectoryAnalysisSettings *settings);
        /*! \brief
         * Initializes the analysis.
         *
         * \param[in]    settings Settings to pass to and from the module.
         * \param[in]    top      Topology information.
         *
         * When this function is called, selections have been initialized based
         * on user input, and a topology has been loaded if provided by the
         * user.  For dynamic selections, the selections have been evaluated to
         * the largest possible selection, i.e., the selections passed to
         * analyzeFrame() are always a subset of the selections provided here.
         */
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top) = 0;
        /*! \brief
         * Performs additional initialization after reading the first frame.
         *
         * When this function is called, selections are the same as in
         * initAnalysis(), i.e., they have not been evaluated for the first
         * frame.
         *
         * It is necessary to override this method only if the module needs to
         * do initialization for which it requires data from the first frame.
         *
         * The default implementation does nothing.
         */
        virtual void initAfterFirstFrame(const TrajectoryAnalysisSettings &settings,
                                         const t_trxframe                 &fr);

        /*! \brief
         * Starts the analysis of frames.
         *
         * \param[in]  opt
         * \param[in]  selections  Frame-local selection collection object.
         * \returns    Data structure for thread-local data.
         *
         * This function is necessary only for threaded parallelization.
         * It is called once for each thread and should initialize a class that
         * contains any required frame-local data in the returned value.
         * The default implementation creates a basic data structure that holds
         * thread-local data handles for all data objects registered with
         * registerAnalysisDataset(), as well as the thread-local selection
         * collection.  These can be accessed in analyzeFrame() using the
         * methods in TrajectoryAnalysisModuleData.
         * If other thread-local data is needed, this function should be
         * overridden and it should create an instance of a class derived from
         * TrajectoryAnalysisModuleData.
         *
         * \see TrajectoryAnalysisModuleData
         */
        virtual TrajectoryAnalysisModuleDataPointer startFrames(
            const AnalysisDataParallelOptions &opt,
            const SelectionCollection         &selections);
        /*! \brief
         * Analyzes a single frame.
         *
         * \param[in]     frnr   Frame number, a zero-based index that
         *      uniquely identifies the frame.
         * \param[in]     fr     Current frame.
         * \param[in]     pbc    Periodic boundary conditions for \p fr.
         * \param[in,out] pdata  Data structure for frame-local data.
         *
         * This method is called once for each frame to be analyzed, and should
         * analyze the positions provided in the selections.  Data handles and
         * selections should be obtained from the \p pdata structure.
         *
         * For threaded analysis, this method is called asynchronously in
         * different threads to analyze different frames.  The \p pdata
         * structure is one of the structures created with startFrames(),
         * but no assumptions should be made about which of these data
         * structures is used.  It is guaranteed that two instances of
         * analyzeFrame() are not running concurrently with the same \p pdata
         * data structure.
         * Any access to data structures not stored in \p pdata should be
         * designed to be thread-safe.
         */
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata) = 0;
        /*! \brief
         * Finishes the analysis of frames.
         *
         * \param[in]  pdata    Data structure for thread-local data.
         *
         * This method is called once for each call of startFrames(), with the
         * data structure returned by the corresponding startFrames().
         * The \p pdata object should be destroyed by the caller after this
         * function has been called.
         *
         * You only need to override this method if you need custom
         * operations to combine data from the frame-local data structures
         * to get the final result.  In such cases, the data should be
         * aggregated in this function and stored in a member attribute.
         *
         * The default implementation does nothing.
         *
         * \see startFrames()
         */
        virtual void finishFrames(TrajectoryAnalysisModuleData *pdata);

        /*! \brief
         * Postprocesses data after frames have been read.
         *
         * \param[in]  nframes  Total number of frames processed.
         *
         * This function is called after all finishFrames() calls have been
         * called.
         * \p nframes will equal the number of calls to analyzeFrame() that
         * have occurred.
         */
        virtual void finishAnalysis(int nframes) = 0;
        /*! \brief
         * Writes output into files and/or standard output/error.
         *
         * All output from the module, excluding data written out for each
         * frame during analyzeFrame(), should be confined into this function.
         * This function is guaranteed to be called only after
         * finishAnalysis().
         */
        virtual void writeOutput() = 0;

        /*! \brief
         * Returns the name of the analysis module.
         *
         * Does not throw.
         */
        const char *name() const;
        /*! \brief
         * Returns short description for the analysis module.
         *
         * Does not throw.
         */
        const char *description() const;
        /*! \brief
         * Returns the number of datasets provided by the module.
         *
         * Does not throw.
         */
        int datasetCount() const;
        /*! \brief
         * Returns a vector with the names of datasets provided by the module.
         *
         * Does not throw.
         */
        const std::vector<std::string> &datasetNames() const;
        /*! \brief
         * Returns a pointer to the data set \p index.
         *
         * \param[in] index  Data set to query for.
         * \returns   Reference to the requested data set.
         * \throws    APIError if \p index is not valid.
         *
         * \p index should be >= 0 and < datasetCount().
         *
         * The return value is not const to allow callers to add modules to the
         * data sets. However, the AbstractAnalysisData interface does not
         * provide any means to alter the data, so the module does not need to
         * care about external modifications.
         */
        AbstractAnalysisData &datasetFromIndex(int index) const;
        /*! \brief
         * Returns a pointer to the data set with name \p name
         *
         * \param[in] name  Data set to query for.
         * \returns   Reference to the requested data set.
         * \throws    APIError if \p name is not valid.
         *
         * \p name should be one of the names returned by datasetNames().
         *
         * The return value is not const to allow callers to add modules to the
         * data sets. However, the AbstractAnalysisData interface does not
         * provide any means to alter the data, so the module does not need to
         * care about external modifications.
         */
        AbstractAnalysisData &datasetFromName(const char *name) const;
        /*! \brief
         * Processes data in AnalysisData objects in serial for each frame.
         *
         * \param[in] frameIndex  Index of the frame that has been finished.
         *
         * This method is called by the framework in order for each frame,
         * after the analysis for that frame has been finished.  These calls
         * always execute in serial and in sequential frame order, even during
         * parallel analysis where multiple analyzeFrame() calls may be
         * executing concurrently.
         *
         * \see AnalysisData::finishFrameSerial()
         */
        void finishFrameSerial(int frameIndex);

    protected:
        /*! \brief
         * Initializes the dataset registration mechanism.
         *
         * \param[in] name         Name for the module.
         * \param[in] description  One-line description for the module.
         * \throws    std::bad_alloc if out of memory.
         */
        TrajectoryAnalysisModule(const char *name, const char *description);

        /*! \brief
         * Registers a dataset that exports data.
         *
         * \param     data  Data object to register.
         * \param[in] name  Name to register the dataset with.
         * \throws    std::bad_alloc if out of memory.
         *
         * Registers \p data as a dataset that provides output from the
         * analysis module.  Callers for the module can access the dataset
         * with datasetFromName() using \p name as an AbstractAnalysisData
         * object.  This allows them to add their own data modules to do extra
         * processing.
         *
         * \p name must be unique across all calls within the same
         * TrajectoryAnalysisModule instance.
         */
        void registerBasicDataset(AbstractAnalysisData *data, const char *name);
        /*! \brief
         * Registers a parallelized dataset that exports data.
         *
         * \param     data  AnalysisData object to register.
         * \param[in] name  Name to register the dataset with.
         * \throws    std::bad_alloc if out of memory.
         *
         * This method works as registerBasicDataset(), but additionally allows
         * data handles for \p data to be accessed using
         * TrajectoryAnalysisData.
         *
         * \see registerBasicDataset()
         */
        void registerAnalysisDataset(AnalysisData *data, const char *name);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        /*! \brief
         * Needed to access the registered analysis data sets.
         */
        friend class TrajectoryAnalysisModuleData;
};

//! Smart pointer to manage a TrajectoryAnalysisModule.
typedef boost::shared_ptr<TrajectoryAnalysisModule>
    TrajectoryAnalysisModulePointer;

} // namespace gmx

#endif
