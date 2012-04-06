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
 * Declares gmx::TrajectoryAnalysisModule and
 * gmx::TrajectoryAnalysisModuleData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H

#include <string>
#include <vector>

#include "../legacyheaders/typedefs.h"

#include "../selection/selection.h" // For gmx::SelectionList
#include "../utility/common.h"
#include "../utility/uniqueptr.h"

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
 * using dataHandle() and parallelSelection().
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
         * This function is called immediately before the destructor.
         * All implementations should call finishDataHandles().
         */
        virtual void finish() = 0;

        /*! \brief
         * Returns a data handle for a given dataset.
         *
         * Allowed data sets are those that have been registered with
         * TrajectoryAnalysisModule::registerAnalysisDataset().
         */
        AnalysisDataHandle dataHandle(const AnalysisData &data);
        /*! \brief
         * Returns a selection that corresponds to the given selection.
         */
        Selection parallelSelection(const Selection &selection);
        /*! \brief
         * Returns a set of selection that corresponds to the given selections.
         */
        SelectionList parallelSelections(const SelectionList &selections);

    protected:
        /*! \brief
         * Initializes thread-local storage for data handles and selections.
         *
         * \param[in] module     Analysis module to use for data objects.
         * \param[in] opt        Data parallelization options.
         * \param[in] selections Thread-local selection collection.
         *
         * Calls AnalysisData::startData() on all data objects registered with
         * TrajectoryAnalysisModule::registerAnalysisDataset() in \p module.
         * The handles are accessible through dataHandle().
         */
        TrajectoryAnalysisModuleData(TrajectoryAnalysisModule *module,
                                     const AnalysisDataParallelOptions &opt,
                                     const SelectionCollection &selections);

        /*! \brief
         * Calls finishData() on all data handles.
         *
         * This function should be called from the implementation of finish()
         * in all subclasses.
         */
        void finishDataHandles();

    private:
        class Impl;

        PrivateImplPointer<Impl> _impl;
};

//! Smart pointer to manage a TrajectoryAnalysisModuleData object.
typedef gmx_unique_ptr<TrajectoryAnalysisModuleData>::type
        TrajectoryAnalysisModuleDataPointer;

/*! \brief
 * Base class for trajectory analysis methods.
 *
 * Trajectory analysis methods should derive from this class and override the
 * necessary virtual functions to implement initialization (initOptions(),
 * initOptionsDone(), initAnalysis(), initAfterFirstFrame()), per-frame analysis
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
         * In addition to initializing the options, this function can also
         * provide information about its requirements using the \p settings
         * object; see TrajectoryAnalysisSettings for more details.
         *
         * If settings depend on the option values provided by the user, see
         * initOptionsDone().
         */
        virtual Options &initOptions(TrajectoryAnalysisSettings *settings) = 0;
        /*! \brief
         * Called after all option values have been set.
         *
         * If the module needs to change settings that affect topology loading
         * or selection initialization based on option values, this function
         * has to be overridden.
         *
         * The default implementation does nothing.
         */
        virtual void initOptionsDone(TrajectoryAnalysisSettings *settings);
        /*! \brief
         * Initializes the analysis.
         *
         * When this function is called, selections have been initialized based
         * on user input, and a topology has been loaded if provided by the
         * user.  For dynamic selections, the selections have been evaluated to
         * the largest possible selection, i.e., the selections passed to
         * analyzeFrame() are always a subset of the selections provided here.
         */
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation &top) = 0;
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
        virtual void initAfterFirstFrame(const t_trxframe &fr);

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
                const SelectionCollection &selections);
        /*! \brief
         * Analyzes a single frame.
         *
         * \param[in]     frnr   Frame number, a zero-based index that
         *      uniquely identifies the frame.
         * \param[in]     fr     Current frame.
         * \param[in]     pbc    Periodic boundary conditions for \p fr.
         * \param[in,out] pdata  Data structure for frame-local data.
         *
         * This function is called once for each frame to be analyzed,
         * and should analyze the positions provided in \p sel.
         *
         * For threaded analysis, this function is called asynchronously in
         * different threads to analyze different frames. The \p pdata
         * structure is one of the structures created with startFrames(),
         * but no assumptions should be made about which of these data
         * structures is used. It is guaranteed that two instances of
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
         * This function is called once for each call of startFrames(),
         * with the data structure returned by the corresponding startFrames().
         * The \p pdata object should be destroyed by the caller after this
         * function has been called.
         *
         * You only need to override this method if you need custom
         * operations to combine data from the frame-local data structures
         * to get the final result. In such cases, the data should be
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
         * This function is called after all finishFrames() calls have been
         * called.
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
         * Returns the number of datasets provided by the module.
         */
        int datasetCount() const;
        /*! \brief
         * Returns a vector with the names of the datasets.
         */
        const std::vector<std::string> &datasetNames() const;
        /*! \brief
         * Returns a pointer to the data set \p index.
         *
         * \param[in] index  Data set to query for.
         * \returns   A pointer to the data set, or NULL if \p index is not
         *      valid.
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
         * \returns   A pointer to the data set, or NULL if \p name is not
         *      recognized.
         *
         * The return value is not const to allow callers to add modules to the
         * data sets. However, the AbstractAnalysisData interface does not
         * provide any means to alter the data, so the module does not need to
         * care about external modifications.
         */
        AbstractAnalysisData &datasetFromName(const char *name) const;

    protected:
        //! Initializes the dataset registration mechanism.
        TrajectoryAnalysisModule();

        /*! \brief
         * Registers a dataset that exports data.
         */
        void registerBasicDataset(AbstractAnalysisData *data, const char *name);
        /*! \brief
         * Registers a parallelized dataset that exports data.
         */
        void registerAnalysisDataset(AnalysisData *data, const char *name);

    private:
        class Impl;

        PrivateImplPointer<Impl> _impl;

        /*! \brief
         * Needed to access the registered analysis data sets.
         */
        friend class TrajectoryAnalysisModuleData;
};

//! Smart pointer to manage a TrajectoryAnalysisModule.
typedef gmx_unique_ptr<TrajectoryAnalysisModule>::type
        TrajectoryAnalysisModulePointer;

} // namespace gmx

#endif
