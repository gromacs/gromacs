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
/*! \defgroup module_trajectoryanalysis Framework for Trajectory Analysis (trajectoryanalysis)
 * \ingroup group_analysismodules
 * \brief
 * Provides functionality for implementing trajectory analysis modules.
 *
 * This module implements a framework for implementing flexible trajectory
 * analysis routines.  It provides a base class for implementing analysis as
 * reusable modules that can be used from different contexts and can also
 * support per-frame parallelization.  It integrally uses functionality from the
 * following modules:
 *  - \ref module_options
 *  - \ref module_analysisdata
 *  - \ref module_selection
 *
 * The main interface of this module is the gmx::TrajectoryAnalysisModule class.
 * Analysis modules should derive from this class, and override the necessary
 * virtual methods to provide the actual initialization and analysis routines.
 * Classes gmx::TrajectoryAnalysisSettings and gmx::TopologyInformation (in
 * addition to classes declared in the above-mentioned modules) are used to pass
 * information to and from these methods.  gmx::TrajectoryAnalysisModuleData can
 * be used in advanced scenarios where the tool requires custom thread-local
 * data for parallel analysis.
 *
 * The sequence charts below provides an overview of how the trajectory
 * analysis modules typically interact with other components.
 * The first chart provides an overview of the call sequence of the most
 * important methods in gmx::TrajectoryAnalysisModule.
 * There is a runner, which is responsible for doing the work that is shared
 * between all trajectory analysis (such as reading the trajectory and
 * processing selections).  The runner then calls different methods in the
 * analysis module at appropriate points to perform the module-specific tasks.
 * The analysis module is responsible for creating and managing
 * gmx::AnalysisData objects, and the chart shows the most important
 * interactions with this module as well.  However, the runner takes
 * responsibility of calling gmx::AnalysisData::finishFrameSerial().
 * Interactions with options (for command-line option processing) and
 * selections is not shown for brevity: see \ref module_options for an overview
 * of how options work, and the second chart for a more detailed view of how
 * selections are accessed from an analysis module.
 * \msc
 *     runner,
 *     module [ URL="\ref gmx::TrajectoryAnalysisModule" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     runner box module [ label="caller owns runner and module objects" ];
 *     module => data [ label="create (in constructor)" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     runner => runner [ label="parse user input" ];
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     runner => runner [ label="initialize topology\nand selections" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     module => data [ label="initialize" ];
 *     runner => runner [ label="read frame 0" ];
 *     runner => module [ label="initAfterFirstFrame()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAfterFirstFrame()" ];
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => module [ label="analyzeFrame(0)",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     module => data [ label="add data",
 *                      URL="\ref gmx::AnalysisDataHandle" ];
 *     module => data [ label="finishFrame()",
 *                      URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner => data [ label="finishFrameSerial()",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => module [ label="analyzeFrame(1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     ...;
 *     --- [ label="loop over frames ends" ];
 *     runner => module [ label="finishAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::finishAnalysis()" ];
 *     module => data [ label="post-process data" ];
 *     runner => module [ label="writeOutput()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::writeOutput()" ];
 * \endmsc
 *
 * The second chart below shows the interaction with selections and options
 * with focus on selection options.  The gmx::TrajectoryAnalysisModule object
 * creates one or more gmx::Selection variables, and uses gmx::SelectionOption
 * to indicate them as the destination for selections.  This happens in
 * gmx::TrajectoryAnalysisModule::initOptions().  After the options have been
 * parsed (includes parsing any options present on the command-line or read
 * from files, but not those provided interactively),
 * gmx::TrajectoryAnalysisModule::optionsFinished() can adjust the selections
 * using gmx::SelectionOptionInfo.  This is done like this to allow the
 * analysis module to influence the interactive prompt of selections based on
 * what command-line options were given.  After optionsFinished() returns, the
 * interactive selection prompt is presented if necessary.  After this point,
 * all access to selections from the analysis module is through the
 * gmx::Selection variables: the runner is responsible for calling methods in
 * the selection library, and these methods update the content referenced by
 * the gmx::Selection variables.  See documentation of
 * gmx::TrajectoryAnalysisModule for details of what the selections contain at
 * each point.
 * \msc
 *     runner,
 *     options [ label="Options", URL="\ref module_options" ],
 *     selection [ label="selections", URL="\ref module_selection" ],
 *     module [ label="module", URL="\ref gmx::TrajectoryAnalysisModule" ];
 *
 *     runner box selection [ label="all these objects are owned by the framework" ];
 *     runner => module [ label="initOptions()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initOptions()" ];
 *     module => options [ label="addOption(SelectionOption)",
 *                         URL="\ref gmx::SelectionOption" ];
 *     module => options [ label="addOption() (other options)",
 *                         URL="\ref gmx::Options::addOption()" ];
 *     ...;
 *     runner << module;
 *     runner => options [ label="parse command-line parameters" ];
 *     options => selection [ label="parse selections" ];
 *     selection -> module [ label="initialize Selection variables",
 *                           URL="\ref gmx::Selection" ];
 *     runner << options;
 *     runner => module [ label="optionsFinished()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::optionsFinished()" ];
 *     module => selection [ label="adjust SelectionOptions",
 *                           URL="\ref gmx::SelectionOptionInfo" ];
 *     runner << module;
 *     runner => selection [ label="prompt missing selections" ];
 *     selection -> module [ label="initialize Selection variables",
 *                         URL="\ref gmx::Selection" ];
 *     runner => selection [ label="compile selections" ];
 *     selection -> module [ label="change content referenced\nby Selection variables" ];
 *     runner => module [ label="initAnalysis()",
 *                        URL="\ref gmx::TrajectoryAnalysisModule::initAnalysis()" ];
 *     ...;
 *     --- [ label="loop over frames starts" ];
 *     runner => runner [ label="read and initialize frame 0" ];
 *     runner => selection [ label="evaluate selections for frame 0" ];
 *     selection -> module [ label="change content referenced\nby Selection variables" ];
 *     ...;
 * \endmsc
 *
 * The final chart shows the flow within the frame loop in the case of parallel
 * (threaded) execution and the interaction with the \ref module_analysisdata
 * module in this case.  Although parallelization has not yet been implemented,
 * it has influenced the design and needs to be understood if one wants to
 * write modules that can take advantage of the parallelization once it gets
 * implemented.  The parallelization takes part over frames: analyzing a single
 * frame is one unit of work.  When the frame loop is started,
 * gmx::TrajectoryAnalysisModule::startFrames() is called for each thread, and
 * initializes an object that contains thread-local data needed during the
 * analysis.  This includes selection information, gmx::AnalysisDataHandle
 * objects, and possibly other module-specific variables.  Then, the runner
 * reads the frames in sequence and passes the work into the different threads,
 * together with the appropriate thread-local data object.
 * The gmx::TrajectoryAnalysisModule::analyzeFrame() calls are only allowed to modify
 * the thread-local data object; everything else is read-only.  For any output,
 * they pass the information to gmx::AnalysisData, which together with the
 * runner takes care of ordering the data from different frames such that it
 * gets processed in the right order.
 * When all frames are analyzed, gmx::TrajectoryAnalysisModule::finishFrames()
 * is called for each thread-local data object to destroy them and to
 * accumulate possible results from them into the main
 * gmx::TrajectoryAnalysisModule object.
 * Note that in the diagram, some part of the work attributed for the runner
 * (e.g., evaluating selections) will actually be carried out in the analysis
 * threads before gmx::TrajectoryAnalysisModule::analyzeFrame() gets called.
 * \msc
 *     runner,
 *     module [ label="module object" ],
 *     thread1 [ label="analysis\nthread 1" ],
 *     thread2 [ label="analysis\nthread 2" ],
 *     data [ label="analysis data", URL="\ref module_analysisdata" ];
 *
 *     module box thread2 [ label="single TrajectoryAnalysisModule object",
 *                          URL="\ref gmx::TrajectoryAnalysisModule" ];
 *     ...;
 *     --- [ label="loop over frames starts" ];
 *     runner => thread1 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread1 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread1 [ label="pdata1" ];
 *     runner => thread2 [ label="startFrames()",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::startFrames()" ];
 *     thread2 => data [ label="startData()",
 *                       URL="\ref gmx::AnalysisData::startData()" ];
 *     runner << thread2 [ label="pdata2" ];
 *     |||;
 *     runner => runner [ label="initialize frame 0" ];
 *     runner => thread1 [ label="analyzeFrame(0, pdata1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     runner => runner [ label="read and initialize frame 1" ];
 *     runner => thread2 [ label="analyzeFrame(1, pdata2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="add data",
 *                       URL="\ref gmx::AnalysisDataHandle" ];
 *     thread2 => data [ label="finishFrame(1)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner << thread2 [ label="analyzeFrame() (frame 1)" ];
 *     runner => runner [ label="read and initialize frame 2" ];
 *     runner => thread2 [ label="analyzeFrame(2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::analyzeFrame()" ];
 *     thread1 => data [ label="finishFrame(0)",
 *                       URL="\ref gmx::AnalysisDataHandle::finishFrame()" ];
 *     runner << thread1 [ label="analyzeFrame() (frame 0)" ];
 *     runner => data [ label="finishFrameSerial() (frame 0)",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     runner => data [ label="finishFrameSerial() (frame 1)",
 *                      URL="\ref gmx::AnalysisData::finishFrameSerial()" ];
 *     ...;
 *     runner => thread1 [ label="finishFrames(pdata1)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread1 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread1 -> module [ label="accumulate results" ];
 *     runner << thread1;
 *     runner => thread2 [ label="finishFrames(pdata2)",
 *                         URL="\ref gmx::TrajectoryAnalysisModule::finishFrames()" ];
 *     thread2 => data [ label="finishData()",
 *                       URL="\ref gmx::AnalysisData::finishData()" ];
 *     thread2 -> module [ label="accumulate results" ];
 *     runner << thread2;
 *     --- [ label="loop over frames ends" ];
 *     ...;
 * \endmsc
 *
 * In addition to the framework for defining analysis modules, this module also
 * provides gmx::TrajectoryAnalysisCommandLineRunner, which implements a
 * command-line program that runs a certain analysis module.
 *
 * Internally, the module also defines a set of trajectory analysis modules that
 * can currently be accessed only through gmx::registerTrajectoryAnalysisModules.
 *
 * For an example of how to implement an analysis tool using the framework, see
 * \ref template.cpp.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
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

#include <memory>
#include <string>
#include <vector>

#include "gromacs/selection/selection.h" // For gmx::SelectionList

struct t_pbc;
struct t_trxframe;

namespace gmx
{

class AbstractAnalysisData;
class AnalysisData;
class AnalysisDataHandle;
class AnalysisDataParallelOptions;
class IOptionsContainer;
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
    AnalysisDataHandle dataHandle(const AnalysisData& data);
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
    static Selection parallelSelection(const Selection& selection);
    /*! \brief
     * Returns a set of selection that corresponds to the given selections.
     *
     * \throws std::bad_alloc if out of memory.
     *
     * Works as parallelSelection(), but for a list of selections at once.
     *
     * \see parallelSelection()
     */
    static SelectionList parallelSelections(const SelectionList& selections);

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
    TrajectoryAnalysisModuleData(TrajectoryAnalysisModule*          module,
                                 const AnalysisDataParallelOptions& opt,
                                 const SelectionCollection&         selections);

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

    std::unique_ptr<Impl> impl_;
};

//! Smart pointer to manage a TrajectoryAnalysisModuleData object.
typedef std::unique_ptr<TrajectoryAnalysisModuleData> TrajectoryAnalysisModuleDataPointer;

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
    virtual void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) = 0;
    /*! \brief
     * Called after all option values have been set.
     *
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
    virtual void optionsFinished(TrajectoryAnalysisSettings* settings);
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
    virtual void initAnalysis(const TrajectoryAnalysisSettings& settings,
                              const TopologyInformation&        top) = 0;
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
    virtual void initAfterFirstFrame(const TrajectoryAnalysisSettings& settings, const t_trxframe& fr);

    /*! \brief
     * Starts the analysis of frames.
     *
     * \param[in]  opt         Parallel options
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
    virtual TrajectoryAnalysisModuleDataPointer startFrames(const AnalysisDataParallelOptions& opt,
                                                            const SelectionCollection& selections);
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
    virtual void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) = 0;
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
    virtual void finishFrames(TrajectoryAnalysisModuleData* pdata);

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
    const std::vector<std::string>& datasetNames() const;
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
    AbstractAnalysisData& datasetFromIndex(int index) const;
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
    AbstractAnalysisData& datasetFromName(const char* name) const;
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
     * \throws    std::bad_alloc if out of memory.
     */
    TrajectoryAnalysisModule();

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
    void registerBasicDataset(AbstractAnalysisData* data, const char* name);
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
    void registerAnalysisDataset(AnalysisData* data, const char* name);

private:
    class Impl;

    std::unique_ptr<Impl> impl_;

    /*! \brief
     * Needed to access the registered analysis data sets.
     */
    friend class TrajectoryAnalysisModuleData;
};

//! Smart pointer to manage a TrajectoryAnalysisModule.
typedef std::unique_ptr<TrajectoryAnalysisModule> TrajectoryAnalysisModulePointer;

} // namespace gmx

#endif
