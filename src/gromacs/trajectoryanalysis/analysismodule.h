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
 * Declaration of TrajanaModule and integrally related classes.
 */
#ifndef GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H
#define GMX_TRAJECTORYANALYSIS_ANALYSISMODULE_H

#include <typedefs.h>

namespace gmx
{

class AnalysisData;
class AnalysisDataHandle;
class Options;
class SelectionCollection;
class TrajectoryAnalysisSettings;
class TopologyInformation;

/*! \brief
 * Base class for thread-local data storage during trajectory analysis.
 *
 * \see TrajectoryAnalysisModule::startFrames()
 */
class TrajectoryAnalysisModuleData
{
    public:
        virtual ~TrajectoryAnalysisModuleData() {};

        /*! \brief
         * Performs any finishing actions after all frames have been processed.
         *
         * \returns  0 on success, a non-zero error code on error.
         *
         * This function is called immediately before the destructor.
         */
        virtual int finish() = 0;

    protected:
        /*! \brief
         * Convenience function for finishing data handles.
         *
         * \param[in,out]  dhp  Pointer to data handle to finish
         *      (\p *dhp can be NULL).
         * \returns The \p oldrc if non-zero, return value of
         *      \c (*dhp)->finishData() otherwise
         *
         * Calls \p (*dhp)->finishData() if \p dhp is not NULL.  \p oldrc is
         * useful for chaining several calls on different handles so that all
         * the handles get finished even if there is an error in one of them,
         * and still return the error code of the first error.
         */
        static int finishHandle(AnalysisDataHandle **dhp, int oldrc = 0);
};

/*! \brief
 * Simple data storage class for a single data handle.
 *
 * Most simple tools should only require a single data handle to be
 * thread-local, so this class implements just that.
 */
class TrajectoryAnalysisModuleDataBasic : public TrajectoryAnalysisModuleData
{
    public:
        static int create(AnalysisData *data, /*AnalysisDataParallelOptions*/ void* opt,
                          TrajectoryAnalysisModuleData **pdatap);

        virtual int finish();

        AnalysisDataHandle *dataHandle() { return _dh; }

    private:
        TrajectoryAnalysisModuleDataBasic();

        AnalysisDataHandle     *_dh;
};


/*! \brief
 * Trajectory analysis method.
 *
 * For parallel analysis using threads, only a single object is constructed,
 * but the methods startFrames(), analyzeFrame() and finishFrames() are
 * called in each thread. Frame-local data should be initialized in
 * startFrames() and stored in a class derived from TrajanaModuleData that
 * is passed to the other methods.
 */
class TrajectoryAnalysisModule
{
    public:
        virtual ~TrajectoryAnalysisModule() {}

        /*! \brief
         * Sets the options object.
         *
         * Subclasses can access the options object through the options()
         * method, and as a short-hand, the selections through the selections()
         * method.
         *
         * There is no need for subclasses to ever use this function or even to
         * know of its existence.
         * Functions that create modules should call this function before
         * calling any other function in the interface.
         * The module does not take ownership of the options object: it is the
         * responsibility of the caller to ensure that the object remains valid
         * for the lifetime of the module.
         *
         * \see TrajanaOptions
         */

        /*! \brief
         * Initializes parameters understood by the module.
         *
         * In addition to initializing the accepted parameters, this function
         * should also set any required options using the options() object,
         * see TrajanaOptions for more details.
         * It can also customize the acceptable selections using selections(),
         * see Selections.
         *
         * If option values depend on the parameter values provided by the
         * user, see initParamsDone().
         */
        virtual Options *initOptions(TrajectoryAnalysisSettings *settings) = 0;
        /*! \brief
         * Called after all parameter values have been set.
         *
         * If the module needs to set options that affect topology loading or
         * selection initialization based on parameters values, this function
         * has to be overridden.
         */
        virtual int initOptionsDone(TrajectoryAnalysisSettings * /*settings*/)
        { return 0; }
        /*! \brief
         * Initializes the analysis.
         *
         * \returns  Zero on success, a non-zero error code on error.
         *
         * When this function is called, selections have been initialized based
         * on user input, and a topology has been loaded if provided by the
         * user. For dynamic selections, the selections have been evaluated to
         * the largest possible selection, i.e., the selections passed to
         * analyzeFrame() are always a subset of the selections provided here.
         */
        virtual int initAnalysis(const TopologyInformation &top) = 0;
        /*! \brief
         * Performs additional initialization after reading the first frame.
         *
         * \returns  Zero on success, a non-zero error code on error.
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
        virtual int initAfterFirstFrame(const t_trxframe &/*fr*/) { return 0; }

        /*! \brief
         * Starts the analysis of frames.
         *
         * \param[in]  opt
         * \param[in]  sel     Frame-local selection object.
         * \param[out] pdatap  Data structure for thread-local data.
         *
         * It is only necessary to override this function if the module
         * needs to support threaded parallelization. Otherwise, the
         * default implementation (which does nothing) can be used, and
         * all data stored as attributes in the subclass. With the default
         * implementation, the \p pdata pointer passed to analyzeFrame()
         * is always NULL.
         *
         * For threaded analysis, this function is called once for each
         * thread, and should initialize a class that contains any
         * required frame-local data in \p *pdatap.
         */
        virtual int startFrames(/*AnalysisDataParallelOptions*/ void* /*opt*/,
                                const SelectionCollection & /*sel*/,
                                TrajectoryAnalysisModuleData **pdatap)
        {
            *pdatap = NULL;
            return 0;
        }
        /*! \brief
         * Analyzes a single frame.
         *
         * \param[in]     frnr   Frame number, a zero-based index that
         *      uniquely identifies the frame.
         * \param[in]     fr     Current frame.
         * \param[in]     pbc    Periodic boundary conditions for \p fr.
         * \param[in,out] pdata  Data structure for frame-local data.
         * \return  0 on success, a non-zero error code or error.
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
         * data structure. It is also guaranteed that the \p sel object is
         * the same as was passed to startFrames() that created \p pdata.
         * Any access to data structures not stored in \p pdata
         * should be designed to be thread-safe.
         */
        virtual int analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
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
        virtual int finishFrames(TrajectoryAnalysisModuleData * /*pdata*/)
        {
            return 0;
        }

        /*! \brief
         * Postprocesses data after frames have been read.
         *
         * This function is called after all finishFrames() calls have been
         * called.
         */
        virtual int finishAnalysis(int nframes) = 0;
        /*! \brief
         * Writes output into files and/or standard output/error.
         *
         * All output from the module, excluding data written out for each
         * frame during analyzeFrame(), should be confined into this function.
         * This function is guaranteed to be called only after
         * finishAnalysis().
         */
        virtual int writeOutput() = 0;

        /*! \brief
         * Returns the number of datasets provided by the module.
         *
         * The default implementation returns 0.
         */
        //virtual int datasetCount() const { return 0; }
        /*! \brief
         * Returns a pointer to the data set \p index.
         *
         * \param[in] index  Data set to query for.
         * \returns   A pointer to the data set, or NULL if \p index is not
         *      valid.
         *
         * The return value is not const to allow callers to add modules to the
         * data sets. However, the AbstractData interface does not provide any
         * means to alter the data, so the module does not need to care about
         * external modifications.
         *
         * The default implementation always returns NULL.
         */
        //virtual AbstractAnalysisData *dataset(int index) const { return NULL; }

    protected:
        TrajectoryAnalysisModule() {}

        // Disallow copy and assign.
        TrajectoryAnalysisModule(const TrajectoryAnalysisModule &);
        void operator =(const TrajectoryAnalysisModule &);
};

} // namespace gmx

#endif
