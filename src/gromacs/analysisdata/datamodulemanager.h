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
/*! \libinternal \file
 * \brief
 * Declares gmx::AnalysisDataModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAMODULEMANAGER_H
#define GMX_ANALYSISDATA_DATAMODULEMANAGER_H

#include <memory>

#include "gromacs/analysisdata/abstractdata.h"

namespace gmx
{

class AnalysisDataParallelOptions;
class AnalysisDataFrameHeader;
class AnalysisDataPointSetRef;
class IAnalysisDataModule;

/*! \libinternal \brief
 * Encapsulates handling of data modules attached to AbstractAnalysisData.
 *
 * See IAnalysisDataModule and \ref module_analysisdata for more
 * details on the notifications and the order in which they should be raised.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataModuleManager
{
public:
    /*! \brief
     * Identifies data properties to check with data modules.
     *
     * \see IAnalysisDataModule::Flag
     */
    enum DataProperty
    {
        eMultipleDataSets, //!< Data has multiple data sets.
        eMultipleColumns,  //!< Data has multiple columns.
        eMultipoint,       //!< Data is multipoint.
        eDataPropertyNR    //!< Number of properties; for internal use only.
    };

    AnalysisDataModuleManager();
    ~AnalysisDataModuleManager();

    /*! \brief
     * Allows the manager to check modules for compatibility with the data.
     *
     * \throws  APIError if any data module already added is not compatible
     *      with the new setting.
     *
     * Does two things: checks any modules already attached to the data and
     * throws if any of them is not compatible, and stores the property
     * to check modules attached in the future.
     *
     * Strong exception safety.
     */
    void dataPropertyAboutToChange(DataProperty property, bool bSet);

    /*! \brief
     * Whether there are modules that do not support parallel processing.
     *
     * Must not be called before notifyDataStart()/notifyParallelDataStart().
     * If notifyDataStart() has been called, returns true if there are any
     * modules (all modules are treated as serial).
     *
     * Does not throw.
     */
    bool hasSerialModules() const;

    /*! \brief
     * Adds a module to process the data.
     *
     * \param     data    Data object to add the module to.
     * \param     module  Module to add.
     * \throws    std::bad_alloc if out of memory.
     * \throws    APIError if
     *      - \p module is not compatible with the data object
     *      - data has already been added to the data object and everything
     *        is not available through getDataFrame().
     * \throws    unspecified Any exception thrown by \p module in its
     *      notification methods (if data has been added).
     *
     * \see AbstractAnalysisData::addModule()
     */
    void addModule(AbstractAnalysisData* data, const AnalysisDataModulePointer& module);
    /*! \brief
     * Applies a module to process data that is ready.
     *
     * \param     data    Data object to apply the module to.
     * \param     module  Module to apply.
     * \throws    APIError in same situations as addModule().
     * \throws    unspecified Any exception thrown by \p module in its
     *      notification methods.
     *
     * \see AbstractAnalysisData::applyModule()
     */
    void applyModule(AbstractAnalysisData* data, IAnalysisDataModule* module);

    /*! \brief
     * Notifies attached modules of the start of serial data.
     *
     * \param   data  Data object that is starting.
     * \throws  APIError if any attached data module is not compatible.
     * \throws  unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::dataStarted().
     *
     * Should be called once, after data properties have been set with
     * the methods in AbstractAnalysisData, and before any other
     * notification methods.
     * The caller should be prepared for requestStorage() calls to \p data
     * from the attached modules.
     *
     * \p data should typically be \c this when calling from a class
     * derived from AbstractAnalysisData.
     *
     * This method initializes all modules for serial processing by calling
     * IAnalysisDataModule::dataStarted().
     */
    void notifyDataStart(AbstractAnalysisData* data);
    /*! \brief
     * Notifies attached modules of the start of parallel data.
     *
     * \param     data    Data object that is starting.
     * \param[in] options Parallelization properties of the input data.
     * \throws  APIError if any attached data module is not compatible.
     * \throws  unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::parallelDataStarted().
     *
     * Can be called instead of notifyDataStart() if \p data supports
     * non-sequential creation of frames.  Works as notifyDataStart(),
     * but instead calls IAnalysisDataModule::parallelDataStarted()
     * and records whether the module supports the parallel mode.
     * Subsequent notification calls then notify the modules according to
     * the mode they accept.
     *
     * See notifyDataStart() for general constraints.
     */
    void notifyParallelDataStart(AbstractAnalysisData* data, const AnalysisDataParallelOptions& options);
    /*! \brief
     * Notifies attached serial modules of the start of a frame.
     *
     * \param[in] header  Header information for the frame that is starting.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::frameStarted().
     *
     * Should be called once for each frame, before notifyPointsAdd() calls
     * for that frame.
     */
    void notifyFrameStart(const AnalysisDataFrameHeader& header) const;
    /*! \brief
     * Notifies attached parallel modules of the start of a frame.
     *
     * \param[in] header  Header information for the frame that is starting.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::frameStarted().
     *
     * If notifyParallelDataStart() has been called, should be called once
     * for each frame, before notifyParallelPointsAdd() calls for that
     * frame.
     * It is allowed to call this method in any order for the frames, but
     * should be called exactly once for each frame.
     */
    void notifyParallelFrameStart(const AnalysisDataFrameHeader& header) const;
    /*! \brief
     * Notifies attached serial modules of the addition of points to the
     * current frame.
     *
     * \param[in] points  Set of points added (also provides access to
     *      frame-level data).
     * \throws    APIError if any attached data module is not compatible.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::pointsAdded().
     *
     * Can be called zero or more times for each frame.
     * The caller should ensure that any column occurs at most once in the
     * calls, unless the data is multipoint.
     * For efficiency reasons, calls to this method should be aggregated
     * whenever possible, i.e., it's better to handle multiple columns or
     * even the whole frame in a single call rather than calling the method
     * for each column separately.
     */
    void notifyPointsAdd(const AnalysisDataPointSetRef& points) const;
    /*! \brief
     * Notifies attached parallel modules of the addition of points to a frame.
     *
     * \param[in] points  Set of points added (also provides access to
     *      frame-level data).
     * \throws    APIError if any attached data module is not compatible.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::pointsAdded().
     *
     * See notifyPointsAdd() for information on the structure of the point
     * sets.
     */
    void notifyParallelPointsAdd(const AnalysisDataPointSetRef& points) const;
    /*! \brief
     * Notifies attached serial modules of the end of a frame.
     *
     * \param[in] header  Header information for the frame that is ending.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::frameFinished().
     *
     * Should be called once for each call of notifyFrameStart(), after any
     * notifyPointsAdd() calls for the frame.
     * \p header should be identical to that used in the corresponding
     * notifyFrameStart() call.
     *
     * This method also notifies parallel modules about serial end of frame.
     */
    void notifyFrameFinish(const AnalysisDataFrameHeader& header) const;
    /*! \brief
     * Notifies attached parallel modules of the end of a frame.
     *
     * \param[in] header  Header information for the frame that is ending.
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::frameFinished().
     *
     * Should be called once for each call of notifyParallelFrameStart(),
     * after any notifyParallelPointsAdd() calls for the frame.
     * \p header should be identical to that used in the corresponding
     * notifyParallelFrameStart() call.
     */
    void notifyParallelFrameFinish(const AnalysisDataFrameHeader& header) const;
    /*! \brief
     * Notifies attached modules of the end of data.
     *
     * \throws    unspecified Any exception thrown by attached data modules
     *      in IAnalysisDataModule::dataFinished().
     *
     * Should be called once, after all the other notification calls.
     */
    void notifyDataFinish() const;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
