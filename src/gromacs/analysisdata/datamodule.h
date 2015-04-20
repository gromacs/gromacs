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
 * Declares gmx::AnalysisDataModuleInterface and related convenience classes.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAMODULE_H
#define GMX_ANALYSISDATA_DATAMODULE_H

namespace gmx
{

class AbstractAnalysisData;
class AnalysisDataFrameHeader;
class AnalysisDataParallelOptions;
class AnalysisDataPointSetRef;

/*! \brief
 * Interface for a module that gets notified whenever data is added.
 *
 * The interface provides one method (flags()) that describes features of
 * data objects the module supports.  Only most common features are included
 * in the flags; custom checks can be implemented in the dataStarted() and/or
 * parallelDataStarted() methods (see below).
 * All other methods in the interface are callbacks that are called by the
 * data object to which the module is attached to describe the data.
 * See \ref module_analysisdata for an overview of the notifications the
 * modules receive, and \ref page_analysisdata for overview of the terminology.
 *
 * Concrete modules typically do not directly derive from this interface, but
 * from either AnalysisDataModuleSerial or AnalysisDataModuleParallel.
 * Both these classes implement one of dataStarted()/parallelDataStarted() by
 * forwarding the calls to the other method of this pair.  This allows the
 * module to only implement the initialization once, without needing to worry
 * how to correctly handle both cases.
 *
 * Currently, if the module throws an exception, it requires the analysis tool
 * to terminate, since AbstractAnalysisData will be left in a state where it
 * is not possible to continue processing.  See a related todo item in
 * AbstractAnalysisData.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataModuleInterface
{
    public:
        /*! \brief
         * Possible flags for flags().
         */
        enum Flag
        {
            //! The module can process multipoint data.
            efAllowMultipoint           = 1<<0,
            //! The module does not make sense for non-multipoint data.
            efOnlyMultipoint            = 1<<1,
            //! The module can process data with more than one column.
            efAllowMulticolumn          = 1<<2,
            //! The module can process data with missing points.
            efAllowMissing              = 1<<3,
            //! The module can process data with multiple data sets.
            efAllowMultipleDataSets     = 1<<4
        };

        virtual ~AnalysisDataModuleInterface() {};

        /*! \brief
         * Returns properties supported by the module.
         *
         * The return value of this method should not change after the module
         * has been added to a data (this responsibility can, and in most cases
         * must, be delegated to the user of the module).
         *
         * The purpose of this method is to remove the need for common checks
         * for data compatibility in the classes that implement the interface.
         * Instead, AbstractAnalysisData performs these checks based on the
         * flags provided.
         *
         * Does not throw.
         */
        virtual int flags() const = 0;

        /*! \brief
         * Called (once) when the data has been set up properly.
         *
         * \param[in] data  Data object to which the module is added.
         * \throws    APIError if the provided data is not compatible.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         *
         * When the data is ready, either this method or parallelDataStarted()
         * is called, depending on the nature of the input data.  If this
         * method is called, the input data will always present the frames in
         * sequential order.
         *
         * The data to which the module is attached is passed as an argument
         * to provide access to properties of the data for initialization
         * and/or validation.  The module can also call
         * AbstractAnalysisData::requestStorage() if needed.
         *
         * This is the only place where the module gets access to the data;
         * if properties of the data are required later, the module should
         * store them internally.  It is guaranteed that the data properties
         * (column count, whether it's multipoint) do not change once this
         * method has been called.
         *
         * Notice that \p data will be a proxy object if the module is added as
         * a column module, not the data object for which
         * AbstractAnalysisData::addColumnModule() was called.
         */
        virtual void dataStarted(AbstractAnalysisData *data) = 0;
        /*! \brief
         * Called (once) for parallel data when the data has been set up.
         *
         * \param[in] data     Data object to which the module is added.
         * \param[in] options  Parallelization properties of the input data.
         * \returns   true if the module can process the input in
         *      non-sequential order.
         * \throws    APIError if the provided data is not compatible.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         *
         * This method is called instead of dataStarted() if the input data has
         * the capability to present data in non-sequential order.
         * If the method returns true, then the module accepts this and frame
         * notification methods may be called in that non-sequential order.
         * If the method returns false, then the frame notification methods are
         * called in sequential order, as if dataStarted() had been called.
         *
         * See dataStarted() for general information on initializing the data.
         * That applies to this method as well, with the exception that calling
         * AbstractAnalysisData::requestStorage() is currently not very well
         * supported (or rather, accessing the requested storage doesn't work).
         */
        virtual bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options) = 0;
        /*! \brief
         * Called at the start of each data frame.
         *
         * \param[in] frame  Header information for the frame that is starting.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void frameStarted(const AnalysisDataFrameHeader &frame) = 0;
        /*! \brief
         * Called one or more times during each data frame.
         *
         * \param[in] points  Set of points added (also provides access to
         *      frame-level data).
         * \throws    APIError if the provided data is not compatible.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         *
         * Can be called once or multiple times for a frame.  For all data
         * objects currently implemented in the library (and all objects that
         * will use AnalysisDataStorage for internal implementation), it is
         * called exactly once for each frame if the data is not multipoint,
         * but currently this restriction is not enforced.
         */
        virtual void pointsAdded(const AnalysisDataPointSetRef &points) = 0;
        /*! \brief
         * Called when a data frame is finished.
         *
         * \param[in] header  Header information for the frame that is ending.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void frameFinished(const AnalysisDataFrameHeader &header) = 0;
        /*! \brief
         * Called in sequential order for each frame after they are finished.
         *
         * \param[in] frameIndex   Index of the next finished frame.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         *
         * This method is called after frameFinished(), but with an additional
         * constraint that it is always called in serial and with an increasing
         * \p frameIndex.  Parallel data modules need this to serialize their
         * data for downsteam serial modules; AnalysisDataModuleSerial provides
         * an empty implementation, as there frameFinished() can be used for
         * the same purpose.
         */
        virtual void frameFinishedSerial(int frameIndex) = 0;
        /*! \brief
         * Called (once) when no more data is available.
         *
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void dataFinished() = 0;
};

/*! \brief
 * Convenience base class for serial analysis data modules.
 *
 * Implements the parallelDataStarted() method such that initialization is
 * always forwarded to dataStarted(), and the module always behaves as serial
 * (parallelDataStarted() returns false).
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataModuleSerial : public AnalysisDataModuleInterface
{
    public:
        virtual ~AnalysisDataModuleSerial() {}

        virtual int flags() const = 0;

        virtual void dataStarted(AbstractAnalysisData *data)              = 0;
        virtual void frameStarted(const AnalysisDataFrameHeader &frame)   = 0;
        virtual void pointsAdded(const AnalysisDataPointSetRef &points)   = 0;
        virtual void frameFinished(const AnalysisDataFrameHeader &header) = 0;
        virtual void dataFinished() = 0;

    private:
        virtual bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options);
        virtual void frameFinishedSerial(int /*frameIndex*/) {}
};

/*! \brief
 * Convenience base class for parallel analysis data modules.
 *
 * Implements the dataStarted() method such that initialization is always done
 * in parallelDataStarted().  dataStarted() calls are forwarded to
 * parallelDataStarted() using a dummy serial AnalysisDataParallelOptions.
 *
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
class AnalysisDataModuleParallel : public AnalysisDataModuleInterface
{
    public:
        virtual ~AnalysisDataModuleParallel() {}

        virtual int flags() const = 0;

        virtual bool parallelDataStarted(
            AbstractAnalysisData              *data,
            const AnalysisDataParallelOptions &options)                   = 0;
        virtual void frameStarted(const AnalysisDataFrameHeader &frame)   = 0;
        virtual void pointsAdded(const AnalysisDataPointSetRef &points)   = 0;
        virtual void frameFinished(const AnalysisDataFrameHeader &header) = 0;
        virtual void frameFinishedSerial(int index) = 0;
        virtual void dataFinished()                 = 0;

    private:
        virtual void dataStarted(AbstractAnalysisData *data);
};

} // namespace gmx

#endif
