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
 * Declares gmx::AnalysisDataModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAMODULE_H
#define GMX_ANALYSISDATA_DATAMODULE_H

#include "../legacyheaders/types/simple.h"

namespace gmx
{

class AbstractAnalysisData;
class AnalysisDataFrameHeader;
class AnalysisDataPointSetRef;

/*! \brief
 * Interface for a module that gets notified whenever data is added.
 *
 * The interface provides one method (flags()) that describes features of
 * data objects the module supports.  Only most common features are included
 * in the flags; custom checks can be implemented in the dataStarted() method
 * (see below).
 * All other methods in the interface are callbacks that are called by the
 * data object to which the module is attached to describe the data.
 *
 * The frames are presented to the module always in the order of increasing
 * indices, even if they become ready in a different order in the attached
 * data.
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
        enum {
            //! The module can process multipoint data.
            efAllowMultipoint  = 0x01,
            //! The module does not make sense for non-multipoint data.
            efOnlyMultipoint   = 0x02,
            //! The module can process data with more than one column.
            efAllowMulticolumn = 0x04,
            //! The module can process data with missing points.
            efAllowMissing     = 0x08,
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
         * Called at the start of each data frame.
         *
         * \param[in] frame  Header information for the frame that is starting.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void frameStarted(const AnalysisDataFrameHeader &frame)   = 0;
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
        virtual void pointsAdded(const AnalysisDataPointSetRef &points)   = 0;
        /*! \brief
         * Called when a data frame is finished.
         *
         * \param[in] header  Header information for the frame that is ending.
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void frameFinished(const AnalysisDataFrameHeader &header) = 0;
        /*! \brief
         * Called (once) when no more data is available.
         *
         * \throws    unspecified  Can throw any exception required by the
         *      implementing class to report errors.
         */
        virtual void dataFinished() = 0;

    protected:
        AnalysisDataModuleInterface() {}
};

} // namespace gmx

#endif
