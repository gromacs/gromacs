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
 * Declares gmx::AnalysisDataModuleInterface.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_DATAMODULE_H
#define GMX_ANALYSISDATA_DATAMODULE_H

#include "types/simple.h"

namespace gmx
{

class AbstractAnalysisData;

/*! \brief
 * Interface for a module that gets notified whenever data is added.
 *
 * The interface provides one method (flags()) that describes features of
 * data objects the module supports. Only most common features are included
 * in the flags; custom checks can be implemented in the dataStarted() method
 * (see below).
 * All other methods in the interface are callbacks that are called by the
 * data object to which the module is attached to describe the data.
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
            efAllowMultipoint    = 0x01,
            //! The module does not make sense for non-multipoint data.
            efOnlyMultipoint     = 0x02,
            //! The module can process data with more than one column.
            efAllowMulticolumn   = 0x04,
            //! The module can process data with missing points.
            efAllowMissing       = 0x08,
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
         * Instead, AbstractData performs these checks based on the flags
         * provided.
         */
        virtual int flags() const = 0;

        /*! \brief
         * Called (once) when the data has been set up properly.
         *
         * The data to which the module is attached is passed as an argument
         * to provide access to properties of the data for initialization
         * and/or validation.
         * This is the only place where the module gets access to the data;
         * if properties of the data are required later, the module should
         * store them internally. It is guaranteed that the data properties
         * (column count, whether it's multipoint) do not change once this
         * method has been called.
         */
        virtual int dataStarted(AbstractAnalysisData *data) = 0;
        /*! \brief
         * Called at the start of each data frame.
         */
        virtual int frameStarted(real x, real dx) = 0;
        /*! \brief
         * Called one or more times during each data frame.
         *
         * For convenience, the \p x and \p dx values for the frame are
         * passed to each call of this function.
         */
        virtual int pointsAdded(real x, real dx, int firstcol, int n,
                                const real *y, const real *dy,
                                const bool *present) = 0;
        /*! \brief
         * Called when a data frame is finished.
         */
        virtual int frameFinished() = 0;
        /*! \brief
         * Called (once) when no more data is available.
         */
        virtual int dataFinished() = 0;

    protected:
        AnalysisDataModuleInterface() {}
};

} // namespace gmx

#endif
