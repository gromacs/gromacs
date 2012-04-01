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
/*! \internal \file
 * \brief
 * Declares internal implementation class for gmx::AbstractAnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ABSTRACTDATA_IMPL_H
#define GMX_ANALYSISDATA_ABSTRACTDATA_IMPL_H

#include <vector>

#include "../legacyheaders/types/simple.h"

#include "gromacs/utility/uniqueptr.h"

#include "abstractdata.h"
#include "dataframe.h"

namespace gmx
{

/*! \internal \brief
 * Private implementation class for AbstractAnalysisData.
 *
 * \ingroup module_analysisdata
 */
class AbstractAnalysisData::Impl
{
    public:
        //! Shorthand for list of modules added to the data.
        typedef std::vector<AnalysisDataModulePointer> ModuleList;

        Impl();
        ~Impl();

        /*! \brief
         * Present data already added to the data object to a module.
         *
         * \param[in] data   Data object to read data from.
         * \param[in] module Module to present the data to.
         * \throws    APIError if \p module is not compatible with the data
         *      object.
         * \throws    APIError if all data is not available through
         *      getDataFrame().
         * \throws    unspecified Any exception thrown by \p module in its data
         *      notification methods.
         *
         * Uses getDataFrame() in \p data to access all data in the object, and
         * calls the notification functions in \p module as if the module had
         * been registered to the data object when the data was added.
         */
        void presentData(AbstractAnalysisData *data,
                         AnalysisDataModuleInterface *module);

        //! List of modules added to the data.
        ModuleList              _modules;
        //! true if all modules support missing data.
        bool                    _bAllowMissing;
        //! Whether notifyDataStart() has been called.
        mutable bool            _bDataStart;
        //! Whether new data is being added.
        mutable bool            _bInData;
        //! Whether data for a frame is being added.
        mutable bool            _bInFrame;
        //! Index of the currently active frame.
        mutable int             _currIndex;
        /*! \brief
         * Total number of frames in the data.
         *
         * The counter is incremented in notifyFrameFinish().
         */
        int                     _nframes;
};

} // namespace gmx

#endif
