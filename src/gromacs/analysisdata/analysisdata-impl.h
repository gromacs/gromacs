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
 * Declares private implementation class for gmx::AnalysisData.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_ANALYSISDATA_IMPL_H
#define GMX_ANALYSISDATA_ANALYSISDATA_IMPL_H

#include <vector>

#include "analysisdata.h"

#include "abstractdata-impl.h"

namespace gmx
{

/*! \internal \brief
 * Private implementation class for AnalysisData.
 *
 * \ingroup module_analysisdata
 */
class AnalysisData::Impl
{
    public:
        //! Shorthand for a list of data handles.
        typedef std::vector<AnalysisDataHandle *> HandleList;
        //! Shorthand for a list of frames.
        typedef std::vector<AnalysisDataFrame *>  FrameList;

        //! Creates an implementation class associated with the given object.
        explicit Impl(AnalysisData *data);
        ~Impl();

        /*! \brief
         * Handles a new frame.
         *
         * If all earlier frames are ready, the data is directly passed to
         * AbstractStoredData.  Otherwise, it is put into the correct location
         * in \a _pending.  Calls processPendingFrame() after processing \p fr.
         */
        int addPendingFrame(AnalysisDataFrame *fr);
        /*! \brief
         * Passes pending frames to base class if all earlier frames are ready.
         */
        int processPendingFrames();
        //! Increments \a _pstart.
        void incrementPStart();

        //! The data object that uses this implementation class.
        AnalysisData           &_data;
        //! List of handles for this data object.
        HandleList              _handles;
        /*! \brief
         * Index into \a _pending that points to the location where the current
         * frame would go.
         */
        size_t                  _pstart;
        /*! \brief
         * Circular buffer for frames that are ready but waiting for earlier
         * frames.
         */
        FrameList               _pending;
};

/*! \internal \brief
 * Private implementation class for AnalysisDataHandle.
 *
 * \ingroup module_analysisdata
 */
class AnalysisDataHandle::Impl
{
    public:
        //! Creates a handle associated with the given data object.
        explicit Impl(AnalysisData *data);
        ~Impl();

        //! The data object that this handle belongs to.
        AnalysisData            &_data;
        //! Frame object where the current frame is being accumulated.
        AnalysisDataFrame       *_frame;
};

} // namespace gmx

#endif
