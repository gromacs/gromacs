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
 * Declares private implementation class for gmx::test::MockAnalysisModule.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_analysisdata
 */
#ifndef GMX_ANALYSISDATA_TESTS_MOCK_MODULE_IMPL_H
#define GMX_ANALYSISDATA_TESTS_MOCK_MODULE_IMPL_H

#include "mock_module.h"

#include <boost/scoped_ptr.hpp>

namespace gmx
{
namespace test
{

/*! \internal \brief
 * Private implementation class for gmx::test::MockAnalysisModule.
 *
 * \ingroup module_analysisdata
 */
class MockAnalysisModule::Impl
{
    public:
        //! Initializes a mock object with the given flags.
        explicit Impl(int flags);

        /*! \brief
         * Callback used to check frame start against reference data.
         *
         * Called to check parameters and order of calls to frameStarted().
         * In addition to reference data checks, this method checks statically
         * that the new frame matches \a frameIndex_.
         */
        void startReferenceFrame(const AnalysisDataFrameHeader &header);
        /*! \brief
         * Callback used to check frame points against reference data.
         *
         * Called to check parameters and order of calls to pointsAdded().
         */
        void checkReferencePoints(const AnalysisDataPointSetRef &points);
        /*! \brief
         * Callback used to check frame finish against reference data.
         *
         * Called to check parameters and order of calls to frameFinished().
         * \a frameIndex_ is incremented here.
         */
        void finishReferenceFrame(const AnalysisDataFrameHeader &header);

        /*! \brief
         * Reference data checker to use for checking frames.
         *
         * Must be non-NULL if startReferenceFrame() is called.
         */
        boost::scoped_ptr<TestReferenceChecker>  rootChecker_;
        /*! \brief
         * Reference data checker to use to check the current frame.
         *
         * Non-NULL between startReferenceFrame() and finishReferenceFrame()
         * calls.
         */
        boost::scoped_ptr<TestReferenceChecker>  frameChecker_;
        //! Flags that will be returned by the mock module.
        int                     flags_;
        //! Index of the current/next frame.
        int                     frameIndex_;
};

} // namespace test
} // namespace gmx

#endif
