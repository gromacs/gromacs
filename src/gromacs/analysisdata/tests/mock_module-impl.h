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
        explicit Impl(int flags);

        void startReferenceFrame(const AnalysisDataFrameHeader &header);
        void checkReferencePoints(const AnalysisDataPointSetRef &points);
        void finishReferenceFrame(const AnalysisDataFrameHeader &header);

        boost::scoped_ptr<TestReferenceChecker>  rootChecker_;
        boost::scoped_ptr<TestReferenceChecker>  frameChecker_;
        int                     flags_;
        int                     frameIndex_;
};

} // namespace test
} // namespace gmx

#endif
