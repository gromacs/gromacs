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
 * Declares mock implementation of gmx::AnalysisDataModuleInterface.
 *
 * Requires Google Mock.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_testutils
 */
#ifndef GMX_TESTUTILS_MOCK_DATAMODULE_H
#define GMX_TESTUTILS_MOCK_DATAMODULE_H

#include <boost/shared_ptr.hpp>
#include <gmock/gmock.h>

#include "gromacs/analysisdata/dataframe.h"
#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/utility/common.h"

namespace gmx
{
namespace test
{

class AnalysisDataTestInput;
class TestReferenceChecker;

class MockAnalysisDataModule : public AnalysisDataModuleInterface
{
    public:
        explicit MockAnalysisDataModule(int flags);
        virtual ~MockAnalysisDataModule();

        virtual int flags() const;

        MOCK_METHOD1(dataStarted, void(AbstractAnalysisData *data));
        MOCK_METHOD1(frameStarted, void(const AnalysisDataFrameHeader &header));
        MOCK_METHOD1(pointsAdded, void(const AnalysisDataPointSetRef &points));
        MOCK_METHOD1(frameFinished, void(const AnalysisDataFrameHeader &header));
        MOCK_METHOD0(dataFinished, void());

        void setupStaticCheck(const AnalysisDataTestInput &data,
                              AbstractAnalysisData *source);
        void setupStaticColumnCheck(const AnalysisDataTestInput &data,
                                    int firstcol, int n,
                                    AbstractAnalysisData *source);
        void setupStaticStorageCheck(const AnalysisDataTestInput &data,
                                     int storageCount,
                                     AbstractAnalysisData *source);
        void setupReferenceCheck(const TestReferenceChecker &checker,
                                 AbstractAnalysisData *source);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

//! Smart pointer to manage an MockAnalysisDataModule object.
typedef boost::shared_ptr<MockAnalysisDataModule>
MockAnalysisDataModulePointer;

} // namespace test
} // namespace gmx

#endif
