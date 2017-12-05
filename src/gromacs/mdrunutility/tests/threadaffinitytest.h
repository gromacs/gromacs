/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
#ifndef GMX_MDRUNUTILITY_TESTS_THREADAFFINITYTEST_H
#define GMX_MDRUNUTILITY_TESTS_THREADAFFINITYTEST_H

#include <initializer_list>
#include <memory>

#include <gmock/gmock.h>

#include "gromacs/hardware/hw_info.h"
#include "gromacs/mdrunutility/threadaffinity.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/loggertest.h"

struct t_commrec;

namespace gmx
{

class HardwareTopology;

namespace test
{

class MockThreadAffinityAccess : public IThreadAffinityAccess
{
    public:
        MockThreadAffinityAccess();
        ~MockThreadAffinityAccess();

        void setSupported(bool supported) { supported_ = supported; }
        void setPhysicalNodeId(int nodeId) { physicalNodeId_ = nodeId; }

        virtual bool isThreadAffinitySupported() const { return supported_; }
        virtual int physicalNodeId() const { return physicalNodeId_; }
        MOCK_METHOD1(setCurrentThreadAffinityToCore, bool(int core));

    private:
        bool supported_;
        int  physicalNodeId_;
};

class ThreadAffinityTestHelper
{
    public:
        ThreadAffinityTestHelper();
        ~ThreadAffinityTestHelper();

        void setAffinitySupported(bool supported)
        {
            affinityAccess_.setSupported(supported);
        }
        void setAffinityOption(int affinityOption)
        {
            hwOpt_.thread_affinity = affinityOption;
        }
        void setOffsetAndStride(int offset, int stride)
        {
            hwOpt_.core_pinning_offset = offset;
            hwOpt_.core_pinning_stride = stride;
        }

        void setPhysicalNodeId(int nodeId)
        {
            affinityAccess_.setPhysicalNodeId(nodeId);
        }

        void setLogicalProcessorCount(int logicalProcessorCount);

        void setTotNumThreadsIsAuto(bool isAuto)
        {
            hwOpt_.totNumThreadsIsAuto = isAuto;
        }

        void expectAffinitySet(int core)
        {
            EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(core));
        }
        void expectAffinitySet(std::initializer_list<int> cores)
        {
            for (int core : cores)
            {
                expectAffinitySet(core);
            }
        }
        void expectAffinitySetThatFails(int core)
        {
            using ::testing::Return;
            EXPECT_CALL(affinityAccess_, setCurrentThreadAffinityToCore(core))
                .WillOnce(Return(false));
        }

        void expectWarningMatchingRegex(const char *re)
        {
            expectWarningMatchingRegexIf(re, true);
        }
        void expectWarningMatchingRegexIf(const char *re, bool condition)
        {
            expectLogMessageMatchingRegexIf(MDLogger::LogLevel::Warning, re, condition);
        }
        void expectInfoMatchingRegex(const char *re)
        {
            expectInfoMatchingRegexIf(re, true);
        }
        void expectInfoMatchingRegexIf(const char *re, bool condition)
        {
            expectLogMessageMatchingRegexIf(MDLogger::LogLevel::Info, re, condition);
        }
        void expectGenericFailureMessage()
        {
            expectGenericFailureMessageIf(true);
        }
        void expectGenericFailureMessageIf(bool condition)
        {
            expectWarningMatchingRegexIf("NOTE: Thread affinity was not set.", condition);
        }
        void expectPinningMessage(bool userSpecifiedStride, int stride)
        {
            std::string pattern = formatString("Pinning threads .* %s.* stride of %d",
                                               userSpecifiedStride ? "user" : "auto",
                                               stride);
            expectInfoMatchingRegex(pattern.c_str());
        }
        void expectLogMessageMatchingRegexIf(MDLogger::LogLevel level,
                                             const char *re, bool condition)
        {
            if (condition)
            {
                logHelper_.expectEntryMatchingRegex(level, re);
            }
        }

        void setAffinity(int nthread_local)
        {
            if (hwTop_ == nullptr)
            {
                setLogicalProcessorCount(1);
            }
            gmx_set_thread_affinity(logHelper_.logger(), cr_, &hwOpt_, *hwTop_,
                                    nthread_local, &affinityAccess_);
        }

    private:
        t_commrec                         *cr_;
        gmx_hw_opt_t                       hwOpt_;
        std::unique_ptr<HardwareTopology>  hwTop_;
        MockThreadAffinityAccess           affinityAccess_;
        LoggerTestHelper                   logHelper_;
};

} // namespace test
} // namespace gmx

#endif
