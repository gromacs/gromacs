/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018,2019, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "gromacs/mdlib/mdebin.h"

#include <cstdio>

#include <gtest/gtest.h>

#include "gromacs/mdlib/ebin.h"
#include "gromacs/mdlib/makeconstraints.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/unique_cptr.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{
namespace test
{
namespace
{

//! Wraps fclose to discard the return value to use it as a deleter with gmx::unique_cptr.
void fcloseWrapper(FILE *fp)
{
    fclose(fp);
}

class MdebinTest : public ::testing::Test
{
    public:
        TestFileManager fileManager_;

        // Objects needed to make t_mdebin
        t_inputrec   inputrec_;
        gmx_mtop_t   mtop_;
        t_mdebin    *mdebin_;
        unique_cptr<t_mdebin, done_mdebin> mdebinGuard_;

        // Objects needed for default energy output behavior.
        t_mdatoms                    mdatoms_;
        std::unique_ptr<Constraints> constraints_;
        matrix box_ = {{10, 0, 0}, {0, 10, 0}, {0, 0, 10}};
        gmx_enerdata_t               enerdata_;
        tensor totalVirial_, pressure_;

        // TODO This will be more elegant (and run faster) when we
        // refactor the output routines to write to a stream
        // interface, which can already be handled in-memory when
        // running tests.
        std::string logFilename_;
        FILE       *log_;
        unique_cptr<FILE, fcloseWrapper> logFileGuard_;

        TestReferenceData                refData_;
        TestReferenceChecker             checker_;

        MdebinTest() :
            logFilename_(fileManager_.getTemporaryFilePath(".log")),
            log_(std::fopen(logFilename_.c_str(), "w")), logFileGuard_(log_),
            checker_(refData_.rootChecker())
        {
            mdebin_ = init_mdebin(nullptr, &mtop_, &inputrec_, nullptr, false);
            mdebinGuard_.reset(mdebin_);
            constraints_ = makeConstraints(mtop_, inputrec_, false, log_, mdatoms_, nullptr,
                                           nullptr, nullptr, nullptr, false);
        }
        //! Helper function to generate synthetic data to output
        void setStepData(real testValue)
        {
            enerdata_.term[F_LJ]      = (testValue += 0.1);
            enerdata_.term[F_COUL_SR] = (testValue += 0.1);
            enerdata_.term[F_EPOT]    = (testValue += 0.1);
            enerdata_.term[F_EKIN]    = (testValue += 0.1);
            enerdata_.term[F_ETOT]    = (testValue += 0.1);
            enerdata_.term[F_TEMP]    = (testValue += 0.1);
            enerdata_.term[F_PRES]    = (testValue += 0.1);
            totalVirial_[XX][XX]      = (testValue += 0.1);
            totalVirial_[XX][YY]      = 0.0;
            totalVirial_[XX][ZZ]      = 0.0;
            totalVirial_[YY][XX]      = 0.0;
            totalVirial_[YY][YY]      = (testValue += 0.1);
            totalVirial_[YY][ZZ]      = 0.0;
            totalVirial_[ZZ][XX]      = 0.0;
            totalVirial_[ZZ][YY]      = 0.0;
            totalVirial_[ZZ][ZZ]      = (testValue += 0.1);
            pressure_[XX][XX]         = (testValue += 0.1);
            pressure_[XX][YY]         = 0.0;
            pressure_[XX][ZZ]         = 0.0;
            pressure_[YY][XX]         = 0.0;
            pressure_[YY][YY]         = (testValue += 0.1);
            pressure_[YY][ZZ]         = 0.0;
            pressure_[ZZ][XX]         = 0.0;
            pressure_[ZZ][YY]         = 0.0;
            pressure_[ZZ][ZZ]         = (testValue += 0.1);
        }

};

TEST_F(MdebinTest, HandlesEmptyAverages)
{
    ASSERT_NE(log_, nullptr);

    // Test printing values
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprNORMAL, mdebin_,
               nullptr, nullptr, nullptr, nullptr);
    // Test printing averages
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprAVER, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkInteger(mdebin_->ebin->nener, "Number of Energy Terms");
    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

TEST_F(MdebinTest, HandlesSingleStep)
{
    ASSERT_NE(log_, nullptr);

    // Add synthetic data for a single step
    real time      = 1.0;
    real testValue = 1.0;
    setStepData(testValue);
    upd_mdebin(mdebin_, false, true, time, 0.0, &enerdata_,
               nullptr, nullptr, nullptr, box_,
               nullptr, nullptr, totalVirial_, pressure_,
               nullptr, nullptr, constraints_.get());

    // Test printing values
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprNORMAL, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // Test printing averages
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprAVER, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkInteger(mdebin_->ebin->nener, "Number of Energy Terms");
    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

TEST_F(MdebinTest, HandlesTwoSteps)
{
    ASSERT_NE(log_, nullptr);

    // Add synthetic data for the first step
    real time      = 1.0;
    real testValue = 1.0;
    setStepData(testValue);
    upd_mdebin(mdebin_, false, true, time, 0.0, &enerdata_,
               nullptr, nullptr, nullptr, box_,
               nullptr, nullptr, totalVirial_, pressure_,
               nullptr, nullptr, constraints_.get());

    // Test printing values
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprNORMAL, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // Add synthetic data for the second step
    time += 0.005;
    setStepData(testValue += 1.0);
    upd_mdebin(mdebin_, false, true, time, 0.0, &enerdata_,
               nullptr, nullptr, nullptr, box_,
               nullptr, nullptr, totalVirial_, pressure_,
               nullptr, nullptr, constraints_.get());

    // Test printing values
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprNORMAL, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // Test printing averages
    print_ebin(nullptr, false, false, false, log_,
               0, 0, eprAVER, mdebin_,
               nullptr, nullptr, nullptr, nullptr);

    // We need to close the file before the contents are available.
    logFileGuard_.reset(nullptr);

    checker_.checkInteger(mdebin_->ebin->nener, "Number of Energy Terms");
    checker_.checkString(TextReader::readFileToString(logFilename_), "log");
}

}  // namespace
}  // namespace test
}  // namespace gmx
