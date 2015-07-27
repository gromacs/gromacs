/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Tests utilities for test reference data.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/xvgtest.h"

#include <vector>

#include <gtest/gtest.h>

#include "gromacs/utility/filestream.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"

namespace xvgtest
{

/*! \brief
 * Helper class to convert static data to a stream
 */
class StringInputStream : public gmx::TextInputStream
{
    public:
        StringInputStream(std::vector<std::string> *input)
        {
            input_ = input;
            si_    = input_->begin();
        }
        virtual ~StringInputStream() {}

        /*! \brief
         * Reads a line (with newline included) from the stream.
         *
         * \param[out] line    String to receive the line.
         * \returns    `false` if nothing was read because the stream ended.
         *
         * On error or when `false` is returned, \p line will be empty.
         */
        virtual bool readLine(std::string *line)
        {
            if (si_ != input_->end())
            {
                line->assign(*si_);
                ++si_;
                return true;
            }
            return false;
        }
        /*! \brief
         * Closes the stream.
         *
         * It is not allowed to read from a stream after it has been closed.
         * See TextOutputStream::close() for rationale for a close() method
         * separate from the destructor.  For input, failures during close
         * should be rare, but it is clearer to keep the interface symmetric.
         */
        virtual void close() {};
    private:
        std::vector<std::string>          *input_;
        std::vector<std::string>::iterator si_;
};

TEST(XvgTests, ReadFile)
{
    std::vector<std::string> input;

    input.push_back("0     2905.86    -410.199\n");
    input.push_back("0.2     6656.67    -430.437\n");
    input.push_back("0.4     5262.44    -409.399\n");
    input.push_back("0.6     5994.69    -405.763\n");
    input.push_back("0.8     5941.37    -408.337\n");
    input.push_back("1     5869.87    -411.124\n");
    input.push_back("1.2     6216.03    -414.901\n");
    input.push_back("1.4      5433.4     -405.81\n");
    input.push_back("1.6     6342.64    -414.192\n");
    input.push_back("1.8     5881.55    -411.311\n");
    input.push_back("2     5974.44    -406.021\n");
    input.push_back("2.2     6162.28    -409.284\n");
    input.push_back("2.4     6285.24    -406.367\n");
    input.push_back("2.6     6355.36    -397.706\n");
    input.push_back("2.8     6150.98    -394.927\n");
    input.push_back("3     5775.47    -399.688\n");
    input.push_back("3.2     6281.27    -402.533\n");
    input.push_back("3.4     6011.39    -399.425\n");
    input.push_back("3.6     6355.96    -400.802\n");
    input.push_back("3.8     5936.91    -399.281\n");
    input.push_back("4     5966.04    -402.363\n");
    input.push_back("4.2     6198.57    -404.514\n");
    input.push_back("4.4     6177.54    -406.713\n");
    input.push_back("4.6     6089.12     -404.75\n");
    input.push_back("4.8     5971.48    -399.051\n");
    input.push_back("5     6207.51    -396.209\n");

    {
        gmx::test::TestReferenceData    data(gmx::test::erefdataUpdateAll);
        gmx::test::TestReferenceChecker checker(data.rootChecker());
        double                          tolerance = 1e-3;
        checker.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, tolerance));
        StringInputStream               sis(&input);
        gmx::test::checkXvgFile(&sis, &checker);
    }
}

} // namespace
