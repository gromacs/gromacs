/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 * Tests for plotfiles.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/plotfile.h"

#include <string>

#include <gtest/gtest.h>

#include "gromacs/fileio/oenv.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/cmdlinetest.h"
#include "testutils/stringtest.h"
namespace
{

class PlotFileTest : public gmx::test::CommandLineTestBase
{
    public:
        PlotFileTest()
        {
            output_env_init_default(&oenv_);
        }
        ~PlotFileTest()
        {
            output_env_done(oenv_);
        }
        void addLegend()
        {
            gmx::PlotFile xvg("test.xvg", "title", "xaxis", "yaxis", oenv_);
            xvg.addLegend("Nice weather today!\n");
            xvg.startPlot();
            xvg.finishPlot();
        }
        void addComment()
        {
            gmx::PlotFile xvg("test.xvg", "title", "xaxis", "yaxis", oenv_);
            xvg.addComment("Nice weather today!\n");
            xvg.startPlot();
            xvg.finishPlot();
        }
        void addRows(int nlegend)
        {
            gmx::PlotFile xvg("test.xvg", "title", "xaxis", "yaxis", oenv_);
            for (int i = 0; i < nlegend; i++)
            {
                xvg.addLegend(gmx::formatString("L%d", i));
            }
            xvg.startPlot();
            std::vector<real> data;
            data.push_back(3);
            data.push_back(5);
            xvg.addRow(1, data);
            xvg.addRow(2, data);
            xvg.finishPlot();
        }
        void createDestroyPtr()
        {
            gmx::PlotFile *xvg = new gmx::PlotFile("test.xvg", "title", "xaxis", "yaxis", oenv_);
            delete xvg;
        }
    private:
        gmx_output_env_t *oenv_;
};

TEST_F(PlotFileTest, addLegend)
{
    addLegend();
}

TEST_F(PlotFileTest, addComment)
{
    addComment();
}

TEST_F(PlotFileTest, addRowsLegend)
{
    addRows(2);
}

TEST_F(PlotFileTest, addRowsNoLegend)
{
    addRows(0);
}

TEST_F(PlotFileTest, addRowsIncorrectLegend)
{
    addRows(1);
}

TEST_F(PlotFileTest, CreateDestroyPtr)
{
    createDestroyPtr();
}

} // namespace
