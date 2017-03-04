/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Code with simple manager for plotfiles.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "plotfile.h"

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textwriter.h"
#include "gromacs/utility/unique_cptr.h"

namespace gmx
{

class PlotFileImpl
{
    public:
        PlotFileImpl(const char *fn, const char *title, const char *xaxis,
                     const char *yaxis, const struct gmx_output_env_t *oenv)
        {
            fp_.reset(xvgropen(fn, title, xaxis, yaxis, oenv));
            writer_ = new gmx::TextWriter(fp_.get());
        }
        PlotFileImpl(const PlotFileImpl &) = delete;
        ~PlotFileImpl()
        {
            delete writer_;
        }
        gmx::TextWriter *writer_;
    private:
        //! File pointer for storing xvg data
        gmx::unique_cptr<FILE, xvgrclose> fp_;
};

PlotFile::PlotFile(const char                    *fn,
                   const char                    *title,
                   const char                    *xaxis,
                   const char                    *yaxis,
                   const struct gmx_output_env_t *oenv)
{
    impl_ = new PlotFileImpl(fn, title, xaxis, yaxis, oenv);
}

PlotFile::~PlotFile()
{
    delete impl_;
}

void PlotFile::setTimeAxis(const std::vector<real> &t)
{
    t_ = t;
}

void PlotFile::setDataSets(const std::vector<std::vector<real> > &data)
{
    // This overrides the previously stored data
    data_ = data;
}

void PlotFile::plot()
{
    // First check that the length of the data series matches the time
    for (auto &d : data_)
    {
        GMX_RELEASE_ASSERT(t_.size() == d.size(), "Time series has different length than data");
    }
    // Now check the legend size is either zero or matches the number of data sets
    GMX_RELEASE_ASSERT(legend_.size() == 0 || legend_.size() == data_.size(),
                       "Legend size should be zero or match number of data sets");
    for (auto &c : comment_)
    {
        impl_->writer_->writeLine(formatString("# %s", c.c_str()));
    }
    int nleg = 0;
    for (auto &l : legend_)
    {
        impl_->writer_->writeLine(formatString("@s%d \"%s\"", nleg++, l.c_str()));
    }
    for (auto &d : data_)
    {
        impl_->writer_->writeLine("@type xy");
        for (size_t i = 0; i < t_.size(); i++)
        {
            std::string s = formatString("%10g  %16g", t_[i], d[i]);
            impl_->writer_->writeLine(s.c_str());
        }
        impl_->writer_->writeLine("&");
    }
}

}
