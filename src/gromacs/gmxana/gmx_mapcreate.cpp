/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "gromacs/gmxana/gmx_ana.h"

#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/cmdlinehelpcontext.h"
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionfiletype.h"

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/fileio/griddataview.h"
#include "gromacs/fileio/xplor.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"

namespace gmx
{


int gmx_mapcreate(int argc, char *argv[])
{
    std::string desc =
        "[THISMODULE] creates a new density map. "
        "The map may be filled with zeros, ones, or gaussian noise with unit amplitude ";

    std::vector<const char *> bugs =
    {"This tool is under construction.", ""};
    Options                   options;

    std::string               fnOutput;
    options.addOption(StringOption("mo").store(&fnOutput).required());

    std::array<int, DIM>   extend {};
    options.addOption(IntegerOption("extend").store(extend.data()).vector().required());
    std::array<int, DIM>   start {};
    options.addOption(IntegerOption("start").store(start.data()).vector());
    float                  spacing {};
    options.addOption(FloatOption("spacing").store(&spacing).required());
    bool                   ones = false;
    options.addOption(BooleanOption("ones").store(&ones));
    bool                   noise = false;
    options.addOption(BooleanOption("noise").store(&noise));

    gmx::CommandLineParser(&options).parse(&argc, argv);
    options.finish();

    if (ones && noise)
    {
        GMX_THROW(InconsistentInputError(
                          "Cannot set map to noise and ones at the same time."));
    }
    const CommandLineHelpContext *context = GlobalCommandLineHelpContext::get();
    if (context != nullptr)
    {
        CommandLineHelpWriter(options).setHelpText(desc)
            .setKnownIssues(arrayRefFromArray(bugs.data(), bugs.size()))
            .writeHelp(*context);
    }

    XplorData xplorData {};
    std::copy(start.begin(), start.end(), std::begin(xplorData.crs_start));
    std::copy(extend.begin(), extend.end(), std::begin(xplorData.extend));

    // evaluate the end of the grid from start and extend; beware off-by-one
    std::transform(
            std::begin(start), std::end(start),
            std::begin(extend),
            std::begin(xplorData.crs_end),
            [](int strt, int extnd){return strt+extnd-1; }
            );

    // calculate the cell length from spacing and number of grid points
    std::transform(
            std::begin(xplorData.extend), std::end(xplorData.extend),
            std::begin(xplorData.cell_length.as_vec()),
            [spacing](int extnd){return extnd*spacing; }
            );

    auto map = gridDataFromXplor(xplorData);

    // initialize the map with data values
    if (ones)
    {
        std::fill(std::begin(map), std::end(map), 1);
    }
    else if (noise)
    {
        ThreeFry2x64<64>               rng;
        gmx::NormalDistribution<real>  normalDist;
        std::generate(
                std::begin(map), std::end(map),
                [&rng, &normalDist](){return normalDist(rng); }
                );
    }
    else
    {
        std::fill(std::begin(map), std::end(map), 0);
    }

    MapConverter(map).to(fnOutput);
    return 0;
}

}
