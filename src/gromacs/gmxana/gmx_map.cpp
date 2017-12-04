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
#include "gromacs/fileio/griddataio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xplor.h"
#include "gromacs/fileio/mrcmetadata.h"

namespace gmx
{



class MapConverter
{
    public:
        void read(std::string filename);
        void write(std::string filename);
        void setGrid(const GridDataReal3D &grid);
    private:
        GridDataReal3D grid_;
};

void MapConverter::read(std::string filename)
{
    const auto &extension = filename.substr(filename.find_last_of(".") + 1);
    if (extension == "ccp4" || extension == "map" || extension == "mrc")
    {
        grid_ = MrcFile().read(filename);
    }
    if (extension == "xplor")
    {
        auto outputfile = gmx_fio_fopen(filename.c_str(), "r");
        grid_ = gridDataFromXplor(XplorData(outputfile));
        gmx_fio_fclose(outputfile);
    }
}

void MapConverter::write(std::string filename)
{
    const auto &extension = filename.substr(filename.find_last_of(".") + 1);
    if (extension == "ccp4" || extension == "map" || extension == "mrc")
    {
        grid_ = MrcFile().read(filename);
    }
    if (extension == "xplor" || extension == "dat")
    {
        auto outputfile = gmx_fio_fopen(filename.c_str(), "w");
        if (extension == "xplor")
        {
            xplorFromGridData(grid_).write(outputfile);
        }
        if (extension == "dat")
        {
            const auto textDump = MrcMetaData().fromGrid(grid_.getGrid()).set_grid_stats(grid_).to_string();
            fprintf(outputfile, "%s", textDump.c_str());
        }
        gmx_fio_fclose(outputfile);
    }
    // if (extension=="df3")
    // {
    //     Df3File().write(filename, grid_).writePovray();
    // }
}

void MapConverter::setGrid(const GridDataReal3D &grid)
{
    grid_ = grid;
}

int gmx_map(int argc, char *argv[])
{
    std::string desc =
        "[THISMODULE] is a tool to read in and write out (electron) density maps.";

    std::vector<const char *> bugs =
    {"This tool is under construction.", ""};
    Options                   options;

    std::string               fnInput;
    bool                      hasInput;
    options.addOption(StringOption("mi").store(&fnInput).storeIsSet(&hasInput));
    std::string               fnOutput;
    options.addOption(StringOption("mo").store(&fnOutput));

    std::array<int, DIM>   extend {};
    bool                   hasExtend = false;
    options.addOption(IntegerOption("extend").store(extend.data()).vector().storeIsSet(&hasExtend));
    std::array<int, DIM>   start {};
    bool                   hasStart = false;
    options.addOption(IntegerOption("start").store(start.data()).vector().storeIsSet(&hasStart));
    std::array<int, DIM>   end {};
    bool                   hasEnd = false;
    options.addOption(IntegerOption("end").store(end.data()).vector().storeIsSet(&hasEnd));
    std::array<float, DIM> length {};
    bool                   hasLength = false;
    options.addOption(FloatOption("length").store(length.data()).vector().storeIsSet(&hasLength));
    float                  spacing {};
    bool                   hasSpacing = false;
    options.addOption(FloatOption("spacing").store(&spacing).storeIsSet(&hasSpacing));

    gmx::CommandLineParser(&options).parse(&argc, argv);
    options.finish();

    const CommandLineHelpContext *context = GlobalCommandLineHelpContext::get();
    if (context != nullptr)
    {
        CommandLineHelpWriter(options).setHelpText(desc)
            .setKnownIssues(arrayRefFromArray(bugs.data(), bugs.size()))
            .writeHelp(*context);
    }

    MapConverter mapConverter;

    if (hasInput)
    {
        mapConverter.read(fnInput);
    }
    bool hasMapParameter = hasExtend || hasStart || hasEnd || hasLength || hasSpacing;
    if (hasMapParameter && hasInput)
    {
        fprintf(stderr, "\nWARNING: ignoring map input parameter and reading all parameters from input map\n");
    }
    if (hasMapParameter && !hasInput)
    {
        if (hasSpacing && hasExtend && hasLength)
        {
            fprintf(stderr, "\nWARNING: ignoring length, because extend * spacing = length\n");
        }
        XplorData xplorData {};
        if (hasLength)
        {
            std::copy(length.begin(), length.end(), std::begin(xplorData.cell_length.as_vec()));
        }
        if (hasStart)
        {
            std::copy(start.begin(), start.end(), std::begin(xplorData.crs_start));
        }
        if (hasEnd)
        {
            std::copy(end.begin(), end.end(), std::begin(xplorData.crs_end));
        }
        if (hasExtend)
        {
            std::copy(extend.begin(), extend.end(), std::begin(xplorData.extend));
        }
        if (hasSpacing)
        {
            std::transform(std::begin(xplorData.extend), std::end(xplorData.extend), std::begin(xplorData.cell_length.as_vec()), [spacing](real l){return l*spacing; });
        }
        mapConverter.setGrid(gridDataFromXplor(xplorData));
    }
    mapConverter.write(fnOutput);

    return 0;
}

}
