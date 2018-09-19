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

#include "gromacs/fileio/griddataview.h"
#include "gromacs/fileio/xplor.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/toolrunner.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/seed.h"
#include "gromacs/random/threefry.h"


namespace gmx
{

namespace
{
enum class MapCreateType
{
    Nothing,
    Ones,
    Noise
};
const char * const  c_mapCreateTypes[] = { "nothing", "ones", "noise" };
}

class MapCreate final : public ToolRunner
{

    public:
        void initOptions(Options* options) override;
        void optionsFinished() override;
        void analyze() override;
    private:
        std::string            fnOutput_ {};
        std::array<int, DIM>   extend_ {};
        std::array<int, DIM>   start_ {};
        float                  spacing_ {};
        MapCreateType          mapCreateType_;

};

void MapCreate::initOptions(Options* options)
{
    setHelpText("[THISMODULE] creates a new density map. The map may be filled"
                " with zeros, ones, or gaussian noise with unit amplitude."
                " The map dimensions are inferred from its [TT]-extend[tt]"
                " in voxels and the [TT]-spacing[tt] in nm. If the map shall"
                " not start at the coordinate system origin, a shift vector."
                " in voxels may be supplied with the [TT]-start[tt] option.");

    setBugDescriptions({"This tool is under construction.", ""});

    options->addOption(StringOption("mo")
                           .store(&fnOutput_)
                           .defaultValue("new.ccp4")
                           .description("The output map file name.")
                           .required());

    options->addOption(IntegerOption("extend")
                           .store(extend_.data())
                           .vector()
                           .description("The extend of the map in voxels.")
                           .defaultValue(10)
                           .required());

    options->addOption(FloatOption("spacing")
                           .store(&spacing_)
                           .description("The spacing of the map in nm.")
                           .defaultValue(0.2)
                           .required());

    options->addOption(IntegerOption("start")
                           .store(start_.data())
                           .description("The map origin voxel.")
                           .defaultValue(0)
                           .vector());

    options->addOption(EnumOption<MapCreateType>("fillwith")
                           .enumValue(c_mapCreateTypes)
                           .store(&mapCreateType_)
                           .description("Map data generator.")
                           .defaultValue(MapCreateType::Nothing));
}

void MapCreate::optionsFinished()
{
}

void MapCreate::analyze()
{
    XplorData xplorData {};
    std::copy(start_.begin(), start_.end(), std::begin(xplorData.crs_start));
    std::copy(extend_.begin(), extend_.end(), std::begin(xplorData.extend));

    // evaluate the end of the grid from start and extend_; beware off-by-one
    std::transform(
            std::begin(start_), std::end(start_),
            std::begin(extend_),
            std::begin(xplorData.crs_end),
            [](int strt, int extnd){return strt+extnd-1; }
            );

    auto spacing = spacing_;     // copy member  variable content for caputre in lambda
    // calculate the cell length from spacing and number of grid points
    std::transform(
            std::begin(xplorData.extend), std::end(xplorData.extend),
            std::begin(xplorData.cell_length.as_vec()),
            [spacing](int extnd){return extnd*spacing; }
            );

    auto map = gridDataFromXplor(xplorData);

    // initialize the map with data values
    switch (mapCreateType_)
    {
        case MapCreateType::Noise:
        {
            ThreeFry2x64<64>               rng(std::rand());
            gmx::NormalDistribution<real>  normalDist;
            std::generate(
                    std::begin(map), std::end(map),
                    [&rng, &normalDist](){return normalDist(rng); }
                    );
        }
        break;
        case MapCreateType::Ones:
            std::fill(std::begin(map), std::end(map), 1.);
            break;
        case MapCreateType::Nothing:
            std::fill(std::begin(map), std::end(map), 0.);
            break;
    }

    MapConverter(map).to(fnOutput_);
}

}

int gmx_mapcreate(int argc, char *argv[])
{
    return gmx::MapCreate().run(argc, argv);
}
