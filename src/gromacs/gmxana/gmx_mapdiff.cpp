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
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/toolrunner.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/operations/gridinterpolator.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/optionfiletype.h"

namespace gmx
{

class MapDiff final : public ToolRunner
{

    public:
        void initOptions(Options* options) override;
        void optionsFinished() override;
        void analyze() override;
    private:
        std::string fnInput_;
        std::string fnCompare_;
        std::string fnOutput_;
        real        refFactor_;
        real        compareFactor_;
        real        offset_;
};

void MapDiff::initOptions(Options* options)
{
    setHelpText("[THISMODULE] calculates the linear combination of density maps."
                " outputmap = reffactor * refmap + comparefactor * compare + offset."
                " Per default, the tool outputs the difference between two density"
                " maps by using [TT]-reffactor[tt] zero, [TT]-comparefactor[tt]"
                " -1 and [TT]-offset[tt] zero.");
    setBugDescriptions({});

    options->addOption(StringOption("refmap").store(&fnInput_).required().defaultValue("refmap.ccp4"));
    options->addOption(StringOption("compare").store(&fnCompare_).defaultValue("compare.ccp4"));
    options->addOption(StringOption("mo").store(&fnOutput_).required().defaultValue("output.ccp4"));
    options->addOption(FloatOption("reffactor").store(&refFactor_).defaultValue(1));
    options->addOption(FloatOption("comparefactor").store(&compareFactor_).defaultValue(-1));
    options->addOption(FloatOption("offset").store(&offset_).defaultValue(0));
}

void MapDiff::optionsFinished()
{ }

void MapDiff::analyze()
{
    auto referenceMap   = MapConverter(fnInput_).map();
    if (refFactor_ != 1)
    {
        containeroperation::multiply(&referenceMap, refFactor_);
    }
    if (!fnCompare_.empty())
    {
        auto compareMap = MapConverter(fnCompare_).map();
        if (!(referenceMap.getGrid() == compareMap.getGrid()))
        {
            const auto &interpolatedMap = interpolateLinearly(compareMap, referenceMap.getGrid());
            compareMap = interpolatedMap;
            fprintf(stderr, "\nWARNING: Map grids don't match. Performing linear interpolation between grids.\n");
        }
        if (compareFactor_ != 1)
        {
            containeroperation::multiply(&compareMap, compareFactor_);
        }
        containeroperation::add(&referenceMap, compareMap);
    }
    if (offset_ != 0)
    {
        containeroperation::addScalar(&referenceMap, offset_);
    }

    MapConverter(referenceMap).to(fnOutput_);
}

}

int gmx_mapdiff(int argc, char *argv[])
{
    return gmx::MapDiff().run(argc, argv);
}
