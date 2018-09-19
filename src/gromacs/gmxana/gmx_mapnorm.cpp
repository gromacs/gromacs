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
#include "gromacs/math/container/mask.h"
#include "gromacs/options/basicoptions.h"

namespace gmx
{

class MapNorm final : public ToolRunner
{

    public:
        void initOptions(Options* options) override;
        void optionsFinished() override;
        void analyze() override;
    private:
        std::string               fnInput;
        std::string               fnOutput;
};

void MapNorm::initOptions(Options *options)
{
    setHelpText( "[THISMODULE] normalizes density maps to zero mean and unit standard deviation");
    setBugDescriptions({});

    options->addOption(StringOption("mi").store(&fnInput));
    options->addOption(StringOption("mo").store(&fnOutput));
}

void MapNorm::optionsFinished()
{
}

void MapNorm::analyze()
{
    auto map   = MapConverter(fnInput).map();
    containeroperation::addScalar(&map, -containermeasure::meanIgnoreNaN(map));
    containeroperation::divide(&map, containermeasure::rmsIgnoreNaN(map));

    MapConverter(map).to(fnOutput);
}
}

int gmx_mapnorm(int argc, char *argv[])
{
    return gmx::MapNorm().run(argc, argv);
}
