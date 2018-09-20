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

#include "mapconvert.h"

#include "gromacs/compat/make_unique.h"
#include "gromacs/fileio/griddataview.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"


namespace gmx
{

class MapConvert final : public ICommandLineOptionsModule
{
    public:
        void init(CommandLineModuleSettings *settings) override;
        void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings) override;
        void optionsFinished() override;
        int run() override;
    private:
        static const char * const helpText_[];
        std::string               fnInput_;
        std::string               fnOutput_;
};

const char * const MapConvert::helpText_[] = {
    "[THISMODULE] reads a density map, guesses its format from the"
    " file extension then outputs according to the to the output"
    " file format. Supported conversions are ccp4 and xplor to"
    " ccp4, xplor and dat."
};
void MapConvert::init(CommandLineModuleSettings * /*settings*/)
{}
void MapConvert::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings)
{
    settings->setHelpText(helpText_);

    options->addOption(StringOption("mi")
                           .defaultValue("input.ccp4")
                           .store(&fnInput_)
                           .required());

    options->addOption(StringOption("mo")
                           .defaultValue("output.ccp4")
                           .store(&fnOutput_)
                           .required());
}

void MapConvert::optionsFinished()
{
}

int MapConvert::run()
{
    MapConverter(fnInput_).to(fnOutput_);
    return EXIT_SUCCESS;
}


const char mapconvertInfo::shortDescription[] = "Convert between density maps";

ICommandLineOptionsModulePointer mapconvertInfo::create()
{
    return compat::make_unique<MapConvert>();
}

const char mapconvertInfo::name[] = "mapconvert";

} //namespace gmx
