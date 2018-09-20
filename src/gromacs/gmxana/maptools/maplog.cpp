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

#include "maplog.h"

#include <map>

#include "gromacs/compat/make_unique.h"
#include "gromacs/fileio/griddataview.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"

namespace gmx
{
namespace
{
enum class FillUndefinedLog
{
    Zero,
    NaN,
    MinFloat
};
const std::map<FillUndefinedLog, float> c_fillUndefinedLogToFloat = {
    {FillUndefinedLog::Zero, 0},
    {FillUndefinedLog::NaN, std::numeric_limits<float>::quiet_NaN()},
    {FillUndefinedLog::MinFloat, -GMX_FLOAT_MAX}
};
const char * const                      c_fillUndefinedLog[] = {  "NaN", "zero", "minfloat" };
}

class MapLog final : public ICommandLineOptionsModule
{
    public:
        void init(CommandLineModuleSettings *settings) override;
        void initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings) override;
        void optionsFinished() override;
        int run() override;
    private:
        static const char *const  helpText_[];
        std::string               fnInput_;
        std::string               fnOutput_;
        FillUndefinedLog          fillUndef_;
};


const char * const MapLog::helpText_[] = {
    "[THISMODULE] calculates the logarithm of a density map. When"
    " negative values are encountered, the output map reports NaN"
    " per default. The [TT]-fillundef[tt] option allows to chose"
    " zero and the smallest float number alternatively."
};

void MapLog::init(CommandLineModuleSettings * /*settings*/)
{}

void MapLog::initOptions(IOptionsContainer* options, ICommandLineOptionsModuleSettings *settings)
{
    settings->setHelpText(helpText_);

    options->addOption(StringOption("mi").store(&fnInput_).required().defaultValue("input.ccp4"));
    options->addOption(StringOption("mo").store(&fnOutput_).required().defaultValue("output.ccp4"));
    options->addOption(EnumOption<FillUndefinedLog>("fillundef")
                           .enumValue(c_fillUndefinedLog)
                           .store(&fillUndef_)
                           .description("Values to report where log of negative number is taken.")
                           .defaultValue(FillUndefinedLog::NaN));
}

void MapLog::optionsFinished()
{ }

int MapLog::run()
{
    auto        map            = MapConverter(fnInput_).map();
    const auto &undefFillValue = c_fillUndefinedLogToFloat.at(fillUndef_);
    std::transform(std::begin(map), std::end(map), std::begin(map),
                   [undefFillValue](decltype(map) ::value_type x){
                       return x > 0 ? log(x) : undefFillValue;
                   });

    MapConverter(map).to(fnOutput_);
    return EXIT_SUCCESS;
}

const char maplogInfo::shortDescription[] = "Log of density maps";
const char maplogInfo::name[]             = "maplog";
ICommandLineOptionsModulePointer maplogInfo::create()
{
    return compat::make_unique<MapLog>();
}

} // namespace gmx
