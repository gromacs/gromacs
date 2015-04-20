/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
 * Implements gmx::TimeUnitManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "timeunitmanager.h"

#include <cstdlib>

#include <algorithm>

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace
{

/*! \brief
 * Enum values for a time unit.
 *
 * These must correspond to the TimeUnit enum in the header!
 */
const char *const g_timeUnits[] = {
    "fs", "ps", "ns", "us", "ms",  "s"
};
/*! \brief
 * Scaling factors from each time unit to internal units (=picoseconds).
 *
 * These must correspond to the TimeUnit enum in the header!
 */
const double g_timeScaleFactors[] = {
    1e-3,    1,  1e3,  1e6,  1e9, 1e12
};

} // namespace

namespace gmx
{

TimeUnitManager::TimeUnitManager()
    : timeUnit_(eTimeUnit_ps)
{
}

TimeUnitManager::TimeUnitManager(TimeUnit unit)
{
    setTimeUnit(unit);
}

void TimeUnitManager::setTimeUnit(TimeUnit unit)
{
    GMX_RELEASE_ASSERT(unit >= 0 && unit <= eTimeUnit_s,
                       "Invalid time unit");
    timeUnit_ = unit;
}

const char *TimeUnitManager::timeUnitAsString() const
{
    GMX_RELEASE_ASSERT(timeUnit_ >= 0 && timeUnit_ <= eTimeUnit_s,
                       "Invalid time unit");
    return g_timeUnits[timeUnit_];
}

double TimeUnitManager::timeScaleFactor() const
{
    GMX_RELEASE_ASSERT(timeUnit_ >= 0
                       && (size_t)timeUnit_ < sizeof(g_timeScaleFactors)/sizeof(g_timeScaleFactors[0]),
                       "Time unit index has become out-of-range");
    return g_timeScaleFactors[timeUnit_];
}

double TimeUnitManager::inverseTimeScaleFactor() const
{
    return 1.0 / timeScaleFactor();
}

void TimeUnitManager::setTimeUnitFromEnvironment()
{
    const char *const value = std::getenv("GMXTIMEUNIT");
    if (value != NULL)
    {
        ConstArrayRef<const char *>                 timeUnits(g_timeUnits);
        ConstArrayRef<const char *>::const_iterator i =
            std::find(timeUnits.begin(), timeUnits.end(), std::string(value));
        if (i == timeUnits.end())
        {
            std::string message = formatString(
                        "Time unit provided with environment variable GMXTIMEUNIT=%s "
                        "is not recognized as a valid time unit.\n"
                        "Possible values are: %s",
                        value, joinStrings(timeUnits, ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        timeUnit_ = i - timeUnits.begin();
    }
}

void TimeUnitManager::addTimeUnitOption(Options *options, const char *name)
{
    options->addOption(StringOption(name).enumValue(g_timeUnits)
                           .defaultValue(g_timeUnits[timeUnit()])
                           .storeEnumIndex(&timeUnit_)
                           .description("Unit for time values"));
}

namespace
{

/*! \internal \brief
 * Option visitor that scales time options.
 *
 * \tparam FloatingPointOptionInfo  OptionInfo type for an option that provides
 *     isTime() and setScaleFactor() methods.
 *
 * \ingroup module_options
 */
template <class FloatingPointOptionInfo>
class TimeOptionScaler : public OptionsModifyingTypeVisitor<FloatingPointOptionInfo>
{
    public:
        //! Initializes a scaler with the given factor.
        explicit TimeOptionScaler(double factor) : factor_(factor) {}

        void visitSubSection(Options *section)
        {
            OptionsModifyingIterator iterator(section);
            iterator.acceptSubSections(this);
            iterator.acceptOptions(this);
        }

        void visitOptionType(FloatingPointOptionInfo *option)
        {
            if (option->isTime())
            {
                option->setScaleFactor(factor_);
            }
        }

    private:
        double                  factor_;
};

}   // namespace

void TimeUnitManager::scaleTimeOptions(Options *options) const
{
    double factor = timeScaleFactor();
    TimeOptionScaler<DoubleOptionInfo>(factor).visitSubSection(options);
    TimeOptionScaler<FloatOptionInfo>(factor).visitSubSection(options);
}

} // namespace gmx
