/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::TimeUnitManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gmxpre.h"

#include "gromacs/options/timeunitmanager.h"

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include "gromacs/fileio/oenv.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class OptionSectionInfo;

namespace
{

/*! \brief
 * Enum values for a time unit.
 *
 * These must correspond to the TimeUnit enum in the header!
 */
const EnumerationArray<TimeUnit, const char*> c_timeUnitNames = {
    { "fs", "ps", "ns", "us", "ms", "s" }
};
/*! \brief
 * Scaling factors from each time unit to internal units (=picoseconds).
 *
 * These must correspond to the TimeUnit enum in the header!
 */
const EnumerationArray<TimeUnit, double> c_timeUnitScaleFactors = { { 1e-3, 1, 1e3, 1e6, 1e9, 1e12 } };

} // namespace

TimeUnitManager::TimeUnitManager() : timeUnit_(TimeUnit::Default) {}

TimeUnitManager::TimeUnitManager(TimeUnit unit) : timeUnit_(unit) {}

void TimeUnitManager::setTimeUnit(TimeUnit unit)
{
    timeUnit_ = unit;
}

const char* TimeUnitManager::timeUnitAsString() const
{
    return c_timeUnitNames[timeUnit_];
}

double TimeUnitManager::timeScaleFactor() const
{
    return c_timeUnitScaleFactors[timeUnit_];
}

double TimeUnitManager::inverseTimeScaleFactor() const
{
    return 1.0 / timeScaleFactor();
}

/********************************************************************
 * TimeUnitBehavior
 */

TimeUnitBehavior::TimeUnitBehavior() : timeUnit_(TimeUnit::Default), timeUnitStore_(nullptr) {}

void TimeUnitBehavior::setTimeUnit(TimeUnit unit)
{
    timeUnit_ = unit;
    if (timeUnitStore_ != nullptr)
    {
        *timeUnitStore_ = unit;
    }
}

void TimeUnitBehavior::setTimeUnitStore(TimeUnit* store)
{
    timeUnitStore_ = store;
    *store         = timeUnit();
}

void TimeUnitBehavior::setTimeUnitFromEnvironment()
{
    const char* const value = std::getenv("GMXTIMEUNIT");
    if (value != nullptr)
    {
        TimeUnit result = TimeUnit::Count;
        for (TimeUnit t : keysOf(c_timeUnitNames))
        {
            if (std::strcmp(value, c_timeUnitNames[t]) == 0)
            {
                result = t;
                break;
            }
        }
        if (result == TimeUnit::Count)
        {
            std::string message = formatString(
                    "Time unit provided with environment variable GMXTIMEUNIT=%s "
                    "is not recognized as a valid time unit.\n"
                    "Possible values are: %s",
                    value,
                    joinStrings(c_timeUnitNames, ", ").c_str());
            GMX_THROW(InvalidInputError(message));
        }
        setTimeUnit(result);
    }
}

void TimeUnitBehavior::addTimeUnitOption(IOptionsContainer* options, const char* name)
{
    options->addOption(
            EnumOption<TimeUnit>(name).enumValue(c_timeUnitNames).store(&timeUnit_).description("Unit for time values"));
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
template<class FloatingPointOptionInfo>
class TimeOptionScaler : public OptionsModifyingTypeVisitor<FloatingPointOptionInfo>
{
public:
    //! Initializes a scaler with the given factor.
    explicit TimeOptionScaler(double factor) : factor_(factor) {}

    void visitSection(OptionSectionInfo* section) override
    {
        OptionsModifyingIterator iterator(section);
        iterator.acceptSections(this);
        iterator.acceptOptions(this);
    }

    void visitOptionType(FloatingPointOptionInfo* option) override
    {
        if (option->isTime())
        {
            option->setScaleFactor(factor_);
        }
    }

private:
    double factor_;
};

} // namespace

void TimeUnitBehavior::optionsFinishing(Options* options)
{
    double factor = TimeUnitManager(timeUnit()).timeScaleFactor();
    TimeOptionScaler<DoubleOptionInfo>(factor).visitSection(&options->rootSection());
    TimeOptionScaler<FloatOptionInfo>(factor).visitSection(&options->rootSection());
    if (timeUnitStore_ != nullptr)
    {
        *timeUnitStore_ = static_cast<TimeUnit>(timeUnit_);
    }
}

} // namespace gmx
