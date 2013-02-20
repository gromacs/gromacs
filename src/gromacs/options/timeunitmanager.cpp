/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \file
 * \brief
 * Implements gmx::TimeUnitManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_options
 */
#include "gromacs/options/timeunitmanager.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/options/optionsvisitor.h"
#include "gromacs/utility/gmxassert.h"

namespace
{

/*! \brief
 * Enum values for a time unit.
 *
 * These must correspond to the TimeUnit enum in the header!
 */
const char *const g_timeUnits[] = {
    "fs", "ps", "ns", "us", "ms",  "s", NULL
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
 * \ingroup module_options
 */
class TimeOptionScaler : public OptionsModifyingTypeVisitor<DoubleOptionInfo>
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

        void visitOptionType(DoubleOptionInfo *option)
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
    TimeOptionScaler(factor).visitSubSection(options);
}

} // namespace gmx
