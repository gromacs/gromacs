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
 * Implements classes in abstractoption.h and abstractoptionstorage.h.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#include "gromacs/options/abstractoption.h"

#include "gromacs/options/abstractoptionstorage.h"
#include "gromacs/options/optionflags.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"

#include "basicoptionstorage.h"

namespace gmx
{

/********************************************************************
 * AbstractOptionStorage
 */

AbstractOptionStorage::AbstractOptionStorage(const AbstractOption &settings,
                                             OptionFlags staticFlags)
    : flags_(settings.flags_ | staticFlags),
      minValueCount_(settings.minValueCount_),
      maxValueCount_(settings.maxValueCount_),
      bInSet_(false), bSetValuesHadErrors_(false)
{
    // Check that user has not provided incorrect values for vectors.
    if (hasFlag(efOption_Vector) && (minValueCount_ > 1 || maxValueCount_ < 1))
    {
        GMX_THROW(APIError("Inconsistent value counts for vector values"));
    }

    if (settings.name_ != NULL)
    {
        name_  = settings.name_;
    }
    descr_ = settings.createDescription();
    setFlag(efOption_ClearOnNextSet);
}

AbstractOptionStorage::~AbstractOptionStorage()
{
}

bool AbstractOptionStorage::isBoolean() const
{
    return dynamic_cast<const BooleanOptionStorage *>(this) != NULL;
}

void AbstractOptionStorage::startSource()
{
    setFlag(efOption_ClearOnNextSet);
}

void AbstractOptionStorage::startSet()
{
    GMX_RELEASE_ASSERT(!bInSet_, "finishSet() not called");
    // The last condition takes care of the situation where multiple
    // sources are used, and a later source should be able to reassign
    // the value even though the option is already set.
    if (isSet() && !hasFlag(efOption_MultipleTimes)
        && !hasFlag(efOption_ClearOnNextSet))
    {
        GMX_THROW(InvalidInputError("Option specified multiple times"));
    }
    clearSet();
    bInSet_ = true;
    bSetValuesHadErrors_ = false;
}

void AbstractOptionStorage::appendValue(const std::string &value)
{
    GMX_RELEASE_ASSERT(bInSet_, "startSet() not called");
    try
    {
        convertValue(value);
    }
    catch (...)
    {
        bSetValuesHadErrors_ = true;
        throw;
    }
}

void AbstractOptionStorage::finishSet()
{
    GMX_RELEASE_ASSERT(bInSet_, "startSet() not called");
    bInSet_ = false;
    // We mark the option as set even when there are errors to avoid additional
    // errors from required options not set.
    // TODO: There could be a separate flag for this purpose.
    setFlag(efOption_Set);
    if (!bSetValuesHadErrors_)
    {
        // TODO: Correct handling of the efOption_ClearOnNextSet requires
        // processSet() and/or convertValue() to check it internally.
        // OptionStorageTemplate takes care of it, but it's error-prone if
        // a custom option is implemented that doesn't use it.
        processSet();
    }
    bSetValuesHadErrors_ = false;
    clearFlag(efOption_ClearOnNextSet);
    clearSet();
}

void AbstractOptionStorage::finish()
{
    GMX_RELEASE_ASSERT(!bInSet_, "finishSet() not called");
    processAll();
    if (isRequired() && !(isSet() || hasFlag(efOption_ExplicitDefaultValue)))
    {
        GMX_THROW(InvalidInputError("Option is required, but not set"));
    }
}

void AbstractOptionStorage::setMinValueCount(int count)
{
    GMX_RELEASE_ASSERT(!hasFlag(efOption_MultipleTimes),
                       "setMinValueCount() not supported with efOption_MultipleTimes");
    GMX_RELEASE_ASSERT(count >= 0, "Invalid value count");
    minValueCount_ = count;
    if (isSet() && !hasFlag(efOption_DontCheckMinimumCount)
        && valueCount() < minValueCount_)
    {
        GMX_THROW(InvalidInputError("Too few values"));
    }
}

void AbstractOptionStorage::setMaxValueCount(int count)
{
    GMX_RELEASE_ASSERT(!hasFlag(efOption_MultipleTimes),
                       "setMaxValueCount() not supported with efOption_MultipleTimes");
    GMX_RELEASE_ASSERT(count >= -1, "Invalid value count");
    maxValueCount_ = count;
    if (isSet() && maxValueCount_ >= 0 && valueCount() > maxValueCount_)
    {
        GMX_THROW(InvalidInputError("Too many values"));
    }
}

/********************************************************************
 * OptionInfo
 */

/*! \cond libapi */
OptionInfo::OptionInfo(AbstractOptionStorage *option)
    : option_(*option)
{
}
//! \endcond

OptionInfo::~OptionInfo()
{
}

bool OptionInfo::isSet() const
{
    return option().isSet();
}

bool OptionInfo::isHidden() const
{
    return option().isHidden();
}

bool OptionInfo::isRequired() const
{
    return option().isRequired();
}

const std::string &OptionInfo::name() const
{
    return option().name();
}

const std::string &OptionInfo::description() const
{
    return option().description();
}

const char *OptionInfo::type() const
{
    return option().typeString();
}

int OptionInfo::valueCount() const
{
    return option().valueCount();
}

std::string OptionInfo::formatValue(int i) const
{
    return option().formatValue(i);
}

std::string OptionInfo::formatDefaultValueIfSet() const
{
    return option().formatDefaultValueIfSet();
}

} // namespace gmx
