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
/*! \file
 * \brief
 * Declares gmx::TimeUnitManager.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_TIMEUNITMANAGER_H
#define GMX_OPTIONS_TIMEUNITMANAGER_H

#include <memory>

#include "gromacs/fileio/oenv.h"
#include "gromacs/options/ioptionsbehavior.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class IOptionsContainer;
class Options;

/*! \brief
 * Provides common functionality for time unit conversions.
 *
 * Methods/objects that need to deal with time units can either take a
 * TimeUnitManager object, or they can take a TimeUnit value and construct a
 * TimeUnitManager object internally.
 *
 * Default copy constructor and assignment are used: the copy is an independent
 * object that is initialized with the same time unit as the original.
 *
 * \if internal
 * \todo
 * This class is independent of the options implementation.
 * To ease reuse, it could be moved to the utility module, and only
 * TimeUnitBehavior left here.
 * \endif
 *
 * \inpublicapi
 * \ingroup module_options
 */
class TimeUnitManager
{
public:
    //! Creates a time unit manager with the default (ps) time unit.
    TimeUnitManager();
    //! Creates a time unit manager with the given time unit.
    explicit TimeUnitManager(TimeUnit unit);

    //! Returns the currently selected time unit.
    TimeUnit timeUnit() const { return timeUnit_; }
    //! Set a new time unit for the manager.
    void setTimeUnit(TimeUnit unit);

    //! Returns a string constant corresponding to the current time unit.
    const char* timeUnitAsString() const;

    //! Returns the scaling factor to convert times to ps.
    double timeScaleFactor() const;
    //! Returns the scaling factor to convert times from ps.
    double inverseTimeScaleFactor() const;

private:
    //! Currently set time unit for this manager.
    TimeUnit timeUnit_;
};

/*! \brief
 * Options behavior to add a time unit option.
 *
 * This class provides functionality to add a time unit option that affects the
 * input unit for time options (specified with FloatOption::timeValue() or
 * DoubleOption::timeValue()).  When options are finished, it scales each time
 * option such that any user-given values are interpreted as given in the time
 * unit specified by the user, and scaled to picoseconds.  Programmatically
 * given values (e.g., as default values for the options) are not scaled.
 *
 * \inpublicapi
 * \ingroup module_options
 */
class TimeUnitBehavior : public IOptionsBehavior
{
public:
    TimeUnitBehavior();

    //! Returns the current time unit.
    TimeUnit timeUnit() const { return timeUnit_; }
    //! Sets the time unit.
    void setTimeUnit(TimeUnit unit);

    /*! \brief
     * Sets a storage location for the selected time unit.
     *
     * \param[in] store  Location that will receive the selected time unit.
     *
     * \p *store will be set to the time unit selected by the user (or
     * programmatically).  The value is guaranteed to be set once the
     * options have been finished.
     */
    void setTimeUnitStore(TimeUnit* store);

    /*! \brief
     * Sets the default time unit from an environment variable.
     *
     * This should be called before addTimeUnitOption() for consistent
     * behavior.
     */
    void setTimeUnitFromEnvironment();
    /*! \brief
     * Adds a common option for selecting the time unit.
     *
     * \param[in,out] options Options to which the common option is added.
     * \param[in]     name    Name of the option to add.
     *
     * Adds an enum option to \p options to select the time unit for this
     * behavior.
     */
    void addTimeUnitOption(IOptionsContainer* options, const char* name);

    // From IOptionsBehavior
    void initBehavior(Options* /*options*/) override {}
    void optionsFinishing(Options* options) override;
    void optionsFinished() override {}

private:
    TimeUnit  timeUnit_;
    TimeUnit* timeUnitStore_;

    GMX_DISALLOW_COPY_AND_ASSIGN(TimeUnitBehavior);
};

} // namespace gmx

#endif
