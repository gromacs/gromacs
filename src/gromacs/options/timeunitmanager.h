/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#include "gromacs/utility/gmxassert.h"

namespace gmx
{

class Options;

/*! \brief
 * Time values for TimeUnitManager.
 *
 * \if internal
 * Currently, this should match with the time_unit_t enum defined in oenv.h
 * except that there is no NULL first item in this enum.
 * \endif
 *
 * \inpublicapi
 */
enum TimeUnit
{
    eTimeUnit_fs, //!< Femtoseconds.
    eTimeUnit_ps, //!< Picoseconds.
    eTimeUnit_ns, //!< Nanoseconds.
    eTimeUnit_us, //!< Microseconds.
    eTimeUnit_ms, //!< Milliseconds.
    eTimeUnit_s   //!< Seconds.
};

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
 * Most of this class is independent of the options implementation.
 * To ease reuse, it could be split such that the generic part is moved to the
 * utility module, and only the options-specific parts left in the options
 * module.
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
        TimeUnit timeUnit() const
        {
            GMX_ASSERT(timeUnit_ >= 0 && timeUnit_ <= eTimeUnit_s,
                       "Time unit index has become out-of-range");
            return static_cast<TimeUnit>(timeUnit_);
        }
        //! Set a new time unit for the manager.
        void setTimeUnit(TimeUnit unit);

        //! Returns a string constant corresponding to the current time unit.
        const char *timeUnitAsString() const;

        //! Returns the scaling factor to convert times to ps.
        double timeScaleFactor() const;
        //! Returns the scaling factor to convert times from ps.
        double inverseTimeScaleFactor() const;

        /*! \brief
         * Sets the time unit in this manager from an environment variable.
         */
        void setTimeUnitFromEnvironment();
        /*! \brief
         * Adds a common option for selecting the time unit.
         *
         * \param[in,out] options Options to which the common option is added.
         * \param[in]     name    Name of the option to add.
         *
         * Adds an enum option to \p options to select the time unit for this
         * manager.
         */
        void addTimeUnitOption(Options *options, const char *name);
        /*! \brief
         * Scales user input values given to time options.
         *
         * \param[in,out] options Options in which to scale times.
         *
         * Scales each time option (see DoubleOption::timeValue()) in
         * \p options such that any user-given values are interpreted as given
         * in the time unit specified by this manager, and scaled to
         * picoseconds.  Programmatically given values (e.g., as default values
         * for the options) are not scaled.
         */
        void scaleTimeOptions(Options *options) const;

    private:
        /*! \brief
         * Currently set time unit for this manager.
         *
         * Type is int to make it possible to use it with
         * StringOption::storeEnumIndex(), but it should always one of the
         * allowed values for TimeUnit.
         */
        int                     timeUnit_;
};

} // namespace gmx

#endif
