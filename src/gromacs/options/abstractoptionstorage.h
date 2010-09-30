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
/*! \libinternal \file
 * \brief
 * Declares gmx::AbstractOptionStorage.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_ABSTRACTOPTIONSTORAGE_H
#define GMX_OPTIONS_ABSTRACTOPTIONSTORAGE_H

#include <string>

#include "optionflags.h"

namespace gmx
{

class AbstractErrorReporter;
class AbstractOption;
class Option;
class Options;

/*! \libinternal \brief
 * Pure interface for converting, validating, and storing option values.
 *
 * The OptionStorageTemplate template provides basic functionality for this
 * interface, only leaving some of the functions to implement in derived
 * classes.
 * It also integrates well with option settings objects derived from
 * OptionTemplate.
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class AbstractOptionStorage
{
    public:
        virtual ~AbstractOptionStorage();

        /*! \brief
         * Initializes the storage object from the settings object.
         *
         * \param[in] settings  Option settings.
         * \param[in] options   Option collection that will contain the
         *     option.
         * \retval 0 on success.
         */
        int init(const AbstractOption &settings, Options *options);

        /*! \brief
         * Returns a short string describing the type of the option.
         *
         * The caller is free to discard the returned string.
         */
        virtual const char *typeString() const = 0;
        /*! \brief
         * Removes all values from the storage.
         *
         * This function is called also before the first value is added,
         * allowing the storage to set a default value in initialization.
         */
        virtual void clear() = 0;
        /*! \brief
         * Adds a new value, converting it from a string.
         *
         * \param[in] value  String value to convert.
         * \param[in] errors Error reporter object.
         * \retval 0 if the value was successfully converted.
         *
         * This function may be called multiple times if the underlying
         * option is defined to accept multiple values.
         */
        virtual int appendValue(const std::string &value,
                                AbstractErrorReporter *errors) = 0;
        /*! \brief
         * Performs validation and/or actions once a set of values has been
         * added.
         *
         * \param[in] nvalues  Number of values added since the previous call
         *      to finishSet().
         * \param[in] errors Error reporter object.
         * \retval 0 on success.
         *
         * This function may be called multiple times of the underlying option
         * can be specified multiple times.
         */
        virtual int finishSet(int nvalues, AbstractErrorReporter *errors) = 0;
        /*! \brief
         * Performs validation and/or actions once all values have been added.
         *
         * \param[in] errors Error reporter object.
         * \retval 0 if final option values are valid.
         *
         * This function is always called once.
         */
        virtual int finish(AbstractErrorReporter *errors) = 0;
        /*! \brief
         * Returns the number of option values added so far.
         */
        virtual int valueCount() const = 0;
        /*! \brief
         * Returns the i'th value formatted as a string.
         */
        virtual std::string formatValue(int i) const = 0;

    protected:
        //! Creates an uninitialized storage object.
        AbstractOptionStorage();

        //! Returns true if the given flag is set.
        bool hasFlag(OptionFlag flag) const { return _flags.test(flag); }
        //! Sets the given flag.
        void setFlag(OptionFlag flag) { return _flags.set(flag); }
        //! Clears the given flag.
        void clearFlag(OptionFlag flag) { return _flags.clear(flag); }

        //! Returns the minimum number of values required in one set.
        int minValueCount() const { return _minValueCount; }
        //! Returns the maximum allowed number of values in one set (-1 = no limit).
        int maxValueCount() const { return _maxValueCount; }

        //! Returns the Options object that houses the option.
        Options &hostOptions() { return *_options; }
        //! \copydoc hostOptions()
        const Options &hostOptions() const { return *_options; }

        /*! \brief
         * Increments the number of values for the current set.
         *
         * \retval 0 on success.
         * \retval ::eeInvalidInput if the maximum value count has been reached.
         */
        int incrementValueCount();

    private:
        //! Flags for the option.
        OptionFlags             _flags;
        //! Minimum number of values required (in one set).
        int                     _minValueCount;
        //! Maximum allowed number of values (in one set), or -1 if no limit.
        int                     _maxValueCount;
        //! Number of values added so far to the current set, or -1 if not in one.
        int                     _currentValueCount;
        //! Parent Options object.
        Options                *_options;

        /*! \brief
         * Needed to access value count attributes from the parent Option object.
         *
         * There is no clear division between the functionality in the Option
         * object and the associated AbstractOptionStorage, and both need to
         * track the number of values, so they are stored here to make all
         * dependencies one-way.
         */
        friend class Option;

        // Disallow copy and assign.
        AbstractOptionStorage(const AbstractOptionStorage &);
        void operator =(const AbstractOptionStorage &);
};

} // namespace gmx

#endif
