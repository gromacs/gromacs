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
 * Declares gmx::Option.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTION_H
#define GMX_OPTIONS_OPTION_H

#include <string>
#include <vector>

#include "optionflags.h"

namespace gmx
{

class AbstractErrorReporter;

class AbstractOption;
class AbstractOptionStorage;
class Options;

/*! \internal \brief
 * Describes a single option.
 *
 * This class is used internally by the option module to store properties
 * and values of a single option.
 *
 * \ingroup module_options
 */
class Option
{
    public:
        Option();
        ~Option();

        /*! \brief
         * Initializes the option object with the properties set in the
         * settings object.
         *
         * \param[in] settings  Option settings.
         * \param[in] options   Option collection that will contain the
         *     option.
         * \retval 0 on success.
         *
         * Should only be called by Options::addOption().
         */
        int init(const AbstractOption &settings, Options *options);

        //! Returns true if the option is a boolean option.
        bool isBoolean() const { return hasFlag(efBoolean); }
        //! Returns true if the option is a file option.
        bool isFile() const { return hasFlag(efFile); }
        //! Returns true if the option is a hidden option.
        bool isHidden() const { return hasFlag(efHidden); }
        //! Returns the name of the option.
        const std::string &name() const { return _name; }
        //! Returns the description of the option.
        const std::string &description() const { return _descr; }
        //! Returns the type of the option as a string.
        const char *type() const;
        //! Returns the number of values given for the option.
        int valueCount() const;
        //! Returns the i'th value of the option as a string.
        std::string formatValue(int i) const;

        //! Returns true if the option has been set.
        bool isSet() const { return hasFlag(efSet); }

        /*! \brief
         * Starts adding values from a new source for the option.
         *
         * \retval 0 on success.
         *
         * This marks the vurrent value of the option as a default value,
         * causing next call to startSet() to clear it.  This allows values
         * from the new source to overwrite old values.
         */
        int startSource();
        /*! \brief
         * Starts adding a new set of values for the option.
         *
         * \param[in] errors Error reporter for errors.
         * \retval 0 on success.
         *
         * If the parameter is specified multiple times, startSet() should be
         * called before the values for each instance.
         */
        int startSet(AbstractErrorReporter *errors);
        /*! \brief
         * Adds a new value for the option, converting it from a string.
         *
         * \param[in] value  String value to convert.
         * \param[in] errors Error reporter for errors.
         * \retval 0 if the value was successfully converted and added.
         */
        int appendValue(const std::string &value,
                        AbstractErrorReporter *errors);
        /*! \brief
         * Performs validation and/or actions once a set of values has been
         * added.
         *
         * \param[in] errors Error reporter for errors.
         * \retval 0 on success.
         *
         * If the parameter is specified multiple times, finishSet() should be
         * called after the values for each instance.
         */
        int finishSet(AbstractErrorReporter *errors);
        /*! \brief
         * Performs validation and/or actions once all values have been added.
         *
         * \param[in] errors Error reporter for errors.
         * \retval 0 on success.
         *
         * This function should be called after all values have been provided
         * with appendValue().  Calls finishSet() if needed.
         */
        int finish(AbstractErrorReporter *errors);

    private:
        //! Returns true if the given flag is set for the option.
        bool hasFlag(OptionFlag flag) const { return _flags.test(flag); }

        std::string             _name;
        std::string             _descr;
        OptionFlags             _flags;
        /*! \brief
         * Pointer to the object that validates and stores values.
         *
         * The Option object is always responsible for freeing the storage
         * object once it is no longer needed.
         */
        AbstractOptionStorage  *_storage;

        // Disallow copy and assign.
        Option(const Option &);
        void operator =(const Option &);
};

} // namespace gmx

#endif
