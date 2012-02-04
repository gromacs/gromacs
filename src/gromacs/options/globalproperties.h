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
 * Declares gmx::OptionsGlobalProperties.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_GLOBALPROPERTIES_H
#define GMX_OPTIONS_GLOBALPROPERTIES_H

#include <typedefs.h>

namespace gmx
{

class Options;

/*! \libinternal \brief
 * ID numbers for global properties.
 */
enum OptionGlobalPropertyId
{
    eogpTimeScaleFactor,
    eogpPlotFormat,
};

/*! \libinternal \brief
 * Describes global properties of an Options collection.
 *
 * These properties are used to implement features that require all options of
 * a certain type to access some global data.
 * For example, if there are options that specify times, and in addition an
 * option that specifies the unit for these times, all the time options need to
 * know the scaling factor to get the time in internal units.
 *
 * \todo
 * There are things in this class that would be
 * better placed in the analysisdata module (for selecting plot formats).
 * It should be considered whether this should be implemented in some other way
 * (see Redmine issue #839).
 *
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsGlobalProperties
{
    public:
        ~OptionsGlobalProperties();

        //! Request for a global property to be used.
        void request(OptionGlobalPropertyId id)
        {
            _usedProperties |= (1<<id);
        }

        //! Returns the scaling factor to get times in ps.
        double timeScaleFactor() const;
        /*! \brief
         * Returns an output environment structure for interfacing with old
         * code.
         *
         * Currently, the returned structure is always filled with default
         * values for most fields.
         *
         * \deprecated
         */
        output_env_t output_env() const
        {
            return _oenv;
        }

    private:
        OptionsGlobalProperties();

        //! Returns true if request() has been called for the given property.
        bool isPropertyUsed(OptionGlobalPropertyId id) const
        {
            return _usedProperties & (1<<id);
        }
        /*! \brief
         * Adds options for setting requested global properties.
         *
         * \param[in,out] options Options to which the global property options
         *      are added.
         *
         * If a global property has been requested and it can be set/customized
         * by the user, this method adds the necessary option to \p options.
         *
         * This method performs the real work of Options::addDefaultOptions().
         */
        void addDefaultOptions(Options *options);
        /*! \brief
         * Initializes variables dependent on global properties.
         *
         * This method should be called after the values for the options
         * generated with addDefaultOptions() have been set.
         */
        void finish();

        unsigned long           _usedProperties;
        int                     _timeUnit;
        int                     _plotFormat;
        output_env_t            _oenv;

        friend class Options;

        // Disallow copy and assign.
        OptionsGlobalProperties(const OptionsGlobalProperties &);
        void operator =(const OptionsGlobalProperties &);
};

} // namespace gmx

#endif
