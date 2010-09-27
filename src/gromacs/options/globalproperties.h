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

namespace gmx
{

class Options;
class SelectionCollection;

/*! \libinternal \brief
 * ID numbers for global properties.
 */
enum OptionGlobalPropertyId
{
    eogpTimeScaleFactor,
    eogpSelectionCollection,
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
 * \inlibraryapi
 * \ingroup module_options
 */
class OptionsGlobalProperties
{
    public:
        //! Request for a global property to be used.
        void request(OptionGlobalPropertyId id)
        {
            _usedProperties |= (1<<id);
        }

        //! Set the selection collection for selection option output.
        void setSelectionCollection(SelectionCollection *sc)
        {
            _selectionCollection = sc;
        }

        //! Returns the scaling factor to get times in ps.
        double timeScaleFactor() const;
        //! Returns the selection collection.
        SelectionCollection *selectionCollection() const
        {
            return _selectionCollection;
        }

    private:
        OptionsGlobalProperties();

        //! Returns true if request() has been called for the given property.
        bool isPropertyUsed(OptionGlobalPropertyId id) const
        {
            return _usedProperties & (1<<id);
        }
        void addDefaultOptions(Options *options);

        unsigned long           _usedProperties;
        int                     _timeUnit;
        SelectionCollection    *_selectionCollection;

        friend class Options;

        // Disallow copy and assign.
        OptionsGlobalProperties(const OptionsGlobalProperties &);
        void operator =(const OptionsGlobalProperties &);
};

} // namespace gmx

#endif
