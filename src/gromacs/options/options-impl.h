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
 * Declares private implementation class for gmx::Options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONS_IMPL_H
#define GMX_OPTIONS_OPTIONS_IMPL_H

#include <string>
#include <vector>

#include "options.h"

namespace gmx
{

class OptionsGlobalProperties;
class Option;

/*! \internal \brief
 * Private implementation class for Options.
 *
 * Note that in addition to Options, the OptionsAssigner and OptionsIterator
 * classes also directly access this class.
 *
 * \ingroup module_options
 */
class Options::Impl
{
    public:
        //! Internal flags for the implementation.
        enum Flag {
            efHasFileOptions    = 1<<0,
            efHasNonFileOptions = 1<<1,
        };
        //! Type to hold flags.
        typedef unsigned int Flags;

        //! Convenience type for list of sections.
        typedef std::vector<Options *> SubSectionList;
        //! Convenience type for list of options.
        typedef std::vector<Option *> OptionList;

        //! Sets the name and title.
        Impl(const char *name, const char *title);
        ~Impl();

        /*! \brief
         * Finds a subsection by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found subsection, or NULL if not found.
         */
        Options *findSubSection(const char *name) const;
        /*! \brief
         * Finds an option by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found option, or NULL if not found.
         */
        Option *findOption(const char *name) const;

        /*! \brief
         * Calls Option::startSource() for all options, including subsections.
         *
         * \returns 0 on success, or the first non-zero return value of the
         *      called functions.
         */
        int startSource();

        //! Name for the Options object.
        std::string             _name;
        //! Description title for the Options object.
        std::string             _title;
        //! Full description for the Options object.
        std::string             _description;
        /*! \brief
         * List of subsections, in insertion order.
         *
         * This container contains only references to external objects; memory
         * management is performed elsewhere.
         */
        SubSectionList          _subSections;
        /*! \brief
         * List of options, in insertion order.
         *
         * All objects in this container are owned by this object, and are
         * freed in the destructor.
         */
        OptionList              _options;
        //! Options object that contains this object as a subsection, or NULL.
        Options                *_parent;
        /*! \brief
         * Object that contains global properties, or NULL if \a _parent != NULL.
         *
         * This object is always owned by the Options object.
         * For subsections, the global properties are kept in the parent, and
         * this pointer is NULL.
         */
        OptionsGlobalProperties *_globalProperties;
        /*! \brief
         * Flags for extra information.
         */
        Flags                   _flags;
};

} // namespace gmx

#endif
