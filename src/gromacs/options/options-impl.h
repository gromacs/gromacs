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

#include "abstractoption.h"
#include "options.h"

namespace gmx
{

class AbstractOptionStorage;

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
        //! Convenience type for list of sections.
        typedef std::vector<Options *> SubSectionList;
        //! Convenience type for list of options.
        typedef std::vector<AbstractOptionStoragePointer> OptionList;

        //! Sets the name and title.
        Impl(const char *name, const char *title);
        ~Impl();

        /*! \brief
         * Finds a subsection by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found subsection, or NULL if not found.
         *
         * Does not throw.
         */
        Options *findSubSection(const char *name) const;
        /*! \brief
         * Finds an option by name.
         *
         * \param[in] name  Name to search for.
         * \returns Pointer to the found option, or NULL if not found.
         *
         * Does not throw.
         */
        AbstractOptionStorage *findOption(const char *name) const;

        /*! \brief
         * Calls AbstractOptionStorage::startSource() for all options,
         * including subsections.
         *
         * Does not throw.
         */
        void startSource();

        //! Name for the Options object.
        std::string name_;
        //! Description title for the Options object.
        std::string title_;
        //! Full description for the Options object.
        std::string description_;
        /*! \brief
         * List of subsections, in insertion order.
         *
         * This container contains only references to external objects; memory
         * management is performed elsewhere.
         */
        SubSectionList subSections_;
        /*! \brief
         * List of options, in insertion order.
         *
         * All objects in this container are owned by this object, and are
         * freed in the destructor.
         */
        OptionList options_;
        //! Options object that contains this object as a subsection, or NULL.
        Options *parent_;
};

} // namespace gmx

#endif
