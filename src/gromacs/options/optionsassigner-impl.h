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
 * Declares private implementation class for gmx::OptionsAssigner.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_OPTIONSASSIGNER_IMPL_H
#define GMX_OPTIONS_OPTIONSASSIGNER_IMPL_H

#include <vector>

#include "optionsassigner.h"

namespace gmx
{

class AbstractOptionStorage;
class Options;

/*! \internal \brief
 * Private implementation class for OptionsAssigner.
 *
 * \ingroup module_options
 */
class OptionsAssigner::Impl
{
    public:
        //! Possible flags for controlling assignment behavior.
        enum Flag
        {
            //! Recognize boolean option "name" also as "noname".
            efAcceptBooleanNoPrefix     = 1<<0,
            //! Look for options in all sections, not just the current one.
            efNoStrictSectioning        = 1<<1,
        };
        //! Sets the option object to assign to.
        Impl(Options *options);

        //! Sets or clears the given flag.
        void setFlag(Flag flag, bool bSet);

        //! Returns true if the given flag is set.
        bool hasFlag(Flag flag) const { return _flags & flag; }
        //! Returns true if a subsection has been set.
        bool inSubSection() const { return _sectionStack.size() > 1; }
        //! Returns the Options object for the current section.
        Options &currentSection() const { return *_sectionStack.back(); }
        /*! \brief
         * Finds an option by the given name.
         *
         * \param[in] name  Name of the option to look for.
         * \returns Pointer to the found option, or NULL if none found.
         *
         * This function takes into account the flags specified, and may change
         * the internal state of the assigner to match the option found.
         * If no option is found, the internal state is not modified.
         */
        AbstractOptionStorage *findOption(const char *name);

        //! Options object to assign to.
        Options                &_options;
        //! Flags that control assignment behavior.
        unsigned long           _flags;
        /*! \brief
         * List of (sub)sections being assigned to.
         *
         * The first element always points to \a _options.
         */
        std::vector<Options *>  _sectionStack;
        //! Current option being assigned to, or NULL if none.
        AbstractOptionStorage  *_currentOption;
        //! Number of values assigned so far to the current option.
        int                     _currentValueCount;
        //! If true, a "no" prefix was given for the current boolean option.
        bool                    _reverseBoolean;
};

} // namespace gmx

#endif
