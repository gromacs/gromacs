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
/*! \file
 * \brief
 * Declares gmx::SelectionOptionManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTIONMANAGER_H
#define GMX_SELECTION_SELECTIONOPTIONMANAGER_H

#include <string>

#include "../utility/common.h"

namespace gmx
{

class Options;
class SelectionCollection;
class SelectionOptionStorage;

/*! \brief
 * Handles interaction of selection options with other options and user input.
 *
 * This class implements interaction of SelectionOption with
 * SelectionCollection, and also implements features of SelectionOption that
 * require actions outside options parsing.
 * It also implements the coupling between SelectionOption and
 * SelectionFileOption.
 *
 * The main features of this class are:
 *  - convertOptionValue(), which is used to convert string values into
 *    selections for options.
 *  - requestOptionDelayedParsing(), which is called by the internal
 *    implementation of selection options when an option is provided on the
 *    command line without a value.  Such calls are remembered, and the value
 *    for all requested options can be later provided by calling one of
 *    parseRequestedFromStdin(), parseRequestedFromFile() or
 *    parseRequstedFromString().
 *
 * \see setManagerForSelectionOptions()
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOptionManager
{
    public:
        /*! \brief
         * Creates a manager for selection options.
         *
         * \throws  std::bad_alloc if out of memory.
         */
        explicit SelectionOptionManager(SelectionCollection *selections);
        ~SelectionOptionManager();

        /*! \brief
         * Adds a selection option to be managed.
         *
         * \param     storage  Storage object for the option to register.
         * \throws    std::bad_alloc if out of memory.
         *
         * This is only for internal use by the selection module.
         * It is not possible to obtain a SelectionOptionStorage pointer
         * through any public or library API.
         *
         * Strong exception safety.
         */
        void registerOption(SelectionOptionStorage *storage);
        /*! \brief
         * Converts a string value to selections for an option.
         *
         * \param     storage  Storage object to receive the selections.
         * \param[in] value    Value to convert.
         * \throws    std::bad_alloc if out of memory.
         * \throws    InvalidInputError if the selection string is not valid,
         *      or uses a feature not supported by the option.
         *
         * This is only for internal use by the selection module.
         * It is not possible to obtain a SelectionOptionStorage pointer
         * through any public or library API.
         */
        void convertOptionValue(SelectionOptionStorage *storage,
                                const std::string &value);
        /*! \brief
         * Adds a selection option for delayed user input.
         *
         * \param     storage  Storage object for the option to request.
         * \throws    std::bad_alloc if out of memory.
         *
         * This is only for internal use by the selection module.
         * It is not possible to obtain a SelectionOptionStorage pointer
         * through any public or library API.
         *
         * Strong exception safety.
         */
        void requestOptionDelayedParsing(SelectionOptionStorage *storage);

        /*! \brief
         * Parses selection(s) from standard input for options not yet
         * provided.
         *
         * \param[in]  bInteractive Whether the parser should behave
         *      interactively.
         * \throws     unspecified  Can throw any exception thrown by
         *      SelectionCollection::parseFromStdin().
         * \throws     std::bad_alloc if out of memory.
         *
         * This method cooperates with SelectionOption to allow interactive
         * input of requested selections after all options have been processed.
         * It should be called after the Options::finish() method has been
         * called on all options that add selections to this collection.
         * For each required selection option that has not been given, as well
         * as for optional selection options that have been specified without
         * values, it will prompt the user to input the necessary selections.
         */
        void parseRequestedFromStdin(bool bInteractive);
        /*! \brief
         * Parses selection(s) from a file for options not yet provided.
         *
         * \param[in]  filename Name of the file to parse selections from.
         * \throws     unspecified  Can throw any exception thrown by
         *      SelectionCollection::parseFromFile().
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError if
         *      - the number of selections in \p filename doesn't match the
         *        number requested.
         *      - any selection uses a feature that is not allowed for the
         *        corresponding option.
         *      - if there is a request for any number of selections that is
         *        not the last (in which case it is not possible to determine
         *        which selections belong to which request).
         *
         * This method behaves as parseRequestedFromStdin(), with two
         * exceptions:
         *  -# It reads the selections from a file instead of standard input.
         *  -# If no requests are pending, assigns values to all required
         *     options that have not yet been set.
         *
         * This method used to implement SelectionFileOption.
         *
         * \see parseRequestedFromStdin()
         */
        void parseRequestedFromFile(const std::string &filename);
        /*! \brief
         * Parses selection(s) from a string for options not yet provided.
         *
         * \param[in]  str     String to parse.
         * \throws     unspecified  Can throw any exception thrown by
         *      SelectionCollection::parseFromString().
         * \throws     std::bad_alloc if out of memory.
         * \throws     InvalidInputError in same conditions as
         *      parseRequestedFromFile().
         *
         * This method behaves as parseRequestedFromFile(), but reads the
         * selections from a string instead of a file.
         * This method is mainly used for testing.
         *
         * \see parseRequestedFromFile()
         */
        void parseRequestedFromString(const std::string &str);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;

        /*! \brief
         * Needed for handling delayed selection parsing requests.
         */
        friend class SelectionOptionStorage;
};

/*! \brief
 * Set manager for all selection options.
 *
 * Recursively sets the manager to \p manager for all selection options in
 * \p options.
 * Must be called before value assignment starts for \p options.
 *
 * Does not throw.
 *
 * \inpublicapi
 */
void setManagerForSelectionOptions(Options *options,
                                   SelectionOptionManager *manager);

} // namespace gmx

#endif
