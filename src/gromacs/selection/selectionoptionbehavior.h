/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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
 * Declares gmx::SelectionOptionBehavior.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_selection
 */
#ifndef GMX_SELECTION_SELECTIONOPTIONBEHAVIOR_H
#define GMX_SELECTION_SELECTIONOPTIONBEHAVIOR_H

#include <string>

#include "gromacs/options/ioptionsbehavior.h"
#include "gromacs/utility/classhelpers.h"

struct gmx_mtop_t;

namespace gmx
{

class IOptionsContainer;
class Options;
class SelectionCollection;

/*! \brief
 * Provides topology information to SelectionOptionBehavior.
 *
 * Modules that use SelectionOptionBehavior need to implement this interface
 * to provide functionality to load topology information for use with the
 * selections.
 *
 * If future need arises to use similar information elsewhere, this can be
 * moved to, e.g., the topology module, but for now it is here for simplicity.
 * Also, if/when there will be more modules that use this, we can identify
 * common code from those users and possibly provide a shared implementation
 * (e.g., in the form of another IOptionsBehavior), but currently there are too
 * few users to identify any useful reusable functionality from the callers.
 *
 * \see SelectionCollection::setTopology().
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class ITopologyProvider
{
    public:
        /*! \brief
         * Returns the topology to use.
         *
         * \param[in] required  Whether the topology is required by the caller.
         *
         * Can return NULL if \p required is `false` and the topology is not
         * provided by the user.
         * If \p required is `true`, should throw an error if the topology
         * cannot be loaded.
         *
         * This method may get called multiple times, potentially with
         * different values of \p required.  Subsequent calls should just
         * return the same topology that was loaded in the first call.
         */
        virtual gmx_mtop_t *getTopology(bool required) = 0;
        /*! \brief
         * Returns the number of atoms.
         *
         * This method is only called if getTopology() returns NULL.
         * It should return the number of atoms that at most need to be
         * selected by the selections.
         */
        virtual int getAtomCount()                     = 0;

    protected:
        virtual ~ITopologyProvider();
};

/*! \brief
 * Options behavior to allow using SelectionOptions.
 *
 * This behavior wraps SelectionOptionManager, as well as all calls to
 * SelectionCollection up to and including selection compilation.
 *
 * To use selections through SelectionOption options in a
 * ICommandLineOptionsModule, you need to
 *  * create a SelectionCollection object,
 *  * implement ITopologyProvider to return the topology or number of atoms
 *    to be used with the selections,
 *  * create and add a SelectionOptionBehavior, and call initOptions() to add
 *    common options for selection control,
 *  * use SelectionOption options to specify the selections used, and
 *  * evaluate the selections when running the module for the desired sets of
 *    coordinates (see SelectionCollection::evaluate()).
 *
 * The SelectionOptionBehavior provides the following functionalities to the
 * module:
 *  * Creates SelectionOptionManager and manages all calls to it.
 *  * Creates an option to provide an `ndx` file, and loads it when necessary.
 *  * Creates options to control general aspects of selections (see
 *    SelectionCollection::initOptions()), as well as a generic option to
 *    provide selections from a file.
 *  * After all options have been processed, provides an interactive
 *    command-line prompt for any missing selections.
 *  * Compiles the selections.
 *
 * The behavior needs to be added before any options are created.
 *
 * \inpublicapi
 * \ingroup module_selection
 */
class SelectionOptionBehavior : public IOptionsBehavior
{
    public:
        /*! \brief
         * Creates a behavior to use selections.
         *
         * \param[in,out] selections  Selection collection to use.
         * \param[in]     topologyProvider  Callback to load/provide topology
         *     information to selections when required.
         *
         * The methods in \p topologyProvider are called after all options have
         * been parsed and finished, so the caller can, e.g., load the topology
         * from a file specified by a file option.
         */
        SelectionOptionBehavior(SelectionCollection *selections,
                                ITopologyProvider   *topologyProvider);
        ~SelectionOptionBehavior();

        /*! \brief
         * Add common options for controlling selections.
         *
         * This method is separate from the constructor so that the caller can
         * control the order of options better.
         */
        void initOptions(IOptionsContainer *options);

        // From IOptionsBehavior
        virtual void initBehavior(Options *options);
        virtual void optionsFinishing(Options * /*options*/) {}
        virtual void optionsFinished();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};


} // namespace gmx

#endif
