/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::MDModules.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_MDMODULES_H
#define GMX_MDRUNUTILITY_MDMODULES_H

#include "gromacs/utility/classhelpers.h"

struct t_inputrec;

namespace gmx
{

class KeyValueTreeObject;
class IKeyValueTreeErrorHandler;
class IKeyValueTreeTransformRules;

/*! \libinternal \brief
 * Factory for t_inputrec.
 *
 * This class acts as a central place for constructing t_inputrec (and possibly
 * other mdrun structures in the future) and wiring up dependencies between
 * modules that are referenced from these structures.  This class owns all such
 * modules, and needs to remain in existence as long as the returned data
 * structures are in use.  Ideally, it is also the only place that creates
 * instances of these modules (outside test code).
 *
 * The general idea is that each module takes care of its own data rather than
 * mdrun having to know about all the details of each type of force calculation.
 * Initially this is applied for simple things like electric field calculations
 * but later more complex forces will be supported too.
 *
 * The current approach uses t_inputrec and IInputRecExtension to pass
 * references to the modules to other code to avoid changing many function
 * signatures.  Also, the current usage means that nearly every use of
 * t_inputrec (in particular, reading it from mdp or tpr files) needs to be
 * initialized through MDModules for correct functionality.  For the future, a
 * better approach would be to pass around a reference to MDModules instead and
 * call it directly for cases that are not related to t_inputrec functionality.
 * This (and other refactoring) would allow simplifying IInputRecExtension.
 * IForceProvider is the other interface currently used to interact with these
 * modules.  Also, all the places where these interfaces are used should become
 * loops over a container of these interfaces, and/or groups of them (e.g.
 * applied forces), instead of the current single pointer.
 *
 * The assignOptionsToModules() and
 * assignOptionsToModulesFromInputrec() methods of this class also
 * take responsibility for wiring up the options (and their defaults)
 * for each module, respectively for mdp- and tpr-style input of those
 * options.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
class MDModules
{
    public:
        MDModules();
        ~MDModules();

        /*! \brief
         * Returns an initialized t_inputrec structure.
         *
         * The inputrec structure is owned by MDModules and will be destroyed
         * with it.
         */
        t_inputrec *inputrec();
        //! \copydoc t_inputrec *inputrec()
        const t_inputrec *inputrec() const;

        /*! \brief Initializes a transform from mdp values to
         * sectioned options.
         *
         * The transform is specified from a flat KeyValueTreeObject that
         * contains each mdp value as a property, to a structure which is then
         * assigned to the options defined with initMdpOptions().
         *
         * Once the transition from mdp to key-value input is
         * complete, this method will probably not exist.
         */
        void initMdpTransform(IKeyValueTreeTransformRules *rules);

        /*! \brief Use \c mdpOptionValues to set the options (e.g.read
         * from mdp input) for each module.
         *
         * \param[in] mdpOptionValues Contains keys and values from user
         *     input (and defaults) to configure modules that have
         *     registered options with those keys.
         * \param[out] errorHandler  Called to report errors. */
        void assignOptionsToModulesFromMdp(const KeyValueTreeObject  &mdpOptionValues,
                                           IKeyValueTreeErrorHandler *errorHandler);

        /*! \brief
         * Initializes modules based on inputrec values read from tpr file.
         *
         * This needs to be called after read_tpx_state() if the modules need
         * to be accessed.
         */
        void assignOptionsToModulesFromTpr();

        /*! \brief
         * Normalizes inputrec parameters to match current code version.
         *
         * This orders the parameters in inputrec to match the current code and
         * adds any missing defaults.
         */
        void adjustInputrecBasedOnModules();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
