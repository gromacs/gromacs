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

struct ForceProviders;

struct t_inputrec;

namespace gmx
{

class KeyValueTreeObjectBuilder;
class KeyValueTreeObject;
class IKeyValueTreeErrorHandler;
class IKeyValueTreeTransformRules;
class IMDOutputProvider;
class KeyValueTreeObject;

/*! \libinternal \brief
 * Manages the collection of all modules used for mdrun.
 *
 * This class acts as a central place for constructing modules for mdrun
 * and wiring up dependencies between them.  This class should be the only
 * place that contains the full list of modules, although in the future, some
 * code (e.g., in tools) may benefit from the ability to only create one or a
 * few modules and use them.
 *
 * The general idea is that each module takes care of its own data rather than
 * mdrun having to know about all the details of each type of force calculation.
 * Initially this is applied for simple things like electric field calculations
 * but later more complex forces will be supported too.
 *
 * Currently, where the set of modules needs to be accessed, either a pointer
 * to MDModules is passed around, or an instance of IMDOutputProvider or
 * ForceProviders returned from MDModules.  These objects returned from
 * MDModules call the corresponding methods in the relevant modules.
 * In the future, some additional logic may need to be introduced at
 * the call sites that can also influence the signature of the methods,
 * similar to what ForceProviders already does for force computation.
 *
 * The assignOptionsToModules() and adjustInputrecBasedOnModules() methods of
 * this class also take responsibility for wiring up the options (and their
 * defaults) for each module.
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
         * Initializes a transform from mdp values to sectioned options.
         *
         * \see IMdpOptionProvider::initMdpTransform()
         *
         * Initializes the combined transform from all modules.
         */
        void initMdpTransform(IKeyValueTreeTransformRules *rules);

        /*! \brief Initializes a builder of flat mdp-style key-value pairs
         * suitable for output.
         *
         * If used as input to initMdpTransform(), the key-value pairs
         * resulting from this function would leave the module
         * settings unchanged.
         *
         * Once the transition from mdp to key-value input is
         * complete, this method will probably not exist.
         */
        void buildMdpOutput(KeyValueTreeObjectBuilder *builder);

        /*! \brief
         * Sets input parameters from `params` for each module.
         *
         * \param[in]  params  Contains keys and values from user
         *     input (and defaults) to configure modules that have
         *     registered options with those keys.
         * \param[out] errorHandler  Called to report errors.
         */
        void assignOptionsToModules(const KeyValueTreeObject  &params,
                                    IKeyValueTreeErrorHandler *errorHandler);

        /*! \brief
         * Normalizes inputrec parameters to match current code version.
         *
         * This orders the parameters in `ir->param` to match the current code
         * and adds any missing defaults.  It also throws an error if the
         * inputrec contains parameters that are not recognized by any module.
         */
        void adjustInputrecBasedOnModules(t_inputrec *ir);

        /*! \brief
         * Returns an interface for initializing and finalizing output for modules.
         */
        IMDOutputProvider *outputProvider();
        /*! \brief
         * Returns an object for computing forces from the modules.
         */
        ForceProviders *initForceProviders();

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
