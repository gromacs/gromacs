/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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

#include "gromacs/math/paddedvector.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

#include <cstdio>

struct t_commrec;
struct gmx_output_env_t;
struct t_filenm;
struct t_forcerec;
struct t_inputrec;
struct t_mdatoms;

namespace gmx
{

class KeyValueTreeObject;
class IKeyValueTreeErrorHandler;
class IKeyValueTreeTransformRules;

/*! \libinternal \brief
 * Owner of modules used by mdrun.
 *
 * This class acts as a central place for owning and managing data
 * structures for mdrun modules, and wiring up dependencies between
 * them. (It also has responsibility for constructing t_inputrec (and
 * possibly other mdrun structures in the future) while refactoring to
 * eliminate it takes place.) This class owns all such modules, and
 * needs to remain in existence as long as the returned data
 * structures are in use.  Ideally, it is also the only place that
 * creates instances of these modules (outside test code).
 *
 * The general idea is that each module takes care of its own data rather than
 * mdrun having to know about all the details of each type of force calculation.
 * Initially this is applied for simple things like electric field calculations
 * but later more complex forces will be supported too.
 *
 * Functionality common to (what will eventually be) the set of mdrun
 * modules should be invoked by calling methods of MDModules.  The
 * current usage means that nearly every use of t_inputrec (in
 * particular, reading it from mdp or tpr files) needs to be
 * initialized through MDModules for correct functionality.
 * IForceProvider is the other interface currently used to interact
 * with these modules.
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
class MDModules
{
    public:
        MDModules();
        ~MDModules();

        /*! \brief Returns an initialized t_inputrec structure.
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
        /*! \brief Use \c tree to set the module options.
         *
         * \param[in] tree Contains keys and values from user input
         *                 (and defaults) to configure modules that
         *                 have registered options with those keys
         * \param[out] errorHandler  Called to report errors. */
        void assignOptionsToModules(KeyValueTreeObject &&tree,
                                    IKeyValueTreeErrorHandler *errorHandler);

        /*! \brief Print parameters belonging to modules to \c fp.
         *
         * Used in writing mdrun log files and gmx dump.
         *
         * \param[in] fp     File pointer
         * \param[in] indent Initial indentation level for printing
         */
        void printParameters(FILE *fp, int indent) const;

        /*! \brief Initiate any module-specific handling for writing
         * mdrun output files.
         *
         * \param[in] fplog File pointer for log messages
         * \param[in] nfile Number of files
         * \param[in] fnm   Array of filenames and properties
         * \param[in] bAppendFiles Whether or not we should append to files
         * \param[in] oenv  The output environment for xvg files
         */
        void initOutput(FILE *fplog,
                        int nfile,
                        const t_filenm fnm[],
                        bool bAppendFiles,
                        const gmx_output_env_t *oenv);

        //! Finalize module-specific handling for mdrun output.
        void finishOutput();

        /*! \brief Compare matching modules from two MDModules.
         *
         * Used in gmx check.
         *
         * \todo The current implementation only makes it sensible to
         * compare matching tpr versions, ie aware of the same
         * modules.
         *
         * \param[in]    fp     File pointer
         * \param[in]    other  Another MDModules
         * \param[in]    reltol Relative tolerance
         * \param[in]    abstol Absolute tolerance
         */
        void compare(FILE *fp,
                     const MDModules *other,
                     real reltol,
                     real abstol) const;

        /*! \brief Broadcast parameters from each module from master to all ranks.
         *
         * \param[in] cr  Communication record for parallel operations
         */
        void broadCast(const t_commrec *cr);

        /*! \brief Initialize data still stored in forcerec for all
         * modules that also implement IForceProvider.
         *
         * \todo Eventually that data will be stored in the modules
         * themselves, and this method may no longer be required.
         */
        void initForcerec(t_forcerec *fr);

        /*! \brief Compute forces for modules.
         *
         * \todo The function signature is specific for electric
         * fields, and needs to be generalized.
         *
         * \param[in]    cr      Communication record for parallel operations
         * \param[in]    mdatoms Atom information
         * \param[inout] force   The forces
         * \param[in]    t       The actual time in the simulation (ps)
         */
        void calculateForces(const t_commrec  *cr,
                             const t_mdatoms  *mdatoms,
                             PaddedRVecVector *force,
                             double            t);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
