/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
 * \brief Declares classes for components that collaborate to
 * accumulate simulation variables across PP ranks during an MD step.
 *
 * The eventual form of collaboration is described in the following
 * interaction diagram, though the present implementation has to also
 * collaborate with compute_globals(), which is not shown.
 *
 * \msc
   wordwraparcs=true;

   runner,
   builder,
   simulationLoop,
   clients,
   accumulateGlobals;

   builder box builder [ label="builds the accumulateGlobals" ];

   runner => builder [ label="constructs" ];
   runner => clients [ label="constructs many" ];
   clients => builder [ label="register themselves with" ];
   --- [ label="module setup phase completes" ];
   runner => builder [ label="prepares for\n reduction in\n simulationLoop" ];
   builder => clients [ label="queries for requirements" ];
   clients => builder [ label="return requirements" ];
   builder box builder [ label="allocates contiguous storage" ];
   builder => clients [ label="notifies of view matching requirements" ];
   builder => accumulateGlobals [ label="builds" ];
   accumulateGlobals => simulationLoop [ label="is transferred to" ];
   builder => runner [ label="returns" ];
   builder box builder [ label="inactive hereafter" ];
   --- [ label="loop setup phase completes" ];
   runner => simulationLoop [ label="starts" ];
   simulationLoop => clients [ label="calls all" ];
   clients box clients [ label="do work, generating reduction inputs" ];
   clients => accumulateGlobals [ label="fill view\n with reduction\n inputs" ];
   clients => accumulateGlobals [ label="notify when\n reduction is\n required" ];
   clients => simulationLoop [ label="returns" ];
   simulationLoop => accumulateGlobals [ label="calls" ];
   accumulateGlobals box accumulateGlobals [ label="do MPI reduction if required" ];
   accumulateGlobals => clients [ label="in debug mode,\n notify reduction\n is complete" ];
   accumulateGlobals => simulationLoop [ label="returns" ];
   simulationLoop => clients [ label="calls all" ];
   clients box clients [ label="do work using reduction outputs" ];
   clients box clients [ label="clear reduction inputs" ];
   clients => simulationLoop [ label="returns" ];
   simulationLoop box simulationLoop [ label="iterates" ];

 * \endmsc
 *
 * Client objects must implement IAccumulateGlobalsClient and then be
 * registered with an instance of AccumulateGlobalsBuilder. Once
 * setup is complete, all the registered client objects will be
 * queried by \c getNumGlobalsRequired() for the number of values they
 * wish to accumulate. In general, that number will vary when
 * different modules or their options are active, but it is expected
 * to be constant for the lifetime of the simulation. The
 * AccumulateGlobalsBuilder will allocate enough contiguous memory for
 * all clients, which is thus efficient for in-place MPI reduction.
 *
 * In debug mode, the clients are notified after communication
 * completes, so that correctness can be verified.
 *
 * In principle, a simulation could have more than one
 * AccumulateGlobals object, which could be useful if there is need to
 * accumulate at multiple different points during a step. This option
 * should be exploited sparingly, as reduction across all PP ranks
 * requires expensive synchronization.
 *
 * \todo Consider what options may exist for improving the ability of
 * the code to be self-checking, at least in debug mode.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_ACCUMULATEGLOBALS_H
#define GMX_MDRUNUTILITY_ACCUMULATEGLOBALS_H

#include <vector>

#include "gromacs/compat/pointers.h"
#include "gromacs/mdrunutility/iaccumulateglobalsclient.h"
#include "gromacs/utility/arrayref.h"

namespace gmx
{

class AccumulateGlobalsBuilder;

/*! \libinternal
 * \brief Class to manage the buffer used for MPI reduction of simulation globals.
 *
 * \todo Once all/enough modules have been converted to use this
 * framework, then this class could also take the responsibility of
 * holding a communicator and doing the MPI across collaborating
 * ranks, e.g. with PP duty. Then several methods including
 * getReductionView() will no longer be needed.
 */
class AccumulateGlobals
{
    // Permit the builder access to do its job.
    friend class AccumulateGlobalsBuilder;
    public:
        /*! \brief Permits clients to notify that MPI reduction is required this step.
         *
         * \param[in]  clientNotifying  The client making this notification passes its
         *                              this pointer, so in a debug build it will be
         *                              notified after accumulate has completed.
         */
        void notifyReductionRequired(compat::not_null<IAccumulateGlobalsClient *> clientNotifying);
        //! Whether any client requires that the current contents are reduced.
        bool isReductionRequired() const;
        //! Getter for the view to be reduced over MPI.
        ArrayRef<double> getReductionView();
        /*! \brief Called after MPI reduction is complete to permit
         * clients to check their logic.
         *
         * Also clears the flag that requires reduction, in preparation
         * for next MD step. */
        void notifyClientsAfterCommunication();
    private:
        /*! \brief The buffer used for MPI reduction to accumulate global values. */
        std::vector<double>                     values_;
        //! The subset of clients to notify after this reduction completes.
        std::vector < compat::not_null < IAccumulateGlobalsClient * >> clientsToNotifyThisStep_;
        //! Whether reduction is required.
        bool reductionRequired_;
};

/*! \libinternal
 * \brief Class to manage building the buffer used for MPI reduction
 * of simulation globals.
 *
 * Permits multiple modules to register themselves as clients of the
 * object that will be built, and later notifies them of the view to
 * the memory they should use for accumulating a global value.
 *
 * The lifetime of this builder should exceed that of the
 * AccumulateGlobals built by it, if any use of the clients_ view in
 * that AccumulateGlobals occurs. Likewise the lifetime of the clients
 * must exceed the lifetime of both AccumulateGlobals and
 * AccumulateGlobalsBuilder.
 */
class AccumulateGlobalsBuilder
{
    public:
        //! Registers a client that may need to reduce global values.
        void registerClient(compat::not_null<IAccumulateGlobalsClient *> newClient);
        /*! \brief Determine and fulfil client requirements for an
         * AccumulateGlobals.
         *
         * Calls the registered clients to learn how much memory is required,
         * allocates memory, and notifies each client of the view they can
         * later use. */
        AccumulateGlobals build() const;
    private:
        //! The registered clients.
        std::vector < compat::not_null < IAccumulateGlobalsClient * >> clients_;
};

} // namespace gmx

#endif
