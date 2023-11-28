/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 *
 * \brief This file declares functions for inter-rank signalling by mdrun
 *
 * This handles details of responding to termination conditions,
 * coordinating checkpoints, and coordinating multi-simulations.
 *
 * \todo Move this to mdrunutility module alongside gathering
 * multi-simulation communication infrastructure there.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdlib
 */
#ifndef GMX_MDLIB_SIMULATIONSIGNAL_H
#define GMX_MDLIB_SIMULATIONSIGNAL_H

#include <array>

#include "gromacs/utility/real.h"

struct gmx_multisim_t;
struct t_commrec;

//! Kinds of simulation conditions to signal about.
enum
{
    eglsCHKPT,
    eglsSTOPCOND,
    eglsRESETCOUNTERS,
    eglsNR
};

namespace gmx
{

template<typename T>
class ArrayRef;

/*!
 * \brief POD-style object used by mdrun ranks to set and
 * receive signals within and between simulations.
 *
 * Keep in mind that the values of signals are transmitted to other
 * ranks through an MPI_Reduce after casting them to a real (so the
 * signals can be sent together with other data). This means that the
 * only meaningful values are positive, negative or zero.
 *
 * isLocal permits (for example) replica-exchange to require that any
 * checkpointing is synchronized across all simulations, by setting
 * isLocal to false, so that the trigger for action is set only when
 * inter-simulation signalling happens. Replica-exchange can
 * coordinate this at run time when a SimulationSignaller is made. */
class SimulationSignal
{
public:
    //! Constructor
    SimulationSignal(bool isSignalLocal = true) : sig(0), set(0), isLocal(isSignalLocal) {}
    //! The signal set by this rank in do_md().
    signed char sig;
    //! The communicated signal that triggers action, which will be equal for all ranks, once communication has occurred.
    signed char set;
    //! Is the signal in one simulation independent of other simulations?
    bool isLocal;
};

//! Convenience typedef for the group of signals used.
typedef std::array<SimulationSignal, eglsNR> SimulationSignals;

/*!
 * \brief Object used by mdrun ranks to signal to each other at this step.
 *
 * This object has responsibility to read signal values from \c gs,
 * coordinate communication within and perhaps between simulations,
 * and set result signal values in \c gs as appropriate.
 *
 * It is intended to have a very short lifetime, so should remain easy
 * to construct and destruct on the stack just when the global
 * communication occurs. */
class SimulationSignaller
{
public:
    //! Constructor
    SimulationSignaller(SimulationSignals*    signals,
                        const t_commrec*      cr,
                        const gmx_multisim_t* ms,
                        bool                  doInterSim,
                        bool                  doIntraSim);
    /*! \brief Return a reference to an array of signal values to communicate.
     *
     * \return If intra-sim signalling will take place, fill and
     * return a reference to the array of reals in which signals
     * will be communicated with the signal values to be
     * sent. Otherwise return a EmptyArrayRef. */
    gmx::ArrayRef<real> getCommunicationBuffer();
    /*! \brief Handle inter-simulation signal communication.
     *
     * If an inter-simulation signal should be handled, communicate between
     * simulation-main ranks, then propagate from the mains to the
     * rest of the ranks for each simulation. It is the responsibility of
     * the calling code to ensure that any necessary intra-simulation
     * signalling has already occurred, e.g. in global_stat(). */
    void signalInterSim();
    /*! \brief Propagate signals when appropriate.
     *
     * Always propagate an mdrun signal value when doing
     * inter-simulation signalling; otherwise, propagate it only
     * if should be propagated within this simulation,
     * ie. locally. See documentation of SimulationSignal for
     * details. */
    void setSignals();
    //! Convenience wrapper that calls signalInterSim() then setSignals().
    void finalizeSignals();

private:
    //! Source and sink for mdrun signals
    SimulationSignals* signals_;
    //! Communication object.
    const t_commrec* cr_;
    //! Multi-sim handler.
    const gmx_multisim_t* ms_;
    //! Do inter-sim communication at this step.
    bool doInterSim_;
    //! Do intra-sim communication at this step.
    bool doIntraSim_;
    //! Buffer for MPI communication.
    std::array<real, eglsNR> mpiBuffer_;
};

} // namespace gmx

#endif
