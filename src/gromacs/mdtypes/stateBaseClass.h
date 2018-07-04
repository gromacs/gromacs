/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
 *
 * \brief
 * This file contains the definition of the microstate of the simulated system
 *
 * History of observables that needs to be checkpointed should be stored
 * in ObservablesHistory.
 * The state of the mdrun machinery that needs to be checkpointed is also
 * stored elsewhere.
 *
 * \author Berk Hess
 *
 * \inlibraryapi
 * \ingroup module_mdtypes
 */

#ifndef GMX_MDTYPES_STATEBASECLASS_H
#define GMX_MDTYPES_STATEBASECLASS_H

#include <array>
#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct df_history_t;
class history_t;
class LocalState;
struct t_inputrec;

namespace gmx
{
struct AwhHistory;
}

/*! \brief The microstate of the system
 *
 * The global state will contain complete data for all used entries.
 * The local state with domain decomposition will have partial entries
 * for which \p stateEntryIsAtomProperty() is true. Some entries that
 * are used in the global state might not be present in the local state.
 * To avoid confusion, two distinct classes for the global and the
 * local state are introduced.
 * \todo Move pure observables history to ObservablesHistory.
 *
 * The base class State is inherited by GlobalState and LocalState and
 * mainly serves the purpose of avoiding code duplication. It is not meant
 * to be used directly; only instances of the two child classes should
 * be created!
 * As GlobalState and LocalState ARE States and just add a conversion
 * constructor to allow copying data from GlobalState to LocalState and back
 * if needed, inheritance has been considered better design than compositiion
 * in this case.
 */
class State
{
    public:
        //! Constructor
        State();

        //! Copy constructor
        State(const State &state);

        //! Assignment operator
        State & operator= (const State &) = default;

        //! Resizes the T- and P-coupling state variables
        void init_gtc_state(int ngtcExtern, int nnhpresExtern, int nhchainlengthExtern);
        //! Change the number of atoms represented by this state, allocating memory as needed.
        void state_change_natoms(int natomsExtern);

        // All things public
        int                        natoms;         //!< Number of atoms, local + non-local; this is the size of \p x, \p v and \p cg_p, when used
        int                        ngtc;           //!< The number of temperature coupling groups
        int                        nnhpres;        //!< The NH-chain length for the MTTK barostat
        int                        nhchainlength;  //!< The NH-chain length for temperature coupling
        int                        flags;          //!< Set of bit-flags telling which entries are present, see enum at the top of the file
        int                        fep_state;      //!< indicates which of the alchemical states we are in
        std::array<real, efptNR>   lambda;         //!< Free-energy lambda vector
        matrix                     box;            //!< Matrix of box vectors
        matrix                     box_rel;        //!< Relative box vectors to preserve box shape
        matrix                     boxv;           //!< Box velocities for Parrinello-Rahman P-coupling
        matrix                     pres_prev;      //!< Pressure of the previous step for pcoupl
        matrix                     svir_prev;      //!< Shake virial for previous step for pcoupl
        matrix                     fvir_prev;      //!< Force virial of the previous step for pcoupl
        std::vector<double>        nosehoover_xi;  //!< Nose-Hoover coordinates (ngtc)
        std::vector<double>        nosehoover_vxi; //!< Nose-Hoover velocities (ngtc)
        std::vector<double>        nhpres_xi;      //!< Pressure Nose-Hoover coordinates
        std::vector<double>        nhpres_vxi;     //!< Pressure Nose-Hoover velocities
        std::vector<double>        therm_integral; //!< Work exterted N-H/V-rescale T-coupling (ngtc)
        double                     baros_integral; //!< For Berendsen P-coupling conserved quantity
        real                       veta;           //!< Trotter based isotropic P-coupling
        real                       vol0;           //!< Initial volume,required for computing MTTK conserved quantity
        gmx::HostVector<gmx::RVec> x;              //!< The coordinates (natoms)
        PaddedRVecVector           v;              //!< The velocities (natoms)
        PaddedRVecVector           cg_p;           //!< p vector for conjugate gradient minimization

        ekinstate_t                ekinstate;      //!< The state of the kinetic energy

        /* History for special algorithms, should be moved to a history struct */
        history_t                         hist;            //!< Time history for restraints
        df_history_t                     *dfhist;          //!< Free-energy history for free energy analysis
        std::shared_ptr<gmx::AwhHistory>  awhHistory;      //!< Accelerated weight histogram history

        int                               ddp_count;       //!< The DD partitioning count for this state
        int                               ddp_count_cg_gl; //!< The DD partitioning count for index_gl
        std::vector<int>                  cg_gl;           //!< The global cg number of the local cgs

        //! Destructor
        virtual ~State() {}
};
#endif
