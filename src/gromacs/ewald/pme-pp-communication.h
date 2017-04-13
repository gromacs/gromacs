/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \internal \file
 *
 * \brief This file contains declarations and constants necessary for
 * coordinating the communication for the offload of long-ranged PME
 * work to separate MPI rank, for computing energies and forces
 * (Coulomb and LJ).
 *
 * \author Berk Hess <hess@kth.se>
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_ewald
 */

#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/utility/real.h"

/*! \brief MPI Tags used to separate communication of different types of quantities */
enum {
    eCommType_ChargeA, eCommType_ChargeB, eCommType_SQRTC6A, eCommType_SQRTC6B,
    eCommType_SigmaA, eCommType_SigmaB, eCommType_NR, eCommType_COORD,
    eCommType_CNB
};

//@{
/*! \brief Flags used to coordinate PP-PME communication and computation phases
 *
 * Some parts of the code(gmx_pme_send_q, gmx_pme_recv_q_x) assume
 * that the six first flags are exactly in this order.
 */

#define PP_PME_CHARGE         (1<<0)
#define PP_PME_CHARGEB        (1<<1)
#define PP_PME_SQRTC6         (1<<2)
#define PP_PME_SQRTC6B        (1<<3)
#define PP_PME_SIGMA          (1<<4)
#define PP_PME_SIGMAB         (1<<5)
#define PP_PME_COORD          (1<<6)
#define PP_PME_ENER_VIR       (1<<9)
#define PP_PME_FINISH         (1<<10)
#define PP_PME_SWITCHGRID     (1<<11)
#define PP_PME_RESETCOUNTERS  (1<<12)
//@}

/*! \brief Return values for gmx_pme_recv_q_x */
enum {
    pmerecvqxX,            /* calculate PME mesh interactions for new x    */
    pmerecvqxFINISH,       /* the simulation should finish, we should quit */
    pmerecvqxSWITCHGRID,   /* change the PME grid size                     */
    pmerecvqxRESETCOUNTERS /* reset the cycle and flop counters            */
};

/*! \internal
 * \brief Helper struct for PP-PME communication of parameters.
 *
 * The contents are communicated over MPI in memcpy style, so should
 * remain suitable for that.
 */
struct gmx_pme_comm_n_box_t
{
    int             natoms;     /**< Number of atoms */
    matrix          box;        /**< Box */
    int             maxshift_x; /**< Maximum shift in x direction */
    int             maxshift_y; /**< Maximum shift in y direction */
    real            lambda_q;   /**< Free-energy lambda for electrostatics */
    real            lambda_lj;  /**< Free-energy lambda for Lennard-Jones */
    unsigned int    flags;      /**< Control flags */
    gmx_int64_t     step;       /**< MD integration step number */
    //@{
    /*! \brief Used in PME grid tuning */
    ivec            grid_size;
    real            ewaldcoeff_q;
    real            ewaldcoeff_lj;
    //@}
};

/*! \internal
 * \brief Helper struct for PP-PME communication of virial and energy.
 *
 * The contents are communicated over MPI in memcpy style, so should
 * remain suitable for that.
 */
struct gmx_pme_comm_vir_ene_t
{
    //@{
    /*! \brief Virial, energy, and derivative of potential w.r.t. lambda for charge and Lennard-Jones */
    matrix          vir_q;
    matrix          vir_lj;
    real            energy_q;
    real            energy_lj;
    real            dvdlambda_q;
    real            dvdlambda_lj;
    //@}
    float           cycles;     /**< Counter of CPU cycles used */
    gmx_stop_cond_t stop_cond;  /**< Flag used in responding to an external signal to terminate */
};
