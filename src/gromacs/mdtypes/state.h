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

#ifndef GMX_MDTYPES_STATE_H
#define GMX_MDTYPES_STATE_H

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

struct t_inputrec;

namespace gmx
{
struct AwhHistory;
}

/*
 * The t_state struct should contain all the (possibly) non-static
 * information required to define the state of the system.
 * Currently the random seeds for SD and BD are missing.
 */

/* \brief Enum for all entries in \p t_state
 *
 * These enums are used in flags as (1<<est...).
 * The order of these enums should not be changed,
 * since that affects the checkpoint (.cpt) file format.
 */
enum {
    estLAMBDA,
    estBOX, estBOX_REL, estBOXV, estPRES_PREV, estNH_XI,  estTHERM_INT,
    estX,   estV,       estSDX_NOTSUPPORTED,  estCGP,
    estLD_RNG_NOTSUPPORTED, estLD_RNGI_NOTSUPPORTED,
    estDISRE_INITF, estDISRE_RM3TAV,
    estORIRE_INITF, estORIRE_DTAV,
    estSVIR_PREV, estNH_VXI, estVETA, estVOL0, estNHPRES_XI, estNHPRES_VXI, estFVIR_PREV,
    estFEPSTATE, estMC_RNG_NOTSUPPORTED, estMC_RNGI_NOTSUPPORTED,
    estBAROS_INT,
    estNR
};

//! \brief The names of the state entries, defined in src/gmxlib/checkpoint.c
extern const char *est_names[estNR];

/*! \libinternal \brief History information for NMR distance and orientation restraints
 *
 * Often this is only used for reporting observables, and thus should not
 * actually be part of the microstate. But with time-dependent restraining
 * they are actually part of the (non-Markovian) microstate.
 * \todo Rename this with a more descriptive name.
 */
class history_t
{
    public:
        //! Constructor
        history_t();

        real  disre_initf;  //!< The scaling factor for initializing the time av.
        int   ndisrepairs;  //!< The number of distance restraints
        real *disre_rm3tav; //!< The r^-3 time averaged pair distances
        real  orire_initf;  //!< The scaling factor for initializing the time av.
        int   norire_Dtav;  //!< The number of matrix element in dtav (npair*5)
        real *orire_Dtav;   //!< The time averaged orientation tensors
};

/*! \libinternal \brief Struct used for checkpointing only
 *
 * This struct would not be required with unlimited precision.
 * But because of limited precision, the COM motion removal implementation
 * can cause the kinetic energy in the MD loop to differ by a few bits from
 * the kinetic energy one would determine from state.v.
 */
class ekinstate_t
{
    public:
        ekinstate_t();

        gmx_bool             bUpToDate;      //!< Test if all data is up to date
        int                  ekin_n;         //!< The number of tensors
        tensor              *ekinh;          //!< Half step Ekin, size \p ekin_n
        tensor              *ekinf;          //!< Full step Ekin, size \p ekin_n
        tensor              *ekinh_old;      //!< Half step Ekin of the previous step, size \p ekin_n
        tensor               ekin_total;     //!< Total kinetic energy
        std::vector<double>  ekinscalef_nhc; //!< Nose-Hoover Ekin scaling factors for full step Ekin
        std::vector<double>  ekinscaleh_nhc; //!< Nose-Hoover Ekin scaling factors for half step Ekin
        std::vector<double>  vscale_nhc;     //!< Nose-Hoover velocity scaling factors
        real                 dekindl;        //!< dEkin/dlambda, with free-energy
        real                 mvcos;          //!< Cosine(z) component of the momentum, for viscosity calculations
};

/*! \brief Free-energy sampling history struct
 *
 * \todo Split out into microstate and observables history.
 */
typedef struct df_history_t
{
    int      nlambda;        //!< total number of lambda states - for history

    gmx_bool bEquil;         //!< Have we reached equilibration
    int     *n_at_lam;       //!< number of points observed at each lambda
    real    *wl_histo;       //!< histogram for WL flatness determination
    real     wl_delta;       //!< current wang-landau delta

    real    *sum_weights;    //!< weights of the states
    real    *sum_dg;         //!< free energies of the states -- not actually used for weighting, but informational
    real    *sum_minvar;     //!< corrections to weights for minimum variance
    real    *sum_variance;   //!< variances of the states

    real   **accum_p;        //!< accumulated bennett weights for n+1
    real   **accum_m;        //!< accumulated bennett weights for n-1
    real   **accum_p2;       //!< accumulated squared bennett weights for n+1
    real   **accum_m2;       //!< accumulated squared bennett weights for n-1

    real   **Tij;            //!< transition matrix
    real   **Tij_empirical;  //!< Empirical transition matrix

} df_history_t;


/*! \brief The microstate of the system
 *
 * The global state will contain complete data for all used entries.
 * The local state with domain decomposition will have partial entries
 * for which \p stateEntryIsAtomProperty() is true. Some entries that
 * are used in the global state might not be present in the local state.
 * \todo Move pure observables history to ObservablesHistory.
 */
class t_state
{
    public:
        //! Constructor
        t_state();

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
};

#ifndef DOXYGEN
/* We don't document the structs below, as they don't belong here.
 * TODO: Move the next two structs out of state.h.
 */

typedef struct t_extmass
{
    double *Qinv;  /* inverse mass of thermostat -- computed from inputs, but a good place to store */
    double *QPinv; /* inverse mass of thermostat for barostat -- computed from inputs, but a good place to store */
    double  Winv;  /* Pressure mass inverse -- computed, not input, but a good place to store. Need to make a matrix later */
    tensor  Winvm; /* inverse pressure mass tensor, computed       */
} t_extmass;


typedef struct
{
    real    veta;
    double  rscale;
    double  vscale;
    double  rvscale;
    double  alpha;
    double *vscale_nhc;
} t_vetavars;

#endif // DOXYGEN

//! Resizes the T- and P-coupling state variables
void init_gtc_state(t_state *state, int ngtc, int nnhpres, int nhchainlength);

//! Change the number of atoms represented by this state, allocating memory as needed.
void state_change_natoms(t_state *state, int natoms);

//! Allocates memory for free-energy history
void init_dfhist_state(t_state *state, int dfhistNumLambda);

/*! \brief Compares two states, write the differences to stdout */
void comp_state(const t_state *st1, const t_state *st2, gmx_bool bRMSD, real ftol, real abstol);

/*! \brief Allocates an rvec pointer and copy the contents of v to it */
rvec *makeRvecArray(gmx::ArrayRef<const gmx::RVec> v,
                    unsigned int                   n);

/*! \brief Determine the relative box components
 *
 * Set box_rel e.g. used in mdrun state, used to preserve the box shape
 * \param[in]    ir      Input record
 * \param[inout] state   State
 */
void set_box_rel(const t_inputrec *ir, t_state *state);

/*! \brief Make sure the relative box shape remains the same
 *
 * This function ensures that the relative box dimensions are
 * preserved, which otherwise might diffuse away due to rounding
 * errors in pressure coupling or the deform option.
 *
 * \param[in]    ir      Input record
 * \param[in]    box_rel Relative box dimensions
 * \param[inout] box     The corrected actual box dimensions
 */
void preserve_box_shape(const t_inputrec *ir, matrix box_rel, matrix box);

#endif
