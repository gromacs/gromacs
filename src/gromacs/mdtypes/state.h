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

#include <cstdio>

#include <array>
#include <limits>
#include <memory>
#include <vector>

#include "gromacs/gpu_utils/hostallocator.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class gmx_ekindata_t;
struct t_inputrec;
struct t_lambda;
struct PressureCouplingOptions;
enum class FreeEnergyPerturbationType;

namespace gmx
{
struct AwhHistory;
enum class CheckpointDataOperation;
template<CheckpointDataOperation operation>
class CheckpointData;
} // namespace gmx

//! Convenience alias for until all is moved in the gmx namespace
template<class T>
using PaddedHostVector = gmx::PaddedHostVector<T>;

/*
 * The t_state struct should contain all the (possibly) non-static
 * information required to define the state of the system.
 * Currently the random seeds for SD and BD are missing.
 */

/*! \brief Enum for all entries in \p t_state
 *
 * These enums are used in flags as (1<<est...).
 * The order of these enums should not be changed,
 * since that affects the checkpoint (.cpt) file format.
 */
enum class StateEntry : int
{
    Lambda,
    Box,
    BoxRel,
    BoxV,
    PressurePrevious,
    Nhxi,
    ThermInt,
    X,
    V,
    SDxNotSupported,
    Cgp,
    LDRngNotSupported,
    LDRngINotSupported,
    DisreInitF,
    DisreRm3Tav,
    OrireInitF,
    OrireDtav,
    SVirPrev,
    Nhvxi,
    Veta,
    Vol0,
    Nhpresxi,
    Nhpresvxi,
    FVirPrev,
    FepState,
    MCRngNotSupported,
    MCRngINotSupported,
    BarosInt,
    PullComPrevStep,
    Count
};

//! \brief The names of the state entries, defined in src/gromacs/fileio/checkpoint.cpp
const char* enumValueToString(StateEntry enumValue);
/*! \brief Convert enum to bitmask value.
 *
 * Used for setting flags in checkpoint header and verifying which flags are set.
 */
template<typename Enum>
inline int enumValueToBitMask(Enum enumValue)
{
    static_assert(static_cast<int>(Enum::Count) <= std::numeric_limits<int>::digits);
    return 1 << static_cast<int>(enumValue);
}

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
    history_t();

    real              disre_initf;  //!< The scaling factor for initializing the time av.
    std::vector<real> disre_rm3tav; //!< The r^-3 time averaged pair distances
    real              orire_initf;  //!< The scaling factor for initializing the time av.
    std::vector<real> orire_Dtav;   //!< The time averaged orientation tensors
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

    bool                bUpToDate;      //!< Test if all data is up to date
    int                 ekin_n;         //!< The number of tensors
    tensor*             ekinh;          //!< Half step Ekin, size \p ekin_n
    tensor*             ekinf;          //!< Full step Ekin, size \p ekin_n
    tensor*             ekinh_old;      //!< Half step Ekin of the previous step, size \p ekin_n
    tensor              ekin_total;     //!< Total kinetic energy
    std::vector<double> ekinscalef_nhc; //!< Nose-Hoover Ekin scaling factors for full step Ekin
    std::vector<double> ekinscaleh_nhc; //!< Nose-Hoover Ekin scaling factors for half step Ekin
    std::vector<double> vscale_nhc;     //!< Nose-Hoover velocity scaling factors
    real                dekindl;        //!< dEkin/dlambda, with free-energy
    real                mvcos; //!< Cosine(z) component of the momentum, for viscosity calculations
    /*! \brief Whether KE terms have been read from the checkpoint.
     *
     * Only used for managing whether the call to compute_globals
     * before we enter the MD loop should compute these quantities
     * fresh, or not. */
    bool hasReadEkinState;

    /*!
     * \brief Allows to read and write checkpoint within modular simulator
     * \tparam operation  Whether we're reading or writing
     * \param checkpointData  The CheckpointData object
     */
    template<gmx::CheckpointDataOperation operation>
    void doCheckpoint(gmx::CheckpointData<operation> checkpointData);
};

/*! \brief Free-energy sampling history struct
 *
 * \todo Split out into microstate and observables history.
 */
struct df_history_t
{
    int nlambda; //!< total number of lambda states - for history

    bool  bEquil;   //!< Have we reached equilibration
    int*  n_at_lam; //!< number of points observed at each lambda
    real* wl_histo; //!< histogram for WL flatness determination
    real  wl_delta; //!< current wang-landau delta

    real* sum_weights; //!< weights of the states
    real* sum_dg; //!< free energies of the states -- not actually used for weighting, but informational
    real* sum_minvar;   //!< corrections to weights for minimum variance
    real* sum_variance; //!< variances of the states

    real** accum_p;  //!< accumulated bennett weights for n+1
    real** accum_m;  //!< accumulated bennett weights for n-1
    real** accum_p2; //!< accumulated squared bennett weights for n+1
    real** accum_m2; //!< accumulated squared bennett weights for n-1

    real** Tij;           //!< transition matrix
    real** Tij_empirical; //!< Empirical transition matrix

    /*! \brief Allows to read and write checkpoint within modular simulator
     *
     * \tparam operation  Whether we're reading or writing
     * \param checkpointData  The CheckpointData object
     * \param elamstats  How the lambda weights are calculated
     */
    template<gmx::CheckpointDataOperation operation>
    void doCheckpoint(gmx::CheckpointData<operation> checkpointData, LambdaWeightCalculation elamstats);
};


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
    t_state();

    //! Returns the number of atoms represented by this state.
    int numAtoms() const { return numAtoms_; }

    //! Change the number of atoms represented by this state, allocating memory as needed.
    void changeNumAtoms(int numAtoms);

    //! Returns whether entry \p entry is present
    bool hasEntry(StateEntry entry) const { return (flags_ & enumValueToBitMask(entry)) != 0; }

    //! Set entry \p entry to present, resizes corresponding vector to numAtoms() when relevant
    void addEntry(StateEntry entry);

    //! Return an integer with bits set for entries that are present
    int flags() const { return flags_; }

    //! Sets the present entries to the ones set in \p flags
    void setFlags(int flags);

private:
    int numAtoms_; //!< Number of atoms, local + non-local; this is the size of \p x, \p v and \p cg_p, when used
    int flags_; //!< Set of bit-flags telling which entries are present, see enum at the top of the file

    // The rest is still public
public:
    int ngtc;          //!< The number of temperature coupling groups
    int nnhpres;       //!< The number of NH-chains for the MTTK barostat (always 1 or 0)
    int nhchainlength; //!< The NH-chain length for temperature coupling and MTTK barostat
    int fep_state;     //!< indicates which of the alchemical states we are in
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real> lambda; //!< Free-energy lambda vector
    matrix                                                          box; //!< Matrix of box vectors
    //! Relative box vectors characteristic of the box shape, used to to preserve that box shape
    matrix              box_rel;
    matrix              boxv;           //!< Box velocities for Parrinello-Rahman P-coupling
    matrix              pres_prev;      //!< Pressure of the previous step for pcoupl
    matrix              svir_prev;      //!< Shake virial for previous step for pcoupl
    matrix              fvir_prev;      //!< Force virial of the previous step for pcoupl
    std::vector<double> nosehoover_xi;  //!< Nose-Hoover coordinates (ngtc)
    std::vector<double> nosehoover_vxi; //!< Nose-Hoover velocities (ngtc)
    std::vector<double> nhpres_xi;      //!< Pressure Nose-Hoover coordinates
    std::vector<double> nhpres_vxi;     //!< Pressure Nose-Hoover velocities
    std::vector<double> therm_integral; //!< Work exterted N-H/V-rescale T-coupling (ngtc)
    double              baros_integral; //!< For Berendsen P-coupling conserved quantity
    real                veta;           //!< Trotter based isotropic P-coupling
    real                vol0; //!< Initial volume,required for computing MTTK conserved quantity
    PaddedHostVector<gmx::RVec> x;    //!< The coordinates (numAtoms_)
    PaddedHostVector<gmx::RVec> v;    //!< The velocities (numAtoms_)
    PaddedHostVector<gmx::RVec> cg_p; //!< p vector for conjugate gradient minimization

    ekinstate_t ekinstate; //!< The state of the kinetic energy

    /* History for special algorithms, should be moved to a history struct */
    history_t                        hist;       //!< Time history for restraints
    df_history_t*                    dfhist;     //!< Free-energy history for free energy analysis
    std::shared_ptr<gmx::AwhHistory> awhHistory; //!< Accelerated weight histogram history

    int              ddp_count;       //!< The DD partitioning count for this state
    int              ddp_count_cg_gl; //!< The DD partitioning count for index_gl
    std::vector<int> cg_gl;           //!< The global cg number of the local cgs

    std::vector<double> pull_com_prev_step; //!< The COM of the previous step of each pull group
};

#ifndef DOXYGEN
/* We don't document the structs below, as they don't belong here.
 * TODO: Move the next two structs out of state.h.
 */

struct t_extmass
{
    std::vector<double> Qinv; /* inverse mass of thermostat -- computed from inputs, but a good place to store */
    std::vector<double> QPinv; /* inverse mass of thermostat for barostat -- computed from inputs, but a good place to store */
    double              Winv; /* Pressure mass inverse -- computed, not input, but a good place to store. Need to make a matrix later */
};

#endif // DOXYGEN

//! Resizes the T- and P-coupling state variables
void init_gtc_state(t_state* state, int ngtc, int nnhpres, int nhchainlength);

//! Allocates memory for free-energy history
void init_dfhist_state(t_state* state, int dfhistNumLambda);

/*! \brief Compares two states, write the differences to stdout */
void comp_state(const t_state* st1, const t_state* st2, bool bRMSD, real ftol, real abstol);

/*! \brief Allocates an rvec pointer and copy the contents of v to it */
rvec* makeRvecArray(gmx::ArrayRef<const gmx::RVec> v, gmx::Index n);

/*! \brief Determine the relative box components
 *
 * Set box_rel e.g. used in mdrun state, used to preserve the box shape
 * \param[in]    ir      Input record
 * \param[inout] state   State
 */
void set_box_rel(const t_inputrec* ir, t_state* state);

/*! \brief Make sure the relative box shape remains the same
 *
 * This function ensures that the relative box dimensions are
 * preserved, which otherwise might diffuse away due to rounding
 * errors in pressure coupling or the deform option.
 *
 * \param[in]    pressureCoupling  The pressure-coupling options
 * \param[in]    deform  The box-deformation tensor
 * \param[in]    box_rel Relative box dimensions
 * \param[inout] box     The corrected actual box dimensions
 */
void preserveBoxShape(const PressureCouplingOptions& pressureCoupling,
                      const tensor                   deform,
                      matrix                         box_rel,
                      matrix                         box);

/*! \brief Returns an arrayRef to the positions in \p state when \p state!=null
 *
 * When \p state=nullptr, returns an empty arrayRef.
 *
 * \note The size returned is the number of atoms, without padding.
 *
 * \param[in] state  The state, can be nullptr
 */
static inline gmx::ArrayRef<const gmx::RVec> positionsFromStatePointer(const t_state* state)
{
    if (state)
    {
        return gmx::makeConstArrayRef(state->x).subArray(0, state->numAtoms());
    }
    else
    {
        return {};
    }
};

/*! \brief Prints the current lambda state to the log file.
 *
 * \param[in] fplog  The log file. If fplog == nullptr there will be no output.
 * \param[in] lambda The array of lambda values.
 * \param[in] isInitialOutput Whether this output is the initial lambda state or not.
 */
void printLambdaStateToLog(FILE* fplog, gmx::ArrayRef<const real> lambda, bool isInitialOutput);


/*! \brief Fills fep_state and lambda if needed
 *
 * If FEP or simulated tempering is in use, fills \p fep_state
 * and \p lambda on the main rank and sets the reference temperatures
 * in \p ekind on all ranks.
 *
 * Reports the initial lambda state to the log file. */
void initialize_lambdas(FILE*                      fplog,
                        FreeEnergyPerturbationType freeEnergyPerturbationType,
                        bool                       haveSimulatedTempering,
                        const t_lambda&            fep,
                        gmx::ArrayRef<const real>  simulatedTemperingTemps,
                        gmx_ekindata_t*            ekind,
                        bool                       isMain,
                        int*                       fep_state,
                        gmx::ArrayRef<real>        lambda);

#endif
