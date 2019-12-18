/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
 * \brief Header for the code that writes energy-like quantities.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \author Artem Zhmurov <zhmurov@gmail.com>
 *
 * \ingroup module_mdlib
 * \inlibraryapi
 */
#ifndef GMX_MDLIB_ENERGYOUTPUT_H
#define GMX_MDLIB_ENERGYOUTPUT_H

#include <cstdio>

#include "gromacs/mdtypes/enerdata.h"

class energyhistory_t;
struct ener_file;
struct gmx_ekindata_t;
struct gmx_enerdata_t;
struct SimulationGroups;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_t;
struct t_ebin;
struct t_expanded;
struct t_fcdata;
struct t_grpopts;
struct t_inputrec;
struct t_lambda;
class t_state;

struct t_mde_delta_h_coll;

namespace gmx
{
class Awh;
class Constraints;
struct MdModulesNotifier;
enum class StartingBehavior;
} // namespace gmx

//! \brief Printed names for intergroup energies
extern const char* egrp_nm[egNR + 1];

/* \brief delta_h block type enum: the kinds of energies written out. */
enum
{
    //! Delta H BAR energy difference
    dhbtDH = 0,
    //! dH/dlambda derivative
    dhbtDHDL = 1,
    //! System energy
    dhbtEN,
    //! pV term
    dhbtPV,
    //! Expanded ensemble statistics
    dhbtEXPANDED,
    //! Total number of energy types in this enum
    dhbtNR
};

namespace gmx
{

/* Energy output class
 *
 * This is the collection of energy averages collected during mdrun, and to
 * be written out to the .edr file.
 *
 * \todo Use more std containers.
 * \todo Remove GMX_CONSTRAINTVIR
 * \todo Write free-energy output also to energy file (after adding more tests)
 */
class EnergyOutput
{
public:
    /*! \brief Initiate MD energy bin
     *
     * \param[in] fp_ene     Energy output file.
     * \param[in] mtop       Topology.
     * \param[in] ir         Input parameters.
     * \param[in] pull_work  Pulling simulations data
     * \param[in] fp_dhdl    FEP file.
     * \param[in] isRerun    Is this is a rerun instead of the simulations.
     * \param[in] startingBehavior  Run starting behavior.
     * \param[in] mdModulesNotifier Notifications to MD modules.
     */
    EnergyOutput(ener_file*               fp_ene,
                 const gmx_mtop_t*        mtop,
                 const t_inputrec*        ir,
                 const pull_t*            pull_work,
                 FILE*                    fp_dhdl,
                 bool                     isRerun,
                 StartingBehavior         startingBehavior,
                 const MdModulesNotifier& mdModulesNotifier);

    ~EnergyOutput();

    /*! \brief Update the averaging structures.
     *
     * Called every step on which the thermodynamic values are evaluated.
     *
     * \param[in] bDoDHDL  Whether the FEP is enabled.
     * \param[in] bSum     If this stepshould be recorded to compute sums and averaes.
     * \param[in] time     Current simulation time.
     * \param[in] tmass    Total mass
     * \param[in] enerd    Energy data object.
     * \param[in] state    System state.
     * \param[in] fep      FEP data.
     * \param[in] expand   Expanded ensemble (for FEP).
     * \param[in] lastbox  PBC data.
     * \param[in] svir     Constraint virial.
     * \param[in] fvir     Force virial.
     * \param[in] vir      Total virial.
     * \param[in] pres     Pressure.
     * \param[in] ekind    Kinetic energy data.
     * \param[in] mu_tot   Total dipole.
     * \param[in] constr   Constraints object to get RMSD from (for LINCS only).
     */
    void addDataAtEnergyStep(bool                    bDoDHDL,
                             bool                    bSum,
                             double                  time,
                             real                    tmass,
                             const gmx_enerdata_t*   enerd,
                             const t_state*          state,
                             const t_lambda*         fep,
                             const t_expanded*       expand,
                             const matrix            lastbox,
                             const tensor            svir,
                             const tensor            fvir,
                             const tensor            vir,
                             const tensor            pres,
                             const gmx_ekindata_t*   ekind,
                             const rvec              mu_tot,
                             const gmx::Constraints* constr);

    /*! \brief Update the data averaging structure counts.
     *
     * Updates the number of steps, the values have not being computed.
     */
    void recordNonEnergyStep();

    /*! \brief Writes current quantites to log and energy files.
     *
     * Prints current values of energies, pressure, temperature, restraint
     * data, etc. to energy output file and to the log file (if not nullptr).
     *
     * This function only does something useful when bEne || bDR || bOR || log.
     *
     * \todo Perhaps this responsibility should involve some other
     *       object visiting all the contributing objects.
     *
     * \param[in] fp_ene   Energy file for the output.
     * \param[in] bEne     If it is a step for energy output or last step.
     * \param[in] bDR      If it is a step of writing distance restraints.
     * \param[in] bOR      If it is a step of writing orientation restraints.
     * \param[in] log      Pointer to the log file.
     * \param[in] step     Current step.
     * \param[in] time     Current simulation time.
     * \param[in] fcd      Bonded force computation data,
     *                     including orientation and distance restraints.
     * \param[in] awh      AWH data.
     */
    void printStepToEnergyFile(ener_file* fp_ene,
                               bool       bEne,
                               bool       bDR,
                               bool       bOR,
                               FILE*      log,
                               int64_t    step,
                               double     time,
                               t_fcdata*  fcd,
                               gmx::Awh*  awh);

    /*! \brief Print reference temperatures for annealing groups.
     *
     * Nothing is done if log is nullptr.
     *
     * \param[in] log     Log file to print to.
     * \param[in] groups  Information on atom groups.
     * \param[in] opts    Atom temperature coupling groups options
     *                    (annealing is done by groups).
     */
    void printAnnealingTemperatures(FILE* log, SimulationGroups* groups, t_grpopts* opts);

    /*! \brief Prints average values to log file.
     *
     * This is called at the end of the simulation run to print accumulated average values.
     * Nothing it done if log is nullptr.
     *
     * \param[in]   log      Where to print.
     * \param[in]   groups   Atom groups.
     */
    void printAverages(FILE* log, const SimulationGroups* groups);

    /*! \brief Get the number of thermodynamic terms recorded.
     *
     * \todo Refactor this to return the expected output size,
     *       rather than exposing the implementation details about
     *       thermodynamic terms.
     *
     * \returns Number of thermodynamic terms.
     */
    int numEnergyTerms() const;

    /*! \brief Fill the energyhistory_t data.
     *
     * Between .edr writes, the averages are history dependent,
     * and that history needs to be retained in checkpoints.
     * These functions set/read the energyhistory_t class
     * that is written to checkpoints.
     *
     * \param[out] enerhist  Energy history data structure.
     */
    void fillEnergyHistory(energyhistory_t* enerhist) const;

    /*! \brief Restore from energyhistory_t data.
     *
     * \param[in] enerhist  Energy history data structure.
     */
    void restoreFromEnergyHistory(const energyhistory_t& enerhist);

    //! Print an output header to the log file.
    void printHeader(FILE* log, int64_t steps, double time);

private:
    //! Timestep
    double delta_t_ = 0;

    //! Structure to store energy components and their running averages
    t_ebin* ebin_ = nullptr;

    //! Is the periodic box triclinic
    bool bTricl_ = false;
    //! NHC trotter is used
    bool bNHC_trotter_ = false;
    //! If Nose-Hoover chains should be printed
    bool bPrintNHChains_ = false;
    //! If MTTK integrator was used
    bool bMTTK_ = false;

    //! Temperature control scheme
    int etc_ = 0;

    //! Which of the main energy terms should be printed
    bool bEner_[F_NRE] = { false };
    //! Index for main energy terms
    int ie_ = 0;
    //! Number of energy terms from F_NRE list to be saved (i.e. number of 'true' in bEner)
    int f_nre_ = 0;

    //! Index for constraints RMSD
    int iconrmsd_ = 0;
    /* !\brief How many constraints RMSD terms there are.
     * Usually 1 if LINCS is used and 0 otherwise)
     * nCrmsd > 0 indicates when constraints RMSD is saves, hence no boolean
     */
    int nCrmsd_ = 0;

    //! Is the periodic box dynamic
    bool bDynBox_ = false;
    //! Index for box dimensions
    int ib_ = 0;
    //! Index for box volume
    int ivol_ = 0;
    //! Index for density
    int idens_ = 0;
    //! Triclinic box and not a rerun
    bool bDiagPres_ = false;
    //! Reference pressure, averaged over dimensions
    real ref_p_ = 0.0;
    //! Index for thermodynamic work (pV)
    int ipv_ = 0;
    //! Index for entalpy (pV + total energy)
    int ienthalpy_ = 0;

    /*! \brief If the constraints virial should be printed.
     * Can only be true if "GMX_CONSTRAINTVIR" environmental variable is set */
    bool bConstrVir_ = false;
    //! Index for constrains virial
    int isvir_ = 0;
    //! Index for force virial
    int ifvir_ = 0;

    //! If we have pressure computed
    bool bPres_ = false;
    //! Index for total virial
    int ivir_ = 0;
    //! Index for pressure
    int ipres_ = 0;
    /*! \brief Index for surface tension
     * [(pres[ZZ][ZZ]-(pres[XX][XX]+pres[YY][YY])*0.5)*box[ZZ][ZZ]]*/
    int isurft_ = 0;

    //! Pressure control scheme
    int epc_ = 0;
    //! Index for velocity of the box borders
    int ipc_ = 0;

    //! If dipole was calculated
    bool bMu_ = false;
    //! Index for the dipole
    int imu_ = 0;

    //! Index for coseine acceleration used for viscocity calculation
    int ivcos_ = 0;
    //! Index for viscocity
    int ivisc_ = 0;

    //! Which energy terms from egNR list should be printed in group-to-group block
    bool bEInd_[egNR] = { false };
    //! Number of energy terms to be printed (these, for which bEInd[] == true)
    int nEc_ = 0;
    //! Number of energy output groups
    int nEg_ = 0;
    //! Number of intergroup energy sets to be printed for each energy term (nE = (nEg*(nEg+1))/2)
    int nE_ = 0;
    //! Indexes for integroup energy sets (each set with nEc energies)
    int* igrp_ = nullptr;

    //! Number of temperature coupling groups
    int nTC_ = 0;
    //! Index for temperature
    int itemp_ = 0;

    //! Number of Nose-Hoover chains
    int nNHC_ = 0;
    //! Number of temperature coupling coefficients in case of NH Chains
    int mde_n_ = 0;
    //! Index for temperature coupling coefficient in case of NH chains
    int itc_ = 0;

    //! Number of temperature coupling terms if the temperature control is dealt by barostat (MTTK case)
    int nTCP_ = 0;
    //! Scalling factors for temperaturs control in MTTK
    int mdeb_n_ = 0;
    //! Index for scalling factor of MTTK
    int itcb_ = 0;

    //! Number of acceleration groups
    int nU_ = 0;
    //! Index for group velocities
    int iu_ = 0;

    //! Array to accumulate values during update
    real* tmp_r_ = nullptr;
    //! Array to accumulate values during update
    rvec* tmp_v_ = nullptr;

    //! The dhdl.xvg output file
    FILE* fp_dhdl_ = nullptr;
    //! Energy components for dhdl.xvg output
    double* dE_ = nullptr;
    //! The delta U components (raw data + histogram)
    t_mde_delta_h_coll* dhc_ = nullptr;
    //! Temperatures for simulated tempering groups
    real* temperatures_ = nullptr;
    //! Number of temperatures actually saved
    int numTemperatures_ = 0;
};

} // namespace gmx

//! Open the dhdl file for output
FILE* open_dhdl(const char* filename, const t_inputrec* ir, const gmx_output_env_t* oenv);

#endif
