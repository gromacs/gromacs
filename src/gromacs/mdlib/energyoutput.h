/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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

#include <cstdint>
#include <cstdio>

#include <array>
#include <memory>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

class energyhistory_t;
struct ener_file;
class gmx_ekindata_t;
struct gmx_enerdata_t;
struct SimulationGroups;
struct gmx_mtop_t;
struct gmx_output_env_t;
struct pull_t;
struct t_ebin;
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
struct MDModulesNotifiers;
enum class StartingBehavior;
} // namespace gmx

extern const char* const pvEnergyFieldName;

extern const char* const enthalpyEnergyFieldName;

extern const std::array<const char*, 9> virialEnergyFieldNames;

//! \brief Printed names for intergroup energies
const char* enumValueToString(NonBondedEnergyTerms enumValue);

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
class EnergyDriftTracker;

/*! \internal
 * \brief Arrays connected to Pressure and Temperature coupling
 */
struct PTCouplingArrays
{
    //! Box velocities for Parrinello-Rahman P-coupling.
    const rvec* boxv;
    //! Nose-Hoover coordinates.
    ArrayRef<const double> nosehoover_xi;
    //! Nose-Hoover velocities.
    ArrayRef<const double> nosehoover_vxi;
    //! Pressure Nose-Hoover coordinates.
    ArrayRef<const double> nhpres_xi;
    //! Pressure Nose-Hoover velocities.
    ArrayRef<const double> nhpres_vxi;
};

/* Energy output class
 *
 * This is the collection of energy averages collected during mdrun, and to
 * be written out to the .edr file.
 *
 * \todo Use more std containers.
 * \todo Write free-energy output also to energy file (after adding more tests)
 */
class EnergyOutput
{
public:
    /*! \brief Initiate MD energy bin
     *
     * \param[in] fp_ene     Energy output file.
     * \param[in] mtop       Topology.
     * \param[in] inputrec   Input parameters.
     * \param[in] pull_work  Pulling simulations data
     * \param[in] fp_dhdl    FEP file.
     * \param[in] isRerun    Is this is a rerun instead of the simulations.
     * \param[in] startingBehavior  Run starting behavior.
     * \param[in] simulationsShareState  Tells whether the physical state is shared over simulations
     * \param[in] mdModulesNotifiers Notifications to MD modules.
     */
    EnergyOutput(ener_file*                fp_ene,
                 const gmx_mtop_t&         mtop,
                 const t_inputrec&         inputrec,
                 const pull_t*             pull_work,
                 FILE*                     fp_dhdl,
                 bool                      isRerun,
                 StartingBehavior          startingBehavior,
                 bool                      simulationsShareState,
                 const MDModulesNotifiers& mdModulesNotifiers);

    ~EnergyOutput();

    /*! \brief Update the averaging structures.
     *
     * Called every step on which the thermodynamic values are evaluated.
     *
     * \param[in] bDoDHDL           Whether the FEP is enabled.
     * \param[in] bSum              If this stepshould be recorded to compute sums and averages.
     * \param[in] time              Current simulation time.
     * \param[in] tmass             Total mass
     * \param[in] enerd             Energy data object.
     * \param[in] fep               FEP data.
     * \param[in] lastbox           PBC data.
     * \param[in] ptCouplingArrays  Arrays connected to pressure and temperature coupling.
     * \param[in] fep_state         The current alchemical state we are in.
     * \param[in] vir               Total virial.
     * \param[in] pres              Pressure.
     * \param[in] ekind             Kinetic energy data.
     * \param[in] mu_tot            Total dipole.
     * \param[in] constr            Constraints object to get RMSD from (for LINCS only).
     */
    void addDataAtEnergyStep(bool                    bDoDHDL,
                             bool                    bSum,
                             double                  time,
                             real                    tmass,
                             const gmx_enerdata_t*   enerd,
                             const t_lambda*         fep,
                             const matrix            lastbox,
                             PTCouplingArrays        ptCouplingArrays,
                             int                     fep_state,
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
     * \param[in] ekind   Kinetic energy data, including current reference temperatures
     */
    static void printAnnealingTemperatures(FILE*                   log,
                                           const SimulationGroups& groups,
                                           const t_grpopts&        opts,
                                           const gmx_ekindata_t&   ekind);

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
    static void printHeader(FILE* log, int64_t steps, double time);

    /*! \brief Print conserved energy drift message to \p fplog
     *
     * Note that this is only over the current run (times are printed),
     * this is not from the original start time for runs with continuation.
     * This has the advantage that one can find if conservation issues are
     * from the current run with the current settings on the current hardware.
     */
    void printEnergyConservation(FILE* fplog, int simulationPart, bool usingMdIntegrator) const;

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
    TemperatureCoupling etc_ = TemperatureCoupling::No;

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
    PressureCoupling epc_ = PressureCoupling::No;
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

    //! Which energy terms from NonBondedEnergyTerms list should be printed in group-to-group block
    gmx::EnumerationArray<NonBondedEnergyTerms, bool> bEInd_;
    //! Number of energy terms to be printed (these, for which bEInd[] == true)
    int nEc_ = 0;
    //! Number of energy output groups
    int nEg_ = 0;
    //! Number of intergroup energy sets to be printed for each energy term (nE = (nEg*(nEg+1))/2)
    int nE_ = 0;
    //! Indexes for integroup energy sets (each set with nEc energies)
    std::vector<int> igrp_;

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

    //! Array to accumulate values during update
    std::vector<real> tmp_r_;

    //! The dhdl.xvg output file
    FILE* fp_dhdl_ = nullptr;
    //! Whether the free-energy lambda moves dynamically between lambda states
    bool haveFepLambdaMoves_;
    //! Energy components for dhdl.xvg output
    std::vector<double> dE_;
    //! The delta U components (raw data + histogram)
    std::unique_ptr<t_mde_delta_h_coll> dhc_;
    //! Temperatures for simulated tempering groups
    std::vector<real> temperatures_;

    //! For tracking the conserved or total energy
    std::unique_ptr<EnergyDriftTracker> conservedEnergyTracker_;
};

} // namespace gmx

//! Open the dhdl file for output
FILE* open_dhdl(const char* filename, const t_inputrec* ir, const gmx_output_env_t* oenv);

#endif
