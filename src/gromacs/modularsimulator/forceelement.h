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
/*! \internal \file
 * \brief Declares the force element for the modular simulator
 *
 * This element calculates the forces, with or without shells or
 * flexible constraints.
 *
 * \author Pascal Merz <pascal.merz@me.com>
 * \ingroup module_modularsimulator
 *
 * This header is only used within the modular simulator module
 */

#ifndef GMX_MODULARSIMULATOR_FORCEELEMENT_H
#define GMX_MODULARSIMULATOR_FORCEELEMENT_H

#include <array>

#include "gromacs/domdec/dlbtiming.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/real.h"

#include "modularsimulatorinterfaces.h"
#include "topologyholder.h"

struct gmx_enfrot;
struct gmx_shellfc_t;
struct gmx_wallcycle;
class CpuPpLongRangeNonbondeds;
struct pull_t;
struct t_nrnb;

namespace gmx
{
class Awh;
class EnergyData;
class FreeEnergyPerturbationData;
class GlobalCommunicationHelper;
class ImdSession;
class LegacySimulatorData;
class MDAtoms;
struct MDModulesNotifiers;
class MdrunScheduleWorkload;
class ModularSimulatorAlgorithmBuilderHelper;
class ObservablesReducer;
class StatePropagatorData;
class VirtualSitesHandler;

/*! \internal
 * \ingroup module_modularsimulator
 * \brief Force element
 *
 * The force element manages the call to either
 * do_force(...) or relax_shell_flexcon(...)
 */
class ForceElement final :
    public ISimulatorElement,
    public ITopologyHolderClient,
    public INeighborSearchSignallerClient,
    public IEnergySignallerClient,
    public IDomDecHelperClient
{
public:
    //! Constructor
    ForceElement(StatePropagatorData*        statePropagatorData,
                 EnergyData*                 energyData,
                 FreeEnergyPerturbationData* freeEnergyPerturbationData,
                 bool                        isVerbose,
                 bool                        isDynamicBox,
                 FILE*                       fplog,
                 const t_commrec*            cr,
                 const t_inputrec*           inputrec,
                 const MDModulesNotifiers&   mdModulesNotifiers,
                 const MDAtoms*              mdAtoms,
                 t_nrnb*                     nrnb,
                 t_forcerec*                 fr,
                 gmx_wallcycle*              wcycle,
                 MdrunScheduleWorkload*      runScheduleWork,
                 VirtualSitesHandler*        vsite,
                 ImdSession*                 imdSession,
                 pull_t*                     pull_work,
                 Constraints*                constr,
                 const gmx_mtop_t&           globalTopology,
                 gmx_enfrot*                 enforcedRotation);
    //! Destructor
    ~ForceElement();

    /*! \brief Register force calculation for step / time
     *
     * \param step                 The step number
     * \param time                 The time
     * \param registerRunFunction  Function allowing to register a run function
     */
    void scheduleTask(Step step, Time time, const RegisterRunFunction& registerRunFunction) override;

    //! Check that we got the local topology
    void elementSetup() override;
    //! Print some final output
    void elementTeardown() override;

    /*! \brief Factory method implementation
     *
     * \param legacySimulatorData  Pointer allowing access to simulator level data
     * \param builderHelper  ModularSimulatorAlgorithmBuilder helper object
     * \param statePropagatorData  Pointer to the \c StatePropagatorData object
     * \param energyData  Pointer to the \c EnergyData object
     * \param freeEnergyPerturbationData  Pointer to the \c FreeEnergyPerturbationData object
     * \param globalCommunicationHelper   Pointer to the \c GlobalCommunicationHelper object
     * \param observablesReducer          Pointer to the \c ObservablesReducer object
     *
     * \return  Pointer to the element to be added. Element needs to have been stored using \c storeElement
     */
    static ISimulatorElement* getElementPointerImpl(LegacySimulatorData* legacySimulatorData,
                                                    ModularSimulatorAlgorithmBuilderHelper* builderHelper,
                                                    StatePropagatorData*        statePropagatorData,
                                                    EnergyData*                 energyData,
                                                    FreeEnergyPerturbationData* freeEnergyPerturbationData,
                                                    GlobalCommunicationHelper* globalCommunicationHelper,
                                                    ObservablesReducer*        observablesReducer);

    //! Callback on domain decomposition repartitioning
    DomDecCallback registerDomDecCallback() override;

private:
    //! ITopologyHolderClient implementation
    void setTopology(const gmx_localtop_t* top) override;
    //! INeighborSearchSignallerClient implementation
    std::optional<SignallerCallback> registerNSCallback() override;
    //! IEnergySignallerClient implementation
    std::optional<SignallerCallback> registerEnergyCallback(EnergySignallerEvent event) override;
    //! The actual do_force call
    template<bool doShellFC>
    void run(Step step, Time time, unsigned int flags);

    //! The shell / FC helper struct
    gmx_shellfc_t* shellfc_;
    //! Whether shells or flexible constraints are present
    const bool doShellFC_;

    //! The next NS step
    Step nextNSStep_;
    //! The next energy calculation step
    Step nextEnergyCalculationStep_;
    //! The next energy calculation step
    Step nextVirialCalculationStep_;
    //! The next free energy calculation step
    Step nextFreeEnergyCalculationStep_;

    // TODO: Clarify relationship to data objects and find a more robust alternative to raw pointers (#3583)
    //! Pointer to the micro state
    StatePropagatorData* statePropagatorData_;
    //! Pointer to the energy data
    EnergyData* energyData_;
    //! Pointer to the free energy perturbation data
    FreeEnergyPerturbationData* freeEnergyPerturbationData_;

    //! The local topology - updated by Topology via Client system
    const gmx_localtop_t* localTopology_;

    //! Whether we're having a dynamic box
    const bool isDynamicBox_;
    //! Whether we're being verbose
    const bool isVerbose_;
    //! The number of shell relaxation steps we did
    Step nShellRelaxationSteps_;

    //! DD / DLB helper object
    const DDBalanceRegionHandler ddBalanceRegionHandler_;
    //! Long range force calculator
    std::unique_ptr<CpuPpLongRangeNonbondeds> longRangeNonbondeds_;

    /* \brief The FEP lambda vector
     *
     * Used if FEP is off, since do_force
     * requires lambda to be allocated anyway
     */
    gmx::EnumerationArray<FreeEnergyPerturbationCouplingType, real> lambda_;

    // Access to ISimulator data
    //! Handles logging.
    FILE* fplog_;
    //! Handles communication.
    const t_commrec* cr_;
    //! Contains user input mdp options.
    const t_inputrec* inputrec_;
    //! Notifiers for MDModules
    const MDModulesNotifiers& mdModulesNotifiers_;
    //! Atom parameters for this domain.
    const MDAtoms* mdAtoms_;
    //! Manages flop accounting.
    t_nrnb* nrnb_;
    //! Manages wall cycle accounting.
    gmx_wallcycle* wcycle_;
    //! Parameters for force calculations.
    t_forcerec* fr_;
    //! Handles virtual sites.
    VirtualSitesHandler* vsite_;
    //! The Interactive Molecular Dynamics session.
    ImdSession* imdSession_;
    //! The pull work object.
    pull_t* pull_work_;
    //! Schedule of work for each MD step for this task.
    MdrunScheduleWorkload* runScheduleWork_;
    //! Handles constraints.
    Constraints* constr_;
    //! Handles enforced rotation.
    gmx_enfrot* enforcedRotation_;
};

} // namespace gmx

#endif // GMX_MODULARSIMULATOR_FORCEELEMENT_H
