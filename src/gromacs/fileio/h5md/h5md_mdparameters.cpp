/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Defines the description of the MD parameters written to H5MD files.
 *
 * MiMiC counts as molecular dynamics, but CPMD integrates, so the step count, time step,
 * thermostat, barostat and COM motion removal of GROMACS have no effect and are omitted.
 * MiMiC itself is not described, since it is outside the engine-agnostic parameters.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */

#include "gmxpre.h"

#include "h5md_mdparameters.h"

#include <cstdint>

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/vectypes.h"

namespace gmx
{

namespace
{

/*! \brief Add the parameter \p name holding \p value.
 *
 * Integers are written as int64_t and floating-point numbers as double, whatever their
 * type in the input record, so that a reader finds the same type in every build. A
 * floating-point value must be given with a unit, which may be empty; without one it
 * would convert to the integer overload.
 */
void addParameter(KeyValueTreeObjectBuilder* parent, const std::string& name, const std::string& value)
{
    parent->addObject(name).addValue<std::string>("value", value);
}

//! \copydoc addParameter()
void addParameter(KeyValueTreeObjectBuilder* parent, const std::string& name, const int64_t value)
{
    parent->addObject(name).addValue<int64_t>("value", value);
}

//! \brief Add the parameter \p name holding \p value, with \p unit unless it is empty.
void addParameter(KeyValueTreeObjectBuilder* parent,
                  const std::string&         name,
                  const double               value,
                  const std::string&         unit)
{
    KeyValueTreeObjectBuilder parameter = parent->addObject(name);
    parameter.addValue<double>("value", value);
    if (!unit.empty())
    {
        parameter.addValue<std::string>("unit", unit);
    }
}

/*! \brief Add the parameter \p name holding \p values, with \p unit unless it is empty.
 *
 * Nothing is added when there are no values, so that every parameter holds a value.
 */
template<typename T>
void addArrayParameter(KeyValueTreeObjectBuilder* parent,
                       const std::string&         name,
                       const std::vector<T>&      values,
                       const std::string&         unit)
{
    if (values.empty())
    {
        return;
    }
    KeyValueTreeObjectBuilder          parameter = parent->addObject(name);
    KeyValueTreeUniformArrayBuilder<T> array     = parameter.addUniformArray<T>("value");
    for (const T& value : values)
    {
        array.addValue(value);
    }
    if (!unit.empty())
    {
        parameter.addValue<std::string>("unit", unit);
    }
}

//! \brief Add the subgroup \p name for an algorithm described in the publication \p doi.
KeyValueTreeObjectBuilder addAlgorithmSubgroup(KeyValueTreeObjectBuilder* parent,
                                               const std::string&         name,
                                               const std::string&         doi)
{
    KeyValueTreeObjectBuilder subgroup = parent->addObject(name);
    subgroup.addValue<std::string>("doi", doi);
    return subgroup;
}

//! \brief Return \p values as doubles.
std::vector<double> asDoubles(ArrayRef<const real> values)
{
    return std::vector<double>(values.begin(), values.end());
}

//! \brief Add the parameter \p name holding the tensor \p values row by row, with its shape and \p unit.
void addTensorParameter(KeyValueTreeObjectBuilder* parent,
                        const std::string&         name,
                        const tensor&              values,
                        const std::string&         unit)
{
    KeyValueTreeObjectBuilder               parameter = parent->addObject(name);
    KeyValueTreeUniformArrayBuilder<double> array     = parameter.addUniformArray<double>("value");
    for (const real value : gmx::constArrayRefFromArray(&values[0][0], DIM * DIM))
    {
        array.addValue(value);
    }
    parameter.addUniformArray<int64_t>("shape", { DIM, DIM });
    parameter.addValue<std::string>("unit", unit);
}

//! \brief Return the names of the groups of \p type in \p mtop.
std::vector<std::string> groupNames(const gmx_mtop_t& mtop, const SimulationAtomGroupType type)
{
    std::vector<std::string> names;
    for (const int index : mtop.groups.groups[type])
    {
        names.emplace_back(*mtop.groups.groupNames[index]);
    }
    return names;
}

//! \brief Return the name of the simulation type of \p ir.
const char* simulationTypeName(const t_inputrec& ir, const bool isRerun)
{
    // A rerun recomputes the energies of an existing trajectory with the integrator unused
    if (isRerun)
    {
        return "rerun";
    }
    else if (EI_MD(ir.eI) && !EI_MIMIC(ir.eI))
    {
        return "molecular_dynamics";
    }
    else if (EI_RANDOM(ir.eI))
    {
        return "stochastic_dynamics";
    }
    else if (EI_ENERGY_MINIMIZATION(ir.eI))
    {
        return "energy_minimization";
    }
    else
    {
        return "other";
    }
}

//! \brief Return the name of the molecular dynamics \p integrator.
const char* molecularDynamicsIntegratorName(const IntegrationAlgorithm integrator)
{
    if (integrator == IntegrationAlgorithm::MD)
    {
        return "leap_frog";
    }
    else
    {
        GMX_RELEASE_ASSERT(EI_VV(integrator),
                           "Invalid integration algorithm for molecular dynamics");
        return "velocity_verlet";
    }
}

//! \brief Add the stochastic dynamics subgroup of \p runControl.
void addStochasticDynamics(KeyValueTreeObjectBuilder* runControl, const t_inputrec& ir, const gmx_mtop_t& mtop)
{
    const bool isLangevin = (ir.eI == IntegrationAlgorithm::SD1);
    const bool isBrownian = (ir.eI == IntegrationAlgorithm::BD);
    GMX_RELEASE_ASSERT(isLangevin || isBrownian,
                       "Invalid integration algorithm for stochastic dynamics");

    /* The temperature control is part of the stochastic dynamics, and grompp switches off
     * the temperature coupling, so the heat bath is described here instead. */
    KeyValueTreeObjectBuilder stochasticDynamics = runControl->addObject("stochastic_dynamics");
    addParameter(&stochasticDynamics, "dynamics", isLangevin ? "langevin" : "brownian");
    addParameter(&stochasticDynamics, "integrator", isLangevin ? "leap_frog" : "euler");
    addParameter(&stochasticDynamics, "timestep", ir.delta_t, "ps");
    addArrayParameter(&stochasticDynamics,
                      "groups",
                      groupNames(mtop, SimulationAtomGroupType::TemperatureCoupling),
                      "");
    addArrayParameter(&stochasticDynamics,
                      "reference_temperature",
                      asDoubles(constArrayRefFromArray(ir.opts.ref_t, ir.opts.ngtc)),
                      "K");
    addParameter(&stochasticDynamics, "random_seed", ir.ld_seed);

    const std::vector<double> couplingTimes =
            asDoubles(constArrayRefFromArray(ir.opts.tau_t, ir.opts.ngtc));
    if (isLangevin)
    {
        KeyValueTreeObjectBuilder langevin = stochasticDynamics.addObject("langevin");
        addArrayParameter(&langevin, "coupling_time", couplingTimes, "ps");
    }
    else
    {
        // Brownian dynamics takes either one friction coefficient or a friction of mass / tau-t per group
        KeyValueTreeObjectBuilder brownian = stochasticDynamics.addObject("brownian");
        if (ir.bd_fric > 0)
        {
            addParameter(
                    &brownian, "friction_coefficient", static_cast<double>(ir.bd_fric), "u ps-1");
        }
        else if (ir.bd_fric == 0)
        {
            addArrayParameter(&brownian, "coupling_time", couplingTimes, "ps");
        }
    }
}

//! \brief Return the name of the energy minimization algorithm \p integrator.
const char* energyMinimizationAlgorithmName(const IntegrationAlgorithm integrator)
{
    switch (integrator)
    {
        case IntegrationAlgorithm::Steep: return "steepest_descent";
        case IntegrationAlgorithm::CG: return "conjugate_gradient";
        case IntegrationAlgorithm::LBFGS: return "l_bfgs";
        default: GMX_RELEASE_ASSERT(false, "Invalid energy minimization algorithm"); return "";
    }
}

//! \brief Return the name of the center of mass motion removal \p mode.
const char* comRemovalModeName(const ComRemovalAlgorithm mode)
{
    switch (mode)
    {
        case ComRemovalAlgorithm::Linear: return "translation";
        case ComRemovalAlgorithm::Angular: return "translation_and_rotation";
        case ComRemovalAlgorithm::LinearAccelerationCorrection:
            return "translational_acceleration_correction";
        default: GMX_RELEASE_ASSERT(false, "Invalid center of mass removal algorithm"); return "";
    }
}

//! \brief Add the run control category to \p builder.
void addRunControl(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir, const gmx_mtop_t& mtop, const bool isRerun)
{
    KeyValueTreeObjectBuilder runControl = builder->addObject("run_control");
    addParameter(&runControl, "simulation_type", simulationTypeName(ir, isRerun));

    /* A rerun recomputes the energies of the frames of an existing trajectory: nothing is
     * integrated, and the input record need not be that of the run which produced the
     * trajectory. Only the simulation type is required in this category, so it is all
     * that a rerun writes.
     *
     * TODO: Define the rerun as a simulation type of its own, with the parameters that are
     * specific to it, such as the trajectory it reads. */
    if (isRerun)
    {
        return;
    }

    if (!EI_MIMIC(ir.eI))
    {
        addParameter(&runControl, "numsteps", ir.nsteps);
    }

    if (EI_MD(ir.eI) && !EI_MIMIC(ir.eI))
    {
        KeyValueTreeObjectBuilder molecularDynamics = runControl.addObject("molecular_dynamics");
        addParameter(&molecularDynamics, "timestep", ir.delta_t, "ps");
        addParameter(&molecularDynamics, "integrator", molecularDynamicsIntegratorName(ir.eI));
    }

    if (EI_RANDOM(ir.eI))
    {
        addStochasticDynamics(&runControl, ir, mtop);
    }

    if (EI_ENERGY_MINIMIZATION(ir.eI))
    {
        KeyValueTreeObjectBuilder minimization = runControl.addObject("energy_minimization");
        addParameter(&minimization, "algorithm", energyMinimizationAlgorithmName(ir.eI));
        addParameter(&minimization, "tolerance", static_cast<double>(ir.em_tol), "kJ mol-1 nm-1");
        addParameter(&minimization, "step_size", static_cast<double>(ir.em_stepsize), "nm");
    }

    // COM motion is only removed during dynamics integrated by GROMACS, and not with comm-mode = None
    if (EI_DYNAMICS(ir.eI) && !EI_MIMIC(ir.eI) && ir.comm_mode != ComRemovalAlgorithm::No)
    {
        KeyValueTreeObjectBuilder comRemoval =
                runControl.addObject("center_of_mass_motion_removal");
        addParameter(&comRemoval, "interval", static_cast<int64_t>(ir.nstcomm));
        addParameter(&comRemoval, "mode", comRemovalModeName(ir.comm_mode));
        /* grompp generates the group "rest" for the atoms that no named group covers, which
         * is the whole system when no groups are given, as an absent list of groups means. */
        std::vector<std::string> comGroups =
                groupNames(mtop, SimulationAtomGroupType::MassCenterVelocityRemoval);
        comGroups.erase(std::remove(comGroups.begin(), comGroups.end(), "rest"), comGroups.end());
        addArrayParameter(&comRemoval, "groups", comGroups, "");
    }

    if (ir.massRepartitionFactor > 1.0)
    {
        KeyValueTreeObjectBuilder massRepartitioning = runControl.addObject("mass_repartitioning");
        addParameter(
                &massRepartitioning, "factor", static_cast<double>(ir.massRepartitionFactor), "");
    }
}

//! \brief Add the output control category to \p builder.
void addOutputControl(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir)
{
    KeyValueTreeObjectBuilder outputControl = builder->addObject("output_control");
    addParameter(&outputControl, "position", static_cast<int64_t>(ir.outputControl.nstxout));
    addParameter(&outputControl, "velocity", static_cast<int64_t>(ir.outputControl.nstvout));
    addParameter(&outputControl, "force", static_cast<int64_t>(ir.outputControl.nstfout));
    addParameter(&outputControl, "energy", static_cast<int64_t>(ir.outputControl.nstenergy));
    addParameter(&outputControl, "log", static_cast<int64_t>(ir.outputControl.nstlog));
}

//! \brief Return the name of the temperature coupling \p algorithm.
const char* temperatureCouplingName(const TemperatureCoupling algorithm)
{
    switch (algorithm)
    {
        case TemperatureCoupling::VRescale: return "stochastic_velocity_rescaling";
        case TemperatureCoupling::Andersen: return "andersen";
        case TemperatureCoupling::AndersenMassive: return "andersen_massive";
        case TemperatureCoupling::Berendsen:
        case TemperatureCoupling::Yes: return "berendsen";
        case TemperatureCoupling::NoseHoover: return "nose_hoover";
        default: GMX_RELEASE_ASSERT(false, "Invalid temperature coupling algorithm"); return "";
    }
}

//! \brief Return the publication that describes the temperature coupling \p algorithm.
const char* temperatureCouplingDoi(const TemperatureCoupling algorithm)
{
    switch (algorithm)
    {
        case TemperatureCoupling::VRescale: return "10.1063/1.2408420";
        case TemperatureCoupling::Andersen: return "10.1063/1.439486";
        case TemperatureCoupling::AndersenMassive: return "10.1021/ct400109a";
        case TemperatureCoupling::Berendsen:
        case TemperatureCoupling::Yes: return "10.1063/1.448118";
        case TemperatureCoupling::NoseHoover: return "10.1063/1.447334";
        default: GMX_RELEASE_ASSERT(false, "Invalid temperature coupling algorithm"); return "";
    }
}

//! \brief Add the temperature coupling category to \p builder.
void addTemperatureCoupling(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir, const gmx_mtop_t& mtop)
{
    /* The specification makes every parameter of this category mandatory, so a simulation
     * without a thermostat omits the category rather than filling it with settings that
     * did not act. */
    if (ir.etc == TemperatureCoupling::No || EI_MIMIC(ir.eI))
    {
        return;
    }

    KeyValueTreeObjectBuilder temperatureCoupling = builder->addObject("temperature_coupling");
    addParameter(&temperatureCoupling, "algorithm", temperatureCouplingName(ir.etc));
    addParameter(&temperatureCoupling, "interval", static_cast<int64_t>(ir.nsttcouple));
    addArrayParameter(&temperatureCoupling,
                      "groups",
                      groupNames(mtop, SimulationAtomGroupType::TemperatureCoupling),
                      "");
    addArrayParameter(&temperatureCoupling,
                      "reference_temperature",
                      asDoubles(constArrayRefFromArray(ir.opts.ref_t, ir.opts.ngtc)),
                      "K");

    /* The value of "algorithm" names the subgroup that holds the parameters of the
     * algorithm. Nose-Hoover takes tau-t as the period of the oscillation rather than as
     * a coupling time, so it is named accordingly. */
    KeyValueTreeObjectBuilder algorithm = addAlgorithmSubgroup(
            &temperatureCoupling, temperatureCouplingName(ir.etc), temperatureCouplingDoi(ir.etc));
    const std::vector<double> tauT = asDoubles(constArrayRefFromArray(ir.opts.tau_t, ir.opts.ngtc));
    if (ir.etc == TemperatureCoupling::NoseHoover)
    {
        addArrayParameter(&algorithm, "period_time", tauT, "ps");
        addParameter(&algorithm, "chain_length", static_cast<int64_t>(ir.opts.nhchainlength));
    }
    else
    {
        addArrayParameter(&algorithm, "coupling_time", tauT, "ps");
    }
}

//! \brief Return the name of the pressure coupling \p algorithm.
const char* pressureCouplingName(const PressureCoupling algorithm)
{
    switch (algorithm)
    {
        case PressureCoupling::CRescale: return "stochastic_cell_rescaling";
        case PressureCoupling::ParrinelloRahman: return "parrinello_rahman";
        case PressureCoupling::Berendsen: return "berendsen";
        case PressureCoupling::Mttk: return "mttk";
        default: GMX_RELEASE_ASSERT(false, "Invalid pressure coupling algorithm"); return "";
    }
}

//! \brief Return the publication that describes the pressure coupling \p algorithm.
const char* pressureCouplingDoi(const PressureCoupling algorithm)
{
    switch (algorithm)
    {
        case PressureCoupling::CRescale: return "10.1063/5.0020514";
        case PressureCoupling::ParrinelloRahman: return "10.1063/1.328693";
        case PressureCoupling::Berendsen: return "10.1063/1.448118";
        case PressureCoupling::Mttk: return "10.1080/00268979600100761";
        default: GMX_RELEASE_ASSERT(false, "Invalid pressure coupling algorithm"); return "";
    }
}

//! \brief Add the pressure coupling category to \p builder.
void addPressureCoupling(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir)
{
    const PressureCouplingOptions& options = ir.pressureCouplingOptions;
    // As for temperature coupling, a simulation without a barostat omits the category
    if (options.epc == PressureCoupling::No || EI_MIMIC(ir.eI))
    {
        return;
    }

    KeyValueTreeObjectBuilder pressureCoupling = builder->addObject("pressure_coupling");
    addParameter(&pressureCoupling, "algorithm", pressureCouplingName(options.epc));
    addTensorParameter(&pressureCoupling, "reference_pressure", options.ref_p, "10+5 Pa");
    addTensorParameter(&pressureCoupling, "compressibility", options.compress, "10-5 Pa-1");
    addParameter(&pressureCoupling, "interval", static_cast<int64_t>(options.nstpcouple));

    // Parrinello-Rahman and MTTK take tau-p as the period of the oscillation
    KeyValueTreeObjectBuilder algorithm = addAlgorithmSubgroup(
            &pressureCoupling, pressureCouplingName(options.epc), pressureCouplingDoi(options.epc));
    const bool isPeriod = (options.epc == PressureCoupling::ParrinelloRahman
                           || options.epc == PressureCoupling::Mttk);
    addParameter(&algorithm,
                 isPeriod ? "period_time" : "coupling_time",
                 static_cast<double>(options.tau_p),
                 "ps");
}

//! \brief Return the name of the Coulomb interaction \p type.
const char* coulombTypeName(const CoulombInteractionType type)
{
    switch (type)
    {
        case CoulombInteractionType::Cut: return "cutoff";
        case CoulombInteractionType::RF: return "reaction_field";
        case CoulombInteractionType::Pme: return "pme";
        case CoulombInteractionType::Ewald: return "ewald";
        case CoulombInteractionType::P3mAD: return "pppm";
        case CoulombInteractionType::PmeUser: return "pme_user";
        case CoulombInteractionType::PmeSwitch: return "pme_switch";
        case CoulombInteractionType::PmeUserSwitch: return "pme_user_switch";
        case CoulombInteractionType::RFZero: return "reaction_field_zero";
        case CoulombInteractionType::Switch: return "switch";
        case CoulombInteractionType::Shift: return "shift";
        case CoulombInteractionType::User: return "user";
        case CoulombInteractionType::Fmm: return "fmm";
        default: GMX_RELEASE_ASSERT(false, "Invalid Coulomb interaction type"); return "";
    }
}

//! \brief Return the name of the Coulomb potential \p modifier.
const char* coulombModifierName(const InteractionModifiers modifier)
{
    switch (modifier)
    {
        case InteractionModifiers::PotShift:
        case InteractionModifiers::PotShiftVerletUnsupported: return "potential_shift";
        case InteractionModifiers::None:
        case InteractionModifiers::ExactCutoff: return "none";
        case InteractionModifiers::PotSwitch:
        case InteractionModifiers::ForceSwitch:
            GMX_RELEASE_ASSERT(
                    false,
                    "PotSwitch/ForceSwitch coulomb_modifier should have been rejected by grompp");
            return "";
        default: GMX_RELEASE_ASSERT(false, "Invalid Coulomb interaction modifier"); return "";
    }
}

//! \brief Add the Coulomb interactions category to \p builder.
void addCoulombInteractions(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir)
{
    KeyValueTreeObjectBuilder coulomb = builder->addObject("coulomb_interactions");
    addParameter(&coulomb, "type", coulombTypeName(ir.coulombtype));
    addParameter(&coulomb, "cutoff_distance", static_cast<double>(ir.rcoulomb), "nm");
    addParameter(&coulomb, "relative_permittivity", static_cast<double>(ir.epsilon_r), "");

    if (ir.coulombtype == CoulombInteractionType::Cut)
    {
        KeyValueTreeObjectBuilder cutoff = coulomb.addObject("cutoff");
        addParameter(&cutoff, "potential_modifier", coulombModifierName(ir.coulomb_modifier));
    }
    else if (ir.coulombtype == CoulombInteractionType::RF)
    {
        KeyValueTreeObjectBuilder reactionField = coulomb.addObject("reaction_field");
        addParameter(
                &reactionField, "relative_permittivity", static_cast<double>(ir.epsilon_rf), "");
    }
}

//! \brief Return the name of the van der Waals interaction \p type.
const char* vdwTypeName(const VanDerWaalsType type)
{
    switch (type)
    {
        case VanDerWaalsType::Cut: return "cutoff";
        case VanDerWaalsType::Pme: return "pme";
        case VanDerWaalsType::User: return "user";
        case VanDerWaalsType::Switch:
        case VanDerWaalsType::Shift:
            GMX_RELEASE_ASSERT(false,
                               "The vdwtype of Switch/Shift should have been converted to "
                               "cutoff+modifier by grompp stage.");
            return "";
        default: GMX_RELEASE_ASSERT(false, "Invalid Van der Waals interaction type"); return "";
    }
}

//! \brief Return the name of the dispersion \p correction.
const char* dispersionCorrectionName(const DispersionCorrectionType correction)
{
    switch (correction)
    {
        case DispersionCorrectionType::Ener:
        case DispersionCorrectionType::AllEner: return "energy";
        case DispersionCorrectionType::EnerPres:
        case DispersionCorrectionType::AllEnerPres: return "energy_and_pressure";
        default: GMX_RELEASE_ASSERT(false, "Invalid dispersion correction type"); return "";
    }
}

//! \brief Return the name of the van der Waals potential \p modifier.
const char* vdwModifierName(const InteractionModifiers modifier)
{
    switch (modifier)
    {
        case InteractionModifiers::PotShift:
        case InteractionModifiers::PotShiftVerletUnsupported: return "potential_shift";
        case InteractionModifiers::None:
        case InteractionModifiers::ExactCutoff: return "none";
        case InteractionModifiers::PotSwitch: return "potential_switch";
        case InteractionModifiers::ForceSwitch: return "force_switch";
        default: GMX_RELEASE_ASSERT(false, "Invalid Van der Waals interaction modifier"); return "";
    }
}

//! \brief Add the van der Waals interactions category to \p builder.
void addVanDerWaalsInteractions(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir)
{
    KeyValueTreeObjectBuilder vdw = builder->addObject("vdw_interactions");
    addParameter(&vdw, "type", vdwTypeName(ir.vdwtype));
    addParameter(&vdw, "cutoff_distance", static_cast<double>(ir.rvdw), "nm");
    if (ir.eDispCorr != DispersionCorrectionType::No)
    {
        addParameter(&vdw, "dispersion_correction", dispersionCorrectionName(ir.eDispCorr));
    }

    if (ir.vdwtype == VanDerWaalsType::Cut)
    {
        KeyValueTreeObjectBuilder cutoff = vdw.addObject("cutoff");
        addParameter(&cutoff, "potential_modifier", vdwModifierName(ir.vdw_modifier));
        if (ir.vdw_modifier == InteractionModifiers::PotSwitch
            || ir.vdw_modifier == InteractionModifiers::ForceSwitch)
        {
            addParameter(&cutoff, "switch_distance", static_cast<double>(ir.rvdw_switch), "nm");
        }
    }
}

//! \brief Add the boundary conditions category to \p builder.
void addBoundaryConditions(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir)
{
    KeyValueTreeObjectBuilder boundaryConditions = builder->addObject("boundary_conditions");
    // A screw box is periodic in all dimensions; its rotation across x is not described
    const bool periodicInXYZ = (ir.pbcType == PbcType::Xyz || ir.pbcType == PbcType::Screw);
    const auto periodicInXY  = static_cast<int64_t>(periodicInXYZ || ir.pbcType == PbcType::XY);
    const auto periodicInZ   = static_cast<int64_t>(periodicInXYZ);
    addArrayParameter(&boundaryConditions,
                      "periodicity",
                      std::vector<int64_t>{ periodicInXY, periodicInXY, periodicInZ },
                      "");
}

//! \brief Return the name of the constraint \p algorithm.
const char* constraintAlgorithmName(const ConstraintAlgorithm algorithm)
{
    switch (algorithm)
    {
        case ConstraintAlgorithm::Lincs: return "lincs";
        case ConstraintAlgorithm::Shake: return "shake";
        default: GMX_RELEASE_ASSERT(false, "Invalid constraint algorithm"); return "";
    }
}

//! \brief Add the constraints category to \p builder.
void addConstraints(KeyValueTreeObjectBuilder* builder, const t_inputrec& ir, const gmx_mtop_t& mtop)
{
    /* Whether any bonds are constrained is decided when the topology is built, not by the
     * input record: the mdp "constraints" setting is consumed by grompp and never reaches
     * t_inputrec, while eConstrAlg only names the solver to use if constraints exist and
     * defaults to LINCS regardless. LINCS and SHAKE only act on these two interaction
     * types; rigid water is solved by SETTLE and described with the force field. */
    if (gmx_mtop_ftype_count(mtop, InteractionFunction::Constraints) == 0
        && gmx_mtop_ftype_count(mtop, InteractionFunction::ConstraintsNoCoupling) == 0)
    {
        return;
    }

    KeyValueTreeObjectBuilder constraints = builder->addObject("constraints");
    addParameter(&constraints, "algorithm", constraintAlgorithmName(ir.eConstrAlg));
    if (ir.eConstrAlg == ConstraintAlgorithm::Lincs)
    {
        KeyValueTreeObjectBuilder lincs = addAlgorithmSubgroup(
                &constraints,
                "lincs",
                "10.1002/(SICI)1096-987X(199709)18:12<1463::AID-JCC4>3.0.CO;2-H");
        addParameter(&lincs, "iterations", static_cast<int64_t>(ir.nLincsIter));
        addParameter(&lincs, "order", static_cast<int64_t>(ir.nProjOrder));
    }
    else
    {
        KeyValueTreeObjectBuilder shake =
                addAlgorithmSubgroup(&constraints, "shake", "10.1016/0021-9991(77)90098-5");
        addParameter(&shake, "tolerance", static_cast<double>(ir.shake_tol), "");
    }
}

} // namespace

void buildMdParameterTree(KeyValueTreeObjectBuilder* builder,
                          const t_inputrec&          inputRecord,
                          const gmx_mtop_t&          mtop,
                          const bool                 isRerun)
{
    // A rerun integrates nothing, so no thermostat, barostat or constraint solver acts in it
    const bool integrates = !isRerun;

    addRunControl(builder, inputRecord, mtop, isRerun);
    addOutputControl(builder, inputRecord);
    if (integrates)
    {
        addTemperatureCoupling(builder, inputRecord, mtop);
        addPressureCoupling(builder, inputRecord);
    }
    addCoulombInteractions(builder, inputRecord);
    addVanDerWaalsInteractions(builder, inputRecord);
    addBoundaryConditions(builder, inputRecord);
    if (integrates)
    {
        addConstraints(builder, inputRecord, mtop);
    }
}

} // namespace gmx
