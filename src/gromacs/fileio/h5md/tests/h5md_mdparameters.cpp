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
 * \brief Tests for the description of the MD parameters written to H5MD files.
 *
 * \author Yang Zhang <yang.zhang@scilifelab.se>
 * \ingroup module_fileio
 */
#include "gmxpre.h"

#include "gromacs/fileio/h5md/h5md_mdparameters.h"

#include <gtest/gtest.h>

#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
namespace test
{
namespace
{

//! Describe the MD parameters of \p ir and \p mtop, as for a run unless \p isRerun is set
KeyValueTreeObject buildTree(const t_inputrec& ir, const gmx_mtop_t& mtop, const bool isRerun = false)
{
    KeyValueTreeBuilder       builder;
    KeyValueTreeObjectBuilder root = builder.rootObject();
    buildMdParameterTree(&root, ir, mtop, isRerun);
    return builder.build();
}

//! Describe the MD parameters of \p ir with a topology that holds no groups or constraints
KeyValueTreeObject buildTree(const t_inputrec& ir)
{
    const gmx_mtop_t mtop;
    return buildTree(ir, mtop);
}

//! Return the "shape" of \p parameter
std::vector<int64_t> shapeOf(const KeyValueTreeObject& parameter)
{
    std::vector<int64_t> shape;
    for (const KeyValueTreeValue& dimension : parameter["shape"].asArray().values())
    {
        shape.push_back(dimension.cast<int64_t>());
    }
    return shape;
}

//! Add the groups \p names of \p type to \p mtop, in the order given
void addGroups(gmx_mtop_t* mtop, const SimulationAtomGroupType type, const std::vector<std::string>& names)
{
    for (const std::string& name : names)
    {
        mtop->groups.groups[type].push_back(static_cast<int>(mtop->groups.groupNames.size()));
        mtop->groups.groupNames.push_back(put_symtab(&mtop->symtab, name.c_str()));
    }
}

/*! \brief Give \p ir one heat bath per reference temperature, with the given coupling times.
 *
 * done_inputrec() frees anneal_time[i] and anneal_temp[i] of every heat bath, so these
 * arrays must exist too.
 */
void setHeatBaths(t_inputrec*              ir,
                  const std::vector<real>& referenceTemperatures,
                  const std::vector<real>& couplingTimes)
{
    GMX_RELEASE_ASSERT(couplingTimes.size() == referenceTemperatures.size(),
                       "Each heat bath needs a coupling time");
    const int numBaths = static_cast<int>(referenceTemperatures.size());
    ir->opts.ngtc      = numBaths;
    snew(ir->opts.ref_t, numBaths);
    snew(ir->opts.tau_t, numBaths);
    snew(ir->opts.anneal_time, numBaths);
    snew(ir->opts.anneal_temp, numBaths);
    for (int i = 0; i < numBaths; i++)
    {
        ir->opts.ref_t[i] = referenceTemperatures[i];
        ir->opts.tau_t[i] = couplingTimes[i];
    }
}

//! Give \p ir two heat baths, at 300 K and 310 K with tau-t 1 and 2 ps, named in \p mtop
void setTwoHeatBaths(t_inputrec* ir, gmx_mtop_t* mtop)
{
    setHeatBaths(ir, { 300, 310 }, { 1, 2 });
    addGroups(mtop, SimulationAtomGroupType::TemperatureCoupling, { "Protein", "non-Protein" });
}


TEST(H5mdMdParameters, RunControl)
{
    t_inputrec inputRecord = {};
    inputRecord.eI         = IntegrationAlgorithm::Steep;
    inputRecord.nsteps     = 50000;
    inputRecord.delta_t    = 0.001;
    // Energy minimization subcategory
    inputRecord.em_tol      = 10.0f;
    inputRecord.em_stepsize = 0.01f;
    // Center of mass removal is set but must not be written for energy minimization
    inputRecord.nstcomm = 100;
    // Mass repartitioning subcategory
    inputRecord.massRepartitionFactor = 3.0f;

    KeyValueTreeObject tree = buildTree(inputRecord);

    EXPECT_TRUE(tree.keyExists("run_control"));
    const auto& runControl = tree["run_control"].asObject();

    EXPECT_TRUE(runControl.keyExists("numsteps"));
    EXPECT_EQ(runControl["numsteps"].asObject()["value"].cast<int64_t>(), 50000);

    EXPECT_TRUE(runControl.keyExists("simulation_type"));
    EXPECT_EQ(runControl["simulation_type"].asObject()["value"].cast<std::string>(),
              "energy_minimization");

    EXPECT_TRUE(runControl.keyExists("energy_minimization"));
    const auto& minimization = runControl["energy_minimization"].asObject();
    EXPECT_TRUE(minimization.keyExists("algorithm"));
    EXPECT_EQ(minimization["algorithm"].asObject()["value"].cast<std::string>(), "steepest_descent");
    EXPECT_TRUE(minimization.keyExists("tolerance"));
    EXPECT_NEAR(minimization["tolerance"].asObject()["value"].cast<double>(), 10.0, 1e-6);
    EXPECT_TRUE(minimization.keyExists("step_size"));
    EXPECT_NEAR(minimization["step_size"].asObject()["value"].cast<double>(), 0.01, 1e-6);

    EXPECT_TRUE(runControl.keyExists("mass_repartitioning"));
    const auto& massRepartition = runControl["mass_repartitioning"].asObject();
    EXPECT_TRUE(massRepartition.keyExists("factor"));
    EXPECT_NEAR(massRepartition["factor"].asObject()["value"].cast<double>(), 3.0, 1e-6);

    EXPECT_FALSE(runControl.keyExists("center_of_mass_motion_removal"));
}

TEST(H5mdMdParameters, RunControlComRemovalModes)
{
    const std::vector<std::pair<ComRemovalAlgorithm, std::string>> cases = {
        { ComRemovalAlgorithm::Linear, "translation" },
        { ComRemovalAlgorithm::Angular, "translation_and_rotation" },
        { ComRemovalAlgorithm::LinearAccelerationCorrection,
          "translational_acceleration_correction" }
    };
    for (const auto& [mode, modeName] : cases)
    {
        SCOPED_TRACE(modeName);
        t_inputrec ir           = {};
        ir.eI                   = IntegrationAlgorithm::MD;
        ir.comm_mode            = mode;
        ir.nstcomm              = 100;
        KeyValueTreeObject tree = buildTree(ir);

        const auto& runControl = tree["run_control"].asObject();
        ASSERT_TRUE(runControl.keyExists("center_of_mass_motion_removal"));
        const auto& comRemoval = runControl["center_of_mass_motion_removal"].asObject();
        EXPECT_EQ(comRemoval["interval"].asObject()["value"].cast<int64_t>(), 100);
        EXPECT_EQ(comRemoval["mode"].asObject()["value"].cast<std::string>(), modeName);
        // Group information is not created because it is skipped with an empty return
        EXPECT_FALSE(comRemoval.keyExists("groups"));
    }
}

TEST(H5mdMdParameters, RunControlComRemovalGroups)
{
    // grompp names the whole system "rest" when no groups for COM removal are given
    const std::vector<std::vector<std::string>> cases = { { "rest" }, { "Protein", "SOL" } };
    for (const auto& groupNameStrings : cases)
    {
        SCOPED_TRACE(formatString("groups: %s", joinStrings(groupNameStrings, ", ").c_str()));
        gmx_mtop_t mtop;
        addGroups(&mtop, SimulationAtomGroupType::MassCenterVelocityRemoval, groupNameStrings);

        t_inputrec ir           = {};
        ir.eI                   = IntegrationAlgorithm::MD;
        ir.nstcomm              = 100;
        KeyValueTreeObject tree = buildTree(ir, mtop);

        const auto& comRemoval =
                tree["run_control"].asObject()["center_of_mass_motion_removal"].asObject();
        if (groupNameStrings.front() == "rest")
        {
            // The whole system is meant when no groups are listed
            EXPECT_FALSE(comRemoval.keyExists("groups"));
        }
        else
        {
            const auto& values = comRemoval["groups"].asObject()["value"].asArray().values();
            ASSERT_EQ(values.size(), 2U);
            EXPECT_EQ(values[0].cast<std::string>(), "Protein");
            EXPECT_EQ(values[1].cast<std::string>(), "SOL");
        }
    }
}

TEST(H5mdMdParameters, RunControlRerun)
{
    t_inputrec ir            = {};
    ir.eI                    = IntegrationAlgorithm::MD;
    ir.nsteps                = 1000;
    ir.delta_t               = 0.002;
    ir.comm_mode             = ComRemovalAlgorithm::Linear;
    ir.nstcomm               = 100;
    ir.massRepartitionFactor = 3.0;
    KeyValueTreeObject tree  = buildTree(ir, gmx_mtop_t(), true);

    // A rerun recomputes energies of an existing trajectory, whatever the integrator says
    const auto& runControl = tree["run_control"].asObject();
    EXPECT_EQ(runControl["simulation_type"].asObject()["value"].cast<std::string>(), "rerun");
    // Nothing is integrated, so the input record contributes no other run control setting
    EXPECT_EQ(runControl.properties().size(), 1U);
}

TEST(H5mdMdParameters, RunControlMimic)
{
    // grompp puts the QM atoms in one group and the MM atoms in the generated "rest" group
    gmx_mtop_t mtop;
    addGroups(&mtop, SimulationAtomGroupType::QuantumMechanics, { "QM", "rest" });

    t_inputrec ir           = {};
    ir.eI                   = IntegrationAlgorithm::Mimic;
    ir.nsteps               = 1000;
    ir.delta_t              = 0.002;
    KeyValueTreeObject tree = buildTree(ir, mtop);

    const auto& runControl = tree["run_control"].asObject();
    // MiMiC is outside the engine-agnostic parameters, so it is not described
    EXPECT_EQ(runControl["simulation_type"].asObject()["value"].cast<std::string>(), "other");
    EXPECT_FALSE(runControl.keyExists("mimic"));
    // CPMD integrates, so the integration settings of GROMACS are not written
    EXPECT_FALSE(runControl.keyExists("numsteps"));
    EXPECT_FALSE(runControl.keyExists("molecular_dynamics"));
    EXPECT_FALSE(runControl.keyExists("center_of_mass_motion_removal"));
}

TEST(H5mdMdParameters, RunControlSkipsMinimizationForMD)
{
    t_inputrec inputRecord = {};
    inputRecord.eI         = IntegrationAlgorithm::MD;
    inputRecord.nsteps     = 50000;
    inputRecord.delta_t    = 0.002;

    KeyValueTreeObject tree = buildTree(inputRecord);

    EXPECT_TRUE(tree.keyExists("run_control"));
    const auto& runControl = tree["run_control"].asObject();
    EXPECT_TRUE(runControl.keyExists("numsteps"));
    EXPECT_FALSE(runControl.keyExists("energy_minimization"));
    // Masses are only repartitioned with a factor above 1
    EXPECT_FALSE(runControl.keyExists("mass_repartitioning"));

    ASSERT_TRUE(runControl.keyExists("molecular_dynamics"));
    const auto& md = runControl["molecular_dynamics"].asObject();
    EXPECT_EQ(md["integrator"].asObject()["value"].cast<std::string>(), "leap_frog");
    EXPECT_NEAR(md["timestep"].asObject()["value"].cast<double>(), 0.002, 1e-9);
    EXPECT_EQ(md["timestep"].asObject()["unit"].cast<std::string>(), "ps");
}

TEST(H5mdMdParameters, RunControlLangevinDynamics)
{
    t_inputrec ir = {};
    gmx_mtop_t mtop;
    setTwoHeatBaths(&ir, &mtop);
    ir.eI      = IntegrationAlgorithm::SD1;
    ir.delta_t = 0.002;
    ir.ld_seed = 1993;

    KeyValueTreeObject tree = buildTree(ir, mtop);

    const auto& runControl = tree["run_control"].asObject();
    EXPECT_EQ(runControl["simulation_type"].asObject()["value"].cast<std::string>(),
              "stochastic_dynamics");
    EXPECT_FALSE(runControl.keyExists("molecular_dynamics"));
    ASSERT_TRUE(runControl.keyExists("stochastic_dynamics"));
    const auto& sd = runControl["stochastic_dynamics"].asObject();

    EXPECT_EQ(sd["dynamics"].asObject()["value"].cast<std::string>(), "langevin");
    EXPECT_EQ(sd["integrator"].asObject()["value"].cast<std::string>(), "leap_frog");
    EXPECT_EQ(sd["random_seed"].asObject()["value"].cast<int64_t>(), 1993);
    const auto& groups = sd["groups"].asObject()["value"].asArray().values();
    ASSERT_EQ(groups.size(), 2U);
    EXPECT_EQ(groups[1].cast<std::string>(), "non-Protein");
    const auto& refT = sd["reference_temperature"].asObject()["value"].asArray().values();
    ASSERT_EQ(refT.size(), 2U);
    EXPECT_NEAR(refT[1].cast<double>(), 310.0, 1e-6);

    // The friction is a subgroup of the stochastic dynamics, not a temperature coupling
    ASSERT_TRUE(sd.keyExists("langevin"));
    const auto& langevin = sd["langevin"].asObject();
    const auto& tauT     = langevin["coupling_time"].asObject()["value"].asArray().values();
    ASSERT_EQ(tauT.size(), 2U);
    EXPECT_NEAR(tauT[1].cast<double>(), 2.0, 1e-6);
    EXPECT_EQ(langevin["coupling_time"].asObject()["unit"].cast<std::string>(), "ps");
    EXPECT_FALSE(sd.keyExists("brownian"));
}

TEST(H5mdMdParameters, RunControlBrownianDynamics)
{
    for (const real frictionCoefficient : { real(0), real(5) })
    {
        SCOPED_TRACE(formatString("bd-fric = %g", frictionCoefficient));
        t_inputrec ir = {};
        gmx_mtop_t mtop;
        setTwoHeatBaths(&ir, &mtop);
        ir.eI      = IntegrationAlgorithm::BD;
        ir.delta_t = 0.002;
        ir.bd_fric = frictionCoefficient;

        KeyValueTreeObject tree = buildTree(ir, mtop);

        const auto& sd = tree["run_control"].asObject()["stochastic_dynamics"].asObject();
        EXPECT_EQ(sd["dynamics"].asObject()["value"].cast<std::string>(), "brownian");
        EXPECT_EQ(sd["integrator"].asObject()["value"].cast<std::string>(), "euler");
        EXPECT_FALSE(sd.keyExists("langevin"));
        ASSERT_TRUE(sd.keyExists("brownian"));
        const auto& brownian = sd["brownian"].asObject();

        // Exactly one of the two ways to give the friction is written
        if (frictionCoefficient > 0)
        {
            EXPECT_FALSE(brownian.keyExists("coupling_time"));
            EXPECT_NEAR(brownian["friction_coefficient"].asObject()["value"].cast<double>(), 5.0, 1e-6);
            EXPECT_EQ(brownian["friction_coefficient"].asObject()["unit"].cast<std::string>(),
                      "u ps-1");
        }
        else
        {
            EXPECT_FALSE(brownian.keyExists("friction_coefficient"));
            EXPECT_EQ(brownian["coupling_time"].asObject()["value"].asArray().values().size(), 2U);
        }
    }
}

TEST(H5mdMdParameters, RunControlComRemovalOnlyDuringDynamics)
{
    struct Case
    {
        IntegrationAlgorithm integrator;
        ComRemovalAlgorithm  mode;
        bool                 expectComRemoval;
    };
    const std::vector<Case> cases = { { IntegrationAlgorithm::MD, ComRemovalAlgorithm::Linear, true },
                                      { IntegrationAlgorithm::SD1, ComRemovalAlgorithm::Linear, true },
                                      { IntegrationAlgorithm::MD, ComRemovalAlgorithm::No, false },
                                      { IntegrationAlgorithm::Steep, ComRemovalAlgorithm::Linear, false },
                                      { IntegrationAlgorithm::Mimic, ComRemovalAlgorithm::Linear, false } };
    for (const auto& c : cases)
    {
        SCOPED_TRACE(formatString(
                "integrator %s, mode %s", enumValueToString(c.integrator), enumValueToString(c.mode)));
        t_inputrec ir = {};
        gmx_mtop_t mtop;
        setTwoHeatBaths(&ir, &mtop);
        ir.eI                   = c.integrator;
        ir.comm_mode            = c.mode;
        ir.nstcomm              = 100;
        KeyValueTreeObject tree = buildTree(ir, mtop);

        EXPECT_EQ(tree["run_control"].asObject().keyExists("center_of_mass_motion_removal"),
                  c.expectComRemoval);
    }
}

TEST(H5mdMdParameters, OutputControl)
{
    t_inputrec inputRecord                            = {};
    inputRecord.outputControl.nstxout                 = 50000;
    inputRecord.outputControl.nstvout                 = 100000;
    inputRecord.outputControl.nstfout                 = 200;
    inputRecord.outputControl.nstenergy               = 500;
    inputRecord.outputControl.nstlog                  = 1000;
    inputRecord.outputControl.nstcalcenergy           = 100;
    inputRecord.outputControl.nstxout_compressed      = 10000;
    inputRecord.outputControl.x_compression_precision = 1000.0f;

    KeyValueTreeObject tree = buildTree(inputRecord);

    EXPECT_TRUE(tree.keyExists("output_control"));
    const auto& output = tree["output_control"].asObject();

    EXPECT_TRUE(output.keyExists("position"));
    EXPECT_EQ(output["position"].asObject()["value"].cast<int64_t>(), 50000);

    EXPECT_TRUE(output.keyExists("velocity"));
    EXPECT_EQ(output["velocity"].asObject()["value"].cast<int64_t>(), 100000);

    EXPECT_TRUE(output.keyExists("force"));
    EXPECT_EQ(output["force"].asObject()["value"].cast<int64_t>(), 200);

    EXPECT_TRUE(output.keyExists("energy"));
    EXPECT_EQ(output["energy"].asObject()["value"].cast<int64_t>(), 500);

    EXPECT_TRUE(output.keyExists("log"));
    EXPECT_EQ(output["log"].asObject()["value"].cast<int64_t>(), 1000);
}

TEST(H5mdMdParameters, TemperatureCouplingVRescale)
{
    t_inputrec ir = {};
    gmx_mtop_t mtop;
    setTwoHeatBaths(&ir, &mtop);
    ir.etc        = TemperatureCoupling::VRescale;
    ir.nsttcouple = 10;

    KeyValueTreeObject tree = buildTree(ir, mtop);

    ASSERT_TRUE(tree.keyExists("temperature_coupling"));
    const auto& tc = tree["temperature_coupling"].asObject();

    EXPECT_EQ(tc["algorithm"].asObject()["value"].cast<std::string>(),
              "stochastic_velocity_rescaling");
    EXPECT_EQ(tc["interval"].asObject()["value"].cast<int64_t>(), 10);

    ASSERT_TRUE(tc.keyExists("groups"));
    const auto& groups = tc["groups"].asObject()["value"].asArray().values();
    ASSERT_EQ(groups.size(), 2U);
    EXPECT_EQ(groups[0].cast<std::string>(), "Protein");
    EXPECT_EQ(groups[1].cast<std::string>(), "non-Protein");

    ASSERT_TRUE(tc.keyExists("reference_temperature"));
    const auto& refTempObj = tc["reference_temperature"].asObject();
    EXPECT_EQ(refTempObj["unit"].cast<std::string>(), "K");
    const auto& refTempVals = refTempObj["value"].asArray().values();
    ASSERT_EQ(refTempVals.size(), 2U);
    EXPECT_NEAR(refTempVals[0].cast<double>(), 300.0, 1e-6);
    EXPECT_NEAR(refTempVals[1].cast<double>(), 310.0, 1e-6);

    ASSERT_TRUE(tc.keyExists("stochastic_velocity_rescaling"));
    const auto& svr = tc["stochastic_velocity_rescaling"].asObject();
    ASSERT_TRUE(svr.keyExists("coupling_time"));
    EXPECT_EQ(svr["coupling_time"].asObject()["unit"].cast<std::string>(), "ps");
    ASSERT_TRUE(svr.keyExists("doi"));
    EXPECT_EQ(svr["doi"].cast<std::string>(), "10.1063/1.2408420");
    const auto& tauVals = svr["coupling_time"].asObject()["value"].asArray().values();
    ASSERT_EQ(tauVals.size(), 2U);
    EXPECT_NEAR(tauVals[0].cast<double>(), 1.0, 1e-6);
    EXPECT_NEAR(tauVals[1].cast<double>(), 2.0, 1e-6);

    // Other thermostat subcategories must be absent
    EXPECT_FALSE(tc.keyExists("berendsen"));
    EXPECT_FALSE(tc.keyExists("andersen"));
    EXPECT_FALSE(tc.keyExists("nose_hoover"));
}

TEST(H5mdMdParameters, TemperatureCouplingNoseHoover)
{
    t_inputrec ir = {};
    setHeatBaths(&ir, { 300, 300 }, { 0.5, 0.5 });
    ir.etc                = TemperatureCoupling::NoseHoover;
    ir.nsttcouple         = 1;
    ir.opts.nhchainlength = 10;

    KeyValueTreeObject tree = buildTree(ir);

    const auto& tc = tree["temperature_coupling"].asObject();
    EXPECT_EQ(tc["algorithm"].asObject()["value"].cast<std::string>(), "nose_hoover");

    ASSERT_TRUE(tc.keyExists("nose_hoover"));
    const auto& nh = tc["nose_hoover"].asObject();

    ASSERT_TRUE(nh.keyExists("period_time"));
    const auto& periodVals = nh["period_time"].asObject()["value"].asArray().values();
    ASSERT_EQ(periodVals.size(), 2U);
    EXPECT_DOUBLE_EQ(periodVals[0].cast<double>(), 0.5);

    ASSERT_TRUE(nh.keyExists("chain_length"));
    EXPECT_EQ(nh["chain_length"].asObject()["value"].cast<int64_t>(), 10);
    ASSERT_TRUE(nh.keyExists("doi"));
    EXPECT_EQ(nh["doi"].cast<std::string>(), "10.1063/1.447334");

    EXPECT_FALSE(tc.keyExists("stochastic_velocity_rescaling"));
    EXPECT_FALSE(tc.keyExists("berendsen"));
    EXPECT_FALSE(tc.keyExists("andersen"));
}

TEST(H5mdMdParameters, TemperatureCouplingAndersen)
{
    t_inputrec ir = {};
    ir.etc        = TemperatureCoupling::Andersen;
    ir.nsttcouple = 1;
    setHeatBaths(&ir, { 300 }, { 0.5 });

    KeyValueTreeObject tree = buildTree(ir);

    const auto& tc = tree["temperature_coupling"].asObject();
    EXPECT_EQ(tc["algorithm"].asObject()["value"].cast<std::string>(), "andersen");
    ASSERT_TRUE(tc.keyExists("andersen")) << "The algorithm value must name an existing subgroup";
    EXPECT_EQ(tc["andersen"].asObject()["doi"].cast<std::string>(), "10.1063/1.439486");
    EXPECT_FALSE(tc.keyExists("andersen_massive"));
}

TEST(H5mdMdParameters, TemperatureCouplingAndersenMassive)
{
    t_inputrec ir = {};
    ir.etc        = TemperatureCoupling::AndersenMassive;
    ir.nsttcouple = 1;
    setHeatBaths(&ir, { 300 }, { 0.5 });

    KeyValueTreeObject tree = buildTree(ir);

    const auto& tc = tree["temperature_coupling"].asObject();
    EXPECT_EQ(tc["algorithm"].asObject()["value"].cast<std::string>(), "andersen_massive");
    ASSERT_TRUE(tc.keyExists("andersen_massive"))
            << "The algorithm value must name an existing subgroup";
    EXPECT_EQ(tc["andersen_massive"].asObject()["doi"].cast<std::string>(), "10.1021/ct400109a");
    EXPECT_FALSE(tc.keyExists("andersen"))
            << "The plain Andersen subgroup must not be used for the massive variant";
}

TEST(H5mdMdParameters, PressureCouplingCRescale)
{
    t_inputrec ir                               = {};
    ir.pressureCouplingOptions.epc              = PressureCoupling::CRescale;
    ir.pressureCouplingOptions.epct             = PressureCouplingType::Isotropic;
    ir.pressureCouplingOptions.nstpcouple       = 100;
    ir.pressureCouplingOptions.tau_p            = 5.0;
    ir.pressureCouplingOptions.refcoord_scaling = RefCoordScaling::No;
    ir.pressureCouplingOptions.ref_p[XX][XX]    = 1.0;
    ir.pressureCouplingOptions.ref_p[YY][YY]    = 1.0;
    ir.pressureCouplingOptions.ref_p[ZZ][ZZ]    = 1.0;
    ir.pressureCouplingOptions.compress[XX][XX] = 4.5e-5;
    ir.pressureCouplingOptions.compress[YY][YY] = 4.5e-5;
    ir.pressureCouplingOptions.compress[ZZ][ZZ] = 4.5e-5;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("pressure_coupling"));
    const auto& pc = tree["pressure_coupling"].asObject();

    EXPECT_EQ(pc["algorithm"].asObject()["value"].cast<std::string>(), "stochastic_cell_rescaling");
    EXPECT_EQ(pc["interval"].asObject()["value"].cast<int64_t>(), 100);
    // The scaling of position restraint references is GROMACS-specific and not written
    EXPECT_FALSE(pc.keyExists("reference_coordinate_scaling"));

    ASSERT_TRUE(pc.keyExists("reference_pressure"));
    ASSERT_EQ(pc["reference_pressure"].asObject()["value"].asArray().values().size(), 9U);
    EXPECT_EQ(pc["reference_pressure"].asObject()["unit"].cast<std::string>(), "10+5 Pa");
    EXPECT_EQ(shapeOf(pc["reference_pressure"].asObject()), (std::vector<int64_t>{ 3, 3 }));

    ASSERT_TRUE(pc.keyExists("compressibility"));
    ASSERT_EQ(pc["compressibility"].asObject()["value"].asArray().values().size(), 9U);
    EXPECT_EQ(pc["compressibility"].asObject()["unit"].cast<std::string>(), "10-5 Pa-1");
    EXPECT_EQ(shapeOf(pc["compressibility"].asObject()), (std::vector<int64_t>{ 3, 3 }));

    ASSERT_TRUE(pc.keyExists("stochastic_cell_rescaling"));
    const auto& scr = pc["stochastic_cell_rescaling"].asObject();
    ASSERT_TRUE(scr.keyExists("coupling_time"));
    EXPECT_NEAR(scr["coupling_time"].asObject()["value"].cast<double>(), 5.0f, 1e-6);
    EXPECT_EQ(scr["coupling_time"].asObject()["unit"].cast<std::string>(), "ps");

    EXPECT_FALSE(pc.keyExists("parrinello_rahman"));
    EXPECT_FALSE(pc.keyExists("berendsen"));
}

TEST(H5mdMdParameters, PressureCouplingParrinelloRahman)
{
    t_inputrec ir                               = {};
    ir.pressureCouplingOptions.epc              = PressureCoupling::ParrinelloRahman;
    ir.pressureCouplingOptions.epct             = PressureCouplingType::Isotropic;
    ir.pressureCouplingOptions.nstpcouple       = 50;
    ir.pressureCouplingOptions.tau_p            = 2.0;
    ir.pressureCouplingOptions.refcoord_scaling = RefCoordScaling::All;

    KeyValueTreeObject tree = buildTree(ir);

    const auto& pc = tree["pressure_coupling"].asObject();
    EXPECT_EQ(pc["algorithm"].asObject()["value"].cast<std::string>(), "parrinello_rahman");
    EXPECT_FALSE(pc.keyExists("reference_coordinate_scaling"));

    ASSERT_TRUE(pc.keyExists("parrinello_rahman"));
    const auto& pr = pc["parrinello_rahman"].asObject();
    ASSERT_TRUE(pr.keyExists("period_time"));
    EXPECT_NEAR(pr["period_time"].asObject()["value"].cast<double>(), 2.0, 1e-6);

    EXPECT_FALSE(pc.keyExists("stochastic_cell_rescaling"));
    EXPECT_FALSE(pc.keyExists("berendsen"));
    EXPECT_FALSE(pc.keyExists("mttk"));
}

/*! \brief MTTK is recorded like any other barostat.
 *
 * It is deprecated but still runnable, and grompp accepts it, so refusing to describe it
 * would abort a simulation that GROMACS is willing to perform.
 */
TEST(H5mdMdParameters, PressureCouplingMttk)
{
    t_inputrec ir                               = {};
    ir.pressureCouplingOptions.epc              = PressureCoupling::Mttk;
    ir.pressureCouplingOptions.epct             = PressureCouplingType::Isotropic;
    ir.pressureCouplingOptions.nstpcouple       = 10;
    ir.pressureCouplingOptions.tau_p            = 5.0;
    ir.pressureCouplingOptions.refcoord_scaling = RefCoordScaling::No;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("pressure_coupling"));
    const auto& pc = tree["pressure_coupling"].asObject();
    EXPECT_EQ(pc["algorithm"].asObject()["value"].cast<std::string>(), "mttk");

    ASSERT_TRUE(pc.keyExists("mttk"));
    const auto& mttk = pc["mttk"].asObject();
    ASSERT_TRUE(mttk.keyExists("period_time"));
    EXPECT_NEAR(mttk["period_time"].asObject()["value"].cast<double>(), 5.0, 1e-6);
    EXPECT_EQ(mttk["period_time"].asObject()["unit"].cast<std::string>(), "ps");

    // The specification does not predefine this algorithm, so the subgroup carries a reference
    EXPECT_EQ(mttk["doi"].cast<std::string>(), "10.1080/00268979600100761");

    EXPECT_FALSE(pc.keyExists("parrinello_rahman"));
    EXPECT_FALSE(pc.keyExists("stochastic_cell_rescaling"));
    EXPECT_FALSE(pc.keyExists("berendsen"));
}

TEST(H5mdMdParameters, CoulombInteraction)
{
    t_inputrec ir       = {};
    ir.coulombtype      = CoulombInteractionType::Pme;
    ir.coulomb_modifier = InteractionModifiers::PotShift;
    ir.rcoulomb         = 1.2;
    ir.rcoulomb_switch  = 0.0;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("coulomb_interactions"));
    const auto& coulomb = tree["coulomb_interactions"].asObject();

    EXPECT_EQ(coulomb["type"].asObject()["value"].cast<std::string>(), "pme");
    EXPECT_NEAR(coulomb["cutoff_distance"].asObject()["value"].cast<double>(), 1.2, 1e-6);
    EXPECT_EQ(coulomb["cutoff_distance"].asObject()["unit"].cast<std::string>(), "nm");
}

TEST(H5mdMdParameters, CoulombInteractionP3mIsNamedPppm)
{
    t_inputrec ir  = {};
    ir.coulombtype = CoulombInteractionType::P3mAD;

    KeyValueTreeObject tree = buildTree(ir);

    // The engine-agnostic name of the particle-particle particle-mesh method
    EXPECT_EQ(tree["coulomb_interactions"].asObject()["type"].asObject()["value"].cast<std::string>(),
              "pppm");
}

TEST(H5mdMdParameters, VanDerWaalsInteraction)
{
    t_inputrec ir   = {};
    ir.vdwtype      = VanDerWaalsType::Cut;
    ir.vdw_modifier = InteractionModifiers::PotShift;
    ir.rvdw         = 1.0;
    ir.rvdw_switch  = 0.0;
    ir.eDispCorr    = DispersionCorrectionType::EnerPres;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("vdw_interactions"));
    const auto& vdw = tree["vdw_interactions"].asObject();

    EXPECT_EQ(vdw["type"].asObject()["value"].cast<std::string>(), "cutoff");
    EXPECT_NEAR(vdw["cutoff_distance"].asObject()["value"].cast<double>(), 1.0, 1e-6);
    EXPECT_EQ(vdw["cutoff_distance"].asObject()["unit"].cast<std::string>(), "nm");
    EXPECT_EQ(vdw["dispersion_correction"].asObject()["value"].cast<std::string>(),
              "energy_and_pressure");

    ASSERT_TRUE(vdw.keyExists("cutoff"));
    const auto& cutoff = vdw["cutoff"].asObject();
    EXPECT_EQ(cutoff["potential_modifier"].asObject()["value"].cast<std::string>(),
              "potential_shift");
    // A shifted potential has no switching region
    EXPECT_FALSE(cutoff.keyExists("switch_distance"));
}

TEST(H5mdMdParameters, VanDerWaalsInteractionSwitchDistance)
{
    for (const auto modifier : { InteractionModifiers::PotSwitch, InteractionModifiers::ForceSwitch })
    {
        SCOPED_TRACE(enumValueToString(modifier));
        t_inputrec ir   = {};
        ir.vdwtype      = VanDerWaalsType::Cut;
        ir.vdw_modifier = modifier;
        ir.rvdw         = 1.2;
        ir.rvdw_switch  = 1.0;

        KeyValueTreeObject tree = buildTree(ir);

        const auto& cutoff = tree["vdw_interactions"].asObject()["cutoff"].asObject();
        ASSERT_TRUE(cutoff.keyExists("switch_distance"));
        EXPECT_NEAR(cutoff["switch_distance"].asObject()["value"].cast<double>(), 1.0, 1e-6);
        EXPECT_EQ(cutoff["switch_distance"].asObject()["unit"].cast<std::string>(), "nm");
    }
}

TEST(H5mdMdParameters, CoulombInteractionCutoffAndReactionField)
{
    {
        SCOPED_TRACE("Plain cut-off");
        t_inputrec ir       = {};
        ir.coulombtype      = CoulombInteractionType::Cut;
        ir.coulomb_modifier = InteractionModifiers::PotShift;

        KeyValueTreeObject tree    = buildTree(ir);
        const auto&        coulomb = tree["coulomb_interactions"].asObject();
        ASSERT_TRUE(coulomb.keyExists("cutoff"));
        EXPECT_EQ(coulomb["cutoff"].asObject()["potential_modifier"].asObject()["value"].cast<std::string>(),
                  "potential_shift");
        EXPECT_FALSE(coulomb.keyExists("reaction_field"));
    }
    {
        SCOPED_TRACE("Reaction field");
        t_inputrec ir  = {};
        ir.coulombtype = CoulombInteractionType::RF;
        ir.epsilon_rf  = 62.0;

        KeyValueTreeObject tree    = buildTree(ir);
        const auto&        coulomb = tree["coulomb_interactions"].asObject();
        ASSERT_TRUE(coulomb.keyExists("reaction_field"));
        EXPECT_NEAR(
                coulomb["reaction_field"].asObject()["relative_permittivity"].asObject()["value"].cast<double>(),
                62.0,
                1e-6);
        EXPECT_FALSE(coulomb.keyExists("cutoff"));
    }
}

TEST(H5mdMdParameters, PeriodicBoundaryConditionXYZ)
{
    t_inputrec ir = {};
    ir.pbcType    = PbcType::Xyz;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("boundary_conditions"));
    const auto& pbc = tree["boundary_conditions"].asObject();
    ASSERT_TRUE(pbc.keyExists("periodicity"));
    const auto& periodicity = pbc["periodicity"].asObject()["value"].asArray().values();

    ASSERT_EQ(periodicity.size(), 3U);
    EXPECT_EQ(periodicity[0].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[1].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[2].cast<int64_t>(), 1);
}

TEST(H5mdMdParameters, PeriodicBoundaryConditionXY)
{
    t_inputrec ir = {};
    ir.pbcType    = PbcType::XY;

    KeyValueTreeObject tree = buildTree(ir);

    const auto& pbc = tree["boundary_conditions"].asObject();
    ASSERT_TRUE(pbc.keyExists("periodicity"));
    const auto& periodicity = pbc["periodicity"].asObject()["value"].asArray().values();
    ASSERT_EQ(periodicity.size(), 3U);
    EXPECT_EQ(periodicity[0].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[1].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[2].cast<int64_t>(), 0);
}

//! A screw box is periodic in all dimensions, with a rotation when crossing the x boundary
TEST(H5mdMdParameters, PeriodicBoundaryConditionScrew)
{
    t_inputrec ir = {};
    ir.pbcType    = PbcType::Screw;

    KeyValueTreeObject tree = buildTree(ir);

    const auto& pbc = tree["boundary_conditions"].asObject();
    ASSERT_TRUE(pbc.keyExists("periodicity"));
    const auto& periodicity = pbc["periodicity"].asObject()["value"].asArray().values();
    ASSERT_EQ(periodicity.size(), 3U);
    EXPECT_EQ(periodicity[0].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[1].cast<int64_t>(), 1);
    EXPECT_EQ(periodicity[2].cast<int64_t>(), 1);
}

/*! \brief Give \p mtop a single constrained bond.
 *
 * Whether a run has constraints is a property of the topology, not of the input record, so
 * the constraints category needs one to apply at all. The topology is filled in place
 * because gmx_mtop_t is neither copyable nor movable.
 *
 * \param[in,out] mtop  Topology to fill.
 */
void addOneConstraint(gmx_mtop_t* mtop)
{
    gmx_moltype_t moltype;
    moltype.atoms.nr                                       = 2;
    moltype.ilist[InteractionFunction::Constraints].iatoms = { 0, 0, 1 };

    gmx_molblock_t molblock;
    molblock.type = 0;
    molblock.nmol = 1;

    mtop->moltype.push_back(moltype);
    mtop->molblock.push_back(molblock);
    mtop->natoms = 2;
    mtop->finalize();
}

TEST(H5mdMdParameters, ConstraintsLINCS)
{
    t_inputrec ir = {};
    ir.eConstrAlg = ConstraintAlgorithm::Lincs;
    ir.nLincsIter = 2;
    ir.nProjOrder = 6;

    gmx_mtop_t mtop;
    addOneConstraint(&mtop);

    KeyValueTreeObject tree = buildTree(ir, mtop);

    ASSERT_TRUE(tree.keyExists("constraints"));
    const auto& c = tree["constraints"].asObject();

    EXPECT_EQ(c["algorithm"].asObject()["value"].cast<std::string>(), "lincs");

    ASSERT_TRUE(c.keyExists("lincs"));
    const auto& lincs = c["lincs"].asObject();
    EXPECT_EQ(lincs["iterations"].asObject()["value"].cast<int64_t>(), 2);
    EXPECT_EQ(lincs["order"].asObject()["value"].cast<int64_t>(), 6);

    EXPECT_FALSE(c.keyExists("shake"));
}

TEST(H5mdMdParameters, ConstraintsShake)
{
    t_inputrec ir = {};
    ir.eConstrAlg = ConstraintAlgorithm::Shake;
    ir.shake_tol  = 1e-4;

    gmx_mtop_t mtop;
    addOneConstraint(&mtop);

    KeyValueTreeObject tree = buildTree(ir, mtop);

    ASSERT_TRUE(tree.keyExists("constraints"));
    const auto& c = tree["constraints"].asObject();
    EXPECT_EQ(c["algorithm"].asObject()["value"].cast<std::string>(), "shake");

    ASSERT_TRUE(c.keyExists("shake"));
    EXPECT_NEAR(c["shake"].asObject()["tolerance"].asObject()["value"].cast<double>(), 1e-4, 1e-6);

    EXPECT_FALSE(c.keyExists("lincs"));
}

/*! \brief A category that does not apply to the simulation is omitted whole.
 *
 * The specification makes every parameter of these categories mandatory and the categories
 * themselves optional, so they are written complete or not at all. Writing the category
 * with the settings of an algorithm that never acted would misreport the simulation.
 */
TEST(H5mdMdParameters, CategoriesThatDoNotApplyAreOmitted)
{
    t_inputrec ir                  = {};
    ir.etc                         = TemperatureCoupling::No;
    ir.pressureCouplingOptions.epc = PressureCoupling::No;
    ir.eConstrAlg                  = ConstraintAlgorithm::Lincs;

    {
        SCOPED_TRACE("A topology without constrained bonds leaves out all three categories.");
        const gmx_mtop_t         unconstrained;
        const KeyValueTreeObject tree = buildTree(ir, unconstrained);

        EXPECT_FALSE(tree.keyExists("temperature_coupling"))
                << "Reference temperatures of a run without a thermostat never acted.";
        EXPECT_FALSE(tree.keyExists("pressure_coupling"))
                << "Reference pressure of a run without a barostat never acted.";
        EXPECT_FALSE(tree.keyExists("constraints"))
                << "eConstrAlg names a solver even when no bond is constrained.";
    }

    {
        SCOPED_TRACE("Constraints come back once the topology holds a constrained bond.");
        gmx_mtop_t constrained;
        addOneConstraint(&constrained);
        const KeyValueTreeObject tree = buildTree(ir, constrained);

        EXPECT_TRUE(tree.keyExists("constraints"));
        EXPECT_FALSE(tree.keyExists("temperature_coupling"));
        EXPECT_FALSE(tree.keyExists("pressure_coupling"));
    }
}

//! With MiMiC, CPMD integrates, so the thermostat and barostat of the input record never act.
TEST(H5mdMdParameters, CouplingCategoriesAreOmittedForMimic)
{
    t_inputrec ir                  = {};
    ir.etc                         = TemperatureCoupling::VRescale;
    ir.pressureCouplingOptions.epc = PressureCoupling::CRescale;

    {
        SCOPED_TRACE("Molecular dynamics writes both categories.");
        ir.eI                         = IntegrationAlgorithm::MD;
        const KeyValueTreeObject tree = buildTree(ir);
        EXPECT_TRUE(tree.keyExists("temperature_coupling"));
        EXPECT_TRUE(tree.keyExists("pressure_coupling"));
    }

    {
        SCOPED_TRACE("MiMiC writes neither.");
        ir.eI                         = IntegrationAlgorithm::Mimic;
        const KeyValueTreeObject tree = buildTree(ir);
        EXPECT_FALSE(tree.keyExists("temperature_coupling"));
        EXPECT_FALSE(tree.keyExists("pressure_coupling"));
    }
}

/*! \brief A rerun only recomputes energies, so the settings of the integration are omitted.
 *
 * The thermostat, barostat and constraint solver of the input record never act in a rerun,
 * while the interactions and the output still do.
 */
TEST(H5mdMdParameters, RerunOmitsCategoriesOfTheIntegration)
{
    t_inputrec ir                  = {};
    ir.eI                          = IntegrationAlgorithm::MD;
    ir.etc                         = TemperatureCoupling::VRescale;
    ir.pressureCouplingOptions.epc = PressureCoupling::CRescale;
    ir.eConstrAlg                  = ConstraintAlgorithm::Lincs;
    gmx_mtop_t mtop;
    addOneConstraint(&mtop);

    {
        SCOPED_TRACE("A run writes the categories of the integration.");
        const KeyValueTreeObject tree = buildTree(ir, mtop, false);
        EXPECT_TRUE(tree.keyExists("temperature_coupling"));
        EXPECT_TRUE(tree.keyExists("pressure_coupling"));
        EXPECT_TRUE(tree.keyExists("constraints"));
    }

    {
        SCOPED_TRACE("A rerun of the same input omits them.");
        const KeyValueTreeObject tree = buildTree(ir, mtop, true);
        EXPECT_FALSE(tree.keyExists("temperature_coupling"));
        EXPECT_FALSE(tree.keyExists("pressure_coupling"));
        EXPECT_FALSE(tree.keyExists("constraints"));
        EXPECT_TRUE(tree.keyExists("run_control"));
        EXPECT_TRUE(tree.keyExists("output_control"));
        EXPECT_TRUE(tree.keyExists("coulomb_interactions"));
        EXPECT_TRUE(tree.keyExists("vdw_interactions"));
        EXPECT_TRUE(tree.keyExists("boundary_conditions"));
    }
}

/*! \brief A parameter without values is left out, so that every parameter holds a value.
 *
 * A thermostat without any coupling group in the input record, described without a
 * topology, has no reference temperatures, coupling times or group names to write.
 */
TEST(H5mdMdParameters, ParametersWithoutValuesAreOmitted)
{
    t_inputrec ir = {};
    ir.etc        = TemperatureCoupling::VRescale;
    ir.nsttcouple = 10;

    KeyValueTreeObject tree = buildTree(ir);

    ASSERT_TRUE(tree.keyExists("temperature_coupling"));
    const auto& tc = tree["temperature_coupling"].asObject();
    EXPECT_TRUE(tc.keyExists("algorithm"));
    EXPECT_TRUE(tc.keyExists("interval"));
    EXPECT_FALSE(tc.keyExists("groups"));
    EXPECT_FALSE(tc.keyExists("reference_temperature"));
    ASSERT_TRUE(tc.keyExists("stochastic_velocity_rescaling"));
    EXPECT_FALSE(tc["stochastic_velocity_rescaling"].asObject().keyExists("coupling_time"));
}

} // namespace

} // namespace test
} // namespace gmx
