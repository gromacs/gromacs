/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \internal \brief
 * Registers command-line modules for pre-5.0 binaries.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "gmxpre.h"

#include "legacymodules.h"

#include <cstdio>

#include "mdrun/mdrun_main.h"
#include "mdrun/nonbonded_bench.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxpreprocess/editconf.h"
#include "gromacs/gmxpreprocess/genconf.h"
#include "gromacs/gmxpreprocess/genion.h"
#include "gromacs/gmxpreprocess/genrestr.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxpreprocess/insert_molecules.h"
#include "gromacs/gmxpreprocess/pdb2gmx.h"
#include "gromacs/gmxpreprocess/solvate.h"
#include "gromacs/gmxpreprocess/x2top.h"
#include "gromacs/tools/check.h"
#include "gromacs/tools/convert_tpr.h"
#include "gromacs/tools/dump.h"
#include "gromacs/tools/eneconv.h"
#include "gromacs/tools/make_ndx.h"
#include "gromacs/tools/mk_angndx.h"
#include "gromacs/tools/pme_error.h"
#include "gromacs/tools/report_methods.h"
#include "gromacs/tools/trjcat.h"
#include "gromacs/tools/trjconv.h"
#include "gromacs/tools/tune_pme.h"

namespace
{

/*! \brief
 * Command line module that provides information about obsolescence.
 *
 * Prints a message directing the user to a wiki page describing replacement
 * options.
 */
class ObsoleteToolModule : public gmx::ICommandLineModule
{
public:
    //! Creates an obsolete tool module for a tool with the given name.
    explicit ObsoleteToolModule(const char* name) : name_(name) {}

    const char* name() const override { return name_; }
    const char* shortDescription() const override { return nullptr; }

    void init(gmx::CommandLineModuleSettings* /*settings*/) override {}
    int  run(int /*argc*/, char* /*argv*/[]) override
    {
        printMessage();
        return 0;
    }
    void writeHelp(const gmx::CommandLineHelpContext& /*context*/) const override
    {
        printMessage();
    }

private:
    static void printMessage()
    {
        std::fprintf(stderr,
                     "This tool is no longer present in GROMACS. Please see\n"
                     "  "
                     "https://manual.gromacs.org/current/user-guide/"
                     "cmdline.html#command-changes-between-versions\n"
                     "for ideas how to perform the same tasks with the "
                     "new tools.\n");
    }

    const char* name_;
};

//! Initializer for a module that defaults to nice level zero.
void initSettingsNoNice(gmx::CommandLineModuleSettings* settings)
{
    settings->setDefaultNiceLevel(0);
}

/*! \brief
 * Convenience function for creating and registering a module.
 *
 * \param[in] manager          Module manager to which to register the module.
 * \param[in] mainFunction     Main function to wrap.
 * \param[in] name             Name for the new module.
 * \param[in] shortDescription One-line description for the new module.
 */
void registerModule(gmx::CommandLineModuleManager*               manager,
                    gmx::CommandLineModuleManager::CMainFunction mainFunction,
                    const char*                                  name,
                    const char*                                  shortDescription)
{
    manager->addModuleCMain(name, shortDescription, mainFunction);
}

/*! \brief
 * Convenience function for creating and registering a module that defaults to
 * -nice 0.
 *
 * \param[in] manager          Module manager to which to register the module.
 * \param[in] mainFunction     Main function to wrap.
 * \param[in] name             Name for the new module.
 * \param[in] shortDescription One-line description for the new module.
 */
void registerModuleNoNice(gmx::CommandLineModuleManager*               manager,
                          gmx::CommandLineModuleManager::CMainFunction mainFunction,
                          const char*                                  name,
                          const char*                                  shortDescription)
{
    manager->addModuleCMainWithSettings(name, shortDescription, mainFunction, &initSettingsNoNice);
}

/*! \brief
 * Convenience function for registering a module for an obsolete tool.
 *
 * \param[in] manager          Module manager to which to register the module.
 * \param[in] name             Name for the obsolete tool.
 */
void registerObsoleteTool(gmx::CommandLineModuleManager* manager, const char* name)
{
    gmx::CommandLineModulePointer module(new ObsoleteToolModule(name));
    manager->addModule(std::move(module));
}

} // namespace

void registerLegacyModules(gmx::CommandLineModuleManager* manager)
{
    registerModule(manager, &gmx_check, "check", "Check and compare files");
    gmx::ICommandLineOptionsModule::registerModuleFactory(
            manager, gmx::DumpInfo::name, gmx::DumpInfo::shortDescription, &gmx::DumpInfo::create);
    registerModule(manager, &gmx_grompp, "grompp", "Make a run input file");
    gmx::ICommandLineOptionsModule::registerModuleFactory(manager,
                                                          gmx::ConvertTprInfo::name,
                                                          gmx::ConvertTprInfo::shortDescription,
                                                          &gmx::ConvertTprInfo::create);
    registerObsoleteTool(manager, "tpbconv");
    registerModule(manager, &gmx_x2top, "x2top", "Generate a primitive topology from coordinates");

    registerModuleNoNice(
            manager,
            &gmx::gmx_mdrun,
            "mdrun",
            "Perform a simulation, do a normal mode analysis or an energy minimization");

    gmx::ICommandLineOptionsModule::registerModuleFactory(manager,
                                                          gmx::NonbondedBenchmarkInfo::name,
                                                          gmx::NonbondedBenchmarkInfo::shortDescription,
                                                          &gmx::NonbondedBenchmarkInfo::create);

    gmx::ICommandLineOptionsModule::registerModuleFactory(manager,
                                                          gmx::InsertMoleculesInfo::name(),
                                                          gmx::InsertMoleculesInfo::shortDescription(),
                                                          &gmx::InsertMoleculesInfo::create);

    gmx::ICommandLineOptionsModule::registerModuleFactory(manager,
                                                          gmx::ReportMethodsInfo::name,
                                                          gmx::ReportMethodsInfo::shortDescription,
                                                          &gmx::ReportMethodsInfo::create);

    gmx::ICommandLineOptionsModule::registerModuleFactory(
            manager, gmx::pdb2gmxInfo::name, gmx::pdb2gmxInfo::shortDescription, &gmx::pdb2gmxInfo::create);

    // Modules from gmx_ana.h.
    registerModule(manager, &gmx_editconf, "editconf", "Convert and manipulates structure files");
    registerModule(manager, &gmx_eneconv, "eneconv", "Convert energy files");
    registerModule(manager, &gmx_solvate, "solvate", "Solvate a system");
    registerObsoleteTool(manager, "genbox");
    registerModule(
            manager, &gmx_genconf, "genconf", "Multiply a conformation in 'random' orientations");
    registerModule(manager,
                   &gmx_genion,
                   "genion",
                   "Generate monoatomic ions on energetically favorable positions");
    registerModule(manager,
                   &gmx_genrestr,
                   "genrestr",
                   "Generate position restraints or distance restraints for index groups");
    registerModule(manager,
                   &gmx_make_edi,
                   "make_edi",
                   "Generate input files for essential dynamics sampling");
    registerModule(manager, &gmx_make_ndx, "make_ndx", "Make index files");
    registerModule(manager, &gmx_mk_angndx, "mk_angndx", "Generate index files for 'gmx angle'");
    registerModule(manager, &gmx_trjcat, "trjcat", "Concatenate trajectory files");
    registerModule(manager, &gmx_trjconv, "trjconv", "Convert and manipulates trajectory files");
    registerModule(manager,
                   &gmx_trjorder,
                   "trjorder",
                   "Order molecules according to their distance to a group");
    registerModule(
            manager, &gmx_xpm2ps, "xpm2ps", "Convert XPM (XPixelMap) matrices to postscript or XPM");

    registerModule(manager, &gmx_anaeig, "anaeig", "Analyze eigenvectors/normal modes");
    registerModule(manager, &gmx_analyze, "analyze", "Analyze data sets");
    registerModule(manager,
                   &gmx_g_angle,
                   "angle",
                   "Calculate distributions and correlations for angles and dihedrals");
    registerModule(
            manager, &gmx_awh, "awh", "Extract data from an accelerated weight histogram (AWH) run");
    registerModule(manager,
                   &gmx_bar,
                   "bar",
                   "Calculate free energy difference estimates through Bennett's acceptance ratio");
    registerObsoleteTool(manager, "bond");
    registerObsoleteTool(manager, "dist");
    registerObsoleteTool(manager, "sas");
    registerObsoleteTool(manager, "sgangle");

    registerModule(manager, &gmx_bundle, "bundle", "Analyze bundles of axes, e.g., helices");
    registerModule(manager,
                   &gmx_chi,
                   "chi",
                   "Calculate everything you want to know about chi and other dihedrals");
    registerModule(manager, &gmx_cluster, "cluster", "Cluster structures");
    registerModule(
            manager, &gmx_clustsize, "clustsize", "Calculate size distributions of atomic clusters");
    registerModule(manager, &gmx_confrms, "confrms", "Fit two structures and calculates the RMSD");
    registerModule(manager, &gmx_covar, "covar", "Calculate and diagonalize the covariance matrix");
    registerModule(manager,
                   &gmx_current,
                   "current",
                   "Calculate dielectric constants and current autocorrelation function");
    registerModule(manager, &gmx_density, "density", "Calculate the density of the system");
    registerModule(
            manager, &gmx_densmap, "densmap", "Calculate 2D planar or axial-radial density maps");
    registerModule(manager, &gmx_densorder, "densorder", "Calculate surface fluctuations");
    registerModule(manager,
                   &gmx_dielectric,
                   "dielectric",
                   "Calculate frequency dependent dielectric constants");
    registerModule(manager, &gmx_dipoles, "dipoles", "Compute the total dipole plus fluctuations");
    registerModule(manager, &gmx_disre, "disre", "Analyze distance restraints");
    registerModule(
            manager, &gmx_dos, "dos", "Analyze density of states and properties based on that");
    registerModule(manager, &gmx_dyecoupl, "dyecoupl", "Extract dye dynamics from trajectories");
    registerModule(manager, &gmx_enemat, "enemat", "Extract an energy matrix from an energy file");
    registerModule(
            manager, &gmx_energy, "energy", "Writes energies to xvg files and display averages");
    registerModule(manager,
                   &gmx_filter,
                   "filter",
                   "Frequency filter trajectories, useful for making smooth movies");
    registerModule(manager, &gmx_gyrate, "gyrate-legacy", "Calculate the radius of gyration");
    registerModule(manager, &gmx_h2order, "h2order", "Compute the orientation of water molecules");
    registerModule(manager, &gmx_hbond, "hbond-legacy", "Compute and analyze hydrogen bonds");
    registerModule(manager, &gmx_helix, "helix", "Calculate basic properties of alpha helices");
    registerModule(manager,
                   &gmx_helixorient,
                   "helixorient",
                   "Calculate local pitch/bending/rotation/orientation inside helices");
    registerModule(manager,
                   &gmx_hydorder,
                   "hydorder",
                   "Compute tetrahedrality parameters around a given atom");
    registerModule(manager, &gmx_lie, "lie", "Estimate free energy from linear combinations");
    registerModule(manager, &gmx_mdmat, "mdmat", "Calculate residue contact maps");
    registerModule(
            manager, &gmx_mindist, "mindist", "Calculate the minimum distance between two groups");
    registerModule(manager, &gmx_nmeig, "nmeig", "Diagonalize the Hessian for normal mode analysis");
    registerModule(manager,
                   &gmx_nmens,
                   "nmens",
                   "Generate an ensemble of structures from the normal modes");
    registerModule(manager,
                   &gmx_nmr,
                   "nmr",
                   "Analyze nuclear magnetic resonance properties from an energy file");
    registerModule(manager,
                   &gmx_nmtraj,
                   "nmtraj",
                   "Generate a virtual oscillating trajectory from an eigenvector");
    registerModule(
            manager, &gmx_order, "order", "Compute the order parameter per atom for carbon tails");
    registerModule(manager,
                   &gmx_pme_error,
                   "pme_error",
                   "Estimate the error of using PME with a given input file");
    registerModule(manager, &gmx_polystat, "polystat", "Calculate static properties of polymers");
    registerModule(manager,
                   &gmx_potential,
                   "potential",
                   "Calculate the electrostatic potential across the box");
    registerModule(manager,
                   &gmx_principal,
                   "principal",
                   "Calculate principal axes of inertia for a group of atoms");
    registerModule(manager, &gmx_rama, "rama", "Compute Ramachandran plots");
    registerModule(manager,
                   &gmx_rms,
                   "rms",
                   "Calculate RMSDs with a reference structure and RMSD matrices");
    registerModule(manager,
                   &gmx_rmsdist,
                   "rmsdist",
                   "Calculate atom pair distances averaged with power -2, -3 or -6");
    registerModule(manager, &gmx_rmsf, "rmsf", "Calculate atomic fluctuations");
    registerModule(manager,
                   &gmx_rotacf,
                   "rotacf",
                   "Calculate the rotational correlation function for molecules");
    registerModule(manager,
                   &gmx_rotmat,
                   "rotmat",
                   "Plot the rotation matrix for fitting to a reference structure");
    registerModule(manager, &gmx_saltbr, "saltbr", "Compute salt bridges");
    registerModule(
            manager, &gmx_sans, "sans-legacy", "Compute small angle neutron scattering spectra");
    registerModule(
            manager, &gmx_saxs, "saxs-legacy", "Compute small angle X-ray scattering spectra");
    registerModule(
            manager, &gmx_sham, "sham", "Compute free energies or other histograms from histograms");
    registerModule(manager,
                   &gmx_sigeps,
                   "sigeps",
                   "Convert c6/12 or c6/cn combinations to and from sigma/epsilon");
    registerModule(manager, &gmx_sorient, "sorient", "Analyze solvent orientation around solutes");
    registerModule(manager, &gmx_spatial, "spatial", "Calculate the spatial distribution function");
    registerModule(manager,
                   &gmx_spol,
                   "spol",
                   "Analyze solvent dipole orientation and polarization around solutes");
    registerModule(manager, &gmx_tcaf, "tcaf", "Calculate viscosities of liquids");
    registerModule(manager,
                   &gmx_traj,
                   "traj",
                   "Plot x, v, f, box, temperature and rotational energy from trajectories");
    registerModule(manager,
                   &gmx_tune_pme,
                   "tune_pme",
                   "Time mdrun as a function of PME ranks to optimize settings");
    registerModule(manager,
                   &gmx_vanhove,
                   "vanhove",
                   "Compute Van Hove displacement and correlation functions");
    registerModule(manager, &gmx_velacc, "velacc", "Calculate velocity autocorrelation functions");
    registerModule(manager,
                   &gmx_wham,
                   "wham",
                   "Perform weighted histogram analysis after umbrella sampling");
    registerModule(manager, &gmx_wheel, "wheel", "Plot helical wheels");

    {
        gmx::CommandLineModuleGroup group =
                manager->addModuleGroup("Generating topologies and coordinates");
        group.addModuleWithDescription("editconf", "Edit the box and write subgroups");
        group.addModule("x2top");
        group.addModule("solvate");
        group.addModule("insert-molecules");
        group.addModule("genconf");
        group.addModule("genion");
        group.addModule("genrestr");
        group.addModule("pdb2gmx");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Running a simulation");
        group.addModule("grompp");
        group.addModule("mdrun");
        group.addModule("convert-tpr");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Viewing trajectories");
        group.addModule("nmtraj");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Processing energies");
        group.addModule("enemat");
        group.addModule("energy");
        group.addModuleWithDescription("mdrun",
                                       "(Re)calculate energies for trajectory frames with -rerun");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Converting files");
        group.addModule("editconf");
        group.addModule("eneconv");
        group.addModule("sigeps");
        group.addModule("trjcat");
        group.addModule("trjconv");
        group.addModule("xpm2ps");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Tools");
        group.addModule("analyze");
        group.addModule("awh");
        group.addModule("filter");
        group.addModule("lie");
        group.addModule("pme_error");
        group.addModule("sham");
        group.addModule("spatial");
        group.addModule("traj");
        group.addModule("tune_pme");
        group.addModule("wham");
        group.addModule("check");
        group.addModule("dump");
        group.addModule("make_ndx");
        group.addModule("mk_angndx");
        group.addModule("trjorder");
        group.addModule("xpm2ps");
        group.addModule("report-methods");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Distances between structures");
        group.addModule("cluster");
        group.addModule("confrms");
        group.addModule("rms");
        group.addModule("rmsf");
    }
    {
        gmx::CommandLineModuleGroup group =
                manager->addModuleGroup("Distances in structures over time");
        group.addModule("mindist");
        group.addModule("mdmat");
        group.addModule("polystat");
        group.addModule("rmsdist");
    }
    {
        gmx::CommandLineModuleGroup group =
                manager->addModuleGroup("Mass distribution properties over time");
        group.addModule("gyrate-legacy");
        group.addModule("polystat");
        group.addModule("rdf");
        group.addModule("rotacf");
        group.addModule("rotmat");
        group.addModule("sans-legacy");
        group.addModule("saxs-legacy");
        group.addModule("traj");
        group.addModule("vanhove");
    }
    {
        gmx::CommandLineModuleGroup group =
                manager->addModuleGroup("Analyzing bonded interactions");
        group.addModule("angle");
        group.addModule("mk_angndx");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Structural properties");
        group.addModule("bundle");
        group.addModule("clustsize");
        group.addModule("disre");
        group.addModule("hbond-legacy");
        group.addModule("order");
        group.addModule("principal");
        group.addModule("rdf");
        group.addModule("saltbr");
        group.addModule("sorient");
        group.addModule("spol");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Kinetic properties");
        group.addModule("bar");
        group.addModule("current");
        group.addModule("dos");
        group.addModule("dyecoupl");
        group.addModule("principal");
        group.addModule("tcaf");
        group.addModule("traj");
        group.addModule("vanhove");
        group.addModule("velacc");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Electrostatic properties");
        group.addModule("current");
        group.addModule("dielectric");
        group.addModule("dipoles");
        group.addModule("potential");
        group.addModule("spol");
        group.addModule("genion");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Protein-specific analysis");
        group.addModule("chi");
        group.addModule("helix");
        group.addModule("helixorient");
        group.addModule("rama");
        group.addModule("wheel");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Interfaces");
        group.addModule("bundle");
        group.addModule("density");
        group.addModule("densmap");
        group.addModule("densorder");
        group.addModule("h2order");
        group.addModule("hydorder");
        group.addModule("order");
        group.addModule("potential");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Covariance analysis");
        group.addModuleWithDescription("anaeig", "Analyze the eigenvectors");
        group.addModule("covar");
        group.addModule("make_edi");
    }
    {
        gmx::CommandLineModuleGroup group = manager->addModuleGroup("Normal modes");
        group.addModuleWithDescription("anaeig", "Analyze the normal modes");
        group.addModule("nmeig");
        group.addModule("nmtraj");
        group.addModule("nmens");
        group.addModule("grompp");
        group.addModuleWithDescription("mdrun",
                                       "Find a potential energy minimum and calculate the Hessian");
    }
}
