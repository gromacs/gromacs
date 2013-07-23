/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \internal \brief
 * Registers command-line modules for pre-5.0 binaries.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "legacymodules.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"

#include "gromacs/gmxana/gmx_ana.h"

extern "C"
{

int gmx_gmxcheck(int argc, char *argv[]);
int gmx_gmxdump(int argc, char *argv[]);
int gmx_grompp(int argc, char *argv[]);
int gmx_pdb2gmx(int argc, char *argv[]);
int gmx_protonate(int argc, char *argv[]);
int gmx_tpbconv(int argc, char *argv[]);
int gmx_x2top(int argc, char *argv[]);

}

namespace
{

/*! \brief
 * Convenience function for creating and registering a module.
 *
 * \param[in] manager          Module manager to which to register the module.
 * \param[in] mainFunction     Main function to wrap.
 * \param[in] name             Name for the new module.
 * \param[in] shortDescription One-line description for the new module.
 */
void registerModule(gmx::CommandLineModuleManager                *manager,
                    gmx::CommandLineModuleManager::CMainFunction  mainFunction,
                    const char *name, const char *shortDescription)
{
    manager->addModuleCMain(name, shortDescription, mainFunction);
}

} // namespace

void registerLegacyModules(gmx::CommandLineModuleManager *manager)
{
    // Modules from this directory (were in src/kernel/).
    registerModule(manager, &gmx_gmxcheck, "gmxcheck",
                   "Check and compare files");
    registerModule(manager, &gmx_gmxdump, "gmxdump",
                   "Make binary files human readable");
    registerModule(manager, &gmx_grompp, "grompp",
                   "Make a run input file");
    registerModule(manager, &gmx_pdb2gmx, "pdb2gmx",
                   "Convert coordinate files to topology and FF-compliant coordinate files");
    registerModule(manager, &gmx_tpbconv, "tpbconv",
                   "Make a run input file for restarting a crashed run");

    registerModule(manager, &gmx_protonate, "protonate",
                   "Protonate structures");
    registerModule(manager, &gmx_x2top, "x2top",
                   "Generate a primitive topology from coordinates");

    // Modules from gmx_ana.h.
    registerModule(manager, &gmx_do_dssp, "do_dssp",
                   "Assign secondary structure and calculate solvent accessible surface area");
    registerModule(manager, &gmx_editconf, "editconf",
                   "Convert and manipulates structure files");
    registerModule(manager, &gmx_eneconv, "eneconv",
                   "Convert energy files");
    registerModule(manager, &gmx_genbox, "genbox",
                   "Solvate a system");
    registerModule(manager, &gmx_genconf, "genconf",
                   "Multiply a conformation in 'random' orientations");
    registerModule(manager, &gmx_genion, "genion",
                   "Generate monoatomic ions on energetically favorable positions");
    registerModule(manager, &gmx_genpr, "genrestr",
                   "Generate position restraints or distance restraints for index groups");
    registerModule(manager, &gmx_make_edi, "make_edi",
                   "Generate input files for essential dynamics sampling");
    registerModule(manager, &gmx_make_ndx, "make_ndx",
                   "Make index files");
    registerModule(manager, &gmx_mk_angndx, "mk_angndx",
                   "Generate index files for 'gmx angle'");
    registerModule(manager, &gmx_trjcat, "trjcat",
                   "Concatenate trajectory files");
    registerModule(manager, &gmx_trjconv, "trjconv",
                   "Convert and manipulates trajectory files");
    registerModule(manager, &gmx_trjorder, "trjorder",
                   "Order molecules according to their distance to a group");
    registerModule(manager, &gmx_xpm2ps, "xpm2ps",
                   "Convert XPM (XPixelMap) matrices to postscript or XPM");

    registerModule(manager, &gmx_anadock, "anadock",
                   "Cluster structures from Autodock runs");
    registerModule(manager, &gmx_anaeig, "anaeig",
                   "Analyze eigenvectors/normal modes");
    registerModule(manager, &gmx_analyze, "analyze",
                   "Analyze data sets");
    registerModule(manager, &gmx_g_angle, "angle",
                   "Calculate distributions and correlations for angles and dihedrals");
    registerModule(manager, &gmx_bar, "bar",
                   "Calculate free energy difference estimates through Bennett's acceptance ratio");
    registerModule(manager, &gmx_bond, "bond",
                   "Calculate distances between atoms and bond length distributions");
    registerModule(manager, &gmx_bundle, "bundle",
                   "Analyze bundles of axes, e.g., helices");
    registerModule(manager, &gmx_chi, "chi",
                   "Calculate everything you want to know about chi and other dihedrals");
    registerModule(manager, &gmx_cluster, "cluster",
                   "Cluster structures");
    registerModule(manager, &gmx_clustsize, "clustsize",
                   "Calculate size distributions of atomic clusters");
    registerModule(manager, &gmx_confrms, "confrms",
                   "Fit two structures and calculates the RMSD");
    registerModule(manager, &gmx_covar, "covar",
                   "Calculate and diagonalize the covariance matrix");
    registerModule(manager, &gmx_current, "current",
                   "Calculate dielectric constants and charge autocorrelation function");
    registerModule(manager, &gmx_density, "density",
                   "Calculate the density of the system");
    registerModule(manager, &gmx_densmap, "densmap",
                   "Calculate 2D planar or axial-radial density maps");
    registerModule(manager, &gmx_densorder, "densorder",
                   "Calculate surface fluctuations");
    registerModule(manager, &gmx_dielectric, "dielectric",
                   "Calculate frequency dependent dielectric constants");
    registerModule(manager, &gmx_dipoles, "dipoles",
                   "Compute the total dipole plus fluctuations");
    registerModule(manager, &gmx_disre, "disre",
                   "Analyze distance restraints");
    registerModule(manager, &gmx_dist, "dist",
                   "Calculate distances between centers of mass of two groups");
    registerModule(manager, &gmx_dos, "dos",
                   "Analyze density of states and properties based on that");
    registerModule(manager, &gmx_dyecoupl, "dyecoupl",
                   "Extract dye dynamics from trajectories");
    registerModule(manager, &gmx_dyndom, "dyndom",
                   "Interpolate and extrapolate structure rotations");
    registerModule(manager, &gmx_enemat, "enemat",
                   "Extract an energy matrix from an energy file");
    registerModule(manager, &gmx_energy, "energy",
                   "Writes energies to xvg files and display averages");
    registerModule(manager, &gmx_filter, "filter",
                   "Frequency filter trajectories, useful for making smooth movies");
    registerModule(manager, &gmx_gyrate, "gyrate",
                   "Calculate the radius of gyration");
    registerModule(manager, &gmx_h2order, "h2order",
                   "Compute the orientation of water molecules");
    registerModule(manager, &gmx_hbond, "hbond",
                   "Compute and analyze hydrogen bonds");
    registerModule(manager, &gmx_helix, "helix",
                   "Calculate basic properties of alpha helices");
    registerModule(manager, &gmx_helixorient, "helixorient",
                   "Calculate local pitch/bending/rotation/orientation inside helices");
    registerModule(manager, &gmx_hydorder, "hydorder",
                   "Compute tetrahedrality parameters around a given atom");
    registerModule(manager, &gmx_kinetics, "kinetics",
                   "Analyze kinetic constants from properties based on the Eyring model");
    registerModule(manager, &gmx_lie, "lie",
                   "Estimate free energy from linear combinations");
    registerModule(manager, &gmx_mdmat, "mdmat",
                   "Calculate residue contact maps");
    registerModule(manager, &gmx_mindist, "mindist",
                   "Calculate the minimum distance between two groups");
    registerModule(manager, &gmx_morph, "morph",
                   "Interpolate linearly between conformations");
    registerModule(manager, &gmx_msd, "msd",
                   "Calculates mean square displacements");
    registerModule(manager, &gmx_nmeig, "nmeig",
                   "Diagonalize the Hessian for normal mode analysis");
    registerModule(manager, &gmx_nmens, "nmens",
                   "Generate an ensemble of structures from the normal modes");
    registerModule(manager, &gmx_nmtraj, "nmtraj",
                   "Generate a virtual oscillating trajectory from an eigenvector");
    registerModule(manager, &gmx_options, "options", NULL);
    registerModule(manager, &gmx_order, "order",
                   "Compute the order parameter per atom for carbon tails");
    registerModule(manager, &gmx_pme_error, "pme_error",
                   "Estimate the error of using PME with a given input file");
    registerModule(manager, &gmx_polystat, "polystat",
                   "Calculate static properties of polymers");
    registerModule(manager, &gmx_potential, "potential",
                   "Calculate the electrostatic potential across the box");
    registerModule(manager, &gmx_principal, "principal",
                   "Calculate principal axes of inertia for a group of atoms");
    registerModule(manager, &gmx_rama, "rama",
                   "Compute Ramachandran plots");
    registerModule(manager, &gmx_rdf, "rdf",
                   "Calculate radial distribution functions");
    registerModule(manager, &gmx_rms, "rms",
                   "Calculate RMSDs with a reference structure and RMSD matrices");
    registerModule(manager, &gmx_rmsdist, "rmsdist",
                   "Calculate atom pair distances averaged with power -2, -3 or -6");
    registerModule(manager, &gmx_rmsf, "rmsf",
                   "Calculate atomic fluctuations");
    registerModule(manager, &gmx_rotacf, "rotacf",
                   "Calculate the rotational correlation function for molecules");
    registerModule(manager, &gmx_rotmat, "rotmat",
                   "Plot the rotation matrix for fitting to a reference structure");
    registerModule(manager, &gmx_saltbr, "saltbr",
                   "Compute salt bridges");
    registerModule(manager, &gmx_sans, "sans",
                   "Compute the small angle neutron scattering spectra");
    registerModule(manager, &gmx_sas, "sas",
                   "Compute solvent accessible surface area");
    registerModule(manager, &gmx_saxs, "saxs",
                   "Calculates SAXS structure factors based on Cromer's method");
    registerModule(manager, &gmx_sgangle, "sgangle",
                   "Compute the angle and distance between two groups");
    registerModule(manager, &gmx_sham, "sham",
                   "Compute free energies or other histograms from histograms");
    registerModule(manager, &gmx_sigeps, "sigeps",
                   "Convert c6/12 or c6/cn combinations to and from sigma/epsilon");
    registerModule(manager, &gmx_sorient, "sorient",
                   "Analyze solvent orientation around solutes");
    registerModule(manager, &gmx_spatial, "spatial",
                   "Calculate the spatial distribution function");
    registerModule(manager, &gmx_spol, "spol",
                   "Analyze solvent dipole orientation and polarization around solutes");
    registerModule(manager, &gmx_tcaf, "tcaf",
                   "Calculate viscosities of liquids");
    registerModule(manager, &gmx_traj, "traj",
                   "Plot x, v, f, box, temperature and rotational energy from trajectories");
    registerModule(manager, &gmx_tune_pme, "tune_pme",
                   "Time mdrun as a function of PME nodes to optimize settings");
    registerModule(manager, &gmx_vanhove, "vanhove",
                   "Compute Van Hove correlation functions");
    registerModule(manager, &gmx_velacc, "velacc",
                   "Calculate velocity autocorrelation functions");
    registerModule(manager, &gmx_wham, "wham",
                   "Perform weighted histogram analysis after umbrella sampling");
    registerModule(manager, &gmx_wheel, "wheel",
                   "Plot helical wheels");

    // TODO: Also include binaries from other directories than src/tools/:
    //        "g_xrama|Show animated Ramachandran plots");
    //        "mdrun|finds a potential energy minimum and calculates the Hessian");
    //        "mdrun|performs a simulation, do a normal mode analysis or an energy minimization");
    //        "mdrun|with -rerun (re)calculates energies for trajectory frames");
    //        "ngmx|Display a trajectory");
}
