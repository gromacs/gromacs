/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
 */
/*! \internal \brief
 * Implements command-line modules for pre-5.0 binaries.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
#include "legacymodules.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"

#include "gromacs/legacyheaders/gmx_ana.h"

namespace
{

/*! \internal \brief
 * Command-line module for wrapping pre-5.0 binaries.
 *
 * Implements a gmx::CommandLineModuleInterface, given a function with
 * C/C++ main signature.
 */
class LegacyCmdLineWrapper : public gmx::CommandLineModuleInterface
{
    public:
        //! Function pointer type for the main function of the module.
        typedef int (*MainFunction)(int argc, char *argv[]);

        /*! \brief
         * Convenience function for creating and registering a module.
         *
         * \param[in] manager  Module manager to which to register the module.
         * \param[in] main     Main function to wrap.
         * \param[in] name     Name for the new module.
         * \param[in] shortDescription One-line description for the new module.
         */
        static void registerModule(gmx::CommandLineModuleManager *manager,
                                   MainFunction main, const char *name,
                                   const char *shortDescription)
        {
            gmx::CommandLineModulePointer module(
                    new LegacyCmdLineWrapper(main, name, shortDescription));
            manager->addModule(gmx::move(module));
        }

        /*! \brief
         * Creates a wrapper module for the given main function.
         *
         * \see registerModule()
         */
        LegacyCmdLineWrapper(MainFunction main, const char *name,
                             const char *shortDescription)
            : main_(main), name_(name), shortDescription_(shortDescription)
        {
        }

        virtual const char *name() const
        {
            return name_;
        }
        virtual const char *shortDescription() const
        {
            return shortDescription_;
        }

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const gmx::HelpWriterContext &context) const;

    private:
        MainFunction            main_;
        const char             *name_;
        const char             *shortDescription_;

};

int LegacyCmdLineWrapper::run(int argc, char *argv[])
{
    return main_(argc, argv);
}

void LegacyCmdLineWrapper::writeHelp(const gmx::HelpWriterContext &context) const
{
    if (context.outputFormat() != gmx::eHelpOutputFormat_Console)
    {
        GMX_THROW(gmx::NotImplementedError(
                          "Command-line help is not implemented for this output format"));
    }
    char *argv[2];
    // TODO: The constness should not be cast away.
    argv[0] = const_cast<char *>(name_);
    argv[1] = const_cast<char *>("-h");
    main_(2, argv);
}

} // namespace

void registerLegacyModules(gmx::CommandLineModuleManager *manager)
{
    LegacyCmdLineWrapper::registerModule(manager, &gmx_do_dssp, "do_dssp",
                                         "Assign secondary structure and calculate solvent accessible surface area");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_editconf, "editconf",
                                         "Convert and manipulates structure files");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_eneconv, "eneconv",
                                         "Convert energy files");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_genbox, "genbox",
                                         "Solvate a system");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_genconf, "genconf",
                                         "Multiply a conformation in 'random' orientations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_genion, "genion",
                                         "Generate monoatomic ions on energetically favorable positions");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_genpr, "genrestr",
                                         "Generate position restraints or distance restraints for index groups");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_make_edi, "make_edi",
                                         "Generate input files for essential dynamics sampling");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_make_ndx, "make_ndx",
                                         "Make index files");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_mk_angndx, "mk_angndx",
                                         "Generate index files for 'gmx angle'");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_trjcat, "trjcat",
                                         "Concatenate trajectory files");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_trjconv, "trjconv",
                                         "Convert and manipulates trajectory files");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_trjorder, "trjorder",
                                         "Order molecules according to their distance to a group");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_xpm2ps, "xpm2ps",
                                         "Convert XPM (XPixelMap) matrices to postscript or XPM");

    // TODO: Include remaining binaries from src/tools/.
    // These are commented out below, and have some issues to consider how to
    // best handle them.
    LegacyCmdLineWrapper::registerModule(manager, &gmx_anadock, "anadock",
                                         "Cluster structures from Autodock runs");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_anaeig, "anaeig",
                                         "Analyze eigenvectors/normal modes");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_analyze, "analyze",
                                         "Analyze data sets");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_g_angle, "angle",
                                         "Calculate distributions and correlations for angles and dihedrals");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_bar, "bar",
                                         "Calculate free energy difference estimates through Bennett's acceptance ratio");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_bond, "bond",
                                         "Calculate distances between atoms and bond length distributions");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_bundle, "bundle",
                                         "Analyze bundles of axes, e.g., helices");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_chi, "chi",
                                         "Calculate everything you want to know about chi and other dihedrals");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_cluster, "cluster",
                                         "Cluster structures");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_clustsize, "clustsize",
                                         "Calculate size distributions of atomic clusters");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_confrms, "confrms",
                                         "Fit two structures and calculates the RMSD");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_covar, "covar",
                                         "Calculate and diagonalize the covariance matrix");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_current, "current",
                                         "Calculate dielectric constants and charge autocorrelation function");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_density, "density",
                                         "Calculate the density of the system");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_densmap, "densmap",
                                         "Calculate 2D planar or axial-radial density maps");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_densorder, "densorder",
                                         "Calculate surface fluctuations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dielectric, "dielectric",
                                         "Calculate frequency dependent dielectric constants");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dipoles, "dipoles",
                                         "Compute the total dipole plus fluctuations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_disre, "disre",
                                         "Analyze distance restraints");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dist, "dist",
                                         "Calculate distances between centers of mass of two groups");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dos, "dos",
                                         "Analyze density of states and properties based on that");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dyecoupl, "dyecoupl",
                                         "Extract dye dynamics from trajectories");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dyndom, "dyndom",
                                         "Interpolate and extrapolate structure rotations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_enemat, "enemat",
                                         "Extract an energy matrix from an energy file");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_energy, "energy",
                                         "Writes energies to xvg files and display averages");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_filter, "filter",
                                         "Frequency filter trajectories, useful for making smooth movies");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_gyrate, "gyrate",
                                         "Calculate the radius of gyration");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_h2order, "h2order",
                                         "Compute the orientation of water molecules");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_hbond, "hbond",
                                         "Compute and analyze hydrogen bonds");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_helix, "helix",
                                         "Calculate basic properties of alpha helices");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_helixorient, "helixorient",
                                         "Calculate local pitch/bending/rotation/orientation inside helices");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_hydorder, "hydorder",
                                         "Compute tetrahedrality parameters around a given atom");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_kinetics, "kinetics",
                                         "Analyze kinetic constants from properties based on the Eyring model");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_lie, "lie",
                                         "Estimate free energy from linear combinations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_mdmat, "mdmat",
                                         "Calculate residue contact maps");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_mindist, "mindist",
                                         "Calculate the minimum distance between two groups");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_morph, "morph",
                                         "Interpolate linearly between conformations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_msd, "msd",
                                         "Calculates mean square displacements");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_nmeig, "nmeig",
                                         "Diagonalize the Hessian for normal mode analysis");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_nmens, "nmens",
                                         "Generate an ensemble of structures from the normal modes");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_nmtraj, "nmtraj",
                                         "Generate a virtual oscillating trajectory from an eigenvector");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_options, "options",
                                         NULL);
    LegacyCmdLineWrapper::registerModule(manager, &gmx_order, "order",
                                         "Compute the order parameter per atom for carbon tails");
    //LegacyCmdLineWrapper::registerModule(manager, &gmx_pme_error, "pme_error",
    //        "Estimate the error of using PME with a given input file");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_polystat, "polystat",
                                         "Calculate static properties of polymers");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_potential, "potential",
                                         "Calculate the electrostatic potential across the box");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_principal, "principal",
                                         "Calculate principal axes of inertia for a group of atoms");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rama, "rama",
                                         "Compute Ramachandran plots");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rdf, "rdf",
                                         "Calculate radial distribution functions");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rms, "rms",
                                         "Calculate RMSDs with a reference structure and RMSD matrices");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rmsdist, "rmsdist",
                                         "Calculate atom pair distances averaged with power -2, -3 or -6");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rmsf, "rmsf",
                                         "Calculate atomic fluctuations");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rotacf, "rotacf",
                                         "Calculate the rotational correlation function for molecules");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_rotmat, "rotmat",
                                         "Plot the rotation matrix for fitting to a reference structure");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_saltbr, "saltbr",
                                         "Compute salt bridges");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sans, "sans",
                                         "Compute the small angle neutron scattering spectra");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sas, "sas",
                                         "Compute solvent accessible surface area");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_saxs, "saxs",
                                         "Calculates SAXS structure factors based on Cromer's method");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sgangle, "sgangle",
                                         "Compute the angle and distance between two groups");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sham, "sham",
                                         "Compute free energies or other histograms from histograms");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sigeps, "sigeps",
                                         "Convert c6/12 or c6/cn combinations to and from sigma/epsilon");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_sorient, "sorient",
                                         "Analyze solvent orientation around solutes");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_spatial, "spatial",
                                         "Calculate the spatial distribution function");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_spol, "spol",
                                         "Analyze solvent dipole orientation and polarization around solutes");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_tcaf, "tcaf",
                                         "Calculate viscosities of liquids");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_traj, "traj",
                                         "Plot x, v, f, box, temperature and rotational energy from trajectories");
    //LegacyCmdLineWrapper::registerModule(manager, &gmx_tune_pme, "tune_pme",
    //        "Time mdrun as a function of PME nodes to optimize settings");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_vanhove, "vanhove",
                                         "Compute Van Hove correlation functions");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_velacc, "velacc",
                                         "Calculate velocity autocorrelation functions");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_wham, "wham",
                                         "Perform weighted histogram analysis after umbrella sampling");
    LegacyCmdLineWrapper::registerModule(manager, &gmx_wheel, "wheel",
                                         "Plot helical wheels");

    // TODO: Also include binaries from other directories than src/tools/:
    //        "g_protonate|Protonate structures");
    //        "g_x2top|Generate a primitive topology from coordinates ");
    //        "g_xrama|Show animated Ramachandran plots");
    //        "gmxcheck|Check and compare files");
    //        "gmxdump|Make binary files human readable");
    //        "grompp|Make a run input file");
    //        "mdrun|finds a potential energy minimum and calculates the Hessian");
    //        "mdrun|performs a simulation, do a normal mode analysis or an energy minimization");
    //        "mdrun|with -rerun (re)calculates energies for trajectory frames");
    //        "ngmx|Display a trajectory");
    //        "pdb2gmx|Convert coordinate files to topology and FF-compliant coordinate files");
    //        "tpbconv|Make a run input file for restarting a crashed run");
}
