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
#include "gmxpre.h"

#include "gromacs/utility/pleasecite.h"

#include <cstring>

#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

typedef struct
{
    const char* key;
    const char* author;
    const char* title;
    const char* journal;
    int         volume, year;
    const char* pages;
} t_citerec;

static constexpr int sc_lineWidth = 79;

void please_cite(FILE* fp, const char* key)
{
    static const t_citerec citedb[] = {
        { "Allen1987a",
          "M. P. Allen and D. J. Tildesley",
          "Computer simulation of liquids",
          "Oxford Science Publications",
          1,
          1987,
          "1" },
        { "Berendsen95a",
          "H. J. C. Berendsen, D. van der Spoel and R. van Drunen",
          "GROMACS: A message-passing parallel molecular dynamics implementation",
          "Comp. Phys. Comm.",
          91,
          1995,
          "43-56" },
        { "Berendsen84a",
          "H. J. C. Berendsen, J. P. M. Postma, A. DiNola and J. R. Haak",
          "Molecular dynamics with coupling to an external bath",
          "J. Chem. Phys.",
          81,
          1984,
          "3684-3690" },
        { "Ryckaert77a",
          "J. P. Ryckaert and G. Ciccotti and H. J. C. Berendsen",
          "Numerical Integration of the Cartesian Equations of Motion of a System with "
          "Constraints; Molecular Dynamics of n-Alkanes",
          "J. Comp. Phys.",
          23,
          1977,
          "327-341" },
        { "Miyamoto92a",
          "S. Miyamoto and P. A. Kollman",
          "SETTLE: An Analytical Version of the SHAKE and RATTLE Algorithms for Rigid Water Models",
          "J. Comp. Chem.",
          13,
          1992,
          "952-962" },
        { "Cromer1968a",
          "D. T. Cromer & J. B. Mann",
          "X-ray scattering factors computed from numerical Hartree-Fock wave functions",
          "Acta Cryst. A",
          24,
          1968,
          "321" },
        { "Barth95a",
          "E. Barth and K. Kuczera and B. Leimkuhler and R. D. Skeel",
          "Algorithms for Constrained Molecular Dynamics",
          "J. Comp. Chem.",
          16,
          1995,
          "1192-1209" },
        { "Essmann95a",
          "U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee and L. G. Pedersen ",
          "A smooth particle mesh Ewald method",
          "J. Chem. Phys.",
          103,
          1995,
          "8577-8592" },
        { "Torda89a",
          "A. E. Torda and R. M. Scheek and W. F. van Gunsteren",
          "Time-dependent distance restraints in molecular dynamics simulations",
          "Chem. Phys. Lett.",
          157,
          1989,
          "289-294" },
        { "Tironi95a",
          "I. G. Tironi and R. Sperb and P. E. Smith and W. F. van Gunsteren",
          "Generalized reaction field method for molecular dynamics simulations",
          "J. Chem. Phys",
          102,
          1995,
          "5451-5459" },
        { "Hess97a",
          "B. Hess and H. Bekker and H. J. C. Berendsen and J. G. E. M. Fraaije",
          "LINCS: A Linear Constraint Solver for molecular simulations",
          "J. Comp. Chem.",
          18,
          1997,
          "1463-1472" },
        { "Hess2008a",
          "B. Hess",
          "P-LINCS: A Parallel Linear Constraint Solver for molecular simulation",
          "J. Chem. Theory Comput.",
          4,
          2008,
          "116-122" },
        { "Hess2008b",
          "B. Hess and C. Kutzner and D. van der Spoel and E. Lindahl",
          "GROMACS 4: Algorithms for highly efficient, load-balanced, and scalable molecular "
          "simulation",
          "J. Chem. Theory Comput.",
          4,
          2008,
          "435-447" },
        { "Hub2010",
          "J. S. Hub, B. L. de Groot and D. van der Spoel",
          "g_wham - A free weighted histogram analysis implementation including robust error and "
          "autocorrelation estimates",
          "J. Chem. Theory Comput.",
          6,
          2010,
          "3713-3720" },
        { "In-Chul99a",
          "Y. In-Chul and M. L. Berkowitz",
          "Ewald summation for systems with slab geometry",
          "J. Chem. Phys.",
          111,
          1999,
          "3155-3162" },
        { "DeGroot97a",
          "B. L. de Groot and D. M. F. van Aalten and R. M. Scheek and A. Amadei and G. Vriend and "
          "H. J. C. Berendsen",
          "Prediction of Protein Conformational Freedom From Distance Constrains",
          "Proteins",
          29,
          1997,
          "240-251" },
        { "Spoel98a",
          "D. van der Spoel and P. J. van Maaren and H. J. C. Berendsen",
          "A systematic study of water models for molecular simulation. Derivation of models "
          "optimized for use with a reaction-field.",
          "J. Chem. Phys.",
          108,
          1998,
          "10220-10230" },
        { "Wishart98a",
          "D. S. Wishart and A. M. Nip",
          "Protein Chemical Shift Analysis: A Practical Guide",
          "Biochem. Cell Biol.",
          76,
          1998,
          "153-163" },
        { "Maiorov95",
          "V. N. Maiorov and G. M. Crippen",
          "Size-Independent Comparison of Protein Three-Dimensional Structures",
          "PROTEINS: Struct. Funct. Gen.",
          22,
          1995,
          "273-283" },
        { "Feenstra99",
          "K. A. Feenstra and B. Hess and H. J. C. Berendsen",
          "Improving Efficiency of Large Time-scale Molecular Dynamics Simulations of "
          "Hydrogen-rich Systems",
          "J. Comput. Chem.",
          20,
          1999,
          "786-798" },
        { "Lourenco2013a",
          "Tuanan C. Lourenco and Mariny F. C. Coelho and Teodorico C. Ramalho and David van der "
          "Spoel and Luciano T. Costa",
          "Insights on the Solubility of CO2 in 1-Ethyl-3-methylimidazolium "
          "Bis(trifluoromethylsulfonyl)imide from the Microscopic Point of View",
          "Environ. Sci. Technol.",
          47,
          2013,
          "7421-7429" },
        { "Timneanu2004a",
          "N. Timneanu and C. Caleman and J. Hajdu and D. van der Spoel",
          "Auger Electron Cascades in Water and Ice",
          "Chem. Phys.",
          299,
          2004,
          "277-283" },
        { "Pascal2011a",
          "T. A. Pascal and S. T. Lin and W. A. Goddard III",
          "Thermodynamics of liquids: standard molar entropies and heat capacities of common "
          "solvents from 2PT molecular dynamics",
          "Phys. Chem. Chem. Phys.",
          13,
          2011,
          "169-181" },
        { "Caleman2008a",
          "C. Caleman and D. van der Spoel",
          "Picosecond Melting of Ice by an Infrared Laser Pulse: A Simulation Study",
          "Angew. Chem. Int. Ed",
          47,
          2008,
          "1417-1420" },
        { "Caleman2011b",
          "C. Caleman and P. J. van Maaren and M. Hong and J. S. Hub and L. T. da Costa and D. van "
          "der Spoel",
          "Force Field Benchmark of Organic Liquids: Density, Enthalpy of Vaporization, Heat "
          "Capacities, Surface Tension, Isothermal Compressibility, Volumetric Expansion "
          "Coefficient, and Dielectric Constant",
          "J. Chem. Theo. Comp.",
          8,
          2012,
          "61" },
        { "Lindahl2001a",
          "E. Lindahl and B. Hess and D. van der Spoel",
          "GROMACS 3.0: A package for molecular simulation and trajectory analysis",
          "J. Mol. Mod.",
          7,
          2001,
          "306-317" },
        { "Wang2001a",
          "J. Wang and W. Wang and S. Huo and M. Lee and P. A. Kollman",
          "Solvation model based on weighted solvent accessible surface area",
          "J. Phys. Chem. B",
          105,
          2001,
          "5055-5067" },
        { "Eisenberg86a",
          "D. Eisenberg and A. D. McLachlan",
          "Solvation energy in protein folding and binding",
          "Nature",
          319,
          1986,
          "199-203" },
        { "Bondi1964a",
          "A. Bondi",
          "van der Waals Volumes and Radii",
          "J. Phys. Chem.",
          68,
          1964,
          "441-451" },
        { "Eisenhaber95",
          "Frank Eisenhaber and Philip Lijnzaad and Patrick Argos and Chris Sander and Michael "
          "Scharf",
          "The Double Cube Lattice Method: Efficient Approaches to Numerical Integration of "
          "Surface Area and Volume and to Dot Surface Contouring of Molecular Assemblies",
          "J. Comp. Chem.",
          16,
          1995,
          "273-284" },
        { "Hess2002",
          "B. Hess, H. Saint-Martin and H.J.C. Berendsen",
          "Flexible constraints: an adiabatic treatment of quantum degrees of freedom, with "
          "application to the flexible and polarizable MCDHO model for water",
          "J. Chem. Phys.",
          116,
          2002,
          "9602-9610" },
        { "Hess2003",
          "B. Hess and R.M. Scheek",
          "Orientation restraints in molecular dynamics simulations using time and ensemble "
          "averaging",
          "J. Magn. Res.",
          164,
          2003,
          "19-27" },
        { "Rappe1991a",
          "A. K. Rappe and W. A. Goddard III",
          "Charge Equillibration for Molecular Dynamics Simulations",
          "J. Phys. Chem.",
          95,
          1991,
          "3358-3363" },
        { "Mu2005a",
          "Y. Mu, P. H. Nguyen and G. Stock",
          "Energy landscape of a small peptide revealed by dihedral angle principal component "
          "analysis",
          "Prot. Struct. Funct. Bioinf.",
          58,
          2005,
          "45-52" },
        { "Okabe2001a",
          "T. Okabe and M. Kawata and Y. Okamoto and M. Mikami",
          "Replica-exchange {M}onte {C}arlo method for the isobaric-isothermal ensemble",
          "Chem. Phys. Lett.",
          335,
          2001,
          "435-439" },
        { "Hukushima96a",
          "K. Hukushima and K. Nemoto",
          "Exchange Monte Carlo Method and Application to Spin Glass Simulations",
          "J. Phys. Soc. Jpn.",
          65,
          1996,
          "1604-1608" },
        { "Tropp80a",
          "J. Tropp",
          "Dipolar Relaxation and Nuclear Overhauser effects in nonrigid molecules: The effect of "
          "fluctuating internuclear distances",
          "J. Chem. Phys.",
          72,
          1980,
          "6035-6043" },
        { "Bultinck2002a",
          "P. Bultinck and W. Langenaeker and P. Lahorte and F. De Proft and P. Geerlings and M. "
          "Waroquier and J. P. Tollenaere",
          "The electronegativity equalization method I: Parametrization and validation for atomic "
          "charge calculations",
          "J. Phys. Chem. A",
          106,
          2002,
          "7887-7894" },
        { "Yang2006b",
          "Q. Y. Yang and K. A. Sharp",
          "Atomic charge parameters for the finite difference Poisson-Boltzmann method using "
          "electronegativity neutralization",
          "J. Chem. Theory Comput.",
          2,
          2006,
          "1152-1167" },
        { "Spoel2005a",
          "D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A. E. Mark and H. J. C. Berendsen",
          "GROMACS: Fast, Flexible and Free",
          "J. Comp. Chem.",
          26,
          2005,
          "1701-1719" },
        { "Spoel2006b",
          "D. van der Spoel, P. J. van Maaren, P. Larsson and N. Timneanu",
          "Thermodynamics of hydrogen bonding in hydrophilic and hydrophobic media",
          "J. Phys. Chem. B",
          110,
          2006,
          "4393-4398" },
        { "Spoel2006d",
          "D. van der Spoel and M. M. Seibert",
          "Protein folding kinetics and thermodynamics from atomistic simulations",
          "Phys. Rev. Letters",
          96,
          2006,
          "238102" },
        { "Palmer94a",
          "B. J. Palmer",
          "Transverse-current autocorrelation-function calculations of the shear viscosity for "
          "molecular liquids",
          "Phys. Rev. E",
          49,
          1994,
          "359-366" },
        { "Bussi2007a",
          "G. Bussi, D. Donadio and M. Parrinello",
          "Canonical sampling through velocity rescaling",
          "J. Chem. Phys.",
          126,
          2007,
          "014101" },
        { "Hub2006",
          "J. S. Hub and B. L. de Groot",
          "Does CO2 permeate through Aquaporin-1?",
          "Biophys. J.",
          91,
          2006,
          "842-848" },
        { "Hub2008",
          "J. S. Hub and B. L. de Groot",
          "Mechanism of selectivity in aquaporins and aquaglyceroporins",
          "PNAS",
          105,
          2008,
          "1198-1203" },
        { "Friedrich2009",
          "M. S. Friedrichs, P. Eastman, V. Vaidyanathan, M. Houston, S. LeGrand, A. L. Beberg, D. "
          "L. Ensign, C. M. Bruns, and V. S. Pande",
          "Accelerating Molecular Dynamic Simulation on Graphics Processing Units",
          "J. Comp. Chem.",
          30,
          2009,
          "864-872" },
        { "Engin2010",
          "O. Engin, A. Villa, M. Sayar and B. Hess",
          "Driving Forces for Adsorption of Amphiphilic Peptides to Air-Water Interface",
          "J. Phys. Chem. B",
          114,
          2010,
          "11093" },
        { "Wang2010",
          "H. Wang, F. Dommert, C.Holm",
          "Optimizing working parameters of the smooth particle mesh Ewald algorithm in terms of "
          "accuracy and efficiency",
          "J. Chem. Phys. B",
          133,
          2010,
          "034117" },
        { "Sugita1999a",
          "Y. Sugita, Y. Okamoto",
          "Replica-exchange molecular dynamics method for protein folding",
          "Chem. Phys. Lett.",
          314,
          1999,
          "141-151" },
        { "Kutzner2011",
          "C. Kutzner and J. Czub and H. Grubmuller",
          "Keep it Flexible: Driving Macromolecular Rotary Motions in Atomistic Simulations with "
          "GROMACS",
          "J. Chem. Theory Comput.",
          7,
          2011,
          "1381-1393" },
        { "Hoefling2011",
          "M. Hoefling, N. Lima, D. Haenni, C.A.M. Seidel, B. Schuler, H. Grubmuller",
          "Structural Heterogeneity and Quantitative FRET Efficiency Distributions of Polyprolines "
          "through a Hybrid Atomistic Simulation and Monte Carlo Approach",
          "PLoS ONE",
          6,
          2011,
          "e19791" },
        { "Hockney1988",
          "R. W. Hockney and J. W. Eastwood",
          "Computer simulation using particles",
          "IOP, Bristol",
          1,
          1988,
          "1" },
        { "Ballenegger2012",
          "V. Ballenegger, J.J. Cerda, and C. Holm",
          "How to Convert SPME to P3M: Influence Functions and Error Estimates",
          "J. Chem. Theory Comput.",
          8,
          2012,
          "936-947" },
        { "Garmay2012",
          "Garmay Yu, Shvetsov A, Karelov D, Lebedev D, Radulescu A, Petukhov M, Isaev-Ivanov V",
          "Correlated motion of protein subdomains and large-scale conformational flexibility of "
          "RecA protein filament",
          "Journal of Physics: Conference Series",
          340,
          2012,
          "012094" },
        { "Kutzner2011b",
          "C. Kutzner, H. Grubmuller, B. L. de Groot, and U. Zachariae",
          "Computational Electrophysiology: The Molecular Dynamics of Ion Channel Permeation and "
          "Selectivity in Atomistic Detail",
          "Biophys. J.",
          101,
          2011,
          "809-817" },
        { "Lundborg2014",
          "M. Lundborg, R. Apostolov, D. Spangberg, A. Gardenas, D. van der Spoel and E. Lindahl",
          "An efficient and extensible format, library, and API for binary trajectory data from "
          "molecular simulations",
          "J. Comput. Chem.",
          35,
          2014,
          "260-269" },
        { "Goga2012",
          "N. Goga and A. J. Rzepiela and A. H. de Vries and S. J. Marrink and H. J. C. Berendsen",
          "Efficient Algorithms for Langevin and DPD Dynamics",
          "J. Chem. Theory Comput.",
          8,
          2012,
          "3637--3649" },
        { "Pronk2013",
          "S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M. R. Shirts, J. "
          "C. Smith, P. M. Kasson, D. van der Spoel, B. Hess, and E. Lindahl",
          "GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation "
          "toolkit",
          "Bioinformatics",
          29,
          2013,
          "845-54" },
        { "Pall2015",
          "S. Páll, M. J. Abraham, C. Kutzner, B. Hess, E. Lindahl",
          "Tackling Exascale Software Challenges in Molecular Dynamics Simulations with GROMACS",
          "In S. Markidis & E. Laure (Eds.), Solving Software Challenges for Exascale",
          8759,
          2015,
          "3-27" },
        { "Abraham2015",
          "M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl",
          "GROMACS: High performance molecular simulations through multi-level parallelism from "
          "laptops to supercomputers",
          "SoftwareX",
          1,
          2015,
          "19-25" },
        { "Ballenegger2009",
          "V. Ballenegger, A. Arnold, J. J. Cerdà",
          "Simulations of non-neutral slab systems with long-range electrostatic interactions in "
          "two-dimensional periodic boundary conditions",
          "J. Chem. Phys",
          131,
          2009,
          "094107" },
        { "Hub2014a",
          "J. S. Hub, B. L. de Groot, H. Grubmueller, G. Groenhof",
          "Quantifying Artifacts in Ewald Simulations of Inhomogeneous Systems with a Net Charge",
          "J. Chem. Theory Comput.",
          10,
          2014,
          "381-393" },
        { "Spoel2018a",
          "D. van der Spoel, M. M. Ghahremanpour, J. Lemkul",
          "Small Molecule Thermochemistry: A Tool For Empirical Force Field Development",
          "J. Phys. Chem. A",
          122,
          2018,
          "8982-8988" },
        { "Lindahl2014",
          "V. Lindahl, J. Lidmar, B. Hess",
          "Accelerated weight histogram method for exploring free energy landscapes",
          "J. Chem. Phys.",
          141,
          2014,
          "044110" },
        { "Spoel2020",
          "D. van der Spoel, H. Henschel, P. J. van Maaren, M. M. Ghahremanpour, L. T. Costa",
          "A potential for molecular simulation of compounds with linear moieties",
          "J. Chem. Phys.",
          153,
          2020,
          "084503" },
        { "Bernetti2020",
          "M. Bernetti, G. Bussi",
          "Pressure control using stochastic cell rescaling",
          "J. Chem. Phys.",
          153,
          2020,
          "114107" },
        { "Lundborg2021",
          "M. Lundborg, J. Lidmar, B. Hess",
          "The accelerated weight histogram method for alchemical free energy calculations",
          "J. Chem. Phys.",
          154,
          2021,
          "204103" },
        { "Kabsch1983",
          "W. Kabsch, C. Sander",
          "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and "
          "geometrical features.",
          "Biopolymers",
          22,
          1983,
          "2577-2637" },
        { "Shvetsov2013",
          "A. V. Shvetsov, A. E. Schmidt, D. V. Lebedev & V. V. Isaev-Ivanov",
          "Method for calculating small-angle neutron scattering spectra using all-atom molecular "
          "dynamics trajectories",
          "Journal of Surface Investigation. X-ray, Synchrotron and Neutron Techniques",
          7,
          2013,
          "1124–1127" },
        { "Lundborg2023",
          "M. Lundborg, J. Lidmar, B. Hess",
          "On the Path to Optimal Alchemistry",
          "Protein J.",
          42,
          2023,
          "477-489" },
        { "Gorelov2024",
          "S. Gorelov, A. Titov, O. Tolicheva, A. Konevega, A. Shvetsov",
          "DSSP in GROMACS: Tool for Defining Secondary Structures of Proteins in Trajectories",
          "Journal of Chemical Information and Modeling",
          0,
          2024,
          "0" },
    };
#define NSTR static_cast<int>(asize(citedb))

    if (fp == nullptr)
    {
        return;
    }

    int index = 0;
    for (; index < NSTR && (strcmp(citedb[index].key, key) != 0); index++) {}

    fprintf(fp, "\n++++ PLEASE READ AND CITE THE FOLLOWING REFERENCE ++++\n");
    if (index < NSTR)
    {
        /* Insert newlines */
        char* author = wrap_lines(citedb[index].author, sc_lineWidth, 0, FALSE);
        char* title  = wrap_lines(citedb[index].title, sc_lineWidth, 0, FALSE);
        fprintf(fp,
                "%s\n%s\n%s %d (%d) pp. %s\n",
                author,
                title,
                citedb[index].journal,
                citedb[index].volume,
                citedb[index].year,
                citedb[index].pages);
        sfree(author);
        sfree(title);
    }
    else
    {
        fprintf(fp, "Entry %s not found in citation database\n", key);
    }
    fprintf(fp, "-------- -------- --- Thank You --- -------- --------\n\n");
    fflush(fp);
}

namespace
{

//! Write a message to \c fp to request citation also of the source-code DOI.
void writeSourceDoi(FILE* fp)
{
    /* Check if we are in release mode or not.
     * TODO The check should properly target something else than
     * the string being empty
     */
    if (strlen(gmxDOI()) == 0)
    {
        /* Not a release build, return without printing anything */
        return;
    }
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setLineLength(sc_lineWidth);
    wrapper.settings().setFirstLineIndent(0);
    const std::string doiString = wrapper.wrapToString(gmxDOI());

    if (fp == nullptr)
    {
        return;
    }
    fprintf(fp, "\n++++ PLEASE CITE THE DOI FOR THIS VERSION OF GROMACS ++++\n");
    fprintf(fp, "%s%s\n", "https://doi.org/", doiString.c_str());
    fprintf(fp, "-------- -------- --- Thank You --- -------- --------\n\n");
    fflush(fp);
}

} // namespace

void pleaseCiteGromacs(FILE* fplog)
{
    if (fplog == nullptr)
    {
        return;
    }

    please_cite(fplog, "Abraham2015");
    please_cite(fplog, "Pall2015");
    please_cite(fplog, "Pronk2013");
    please_cite(fplog, "Hess2008b");
    please_cite(fplog, "Spoel2005a");
    please_cite(fplog, "Lindahl2001a");
    please_cite(fplog, "Berendsen95a");
    writeSourceDoi(fplog);
}
