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

#include <cstdio>
#include <cstring>

#include <string>

#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/basedefinitions.h"
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
    int         year;
    const char* doi;
} t_citerec;

static constexpr int sc_lineWidth = 79;

void please_cite(FILE* fp, const char* key)
{
    static const t_citerec citedb[] = {
        { "Allen2017",
          "M. P. Allen, D. J. Tildesley",
          "Computer simulation of liquids",
          "Oxford Science Publications",
          2017,
          "10.1093/oso/9780198803195.001.0001" },
        { "Berendsen95a",
          "H. J. C. Berendsen, D. van der Spoel and R. van Drunen",
          "GROMACS: A message-passing parallel molecular dynamics implementation",
          "Comp. Phys. Comm.",
          1995,
          "10.1016/0010-4655(95)00042-E" },
        { "Berendsen84a",
          "H. J. C. Berendsen, J. P. M. Postma, A. DiNola and J. R. Haak",
          "Molecular dynamics with coupling to an external bath",
          "J. Chem. Phys.",
          1984,
          "10.1063/1.448118" },
        { "Ryckaert77a",
          "J. P. Ryckaert, G. Ciccotti, H. J. C. Berendsen",
          "Numerical Integration of the Cartesian Equations of Motion of a System with "
          "Constraints; Molecular Dynamics of n-Alkanes",
          "J. Comp. Phys.",
          1977,
          "10.1016/0021-9991(77)90098-5" },
        { "Miyamoto92a",
          "S. Miyamoto, P. A. Kollman",
          "SETTLE: An Analytical Version of the SHAKE and RATTLE Algorithms for Rigid Water Models",
          "J. Comp. Chem.",
          1992,
          "10.1002/jcc.540130805" },
        { "Cromer1968a",
          "D. T. Cromer, J. B. Mann",
          "X-ray scattering factors computed from numerical Hartree-Fock wave functions",
          "Acta Cryst. A",
          1968,
          "10.1107/S0567739468000550" },
        { "Barth95a",
          "E. Barth, K. Kuczera, B. Leimkuhler, R. D. Skeel",
          "Algorithms for Constrained Molecular Dynamics",
          "J. Comp. Chem.",
          1995,
          "10.1002/jcc.540161003" },
        { "Essmann95a",
          "U. Essmann, L. Perera, M. L. Berkowitz, T. Darden, H. Lee, L. G. Pedersen ",
          "A smooth particle mesh Ewald method",
          "J. Chem. Phys.",
          1995,
          "10.1063/1.470117" },
        { "Torda89a",
          "A. E. Torda, R. M. Scheek,  W. F. van Gunsteren",
          "Time-dependent distance restraints in molecular dynamics simulations",
          "Chem. Phys. Lett.",
          1989,
          "10.1016/0009-2614(89)87249-5" },
        { "Hess97a",
          "B. Hess, H. Bekker, H. J. C. Berendsen, J. G. E. M. Fraaije",
          "LINCS: A Linear Constraint Solver for molecular simulations",
          "J. Comp. Chem.",
          1997,
          "10.1002/(sici)1096-987x(199709)18:12<1463::aid-jcc4>3.0.co;2-h" },
        { "Hess2008a",
          "B. Hess",
          "P-LINCS: A Parallel Linear Constraint Solver for molecular simulation",
          "J. Chem. Theory Comput.",
          2008,
          "10.1021/ct700200b" },
        { "Hess2008b",
          "B. Hess, C. Kutzner, D. van der Spoel, E. Lindahl",
          "GROMACS 4: Algorithms for highly efficient, load-balanced, and scalable molecular "
          "simulation",
          "J. Chem. Theory Comput.",
          2008,
          "10.1021/ct700301q" },
        { "Hub2010",
          "J. S. Hub, B. L. de Groot, D. van der Spoel",
          "g_wham - A free weighted histogram analysis implementation including robust error and "
          "autocorrelation estimates",
          "J. Chem. Theory Comput.",
          2010,
          "10.1021/ct100494z" },
        { "In-Chul99a",
          "Y. In-Chul, M. L. Berkowitz",
          "Ewald summation for systems with slab geometry",
          "J. Chem. Phys.",
          1999,
          "10.1063/1.479595" },
        { "Spoel98a",
          "D. van der Spoel, P. J. van Maaren, H. J. C. Berendsen",
          "A systematic study of water models for molecular simulation. Derivation of models "
          "optimized for use with a reaction-field.",
          "J. Chem. Phys.",
          1998,
          "10.1063/1.476482" },
        { "Wishart98a",
          "D. S. Wishart, A. M. Nip",
          "Protein Chemical Shift Analysis: A Practical Guide",
          "Biochem. Cell Biol.",
          1998,
          "10.1139/bcb-76-2-3-153" },
        { "Maiorov95",
          "V. N. Maiorov, G. M. Crippen",
          "Size-Independent Comparison of Protein Three-Dimensional Structures",
          "PROTEINS: Struct. Funct. Gen.",
          1995,
          "10.1002/prot.340220308" },
        { "Lourenco2013a",
          "Tuanan C. Lourenco, Mariny F. C. Coelho, Teodorico C. Ramalho, David van der "
          "Spoel, Luciano T. Costa",
          "Insights on the Solubility of CO2 in 1-Ethyl-3-methylimidazolium "
          "Bis(trifluoromethylsulfonyl)imide from the Microscopic Point of View",
          "Environ. Sci. Technol.",
          2013,
          "10.1021/es4020986" },
        { "Pascal2011a",
          "T. A. Pascal, S. T. Lin, W. A. Goddard III",
          "Thermodynamics of liquids: standard molar entropies and heat capacities of common "
          "solvents from 2PT molecular dynamics",
          "Phys. Chem. Chem. Phys.",
          2011,
          "10.1039/C0CP01549K" },
        { "Caleman2008a",
          "C. Caleman, D. van der Spoel",
          "Picosecond Melting of Ice by an Infrared Laser Pulse: A Simulation Study",
          "Angew. Chem. Int. Ed",
          2008,
          "10.1002/anie.200703987" },
        { "Caleman2011b",
          "C. Caleman, P. J. van Maaren, M. Hong, J. S. Hub, L. T. da Costa, D. van "
          "der Spoel",
          "Force Field Benchmark of Organic Liquids: Density, Enthalpy of Vaporization, Heat "
          "Capacities, Surface Tension, Isothermal Compressibility, Volumetric Expansion "
          "Coefficient, and Dielectric Constant",
          "J. Chem. Theo. Comp.",
          2012,
          "10.1021/ct200731v" },
        { "Lindahl2001a",
          "E. Lindahl, B. Hess, D. van der Spoel",
          "GROMACS 3.0: A package for molecular simulation and trajectory analysis",
          "J. Mol. Mod.",
          2001,
          "10.1007/s008940100045" },
        { "Eisenberg86a",
          "D. Eisenberg, A. D. McLachlan",
          "Solvation energy in protein folding and binding",
          "Nature",
          1986,
          "10.1038/319199a0" },
        { "Bondi1964a",
          "A. Bondi",
          "van der Waals Volumes and Radii",
          "J. Phys. Chem.",
          1964,
          "10.1021/j100785a001" },
        { "Eisenhaber95",
          "Frank Eisenhaber, Philip Lijnzaad, Patrick Argos, Chris Sander, Michael Scharf",
          "The Double Cube Lattice Method: Efficient Approaches to Numerical Integration of "
          "Surface Area and Volume and to Dot Surface Contouring of Molecular Assemblies",
          "J. Comp. Chem.",
          1995,
          "10.1002/jcc.540160303" },
        { "Hess2002",
          "B. Hess, H. Saint-Martin, H.J.C. Berendsen",
          "Flexible constraints: an adiabatic treatment of quantum degrees of freedom, with "
          "application to the flexible and polarizable MCDHO model for water",
          "J. Chem. Phys.",
          2002,
          "10.1063/1.1478056" },
        { "Hess2003",
          "B. Hess, R.M. Scheek",
          "Orientation restraints in molecular dynamics simulations using time and ensemble "
          "averaging",
          "J. Magn. Res.",
          2003,
          "10.1016/S1090-7807(03)00178-2" },
        { "Mu2005a",
          "Y. Mu, P. H. Nguyen, G. Stock",
          "Energy landscape of a small peptide revealed by dihedral angle principal component "
          "analysis",
          "Prot. Struct. Funct. Bioinf.",
          2005,
          "10.1002/prot.20310" },
        { "Okabe2001a",
          "T. Okabe, M. Kawata, Y. Okamoto, M. Mikami",
          "Replica-exchange {M}onte {C}arlo method for the isobaric-isothermal ensemble",
          "Chem. Phys. Lett.",
          2001,
          "10.1016/S0009-2614(01)00055-0" },
        { "Tropp80a",
          "J. Tropp",
          "Dipolar Relaxation and Nuclear Overhauser effects in nonrigid molecules: The effect of "
          "fluctuating internuclear distances",
          "J. Chem. Phys.",
          1980,
          "10.1063/1.439059" },
        { "Spoel2005a",
          "D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A. E. Mark, H. J. C. Berendsen",
          "GROMACS: Fast, Flexible and Free",
          "J. Comp. Chem.",
          2005,
          "10.1002/jcc.20291" },
        { "Spoel2006b",
          "D. van der Spoel, P. J. van Maaren, P. Larsson, N. Timneanu",
          "Thermodynamics of hydrogen bonding in hydrophilic and hydrophobic media",
          "J. Phys. Chem. B",
          2006,
          "10.1021/jp0572535" },
        { "Bussi2007a",
          "G. Bussi, D. Donadio, M. Parrinello",
          "Canonical sampling through velocity rescaling",
          "J. Chem. Phys.",
          2007,
          "10.1063/1.2408420" },
        { "Hub2006",
          "J. S. Hub, B. L. de Groot",
          "Does CO2 permeate through Aquaporin-1?",
          "Biophys. J.",
          2006,
          "10.1529/biophysj.106.081406" },
        { "Engin2010",
          "O. Engin, A. Villa, M. Sayar, B. Hess",
          "Driving Forces for Adsorption of Amphiphilic Peptides to Air-Water Interface",
          "J. Phys. Chem. B",
          2010,
          "10.1021/jp1024922" },
        { "Wang2010",
          "H. Wang, F. Dommert, C.Holm",
          "Optimizing working parameters of the smooth particle mesh Ewald algorithm in terms of "
          "accuracy and efficiency",
          "J. Chem. Phys. B",
          2010,
          "10.1063/1.3446812" },
        { "Sugita1999a",
          "Y. Sugita, Y. Okamoto",
          "Replica-exchange molecular dynamics method for protein folding",
          "Chem. Phys. Lett.",
          1999,
          "10.1016/S0009-2614(99)01123-9" },
        { "Kutzner2011",
          "C. Kutzner, J. Czub, H. Grubmuller",
          "Keep it Flexible: Driving Macromolecular Rotary Motions in Atomistic Simulations with "
          "GROMACS",
          "J. Chem. Theory Comput.",
          2011,
          "10.1021/ct100666v" },
        { "Hoefling2011",
          "M. Hoefling, N. Lima, D. Haenni, C.A.M. Seidel, B. Schuler, H. Grubmuller",
          "Structural Heterogeneity and Quantitative FRET Efficiency Distributions of Polyprolines "
          "through a Hybrid Atomistic Simulation and Monte Carlo Approach",
          "PLoS ONE",
          2011,
          "10.1371/journal.pone.0019791" },
        { "Hockney1988",
          "R. W. Hockney, J. W. Eastwood",
          "Computer simulation using particles",
          "IOP, Bristol",
          1988,
          "10.1201/9780367806934" },
        { "Ballenegger2012",
          "V. Ballenegger, J.J. Cerda, C. Holm",
          "How to Convert SPME to P3M: Influence Functions and Error Estimates",
          "J. Chem. Theory Comput.",
          2012,
          "10.1021/ct2001792" },
        { "Garmay2012",
          "Garmay Yu, Shvetsov A, Karelov D, Lebedev D, Radulescu A, Petukhov M, Isaev-Ivanov V",
          "Correlated motion of protein subdomains and large-scale conformational flexibility of "
          "RecA protein filament",
          "Journal of Physics: Conference Series",
          2012,
          "10.1088/1742-6596/340/1/012094" },
        { "Kutzner2011b",
          "C. Kutzner, H. Grubmuller, B. L. de Groot, U. Zachariae",
          "Computational Electrophysiology: The Molecular Dynamics of Ion Channel Permeation and "
          "Selectivity in Atomistic Detail",
          "Biophys. J.",
          2011,
          "10.1016/j.bpj.2011.06.010" },
        { "Lundborg2014",
          "M. Lundborg, R. Apostolov, D. Spangberg, A. Gardenas, D. van der Spoel, E. Lindahl",
          "An efficient and extensible format, library, and API for binary trajectory data from "
          "molecular simulations",
          "J. Comput. Chem.",
          2014,
          "10.1002/jcc.23495" },
        { "Goga2012",
          "N. Goga, A. J. Rzepiela, A. H. de Vries, S. J. Marrink, H. J. C. Berendsen",
          "Efficient Algorithms for Langevin and DPD Dynamics",
          "J. Chem. Theory Comput.",
          2012,
          "10.1021/ct3000876" },
        { "Pronk2013",
          "S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R. Apostolov, M. R. Shirts, J. "
          "C. Smith, P. M. Kasson, D. van der Spoel, B. Hess, E. Lindahl",
          "GROMACS 4.5: a high-throughput and highly parallel open source molecular simulation "
          "toolkit",
          "Bioinformatics",
          2013,
          "10.1093/bioinformatics/btt055" },
        { "Pall2015",
          "S. Páll, M. J. Abraham, C. Kutzner, B. Hess, E. Lindahl",
          "Tackling Exascale Software Challenges in Molecular Dynamics Simulations with GROMACS",
          "In S. Markidis & E. Laure (Eds.), Solving Software Challenges for Exascale",
          2015,
          "10.1007/978-3-319-15976-8_1" },
        { "Abraham2015",
          "M. J. Abraham, T. Murtola, R. Schulz, S. Páll, J. C. Smith, B. Hess, E. Lindahl",
          "GROMACS: High performance molecular simulations through multi-level parallelism from "
          "laptops to supercomputers",
          "SoftwareX",
          2015,
          "10.1016/j.softx.2015.06.001" },
        { "Ballenegger2009",
          "V. Ballenegger, A. Arnold, J. J. Cerdà",
          "Simulations of non-neutral slab systems with long-range electrostatic interactions in "
          "two-dimensional periodic boundary conditions",
          "J. Chem. Phys",
          2009,
          "10.1063/1.3216473" },
        { "Hub2014a",
          "J. S. Hub, B. L. de Groot, H. Grubmueller, G. Groenhof",
          "Quantifying Artifacts in Ewald Simulations of Inhomogeneous Systems with a Net Charge",
          "J. Chem. Theory Comput.",
          2014,
          "10.1021/ct400626b" },
        { "Spoel2018a",
          "D. van der Spoel, M. M. Ghahremanpour, J. Lemkul",
          "Small Molecule Thermochemistry: A Tool For Empirical Force Field Development",
          "J. Phys. Chem. A",
          2018,
          "10.1021/acs.jpca.8b09867" },
        { "Lindahl2014",
          "V. Lindahl, J. Lidmar, B. Hess",
          "Accelerated weight histogram method for exploring free energy landscapes",
          "J. Chem. Phys.",
          2014,
          "10.1063/1.4890371" },
        { "Bernetti2020",
          "M. Bernetti, G. Bussi",
          "Pressure control using stochastic cell rescaling",
          "J. Chem. Phys.",
          2020,
          "10.1063/5.0020514" },
        { "Lundborg2021",
          "M. Lundborg, J. Lidmar, B. Hess",
          "The accelerated weight histogram method for alchemical free energy calculations",
          "J. Chem. Phys.",
          2021,
          "10.1063/5.0044352" },
        { "Kabsch1983",
          "W. Kabsch, C. Sander",
          "Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and "
          "geometrical features",
          "Biopolymers",
          1983,
          "10.1002/bip.360221211" },
        { "Shvetsov2013",
          "A. V. Shvetsov, A. E. Schmidt, D. V. Lebedev, V. V. Isaev-Ivanov",
          "Method for calculating small-angle neutron scattering spectra using all-atom molecular "
          "dynamics trajectories",
          "Journal of Surface Investigation. X-ray, Synchrotron and Neutron Techniques",
          2013,
          "10.1134/S1027451013060372" },
        { "Lundborg2023",
          "M. Lundborg, J. Lidmar, B. Hess",
          "On the Path to Optimal Alchemistry",
          "Protein J.",
          2023,
          "10.1007/s10930-023-10137-1" },
        { "Gorelov2024a",
          "S. Gorelov, A. Titov, O. Tolicheva, A. Konevega, A. Shvetsov",
          "DSSP in GROMACS: Tool for Defining Secondary Structures of Proteins in Trajectories",
          "Journal of Chemical Information and Modeling",
          2024,
          "10.1021/acs.jcim.3c01344" },
        { "Gorelov2024b",
          "S. Gorelov, A. Titov, O. Tolicheva, A. Konevega, A. Shvetsov",
          "Determination of Hydrogen Bonds in GROMACS: A New Implementation to Overcome Memory "
          "Limitation",
          "Journal of Chemical Information and Modeling",
          2024,
          "10.1021/acs.jcim.3c02087" },
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
                "%s\n%s\n%s (%d)\nDOI: %s\n",
                author,
                title,
                citedb[index].journal,
                citedb[index].year,
                citedb[index].doi);
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
