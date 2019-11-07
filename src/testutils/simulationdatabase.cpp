/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018,2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements declarations from in simulationdatabase.h
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "simulationdatabase.h"

#include <algorithm>
#include <map>
#include <string>

#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

namespace
{

struct DatabaseEntry
{
    MdpFieldValues   mdpFieldValues;
    std::vector<int> validPpRankCounts;
};

//! Helper typedef
using MdpFileValues = std::map<std::string, DatabaseEntry>;

//! Database of .mdp strings that supports prepareDefaultMdpValues()
const MdpFileValues mdpFileValueDatabase_g{
    // Simple system with 12 argon atoms, fairly widely separated
    { "argon12",
      { { { "ref-t", "80" }, { "compressibility", "5e-10" }, { "tau-p", "1000" } }, { 1, 2, 3, 4 } } },
    // Simple system with 5 water molecules, fairly widely separated
    { "tip3p5", { { { "compressibility", "5e-10" }, { "tau-p", "1000" } }, { 1, 2, 3, 4, 5, 6, 8, 9 } } },
    // Simple system with 5832 argon atoms, suitable for normal pressure coupling
    { "argon5832",
      { { { "ref-t", "80" } },
        { // TODO This test case is not currently used, so we
          // have not tested which rank counts work.
          1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
    // Simple system with 2 nearby water molecules
    { "spc2",
      { {},
        { // TODO This test case is not currently used, so we
          // have not tested which rank counts work.
          1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
    // Simple system with 216 water molecules, condensed phase
    { "spc216",
      { {},
        {
                // TODO This test case is not currently used, so we
                // have not tested which rank counts work.
                1, 2, 3, 4, 5, 6, 7, 8, 9 // TODO tpi test
        } } },
    // Capped alanine peptide in vacuo with virtual sites
    { "alanine_vsite_vacuo",
      { { { "constraints", "all-bonds" }, { "compressibility", "5e-10" }, { "tau-p", "1000" } },
        { 1, 2, 3, 4, 6, 9 } } },
    // Capped alanine peptide in aqueous condensed phase, with virtual sites
    { "alanine_vsite_solvated",
      { { { "constraints", "all-bonds" }, { "compressibility", "5e-10" }, { "tau-p", "1000" } },
        { // TODO This test case is not currently used, so we
          // have not tested which rank counts work.
          1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
    // Zwitterionic glycine in vacuo
    { "glycine_vacuo", { { { "constraints", "h-bonds" } }, { 1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
    // Zwitterionic glycine in vacuo, without constraints
    { "glycine_no_constraints_vacuo", { { { "constraints", "none" } }, { 1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
    // Simple mdrun tests of energy
    { "angles1", { {}, { 1, 2 } } },
    // Scaled water for NMA
    { "scaled-water", { {}, { 1, 2, 3, 4, 5, 6 } } },
    // Villin for NMA
    { "villin", { {}, { 1, 2, 3, 4, 5, 6 } } },
    // SPC-Dimer for NMA
    { "spc-dimer", { {}, { 1, 2, 3, 4, 5, 6 } } },
    // SW-Dimer for NMA
    { "sw-dimer", { { { "nstcalcenergy", "1" } }, { 1, 2, 3, 4, 5, 6 } } },
    // TIP5P for NMA
    { "one-tip5p", { {}, { 1, 2, 3, 4, 5, 6 } } },
    // ICE-Binding protein for NMA
    { "ice-binding", { {}, { 1, 2, 3, 4, 5, 6 } } },
    // Nonanol molecule in vacuo, topology suitable for testing FEP
    // on KE, angles, dihedral restraints, coulomb and vdw
    { "nonanol_vacuo",
      { { { "nsteps", "16" },
          { "compressibility", "5e-10" },
          { "tau-p", "1000" },
          { "constraints", "h-bonds" },
          { "other",
            R"(free-energy         = yes
                                  sc-alpha            = 0.5
                                  sc-r-power          = 6
                                  mass-lambdas        = 0.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
                                  bonded-lambdas      = 0.0 0.0 0.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
                                  restraint-lambdas   = 0.0 0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0 1.0 1.0
                                  vdw-lambdas         = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0
                                  coul-lambdas        = 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 1.0
                                  ;couple-moltype      = nonanol
                                  ;couple-lambda0      = none
                                  ;couple-lambda1      = vdw-q
                                  ;couple-intramol     = yes)" } },
        { 1, 2, 3, 4, 5, 6, 8, 9 } } }
};

/*! \brief Prepare default .mdp values
 *
 * Insert suitable .mdp defaults, so that \c mdpFileValueDatabase_g
 * does not need to specify repetitive values. This works because
 * std::map.insert() does not over-write elements that already exist.
 *
 * \todo ideally some of these default values should be the same as
 * grompp uses, and sourced from the same place, but that code is a
 * bit of a jungle until we transition to using IMdpOptions more.
 *
 * \throws  std::bad_alloc     if out of memory
 *          std::out_of_range  if \c simulationName is not in the database
 *
 * Note: Any mdp options that are not added here cannot be used
 */
MdpFieldValues prepareDefaultMdpFieldValues(const std::string& simulationName)
{
    using MdpField = MdpFieldValues::value_type;

    auto mdpFieldValues = mdpFileValueDatabase_g.at(simulationName).mdpFieldValues;
    mdpFieldValues.insert(MdpField("nsteps", "16"));
    mdpFieldValues.insert(MdpField("nstenergy", "4"));
    mdpFieldValues.insert(MdpField("nstxout", "4"));
    mdpFieldValues.insert(MdpField("nstvout", "4"));
    mdpFieldValues.insert(MdpField("nstfout", "4"));
    mdpFieldValues.insert(MdpField("nstxout-compressed", "0"));
    mdpFieldValues.insert(MdpField("nstdhdl", "4"));
    mdpFieldValues.insert(MdpField("comm-mode", "linear"));
    mdpFieldValues.insert(MdpField("nstcomm", "4"));
    mdpFieldValues.insert(MdpField("ref-t", "298"));
    mdpFieldValues.insert(MdpField("nsttcouple", "4"));
    mdpFieldValues.insert(MdpField("tau-p", "1"));
    mdpFieldValues.insert(MdpField("nstpcouple", "4"));
    mdpFieldValues.insert(MdpField("compressibility", "5e-5"));
    mdpFieldValues.insert(MdpField("constraints", "none"));
    mdpFieldValues.insert(MdpField("other", ""));
    mdpFieldValues.insert(MdpField("coulombtype", "Cut-off"));
    mdpFieldValues.insert(MdpField("rcoulomb", "0.7"));
    mdpFieldValues.insert(MdpField("vdwtype", "Cut-off"));
    mdpFieldValues.insert(MdpField("rvdw", "0.7"));
    mdpFieldValues.insert(MdpField("nstcalcenergy", "100"));

    return mdpFieldValues;
}

} // namespace

bool isNumberOfPpRanksSupported(const std::string& simulationName, int possibleNumberOfPpRanks)
{
    const auto& possibleNumbers = mdpFileValueDatabase_g.at(simulationName).validPpRankCounts;
    return (std::find(std::begin(possibleNumbers), std::end(possibleNumbers), possibleNumberOfPpRanks)
            != std::end(possibleNumbers));
}

std::string reportNumbersOfPpRanksSupported(const std::string& simulationName)
{
    const auto& possibleNumbers = mdpFileValueDatabase_g.at(simulationName).validPpRankCounts;
    return formatAndJoin(std::begin(possibleNumbers), std::end(possibleNumbers), ",",
                         StringFormatter("%d"));
}

MdpFieldValues prepareMdpFieldValues(const std::string& simulationName,
                                     const std::string& integrator,
                                     const std::string& tcoupl,
                                     const std::string& pcoupl)
{
    using MdpField = MdpFieldValues::value_type;

    auto mdpFieldValues = prepareDefaultMdpFieldValues(simulationName);
    mdpFieldValues.insert(MdpField("integrator", integrator));
    mdpFieldValues.insert(MdpField("tcoupl", tcoupl));
    mdpFieldValues.insert(MdpField("pcoupl", pcoupl));
    return mdpFieldValues;
}

MdpFieldValues prepareMdpFieldValues(const char* simulationName,
                                     const char* integrator,
                                     const char* tcoupl,
                                     const char* pcoupl)
{
    return prepareMdpFieldValues(std::string(simulationName), integrator, tcoupl, pcoupl);
}
std::string prepareMdpFileContents(const MdpFieldValues& mdpFieldValues)
{
    /* Set up an .mdp file that permits a highly reproducible
     * simulation. The format string needs to be configured with
     * values for various .mdp fields to suit the kind of system
     * used and testing needed. It also
     * - writes frames from different kinds of steps: starting, ending, intermediate NS, intermediate non-NS
     * - has other steps between frame-writing steps
     * - has enough buffer that e.g. a rerun will compute the same potential energy even though it does NS every frame
     * - has very slow pressure coupling and weak compressibility (since otherwise the box will shrink too fast)
     * - can have arbitrary chunks of .mdp content appended to it (but it is up to grompp how it deals with duplicate fields)
     * - sets random seeds to known values
     * - uses cutoffs that fit inside boxes even after GPU mdrun scales rlist
     *
     * Note that forces computed in the absence of energy computations
     * generally follow a different code path from those computed with
     * energies. However a rerun always computes energies, so we don't
     * currently have a good way to compare forces at steps where
     * energies were not computed with those from rerun on the same
     * coordinates.
     *
     * Note: Any mdp options that are not printed here cannot be used
     */
    return formatString(
            R"(coulombtype             = %s
                           rcoulomb                = %s
                           vdwtype                 = %s
                           rvdw                    = %s
                           rlist                   = -1
                           bd-fric                 = 1000
                           verlet-buffer-tolerance = 0.000001
                           nsteps                  = %s
                           nstenergy               = %s
                           nstxout                 = %s
                           nstvout                 = %s
                           nstfout                 = %s
                           nstxout-compressed      = %s
                           nstdhdl                 = %s
                           nstlist                 = 8
                           integrator              = %s
                           ld-seed                 = 234262
                           tcoupl                  = %s
                           nsttcouple              = %s
                           ref-t                   = %s
                           tau-t                   = 1
                           tc-grps                 = System
                           pcoupl                  = %s
                           nstpcouple              = %s
                           pcoupltype              = isotropic
                           ref-p                   = 1
                           tau-p                   = %s
                           compressibility         = %s
                           constraints             = %s
                           constraint-algorithm    = lincs
                           lincs-order             = 2
                           lincs-iter              = 5
                           nstcalcenergy           = %s
                           comm-mode               = %s
                           nstcomm                 = %s
                           %s)",
            mdpFieldValues.at("coulombtype").c_str(), mdpFieldValues.at("rcoulomb").c_str(),
            mdpFieldValues.at("vdwtype").c_str(), mdpFieldValues.at("rvdw").c_str(),
            mdpFieldValues.at("nsteps").c_str(), mdpFieldValues.at("nstenergy").c_str(),
            mdpFieldValues.at("nstxout").c_str(), mdpFieldValues.at("nstvout").c_str(),
            mdpFieldValues.at("nstfout").c_str(), mdpFieldValues.at("nstxout-compressed").c_str(),
            mdpFieldValues.at("nstdhdl").c_str(), mdpFieldValues.at("integrator").c_str(),
            mdpFieldValues.at("tcoupl").c_str(), mdpFieldValues.at("nsttcouple").c_str(),
            mdpFieldValues.at("ref-t").c_str(), mdpFieldValues.at("pcoupl").c_str(),
            mdpFieldValues.at("nstpcouple").c_str(), mdpFieldValues.at("tau-p").c_str(),
            mdpFieldValues.at("compressibility").c_str(), mdpFieldValues.at("constraints").c_str(),
            mdpFieldValues.at("nstcalcenergy").c_str(), mdpFieldValues.at("comm-mode").c_str(),
            mdpFieldValues.at("nstcomm").c_str(), mdpFieldValues.at("other").c_str());
}

} // namespace test
} // namespace gmx
