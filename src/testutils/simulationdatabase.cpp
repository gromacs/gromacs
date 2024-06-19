/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * \brief
 * Implements declarations from in simulationdatabase.h
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_testutils
 */
#include "gmxpre.h"

#include "testutils/simulationdatabase.h"

#include <algorithm>
#include <iterator>
#include <map>
#include <string>
#include <vector>

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
      { { { "ref-t", "80" },
          { "compressibility", "5e-10" },
          { "tau-p", "1000" },
          { "ensemble-temperature-setting", "constant" },
          { "ensemble-temperature", "80" } },
        { 1, 2, 3, 4 } } },
    // Simple system with 5 water molecules, fairly widely separated
    { "tip3p5",
      { { { "compressibility", "5e-10" },
          { "tau-p", "1000" },
          { "ensemble-temperature-setting", "constant" },
          { "ensemble-temperature", "298" } },
        { 1, 2, 3, 4, 5, 6, 8, 9 } } },
    // Simple system with 5832 argon atoms, suitable for normal pressure coupling
    { "argon5832",
      { { { "ref-t", "80" } },
        { // TODO This test case is not currently used, so we
          // have not tested which rank counts work.
          1,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9 } } },
    // Simple system with 2 nearby water molecules
    { "spc2",
      { {},
        { // TODO This test case is not currently used, so we
          // have not tested which rank counts work.
          1,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9 } } },
    // Simple system with 216 water molecules, condensed phase
    { "spc216",
      { {},
        {
                // TODO This test case is not currently used, so we
                // have not tested which rank counts work.
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9 // TODO tpi test
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
          1,
          2,
          3,
          4,
          5,
          6,
          7,
          8,
          9 } } },
    // Capped alanine peptide in vacuo (no virtual sites)
    { "alanine_vacuo", { { { "constraints", "h-bonds" } }, { 1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
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
          { "free-energy", "yes" },
          { "sc-alpha", "0.5" },
          { "sc-r-power", "6" },
          { "mass-lambdas", "0.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0" },
          { "bonded-lambdas", "0.0 0.0 0.0 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0" },
          { "restraint-lambdas", "0.0 0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0 1.0 1.0" },
          { "vdw-lambdas", "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0" },
          { "coul-lambdas", "0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.5 1.0" },
          { ";couple-moltype", "nonanol" },
          { ";couple-lambda0", "none" },
          { ";couple-lambda1", "vdw-q" },
          { ";couple-intramol", "yes" } },
        { 1, 2, 3, 4, 5, 6, 8, 9 } } },
    // Artificial test system including all virtual site types
    { "vsite_test", { {}, { 1, 2, 3, 4, 5, 6, 7, 8, 9 } } },
};

//! Helper typedef for mdp database
using MdpDatabase = std::map<MdpParameterDatabase, MdpFieldValues>;
//! Database of additional mdp options used for specific algorithms
const MdpDatabase c_additionalMdpOptions{ { MdpParameterDatabase::Default, {} },
                                          { MdpParameterDatabase::Pull,
                                            { { "coulombtype", "reaction-field" },
                                              { "pull", "yes" },
                                              // Prev step reference is checkpointed - rest of pull is not cpt-dependent
                                              { "pull-pbc-ref-prev-step-com", "yes" },
                                              { "pull-ngroups", "2" },
                                              { "pull-group1-name", "FirstWaterMolecule" },
                                              { "pull-group2-name", "SecondWaterMolecule" },
                                              { "pull-ncoords", "1" },
                                              { "pull-coord1-type", "umbrella" },
                                              { "pull-coord1-geometry", "distance" },
                                              { "pull-coord1-groups", "1 2" },
                                              { "pull-coord1-init", "1" },
                                              { "pull-coord1-k", "10000" } } },
                                          { MdpParameterDatabase::Awh,
                                            { { "pull", "yes" },
                                              { "pull-ngroups", "5" },
                                              { "pull-ncoords", "2" },
                                              { "pull-group1-name", "C_&_r_1" },
                                              { "pull-group2-name", "N_&_r_2" },
                                              { "pull-group3-name", "CA" },
                                              { "pull-group4-name", "C_&_r_2" },
                                              { "pull-group5-name", "N_&_r_3" },
                                              { "pull-coord1-geometry", "dihedral" },
                                              { "pull-coord1-groups", "1 2 2 3 3 4" },
                                              { "pull-coord1-k", "4000" },
                                              { "pull-coord1-kB", "1000" },
                                              { "pull-coord2-geometry", "dihedral" },
                                              { "pull-coord2-groups", "2 3 3 4 4 5" },
                                              { "pull-coord2-k", "4000" },
                                              { "pull-coord2-kB", "1000" },
                                              { "pull-coord1-type", "external-potential" },
                                              { "pull-coord1-potential-provider", "awh" },
                                              { "pull-coord2-type", "external-potential" },
                                              { "pull-coord2-potential-provider", "awh" },
                                              { "awh", "yes" },
                                              { "awh-potential", "convolved" },
                                              { "awh-nstout", "4" },
                                              { "awh-nstsample", "4" },
                                              { "awh-nsamples-update", "1" },
                                              { "awh-share-multisim", "no" },
                                              { "awh-nbias", "2" },
                                              { "awh1-ndim", "1" },
                                              { "awh1-dim1-coord-index", "2" },
                                              { "awh1-dim1-start", "150" },
                                              { "awh1-dim1-end", "180" },
                                              { "awh1-dim1-force-constant", "4000" },
                                              { "awh1-dim1-diffusion", "0.1" },
                                              { "awh2-ndim", "1" },
                                              { "awh2-dim1-coord-index", "1" },
                                              { "awh2-dim1-start", "178" },
                                              { "awh2-dim1-end", "-178" },
                                              { "awh2-dim1-force-constant", "4000" },
                                              { "awh2-dim1-diffusion", "0.1" } } },
                                          { MdpParameterDatabase::ExpandedEnsemble,
                                            { { "free-energy", "expanded" },
                                              { "nstcalcenergy", "1" },
                                              { "nstexpanded", "1" },
                                              { "lmc-stats", "weighted-wang-landau" },
                                              { "lmc-move", "metropolized-gibbs" },
                                              { "lmc-seed", "52435" },
                                              // very large wang-landau incrementor to exaggerate effect
                                              { "init-wl-delta", "10" },
                                              { "nstlog", "1" } } } };

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
    mdpFieldValues.insert(MdpField("coulombtype", "Cut-off"));
    mdpFieldValues.insert(MdpField("rcoulomb", "0.7"));
    mdpFieldValues.insert(MdpField("vdwtype", "Cut-off"));
    mdpFieldValues.insert(MdpField("rvdw", "0.7"));
    mdpFieldValues.insert(MdpField("nstcalcenergy", "100"));
    mdpFieldValues.insert(MdpField("rlist", "-1"));
    mdpFieldValues.insert(MdpField("bd-fric", "1000"));
    mdpFieldValues.insert(MdpField("verlet-buffer-tolerance", "0.000001"));
    mdpFieldValues.insert(MdpField("nstlist", "8"));
    mdpFieldValues.insert(MdpField("ld-seed", "234262"));
    mdpFieldValues.insert(MdpField("tau-t", "1"));
    mdpFieldValues.insert(MdpField("tc-grps", "System"));
    mdpFieldValues.insert(MdpField("pcoupltype", "isotropic"));
    mdpFieldValues.insert(MdpField("ref-p", "1"));
    mdpFieldValues.insert(MdpField("constraint-algorithm", "lincs"));
    mdpFieldValues.insert(MdpField("lincs-order", "2"));
    mdpFieldValues.insert(MdpField("lincs-iter", "5"));

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
    return formatAndJoin(
            std::begin(possibleNumbers), std::end(possibleNumbers), ",", StringFormatter("%d"));
}

MdpFieldValues prepareMdpFieldValues(const std::string&   simulationName,
                                     const std::string&   integrator,
                                     const std::string&   tcoupl,
                                     const std::string&   pcoupl,
                                     MdpParameterDatabase additionalMdpParameters)
{
    using MdpField = MdpFieldValues::value_type;

    auto mdpFieldValues = prepareDefaultMdpFieldValues(simulationName);
    mdpFieldValues.insert(MdpField("integrator", integrator));
    mdpFieldValues.insert(MdpField("tcoupl", tcoupl));
    mdpFieldValues.insert(MdpField("pcoupl", pcoupl));
    for (const auto& mdpField : c_additionalMdpOptions.at(additionalMdpParameters))
    {
        // Here, we are overwriting default values - we assume the additional
        // parameters take precedence over the default parameters
        mdpFieldValues[mdpField.first] = mdpField.second;
    }

    return mdpFieldValues;
}

MdpFieldValues prepareMdpFieldValues(const char*          simulationName,
                                     const char*          integrator,
                                     const char*          tcoupl,
                                     const char*          pcoupl,
                                     MdpParameterDatabase additionalMdpParameters)
{
    return prepareMdpFieldValues(
            std::string(simulationName), integrator, tcoupl, pcoupl, additionalMdpParameters);
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
     */
    std::string optionString;
    for (auto const& [key, value] : mdpFieldValues)
    {
        optionString += formatString("%-24s = %s\n", key.c_str(), value.c_str());
    }
    return optionString;
}

} // namespace test
} // namespace gmx
