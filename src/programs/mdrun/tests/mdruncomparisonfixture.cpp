/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
 * Implements classes in mdruncomparisonfixture.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "mdruncomparisonfixture.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace gmx
{
namespace test
{

MdrunComparisonFixture::~MdrunComparisonFixture()
{
}

namespace
{

//! Helper typedef
typedef std::map<std::string, MdrunComparisonFixture::MdpFieldValues> MdpFileValues;

//! Database of .mdp strings that supports MdrunComparisonFixture::prepareSimulation()
MdpFileValues mdpFileValueDatabase_g
{
    // Simple system with 12 argon atoms, fairly widely separated
    {
        "argon12", { {
                         "ref-t", "80"
                     },
                     {
                         "compressibility", "5e-10"
                     },
                     {
                         "tau-p", "1000"
                     } }
    },
    // Simple system with 5 water molecules, fairly widely separated
    {
        "spc5", { {
                      "compressibility", "5e-10"
                  },
                  {
                      "tau-p", "1000"
                  } }
    },
    // Simple system with 5832 argon atoms, suitable for normal pressure coupling
    {
        "argon5832", { {
                           "ref-t", "80"
                       } }
    },
    // Simple system with 216 water molecules, condensed phase
    {
        "spc216", { }
    },
    // Capped alanine peptide in vacuo with virtual sites
    {
        "alanine_vsite_vacuo", { {
                                     "constraints", "all-bonds"
                                 },
                                 {
                                     "compressibility", "5e-10"
                                 },
                                 {
                                     "tau-p", "1000"
                                 } }
    },
    // Capped alanine peptide in aqueous condensed phase, with virtual sites
    {
        "alanine_vsite_solvated", { {
                                        "constraints", "all-bonds"
                                    } }
    },
    // Nonanol molecule in vacuo, topology suitable for FEP testing
    {
        "nonanol_vacuo", { {
                               "nsteps", "16"
                           },
                           {
                               "compressibility", "5e-10"
                           },
                           {
                               "tau-p", "1000"
                           },
                           {
                               "constraints", "h-bonds"
                           },
                           {
                               "other",
                               "free-energy       = yes\n"
                               "sc-alpha          = 0.5\n"
                               "sc-r-power        = 6\n"
                               "nstdhdl           = 4\n"
                               "init-lambda-state = 3\n"
                               "fep_lambdas       = 0.00 0.50 1.00 1.00 1.00\n"
                               "vdw_lambdas       = 0.00 0.00 0.00 0.50 1.00\n"
                               "couple-moltype    = nonanol\n"
                               "couple-lambda0    = vdw-q\n"
                               "couple-lambda1    = none\n"
                               "couple-intramol   = yes\n"
                           } }
    }
};

}       // namespace

MdrunComparisonFixture::MdpFieldValues MdrunComparisonFixture::prepareMdpFieldValues(const char *simulationName)
{
    /* Insert suitable .mdp defaults, so that the database set up
     * above does not need to specify redundant values. This works
     * because std::map.insert() does not over-write elements that
     * already exist.
     *
     * TODO ideally some of these default values should be the same as
     * grompp uses, and sourced from the same place, but that code is
     * a bit of a jungle. */

    typedef MdpFieldValues::value_type MdpField;

    auto &mdpFieldValues = mdpFileValueDatabase_g.at(simulationName);
    mdpFieldValues.insert(MdpField("nsteps", "16"));
    mdpFieldValues.insert(MdpField("ref-t", "298"));
    mdpFieldValues.insert(MdpField("tau-p", "1"));
    mdpFieldValues.insert(MdpField("compressibility", "5e-5"));
    mdpFieldValues.insert(MdpField("constraints", "none"));
    mdpFieldValues.insert(MdpField("other", ""));

    return mdpFieldValues;
}

void MdrunComparisonFixture::prepareMdpFile(const MdpFieldValues &mdpFieldValues,
                                            const char           *integrator,
                                            const char           *tcoupl,
                                            const char           *pcoupl)
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
    runner_.useStringAsMdpFile(formatString("rcoulomb                = 0.7\n"
                                            "rvdw                    = 0.7\n"
                                            "rlist                   = -1\n"
                                            "bd-fric                 = 1000\n"
                                            "verlet-buffer-tolerance = 0.000001\n"
                                            "nsteps                  = %s\n"
                                            "nstenergy               = 4\n"
                                            "nstlist                 = 8\n"
                                            "nstxout                 = 4\n"
                                            "nstvout                 = 4\n"
                                            "nstfout                 = 4\n"
                                            "integrator              = %s\n"
                                            "ld-seed                 = 234262\n"
                                            "tcoupl                  = %s\n"
                                            "ref-t                   = %s\n"
                                            "tau-t                   = 1\n"
                                            "tc-grps                 = System\n"
                                            "pcoupl                  = %s\n"
                                            "pcoupltype              = isotropic\n"
                                            "ref-p                   = 1\n"
                                            "tau-p                   = %s\n"
                                            "compressibility         = %s\n"
                                            "constraints             = %s\n"
                                            "constraint-algorithm    = lincs\n"
                                            "lincs-order             = 2\n"
                                            "lincs-iter              = 5\n"
                                            "%s",
                                            mdpFieldValues.at("nsteps").c_str(),
                                            integrator, tcoupl,
                                            mdpFieldValues.at("ref-t").c_str(),
                                            pcoupl,
                                            mdpFieldValues.at("tau-p").c_str(),
                                            mdpFieldValues.at("compressibility").c_str(),
                                            mdpFieldValues.at("constraints").c_str(),
                                            mdpFieldValues.at("other").c_str()));
}

void MdrunComparisonFixture::runTest(const char            *simulationName,
                                     const char            *integrator,
                                     const char            *tcoupl,
                                     const char            *pcoupl,
                                     FloatingPointTolerance tolerance)
{
    CommandLine caller;
    caller.append("grompp");
    runTest(caller, simulationName, integrator, tcoupl, pcoupl, tolerance);
}

} // namespace test
} // namespace gmx
