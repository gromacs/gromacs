/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2018, by the GROMACS development team, led by
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
 * Implements declarations from in mdruncomparison.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun_integration_tests
 */
#include "gmxpre.h"

#include "mdruncomparison.h"

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

namespace
{

//! Helper typedef
using MdpFileValues = std::map<std::string, MdpFieldValues>;

//! Database of .mdp strings that supports prepareDefaultMdpValues()
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
                               R"(free-energy       = yes
                                  sc-alpha          = 0.5
                                  sc-r-power        = 6
                                  nstdhdl           = 4
                                  init-lambda-state = 3
                                  fep_lambdas       = 0.00 0.50 1.00 1.00 1.00
                                  vdw_lambdas       = 0.00 0.00 0.00 0.50 1.00
                                  couple-moltype    = nonanol
                                  couple-lambda0    = vdw-q
                                  couple-lambda1    = none
                                  couple-intramol   = yes)"
                           } }
    }
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
 *          std::out_of_range  if \c simulationName is not in the database */
MdpFieldValues prepareDefaultMdpFieldValues(const char *simulationName)
{
    using MdpField = MdpFieldValues::value_type;

    auto &mdpFieldValues = mdpFileValueDatabase_g.at(simulationName);
    mdpFieldValues.insert(MdpField("nsteps", "16"));
    mdpFieldValues.insert(MdpField("ref-t", "298"));
    mdpFieldValues.insert(MdpField("tau-p", "1"));
    mdpFieldValues.insert(MdpField("compressibility", "5e-5"));
    mdpFieldValues.insert(MdpField("constraints", "none"));
    mdpFieldValues.insert(MdpField("other", ""));

    return mdpFieldValues;
}

}       // namespace

MdpFieldValues
prepareMdpFieldValues(const char *simulationName,
                      const char *integrator,
                      const char *tcoupl,
                      const char *pcoupl)
{
    using MdpField = MdpFieldValues::value_type;

    auto mdpFieldValues = prepareDefaultMdpFieldValues(simulationName);
    mdpFieldValues.insert(MdpField("integrator", integrator));
    mdpFieldValues.insert(MdpField("tcoupl", tcoupl));
    mdpFieldValues.insert(MdpField("pcoupl", pcoupl));
    return mdpFieldValues;
}

std::string
prepareMdpFileContents(const MdpFieldValues &mdpFieldValues)
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
    return formatString(R"(rcoulomb                = 0.7
                           rvdw                    = 0.7
                           rlist                   = -1
                           bd-fric                 = 1000
                           verlet-buffer-tolerance = 0.000001
                           nsteps                  = %s
                           nstenergy               = 4
                           nstlist                 = 8
                           nstxout                 = 4
                           nstvout                 = 4
                           nstfout                 = 4
                           integrator              = %s
                           ld-seed                 = 234262
                           tcoupl                  = %s
                           ref-t                   = %s
                           tau-t                   = 1
                           tc-grps                 = System
                           pcoupl                  = %s
                           pcoupltype              = isotropic
                           ref-p                   = 1
                           tau-p                   = %s
                           compressibility         = %s
                           constraints             = %s
                           constraint-algorithm    = lincs
                           lincs-order             = 2
                           lincs-iter              = 5
                           %s)",
                        mdpFieldValues.at("nsteps").c_str(),
                        mdpFieldValues.at("integrator").c_str(),
                        mdpFieldValues.at("tcoupl").c_str(),
                        mdpFieldValues.at("ref-t").c_str(),
                        mdpFieldValues.at("pcoupl").c_str(),
                        mdpFieldValues.at("tau-p").c_str(),
                        mdpFieldValues.at("compressibility").c_str(),
                        mdpFieldValues.at("constraints").c_str(),
                        mdpFieldValues.at("other").c_str());
}

} // namespace test
} // namespace gmx
