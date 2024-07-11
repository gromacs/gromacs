/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright None- The GROMACS Authors
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
/**
 * \internal \file
 *
 * \brief
 * Implements test of RAMD
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 *
 */

#include "gmxpre.h"
#include "gromacs/mdlib/sighandler.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "programs/mdrun/tests/moduletest.h"

namespace gmx {
namespace test {

class RAMDTestFixture : public MdrunTestFixture
{
    protected:
        RAMDTestFixture();
        ~RAMDTestFixture() override;
};


RAMDTestFixture::RAMDTestFixture()
{}

RAMDTestFixture::~RAMDTestFixture()
{}


//! Test fixture for mdrun with RAMD settings
typedef gmx::test::RAMDTestFixture RAMDTest;

TEST_F(RAMDTest, RAMD_1whhi)
{
    runner_.useTopGroAndNdxFromDatabase("1wdhi");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 4
        tcoupl                   = V-rescale
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
        ramd                     = yes
        ramd-seed                = 9876
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = INH
        ramd-group1-force        = 600.0
        ramd-group1-r-min-dist   = 0.025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 2
        ramd-force-out-freq      = 2
        ramd-group1-ligand-pbcatom   = 1631
        ramd-group1-receptor-pbcatom = 3289
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

TEST_F(RAMDTest, RAMD_membrane)
{
    runner_.useTopGroAndNdxFromDatabase("membrane");
    const std::string mdpContents = R"(
        dt                       = 0.004
        nsteps                   = 4
        tcoupl                   = V-rescale
        tc-grps                  = System
        tau-t                    = 0.5
        ref-t                    = 300
        constraints              = all-bonds
        cutoff-scheme            = Verlet
        ramd                     = yes
        ramd-seed                = 9876
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = IXO
        ramd-group1-force        = 600.0
        ramd-group1-r-min-dist   = 0.025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 2
        ramd-force-out-freq      = 2
        ramd-group1-ligand-pbcatom   = 3255
        ramd-group1-receptor-pbcatom = 6580
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());

    // Do an mdrun with RAMD enabled
    ASSERT_EQ(0, runner_.callMdrun());
}

// heat shock protein 90
TEST_F(RAMDTest, RAMD_hsp90)
{
    runner_.useTopGroAndNdxFromDatabase("hsp90");
    const std::string mdpContents = R"(
        integrator               = md
        dt                       = 0.002
        nsteps                   = 4
        nstlog                   = 1
        nstenergy                = 1
        nstxout-compressed       = 1
        continuation             = yes
        constraints              = h-bonds
        constraint-algorithm     = lincs
        cutoff-scheme            = Verlet
        coulombtype              = PME
        rcoulomb                 = 1.2
        vdwtype                  = Cut-off
        DispCorr                 = EnerPres
        comm-mode                = Linear
        nstcomm                  = 100
        comm_grps                = System
        tinit                    = 0.000
        nstxout                  = 25000
        nstvout                  = 50000
        compressed-x-precision   = 1000
        compressed-x-grps        = SYSTEM
        pbc                      = xyz
        rlist                    = 1.10
        fourierspacing           = 0.12
        pme_order                = 4
        ewald_geometry           = 3d
        ewald-rtol               = 1e-5
        ewald-rtol-lj            = 1e-5
        vdw-modifier             = Potential-shift
        rvdw-switch              = 0.0
        rvdw                     = 1.2
        tcoupl                   = Nose-Hoover
        tc_grps                  = Protein_INH Water_and_ions
        tau_t                    = 1.0  1.0
        ref_t                    = 300  300
        Pcoupl                   = Parrinello-Rahman
        pcoupltype               = semiisotropic
        tau_p                    = 5.0
        compressibility          = 4.5e-5  4.5e-5
        ref_p                    = 1.0     1.0
        gen_vel                  = no
        ramd                     = yes
        ramd-seed                = 989
        ramd-ngroups             = 1
        ramd-group1-receptor     = Protein
        ramd-group1-ligand       = INH
        ramd-group1-force        = 585.0
        ramd-group1-r-min-dist   = 0.0025
        ramd-group1-max-dist     = 4.0
        ramd-eval-freq           = 2
        ramd-force-out-freq      = 0
        ramd-group1-ligand-pbcatom   = 1
        ramd-group1-receptor-pbcatom = 1
    )";
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun());
    gmx_reset_stop_condition();
}

const std::string glyr_mdp_base = R"(
    integrator               = md
    dt                       = 0.002
    nsteps                   = 4
    nstlog                   = 1
    nstenergy                = 1
    nstxout-compressed       = 1
    continuation             = yes
    constraints              = h-bonds
    constraint-algorithm     = lincs
    cutoff-scheme            = Verlet
    coulombtype              = PME
    rcoulomb                 = 1.2
    vdwtype                  = Cut-off
    rvdw                     = 1.2
    DispCorr                 = EnerPres
    tcoupl                   = Nose-Hoover
    tc-grps                  = System
    tau-t                    = 1.0
    ref-t                    = 300
    pcoupl                   = Parrinello-Rahman
    tau-p                    = 5.0
    compressibility          = 4.5e-5
    ref-p                    = 1.0
    ramd                     = yes
    ramd-seed                = 1234
    ramd-ngroups             = 5
    ramd-group1-receptor     = Protein1
    ramd-group1-ligand       = INH1
    ramd-group1-force        = 600
    ramd-group1-max_dist     = 40.0
    ramd-group1-r_min_dist   = 0.0025
    ramd-group2-receptor     = Protein2
    ramd-group2-ligand       = INH2
    ramd-group2-force        = 600
    ramd-group2-max_dist     = 40.0
    ramd-group2-r_min_dist   = 0.0025
    ramd-group3-receptor     = Protein3
    ramd-group3-ligand       = INH3
    ramd-group3-force        = 600
    ramd-group3-max_dist     = 40.0
    ramd-group3-r_min_dist   = 0.0025
    ramd-group4-receptor     = Protein4
    ramd-group4-ligand       = INH4
    ramd-group4-force        = 600
    ramd-group4-max_dist     = 40.0
    ramd-group4-r_min_dist   = 0.0025
    ramd-group5-receptor     = Protein5
    ramd-group5-ligand       = INH5
    ramd-group5-force        = 600
    ramd-group5-max_dist     = 40.0
    ramd-group5-r_min_dist   = 0.0025
    ramd-eval_freq           = 2
    ramd-force_out_freq      = 2
    ramd-old-angle-dist      = no
    ramd-group1-ligand-pbcatom   = 0
    ramd-group1-receptor-pbcatom = 2692
    ramd-group2-ligand-pbcatom   = 0
    ramd-group2-receptor-pbcatom = 8075
    ramd-group3-ligand-pbcatom   = 0
    ramd-group3-receptor-pbcatom = 13458
    ramd-group4-ligand-pbcatom   = 0
    ramd-group4-receptor-pbcatom = 18841
    ramd-group5-ligand-pbcatom   = 0
    ramd-group5-receptor-pbcatom = 24224
)";

TEST_F(RAMDTest, RAMD_GlyR)
{
    runner_.useTopGroAndNdxFromDatabase("glyr");
    const std::string mdpContents = glyr_mdp_base;
    runner_.useStringAsMdpFile(mdpContents);

    EXPECT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(0, runner_.callMdrun());
    gmx_reset_stop_condition();
}

const std::string water4_mdp_base = R"(
    integrator               = md
    dt                       = 0.001
    nsteps                   = 100000
    nstlog                   = 10
    rlist                    = 1.0
    coulombtype              = Cut-off
    rcoulomb-switch          = 0
    rcoulomb                 = 1.0
    epsilon-r                = 1
    epsilon-rf               = 1
    vdw-type                 = Cut-off
    rvdw-switch              = 0
    rvdw                     = 1.0
    DispCorr                 = no
    Tcoupl                   = no
    Pcoupl                   = no

    ramd                     = yes
    ramd-seed                = 1234
    ramd-eval-freq           = 10
    ramd-force-out-freq      = 10
    ramd-old-angle-dist      = no
    ramd-ngroups             = 2
    ramd-group1-receptor     = 1SOL
    ramd-group1-ligand       = 2SOL
    ramd-group1-force        = 100
    ramd-group1-max-dist     = 1.0
    ramd-group1-r-min-dist   = 0.0025
    ramd-group2-receptor     = 1SOL
    ramd-group2-ligand       = 3SOL
    ramd-group2-force        = 100
    ramd-group2-max-dist     = 1.0
    ramd-group2-r-min-dist   = 0.0025
)";

TEST_F(RAMDTest, RAMD_connected_ligands)
{
    runner_.useTopGroAndNdxFromDatabase("4water");
    auto mdpContents = water4_mdp_base + R"(
        ramd-connected-ligands = no
    )";
    runner_.useStringAsMdpFile(mdpContents);

    CommandLine caller;
    caller.addOption("-ramd");
    caller.addOption("-reprod");

    EXPECT_EQ(0, runner_.callGrompp());
    ASSERT_EQ(2, runner_.callMdrun(caller));
    gmx_reset_stop_condition();

    TextReader reader_log(runner_.logFileName_);
    std::string line;
    int number_of_steps = -1;
    while (reader_log.readLine(&line)) {
        if (line.find("==== RAMD ==== GROMACS will be stopped after") != std::string::npos) {
            number_of_steps = stoi(gmx::splitString(line)[8]);
        }
    }
    EXPECT_EQ(number_of_steps, 630);

    TextReader reader_pullx(fileManager_.getTemporaryFilePath("state_pullx.xvg"));
    // std::cout << reader_pullx.readAll();
    while (reader_pullx.readLine(&line)) {
        if (line.rfind("0.000", 0) != std::string::npos) {
            EXPECT_EQ(std::string("0.0593702"), gmx::splitString(line)[1]);
        }
    }

    TextReader reader_ramd(fileManager_.getTemporaryFilePath("state.xvg"));
    // std::cout << reader_ramd.readAll();
    while (reader_ramd.readLine(&line)) {
        if (line.rfind("0.000", 0) != std::string::npos) {
            EXPECT_EQ(std::string("0.423961"), gmx::splitString(line)[1]);
        }
    }
}

} // namespace test
} // namespace gmx
