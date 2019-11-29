/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2019
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include <map>

#include <math.h>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "programs/alexandria/babel_io.h"
#include "programs/alexandria/fill_inputrec.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"
#include "programs/alexandria/qgen_acm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "poldata_utils.h"

namespace alexandria
{

namespace
{

enum informat{
    einfLOG = 0,
    einfPDB = 1, 
    einfNR
};

class AtomtypeTest : public gmx::test::CommandLineTestBase
{
    protected:
        void testAtype(const char *molname)
        {                       
            gmx::test::TestReferenceChecker checker_(this->rootChecker());
            int         maxpot   = 100; 
            int         nsymm    = 0;
            const char *conf     = (char *)"minimum";
            const char *basis    = (char *)"";
            const char *jobtype  = (char *)"Opt";
            
            std::string           dataName;
            alexandria::MolProp   molprop;

            dataName = gmx::test::TestFileManager::getInputFilePath(molname);
            
            bool readOK = readBabel(dataName.c_str(), &molprop, molname, molname, 
                                    conf, basis, maxpot, nsymm, jobtype, 0.0);
            EXPECT_TRUE(readOK);
            if (readOK)
            {
                std::vector<std::string> atypes;
                auto exper = molprop.BeginExperiment();
                for(auto ca = exper->BeginAtom(); ca < exper->EndAtom(); ++ca)
                {
                    atypes.push_back(ca->getObtype());
                }
                checker_.checkInteger(static_cast<int>(atypes.size()), molname);
                checker_.checkSequence(atypes.begin(), atypes.end(), "atomtypes");
            }
        }
};

TEST_F (AtomtypeTest, Butanol)
{
    testAtype("1-butanol.pdb");
}

TEST_F (AtomtypeTest, NitroEthyne)
{
    testAtype("1-nitroethyne.pdb");
}

TEST_F (AtomtypeTest, Guanidine)
{
    testAtype("guanidine.pdb");
}

TEST_F (AtomtypeTest, Guanidinium)
{
    testAtype("guanidinium.sdf");
}

TEST_F (AtomtypeTest, DimethylCarbonate)
{
    testAtype("dimethyl-carbonate.pdb");
}

TEST_F (AtomtypeTest, 1EthoxyethylPhosphorylOxyethane)
{
    testAtype("1--ethoxyethylphosphoryl-oxyethane-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, 1Amino1Hydroxyguanidine)
{
    testAtype("1-amino-1-hydroxyguanidine-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, Dimethylguanidinium)
{
    testAtype("1-1-dimethylguanidinium.sdf");
}

TEST_F (AtomtypeTest, Acetate)
{
    testAtype("acetate-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, AceticAcid)
{
    testAtype("acetic-acid-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, acetone)
{
    testAtype("acetone-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, DiethylSulfate)
{
    testAtype("diethyl-sulfate-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, DisulfurMonoxide)
{
    testAtype("disulfur-monoxide-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, EthylSulfate)
{
    testAtype("ethyl-sulfate-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, Glutamate)
{
    testAtype("glutamate-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, GlutamicAcid)
{
    testAtype("glutamic-acid-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, HistidineHdHe)
{
    testAtype("histidine-hdhe-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, nitromethane)
{
    testAtype("nitromethane-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, SulfurDioxide)
{
    testAtype("sulfur-dioxide-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, ThiazylFluoride)
{
    testAtype("thiazyl-fluoride-3-oep.log.pdb");
}

TEST_F (AtomtypeTest, Water)
{
    testAtype("water-3-oep.log.pdb");
}

}

}
