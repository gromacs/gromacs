/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2018
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


#include <math.h>
#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/physicalnodecommunicator.h"
#include "programs/alexandria/fill_inputrec.h"
#include "programs/alexandria/babel_io.h"
#include "programs/alexandria/getmdlogger.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"
#include "programs/alexandria/qgen_acm.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace alexandria
{

namespace
{

enum informat{
    einfLOG = 0,
    einfPDB = 1, 
    einfNR
};

class AcmTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::Poldata             pd_;
        alexandria::MyMol               mp_;
        gmx_atomprop_t                  aps_;

        //init set tolerance
        AcmTest () : checker_(this->rootChecker())
        {
            aps_ = gmx_atomprop_init();

            //read input file for poldata
            std::string dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
            try
            {
                alexandria::readPoldata(dataName.c_str(), pd_, aps_);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;


            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testAcm(ChargeDistributionModel model, informat inputformat)
        {                       
            int        maxpot    = 100; 
            int        nsymm     = 0;
            const char *lot      = nullptr;
            const char *dihopt[] = { nullptr, "No", "Single", "All", nullptr };            
            const char *molnm    = (char *)"1-butanol";
            const char *iupac    = (char *)"1-butanol";
            const char *conf     = (char *)"minimum";
            const char *basis    = (char *)"";
            const char *jobtype  = (char *)"Opt";
            
            std::string           dataName;
            alexandria::MolProp   molprop;
            std::vector<MolProp>  vmp;

            if (inputformat == einfLOG)
            {
                lot      = "B3LYP/Gen";
                dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol-3-oep.log");
            }
            else
            {
                lot      = "AFF/ACM";
                dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol.pdb");
            }
            
            readBabel(dataName.c_str(), 
                      molprop, 
                      molnm, 
                      iupac, 
                      conf, 
                      basis,
                      maxpot, 
                      nsymm, 
                      pd_.getForceField().c_str(), 
                      jobtype);
            
            vmp.push_back(molprop);
            mp_.molProp()->Merge(vmp.begin());
            
            
            // Generate charges and topology
            
            eDih            edih       = (eDih) get_option(dihopt);
            t_inputrec      inputrecInstance;
            t_inputrec     *inputrec   = &inputrecInstance;
            fill_inputrec(inputrec);
            mp_.setInputrec(inputrec);

            mp_.GenerateTopology(aps_, pd_, lot, model, false, false, edih, false, nullptr);

            // Needed for GenerateCharges
            real           hfac                  = 0;
            real           watoms                = 0;
            char          *symm_string           = (char *)"";
            t_commrec     *cr                    = init_commrec();
            auto           pnc                   = gmx::PhysicalNodeCommunicator(MPI_COMM_WORLD, 0);
            auto           mdlog                 = getMdLogger(cr, stdout);
            auto           hwinfo                = gmx_detect_hardware(mdlog, pnc);
            int            qcycle                = 1;
            real           qtol                  = 1e-3;

            mp_.GenerateCharges(pd_, mdlog, aps_, model,
                                eqgACM, watoms, hfac, lot,
                                true, symm_string, cr,
                                nullptr, hwinfo, qcycle,
                                maxpot, qtol, nullptr);

            std::vector<double> qtotValues;
            for (int atom = 0; atom < mp_.topology_->atoms.nr; atom++)
            {
                qtotValues.push_back(mp_.topology_->atoms.atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdAlgorithm_%d", static_cast<int>(ChargeGenerationAlgorithm()));
            checker_.checkInteger(static_cast<int>(qtotValues.size()), "qtotSize");
            checker_.checkSequence(qtotValues.begin(), qtotValues.end(), buf);
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (AcmTest, BultinckLog)
{
    testAcm(eqdBultinck, einfLOG);
}

TEST_F (AcmTest, BultinckPDB)
{
    testAcm(eqdBultinck, einfPDB);
}

TEST_F (AcmTest, RappeLog)
{
    testAcm(eqdRappe, einfLOG);
}

TEST_F (AcmTest, RappePDB)
{
    testAcm(eqdRappe, einfPDB);
}

TEST_F (AcmTest, YangLog)
{
    testAcm(eqdYang, einfLOG);
}

TEST_F (AcmTest, YangPDB)
{
    testAcm(eqdYang, einfPDB);
}

TEST_F (AcmTest, AXpgLOG)
{
    testAcm(eqdAXpg, einfLOG);
}

TEST_F (AcmTest, AXpgPDB)
{
    testAcm(eqdAXpg, einfPDB);
}

TEST_F (AcmTest, AXpsLOG)
{
    testAcm(eqdAXps, einfLOG);
}

TEST_F (AcmTest, AXpsPDB)
{
    testAcm(eqdAXps, einfPDB);
}
 

}

}
