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


#include <cmath>
#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/mdrunutility/mdmodules.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/logger.h"
#include "programs/alexandria/fill_inputrec.h"
#include "programs/alexandria/gauss_io.h"
#include "programs/alexandria/getmdlogger.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"
#include "programs/alexandria/qgen_resp.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

class RespTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::Poldata             pd_;
        alexandria::MyMol               mp_;
        gmx_atomprop_t                  aps_;

        //init set tolecrance
        RespTest () : checker_(this->rootChecker())
        {
            alexandria::MolProp     molprop;
            aps_ = gmx_atomprop_init();

            setenv("GMX_NB_GENERIC", "1", 1);
            //needed for ReadGauss
            const char *molnm    = (char *)"XXX";
            const char *iupac    = (char *)"";
            const char *conf     = (char *)"minimum";
            const char *basis    = (char *)"";
            const char *jobtype  = (char *)"Pop";
            int         maxpot   = 100;
            int         nsymm    = 0;

            //read input file for poldata
            std::string dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
            try
            {
                alexandria::readPoldata(dataName.c_str(), pd_, aps_);
            }
            GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

            //Read input file for molprop
            dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol3-esp.log");
            ReadGauss(dataName.c_str(), molprop, molnm, iupac, conf, basis,
                      maxpot, nsymm, pd_.getForceField().c_str(), jobtype);
            std::vector<MolProp> vmp;
            vmp.push_back(molprop);
            mp_.molProp()->Merge(vmp.begin());

            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 5e-2);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testResp(ChargeDistributionModel qdist, bool bPolar)
        {
            //Generate charges and topology
            const char               *lot        = "B3LYP/aug-cc-pVTZ";

            gmx::MDModules            mdModules;
            t_inputrec               *inputrec   = mdModules.inputrec();
            fill_inputrec(inputrec);
            mp_.setInputrec(inputrec);
            mp_.GenerateTopology(aps_, pd_, lot, qdist, false, false, false, bPolar, false, nullptr);
            
            //Needed for GenerateCharges
            real           hfac        = 0;
            real           watoms      = 0;
            char          *symm_string = (char *)"";
            t_commrec     *cr          = init_commrec();
            gmx::MDLogger  mdlog       = getMdLogger(cr, stdout);
            gmx_hw_info_t *hwinfo      = gmx_detect_hardware(mdlog, cr, false);
            int            qcycle      = 1;
            real           qtol        = 1e-3;

            if(!bPolar)
            {
                mp_.GenerateCharges(pd_, mdlog, aps_, qdist, eqgESP, watoms,
                                    hfac, lot, false, symm_string, cr, nullptr, hwinfo, qcycle, qtol, nullptr);
            }
            else
            {
                if (qdist == eqdAXpg)
                {
                    mp_.GenerateCharges(pd_, mdlog, aps_, qdist, eqgESP, watoms,
                                        hfac, lot, false, symm_string, cr, nullptr, hwinfo, qcycle, qtol, nullptr);
                }
                else if (qdist == eqdAXps)
                {
                    inputrec->coulombtype = eelUSER;
                    mp_.setInputrec(inputrec);
                    std::string tabFile = fileManager().getInputFilePath("table.xvg");
                    mp_.GenerateCharges(pd_, mdlog, aps_, qdist, eqgESP, watoms,
                                        hfac, lot, false, symm_string, cr,
                                        tabFile.c_str(), hwinfo, qcycle, qtol, nullptr);
                }
            }
            std::vector<double> qtotValues;
            for (int atom = 0; atom < mp_.mtop_->moltype[0].atoms.nr; atom++)
            {
                qtotValues.push_back(mp_.mtop_->moltype[0].atoms.atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdModel_%d",
                     static_cast<int>(qdist));
            checker_.checkSequence(qtotValues.begin(),
                                   qtotValues.end(), buf);
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (RespTest, AXpValues)
{
    testResp(eqdAXp, false);
}

TEST_F (RespTest, AXgValues)
{
    testResp(eqdAXg, false);
}

TEST_F (RespTest, AXgPolarValues)
{
    testResp(eqdAXpg, true);
}

TEST_F (RespTest, AXsValues)
{
    testResp(eqdAXs, false);
}

TEST_F (RespTest, AXsPolarValues)
{
    testResp(eqdAXps, true);
}
