/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016, by the GROMACS development team, led by
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
 * Implements test of autocorrelation function routines
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_correlationfunctions
 */


#include <math.h>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
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

            //needed for ReadGauss
            const char *molnm    = (char *)"XXX";
            const char *iupac    = (char *)"";
            const char *conf     = (char *)"minimum";
            const char *basis    = (char *)"";
            const char *jobtype  = (char *)"Pop";
            int         maxpot   = 0;
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

            auto tolerance = gmx::test::relativeToleranceAsFloatingPoint(1.0, 1e-5);
            checker_.setDefaultTolerance(tolerance);
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

        void testResp(ChargeDistributionModel model, bool bPolar)
        {
            //Generate charges and topology
            const char               *lot         = "B3LYP/aug-cc-pVTZ";

            gmx::MDModules            mdModules;
            t_inputrec               *inputrec   = mdModules.inputrec();
            fill_inputrec(inputrec);
            mp_.setInputrec(inputrec);
            mp_.GenerateTopology(aps_, pd_, lot, model,
                                 false, false, false, bPolar);
            //Needed for GenerateCharges
            real           hfac        = 0;
            real           epsr        = 1;
            char          *symm_string = (char *)"";
            t_commrec     *cr          = init_commrec();
            gmx::MDLogger  mdlog       = getMdLogger(cr, stdout);

            mp_.GenerateCharges(pd_, mdlog, aps_, model, eqgESP,
                                hfac, epsr, lot, false, symm_string, cr, NULL,
                                as_rvec_array(mp_.x_->data()));

            std::vector<double> qtotValues;
            for (int atom = 0; atom < mp_.mtop_->moltype[0].atoms.nr; atom++)
            {
                qtotValues.push_back(mp_.mtop_->moltype[0].atoms.atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdModel_%d",
                     static_cast<int>(model));
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
    testResp(eqdAXg, true);
}

TEST_F (RespTest, AXsValues)
{
    testResp(eqdAXs, false);
}

TEST_F (RespTest, AXsPolarValues)
{
    testResp(eqdAXs, true);
}
