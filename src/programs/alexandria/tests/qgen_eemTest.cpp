/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
#include "gromacs/topology/topology.h"
#include "programs/alexandria/gauss_io.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"
#include "programs/alexandria/qgen_eem.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

class EemTest : public gmx::test::CommandLineTestBase
{
    protected:
        gmx::test::TestReferenceChecker checker_;
        alexandria::Poldata             pd_;
        alexandria::MyMol               mp_;
        gmx_atomprop_t                  aps_;

        //init set tolecrance
        EemTest () : checker_(this->rootChecker())
        {
            alexandria::MolProp     molprop;
            aps_ = gmx_atomprop_init();

            //needed for ReadGauss
            char                    * molnm    = (char *)"XXX";
            char                    * iupac    = (char *)"";
            char                    * conf     = (char *)"minimum";
            char                    * basis    = (char *)"";
            int                       maxpot   = 0;
            int                       nsymm    = 0;

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
                      maxpot, nsymm, pd_.getForceField().c_str());
            std::vector<MolProp> vmp;
            vmp.push_back(molprop);
            mp_.molProp()->Merge(vmp.begin());
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
        }

    void testEem(ChargeDistributionModel model)
        {
            //Generate charges and topology
            const char               *lot         = "B3LYP/aug-cc-pVTZ";
            int                       nexcl       = 2;
            const char               *dihopt[]    = { NULL, "No", "Single", "All", NULL };
            eDih                      edih        = (eDih) get_option(dihopt);

            mp_.GenerateTopology(aps_, pd_, lot, model,
                                 nexcl, false, false, edih, false);

            //Needed for GenerateCharges
            real        hfac        = 0;
            real        epsr        = 1;
            char       *symm_string = (char *)"";
            t_commrec  *cr          = init_commrec();
            
            mp_.GenerateCharges(pd_, aps_, model, eqgEEM,
                                hfac, epsr, lot, true, symm_string, cr, NULL);

            std::vector<double> qtotValues;
            for(int atom = 0; atom < mp_.topology_->atoms.nr; atom++)
            {
                qtotValues.push_back(mp_.topology_->atoms.atom[atom].q);
            }
            char buf[256];
            snprintf(buf, sizeof(buf), "qtotValuesEqdAlgorithm_%d",
                     static_cast<int>(ChargeGenerationAlgorithm()));
            checker_.checkSequence(qtotValues.begin(),
                                   qtotValues.end(), buf);
        }

        static void TearDownTestCase()
        {
        }

};

TEST_F (EemTest, Bultinck)
{
    testEem(eqdBultinck);
}

TEST_F (EemTest, Rappe)
{
    testEem(eqdRappe);
}

TEST_F (EemTest, Yang)
{
    testEem(eqdYang);
}

// The tests below are outcommented since we do not have a parameters for the AX? methods.
/*
TEST_F (EemTest, AXp)
{
    testEem(eqdAXp);
}

TEST_F (EemTest, AXg)
{
    testEem(eqdAXg);
}

TEST_F (EemTest, AXs)
{
    testEem(eqdAXs);
}
*/
