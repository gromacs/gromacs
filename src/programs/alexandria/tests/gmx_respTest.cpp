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

#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/state.h"
#include "programs/alexandria/gauss_io.h"
#include "programs/alexandria/gmx_resp.h"
#include "programs/alexandria/mymol.h"
#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"


class RespTest : public ::testing::Test
{

    protected:
        gmx::test::TestReferenceData                     refData_;
        gmx::test::TestReferenceChecker                  checker_;
        static std::vector<alexandria::Resp *>           resp;
        static const int         nrAtoms = 15;


        //init sett tolecrance
        RespTest ( )
            : refData_(gmx::test::erefdataCompare), checker_(refData_.rootChecker())
        {
#ifdef GMX_DOUBLE
            checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-6));
#else
            checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-3));
#endif
        }

        // Static initiation, only run once every test.
        static void SetUpTestCase()
        {
	  std::vector<ChargeDistributionModel> iChargeDistributionModel;
	  iChargeDistributionModel.push_back(eqdAXp);
	  iChargeDistributionModel.push_back(eqdAXg);
	  //	  	  iChargeDistributionModel.push_back(eqdAXs);
          alexandria::MyMol       mp;
          alexandria::MolProp     molprop;
          gmx_atomprop_t          aps = gmx_atomprop_init();


            //needed for ReadGauss
            char                    * molnm    = (char *)"XXX";
            char                    * iupac    = (char *)"";
            char                    * conf     = (char *)"minimum";
            char                    * basis    = (char *)"";
            int                       maxpot   = 0;
            int                       nsymm    = 0;
            int                       nexcl    = 2;
            const char               *dihopt[] = { NULL, "No", "Single", "All", NULL };
            eDih                      edih     = (eDih) get_option(dihopt);

            //Needed for GenerateCharges
            real        hfac        = 0;
            real        epsr        = 1;
            const char *lot         = "B3LYP/aug-cc-pVTZ";
            char       *symm_string = (char *)"";


            ChargeGenerationAlgorithm iChargeGenerationAlgorithm = (ChargeGenerationAlgorithm) eqgRESP;

            //read input file for poldata
            std::string          dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
            alexandria::Poldata *pd       = alexandria::PoldataXml::read(dataName.c_str(), aps);

            //Read input file for molprop
            dataName = gmx::test::TestFileManager::getInputFilePath("1-butanol3-esp.log");
            ReadGauss(dataName.c_str(), molprop, molnm, iupac, conf, basis,
                      maxpot, nsymm, pd->getForceField().c_str());
            mp.molProp()->Merge(molprop);
            for (unsigned int i = 0; i < iChargeDistributionModel.size(); i++)
            {
                //Generate charges and topology
                mp.GenerateTopology(aps, pd, lot, iChargeDistributionModel[i], nexcl, false, false, edih);
                mp.gr_ = new alexandria::Resp(iChargeDistributionModel[i], mp.molProp()->getCharge());
                mp.GenerateCharges(pd, aps, iChargeDistributionModel[i], iChargeGenerationAlgorithm, hfac, epsr, lot, true, symm_string);
                resp.push_back(mp.gr_);
            }
        }

        static void TearDownTestCase()
        {
            for (unsigned int i = 0; i < resp.size(); i++)
            {
                delete resp[i];
            }
        }

};

std::vector<alexandria::Resp *> RespTest::resp;
const int                       RespTest::nrAtoms;



TEST_F (RespTest, qtotValues)
{
    std::vector<double> qtotValues(nrAtoms);
    for (unsigned int i = 0; i < resp.size(); i++)
    {
        for (int atom = 0; atom < nrAtoms; atom++)
        {
            qtotValues[atom] = resp[i]->getQtot(atom);
        }
        std::stringstream sstm;
        sstm << "qtotValuesEqdModel" << i;
        std::string       testName  = sstm.str();
        checker_.checkSequence(qtotValues.begin(), qtotValues.end(), testName.c_str());
    }
}
