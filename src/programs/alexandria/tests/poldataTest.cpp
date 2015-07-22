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
//#include "gmxpre.h"

#include <math.h>
#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "programs/alexandria/poldata.h"
#include "programs/alexandria/poldata_xml.h"

class PoldataTest : public ::testing::Test
{


protected:
  static  alexandria::Poldata * pd;
  gmx::test::TestReferenceData                     refData_;
  gmx::test::TestReferenceChecker                  checker_;
  static const  int numModels = 3;

  //init sett tolecrance
  PoldataTest ( )
    :refData_(gmx::test::erefdataUpdateAll), checker_(refData_.rootChecker())
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
    gmx_atomprop_t aps = gmx_atomprop_init();

    // Reads the file, the file only suport 3 chargedisributionModels
    // eqdAXp,eqdAXg,  eqdAXs,
    std::string dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
    pd = alexandria::PoldataXml::read(dataName.c_str(), aps);  

	  
  }

  static void TearDownTestCase()
  {
  }

   
};


alexandria::Poldata * PoldataTest::pd;
const int PoldataTest::numModels;


TEST_F (PoldataTest, getAtype){

  std::string elem;
  std::string desc;
  std::string name;
  std::string ptype;
  std::string btype;
  std::string vdwparams;
  double                   J0;
  pd->getAtype(&elem,
	       &desc,
	       &name,
	       &ptype,
	       &btype,
	       &vdwparams,
	       &J0);
  checker_.checkString(elem,"elem");
  checker_.checkString(desc,"desc");
  checker_.checkString(name,"name");
  checker_.checkString(ptype,"ptype");
  checker_.checkString(btype,"btype");
  checker_.checkString(vdwparams,"vdwparams");
  checker_.checkDouble(J0,"J0");

}

TEST_F (PoldataTest, chi)
{
  
  
  std::string                    name;
  double                   J0;
  std::vector<double>      values;
  pd->getAtype(NULL,
	       NULL,
	       &name,
	       NULL,
	       NULL,
	       NULL,
	       &J0);
  for (int model = 0; model < numModels; model++)
    {
      values.push_back(pd->getChi0((ChargeDistributionModel)model, name));
    }
  checker_.checkSequence(values.begin(), values.end(), "chi");
}

TEST_F (PoldataTest, row){
  std::string              name;
  double                   J0;
  std::vector<double>      values;
  
  
  pd->getAtype(NULL,
	       NULL,
	       &name,
	       NULL,
	       NULL,
	       NULL,
	       &J0);
  int numAtoms = 3;
  

  for (int atomNr = 0; atomNr < numAtoms; atomNr++)
    {
      for (int model = 0; model <  numModels; model++)
	{
	  values.push_back(pd->getRow((ChargeDistributionModel)model, name, atomNr));
        }
    }
  checker_.checkSequence(values.begin(), values.end(), "row");
}


TEST_F (PoldataTest, zeta)
{
  std::string              name;
  double                   J0;
  std::vector<double>      values;
  pd->getAtype(NULL,
	       NULL,
	       &name,
	       NULL,
	       NULL,
	       NULL,
	       &J0);
  int numAtoms = 3;
  for (int atomNr = 0; atomNr < numAtoms; atomNr++)
    {
      for (int model = 0; model <  numModels; model++)
	{
	  values.push_back(pd->getZeta((ChargeDistributionModel)model, name, atomNr));
        }
    }
  checker_.checkSequence(values.begin(), values.end(), "zeta");
}




TEST_F (PoldataTest, geters)
{
  std::string value = pd->getPolarUnit( );
  checker_.checkString(value,"polarUnit");

  value =  pd->getPolarRef( );
  checker_.checkString(value, "polarRef");

  value =  pd->getForceField( );
  checker_.checkString(value,"forceFiled");

  value =  pd->getLengthUnit( );
  checker_.checkString(value,"lenghtUnit");
}
