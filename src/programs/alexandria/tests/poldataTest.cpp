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

  //using gmx::test::TestReferenceData;
protected:
  static  alexandria::Poldata * pd;
  gmx::test::TestReferenceData                     refData_;
  gmx::test::TestReferenceChecker                  checker_;

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
    std::string dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
	  
    gmx_atomprop_t aps = gmx_atomprop_init();
    //pd = NULL;
    pd = alexandria::PoldataXml::read(dataName.c_str(), aps);  
    std::cout << "e i test" << "\n";         
	  
  }


  void setUp(){


  }
  static void TearDownTestCase()
  {
  }


  void test()
  {
    std::cout << "e i test";
    //ChargeDistributionModel  eqg_model;
    char                    *name;
    double                   J0;
    std::vector<double>      values;
    pd->get_atype(NULL,
		  NULL,
		  &name,
		  NULL,
		  NULL,
		  NULL,
		  &J0);
    //	  pd->get_eemprops( &eqg_model, &name, &J0, &chi0, NULL, NULL, NULL);
    for (int model = 0; model < (int) 2; model++)
      {
	values.push_back(pd->get_chi0((ChargeDistributionModel)model, name));
      }
    checker_.checkSequence(values.begin(), values.end(), "chi");
  }
};


alexandria::Poldata * PoldataTest::pd;


TEST_F (PoldataTest, chi)
{
  test();
}

TEST_F (PoldataTest, row){
  char                    *name;
  double                   J0;
  std::vector<double>      values;
  
  
  pd->get_atype(NULL,
		NULL,
		&name,
		NULL,
		NULL,
		NULL,
		&J0);
  int nrow =3;
  // pd->get_nzeta((ChargeDistributionModel)model, name);

  for (int atomNr = 0; atomNr < std::min(3, nrow); atomNr++)
    {
      for (int model = 0; model < (int) 2; model++)
	{
	  values.push_back(pd->get_row((ChargeDistributionModel)model, name, atomNr));
        }
    }
  checker_.checkSequence(values.begin(), values.end(), "row");
}


TEST_F (PoldataTest, zeta)
{
  char                    *name;
  double                   J0;
  std::vector<double>      values;
  pd->get_atype(NULL,
		NULL,
		&name,
		NULL,
		NULL,
		NULL,
		&J0);
  int nzeta =3;
  // pd->get_nzeta((ChargeDistributionModel)model, name);

  for (int atomNr = 0; atomNr < std::min(3, nzeta); atomNr++)
    {
      for (int model = 0; model < (int) 2; model++)
	{
	  values.push_back(pd->get_zeta((ChargeDistributionModel)model, name, atomNr));
        }
    }
  checker_.checkSequence(values.begin(), values.end(), "zeta");
}

    char * get_polar_unit( );

    char * get_polar_ref( );

    char * get_force_field( );

    char * get_length_unit( );



TEST_F (PoldataTest, geters)
{

  /*  Char * value = NULL;
  value = pd->get_polar_unit( );
  checker_.checkSequence(values.begin(), values.end(), "zeta");

  value = pd->get_polar_ref( );


  value = get_force_field();


  value = get_length_unit();*/

}
