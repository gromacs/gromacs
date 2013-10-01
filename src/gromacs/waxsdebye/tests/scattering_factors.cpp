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
 * Implements tests for reading scattering factors
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \author Alexander Bjorling <alexander.bjorling@chem.gu.se>
 * \author Jonas Ditz <jonas.ditz@icm.uu.se>
 * \ingroup module_waxsdebye
 */
#include "gmxpre.h"

#include "gromacs/waxsdebye/scattering_factors.h"

#include <gtest/gtest.h>

#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

namespace gmx
{

namespace test
{

class ScatteringFactorsTest : public ::testing::Test
{
    protected:

        //! Test variables
        gmx::ScatteringFactorTable                  testTable_;
        double                                      qValues[10];

        test::TestReferenceData                     refData_;
        test::TestReferenceChecker                  checker_;

        ScatteringFactorsTest( )
            : checker_(refData_.rootChecker())
        {
#ifdef GMX_DOUBLE
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-6));
#else
            checker_.setDefaultTolerance(test::relativeToleranceAsFloatingPoint(1, 1e-3));
#endif

            qValues[0] = 0.25;
            qValues[1] = 0.5;
            qValues[2] = 0.75;
            qValues[3] = 1;
            qValues[4] = 1.25;
            qValues[5] = 1.5;
            qValues[6] = 1.75;
            qValues[7] = 2;
            qValues[8] = 5;
            qValues[9] = 10;
        }

        //static init, only run once
        static void SetUpTestCase()
        {
        }

        //static destroy, only run once
        static void TearDownTestCase()
        {
        }


        //! Run the Test
        void test(const char   * filename,
                  esfType        sfType,
                  std::string    forceField,
                  std::string    reference,
                  bool           displacedSolvent,
                  unsigned int   tableSize)
        {
            testTable_.read(filename);

            EXPECT_EQ(sfType, testTable_.getSfType());
            EXPECT_EQ(forceField, testTable_.getForceField());
            EXPECT_EQ(reference, testTable_.getReference());
            EXPECT_EQ(displacedSolvent, testTable_.getDisplacedSolvent());
            EXPECT_EQ(tableSize, testTable_.size());

            int maxType = testTable_.maxType();
            for (int type = 0; type <= maxType; type++)
            {
                for (unsigned int i = 0; i < 10; i++)
                {
                    char   idString[56];
                    std::sprintf(idString, "ScatteringFactorType%dq%2.2f", type, qValues[i]);
                    double computeWithType = testTable_.computeScatteringFactor(type, qValues[i]);
                    checker_.checkReal(computeWithType, idString);
                }
            }
        }
};


TEST_F (ScatteringFactorsTest, SfactorAminoAcidDsFourierTest) {
    const char       * filename           = "sfactor_amino_acid_ds_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Niebling, Björling and Westenhoff, J. Appl. Cryst. 47 (2014) 1190";
    bool               displacedSolvent   = true;
    unsigned int       tableSize          = 21;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}

TEST_F (ScatteringFactorsTest, SfactorAminoAcidNodsFourierTest) {
    const char       * filename           = "sfactor_amino_acid_nods_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Niebling, Björling and Westenhoff, J. Appl. Cryst. 47 (2014) 1190";
    bool               displacedSolvent   = false;
    unsigned int       tableSize          = 21;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}

TEST_F (ScatteringFactorsTest, SfactorAtomDsFourierTest) {
    const char       * filename           = "sfactor_atom_ds_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Doyle and Turner, Acta Cryst. A24 (1968) 390; Fraser, MacRae and Suzuki, J. Appl. Cryst. 11 (1978) 693;  Svergun, Barberato and Koch, J Appl. Cryst. 27 (1995) 768";
    bool               displacedSolvent   = true;
    unsigned int       tableSize          = 7;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}

TEST_F (ScatteringFactorsTest, SfactorAtomNodsFourierTest) {
    const char       * filename           = "sfactor_atom_nods_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Doyle and Turner, Acta Cryst. A24 (1968) 390";
    bool               displacedSolvent   = false;
    unsigned int       tableSize          = 7;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}

TEST_F (ScatteringFactorsTest, SfactorMartiniDsFourierTest) {
    const char       * filename           = "sfactor_martini_ds_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Niebling, Björling and Westenhoff, J. Appl. Cryst. 47 (2014) 1190";
    bool               displacedSolvent   = true;
    unsigned int       tableSize          = 50;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}

TEST_F (ScatteringFactorsTest, SfactorMartiniNodsFourierTest) {
    const char       * filename           = "sfactor_martini_nods_Fourier.xml";
    esfType            sfType             = esfFourier;
    std::string        forceField         = "any";
    std::string        reference          = "Niebling, Björling and Westenhoff, J. Appl. Cryst. 47 (2014) 1190";
    bool               displacedSolvent   = false;
    unsigned int       tableSize          = 50;

    test(filename, sfType, forceField, reference, displacedSolvent, tableSize);
}


}  // namespace test

}  // namespace gmx
