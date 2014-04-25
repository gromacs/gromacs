/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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

#include <cmath>
#include <sstream>

#include <gtest/gtest.h>
#include <gromacs/fileio/xvgr.h>
#include <testutils/refdata.h>
#include <testutils/testasserts.h>
#include <testutils/testfilemanager.h>

#include "gromacs/utility/smalloc.h"
#include "gromacs/energyanalysis/fluctprops.h"
#include "gromacs/energyanalysis/viscosity.h"
#include "gromacs/energyanalysis/freeenergydifference.h"
#include "gromacs/energyanalysis/dhdl.h"
#include "gromacs/energyanalysis/handler.h"

class EnergyanalysisTest : public ::testing::Test
{
    protected:

        EnergyanalysisTest( )
        {
        }

        //static init, only runed onecs every test.
        static void SetUpTestCase()
        {

        }

        void SetUp( )
        {
            //Seting up ref test
        }

        static void TearDownTestCase()
        {

        }

        void TearDown( )
        {
        }


        void test(gmx::EnergyAnalysisPointer energyanalysis, char * option[], int nrOptions, const std::string testName)
        {
            gmx::EnergyHandler handler;
            //run the test
            handler.addAnalysisTool(energyanalysis);
            handler.prepare(&nrOptions, option);
            handler.readFiles();
            testResult(testName);
        }





        void testResult(std::string testName)
        {
            double ** values;
            int       nrColums;

            double    percentBound = 0.01;
            gmx::test::TestReferenceData    data(gmx::test::erefdataCreateMissing);
            gmx::test::TestReferenceChecker checker(data.rootChecker());

            std::string                     fileName = gmx::test::TestFileManager::getInputFilePath((testName + ".xvg").c_str());
            //read result file
            int nrLines =  read_xvg(fileName.c_str(), &values, &nrColums );

            for (int colum = 1; colum  < nrColums; colum++)
            {
                for (int row = 0; row < nrLines; row++)
                {
                    checker.setDefaultTolerance(gmx::test::FloatingPointTolerance(std::abs(values[colum][row]) * percentBound, -1, false));
                    std::string tempName = testName  + convertInt(colum) + convertInt(row);
                    checker.checkReal(values[colum][row], tempName.c_str());
                }
            }
            //alocated in read_xvg
            for (int i = 0; i < nrColums; i++)
            {
                sfree(values[i]);
                values[i] = NULL;
            }
            sfree(values);
            values = NULL;
        }



        //std::to_string does not compile in win
        std::string convertInt(int number)
        {
            std::stringstream ss; //create a stringstream
            ss << number;         //add number to the stream
            return ss.str();      //return a string with the contents of the stream
        }
};



TEST_F (EnergyanalysisTest, fluctProps) {
    int                        nrOptions    = 3;
    int                        curentOption = 1;
    const std::string          testName     = "fluctprops";
    char                     * option[nrOptions*2+1];

    gmx::FluctPropsInfo        info;
    gmx::EnergyAnalysisPointer energyanalysis = info.create();

    std::string                option1 = gmx::test::TestFileManager::getInputFilePath("ener.edr");
    option[curentOption++] = (char *)"-f";
    option[curentOption++] = const_cast<char *>(option1.c_str());

    std::string option2 = gmx::test::TestFileManager::getInputFilePath((testName + ".xvg").c_str());
    option[curentOption++] = (char *)"-convergence";
    option[curentOption++] = const_cast<char *>(option2.c_str());

    std::string option3 = "10";
    option[curentOption++] = (char *)"-nmol";
    option[curentOption++] = const_cast<char *>(option3.c_str());

    test(energyanalysis, option, curentOption, testName);
}


TEST_F (EnergyanalysisTest, viscosity) {
    int                        nrOptions    = 2;
    int                        curentOption = 1;
    std::string                testName     = "viscosity";
    char                     * option[nrOptions*2+1];

    gmx::ViscosityInfo         info;
    gmx::EnergyAnalysisPointer energyanalysis = info.create();

    std::string                option1 = gmx::test::TestFileManager::getInputFilePath("ener.edr");
    option[curentOption++] = (char *)"-f";
    option[curentOption++] = const_cast<char *>(option1.c_str());

    std::string option2 = gmx::test::TestFileManager::getInputFilePath((testName + ".xvg").c_str());
    option[curentOption++] = (char *)"-vis"; //option does not work atm Anders Gärdenäs 13-08-2014
    option[curentOption++] = const_cast<char *>(option2.c_str());

    test(energyanalysis, option, curentOption, testName);
}


//does not calculate might need one more input file
TEST_F (EnergyanalysisTest, dhdl) {
    int                           nrOptions    = 2;
    int                           curentOption = 1;
    std::string                   testName     = "dhdl";
    char                        * option[nrOptions*2+1];

    gmx::FreeEnergyDifferenceInfo info;
    gmx::EnergyAnalysisPointer    energyanalysis = info.create();


    std::string option1 = gmx::test::TestFileManager::getInputFilePath("ener.edr");
    option[curentOption++] = (char *)"-f";
    option[curentOption++] = const_cast<char *>(option1.c_str());

    std::string option2 = gmx::test::TestFileManager::getInputFilePath((testName + ".xvg").c_str());
    option[curentOption++] = (char *)"-outDhdl";
    option[curentOption++] = const_cast<char *>(option2.c_str());

    test(energyanalysis, option, nrOptions, testName);
}




//invalid no output file
TEST_F (EnergyanalysisTest, freeenergydifference) {
    int                           nrOptions    = 1;
    int                           curentOption = 1;
    std::string                   testName     = "freeenergydifference";
    char                        * option[nrOptions*2+1];

    gmx::FreeEnergyDifferenceInfo info;
    gmx::EnergyAnalysisPointer    energyanalysis = info.create();


    std::string option1 = gmx::test::TestFileManager::getInputFilePath("ener.edr");
    option[curentOption++] = (char *)"-f";
    option[curentOption++] = const_cast<char *>(option1.c_str());

    test(energyanalysis, option, curentOption, testName);
}
