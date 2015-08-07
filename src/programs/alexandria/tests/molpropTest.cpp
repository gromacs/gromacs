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
 * Implements test of molprop functionality
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup alexandria
 */
#include "gmxpre.h"

#include <math.h>

#include <cstdlib>

#include <gtest/gtest.h>

#include "gromacs/topology/atomprop.h"
#include "gromacs/utility/snprintf.h"

#include "testutils/cmdlinetest.h"
#include "testutils/refdata.h"
#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"

#include "programs/alexandria/molprop.h"
#include "programs/alexandria/molprop_xml.h"
#include "programs/alexandria/poldata_xml.h"

class MolpropTest : public gmx::test::CommandLineTestBase
{


protected:
    std::vector<alexandria::MolProp>  mp_;
    gmx_atomprop_t                    aps_;
    gmx::test::TestReferenceData      refData_;
    gmx::test::TestReferenceChecker   checker_;

    // Init set tolerance
    MolpropTest()
        :refData_(gmx::test::erefdataUpdateAll), checker_(refData_.rootChecker())
    {
#ifdef GMX_DOUBLE
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-6));
#else
        checker_.setDefaultTolerance(gmx::test::relativeToleranceAsFloatingPoint(1, 1e-3));
#endif
        aps_  = gmx_atomprop_init();

        std::string mpFile = fileManager().getInputFilePath("molprop.dat");
        MolPropRead(mpFile.c_str(), mp_);
    }

    // Static initiation, only run once every test.
    static void SetUpTestCase()
    {	  
    }
  
    void testMolProp ()
    {  
        for(alexandria::MolPropIterator mpi=mp_.begin(); (mpi<mp_.end()); ++mpi)
        {
            checker_.checkString(mpi->getMolname(), "molecule name");
            mpi->GenerateFormula(aps_);
            checker_.checkString(mpi->getFormula(), "formula");
            checker_.checkInteger(mpi->NBond(), "number of bonds");
            int i = 1;
            for(alexandria::BondIterator bi=mpi->BeginBond(); (bi<mpi->EndBond()); ++bi)
            {
                char buf[256];
                snprintf(buf, sizeof(buf), "%d %d %d", 
                         bi->getAi(), bi->getAj(), bi->getBondOrder());
                std::string bond("bond");
                char ibuf[256];
                snprintf(ibuf, sizeof(ibuf), "Bond %d", i++);
                checker_.checkString(buf, ibuf);
            }
        }
    }
    void testCalculations ()
    {  
        for(alexandria::MolPropIterator mpi=mp_.begin(); (mpi<mp_.end()); ++mpi)
        {
            for(alexandria::CalculationIterator ci=mpi->BeginCalculation(); (ci<mpi->EndCalculation()); ++ci)
            {
                checker_.checkString(ci->getProgram(), "program");
                checker_.checkString(ci->getBasisset(), "basisset");
                checker_.checkString(ci->getMethod(), "method");
                for(alexandria::MolecularPolarizabilityIterator poli=ci->BeginPolar(); (poli<ci->EndPolar()); ++ poli)
                {
                    checker_.checkDouble(poli->getXX(), "Polarizability XX");
                    checker_.checkDouble(poli->getYY(), "Polarizability YY");
                    checker_.checkDouble(poli->getZZ(), "Polarizability ZZ");
                    checker_.checkDouble(poli->getXY(), "Polarizability XY");
                    checker_.checkDouble(poli->getXZ(), "Polarizability XZ");
                    checker_.checkDouble(poli->getYZ(), "Polarizability YZ");
                }

                for(alexandria::MolecularDipoleIterator dip=ci->BeginDipole(); (dip<ci->EndDipole()); ++ dip)
                {
                    checker_.checkDouble(dip->getX(), "Dipole X");
                    checker_.checkDouble(dip->getY(), "Dipole Y");
                    checker_.checkDouble(dip->getZ(), "Dipole Z");
                }
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
   
};

TEST_F (MolpropTest, getMolnameFormula){
    testMolProp();
}

TEST_F (MolpropTest, Calculations){
    testCalculations();
}

