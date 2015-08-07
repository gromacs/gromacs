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
        int mol = 1;
        for(alexandria::MolPropIterator mpi=mp_.begin(); (mpi<mp_.end()); ++mpi, ++mol)
        {
            char mbuf[256];
            snprintf(mbuf, sizeof(mbuf), "molecule %d name", mol);
            checker_.checkString(mpi->getMolname(), mbuf);
            mpi->GenerateFormula(aps_);
            snprintf(mbuf, sizeof(mbuf), "molecule %d formula", mol);
            checker_.checkString(mpi->getFormula(), mbuf);
            snprintf(mbuf, sizeof(mbuf), "molecule %d number of bonds", mol);
            checker_.checkInteger(mpi->NBond(), mbuf);
            int i = 1;
            for(alexandria::BondIterator bi=mpi->BeginBond(); (bi<mpi->EndBond()); ++bi)
            {
                char buf[256];
                snprintf(buf, sizeof(buf), "atoms %d %d order %d", 
                         bi->getAi(), bi->getAj(), bi->getBondOrder());
                std::string bond("bond");
                char ibuf[256];
                snprintf(ibuf, sizeof(ibuf), "molecule %d bond %d", mol, i++);
                checker_.checkString(buf, ibuf);
            }
        }
    }
    
    void testExperiments ()
    {  
        int mol = 1;
        for(alexandria::MolPropIterator mpi=mp_.begin(); (mpi<mp_.end()); ++mpi, ++mol)
        {
            char mbuf[256];
            int exp = 1;
            snprintf(mbuf, sizeof(mbuf), "molecule %d number of experiments", mol);
            checker_.checkInteger(mpi->NExperiment(), mbuf);
            for(alexandria::ExperimentIterator expi=mpi->BeginExperiment(); (expi<mpi->EndExperiment()); ++expi, ++exp)
            {
                char cbuf[256];
                snprintf(cbuf, sizeof(cbuf), "molecule %d exper %d", mol, exp);
                int nener = 1;
                for(alexandria::MolecularEnergyIterator ei=expi->BeginEnergy(); (ei<expi->EndEnergy()); ++ei)
                {
                    char ebuf[256];
                    snprintf(mbuf, sizeof(mbuf), "%s energy %d", cbuf, nener++);
                    snprintf(ebuf, sizeof(ebuf), "%s %g +/- %g %s", 
                             ei->getType().c_str(), 
                             ei->getValue(), ei->getError(),
                             ei->getUnit().c_str());
                    checker_.checkString(ebuf, mbuf);
                }
            }
        }
    }
    
    void testCalculations ()
    {  
        int mol = 1;
        for(alexandria::MolPropIterator mpi=mp_.begin(); (mpi<mp_.end()); ++mpi, ++mol)
        {
            char mbuf[256];
            int calc = 1;
            snprintf(mbuf, sizeof(mbuf), "molecule %d number of calcs", mol);
            checker_.checkInteger(mpi->NCalculation(), mbuf);
            for(alexandria::CalculationIterator ci=mpi->BeginCalculation(); (ci<mpi->EndCalculation()); ++ci, ++calc)
            {
                char cbuf[256];
                snprintf(cbuf, sizeof(cbuf), "molecule %d cakc %d", mol, calc);
                snprintf(mbuf, sizeof(mbuf), "%s program", cbuf);
                checker_.checkString(ci->getProgram(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "%s basisset", cbuf);
                checker_.checkString(ci->getBasisset(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "%s method", cbuf);
                checker_.checkString(ci->getMethod(), mbuf);
                snprintf(mbuf, sizeof(mbuf), "%s number of polar", cbuf);
                checker_.checkInteger(ci->NPolar(), mbuf);
                for(alexandria::MolecularPolarizabilityIterator poli=ci->BeginPolar(); (poli<ci->EndPolar()); ++ poli)
                {
                    snprintf(mbuf, sizeof(mbuf), "%s polar XX", cbuf);
                    checker_.checkDouble(poli->getXX(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s polar YY", cbuf);
                    checker_.checkDouble(poli->getYY(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s polar ZZ", cbuf);
                    checker_.checkDouble(poli->getZZ(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s polar XY", cbuf);
                    checker_.checkDouble(poli->getXY(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s polar XZ", cbuf);
                    checker_.checkDouble(poli->getXZ(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s polar YZ", cbuf);
                    checker_.checkDouble(poli->getYZ(), mbuf);
                }
                
                for(alexandria::MolecularQuadrupoleIterator qi=ci->BeginQuadrupole(); (qi<ci->EndQuadrupole()); ++ qi)
                {
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole XX", cbuf);
                    checker_.checkDouble(qi->getXX(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole YY", cbuf);
                    checker_.checkDouble(qi->getYY(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole ZZ", cbuf);
                    checker_.checkDouble(qi->getZZ(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole XY", cbuf);
                    checker_.checkDouble(qi->getXY(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole XZ", cbuf);
                    checker_.checkDouble(qi->getXZ(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s quadrupole YZ", cbuf);
                    checker_.checkDouble(qi->getYZ(), mbuf);
                }

                for(alexandria::MolecularDipoleIterator dip=ci->BeginDipole(); (dip<ci->EndDipole()); ++ dip)
                {
                    snprintf(mbuf, sizeof(mbuf), "%s dipole X", cbuf);
                    checker_.checkDouble(dip->getX(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s dipole Y", cbuf);
                    checker_.checkDouble(dip->getY(), mbuf);
                    snprintf(mbuf, sizeof(mbuf), "%s dipole Z", cbuf);
                    checker_.checkDouble(dip->getZ(), mbuf);
                }
                
                int nener = 1;
                for(alexandria::MolecularEnergyIterator ei=ci->BeginEnergy(); (ei<ci->EndEnergy()); ++ei)
                {
                    char ebuf[256];
                    snprintf(mbuf, sizeof(mbuf), "%s energy %d", cbuf, nener++);
                    snprintf(ebuf, sizeof(ebuf), "%s %g +/- %g %s", 
                             ei->getType().c_str(), 
                             ei->getValue(), ei->getError(),
                             ei->getUnit().c_str());
                    checker_.checkString(ebuf, mbuf);
                }
            }
        }
    }
    
    static void TearDownTestCase()
    {
    }
   
};

TEST_F (MolpropTest, NameFormulaBonds){
    testMolProp();
}

TEST_F (MolpropTest, Calculations){
    testCalculations();
}

TEST_F (MolpropTest, Experiments){
    testExperiments();
}

