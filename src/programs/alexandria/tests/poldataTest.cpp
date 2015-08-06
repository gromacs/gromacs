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
        static  alexandria::Poldata                    * pd;
        gmx::test::TestReferenceData                     refData_;
        gmx::test::TestReferenceChecker                  checker_;
        static const  int numModels = 3;
        static   std::vector<std::string>                atomNames;
        static std::string atomName;

        //init sett tolecrance
        PoldataTest ( )
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
            gmx_atomprop_t aps = gmx_atomprop_init();

            // Reads the file, the file only suport 3 chargedisributionModels
            // eqdAXp,eqdAXg,  eqdAXs,  23/07/15
            std::string dataName = gmx::test::TestFileManager::getInputFilePath("gentop.dat");
            pd = alexandria::PoldataXml::read(dataName.c_str(), aps);

            alexandria::FfatypeIterator iter = pd->getAtypeBegin();
            atomName = iter->type;
            for (; iter != pd->getAtypeEnd(); iter++)
            {
                atomNames.push_back(iter->type);
            }

        }

        static void TearDownTestCase()
        {
        }


};


alexandria::Poldata    * PoldataTest::pd;
const int                PoldataTest::numModels;
std::vector<std::string> PoldataTest::atomNames;
std::string              PoldataTest::atomName;

TEST_F (PoldataTest, getAtype){
    alexandria::FfatypeIterator aType =  pd->getAtypeBegin();

    checker_.checkString(aType->elem, "elem");
    checker_.checkString(aType->desc, "desc");
    checker_.checkString(aType->type, "type");
    checker_.checkString(aType->ptype, "ptype");
    checker_.checkString(aType->btype, "btype");
    checker_.checkString(aType->vdwparams, "vdwparams");
    checker_.checkDouble(aType->refEnthalpy, "refEnthalpy");
}


TEST_F(PoldataTest, searchAtype){
    std::string        elem;
    std::string        desc;
    std::string        atype;
    std::string        ptype;
    std::string        btype;
    std::string        vdwparams;
    pd->searchAtype(atomName,
                    &elem,
                    &desc,
                    &atype,
                    &ptype,
                    &btype,
                    &vdwparams);

    checker_.checkString(elem, "elem");
    checker_.checkString(desc, "desc");
    checker_.checkString(ptype, "ptype");
    checker_.checkString(btype, "btype");
    checker_.checkString(vdwparams, "vdwparams");
    assert(atomName.compare(atype) == 0);
}

TEST_F(PoldataTest, addAtype){
    const std::string        elem         = "elm";
    const std::string        desc         = "temproary test atom";
    const std::string        atype        = "aType";
    const std::string        ptype        = "Type";
    const std::string        btype        = "bType";
    const std::string        vdwparams    = "vdwparams";
    const double             ref_enthalpy = 1000;

    std::string              newElem;
    std::string              newDesc;
    std::string              newAtype;
    std::string              newPtype;
    std::string              newBtype;
    std::string              newVdwparams;

    pd->addAtype( elem,
                  desc,
                  atype,
                  ptype,
                  btype,
                  vdwparams,
                  ref_enthalpy);

    int valid = pd->searchAtype(atype,
                                &newElem,
                                &newDesc,
                                &newAtype,
                                &newPtype,
                                &newBtype,
                                &newVdwparams);
    //will faill if the seartch failed
    assert(valid == 1);

    //test if the extraction where corect
    assert(newElem.compare(elem) == 0);
    assert(newDesc.compare(desc) == 0);
    assert(newAtype.compare(atype) == 0);
    assert(newPtype.compare(ptype) == 0);
    assert(newBtype.compare(btype) == 0);
    assert(newVdwparams.compare(vdwparams) == 0);
}



TEST_F (PoldataTest, Ptype)
{
    alexandria::PtypeIterator ptype = pd->getPtypeBegin();
    checker_.checkString(ptype->type, "type");
    checker_.checkString(ptype->miller, "miller");
    checker_.checkString(ptype->bosque, "bosque");
    checker_.checkDouble(ptype->polarizability, "polarizability");
    checker_.checkDouble(ptype->sigPol, "sigPol");
}

TEST_F (PoldataTest, Miller)
{
    alexandria::MillerIterator miller = pd->getMillerBegin();
    checker_.checkInteger(miller->atomnumber, "atomnumber");
    checker_.checkDouble(miller->tauAhc, "tauAhc");
    checker_.checkDouble(miller->alphaAhp, "alphaAhp");
}


TEST_F (PoldataTest, Bosque)
{
    alexandria::BosqueIterator bosque = pd->getBosqueBegin();
    checker_.checkString(bosque->bosque, "bosque");
    checker_.checkDouble(bosque->polarizability, "polarizability");
}

/*
   TEST_F (PoldataTest, Dihedral)
   {
   alexandria::DihedralIterator dihedral = pd->getDihedralBegin(0);
   checker_.checkString(dihedral->atom1,"atom1");
   checker_.checkString(dihedral->atom2,"atom2");
   checker_.checkString(dihedral->atom3,"atom3");
   checker_.checkString(dihedral->atom4,"atom4");
   checker_.checkString(dihedral->params,"params");
   checker_.checkDouble(dihedral->dihedral,"dihedral");
   checker_.checkDouble(dihedral->sigma,"sigma");
   checker_.checkInteger(dihedral->ntrain,"ntrain");
   }*/


TEST_F (PoldataTest, chi)
{
    std::vector<double>      chi0s;


    for (int model = 0; model < numModels; model++)
    {
        chi0s.push_back(pd->getChi0((ChargeDistributionModel)model, atomName));
    }
    checker_.checkSequence(chi0s.begin(), chi0s.end(), "chi");
}

TEST_F (PoldataTest, row){
    std::vector<double>      rows;
    int numAtoms = 3;

    for (int atomNr = 0; atomNr < numAtoms; atomNr++)
    {
        for (int model = 0; model <  numModels; model++)
        {
            rows.push_back(pd->getRow((ChargeDistributionModel)model, atomName, atomNr));
        }
    }
    checker_.checkSequence(rows.begin(), rows.end(), "row");
}


TEST_F (PoldataTest, zeta)
{
    std::vector<double>      zetas;
    int numAtoms = 3;

    for (int atomNr = 0; atomNr < numAtoms; atomNr++)
    {
        for (int model = 0; model <  numModels; model++)
        {
            zetas.push_back(pd->getZeta((ChargeDistributionModel)model, atomName, atomNr));
        }
    }
    checker_.checkSequence(zetas.begin(), zetas.end(), "zeta");
}

TEST_F (PoldataTest, forceField)
{
    std::string force =  pd->getForceField( );
    checker_.checkString(force, "forceField");
}

TEST_F (PoldataTest, lenghtUnit)
{
    std::string lenght =  pd->getLengthUnit( );
    checker_.checkString(lenght, "lenghtUnit");
}


TEST_F (PoldataTest, polarUnit)
{
    std::string polarUnit = pd->getPolarUnit( );
    checker_.checkString(polarUnit, "polarUnit");
}


TEST_F (PoldataTest, polarRef)
{
    std::string polarRef =  pd->getPolarRef( );
    checker_.checkString(polarRef, "polarRef");
}
