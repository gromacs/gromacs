/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
   \brief
   Tests that check whether \Gromacs can use the underlying FMM implementation

   \author Carsten Kutzner <ckutzne@gwdg.de>
   \author R. Thomas Ullmann <tullman@gwdg.de>

   \ingroup module_fmm
 */

#include "gmxpre.h"

#include "gromacs/fmm/fmm.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/stringutil.h"

#include "testutils/testasserts.h"

namespace
{

/* \brief Test whether we get the correct Coulomb energy and forces for a simple geometric arrangement
 *
 */
class FmmInterfaceTest : public ::testing::Test
{
    public:
        FmmInterfaceTest() {}

        /*! \brief
         * Tests whether FMM yields correct results for a simple two-charge system.
         *
         * We put some values into the energy and forces that are passed to the FMM to also
         * make sure that the FMM does not simply overwrite these values but adds to them.
         *
         * \param[in]    ePBC            Type of periodic boundary condition
         * \param[in]    expectedEnergy  Expected value of the Coulomb energy for this charge configuration
         * \param[in]    expectedForce   Expected value of the Coulomb force in x-direction (others are zero)
         */
        void twoChargesANanometerApart(int ePBC, gmx_unused double expectedEnergy, gmx_unused double expectedForce);

};


void FmmInterfaceTest::twoChargesANanometerApart(gmx_unused int ePBC, gmx_unused double expectedEnergy, gmx_unused double expectedForce)
{
    // Instantiate a fast multipole force provider module
    auto fmm(gmx::createFastMultipoleModule());

    // Prepare MDP input
    gmx::KeyValueTreeBuilder     mdpValues;
    mdpValues.rootObject().addValue("fmm-precision", gmx::formatString("0.001"));
    gmx::KeyValueTreeTransformer transform;
    transform.rules()->addRule().keyMatchType("/", gmx::StringCompareType::CaseAndDashInsensitive);
    auto                         result = transform.transform(mdpValues.build(), nullptr);

    gmx::Options                 moduleOptions;
    fmm->mdpOptionProvider()->initMdpOptions(&moduleOptions);
    gmx::assignOptionsFromKeyValueTree(&moduleOptions, result.object(), nullptr);

#ifndef GMX_WITH_FMM
    // Common message if the tests are run without an underlying FMM implementation
    fprintf(stderr,  "GROMACS compiled without FMM - cannot test FMM energies & forces!\n");
#else
    // Prepare variables needed for the test
    t_mdatoms          *md          = new t_mdatoms;
    md->homenr                      = 2;
    md->chargeA                     = new real[md->homenr];
    md->chargeA[0]                  =  1;
    md->chargeA[1]                  = -1;
    matrix             box          = { { 3.0, 0.0, 0.0 }, { 0.0, 3.0, 0.0 }, { 0.0, 0.0, 3.0 } };
    PaddedRVecVector   positions    = { { 1.0, 1.5, 1.5 }, { 2.0, 1.5, 1.5 } }; // Na, Cl from the NaCl.gro test systems in programs/mdrun/tests
    PaddedRVecVector   preFmmforces = { { 1, 2, 3 }, { 4, 5, 6 } };             // offset to check that GROMACS energy and forces are not overwritten by FMM
    PaddedRVecVector   forces       = preFmmforces;
    t_commrec         *cr           = init_commrec();
    t_inputrec        *ir           = new t_inputrec;
    ir->coulombtype                 = eelFMM;
    ir->ePBC                        = ePBC;
    gmx_mtop_t         mtop;
    init_mtop(&mtop);
    mtop.natoms                     = 2;
    double             preFmmEnergy = 123.456;       // offset to check that GROMACS energy and forces are not overwritten by FMM
    double             energy       = preFmmEnergy;  // apply offset
    gmx_enerdata_t    *enerd        = new gmx_enerdata_t;
    init_enerdata(1, 0, enerd);

    // Test the FMM force provider module
    ForceProviders                forceProvider;
    gmx::ForceProviderInitOptions initOptions(ir->ePBC, ir->coulombtype, &mtop);
    gmx::ForceWithVirial          forceWithVirial(forces, TRUE);

    fmm->initForceProviders(&forceProvider, &initOptions);
    forceProvider.calculateForces(cr, md, box, 0.0, as_rvec_array(positions.data()), &forceWithVirial, enerd);
    energy += enerd->grpp.ener[egCOULSR][0];

    // Subtract pre-Fmm values:
    energy -= preFmmEnergy;
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            forceWithVirial.force_[i][j] -= preFmmforces[i][j];
        }
    }

    // Check the Coulomb energy and forces:
    EXPECT_NEAR(expectedEnergy, energy, 5e-5);

    double forceTolerance = 0.0005;
    EXPECT_NEAR(-expectedForce, forceWithVirial.force_[0][XX], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[0][YY], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[0][ZZ], forceTolerance);

    EXPECT_NEAR(+expectedForce, forceWithVirial.force_[1][XX], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[1][YY], forceTolerance);
    EXPECT_NEAR(           0.0, forceWithVirial.force_[1][ZZ], forceTolerance);

    // Clean up
    done_commrec(cr);
    destroy_enerdata(enerd);
    delete enerd;
    delete ir;
    delete [] md->chargeA;
    delete md;
#endif
}


TEST_F(FmmInterfaceTest, OpenBoundaries)
{
    twoChargesANanometerApart(epbcNONE, -ONE_4PI_EPS0, -ONE_4PI_EPS0);
}

TEST_F(FmmInterfaceTest, PeriodicBoundaries)
{
    twoChargesANanometerApart(epbcXYZ, -151.5503252, -109.858909);
}

} // namespace
