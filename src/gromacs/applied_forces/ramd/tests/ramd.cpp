/*! \internal \file
 * \brief
 * Tests for functionality of the RAMD module.
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "gromacs/applied_forces/ramd/ramd.h"

#include <filesystem>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdrunutility/mdmodulesnotifiers.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider_test_helper.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/options.h"
#include "gromacs/options/treesupport.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringcompare.h"
#include "gromacs/utility/vec.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

class RAMDTest : public ::testing::Test
{
public:
    void addMdpOptionRAMDActive()
    {
        mdpValueBuilder_.rootObject().addValue("ramd-active", std::string("yes"));
    }

    void makeRAMDModuleWithSetOptions()
    {
        KeyValueTreeObject mdpOptionsTree = mdpValueBuilder_.build();

        ramdModule_ = RAMDModuleInfo::create();

        test::fillOptionsFromMdpValues(mdpOptionsTree, ramdModule_->mdpOptionProvider());
    }

    void initializeForceProviders()
    {
        ramdModule_->initForceProviders(&ramdForces_);
    }

protected:
    KeyValueTreeBuilder        mdpValueBuilder_;
    ForceProviders             ramdForces_;
    std::unique_ptr<IMDModule> ramdModule_;
};

TEST_F(RAMDTest, ForceProviderLackingInputThrows)
{
    // Prepare MDP inputs
    addMdpOptionRAMDActive();

    makeRAMDModuleWithSetOptions();

    // Build the force provider, once all input data is gathered
    EXPECT_ANY_THROW(initializeForceProviders());
}

} // namespace

} // namespace gmx
