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
#include "gromacs/utility/logger.h"

#include "testutils/testasserts.h"
#include "testutils/testfilemanager.h"
#include "testutils/testmatchers.h"

namespace gmx
{

namespace
{

TEST(RAMDTest, ForceProviderLackingInputThrows)
{
    KeyValueTreeBuilder mdpValueBuilder;
    mdpValueBuilder.rootObject().addValue("ramd-active", std::string("yes"));
    KeyValueTreeObject mdpOptionsTree = mdpValueBuilder.build();

    std::unique_ptr<IMDModule> ramdModule = RAMDModuleInfo::create();
    test::fillOptionsFromMdpValues(mdpOptionsTree, ramdModule->mdpOptionProvider());

    MDModulesNotifiers notifiers;
    ramdModule->subscribeToSimulationSetupNotifications(&notifiers);
    MDLogger logger;
    notifiers.simulationSetupNotifier_.notify(logger);

    ForceProviders ramdForces;
    EXPECT_ANY_THROW(ramdModule->initForceProviders(&ramdForces));
}

} // namespace

} // namespace gmx
