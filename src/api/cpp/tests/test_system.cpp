//#include "atoms.h"
#include "gmxapi/gmxapi.h"
#include "gmxapi/system.h"
#include "gmxapi/md.h"
#include <gtest/gtest.h>

TEST(ApiSystem, Construction)
{
    {   // Construction
        auto system = gmxapi::System();
    }   // Destruction

    auto system = gmxapi::fromTprFile("topol.tpr");
    ASSERT_TRUE(system != nullptr);
    ASSERT_TRUE(system->md() != nullptr);
    ASSERT_NO_THROW(system->md()->info());
    ASSERT_STREQ("Generic MDEngine object", system->md()->info().c_str());
//    ASSERT_EQ(system->atoms()->x()->size(), 7);
}

TEST(ApiSystem, Accessors)
{
    auto system = gmxapi::fromTprFile("topol.tpr");
//    ASSERT_TRUE(system->atoms() != nullptr);
//    ASSERT_TRUE(system->atoms()->x() != nullptr);
}
