//
// Created by Eric Irrgang on 7/31/17.
//

#include "gmxapi/gmxapi.h"
#include <gtest/gtest.h>

namespace
{

TEST(GmxapiStatusClass, Basic)
{
    {
        auto status = gmxapi::Status();
        ASSERT_FALSE(status.success());
        status = true;
        ASSERT_TRUE(status.success());
    }
    {
        auto status = gmxapi::Status(true);
        ASSERT_TRUE(status.success());
        status = false;
        ASSERT_FALSE(status.success());
    }
    {
        auto status = gmxapi::Status(false);
        ASSERT_FALSE(status.success());
        ASSERT_TRUE(gmxapi::Status(gmxapi::Status(status = true)).success());
    }
}

} // end anonymous namespace
