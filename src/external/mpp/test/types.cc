#include <gtest/gtest.h>

#include <mpp.h>
#include <iostream>
#include <complex>

#include "util/printer.h"

using namespace mpi;

TEST(Type, Char) {
	if(comm::world.rank() == 0) {
		comm::world(1) << 'a';
	} else if (comm::world.rank() == 1) {
		char val;
		comm::world(0) >> val;
		EXPECT_EQ('a', val);
	}
}

/*
TEST(Type, Unsigned) {
	if(comm::world.rank() == 0) {
		comm::world(1) << 10u;
	} else if (comm::world.rank() == 1) {
		unsigned val;
		comm::world(0) >> val;
		EXPECT_EQ(10u, val);
	}
}

TEST(Type, UnsignedLong) {
	if(comm::world.rank() == 0) {
		comm::world(1) << 1000000ul;
	} else if (comm::world.rank() == 1) {
		unsigned long val;
		comm::world(0) >> val;
		EXPECT_EQ(1000000ul, val);
	}
}

TEST(Type, LongDouble) {
	if(comm::world.rank() == 0) {
		long double ld = 0.000000010;
		comm::world(1) << ld;
	} else if (comm::world.rank() == 1) {
		long double val;
		comm::world(0) >> val;
		EXPECT_EQ(0.000000010, val);
	}
}

TEST(Type, Complex) {
	if(comm::world.rank() == 0) {
		comm::world(1) << std::complex<float>(10.3, 2.4);
	} else if (comm::world.rank() == 1) {
		std::complex<float> val;
		comm::world(0) >> val;
		EXPECT_EQ(10.3f, real(val));
		EXPECT_EQ(2.4f, imag(val));
	}
}
*/

int main(int argc, char** argv) {
	init(&argc, &argv);
	// Disables elapsed time by default.
	::testing::GTEST_FLAG(print_time) = true;

	// This allows the user to override the flag on the command line.
	::testing::InitGoogleTest(&argc, argv);
    
    MPIEventForward::Insert();
	size_t errcode = RUN_ALL_TESTS();
	MPIEventForward::Remove();
    finalize();
	return static_cast<int>(errcode);
}
