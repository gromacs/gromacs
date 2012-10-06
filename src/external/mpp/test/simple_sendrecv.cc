#include <gtest/gtest.h>

#include <mpp.h>
#include <iostream>
#include <cmath>
#include <chrono>

#include "util/printer.h"

using namespace mpi;

TEST(SendRecv, Scalar) {
	if(comm::world.rank() == 0) {
		comm::world(1) << 4.2;
		int val;
		auto s = comm::world(1) >> val;
		EXPECT_EQ( 4, val);
		EXPECT_EQ( 1, s.source().rank() );
		EXPECT_EQ( 0, s.tag() );
	} else if (comm::world.rank() == 1) {
		double val;
		auto s = comm::world(0) >> val;
		EXPECT_EQ(4.2, val);
		EXPECT_EQ( 0, s.source().rank() );
		EXPECT_EQ(0, s.tag());
		comm::world(0) << static_cast<int>(floor(val));
	}
}

TEST(SendRecv, Array) {
 	int datav[] = {2, 4, 6, 8};
 	std::vector<int> data(datav, datav+sizeof(datav)/sizeof(int));
 	if(comm::world.rank() == 0) {
 		comm::world(1) << data;
 	} else if (comm::world.rank() == 1) {
 		std::vector<int> vec(4);
 		comm::world(0) >> vec;
 		EXPECT_EQ( static_cast<size_t>(4), vec.size() );
 		// check whether received data is equal to original data
 		EXPECT_TRUE( std::equal(vec.begin(), vec.end(), data.begin(), std::equal_to<int>()) );
 	}
 }
/* 
 TEST(SendRecv, Future) {
 	if ( comm::world.rank() == 0 ) {
 		comm::world(1) << 100;
 	} else if(comm::world.rank() == 1) {
 		int k;
 		request<int> r = comm::world(0) > k;
 		r.get();
 		EXPECT_EQ(100, k);
 	}
 }
*/
 TEST(SendRecv, Tags) {
 
 	if ( comm::world.rank() == 0 ) {
 		comm::world(1) << msg<const int>(100, 11);
 		comm::world(1) << msg<const int>(101, 0);
 	} else if(comm::world.rank() == 1) {
 		int k;
 		comm::world(0) >> msg(k,0);
 		EXPECT_EQ(101, k);
 		comm::world(0) >> msg(k,11);
 		EXPECT_EQ(100, k);
 	}
 }
 
 TEST(SendRecv, PingPong) {
 	int p=0;
 	if(comm::world.rank() == 0) {
 		// start the ping
 		comm::world(1) << p;
 	}
 
 	while ( p <= 10 ) {
 		auto ep = (comm::world(mpi::any) >> p ).source();
 		ep << p+1;
 		EXPECT_TRUE(comm::world.rank()==0?p%2!=0:p%2==0);
 	}
 }
// 
// TEST(SendRecv, Lists) {
// 
// 	if ( comm::world.rank() == 0 ) {
// 		int data[] = {1,2,3,4,5};
// 		std::list<int> l(data,data+sizeof(data)/sizeof(int));
// 		comm::world(1) << l;
// 	} else if(comm::world.rank() == 1) {
// 		std::vector<int> l(5);
// 		comm::world(0) >> l;
// 		EXPECT_EQ(static_cast<size_t>(5), l.size());
// 		EXPECT_EQ(1, l[0]);
// 		EXPECT_EQ(2, l[1]);
// 		EXPECT_EQ(3, l[2]);
// 		EXPECT_EQ(4, l[3]);
// 		EXPECT_EQ(5, l[4]);
// 	}
// }


TEST(Performance, MpiScalar) {
	
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	auto bench = [&rank] () {
		// warmup
		if (rank == 0) {
			for(int i=0; i<100000; ++i) 
				MPI_Send(&i,1,MPI_INT,1,0,MPI_COMM_WORLD);
		} else {
			int val;
			for(int i=0; i<100000; ++i) {
				MPI_Status s;
				MPI_Recv(&val,1,MPI_INT,0,0,MPI_COMM_WORLD,&s);
				EXPECT_EQ(i, val);
			}
		}
	};

	MPI_Barrier(MPI_COMM_WORLD);
	bench();
	MPI_Barrier(MPI_COMM_WORLD);

	auto start = std::chrono::system_clock::now();
	bench();
	auto end = std::chrono::system_clock::now();

	size_t elapsed_time = 
		std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

	if (rank==0)
		std::cout << "MpiScalar: " << elapsed_time << " microsec" << std::endl;

}

TEST(Performance, MppScalar) {

	using mpi::comm;
	auto& world = comm::world;

	auto bench = [&]() {
		// warmup
		if (world.rank() == 0) {
			for(int i=0; i<100000; ++i) 
				world(1) << i;
		} else {
			int val;
			for(int i=0; i<100000; ++i) {
				world(0) >> val;
				EXPECT_EQ(i, val);
			}
		}
	};

	MPI_Barrier(MPI_COMM_WORLD);
	// Warmup
	bench();
	MPI_Barrier(MPI_COMM_WORLD);

	auto start = std::chrono::system_clock::now();
	bench();
	auto end = std::chrono::system_clock::now();
	
	size_t elapsed_time = 
		std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
	if (world.rank() == 0)
		std::cout << "MppScalar: " << elapsed_time << " microsec" << std::endl;
}

int main(int argc, char** argv) {
        mpi::init(&argc, &argv);
	// Disables elapsed time by default.
	::testing::GTEST_FLAG(print_time) = true;

	// This allows the user to override the flag on the command line.
	::testing::InitGoogleTest(&argc, argv);

	MPIEventForward::Insert();
	size_t errcode = RUN_ALL_TESTS();
	MPIEventForward::Remove();

	mpi::finalize();
	return static_cast<int>(errcode);
}
