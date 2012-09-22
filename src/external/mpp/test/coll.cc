#include <gtest/gtest.h>

#define MPP_CXX11_ARRAY
#define MPP_CXX11_RVALREF
#include <mpp.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <array>
#include <list>

#include "util/printer.h"

using namespace mpi;

typedef std::array<int,5> int_array;

//Tests to make sure the collective optimization are actually used 
//(Requires is_static/is_continunous to work correctly)
TEST(IsStatic, Misc) {
    EXPECT_TRUE(mpi_type_traits<const int>::is_static::value);
    EXPECT_TRUE(mpi_type_traits<int_array>::is_static::value);
    EXPECT_FALSE(mpi_type_traits<boost::scoped_ptr<int>>::is_static::value);
    EXPECT_FALSE(mpi_type_traits<std::vector<int>>::is_static::value);
}

TEST(IsContigous, Vector) {
    EXPECT_TRUE(is_contiguous<std::vector<int>::iterator>::value);
    EXPECT_TRUE(is_contiguous<std::vector<int>::const_iterator>::value);
    EXPECT_TRUE(is_contiguous<std::vector<int_array>::iterator>::value);
    EXPECT_FALSE(is_contiguous<std::vector<boost::scoped_ptr<int> >::const_iterator>::value);
}

TEST(IsContigous, MiscContainers) {
    EXPECT_TRUE(is_contiguous<int_array::iterator>::value);
    EXPECT_FALSE(is_contiguous<std::list<int>::iterator>::value);
    EXPECT_TRUE(is_contiguous<int*>::value);
    EXPECT_TRUE(is_contiguous<const int*>::value);
    EXPECT_TRUE(is_contiguous<const int* const>::value);
}

TEST(Gather, ScalarToVector) {
    std::vector<int> recv(comm::world.size());
    int send = comm::world.rank()+1;
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i], i+1);
        }
	}
}

TEST(Gather, VectorToVector) {
    std::vector<int> recv(comm::world.size()*2);
    std::vector<int> send = { comm::world.rank()+1, comm::world.rank()*2};
    comm::world.gather(send.begin(), send.end(), recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i*2], i+1);
            EXPECT_EQ(recv[i*2+1], i*2);
        }
	}
}

TEST(Gather, ArrayToVector) {
    std::vector<int> recv(comm::world.size()*2);
    std::array<int,2> send = {{ comm::world.rank()+1, comm::world.rank()*2}};
    comm::world.gather(send.begin(), send.end(), recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i*2], i+1);
            EXPECT_EQ(recv[i*2+1], i*2);
        }
    }
}

TEST(Gather, VectorToVectorOfVector) {
    std::vector<std::vector<int>> recv(comm::world.size(),std::vector<int>(2));
    std::vector<int> send = { comm::world.rank()+1, comm::world.rank()*2};
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i][0], i+1);
            EXPECT_EQ(recv[i][1], i*2);
        }
    }
}

TEST(Gather, ArrayToVectorOfArray) {
    std::vector<std::array<int,2>> recv(comm::world.size());
    std::array<int,2> send = {{ comm::world.rank()+1, comm::world.rank()*2}};
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i][0], i+1);
            EXPECT_EQ(recv[i][1], i*2);
        }
    }
}

TEST(Gather, PointerToVector) {
    std::vector<int> recv(comm::world.size());
    boost::scoped_ptr<int> send(new int(comm::world.rank()+1));
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i], i+1);
        }
	}
}

TEST(Gather, PointerToPtrVector) {
    std::vector<std::unique_ptr<int> > recv(comm::world.size());
    for (auto it = recv.begin(); it<recv.end(); ++it) {
        it->reset(new int);
    }
    boost::scoped_ptr<int> send(new int(comm::world.rank()+1));
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(*recv[i], i+1);
        }
    }
}

TEST(Gather, ScalarToList) {
    std::list<int> recv(comm::world.size());
    int send = comm::world.rank()+1;
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        std::list<int>::iterator it = recv.begin();
        for (int i=0; i<comm::world.size(); ++i, ++it) {
            EXPECT_EQ(*it, i+1);
        }
	}
}

TEST(Gather, ListToList) {
    std::list<int> recv(comm::world.size()*2);
    std::list<int> send = { comm::world.rank()+1, comm::world.rank()*2};
    comm::world.gather(send.begin(), send.end(), recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        std::list<int>::iterator it = recv.begin();
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(*(it++), i+1);
            EXPECT_EQ(*(it++), i*2);
        }
    }
}

struct S1 {
    char a;
    int b;
    char c;
};
namespace mpi {
template <>
inline MPI_Datatype mpi_type_traits<S1>::get_type(const S1& s) {
    mpi_type_builder builder(s, 2);
    builder.add(s.b);
    builder.add(s.c);
    return builder.build();
}
SET_MPI_STATIC(S1);
}

TEST(Gather, StaticStructVectorToVector) {
    std::vector<S1> recv(comm::world.size()*2);
    std::vector<S1> send(2);
    send[0].b = comm::world.rank()+1;
    send[0].c = 42;
    send[1].b = comm::world.rank()*2;
    send[1].c = 7;
    comm::world.gather(send.begin(), send.end(), recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i*2].b, i+1);
            EXPECT_EQ(recv[i*2].c, 42);
            EXPECT_EQ(recv[i*2+1].b, i*2);
            EXPECT_EQ(recv[i*2+1].c, 7);
        }
    }
}

struct S2 {
    S2() : b(2) {};
    enum E { v1, v2};
    int a;
    std::vector<int> b;
    E c;
};
namespace mpi {
template <>
inline MPI_Datatype mpi_type_traits<S2>::get_type(const S2& s) {
    mpi_type_builder builder(s, 3);
    builder.add(s.a);
    builder.add(s.b);
    builder.add(s.c);
    return builder.build();
}
template <>
inline MPI_Datatype mpi_type_traits<S2::E>::get_type(const S2::E& e) {
    return get_type_enum<S2::E>(e);
}
SET_MPI_STATIC(S2::E);
}

TEST(Gather, StructToVector) {
    std::vector<S2> recv(comm::world.size());
    S2 send;
    send.a = comm::world.rank()+3;
    send.b = {comm::world.rank()+1, comm::world.rank()*2};
    send.c = S2::v2;
    comm::world.gather(&send, (&send)+1, recv.begin(), recv.end(), 0);
    if(comm::world.rank() == 0) {
        for (int i=0; i<comm::world.size(); i++) {
            EXPECT_EQ(recv[i].a, i+3);
            EXPECT_EQ(recv[i].b[0], i+1);
            EXPECT_EQ(recv[i].b[1], i*2);
            EXPECT_EQ(recv[i].c, S2::E::v2);
        }
    }
}

TEST(Reduce, Scalar) {
    int send = comm::world.rank()+1;
    int recv;
    comm::world.reduce(&send, (&send)+1, &recv, (&recv)+1, MPI_SUM, 0);
    if(comm::world.rank() == 0) {
        int n = comm::world.size();
        EXPECT_EQ(recv, n*(n+1)/2);
    }
}

TEST(Reduce, ScalarInPlace) {
    int send = comm::world.rank()+1;
    comm::world.reduce(&send, (&send)+1, MPI_SUM, 0);
    if(comm::world.rank() == 0) {
        int n = comm::world.size();
        EXPECT_EQ(send, n*(n+1)/2);
    }
}

TEST(Reduce, VectorInPlace) {
    std::vector<int> send = {comm::world.rank()+1, 1};
    comm::world.reduce(send.begin(), send.end(), MPI_SUM, 0);
    if(comm::world.rank() == 0) {
        int n = comm::world.size();
        EXPECT_EQ(send[0], n*(n+1)/2);
        EXPECT_EQ(send[1], n);
    }
}

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
