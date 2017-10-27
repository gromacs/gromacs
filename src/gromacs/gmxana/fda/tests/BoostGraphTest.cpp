/*
 * BoostGraphTest.cpp
 *
 *  Created on: Feb 19, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <gtest/gtest.h>
#include "gromacs/gmxana/fda/BoostGraph.h"
#include "testutils/EqualArrays.h"

using namespace fda_analysis;

//! Test for BoostGraph
TEST(BoostGraphTest, Dijkstra1)
{
	const int N = 4;
	std::vector<double> f(N*N);
	f[      1] = 7.0;
	f[      2] = 1.0;
	f[1*N + 2] = 1.0;
	f[1*N + 3] = 1.0;
	f[2*N + 3] = 5.0;

	BoostGraph graph(f);
    BoostGraph::Path shortestPath = graph.findShortestPath(0, 3);

    std::vector<int> refPath{0, 2, 1, 3};
    std::vector<int> shortestPathIndices;
    for (auto const& v : shortestPath) shortestPathIndices.push_back(v);

	EXPECT_TRUE((EqualArrays(refPath, shortestPathIndices)));

	EXPECT_EQ(graph.distance(shortestPath), 3.0);
}

//! Test for BoostGraph
TEST(BoostGraphTest, Dijkstra2)
{
	const int N = 6;
	std::vector<double> f(N*N);
	f[      1] = 3.0;
	f[      2] = 2.0;
	f[1*N + 2] = 1.0;
	f[1*N + 3] = 4.0;
	f[2*N + 3] = 2.0;
	f[2*N + 4] = 3.0;
	f[3*N + 4] = 2.0;
	f[3*N + 5] = 1.0;
	f[4*N + 5] = 2.0;

	BoostGraph graph(f);
    BoostGraph::Path shortestPath = graph.findShortestPath(0, 5);

    std::vector<int> refPath{0, 2, 3, 5};
    std::vector<int> shortestPathIndices;
    for (auto const& v : shortestPath) shortestPathIndices.push_back(v);

	EXPECT_TRUE((EqualArrays(refPath, shortestPathIndices)));

	EXPECT_EQ(graph.distance(shortestPath), 5.0);
}

//! Test for BoostGraph
TEST(BoostGraphTest, Yen1)
{
	const int N = 6;
	std::vector<double> f(N*N);
	f[      1] = 3.0;
	f[      2] = 2.0;
	f[1*N + 2] = 1.0;
	f[1*N + 3] = 4.0;
	f[2*N + 3] = 2.0;
	f[2*N + 4] = 3.0;
	f[3*N + 4] = 2.0;
	f[3*N + 5] = 1.0;
	f[4*N + 5] = 2.0;

	BoostGraph graph(f);
    BoostGraph::PathList kShortestPaths = graph.findKShortestPaths(0, 5, 1);

    std::vector< std::vector<int> > refPath{{0, 2, 3, 5}};

    for (size_t i(0); i != kShortestPaths.size(); ++i)
    {
        std::vector<int> pathIndices;
    	for (auto const& v : kShortestPaths[i]) pathIndices.push_back(v);
    	EXPECT_TRUE((EqualArrays(refPath[i], pathIndices)));
    }
}

//! Test for BoostGraph
TEST(BoostGraphTest, Yen2)
{
	const int N = 6;
	std::vector<double> f(N*N);
	f[      1] = 3.0;
	f[      2] = 2.0;
	f[1*N + 2] = 1.0;
	f[1*N + 3] = 4.0;
	f[2*N + 3] = 2.0;
	f[2*N + 4] = 3.0;
	f[3*N + 4] = 2.0;
	f[3*N + 5] = 1.0;
	f[4*N + 5] = 2.0;

	BoostGraph graph(f);
    BoostGraph::PathList kShortestPaths = graph.findKShortestPaths(0, 5, 3);

    std::vector< std::vector<int> > refPath{{0, 2, 3, 5},
    	                                    {0, 1, 2, 3, 5},
                                            {0, 2, 4, 5}};

    for (size_t i(0); i != kShortestPaths.size(); ++i)
    {
        std::vector<int> pathIndices;
    	for (auto const& v : kShortestPaths[i]) pathIndices.push_back(v);
    	EXPECT_TRUE((EqualArrays(refPath[i], pathIndices)));
    }
}
