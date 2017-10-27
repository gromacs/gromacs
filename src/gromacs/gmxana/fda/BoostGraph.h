/*
 * BoostGraph.h
 *
 *  Created on: Feb 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef BOOSTGRAPH_H_
#define BOOSTGRAPH_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <iostream>
#include <vector>

namespace fda_analysis {

/**
 * Special graph implementation fully compatible to Beifei's R script and the igraph library.
 */
class BoostGraph
{
public:

    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        boost::no_property, boost::property<boost::edge_weight_t, double> > Graph;
    typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
    typedef boost::graph_traits<Graph>::edge_descriptor Edge;

    typedef std::vector<int> Path;
    typedef std::vector<Path> PathList;

	BoostGraph() {};

	//! Build graph by adjacency matrix
	BoostGraph(std::vector<double> const& forceMatrix);

	//! Use Dijkstra algorithm to find the shortest path.
	Path findShortestPath(size_t source, size_t dest) const;

	//! Yen's algorithm to find the k shortest paths.
	PathList findKShortestPaths(size_t source, size_t dest, size_t num) const;

	//! Determine the sum of weights (distance) of a graph.
    double distance(Path const& path) const;

private:

	//! Calculate the spur path from source to dest.
	Path dijkstra(Vertex source, Vertex dest, Graph const& graph) const;

	//! Print graph for debugging
	void print(Graph const& graph) const;

	//! Print path for debugging
	void print(Path const& path) const;

    //! Undirected graph connecting atoms or residues weighted by the force between both.
    Graph graph_;

};

} // namespace fda_analysis

#endif /* BOOSTGRAPH_H_ */
