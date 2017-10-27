/*
 * BoostGraph.cpp
 *
 *  Created on: Feb 13, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <sstream>
#include "BoostGraph.h"
#include "Index.h"
#include "gromacs/utility/fatalerror.h"

namespace fda_analysis {

BoostGraph::BoostGraph(std::vector<double> const& forceMatrix)
{
	int dim = sqrt(forceMatrix.size());
	double force;

    for (int i = 0; i < dim; ++i) {
        for (int j = i; j < dim; ++j) {
        	force = forceMatrix[i*dim+j];
        	if (force != 0.0) {
        		add_edge(i, j, force, graph_);
        	}
        }
    }

	#ifdef PRINT_DEBUG
		std::cout << "Num vertices = " << num_vertices(graph_) << std::endl;
		std::cout << "Num edges = " << num_edges(graph_) << std::endl;
	#endif
}

BoostGraph::Path BoostGraph::findShortestPath(size_t from, size_t to) const
{
	return dijkstra(from, to, graph_);
}

BoostGraph::PathList BoostGraph::findKShortestPaths(size_t from, size_t to, size_t num) const
{
	PathList shortestPaths;
	shortestPaths.push_back(findShortestPath(from, to));

	PathList variants;

	for (size_t k = 1; k < num; ++k)
	{
        Path const& previousPath = shortestPaths[k-1];

	    // The spur node ranges from the first node to the next to last node in the previous k-shortest path.
	    for (size_t i = 0; i < previousPath.size() - 1; ++i)
	    {
	  	    // Spur node is retrieved from the previous k-shortest path, k − 1.
	    	Vertex spurNode = previousPath[i];

	  	    // The sequence of nodes from the source to the spur node of the previous k-shortest path.
	  	    Path rootPath(previousPath.begin(), previousPath.begin() + i);

	  	    Graph tmpGraph(graph_);

	  	    for (auto & p : shortestPaths) {
	  	 	    if (rootPath == Path(p.begin(), p.begin() + i)) {
	  	 	        // Remove the links that are part of the previous shortest paths which share the same root path.
	  	            //std::cout << "remove edge " << p[i] << "  " << p[i+1] << std::endl;
	  	 		    remove_edge(vertex(p[i], tmpGraph), vertex(p[i+1], tmpGraph), tmpGraph);
	  	 	    }
	  	    }

            #ifdef PRINT_DEBUG
                std::cout << "tmp graph after remove_edge" << std::endl;
                print(tmpGraph);
            #endif

	  	    // Calculate the spur path from the spur node to the sink.
	  	    Path spurPath;
	  	    try {
	  	    	spurPath = dijkstra(spurNode, to, tmpGraph);
	  	    } catch ( ... ) {
	  	    	continue;
	  	    }

            #ifdef PRINT_DEBUG
                std::cout << "spurPath:" << std::endl;
                print(spurPath);
            #endif

	  	    // Entire path is made up of the root path and spur path.
	  	    Path totalPath = rootPath;
	  	    totalPath.insert(totalPath.end(), spurPath.begin(), spurPath.end());

            #ifdef PRINT_DEBUG
                std::cout << "totalPath:" << std::endl;
                print(totalPath);
            #endif

	  	    // Add the potential k-shortest path to the heap.
	  	    variants.push_back(totalPath);
	    }

		// This handles the case of there being no spur paths, or no spur paths left.
		// This could happen if the spur paths have already been exhausted (added to A),
		// or there are no spur paths at all - such as when both the source and sink vertices
		// lie along a "dead end".
	    if (variants.empty()) break;

	    // Find the potential k-shortest path.
	    PathList::iterator iterShortestPath;
	    double minDistance = std::numeric_limits<double>::max();
  	    for (auto iter(variants.begin()), iterEnd(variants.end()); iter != iterEnd;  ++iter) {
  	    	if (distance(*iter) < minDistance) {
  	    		minDistance = distance(*iter);
  	    		iterShortestPath = iter;
  	    	}
  	    }

	    // Add the lowest cost path becomes the k-shortest path.
        #ifdef PRINT_DEBUG
            std::cout << k << "-shortest path:" << std::endl;
            print(*iterShortestPath);
        #endif
	    shortestPaths.push_back(*iterShortestPath);
	    variants.erase(iterShortestPath);
	}

    return shortestPaths;
}

double BoostGraph::distance(BoostGraph::Path const& path) const
{
	double dist = 0.0;
	for (size_t i = 0; i != path.size() - 1; ++i) dist +=  boost::get(boost::edge_weight, graph_, boost::edge(path[i], path[i+1], graph_).first);
    return dist;
}

BoostGraph::Path BoostGraph::dijkstra(Vertex source, Vertex dest, Graph const& graph) const
{
    //std::cout << "dijkstra: source = " << source << ", dest = " << dest << std::endl;
	//if (vertex(source, graph)) gmx_fatal(FARGS, "Vertex source is not an element of graph.");
	//if (vertex(dest, graph)) gmx_fatal(FARGS, "Vertex dest is not an element of graph.");

	std::vector<Vertex> predecessor(num_vertices(graph));
	std::vector<double> distance(num_vertices(graph));
	Vertex s = vertex(source, graph);

	dijkstra_shortest_paths(graph, s, boost::predecessor_map(&predecessor[0]).distance_map(&distance[0]));

    #ifdef PRINT_DEBUG
		boost::graph_traits<Graph>::vertex_iterator vi, vend;
		for (boost::tie(vi, vend) = boost::vertices(graph); vi != vend; ++vi) {
			std::cout << "distance(" << *vi << ") = " << distance[*vi] << ", ";
			std::cout << "parent(" << *vi << ") = " << predecessor[*vi] << std::endl;
		}
	#endif

	Path shortestPath;
	Vertex cur = vertex(dest, graph);
    while (cur != s) {
    	shortestPath.push_back(cur);
    	if (cur == predecessor[cur]) throw std::runtime_error("No connection between source and dest.");
        cur = predecessor[cur];
    }
	shortestPath.push_back(source);
	std::reverse(shortestPath.begin(), shortestPath.end());
	return shortestPath;
}

void BoostGraph::print(Graph const& graph) const
{
	boost::graph_traits<Graph>::vertex_iterator vi, vend;
	for (boost::tie(vi, vend) = boost::vertices(graph); vi != vend; ++vi) {
		std::cout << *vi << std::endl;
	}
	boost::graph_traits <Graph>::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = edges(graph); ei != ei_end; ++ei) {
		std::cout << source(*ei, graph) << " -> " << target(*ei, graph) << std::endl;
	}
}

void BoostGraph::print(Path const& path) const
{
	for (auto const& node : path) std::cout << node << " ";
	std::cout << std::endl;
}

} // namespace fda_analysis
