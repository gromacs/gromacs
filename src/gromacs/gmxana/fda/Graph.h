/*
 * Graph.h
 *
 *  Created on: Jan 27, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef GRAPH_H_
#define GRAPH_H_

#include <stddef.h>
#include <iostream>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "Node.h"

namespace fda_analysis {

/**
 * Special graph implementation compatible to Maximes Python script.
 */
class Graph
{
public:

	Graph() {};

	//! Build graph by adjacency matrix
	Graph(std::vector<double> const& forceMatrix, rvec *coord, int *index, int isize);

	void convertInPDBMinGraphOrder(std::string const& outFilename, double threshold,
		size_t minGraphOrder, bool onlyBiggestNetwork, bool append) const;

	void convertInDIMACSMinGraphOrder(std::string const& outFilename, double threshold,
		size_t minGraphOrder, bool onlyBiggestNetwork) const;

	void updateCoordinates(rvec *coord, int *index, int isize);

private:

	typedef std::vector<int> Network;
	typedef std::vector<Network> Networks;

	void createNetworkMinGraphOrder(Networks& networks, double threshold, size_t minGraphOrder, bool onlyBiggestNetwork) const;

	void convertNetworkToPDB(std::string const& outFilename, Networks const& networks, double threshold,
		size_t minGraphOrder, bool append) const;

	void convertNetworkToDIMACS(std::string const& outFilename, Networks const& networks, double threshold,
		size_t minGraphOrder) const;

	void writeAtomToPDB(std::ofstream& os, int num, Node const& node, double force,
        int numNetwork) const;

	friend std::ostream& operator << (std::ostream& os, Graph const& graph);

	std::vector<Node> nodes_;

	std::vector<int> indices_;

    static std::vector<std::string> colors;

    /// Flag for printing warning only one time per run.
    static bool virginValueToLargeForPDB;

};

std::ostream& operator << (std::ostream& os, Graph const& graph);

} // namespace fda_analysis

#endif /* GRAPH_H_ */
