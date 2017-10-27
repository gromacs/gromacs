/*
 * Node.h
 *
 *  Created on: Jan 27, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef NODE_H_
#define NODE_H_

#include <string>
#include <vector>

namespace fda_analysis {

// forwarding for printing node informations
class Graph;

/**
 * Special node implementation compatible to Maximes Python script.
 */
class Node
{
public:

	Node(int index = 0);

	void addConnectedNode(int index, double force);

	int getIndex() const { return index_; }

	std::vector<int> const& getConnectedIndicies() const { return connectedIndicies_; }

private:

	friend class Graph;

	friend std::ostream& operator << (std::ostream& os, Graph const& graph);

	int index_;

	std::vector<int> connectedIndicies_;

	std::vector<double> forces_;

    /// Coordinate in Angstrom.
	double x_;

    /// Coordinate in Angstrom.
	double y_;

    /// Coordinate in Angstrom.
	double z_;

	double ps_;

};

} // namespace fda_analysis

#endif /* NODE_H_ */
