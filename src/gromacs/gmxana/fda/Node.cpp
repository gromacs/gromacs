/*
 * Node.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include "Node.h"
#include <iostream>

namespace fda_analysis {

Node::Node(int index)
 : index_(index), x_(0.0), y_(0.0), z_(0.0), ps_(0.0)
{}

void Node::addConnectedNode(int index, double force)
{
	connectedIndicies_.push_back(index);
	forces_.push_back(force);
}

} // namespace fda_analysis
