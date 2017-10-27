/*
 * Graph.cpp
 *
 *  Created on: Jan 27, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include "Graph.h"
#include "Index.h"
#include "gromacs/utility/fatalerror.h"

namespace fda_analysis {

Graph::Graph(std::vector<double> const& forceMatrix, rvec *coord, int *index, int isize)
{
	double force;
	bool addNode1, addNode2;
	int dim = sqrt(forceMatrix.size());
	std::vector<int>::iterator iterFind;

    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
        	force = forceMatrix[i*dim+j];
            if (force == 0.0) continue;
			Node node1(i), node2(j);

			iterFind = std::find(indices_.begin(), indices_.end(), i);
			if (iterFind != indices_.end()) {
				nodes_[std::distance(indices_.begin(), iterFind)].addConnectedNode(j, force);
				addNode1 = false;
			} else {
				addNode1 = true;
			}

			iterFind = std::find(indices_.begin(), indices_.end(), j);
			if (iterFind != indices_.end()) {
				nodes_[std::distance(indices_.begin(), iterFind)].addConnectedNode(i, force);
				addNode2 = false;
			} else {
				addNode2 = true;
			}

			if (addNode1) {
				node1.addConnectedNode(j, force);
                nodes_.push_back(node1);
                indices_.push_back(i);
			}

			if (addNode2) {
				node2.addConnectedNode(i, force);
                nodes_.push_back(node2);
                indices_.push_back(j);
			}
        }
    }
    for (int i = 0; i < dim; ++i)
    {
    	if (std::find(indices_.begin(), indices_.end(), i) == indices_.end()){
			nodes_.push_back(Node(i));
			indices_.push_back(i);
		}
    }

    updateCoordinates(coord, index, isize);
}

void Graph::updateCoordinates(rvec *coord, int *index, int isize)
{
    std::vector<int>::iterator iterFind;

    // residue-based
    if (index) {
        for (int i = 0; i != isize; ++i)
        {
            iterFind = std::find(indices_.begin(), indices_.end(), i);
            if (iterFind == indices_.end()) gmx_fatal(FARGS, "Error in insertPDBInfo.");
            int pos = std::distance(indices_.begin(), iterFind);

            /// Values will be converted from nm into Angstrom.
            nodes_[pos].x_ = coord[index[i]][0] * 10.0;
            nodes_[pos].y_ = coord[index[i]][1] * 10.0;
            nodes_[pos].z_ = coord[index[i]][2] * 10.0;
        }
    // atomic-based
    } else {
        for (size_t i = 0; i != indices_.size(); ++i)
        {
            /// Values will be converted from nm into Angstrom.
            nodes_[i].x_ = coord[indices_[i]][0] * 10.0;
            nodes_[i].y_ = coord[indices_[i]][1] * 10.0;
            nodes_[i].z_ = coord[indices_[i]][2] * 10.0;
        }
    }
}

void Graph::convertInPDBMinGraphOrder(std::string const& outFilename, double threshold,
    size_t minGraphOrder, bool onlyBiggestNetwork, bool append) const
{
    Networks networks;
    createNetworkMinGraphOrder(networks, threshold, minGraphOrder, onlyBiggestNetwork);
    convertNetworkToPDB(outFilename, networks, threshold, minGraphOrder, append);
}

void Graph::convertInDIMACSMinGraphOrder(std::string const& outFilename, double threshold,
    size_t minGraphOrder, bool onlyBiggestNetwork) const
{
    Networks networks;
    createNetworkMinGraphOrder(networks, threshold, minGraphOrder, onlyBiggestNetwork);
    convertNetworkToDIMACS(outFilename, networks, threshold, minGraphOrder);
}

void Graph::createNetworkMinGraphOrder(Networks& networks, double threshold, size_t minGraphOrder,
    bool onlyBiggestNetwork) const
{
	int iNetwork = 0;
	int idNetwork = 0;

	for (auto node : nodes_)
	{
		Network mergedNetwork;
		bool found = false;

		for (auto network : index(networks))
		{
			if (std::find(network.value.begin(), network.value.end(), node.index_) != network.value.end()) {
				mergedNetwork.push_back(network.index);
				idNetwork = network.index;
				found = true;
			}
		}
		if (!found) {
			mergedNetwork.push_back(networks.size());
			idNetwork = networks.size();
			networks.push_back(Network(1,node.index_));
		}

		bool sup = false;
		for (auto connectedIndex : index(node.connectedIndicies_)) {
			if (std::abs(node.forces_[connectedIndex.index]) >= threshold) {
				sup = true;
				iNetwork = 0;
				bool addNew = true;
				for (auto network : networks) {
					if (std::find(network.begin(), network.end(), connectedIndex.value) != network.end()) {
						mergedNetwork.push_back(iNetwork);
						addNew = false;
						break;
					}
					++iNetwork;
				}
				if (addNew) networks[idNetwork].push_back(connectedIndex.value);
			}
		}
		if (sup) {
			std::sort(mergedNetwork.begin(), mergedNetwork.end());
			mergedNetwork.erase(std::unique(mergedNetwork.begin(), mergedNetwork.end()), mergedNetwork.end());
			if (mergedNetwork.size() > 1) {
				int mergeTimes = 0;
				Network merge;
				for (auto mergeIndex : mergedNetwork) {
					Networks::iterator iter = networks.begin() + (mergeIndex - mergeTimes);
					merge.insert(merge.end(), iter->begin(), iter->end());
					networks.erase(iter);
					++mergeTimes;
				}
				networks.push_back(merge);
			}
		}
	}

	if (onlyBiggestNetwork) {
		size_t count = 0;
		for (auto network : networks) {
			if (network.size() > count) count = network.size();
			if (count > minGraphOrder) minGraphOrder = count;
		}
	}
}

void Graph::convertNetworkToPDB(std::string const& filename, Networks const& networks, double threshold,
    size_t minGraphOrder, bool append) const
{
    std::ofstream pdb;
    if (append) pdb.open(filename, std::ofstream::app);
    else pdb.open(filename);
    if (!pdb) gmx_fatal(FARGS, "Error opening pdb file.");

    std::stringstream connections;
    std::vector< std::pair<int, int> > finishedPairs;
    std::vector<int>::const_iterator iterFind;
    Node node, connectedNode;
    int numNetwork = 0;
    int n = 1;

    real currentForce;
    bool valueToLargeForPDB = false;

    for (auto network : networks) {
        if (network.size() >= minGraphOrder) {
        	for (auto nodeId : network) {
        		iterFind = std::find(indices_.begin(), indices_.end(), nodeId);
        		if (iterFind == indices_.end()) gmx_fatal(FARGS, "Error in insertPDBInfo.");
                node = nodes_[std::distance(indices_.begin(), iterFind)];
                for (auto connectedIndex : index(node.connectedIndicies_)) {
                    if (std::abs(node.forces_[connectedIndex.index]) >= threshold) {
                        int pos = std::distance(indices_.begin(), std::find(indices_.begin(), indices_.end(), connectedIndex.value));
                    	connectedNode = nodes_[pos];
                    	bool add = true;
                    	for (auto finishedPair : finishedPairs) {
                    		if ((node.index_ == finishedPair.first and
                    			 connectedNode.index_ == finishedPair.second) or
                    			(node.index_ == finishedPair.second and
                    		     connectedNode.index_ == finishedPair.first))
                    			add = false;
                    	}
                    	if (add) {
                    	    currentForce = node.forces_[connectedIndex.index];
                    	    if (currentForce > 999.99) valueToLargeForPDB = true;
                    		writeAtomToPDB(pdb, n, node, currentForce, numNetwork);
                    		writeAtomToPDB(pdb, n+1, connectedNode, currentForce, numNetwork);
                            connections << "CONECT" << std::setw(5) << n << std::setw(5) << n+1 << std::endl;
                            finishedPairs.push_back(std::make_pair(node.index_, connectedNode.index_));
                            n += 2;
                    	}
                    }
                }
        	}
        	++numNetwork;
        }
    }

    if (numNetwork > 10) gmx_warning("%d networks are found, which could be difficult to visualize", numNetwork);
    if (virginValueToLargeForPDB and valueToLargeForPDB) {
        gmx_warning("Force values larger than 999.99 are detected. Therefore, the general PDB format of the b-factor column of Real(6.2) is broken. "
                    "It is tested that it works for Pymol and VMD, but it is not guaranteed that it will work for other visualization programs.");
        virginValueToLargeForPDB = false;
    }

    pdb << connections.str() << "ENDMDL" << std::endl;
}

void Graph::convertNetworkToDIMACS(std::string const& outFilename, Networks const& networks, double threshold,
    size_t minGraphOrder) const
{
    std::ofstream dimacs(outFilename);
    std::stringstream edgeSection;
    std::set<int> addedNodes;
    std::map<int, int> corresp;
    std::vector< std::pair<int, int> > finishedPairs;
    std::vector<int>::const_iterator iterFind;
    Node node, connectedNode;
    int numNetwork = 0;
    int n = 1;

    for (auto network : networks) {
        if (network.size() >= minGraphOrder) {
        	for (auto nodeId : network) {
        		iterFind = std::find(indices_.begin(), indices_.end(), nodeId);
        		if (iterFind == indices_.end()) gmx_fatal(FARGS, "Error in insertPDBInfo.");
                node = nodes_[std::distance(indices_.begin(), iterFind)];
                for (auto connectedIndex : index(node.connectedIndicies_)) {
                    if (std::abs(node.forces_[connectedIndex.index]) >= threshold) {
                        int pos = std::distance(indices_.begin(), std::find(indices_.begin(), indices_.end(), connectedIndex.value));
                    	connectedNode = nodes_[pos];
                    	bool add = true;
                    	for (auto finishedPair : finishedPairs) {
                    		if ((node.index_ == finishedPair.first and
                    			 connectedNode.index_ == finishedPair.second) or
                    			(node.index_ == finishedPair.second and
                    		     connectedNode.index_ == finishedPair.first))
                    			add = false;
                    	}
                    	if (add) {
                            if (addedNodes.find(node.index_) == addedNodes.end()) {
                            	addedNodes.insert(node.index_);
                            	corresp[node.index_] = n;
                            	++n;
                            }
                            if (addedNodes.find(connectedNode.index_) == addedNodes.end()) {
                            	addedNodes.insert(connectedNode.index_);
                            	corresp[connectedNode.index_] = n;
                            	++n;
                            }
                            edgeSection << "e " << corresp[node.index_] << " " << corresp[connectedNode.index_] << std::endl;
                            finishedPairs.push_back(std::make_pair(node.index_, connectedNode.index_));
                    	}
                    }
                }
        	}
        	++numNetwork;
        }
    }

    dimacs << "p edge " << addedNodes.size() << " " << finishedPairs.size() << std::endl;
    for (int i = 1; i != n; ++i) dimacs << "n " << i << " 1" << std::endl;
    dimacs << edgeSection.str();
}

void Graph::writeAtomToPDB(std::ofstream& os, int num, Node const& node, double force,
    int numNetwork) const
{
	os << "ATOM"
	   << std::setw(7) << num
	   << std::setw(15) << node.index_
       << std::fixed << std::setprecision(3) << std::setw(12) << node.x_
       << std::setw(8) << node.y_
       << std::setw(8) << node.z_
       << std::setprecision(2) << std::setw(6) << 1.0
       << std::setprecision(2) << std::setw(6) << force
       << std::setw(8) << colors[numNetwork % 32]
       << std::endl;
}

std::vector<std::string> Graph::colors =
    {"AA", "AB", "AC", "AD", "AE", "AF", "AG", "AH", "AI", "AJ",
     "AK", "AL", "AM", "AN", "AO", "AP", "AQ", "AR", "AS", "AT",
     "AU", "AV", "AW", "AX", "AY", "AZ", "BA", "BB", "BC", "BD",
     "BE", "BF"};

bool Graph::virginValueToLargeForPDB = true;

std::ostream& operator << (std::ostream& os, Graph const& graph)
{
	os << "Indices:" << std::endl;
	for (auto index : graph.indices_) {
		os << index << std::endl;
	}
	os << "Nodes:" << std::endl;
	for (auto node : graph.nodes_) {
		os << node.index_ << " "
	       << node.x_ << " "
	       << node.y_ << " "
	       << node.z_ << std::endl;
	}
	return os;
}

} // namespace fda_analysis
