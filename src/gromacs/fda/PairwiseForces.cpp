/*
 * PairwiseForces.cpp
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include "PairwiseForces.h"

namespace fda {

template <typename ForceType>
PairwiseForces<ForceType>::PairwiseForces(std::string const& filename)
{
    int i, j;
    ForceType force;
    PairwiseForceList pairwise_forces;
    std::ifstream is(filename);
    std::string token;
    while (is >> token)
    {
        if (token == "frame") {
            size_t frame;
            is >> frame;
            if (frame) {
                // Sort by i, j, and type
                std::sort(pairwise_forces.begin(), pairwise_forces.end(),
                    [](PairwiseForce const& pf1, PairwiseForce const& pf2) {
                        return pf1.i < pf2.i or
                              (pf1.i == pf2.i and pf1.j < pf2.j) or
                              (pf1.i == pf2.i and pf1.j == pf2.j and pf1.force.type < pf2.force.type);
                    });
                all_pairwise_forces.push_back(pairwise_forces);
                pairwise_forces.clear();
            }
            if (frame != all_pairwise_forces.size()) {
                std::cout << "frame = " << frame << std::endl;
                std::cout << "all_pairwise_forces.size() = " << all_pairwise_forces.size() << std::endl;
                throw std::runtime_error("frame numbers are not consecutively");
            }
            continue;
        }

        try {
           i = std::stoi(token);
        } catch ( ... ) {
            std::cerr << token << std::endl;
            throw;
        }

        is >> j >> force;
        pairwise_forces.push_back(PairwiseForce(i,j,force));
    }

    // Don't forget to push the last frame
    // Sort by i, j, and type
    std::sort(pairwise_forces.begin(), pairwise_forces.end(),
        [](PairwiseForce const& pf1, PairwiseForce const& pf2) {
            return pf1.i < pf2.i or
                  (pf1.i == pf2.i and pf1.j < pf2.j) or
                  (pf1.i == pf2.i and pf1.j == pf2.j and pf1.force.type < pf2.force.type);
        });
    all_pairwise_forces.push_back(pairwise_forces);
}

/// template instantiation
template class PairwiseForces<Force<real>>;
template class PairwiseForces<Force<Vector>>;

} // namespace fda
