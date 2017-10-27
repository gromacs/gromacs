/*
 * PairwiseForces.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_PAIRWISEFORCES_H_
#define SRC_GROMACS_FDA_PAIRWISEFORCES_H_

#include <string>
#include <vector>
#include "gromacs/fda/Force.h"
#include "gromacs/fda/Vector.h"
#include "gromacs/utility/real.h"

namespace fda {

/**
 * Read pairwise forces from file into arrays and compare.
 */
template <typename ForceType>
class PairwiseForces
{
public:

    /// Default constructor
    PairwiseForces() {}

    /// Constructor reading files
    PairwiseForces(std::string const& filename);

    template <class Comparer>
    bool equal(PairwiseForces const& other, Comparer const& comparer) const
    {
        if (all_pairwise_forces.size() != other.all_pairwise_forces.size()) return false;

        // std::zip would be nice; avoid boost::zip as boost was ejected by gromacs
        for (size_t frame = 0; frame !=  all_pairwise_forces.size(); ++frame) {
            auto pairwise_forces = all_pairwise_forces[frame];
            auto other_pairwise_forces = other.all_pairwise_forces[frame];
            if (pairwise_forces.size() != other_pairwise_forces.size()) return false;

            for (size_t i = 0; i !=  pairwise_forces.size(); ++i) {
                auto pairwise_force = pairwise_forces[i];
                auto other_pairwise_force = other_pairwise_forces[i];
                if (!pairwise_force.equal(other_pairwise_force, comparer)) {
                    std::cout << "i (actual, expected) = " << pairwise_force.i << " " << other_pairwise_force.i << std::endl;
                    std::cout << "j (actual, expected) = " << pairwise_force.j << " " << other_pairwise_force.j << std::endl;
                    std::cout << "force (actual, expected) = " << pairwise_force.force << " " << other_pairwise_force.force << std::endl;
                    return false;
                }
            }
        }
        return true;
    }

private:

    friend class PairwiseForcesTest_DefaultConstructor_Test;
    friend class PairwiseForcesTest_ReadFile1_Test;
    friend class PairwiseForcesTest_ReadFile2_Test;
    friend class PairwiseForcesTest_ReadFile3_Test;

    /// Output stream
    friend std::ostream& operator << (std::ostream& os, PairwiseForces const& pf)
    {
      for (size_t i = 0; i != pf.all_pairwise_forces.size(); ++i) {
        os << "frame " << i << std::endl;
        for (auto const& e : pf.all_pairwise_forces[i]) {
            os << e.i << " "
               << e.j << " "
               << e.force << std::endl;
        }
      }
      return os;
    }

    struct PairwiseForce {

        PairwiseForce(int i, int j, ForceType force)
         : i(i), j(j), force(force)
        {}

        bool operator == (PairwiseForce const& other) const {
            return i == other.i and j == other.j and force == other.force;
        }

        bool operator != (PairwiseForce const& other) const {
            return !operator == (other);
        }

        template <class Comparer>
        bool equal(PairwiseForce const& other, Comparer const& comparer) const {
            return i == other.i and j == other.j and force.equal(other.force, comparer);
        }

        int i;
        int j;
        ForceType force;
    };

    typedef typename std::vector<PairwiseForce> PairwiseForceList;

    std::vector<PairwiseForceList> all_pairwise_forces;
};

} // namespace fda

#endif /* SRC_GROMACS_FDA_PAIRWISEFORCES_H_ */
