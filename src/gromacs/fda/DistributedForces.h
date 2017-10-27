/*
 * DistributedForces.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_
#define SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_

#include <vector>
#include "gromacs/math/paddedvector.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "DetailedForce.h"
#include "FDASettings.h"
#include "Force.h"
#include "Vector.h"
#include "Vector2Scalar.h"

/// Forwarding needed for friend declaration
class FDA;

namespace fda {

/**
 * Storage container for distributed forces
 * Same structure for atom and residue based distribution
 *
 * TODO: the bonded interaction type might be possible to be reduced to a
 * single real vector (=eliminate jjnr) as the atoms interacting are known
 * at simulation start and do not change (also their order in the bond list
 * doesn't change). This would however require an extra step of looking up
 * the atom indices when writing out the force matrix. This would also require
 * a change in the way the atom structure is accessed, it makes no sense to
 * keep an array of t_pf_interaction items, but only an array of reals
 * representing forces.
 */
class DistributedForces
{
public:

    /// Constructor
    DistributedForces(int syslen, FDASettings const& fda_settings);

    /// Clear summed/detailed array for the next frame
    void clear();

    /// Clear scalar array
    void clear_scalar();

    void add_summed(int i, int j, Vector const& force, InteractionType type);

    void add_detailed(int i, int j, Vector const& force, PureInteractionType type);

    void write_detailed_vector(std::ostream& os) const;

    void write_detailed_scalar(std::ostream& os, PaddedRVecVector const& x) const;

    void write_summed_vector(std::ostream& os) const;

    void write_summed_scalar(std::ostream& os, PaddedRVecVector const& x) const;

    void write_scalar(std::ostream& os) const;

    void write_total_forces(std::ostream& os, PaddedRVecVector const& x) const;

    void write_scalar_compat_ascii(std::ostream& os) const;

    void write_summed_compat_ascii(std::ostream& os, PaddedRVecVector const& x) const;

    void write_scalar_compat_bin(std::ostream& os) const;

    void write_summed_compat_bin(std::ostream& os, PaddedRVecVector const& x) const;

    /// Divide all scalar forces by the divisor
    void scalar_real_divide(real divisor);

    void summed_merge_to_scalar(PaddedRVecVector const& x);

private:

    friend class ::FDA;
    template <class Base> friend class FDABase;

    /// Total number of atoms/residues in the system
    int syslen;

    /// Indices of second atom (j)
    std::vector<std::vector<int>> indices;

    /// Indices of second atom (j)
    std::vector<std::vector<int>> scalar_indices;

    /// Scalar force pairs
    std::vector<std::vector<Force<real>>> scalar;

    /// Summed vector force pairs
    std::vector<std::vector<Force<Vector>>> summed;

    /// Detailed force pairs
    std::vector<std::vector<DetailedForce>> detailed;

    /// FDA settings
    FDASettings const& fda_settings;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_DISTRIBUTEDFORCES_H_ */
