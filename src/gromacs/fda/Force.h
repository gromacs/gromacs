/*
 * Force.h
 *
 *  Created on: Dec 7, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FORCE_H_
#define SRC_GROMACS_FDA_FORCE_H_

#include <iostream>
#include <iomanip>
#include "gromacs/utility/real.h"
#include "InteractionType.h"
#include "Vector.h"

namespace fda {

template <typename T>
struct Force
{
    Force(T force = 0.0, InteractionType type = 0)
     : force(force), type(type)
    {}

    bool operator == (Force const& other) const {
        return force == other.force and type == other.type;
    }

    bool operator != (Force const& other) const {
        return !operator == (other);
    }

    void operator += (Force const& other) {
        force += other.force;
        type |= other.type;
    }

    template <class Comparer>
    bool equal(Force const& other, Comparer const& comparer) const;

    /// Generic force type: real or Vector
    T force;

    /// Interaction type
    InteractionType type;
};

template <>
template <class Comparer>
bool Force<real>::equal(Force<real> const& other, Comparer const& comparer) const {
    return type == other.type and comparer(force, other.force);
}

template <>
template <class Comparer>
bool Force<Vector>::equal(Force<Vector> const& other, Comparer const& comparer) const {
    return type == other.type and force.equal(other.force, comparer);
}

/// Output stream for ResultType
template <typename T>
std::ostream& operator << (std::ostream& os, Force<T> const& f)
{
    return os << std::scientific << std::setprecision(6)
              << f.force << " " << f.type;
}

/// Input stream for ResultType
template <typename T>
std::istream& operator >> (std::istream& is, Force<T> & f)
{
    return is >> f.force >> f.type;
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_FORCE_H_ */
