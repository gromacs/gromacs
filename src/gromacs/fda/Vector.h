/*
 * Vector.h
 *
 *  Created on: Nov 21, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_VECTOR_H_
#define SRC_GROMACS_FDA_VECTOR_H_

#include <array>
#include <iostream>
#include "gromacs/utility/real.h"

namespace fda {

/// Simple vector class
/// Workaround until GROMACS provides own RAII algebra types
class Vector
{
public:

    Vector(real value = 0.0)
    { v[0] = value; v[1] = value; v[2] = value; }

    Vector(real *v2)
    { v[0] = v2[0]; v[1] = v2[1]; v[2] = v2[2]; }

    bool operator == (Vector const& other) const
    {
      return !operator != (other);
    }

    bool operator != (Vector const& other) const
    {
        return v[0] != other.v[0] or v[1] != other.v[1] or v[2] != other.v[2];
    }

    template <class Comparer>
    bool equal(Vector const& other, Comparer const& comparer) const
    {
      return comparer(v[0], other.v[0]) and comparer(v[1], other.v[1]) and comparer(v[2], other.v[2]);
    }

    void operator += (Vector const& other)
    {
        v[0] += other[0]; v[1] += other[1]; v[2] += other[2];
    }

    void operator += (real *other)
    {
        v[0] += other[0]; v[1] += other[1]; v[2] += other[2];
    }

    real& operator [] (int p) { return v[p]; }
    real const& operator[](int p) const { return v[p]; }

    real* get_pointer() { return &v[0]; }
    real const* get_pointer() const { return &v[0]; }

private:

    /// Static array with dimension = 3
    std::array<real, 3> v;

};

/// Output stream
inline std::ostream& operator << (std::ostream& os, Vector const& v)
{
  return os << v[0] << " " << v[1] << " " << v[2];
}

/// Input stream
inline std::istream& operator >> (std::istream& is, Vector & v)
{
  return is >> v[0] >> v[1] >> v[2];
}

} // namespace fda

#endif /* SRC_GROMACS_FDA_VECTOR_H_ */
