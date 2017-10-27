/*
 * Tensor.h
 *
 *  Created on: Nov 21, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_TENSOR_H_
#define SRC_GROMACS_FDA_TENSOR_H_

#include "gromacs/math/vectypes.h"

namespace fda {

/// As currently no RAII algebra types are available in GROMACS
/// a wrapper class is used around tensor.
class Tensor
{
public:

    /// Default constructor, set elements to zero
    Tensor()
    {
        t[0][0] = 0.0; t[0][1] = 0.0; t[0][2] = 0.0;
        t[1][0] = 0.0; t[1][1] = 0.0; t[1][2] = 0.0;
        t[2][0] = 0.0; t[2][1] = 0.0; t[2][2] = 0.0;
    }

    /// Add and assign Tensor
    void operator += (Tensor const& other)
    {
        t[0][0] += other.t[0][0]; t[0][1] += other.t[0][1]; t[0][2] += other.t[0][2];
        t[1][0] += other.t[1][0]; t[1][1] += other.t[1][1]; t[1][2] += other.t[1][2];
        t[2][0] += other.t[2][0]; t[2][1] += other.t[2][1]; t[2][2] += other.t[2][2];
    }

    /// Add and assign depricated tensor
    void operator += (tensor const& other)
    {
        t[0][0] += other[0][0]; t[0][1] += other[0][1]; t[0][2] += other[0][2];
        t[1][0] += other[1][0]; t[1][1] += other[1][1]; t[1][2] += other[1][2];
        t[2][0] += other[2][0]; t[2][1] += other[2][1]; t[2][2] += other[2][2];
    }

    real& operator () (int x, int y)       { return t[x][y]; }
    real  operator () (int x, int y) const { return t[x][y]; }

private:

    tensor t;

};

} // namespace fda

#endif /* SRC_GROMACS_FDA_TENSOR_H_ */
