//
// Created by sebkelle on 20.11.19.
//

#ifndef GROMACS_ATOMS_H
#define GROMACS_ATOMS_H

#include <tuple>
#include <unordered_map>
#include <string>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "interactions.h"

class TopologyBuilder;

namespace nblib
{

using AtomName = std::string;
using Mass     = real;
using C6       = real;
using C12      = real;

class AtomType
{
public:
    AtomType() noexcept;

    //! Constructor with explicit type specification
    AtomType(AtomName atomName, Mass mass, C6 c6, C12 c12);

    //! Force explicit use of correct types
    template<typename T, typename U, typename V, typename W>
    AtomType(T atomName, U mass, V c6, W c12) = delete;

    bool operator==(const AtomType& b);

    //! Get the name
    AtomName name() const;

    //! Get the mass
    Mass mass() const;

    //! Get the c6 param
    C6 c6() const;

    //! Get the c12 param
    C12 c12() const;

private:
    //! The name
    AtomName name_;
    //! The mass
    Mass mass_;
    //! The c12 param
    C6 c6_;
    //! The c12 param
    C12 c12_;
};

} // namespace nblib
#endif // GROMACS_MOLECULES_H
