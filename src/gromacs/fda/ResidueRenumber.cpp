/*
 * ResiduesRenumber.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include "ResidueRenumber.h"

namespace fda {

std::ostream& operator << (std::ostream& os, ResiduesRenumber r)
{
    switch(r) {
        case ResiduesRenumber::AUTO:
            return os << "auto";
        case ResiduesRenumber::DO:
            return os << "yes";
        case ResiduesRenumber::DONT:
            return os << "no";
        default:
            return os << "invalid";
    }
}

std::istream& operator >> (std::istream& is, ResiduesRenumber& r)
{
    std::string s;
    is >> s;
    std::transform(s.begin(), s.end(), s.begin(), tolower);
    if (s == "auto")
        r = ResiduesRenumber::AUTO;
    else if (s == "yes")
        r = ResiduesRenumber::DO;
    else if (s == "no")
        r = ResiduesRenumber::DONT;
    else
        throw std::runtime_error("Unknown option " + s);
    return is;
}

} // namespace fda
