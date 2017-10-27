/*
 * OnePair.cpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include "OnePair.h"

namespace fda {

std::ostream& operator << (std::ostream& os, OnePair r)
{
    switch(r) {
        case OnePair::DETAILED:
            return os << "detailed";
        case OnePair::SUMMED:
            return os << "summed";
        default:
            return os << "invalid";
    }
}

std::istream& operator >> (std::istream& is, OnePair& r)
{
    std::string s;
    is >> s;
    std::transform(s.begin(), s.end(), s.begin(), tolower);
    if (s == "detailed")
        r = OnePair::DETAILED;
    else if (s == "summed")
        r = OnePair::SUMMED;
    else
        throw std::runtime_error("Unknown option " + s);
    return is;
}

} // namespace fda
