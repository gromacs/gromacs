/*
 * EnumParser.h
 *
 *  Created on: Apr 2, 2015
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef ENUMPARSER_H_
#define ENUMPARSER_H_

#include <map>
#include <string>

#include "gromacs/utility/fatalerror.h"

namespace fda_analysis {

/**
 * Mapping string to enum and reverse.
 * The enums are always case sensitive, whereas the correlated strings can be handled both,
 * depending on the template argument CS, which is by default true.
 */
//template <class T, bool CS = true>
template <class T>
class EnumParser
{
    typedef std::map<std::string, T> MapStringToEnum;
    typedef std::map<T, std::string> MapEnumToString;

public:

    EnumParser() {}

    T operator () (std::string const& value) const
    {
        typename MapStringToEnum::const_iterator iValue = mapStringToEnum.find(value);
        if (iValue == mapStringToEnum.end()) gmx_fatal(FARGS, "EnumParser: enum type not found %s.", value.c_str());
        return iValue->second;
    }

    std::string operator () (T value) const
    {
        typename MapEnumToString::const_iterator iValue = mapEnumToString.find(value);
        if (iValue == mapEnumToString.end()) gmx_fatal(FARGS, "EnumParser: enum type not found %i.", value);
        return iValue->second;
    }

private:

    void map(std::string const& s, T e)
    {
    	mapStringToEnum[s] = e;
    	mapEnumToString[e] = s;
    }

    MapStringToEnum mapStringToEnum;
    MapEnumToString mapEnumToString;

};

//template <class T>
//T EnumParser<T, false>::operator () (std::string const& value) const
//{
//    typename MapStringToEnum::const_iterator iValue = mapStringToEnum.find(value);
//    if (iValue == mapStringToEnum.end()) gmx_fatal(FARGS, "EnumParser: enum type not found.");
//    return iValue->second;
//}

} // namespace fda_analysis

#endif /* ENUMPARSER_H_ */
