//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifExcept.h
**
** \brief Header file for CifExcept class.
*/


#ifndef CIFEXCEPT_H
#define CIFEXCEPT_H


#include <string>


/**
**  \class CifExcept
**
**  \brief Static class that represents some exceptions in CIF files
**    related to data values.
**
*/
class CifExcept
{
  public:
    static bool CanBeUnknown(const std::string& itemName);
    static bool CanBeInapplicable(const std::string& itemName);
    static bool IsBadParentRelation(const std::string& itemName);
    static bool IsBadChildRelation(const std::string& itemName);
};


#endif
