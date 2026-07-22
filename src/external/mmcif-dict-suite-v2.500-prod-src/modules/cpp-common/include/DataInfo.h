//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/**
** \file DataInfo.h
**
** Data info class.
*/


#ifndef DATAINFO_H
#define DATAINFO_H


#include <string>
#include <vector>
#include <iostream>

#include "rcsb_types.h"
#include "GenString.h"


class DataInfo
{
  public:
    DataInfo();
    virtual ~DataInfo();

    virtual void GetVersion(std::string& version) = 0;

    virtual const std::vector<std::string>& GetCatNames() = 0;

    virtual const std::vector<std::string>& GetItemsNames() = 0;

    virtual bool IsCatDefined(const std::string& catName) const = 0;

    virtual bool IsItemDefined(const std::string& itemName) = 0;

    virtual const std::vector<std::string>&
      GetCatKeys(const std::string& catName) = 0;

    virtual const std::vector<std::string>&
      GetCatAttribute(const std::string& catName,
      const std::string& refCatName, const std::string& refAttribName) = 0;

    virtual const std::vector<std::string>&
      GetItemAttribute(const std::string& itemName,
      const std::string& refCatName, const std::string& refAttribName) = 0;

    virtual bool AreAllKeyItems(const std::string& catName,
      const std::vector<std::string>& attribsNames);

    virtual bool IsUnknownValueAllowed(const std::string& catName,
      const std::string& attribName);

    bool AreItemsValuesValid(const std::string& catName,
      const std::vector<std::string>& attribsNames,
      const std::vector<unsigned int>& attribsIndices,
      const std::vector<bool>& allowedNullAttribs,
      const std::vector<std::string>& values,
      const Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    virtual bool IsKeyItem(const std::string& catName,
      const std::string& attribName,
      const Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    virtual bool MustConvertItem(const std::string& catName,
      const std::string& attribName);

    virtual void GetItemsTypes(std::vector<eTypeCode>& attribsTypes,
      const std::string& catName, const std::vector<std::string>& attribsNames);

    virtual void StandardizeEnumItem(std::string& value,
      const std::string& catName,
      const std::string& attribName);

    void GetMandatoryItems(std::vector<std::string>& mandItemsNames,
      const std::string& catName);

    bool IsItemMandatory(const std::string& catName,
      const std::string& attribName);
    virtual bool IsItemMandatory(const std::string& itemName);

    // VLAD - RESOLVE THIS
    virtual bool IsSimpleDataType(const std::string& itemName);
    virtual eTypeCode _GetDataType(const std::string& itemName);

    const std::vector<std::vector<std::string> >&
      GetComboKeys(const std::string& catName);
    std::vector<std::vector<std::vector<std::string> > >&
      GetChildrenKeys(const std::vector<std::string>& parComboKey);
};

#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
extern std::string CIF_ITEM;
#endif

#endif
