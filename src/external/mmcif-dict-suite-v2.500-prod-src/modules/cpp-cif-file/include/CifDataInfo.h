//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/**
** \file CifDataInfo.h
**
** CIF data information class
*/


#ifndef CIFDATAINFO_H
#define CIFDATAINFO_H


#include <string>
#include <vector>

#include "DataInfo.h"
#include "DicFile.h"


class CifDataInfo : public DataInfo
{
  public:
    CifDataInfo(DicFile& dictFile);
    ~CifDataInfo();

    void GetVersion(std::string& version);

    const std::vector<std::string>& GetCatNames();

    const std::vector<std::string>& GetItemsNames();

    bool IsCatDefined(const std::string& catName) const;

    bool IsItemDefined(const std::string& itemName);

    const std::vector<std::string>& GetCatKeys(const std::string& catName);

    const std::vector<std::string>& GetCatAttribute(const std::string& catName,
      const std::string& refCatName, const std::string& refAttrName);

    const std::vector<std::string>&
      GetItemAttribute(const std::string& itemName,
      const std::string& refCatName, const std::string& refAttrName);

    virtual void GetCatItemsNames(std::vector<std::string>& itemsNames,
      const std::string& catName);

  protected:
    DicFile& _dictFile;

  private:
    std::string _version;
    std::vector<std::string> _catsNames;
    std::vector<std::string> _itemsNames;
    std::vector<std::string> _catKeyItems;
    std::vector<std::string> _catAttrib;
    std::vector<std::string> _itemAttrib;
    std::vector<std::string> _itemTypeListAttrib;

    void _GetDictVersion(std::string& dictVer);
    bool _isDictCategory(const std::string& category) const;

    const std::vector<std::string>&
      GetItemAttributeForItemTypeListCat(const std::string& itemName,
      const std::string& refCatName,
      const std::string& refAttrName);
};


#endif
