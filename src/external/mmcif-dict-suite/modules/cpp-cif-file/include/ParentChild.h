//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ParentChild.h
**
** \brief Header file for ParentChild class.
*/


#ifndef PARENTCHILD_H
#define PARENTCHILD_H


#include <string>
#include <vector>
#include <map>

#include "ISTable.h"


class ParentChild
{
  public:
    ParentChild();
    virtual ~ParentChild();

    const std::vector<std::vector<std::string> >&
      GetComboKeys(const std::string& catName);

    std::vector<std::vector<std::vector<std::string> > >&
      GetChildrenKeys(const std::vector<std::string>& parComboKey);

    void GetParents(std::vector<std::vector<std::string> >& parParKeys,
      std::vector<std::vector<std::string> >& comboComboKeys,
      const std::string& childCat);

    void PrintAllParents(const std::string& childCat, std::vector<std::vector<std::string> >& comboComboKeys);

    void GetLinkGroupIdLabel(std::string& linkGroupIdLabel,
      const std::vector<std::string>& parKeys,
      const std::vector<std::string>& childKeys);

    bool IsParKeyPresent(const std::vector<std::string>& parKey,
      const std::string& childCatName);

    bool IsInParentComboKeys(const std::string& itemName);

    ISTable* _groupTableP;
    ISTable* _groupListTableP;

  protected:
    // Maps parent category name to its combo keys.
    std::map<std::string, std::vector<std::vector<std::string> > >
      _parComboKeys;

    // Maps parent combo keys to children combo keys.
    std::map<std::vector<std::string>,
      std::vector<std::vector<std::vector<std::string> > > > _relations;

    void GetComboKeys(const std::string& parCatName,
      const unsigned int maxKeyGroup, ISTable& keysTable, 
      std::vector<std::vector<std::string> >& comboKeys,
      std::vector<std::string>& parKeys);

    virtual void GetParentCifItems(std::vector<std::string>& parCifItems,
      const std::string& cifItemName) = 0;

    void AddParentCategoryToItemLinkedGroup(ISTable& itemLinkedGroup,
      ISTable& itemLinkedGroupList);

    void CreateAllRelations(ISTable& itemLinkedGroup,
      ISTable& itemLinkedGroupList);

    void ISTableFindPairs(std::map<std::string,
      std::vector<std::vector<std::string> > >& childrenKeys,
      const std::vector<std::string>& parKeys, ISTable& itemLinkedGroupList);

    void UpdateMap(std::map<std::string,
      std::vector<std::vector<std::string> > >& childrenKeys,
      const std::string& childCat, std::vector<std::string>& childKeys);

    void UpdateParComboKeys(const std::string& parName,
      std::vector<std::string>& parKeys);

    void UpdateRelations(std::vector<std::string>& parKeys,
      std::vector<std::vector<std::string> >& comboKeys);

    bool KeysMatch(const std::vector<std::string>& firstKey,
      const std::vector<std::string>& secondKey);
};


#endif

