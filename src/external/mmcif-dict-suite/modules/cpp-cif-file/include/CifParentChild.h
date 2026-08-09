//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifParentChild.h
**
** \brief Header file for CifParentChild class.
*/


#ifndef CIFPARENTCHILD_H
#define CIFPARENTCHILD_H


#include <string>
#include <vector>
#include <sstream>

#include "ISTable.h"
#include "TableFile.h"
#include "ParentChild.h"


#define JW_HACK

#ifdef JW_HACK
//  JW_DEBUG invokes significant debug printout associated with parent
//           child checking.  Some hacks have been introduced to improve
//            checking in this module.
// #define JW_DEBUG 1
#endif


class CifParentChild : public ParentChild
{
  public:
    CifParentChild(Block& block);
    CifParentChild(Block& block, ISTable* parChildTableP);

    virtual ~CifParentChild();

    int CheckParentChild(Block& block, ISTable& catTable,
      std::ostringstream& log);

    void WriteGroupTables(Block& block);

  protected:
    void GetParentCifItems(std::vector<std::string>& parCifItems,
      const std::string& cifItemName);

  private:
    ISTable* _parChildTableP;

    ISTable* _inParChildGroupP;
    ISTable* _inParChildGroupListP;

    void Init(Block& block);

    ISTable* CreateKeysTableOld(const std::vector<std::string>& cifItemNames,
      std::map<std::string, unsigned int>& maxKeyGroups);

    void FillKeysTableOld(ISTable& keysTable,
      const std::vector<std::string>& cifItemNames,
      std::map<std::string, unsigned int>& maxKeyGroups);

    void BuildOldTables(const std::vector<std::string>& cats,
      const std::vector<std::vector<std::string> >& items);

    void BuildNewTables(const std::vector<std::string>& cats,
      const std::vector<std::vector<std::string> >& items);

    void FilterMissingItems(std::vector<std::vector<std::string> >& parParKeys,
      std::vector<std::vector<std::string> >& comboComboKeys,
      const std::vector<std::string>& cifItemNames);

    unsigned int LastGroupNum(const std::string& childCat);
};


#endif

