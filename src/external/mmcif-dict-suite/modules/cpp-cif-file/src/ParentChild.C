//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file ParentChild.C
**
** \brief Implementation file for CifFile class.
*/


#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include "Exceptions.h"
#include "GenCont.h"
#include "CifString.h"
#include "ISTable.h"
#include "ParentChild.h"


using std::runtime_error;
using std::find;
using std::string;
using std::pair;
using std::vector;
using std::set;
using std::multiset;
using std::map;
using std::multimap;
using std::cout;
using std::endl;


ParentChild::ParentChild()
{
    _groupTableP = new ISTable("pdbx_item_linked_group");
    _groupTableP->AddColumn("category_id"); // child's category name
    _groupTableP->AddColumn("link_group_id");  // combo key group
    _groupTableP->AddColumn("label");  // atom_label_1
    _groupTableP->AddColumn("context");  // context
    _groupTableP->AddColumn("condition_id");  // context

    _groupListTableP = new ISTable("pdbx_item_linked_group_list");
    _groupListTableP->AddColumn("child_category_id"); // category name
    _groupListTableP->AddColumn("link_group_id");  // combo key group
    _groupListTableP->AddColumn("child_name");  // item name
    _groupListTableP->AddColumn("parent_name");  // item name
    _groupListTableP->AddColumn("parent_category_id");  // category name
}


ParentChild::~ParentChild()
{
    delete (_groupListTableP);
    delete (_groupTableP);
}


const vector<vector<string> >& ParentChild::GetComboKeys(
  const string& catName)
{
    return ((_parComboKeys)[catName]);
}


vector<vector<vector<string> > >& ParentChild::GetChildrenKeys(
  const vector<string>& parComboKey)
{
    return ((_relations)[parComboKey]);
}


void ParentChild::GetComboKeys(const string& parCatName,
  const unsigned int maxKeyGroup, ISTable& keysTable,
  vector<vector<string> >& comboKeys, vector<string>& parKeys)
{
    // First determine the order of keys
    // Look for all keys that have this parent and which group number is 1
    // Return result will have keys that may have single or multipe occurences

    vector<string> keyList;
    keyList.push_back("keyGroup");
    keyList.push_back("parCategory");

    vector<string> keyTarget;
    keyTarget.push_back("1");
    keyTarget.push_back(parCatName);

    vector<string> groupList;
    groupList.push_back("keyGroup");
    groupList.push_back("parKeyCifItem");

    vector<unsigned int> groupOne;

    keysTable.Search(groupOne, keyTarget, keyList);

    vector<string> firstComboKey;

    // Fill in the first combo key and find out the order of parent keys
    for (unsigned int itemI = 0; itemI < groupOne.size(); ++itemI)
    {
        firstComboKey.push_back(keysTable(groupOne[itemI],
          "childKeyCifItem"));

        parKeys.push_back(keysTable(groupOne[itemI],
          "parKeyCifItem"));
    }

    comboKeys.push_back(firstComboKey);

    for (unsigned int groupI = 2; groupI <= maxKeyGroup; ++groupI)
    {
        vector<string> nextComboKey;

        for (unsigned int parKeyI = 0; parKeyI < parKeys.size(); ++parKeyI)
        {
            vector<string> groupTarget;
            groupTarget.push_back(String::IntToString(groupI));
            groupTarget.push_back(parKeys[parKeyI]);

            unsigned int nextGroupIndex;

            nextGroupIndex = keysTable.FindFirst(groupTarget, groupList);
            if (nextGroupIndex == keysTable.GetNumRows())
            {
                // Use it from the first combo key
                nextComboKey.push_back(firstComboKey[parKeyI]);
            }
            else
            {
                nextComboKey.push_back(keysTable(nextGroupIndex,
                  "childKeyCifItem"));
            }
        }
        comboKeys.push_back(nextComboKey);
    }
}


void ParentChild::UpdateParComboKeys(const string& parName,
  vector<string>& parKeys)
{
    if (parName.empty() || parKeys.empty())
    {
        return;
    }

    // Search for the parent
    map<string, vector<vector<string> > >::iterator pos =
      _parComboKeys.find(parName);
    if (pos != _parComboKeys.end())
    {
        // Found the parent
        // Search parent for the keys in order to avoid duplication

        vector<vector<string> >& keys = pos->second;

        bool found = false;
        for (unsigned int vecInd = 0; vecInd < keys.size(); ++vecInd)
        {
            // Order of keys is not important, so use sets for comparison
            set<string> existingKeySet(keys[vecInd].begin(),
              keys[vecInd].end());
            set<string> parKeySet(parKeys.begin(), parKeys.end());
            if (existingKeySet == parKeySet)
            {
                found = true;
                break;
            }
        }

        if (!found)
        {
            // These keys do not exist. Add them to the parent.
            (pos->second).push_back(parKeys);
        }
    }
    else
    {
        vector<vector<string> > tmp;
        tmp.push_back(parKeys);

        // Not found. Insert the parent with keys
        map<string, vector<vector<string> > >::value_type
          valuePair(parName, tmp);

        _parComboKeys.insert(valuePair);
    }
}


void ParentChild::UpdateRelations(vector<string>& parKeys,
  vector<vector<string> >& comboKeys)
{
    if (parKeys.empty() || comboKeys.empty())
    {
        return;
    }

    // Search for the parent
    map<vector<string>,
      vector<vector<vector<string> > > >::iterator pos =
      _relations.find(parKeys);
    if (pos != _relations.end())
    {
        // Found the parent keys
        // Search parent for the keys in order to avoid duplication
        if (find((pos->second).begin(), (pos->second).end(), comboKeys) ==
           (pos->second).end())
        {
            // These keys do not exist. Add them to the parent.
            (pos->second).push_back(comboKeys);
        }
    }
    else
    {
        // Not found. Insert the parent with keys
        vector<vector<vector<string> > > tmp;
        tmp.push_back(comboKeys);

        map<vector<string>,
          vector<vector<vector<string> > > >::value_type
          valuePair(parKeys, tmp);

        _relations.insert(valuePair);
    }
}


void ParentChild::AddParentCategoryToItemLinkedGroup(
  ISTable& itemLinkedGroup, ISTable& itemLinkedGroupList)
{
    vector<string> newCol;

    vector<unsigned int> rowsToDelete;

    for (unsigned int rowI = 0; rowI < itemLinkedGroup.GetNumRows(); ++rowI)
    {
        vector<string> searchCol;
        searchCol.push_back("child_category_id");
        searchCol.push_back("link_group_id");

        vector<string> searchVal;
        searchVal.push_back(itemLinkedGroup(rowI, "category_id"));
        searchVal.push_back(itemLinkedGroup(rowI, "link_group_id"));

        unsigned int found = itemLinkedGroupList.FindFirst(searchVal,
          searchCol);
        if (found == itemLinkedGroupList.GetNumRows())
        {
            cout << "Warning: Link group \"" + searchVal[1] +
              "\" of category \"" + searchVal[0] + "\" not found in table \"" +
              itemLinkedGroupList.GetName() + "\" and this entry will be "\
              "ignored." << endl;

            rowsToDelete.push_back(rowI);

            continue;
        }

        newCol.push_back(itemLinkedGroupList(found, "parent_category_id"));
    }

    itemLinkedGroup.DeleteRows(rowsToDelete);
 
    itemLinkedGroup.AddColumn("parent_category_id", newCol);
}


void ParentChild::CreateAllRelations(
  ISTable& itemLinkedGroup, ISTable& itemLinkedGroupList)
{
    for (unsigned int rowI = 0; rowI < itemLinkedGroup.GetNumRows(); ++rowI)
    {
        vector<string> searchCol;
        searchCol.push_back("child_category_id");
        searchCol.push_back("link_group_id");
        searchCol.push_back("parent_category_id");

        vector<string> searchVal;
        searchVal.push_back(itemLinkedGroup(rowI, "category_id"));
        searchVal.push_back(itemLinkedGroup(rowI, "link_group_id"));
        searchVal.push_back(itemLinkedGroup(rowI, "parent_category_id"));

        vector<unsigned int> found;
        itemLinkedGroupList.Search(found, searchVal, searchCol);

        vector<string> parKeys;
        for (unsigned int foundI = 0; foundI < found.size(); ++foundI)
        {
            parKeys.push_back(itemLinkedGroupList(found[foundI],
              "parent_name")); 
        }

#ifdef VLAD_FIX
    if (parKeys.empty())
    {
        cout << "EMPTY PARENT KEYS !!!" << endl;
        exit(1);
    }
#endif

        map<string, vector<vector<string> > > childrenKeys;

        ISTableFindPairs(childrenKeys, parKeys, itemLinkedGroupList);

        UpdateParComboKeys(itemLinkedGroup(rowI, "parent_category_id"),
          parKeys);

        for (map<string, vector<vector<string> > >::iterator pos =
          childrenKeys.begin();
          pos != childrenKeys.end(); ++pos)
        {
            // pos->first is child name
            // pos->second are vector of childs combo keys
            UpdateRelations(parKeys, pos->second);
        }
    }
}


void ParentChild::ISTableFindPairs(
  map<string, vector<vector<string> > >& childrenKeys,
  const vector<string>& parKeys, ISTable& itemLinkedGroupList)
{
    string parCatName;
    CifString::GetCategoryFromCifItem(parCatName, parKeys[0]);

    // Find all the combo keys in which the first parent key exists.
    vector<string> searchCol;
    searchCol.push_back("parent_name");

    vector<string> searchVal;
    searchVal.push_back(parKeys[0]);

    vector<unsigned int> found1;
    itemLinkedGroupList.Search(found1, searchVal, searchCol);

    for (unsigned int found1I = 0; found1I < found1.size(); ++found1I)
    {
        // See if parent key of this group matches the input parent key.
        vector<string> searchCol2;
        searchCol2.push_back("parent_category_id");
        searchCol2.push_back("link_group_id");
        searchCol2.push_back("child_category_id");

        vector<string> searchVal2;
        searchVal2.push_back(parCatName);
        searchVal2.push_back(itemLinkedGroupList(found1[found1I],
          "link_group_id"));
        searchVal2.push_back(itemLinkedGroupList(found1[found1I],
          "child_category_id"));

        vector<unsigned int> found2;
        itemLinkedGroupList.Search(found2, searchVal2, searchCol2);

        vector<string> foundParKeys;
        for (unsigned int found2I = 0; found2I < found2.size(); ++found2I)
        {
            foundParKeys.push_back(itemLinkedGroupList(found2[found2I],
              "parent_name"));
        }

        if (!KeysMatch(parKeys, foundParKeys))
            continue;

        string currChildCat = itemLinkedGroupList(found1[found1I],
          "child_category_id");

        // The position of child key elements is relevant to the parent keys
        vector<string> childKeys;

        // The multimap key is a parent key and the multimap value is
        // child key.
        multimap<string, string, StringLess> keysMap;

        for (unsigned int found2I = 0; found2I < found2.size(); ++found2I)
        {
            multimap<string, string, StringLess>::value_type
              valuePair(itemLinkedGroupList(found2[found2I],
              "parent_name"), itemLinkedGroupList(found2[found2I],
              "child_name"));
            keysMap.insert(valuePair);
        }

        // This block puts the child key to the appropriate position in
        // the vector of child keys.

        // The map key is a parent key and the map value is its order number.
        // Order number different than 0 can happen for parent keys that have
        // multiple child keys.
        map<string, unsigned int, StringLess> parKeysMap;
        for (unsigned int parKeyI = 0; parKeyI < parKeys.size(); ++parKeyI)
        {
            unsigned int orderNum = 0;

            if (!(parKeysMap.find(parKeys[parKeyI]) == parKeysMap.end()))
            {
                orderNum = parKeysMap[parKeys[parKeyI]];
            }
 
            map<string, unsigned int, StringLess>::value_type
              valuePair(parKeys[parKeyI], orderNum + 1);
            parKeysMap.insert(valuePair);

            pair<multimap<string, string, StringLess>::iterator,
              multimap<string, string, StringLess>::iterator> range;

            range = keysMap.equal_range(parKeys[parKeyI]);
            for (multimap<string, string, StringLess>::iterator it =
              range.first; it != range.second; ++it)
            {
                if (static_cast<unsigned int>(distance(range.first, it)) ==
                  orderNum)
                {
                    childKeys.push_back((*it).second);
                }
            }
        } // For every parent key

        if (parKeys.size() != childKeys.size())
            cout << "PARENT AND CHILD KEYS DIFFERENT SIZE" << endl;

        UpdateMap(childrenKeys, currChildCat, childKeys);
    } // For all combo keys that contain the first parent key item
}


void ParentChild::UpdateMap(
  map<string, vector<vector<string> > >& childrenKeys,
  const string& childCat, vector<string>& childKeys)
{
    // Search for the parent
    map<string,
      vector<vector<string> > >::iterator pos =
      childrenKeys.find(childCat);
    if (pos != childrenKeys.end())
    {
        // Found the child
        // Search child for the keys in order to avoid duplication
        if (find((pos->second).begin(), (pos->second).end(), childKeys) ==
           (pos->second).end())
        {
            // These keys do not exist. Add them to the parent.
            (pos->second).push_back(childKeys);
        }
    }
    else
    {
        // Not found. Insert the parent with keys
        vector<vector<string> > tmp;
        tmp.push_back(childKeys);

        map<string,
          vector<vector<string> > >::value_type
          valuePair(childCat, tmp);

        childrenKeys.insert(valuePair);
    }
}


void ParentChild::GetParents(vector<vector<string> >& parParKeys,
  vector<vector<string > >& comboComboKeys, const string& childCat)
{
    // For every parent with combo keys
    for (map<string, vector<vector<string> > >::iterator pos =
      _parComboKeys.begin(); pos != _parComboKeys.end(); ++pos)
    {
        vector<vector<string> >& parKeys = pos->second;

        // For every combo parent key
        for (unsigned int parKeysI = 0; parKeysI < parKeys.size(); ++parKeysI)
        {
            vector<string>& parKey = parKeys[parKeysI];

            // See if it exists in children

            for (map<vector<string>,
              vector<vector<vector<string> > > >::iterator pos2 =
              _relations.begin(); pos2 != _relations.end(); ++pos2)
            {
                if (pos2->first != parKey)
                    continue;

                vector<vector<vector<string> > >& childrenComboKeys =
                  pos2->second;

                for (unsigned int chComboI = 0; chComboI <
                   childrenComboKeys.size(); ++chComboI)
                {
                    vector<vector<string> >& childComboKeys =
                      childrenComboKeys[chComboI];
 
                    // Find the child category of these childCombo Keys
                    string childCatName;
                    CifString::GetCategoryFromCifItem(childCatName,
                      childComboKeys[0][0]);

                    if (childCatName != childCat)
                        continue;

                    for (unsigned int kI = 0; kI < childComboKeys.size(); ++kI)
                    {
                        parParKeys.push_back(parKey);
                        comboComboKeys.push_back(childComboKeys[kI]);
                    }
                    // WARNIIIIIIIIIIIIIIIIIIINIG
                    // break;
                }
            }
        } 
    }
}


void ParentChild::PrintAllParents(const string& childCat, std::vector<std::vector<std::string> >& comboComboKeys)
{
    cout << "Child: \"" << childCat << "\"" << endl;

    // For every parent with combo keys
    for (map<string, vector<vector<string> > >::iterator pos =
      _parComboKeys.begin(); pos != _parComboKeys.end(); ++pos)
    {
        vector<vector<string> >& parKeys = pos->second;

        cout << "  Parent: \"" << pos->first << "\"" << endl;

        // For every combo parent key
        for (unsigned int parKeysI = 0; parKeysI < parKeys.size(); ++parKeysI)
        {
            vector<string>& parKey = parKeys[parKeysI];
            cout << "  Parent key (";
            
            for (unsigned int keyI = 0; keyI < parKey.size(); ++keyI)
            {
                cout << parKey[keyI];
                if (keyI != (parKey.size() - 1))
                    cout << ", ";
            } 

            cout << ")" << endl;
            // See if it exists in children

            for (map<vector<string>,
              vector<vector<vector<string> > > >::iterator pos2 =
              _relations.begin(); pos2 != _relations.end(); ++pos2)
            {
                if (pos2->first != parKey)
                    continue;

                vector<vector<vector<string> > >& childrenComboKeys =
                  pos2->second;

                vector<vector<string> >& childComboKeys = childrenComboKeys[0];
 
                // Find the child category of these childCombo Keys
                string childCatName;
                CifString::GetCategoryFromCifItem(childCatName,
                  childComboKeys[0][0]);

                if (childCatName != childCat)
                    continue;

                for (unsigned int kI = 0; kI < childComboKeys.size(); ++kI)
                {        
                    comboComboKeys.push_back(childComboKeys[kI]);
                    for (unsigned int sI = 0; sI < childComboKeys[kI].size();
                      ++sI)
                    {
                        cout << "    Child key (";
                        cout << childComboKeys[kI][sI];
                        if (sI != (childComboKeys[kI].size() - 1))
                            cout << ", ";
                    }
                    cout << ")" << endl;

                }
            }
        }
    }
}


void ParentChild::GetLinkGroupIdLabel(string& linkGroupIdLabel,
  const vector<string>& parKeys, const vector<string>& childKeys)
{
    linkGroupIdLabel.clear();

    string parCatName;
    CifString::GetCategoryFromCifItem(parCatName, parKeys[0]);

    string childCatName;
    CifString::GetCategoryFromCifItem(childCatName, childKeys[0]);

    vector<string> searchCol;
    searchCol.push_back("child_category_id");
    searchCol.push_back("child_name");
    searchCol.push_back("parent_name");
    searchCol.push_back("parent_category_id");

    vector<string> searchCol2;
    searchCol2.push_back("child_category_id");
    searchCol2.push_back("link_group_id");
    searchCol2.push_back("child_name");
    searchCol2.push_back("parent_name");
    searchCol2.push_back("parent_category_id");

    vector<string> searchVal;
    searchVal.push_back(childCatName);
    searchVal.push_back(childKeys[0]);
    searchVal.push_back(parKeys[0]);
    searchVal.push_back(parCatName);

    vector<unsigned int> found1;
    _groupListTableP->Search(found1, searchVal, searchCol);

    for (unsigned int found1I = 0; found1I < found1.size(); ++found1I)
    {
        string currLinkGroupId = (*_groupListTableP)(found1[found1I],
          "link_group_id");

        bool matched = true;

        for (unsigned int parKeysI = 1; parKeysI < parKeys.size();
          ++parKeysI)
        {
            vector<string> searchVal2;
            searchVal2.push_back(childCatName);
            searchVal2.push_back(currLinkGroupId);
            searchVal2.push_back(childKeys[parKeysI]);
            searchVal2.push_back(parKeys[parKeysI]);
            searchVal2.push_back(parCatName);

            unsigned int found2 = _groupListTableP->FindFirst(searchVal2,
              searchCol2);
            if (found2 == _groupListTableP->GetNumRows())
            {
                matched = false;
                break;
            }
        }

        if (matched)
        {
            // Find the label for this linkGroupId
            vector<string> searchCol3;
            searchCol3.push_back("category_id");
            searchCol3.push_back("link_group_id");
            searchCol3.push_back("parent_category_id");

            vector<string> searchVal3;
            searchVal3.push_back(childCatName);
            searchVal3.push_back(currLinkGroupId);
            searchVal3.push_back(parCatName);

            unsigned int found3 = _groupTableP->FindFirst(searchVal3,
              searchCol3);

            if (found3 == _groupTableP->GetNumRows())
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "ParentChild::GetLinkGroupIdLabel");
            }

            linkGroupIdLabel = (*_groupTableP)(found3, "label");
            return;
        }
    }
}


bool ParentChild::IsParKeyPresent(const vector<string>& parKey,
  const string& childCatName)
{
    for (map<string, vector<vector<string> > >::iterator pos =
      _parComboKeys.begin(); pos != _parComboKeys.end(); ++pos)
    {
        vector<vector<string> >& second = pos->second;

        for (vector<vector<string> >::iterator pos2 = second.begin();
          pos2 != second.end(); ++pos2)
        {
            if (*pos2 == parKey)
            {
                // Get its children and see if they belong to this child
                // Search for the parent
                for (map<vector<string>,
                  vector<vector<vector<string> > > >::iterator pos3 =
                  _relations.begin(); pos3 != _relations.end(); ++pos3)
                {
                    const vector<string>& pK = pos3->first;
                    if (pK != parKey)
                        continue;
    
                    vector<vector<vector<string> > >& cccK = pos3->second;
    
                    for (unsigned int cccI = 0; cccI < cccK.size(); ++cccI)
                    {
                        string cCat;
                        CifString::GetCategoryFromCifItem(cCat,
                          cccK[cccI][0][0]);

                        if (childCatName == cCat)
                            return (true);
                    }
                }
            }
        }
    }

    return (false);
}


bool ParentChild::IsInParentComboKeys(const std::string& itemName)
{
    string catName;
    CifString::GetCategoryFromCifItem(catName, itemName);

    const vector<vector<string> >& parComboKeys = GetComboKeys(catName);

    for (unsigned int keyI = 0; keyI < parComboKeys.size(); ++keyI)
    {
        if (GenCont::IsInVector(itemName, parComboKeys[keyI]))
        {
            return (true);
        }
    }

    return (false);
}


bool ParentChild::KeysMatch(const vector<string>& firstKey,
  const vector<string>& secondKey)
{
    if (firstKey.size() != secondKey.size())
    {
        return (false);
    }

    multiset<string, StringLess> first(firstKey.begin(), firstKey.end());
    multiset<string, StringLess> second(secondKey.begin(), secondKey.end());
 
    if (first != second)
    {
        return (false);
    }

    return (true);
}
