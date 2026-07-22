//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifParentChild.C
**
** \brief Implementation file for CifFile class.
*/


#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "GenString.h"
#include "GenCont.h"
#include "CifString.h"
#include "ISTable.h"
#include "TableFile.h"
#include "CifParentChild.h"


using std::string;
using std::vector;
using std::map;
using std::ostringstream;
using std::cout;
using std::endl;
using std::sort;


CifParentChild::CifParentChild(Block& block) : ParentChild(),
  _parChildTableP(NULL), _inParChildGroupP(NULL), _inParChildGroupListP(NULL)
{
    _parChildTableP = block.GetTablePtr("item_linked");

    _inParChildGroupP = block.GetTablePtr("pdbx_item_linked_group");
    _inParChildGroupListP = block.GetTablePtr("pdbx_item_linked_group_list");

    if ((_inParChildGroupP == NULL) && (_inParChildGroupListP != NULL))
    {
        throw EmptyValueException("Empty item linked group table",
          "CifParentChild::CifParentChild");
    }

    if ((_inParChildGroupP != NULL) && (_inParChildGroupListP == NULL))
    {
        throw EmptyValueException("Empty item linked group list table",
          "CifParentChild::CifParentChild");
    }

    if (_inParChildGroupP != NULL)
        *_groupTableP = *_inParChildGroupP;
    if (_inParChildGroupListP != NULL)
        *_groupListTableP = *_inParChildGroupListP;

    Init(block);
}


CifParentChild::CifParentChild(Block& block, ISTable* parChildTableP) :
  ParentChild(), _parChildTableP(parChildTableP)
{
    // This method is used when it is needed to build group and grup list
    // tables out of item linked tables.

    ISTable* categoryP = block.GetTablePtr("category");
    ISTable* itemP = block.GetTablePtr("item");

    vector<string> cats;
    for (unsigned int rowI = 0; rowI < categoryP->GetNumRows(); ++rowI)
    {
        cats.push_back((*categoryP)(rowI, "id"));
    }

    sort(cats.begin(), cats.end());

    vector<vector<string> > items;

    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        vector<string> searchValue;
        searchValue.push_back(cats[catI]);
   
        vector<string> searchCols;
        searchCols.push_back("category_id");
  
        vector<unsigned int> found;
        itemP->Search(found, searchValue, searchCols);

        vector<string> cifItemNames;
        for (unsigned int foundI = 0; foundI < found.size(); ++foundI)
        {
           cifItemNames.push_back((*itemP)(found[foundI], "name"));
        } 

        items.push_back(cifItemNames);
    }

    BuildOldTables(cats, items);

    AddParentCategoryToItemLinkedGroup(*_groupTableP, *_groupListTableP);

    CreateAllRelations(*_groupTableP, *_groupListTableP);
}


void CifParentChild::Init(Block& block)
{
    ISTable* categoryP = block.GetTablePtr("category");
    ISTable* itemP = block.GetTablePtr("item");

    vector<string> cats;
    for (unsigned int rowI = 0; rowI < categoryP->GetNumRows(); ++rowI)
    {
        cats.push_back((*categoryP)(rowI, "id"));
    }

    sort(cats.begin(), cats.end());

    vector<vector<string> > items;

    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        vector<string> searchValue;
        searchValue.push_back(cats[catI]);
   
        vector<string> searchCols;
        searchCols.push_back("category_id");
  
        vector<unsigned int> found;
        itemP->Search(found, searchValue, searchCols);

        vector<string> cifItemNames;
        for (unsigned int foundI = 0; foundI < found.size(); ++foundI)
        {
           cifItemNames.push_back((*itemP)(found[foundI], "name"));
        } 

        items.push_back(cifItemNames);
    }

    AddParentCategoryToItemLinkedGroup(*_groupTableP, *_groupListTableP);

    BuildNewTables(cats, items);

    CreateAllRelations(*_groupTableP, *_groupListTableP);
}


CifParentChild::~CifParentChild()
{

}


void CifParentChild::BuildOldTables(const vector<string>& cats, 
  const vector<vector<string> >& items)
{
    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        unsigned int groupNum = 1;

        const vector<string>& cifItemNames = items[catI];

        map<string, unsigned int> maxKeyGroups;
        ISTable* keysTableP = CreateKeysTableOld(cifItemNames, maxKeyGroups);

        if (keysTableP->GetNumRows() == 0)
        {
            // No parents for this child
            delete (keysTableP);
            continue;
        }

        for (map<string, unsigned int>::iterator pos = maxKeyGroups.begin();
          pos != maxKeyGroups.end(); ++pos)
        {
            // pos->first is parent name
            // pos->second is number of parent's groups
            for (unsigned int groupI = 0; groupI < pos->second; ++groupI,
              ++groupNum)
            {
                vector<string> newGroupRow;
                newGroupRow.push_back(cats[catI]);
                newGroupRow.push_back(String::IntToString(groupNum));
                newGroupRow.push_back(cats[catI] + ":" + pos->first + ":" +
                  String::IntToString(groupNum));
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(CifString::InapplicableValue);
                _groupTableP->AddRow(newGroupRow);

                vector<string> searchCol;
                searchCol.push_back("keyGroup");
                searchCol.push_back("parCategory");

                vector<string> searchVal;
                searchVal.push_back(String::IntToString(groupI + 1));
                searchVal.push_back(pos->first);

                vector<unsigned int> found;
                keysTableP->Search(found, searchVal, searchCol);

                for (unsigned int foundI = 0; foundI < found.size(); ++foundI)
                {
                    vector<string> newGroupListRow;
                    newGroupListRow.push_back(cats[catI]);
                    newGroupListRow.\
                      push_back(String::IntToString(groupNum));
                    newGroupListRow.push_back\
                      ((*keysTableP)(found[foundI], "childKeyCifItem"));
                    newGroupListRow.push_back\
                      ((*keysTableP)(found[foundI], "parKeyCifItem"));
                    newGroupListRow.push_back(pos->first);
                    _groupListTableP->AddRow(newGroupListRow);
                } // (for every parent/child pair in a group)
            } // for (all group numbers in a parent)
        } // for (all parents)

        delete (keysTableP); 
    } // for (all categories in the dictionary, acting as child categories)
}


void CifParentChild::BuildNewTables(const vector<string>& cats, 
  const vector<vector<string> >& items)
{
    if (_parChildTableP == NULL)
    {
        return;
    }

    // Maps a pair of child category/parent category to the last group
    // number used.

    map<string, unsigned int> maxGroupNum;

    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        const vector<string>& cifItemNames = items[catI];

        vector<string> childSearchCol;
        childSearchCol.push_back("child_name");

        for (unsigned int childI = 0; childI < cifItemNames.size(); ++childI)
        {
            // Check if item has already been processed
            vector<string> childTarget;
            childTarget.push_back(cifItemNames[childI]);

            if (_groupListTableP != NULL)
            {
                unsigned int row =
                  _groupListTableP->FindFirst(childTarget,
                  childSearchCol);

                if (row != _groupListTableP->GetNumRows())
                    continue;
            }

            vector<unsigned int> parLoc;
            _parChildTableP->Search(parLoc, childTarget, childSearchCol);

            for (unsigned int parI = 0; parI < parLoc.size(); ++parI)
            {
                const string& parentItem = (*_parChildTableP)(parLoc[parI],
                  "parent_name");

                string parCatName;
                CifString::GetCategoryFromCifItem(parCatName, parentItem);

                if (maxGroupNum[cats[catI]] == 0)
                {
                    unsigned int lastGroupNum = LastGroupNum(cats[catI]);
                    maxGroupNum[cats[catI]] = lastGroupNum + 1;
                }
                else
                    maxGroupNum[cats[catI]] += 1;

                cout << "Info: Creating a new group \"" <<
                  String::IntToString(maxGroupNum[cats[catI]]) <<
                  "\" for child category \"" << cats[catI] <<
                  "\" for child item \"" << cifItemNames[childI] <<
                  "\" and parent item \"" << parentItem <<
                  "\", from \"item_linked\" table, as these "\
                  "are not defined in group tables." << endl;

                vector<string> newGroupRow;
                newGroupRow.push_back(cats[catI]);
                newGroupRow.push_back(String::IntToString(
                  maxGroupNum[cats[catI]]));
                newGroupRow.push_back(cats[catI] + ":" + parCatName + ":" +
                  String::IntToString(maxGroupNum[cats[catI]]));
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(parCatName);
                _groupTableP->AddRow(newGroupRow);

                vector<string> row;
                row.push_back(cats[catI]);
                row.push_back(String::IntToString(
                  maxGroupNum[cats[catI]]));
                row.push_back(cifItemNames[childI]);
                row.push_back(parentItem);
                row.push_back(parCatName);
                _groupListTableP->AddRow(row);
            } // for (all child's parent items)
        } // for (all child items)
    }
}


int CifParentChild::CheckParentChild(Block& block, ISTable& catTable,
  ostringstream& log)
{
    int ret = 1;

    const string& childCatName = catTable.GetName();

    vector<vector<string> > parParKeys; 
    vector<vector<string> > comboComboKeys;

    // DEBUG PrintAllParents(parParKeys, comboComboKeys, childCatName);
    GetParents(parParKeys, comboComboKeys, childCatName);


    if (parParKeys.empty())
    {
#ifdef JW_DEBUG
        cout << "No parents to check for category " <<  catTable.GetName() << endl;
#endif
        // No parents. Just return.
        return(ret);
    }

    bool has_multiple_model = false;
    if (childCatName == "pdbx_entity_nonpoly") {
         ISTable* t_ptr = block.GetTablePtr("atom_site");
         if (t_ptr && (t_ptr->GetNumRows() > 1) && t_ptr->IsColumnPresent("pdbx_PDB_model_num")) {
              const string& cifItem_first = (*t_ptr)(0, "pdbx_PDB_model_num");
              const string& cifItem_last = (*t_ptr)(t_ptr->GetNumRows() - 1, "pdbx_PDB_model_num");
              if (cifItem_first != cifItem_last) has_multiple_model = true;
         }
    }

    bool has_linear_polymer = false;
    ISTable* entity_ptr = block.GetTablePtr("entity");
    if (entity_ptr && entity_ptr->IsColumnPresent("type")) {
         for (unsigned int rowI = 0; rowI < entity_ptr->GetNumRows(); ++rowI) {
              const string& cifItem = (*entity_ptr)(rowI, "type");
              if (cifItem == "polymer") {
                   has_linear_polymer = true;
                   break;
              }
         }
    }
    std::string first_model_num = "";
    if ((catTable.GetName() == "atom_site") && catTable.IsColumnPresent("pdbx_PDB_model_num")) {
         first_model_num = catTable(0, "pdbx_PDB_model_num");
    }

    vector<string> cifItemNames;
    CifString::MakeCifItems(cifItemNames, catTable.GetName(),
      catTable.GetColumnNames());

    FilterMissingItems(parParKeys, comboComboKeys, cifItemNames);

    for (unsigned int allParI = 0; allParI < parParKeys.size(); ++allParI)
    {
        vector<string>& parKeys = parParKeys[allParI];
        vector<string>& comboKeys = comboComboKeys[allParI];

#ifdef VLAD_DEBUG
        cout << "  parKeys.size: " << parKeys.size() << endl;
        cout << "  comboComboKeys.size: " << comboComboKeys.size() << endl;
        cout << "  parKeys (";
        for (unsigned int aaaI = 0; aaaI < parKeys.size(); ++aaaI)
        {
            cout << parKeys[aaaI];
            if (aaaI != (parKeys.size() - 1))
                cout << ", ";
        }
        cout << ")" << endl;
#endif

        if (parKeys.empty())
        {
            // No parent keys
            continue;
        }

        string parCatName;
        CifString::GetCategoryFromCifItem(parCatName, parKeys[0]);

        if (has_multiple_model && (catTable.GetName() == "pdbx_entity_nonpoly") && (parCatName == "pdbx_nonpoly_scheme")) continue;

#if JW_DEBUG
        cout << "Current category "    <<  catTable.GetName() << endl;
        cout << "Parent category  "    << parCatName << endl;
        cout << "Max group number "    << pos->second << endl;
        //
        for ( unsigned int ii =0; ii < parKeys.size(); ii++)
        {
            cout << "parent  # ii " << ii << " " << parKeys[ii] << endl;
        }

        for ( unsigned int ii =0; ii < comboKeys.size(); ii++)
        {
            for (unsigned int jj =0; jj < comboKeys[ii].size(); jj++)
            {
                cout << " comboKeys[][] " << ii << " " << jj <<  " child = "
                 << comboKeys[ii][jj] << endl;
            }
        }
#endif

        vector<string> parKeyItems;
        for (unsigned int parKeyItemsI = 0; parKeyItemsI < parKeys.size();
          ++parKeyItemsI)
        {
            string tmpItem;
            CifString::GetItemFromCifItem(tmpItem, parKeys[parKeyItemsI]);
            parKeyItems.push_back(tmpItem);
        }

#ifdef VLAD_PERF
        cout << "Processing parent \"" << parCatName << "\" combo key: " <<
          comboKeyI << " of " << comboKeys.size() << endl;
#endif
#ifdef JW_DEBUG
        cout << "Processing parent \"" << parCatName << "\" combo key: "
             << comboKeyI + 1 << " of " << comboKeys.size() <<  endl;

        cout << "In category " << catTable.GetName() << " row count  "
             <<  catTable.GetNumRows() << endl;
#endif
#ifdef JW_HACK
        int iKeyDif;
        bool allowMissing = false;
#endif
        for (unsigned int l = 0; l < catTable.GetNumRows(); ++l)
        {
#ifdef JW_DEBUG
            cout << "++ In category " << catTable.GetName() << " row: " << l
              << " of " << catTable.GetNumRows() << endl;
#endif
#ifdef VLAD_PERF
            cout << "  Processing combo key row: " << l << " of " <<
              catTable.GetNumRows() << endl;
#endif
            vector<string> parKeyNonEmptyValues;
            vector<string> parKeyAllValues;
            vector<string> parKeyNonEmptyItems;
            vector<string> childKeyNonEmptyItems;
#ifdef JW_HACK
            vector<string> skippedParKeyItems;
            iKeyDif=0;
            allowMissing=false;
#endif
            for (unsigned int keyI = 0; keyI < comboKeys.size(); ++keyI)
            {
                string itemName;
                CifString::GetItemFromCifItem(itemName, comboKeys[keyI]);

                // loop for all rows
                const string& parKeyValue = catTable(l, itemName);
                parKeyAllValues.push_back(parKeyValue);
                if (CifString::IsEmptyValue(parKeyValue))
                {
#ifdef JW_HACK
                    // Skip search for null/na/missing values
                    skippedParKeyItems.push_back(parKeyItems[keyI]);
#endif
                }
                else
                {
                    parKeyNonEmptyValues.push_back(parKeyValue);
                    parKeyNonEmptyItems.push_back(parKeyItems[keyI]);
                    childKeyNonEmptyItems.push_back(comboKeys[keyI]);
                }
            } // For every key item in a combo key

            if (parKeyNonEmptyValues.empty())
                continue;

            bool check_flag = true;
            if ((catTable.GetName() == "atom_site") && ((parCatName == "entity_poly_seq") || (parCatName == "pdbx_poly_seq_scheme")) &&
                !has_linear_polymer) check_flag = false;

#ifndef VLAD_FIX_SINGLE_KEY_UNKNOWN_VALUE_NOT_REPORT_PARENT
            if (!block.IsTablePresent(parCatName))
            {
                // Parent table does not exist. Ignore this fact.
#ifdef JW_HACK
                if ((parCatName != "chem_comp_atom") && (parCatName != "atom_sites_alt") && check_flag)
                {
                   log << "++ERROR - In block \"" << block.GetName() <<
                     "\", parent category \"" << parCatName <<
                     "\", of category \"" << catTable.GetName() <<
                     "\", is missing." << endl;
                }
#endif
                break;
            }

            ISTable* parentCatTableP = block.GetTablePtr(parCatName);
            if (parentCatTableP->GetNumRows() == 0)
            {
#ifdef JW_HACK
                if (check_flag) log << "WARN - In block " << block.GetName() << " category " <<
                             catTable.GetName() << " parent table empty " << parCatName << endl;
#endif
                 // Parent table has no rows. Ignore this fact.
                 break;
            }

            bool check_row_flag = true;
            if ((catTable.GetName() == "atom_site") && catTable.IsColumnPresent("pdbx_PDB_model_num")) {
                 const string& cifItem = catTable(l, "pdbx_PDB_model_num");
                 if (cifItem != first_model_num) check_row_flag = false;
            }

            // Check existence in the parent
            for (unsigned int parKeyI = 0; parKeyI < parKeyItems.size();
              ++parKeyI)
            {
               if (!CifString::IsEmptyValue(parKeyAllValues[parKeyI]) &&
                 !parentCatTableP->IsColumnPresent(parKeyItems[parKeyI]))
               {
                   log << "ERROR - In block \"" << block.GetName() <<
                     "\", for child item \"" <<
                     comboKeys[parKeyI] <<
                     "\", parent key attribute \"" <<
                     parKeyItems[parKeyI] <<
                     "\" is not defined in parent category \"" <<
                     parCatName << "\"" << endl;
                   ret = 0;

                   return(ret);
               }
           }
#endif

#ifdef JW_HACK
            iKeyDif = comboKeys.size() - parKeyNonEmptyValues.size();
            //
            // HACK - to deal with real missing values - 
            //
            if (iKeyDif == 1)
            {
                if (
		      ((catTable.GetName() == "atom_site") && (parCatName == "entity_poly_seq") && (skippedParKeyItems[0] == "num")) ||
		      ((catTable.GetName() == "atom_site") && (parCatName == "pdbx_poly_seq_scheme") && (skippedParKeyItems[0] == "seq_id")) 

		      )
                {
                    allowMissing = true;
		}

#ifdef JW_DEBUG
		if (allowMissing)
                {
                    cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			 << " parent " << parCatName 
			 << " skipped parent item(s) " << skippedParKeyItems[0]
			 << " in child " 
			 << catTable.GetName() << " row: " << l+1 << " of " 
			 << catTable.GetNumRows() << endl;
		  }
#endif

            }

            if (iKeyDif == 2)
            {
                  if ((catTable.GetName() == "pdbx_poly_seq_scheme") && (parCatName == "atom_site") &&
		      ((skippedParKeyItems[0] == "auth_comp_id" &&  skippedParKeyItems[1] == "pdbx_PDB_ins_code") ||
		       (skippedParKeyItems[1] == "auth_comp_id" &&   skippedParKeyItems[0] == "pdbx_PDB_ins_code")))
		    {
		      allowMissing = true;
		    } else if ((catTable.GetName() == "atom_site") && (parCatName == "pdbx_poly_seq_scheme") &&
		      ((skippedParKeyItems[0] == "seq_id" &&   skippedParKeyItems[1] == "pdb_ins_code") ||
		       (skippedParKeyItems[1] == "seq_id" &&   skippedParKeyItems[0] == "pdb_ins_code"))) 
		    {
		      allowMissing = true;
		    }
		  
#ifdef JW_DEBUG
		    if (allowMissing) {
		      cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			   << " parent " << parCatName
			   << " skipped parent item(s) " << skippedParKeyItems[0]
			   << " and " << skippedParKeyItems[1]
			   << " in child " 
			   << catTable.GetName() << " row: " << l+1 << " of " 
			   << catTable.GetNumRows() << endl;
		    }
#endif
            }

            if (iKeyDif == 3)
            { 
                  if ( (catTable.GetName() == "pdbx_poly_seq_scheme") && (parCatName == "atom_site") &&
		       ((skippedParKeyItems[0] == "auth_seq_id" &&
			 skippedParKeyItems[1] == "auth_comp_id" &&
			 skippedParKeyItems[2] == "pdbx_PDB_ins_code") ||
			(skippedParKeyItems[0] == "auth_seq_id"  &&
			 skippedParKeyItems[1] == "auth_comp_id" &&
			 skippedParKeyItems[2] == "pdbx_PDB_ins_code")))  {
		    allowMissing = true;
		  } else if ((catTable.GetName() == "struct_ref_seq_dif") && (parCatName == "pdbx_poly_seq_scheme") &&
		       ((skippedParKeyItems[0] == "mon_id" &&
			 skippedParKeyItems[1] == "seq_id" &&
			 skippedParKeyItems[2] == "pdb_ins_code") ||
			(skippedParKeyItems[1] == "mon_id"  &&
			 skippedParKeyItems[0] == "seq_id" &&
			 skippedParKeyItems[2] == "pdb_ins_code")))  {
		    allowMissing = true;
		  }
#ifdef JW_DEBUG
		  if (allowMissing) {		   
                    cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			 << " parent " << parCatName 
			 << " skipped parent item(s) " << skippedParKeyItems[0]
			 << " and " << skippedParKeyItems[1]
			 << " and " << skippedParKeyItems[2]
			 << " in child " 
			 << catTable.GetName() << " row: " << l+1 << " of " 
			 << catTable.GetNumRows() << endl;
		  }
#endif
            }

            // END HACK -


#ifdef JW_DEBUG
            if (iKeyDif != 0 && !allowMissing)
            {
                cout << "++INFO - child category (keydif=" << iKeyDif << ") " << catTable.GetName() <<
                      " row: " << l << " of " << catTable.GetNumRows() << endl;
            }

            if (iKeyDif > 0)
            {
                if (allowMissing)
                {
                    for ( unsigned int ii =0; ii < skippedParKeyItems.size(); ii++)
                    {
                        cout <<
                          "++INFO - Allowed incomplete parent key in " << parCatName << " "
                          << ii+1 << " of " << skippedParKeyItems.size() <<
                          " -> " << skippedParKeyItems[ii] << endl;
                    }
                }
                else
                {
                    for (unsigned int ii =0;
                      ii < skippedParKeyItems.size(); ii++)
                    {
                        cout << "++INFO - Incomplete parent key in " 
                          << parCatName << " " 
                          << ii+1 << " of " << skippedParKeyItems.size() <<
                          " -> " << skippedParKeyItems[ii] << endl;
                    }
                }
            }
#endif
#endif

#ifdef JW_HACK
            if ((iKeyDif == 0) || ((iKeyDif > 0) && !allowMissing))
#endif

            {
                    unsigned int searchIndex =
                      parentCatTableP->FindFirst(parKeyNonEmptyValues,
                        parKeyNonEmptyItems);
                    if ((searchIndex == parentCatTableP->GetNumRows()) && check_row_flag)
                    {
#ifdef OLD_IMPL
                        log << "ERROR - In block \"" << block.GetName()
                          << " unmatched value in " << catTable.GetName()
                          << " row " << l+1 <<  " in parent category " <<
                          parCatName  << endl;

                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "  \"" << childKeyNonEmptyItems[keyValI]
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
                              parKeyNonEmptyValues[keyValI] << "\"" << endl;
                        }
#else
                        string linkGroupIdLabel;
                        GetLinkGroupIdLabel(linkGroupIdLabel, parKeys,
                          comboKeys);
                        if (linkGroupIdLabel.empty())
                        {
                            log << "BIG TROUBLE; linkGroupId not found!" <<
                              endl;
                        }

                        log << "ERROR PARCHILD \"" << linkGroupIdLabel <<
                          "\" - In block \"" << block.GetName() <<
                          "\", in category \"" << childCatName <<
                          "\", in row# " << l + 1 <<
                          ", unmatched value in the parent \"" << parCatName <<
                          "\"";

                        log << endl;

#ifdef VLAD_MORE_CORRECTIONS
                        log <<  "(";
    
                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "\"" << childKeyNonEmptyItems[keyValI] 
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value = \"" <<
                              parKeyNonEmptyValues[keyValI] << "\"";
                            if (keyValI != parKeyNonEmptyValues.size() - 1)
                                log << ", ";
                        }

                        log << ")" << endl;
#endif

                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "  \"" << childKeyNonEmptyItems[keyValI]
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
                              parKeyNonEmptyValues[keyValI] << "\"" << endl;
                        }
#endif


#ifdef JW_DEBUG
			cout << "ERROR - In block \"" << block.GetName() 
			     << " unsatisfied child value in " << catTable.GetName() 
			     << " row " << l+1 <<  " parent category " << parCatName 
			     << endl;
			
			cout << "ERROR - full key \"" << comboKeys.size()
			     << " realized key size " <<
			  parKeyNonEmptyValues.size() << endl;
			
			for (unsigned int keyValI = 0; keyValI < parKeyNonEmptyValues.size(); ++keyValI)
			  {
			    cout << "  \"" << childKeyNonEmptyItems[keyValI] 
				 << "\"->\"" << parCatName << "." <<
			      parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
			      parKeyNonEmptyValues[keyValI] << "\"" << endl;
			  }

#endif // JW_DEBUG
                    }
                    ret = 0;
            }
        } // For every row in child
    } // For every parent

    return(ret);
}


void CifParentChild::WriteGroupTables(Block& block)
{
    ISTable* tmpGroupTableP = new ISTable();
    *tmpGroupTableP = *_groupTableP;

    ISTable* tmpGroupListTableP = new ISTable();
    *tmpGroupListTableP = *_groupListTableP;

    tmpGroupTableP->DeleteColumn("parent_category_id");

    if (tmpGroupTableP->GetNumRows() != 0)
        block.WriteTable(tmpGroupTableP);
    if (tmpGroupListTableP->GetNumRows() != 0)
        block.WriteTable(tmpGroupListTableP);
}


ISTable* CifParentChild::CreateKeysTableOld(const vector<string>& cifItemNames,
  map<string, unsigned int>& maxKeyGroups)
{
    // Prepare the table that will contain all the needed keys information
    // This table will have four columns: child key CIF item, group number
    // parent key CIF item and parent category name

    ISTable* keysTableP = new ISTable();

    keysTableP->AddColumn("childKeyCifItem");
    keysTableP->AddColumn("keyGroup");
    keysTableP->AddColumn("parKeyCifItem");
    keysTableP->AddColumn("parCategory");

    FillKeysTableOld(*keysTableP, cifItemNames, maxKeyGroups);

    return(keysTableP);
}


void CifParentChild::FillKeysTableOld(ISTable& keysTable,
  const vector<string>& cifItemNames, map<string, unsigned int>& maxKeyGroups)
{
    vector<string> keyList;
    keyList.push_back("parKeyCifItem");

    for (unsigned int cifItemI = 0; cifItemI < cifItemNames.size(); ++cifItemI)
    {
        // cifItemNames are child's items
        vector<string> parCifItems;
        GetParentCifItems(parCifItems, cifItemNames[cifItemI]);

        sort(parCifItems.begin(), parCifItems.end());

        for (unsigned int parI = 0; parI < parCifItems.size(); ++parI)
        {
            string parCatName;
            CifString::GetCategoryFromCifItem(parCatName, parCifItems[parI]);

            vector<string> keyTarget;
            keyTarget.push_back(parCifItems[parI]);

            vector<unsigned int> parents;
            keysTable.Search(parents, keyTarget, keyList);

            unsigned int newKeyGroupInt = 1;
            if (!parents.empty())
            {
                newKeyGroupInt = String::StringToInt(keysTable(
                  parents[parents.size() - 1], "keyGroup"));
                ++newKeyGroupInt;
            }

            if (newKeyGroupInt > maxKeyGroups[parCatName])
                maxKeyGroups[parCatName] = newKeyGroupInt;

            string newKeyGroup = String::IntToString(newKeyGroupInt);

            vector<string> row;

            row.push_back(cifItemNames[cifItemI]);
            row.push_back(newKeyGroup);
            row.push_back(parCifItems[parI]);
            row.push_back(parCatName);

            keysTable.AddRow(row);
        } // For all item's parent items (that can be in different parents)
    } // For all items in the category
}


void CifParentChild::GetParentCifItems(vector<string>& parCifItems,
  const string& cifItemName)
{
    parCifItems.clear();

    vector<string> childCifItem;
    childCifItem.push_back(cifItemName);

    vector<string> childNameCol;
    childNameCol.push_back("child_name");

    vector<unsigned int> parLoc;
    _parChildTableP->Search(parLoc, childCifItem, childNameCol);

    for (unsigned int parLocI = 0; parLocI < parLoc.size(); ++parLocI)
    {
        parCifItems.push_back((*_parChildTableP)(parLoc[parLocI],
          "parent_name"));
    }
}


unsigned int CifParentChild::LastGroupNum(const string& childCat)
{
    vector<string> searchCol;
    searchCol.push_back("category_id");

    vector<string> searchVal;
    searchVal.push_back(childCat);

    vector<unsigned int> found;

    unsigned int lastGroupNum = 0;

    _groupTableP->Search(found, searchVal, searchCol);

    for (unsigned int foundI = 0; foundI < found.size(); ++foundI)
    {
        unsigned int currGroupNum =
          String::StringToInt((*_groupTableP)(found[foundI],
          "link_group_id"));

        if (currGroupNum > lastGroupNum)
            lastGroupNum = currGroupNum;
    }

    return (lastGroupNum);
}


void CifParentChild::FilterMissingItems(vector<vector<string> >& parParKeys,
  vector<vector<string> >& comboComboKeys, const vector<string>& cifItemNames)
{
    // Identify missing items indices

    for (unsigned int allParI = 0; allParI < parParKeys.size(); ++allParI)
    {
        vector<string>& parKeys = parParKeys[allParI];
        vector<string>& comboKeys = comboComboKeys[allParI];

        vector<unsigned int> remIndices;
        for (unsigned int keysI = 0; keysI < comboKeys.size(); ++keysI)
        {
            if (!GenCont::IsInVector(comboKeys[keysI], cifItemNames))
                remIndices.push_back(keysI);
        }
 
        sort (remIndices.begin(), remIndices.end());

        for (unsigned int indI = 0; indI < remIndices.size(); ++indI)
        {
            // Delete elements from end to begin
            parKeys.erase(parKeys.begin() + remIndices[remIndices.size() -
              indI - 1]);
            comboKeys.erase(comboKeys.begin() + remIndices[remIndices.size() -
              indI - 1]);
        }
    }
}

