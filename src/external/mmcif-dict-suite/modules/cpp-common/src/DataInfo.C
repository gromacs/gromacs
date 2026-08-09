//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include <iostream>

#include "Exceptions.h"
#include "rcsb_types.h"
#include "GenCont.h"
#include "RcsbFile.h"
#include "CifString.h"
#include "DataInfo.h"


using std::string;
using std::vector;
using std::cerr;
using std::endl;


#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
string CIF_ITEM;
#endif


DataInfo::DataInfo()
{

}


DataInfo::~DataInfo()
{

}


bool DataInfo::IsUnknownValueAllowed(const string& catName, 
  const string& attribName)
{
    return (true);
}


bool DataInfo::AreItemsValuesValid(const string& catName,
  const vector<string>& columnNames, const vector<unsigned int>& columnIndices,
  const vector<bool>& allowedNullAttribs, const vector<string>& row,
  const Char::eCompareType compareType)
{
    if (row.empty())
    {
        cerr << "  Empty row detected." << endl;
        return (false);
    }

    for (unsigned int i = 0; i < columnNames.size(); ++i)
    {
        if (IsKeyItem(catName, columnNames[i], compareType))
        {
            if (CifString::IsEmptyValue(row[columnIndices[i]]))
            {
                string cifItem;
                CifString::MakeCifItem(cifItem, catName, columnNames[i]);

#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
                if ((cifItem == "_atom_sites_alt.id") &&
                  (row[columnIndices[i]] == CifString::InapplicableValue))
                {
                    CIF_ITEM = cifItem;
                    return (false);
                }
#endif
                cerr << "  Key item \"" << cifItem << "\" has invalid value \""
                  << row[columnIndices[i]] << "\"" << endl;

                return (false);
            }
        }
        else
        {
            if (CifString::IsUnknownValue(row[columnIndices[i]]))
            {
                if (!allowedNullAttribs[i])
                {
                    string cifItem;
                    CifString::MakeCifItem(cifItem, catName, columnNames[i]);

                    cerr << "  Non-key item \"" << cifItem << "\", that must "\
                      "not have unknown value, has invalid value \"" <<
                      row[columnIndices[i]] << "\"" << endl;

                    return (false);
                }
            }
        }
    }

    return (true);
}


bool DataInfo::IsKeyItem(const string& catName,
  const string& attribName, const Char::eCompareType compareType)
{
    if (attribName.empty())
        return(false);

    const vector<string>& keys = GetCatKeys(catName);

    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    return (GenCont::IsInVector(itemName, keys));
}


bool DataInfo::AreAllKeyItems(const string& catName,
  const vector<string>& attributes)
{
    // Return true for cases where all attributes are key.

    if (attributes.empty())
        return(false);

    const vector<string>& keys = GetCatKeys(catName);
 
    unsigned int keyCount = 0;

    for (unsigned int i = 0; i < attributes.size(); ++i)
    {
        if (IsKeyItem(catName, attributes[i]))
        {
            keyCount++;
            continue;
        }
        else
        {
            keyCount = 0;
            break;
        }
    }

    if (!keys.empty() && (keyCount == keys.size()))
    {
        return (true);
    }
    else
    {
        if (keyCount != 0)
        {
            cerr << "Missing keys in category: " << catName << endl;
        }
        return (false);
    }
}


bool DataInfo::MustConvertItem(const string& catName,
  const string& itemName)
{
    string cifItem;
    CifString::MakeCifItem(cifItem, catName, itemName);

    return (IsItemDefined(cifItem));
}


void DataInfo::GetItemsTypes(vector<eTypeCode>& itemsTypes,
  const string& catName, const vector<string>& attribsNames)
{
    itemsTypes.clear();

    for (unsigned int i = 0; i < attribsNames.size(); ++i)
    {
        string cifItem;
        CifString::MakeCifItem(cifItem, catName, attribsNames[i]);

        if (!MustConvertItem(catName, attribsNames[i]))
        {
            itemsTypes.push_back(eTYPE_CODE_NONE);

            continue;
        }

        itemsTypes.push_back(_GetDataType(cifItem));
    }
}


eTypeCode DataInfo::_GetDataType(const string& itemName) 
{
    const vector<string>& primitiveCode = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_TYPE_LIST,
      CifString::CIF_DDL_ITEM_PRIMITIVE_CODE);

    const vector<string>& dataType = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_TYPE,
      CifString::CIF_DDL_ITEM_CODE);

    if (primitiveCode.empty() || dataType.empty())
    {
        cerr << "Warning - ITEM_TYPE_CODE: Item \"" + itemName +
          "\" has no type code defined. Using default string type." << endl;
#ifdef VLAD_TEST
        throw EmptyValueException("Item \"" + itemName + "\" has no type "\
          "code defined.", "DataInfo::_GetDataType");
#endif

        // Missing item primitive or type code, assume string type
        return (eTYPE_CODE_STRING);
    }

    if (primitiveCode[0].empty() || dataType[0].empty())
    {
        cerr << "Warning - ITEM_TYPE_CODE: Item \"" + itemName +
          "\" has no type code defined. Using default string type." << endl;
#ifdef VLAD_TEST
        throw EmptyValueException("Item \"" + itemName + "\" has no type "\
          "code defined.", "DataInfo::_GetDataType");
#endif
        // Missing item primitive or type code, assume string type
        return (eTYPE_CODE_STRING);
    }

    if (String::IsCiEqual(primitiveCode[0], "char") ||
      String::IsCiEqual(primitiveCode[0], "uchar"))
    {
        if (String::IsCiEqual(dataType[0], "yyyy-mm-dd"))
        {
            return (eTYPE_CODE_DATETIME);
        }
        else if (String::IsCiEqual(dataType[0], "text"))
        {
            return (eTYPE_CODE_TEXT);
        }
        else
        {
            // Assume string type
            return (eTYPE_CODE_STRING);
        }
    }
    else if (String::IsCiEqual(primitiveCode[0], "numb"))
    {
        if (String::IsCiEqual(dataType[0], "float"))
        {
            return (eTYPE_CODE_FLOAT);
        }
        else if (String::IsCiEqual(dataType[0], "double"))
        {
            return (eTYPE_CODE_FLOAT);
        }
        else if (String::IsCiEqual(dataType[0], "int"))
        {
            return (eTYPE_CODE_INT);
        }
        else if (String::IsCiEqual(dataType[0], "numb"))
        {
            return (eTYPE_CODE_INT);
        }
        else
        {
            // Assume string type
            return (eTYPE_CODE_STRING);
        }
    }
    else
    {
        throw EmptyValueException("Item \"" + itemName + "\" has invalid "\
          "primitive code \"" + primitiveCode[0] + "\"",
          "DataInfo::_GetDataType");
    }
}


bool DataInfo::IsItemMandatory(const string& catName,
  const string& attribName)
{
    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    return (IsItemMandatory(itemName));
}


void DataInfo::GetMandatoryItems(vector<string>& mandItemsNames,
  const string& catName)
{
    mandItemsNames.clear();

    const vector<string>& itemsNames = GetItemsNames();
    for (unsigned int itemI = 0; itemI < itemsNames.size(); ++itemI)
    {
        const string& itemName = itemsNames[itemI];

        string currCatName;
        CifString::GetCategoryFromCifItem(currCatName, itemName);
        if (currCatName != catName)
        {
            continue;
        }

        if (IsItemMandatory(itemName))
        {
            mandItemsNames.push_back(itemName);
        }
    }
}


bool DataInfo::IsItemMandatory(const string& itemName)
{
    const vector<string>& mCode = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM,
      CifString::CIF_DDL_ITEM_MANDATORY_CODE);

    if (String::IsCiEqual(mCode[0], "Y") || String::IsCiEqual(mCode[0], "YES"))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool DataInfo::IsSimpleDataType(const string& itemName)
{
    bool simpleTyping = false;

    bool iRange = false;
    bool iUnits = false;

    const vector<string>& enums = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
      CifString::CIF_DDL_ITEM_VALUE);

    const vector<string>& units = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_UNITS,
      CifString::CIF_DDL_ITEM_CODE);

    if (!units.empty())
    {
      iUnits = true;
    }

    const vector<string>& rangeMin = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_RANGE,
      CifString::CIF_DDL_ITEM_MINIMUM);

    const vector<string>& rangeMax = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_RANGE,
      CifString::CIF_DDL_ITEM_MAXIMUM);

    if ((!rangeMin.empty()) && (!rangeMax.empty()) &&
      (rangeMin.size() == rangeMax.size()))
    {
        iRange = true;
    }

    if (!iRange && enums.empty() && !iUnits)
    {
        simpleTyping = true;
    }

    return (simpleTyping);
}


void DataInfo::StandardizeEnumItem(string& value,
  const string& catName, const string& attribName)
{
    // If value is empty, no need to go and search. Return.
    if (value.empty())
        return;

    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    const vector<string>& enums = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
      CifString::CIF_DDL_ITEM_VALUE);

    for (unsigned int i = 0; i < enums.size(); ++i)
    {
        if (String::IsCiEqual(value, enums[i]))
        {
            value = enums[i];
            break;
        }
    }
}

