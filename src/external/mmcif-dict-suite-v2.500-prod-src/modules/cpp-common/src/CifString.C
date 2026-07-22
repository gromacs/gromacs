//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#include <string>
#include <vector>

#include "Exceptions.h"
#include "CifString.h"


using std::string;
using std::vector;


const string CifString::CIF_DDL_CATEGORY_BLOCK("datablock");
const string CifString::CIF_DDL_CATEGORY_DATABLOCK("datablock");
const string CifString::CIF_DDL_CATEGORY_DATABLOCK_METHODS("datablock_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM("item");
const string CifString::CIF_DDL_CATEGORY_ITEM_LINKED("item_linked");
const string CifString::CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP("pdbx_item_linked_group");
const string CifString::CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP_LIST("pdbx_item_linked_group_list");
const string CifString::CIF_DDL_CATEGORY_CATEGORY("category");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_EXAMPLES("category_examples");
const string CifString::CIF_DDL_CATEGORY_NDB_CATEGORY_EXAMPLES("ndb_category_examples");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_KEY("category_key");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_GROUP("category_group");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_GROUP_LIST("category_group_list");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_METHODS("category_methods");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY("sub_category");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY_EXAMPLES("sub_category_examples");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY_METHODS("sub_category_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM_SUB_CATEGORY("item_sub_category");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE("item_type");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE_CONDITIONS("item_type_conditions");
const string CifString::CIF_DDL_CATEGORY_ITEM_METHODS("item_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE_LIST("item_type_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_STRUCTURE("item_structure");
const string CifString::CIF_DDL_CATEGORY_ITEM_STRUCTURE_LIST("item_structure_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_DESCRIPTION("item_description");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_DESCRIPTION("ndb_item_description");
const string CifString::CIF_DDL_CATEGORY_NDB_CATEGORY_DESCRIPTION("ndb_category_description");
const string CifString::CIF_DDL_CATEGORY_ITEM_EXAMPLES("item_examples");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_EXAMPLES("ndb_item_examples");
const string CifString::CIF_DDL_CATEGORY_ITEM_DEPENDENT("item_dependent");
const string CifString::CIF_DDL_CATEGORY_ITEM_RELATED("item_related");
const string CifString::CIF_DDL_CATEGORY_ITEM_RANGE("item_range");
const string CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION("item_enumeration");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_ENUMERATION("ndb_item_enumeration");
const string CifString::CIF_DDL_CATEGORY_ITEM_DEFAULT("item_default");
const string CifString::CIF_DDL_CATEGORY_ITEM_ALIASES("item_aliases");
const string CifString::CIF_DDL_CATEGORY_DICTIONARY("dictionary");
const string CifString::CIF_DDL_CATEGORY_DICTIONARY_HISTORY("dictionary_history");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS("item_units");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS_LIST("item_units_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS_CONVERSION("item_units_conversion");
const string CifString::CIF_DDL_CATEGORY_METHOD_LIST("method_list");

const string CifString::CIF_DDL_ITEM_ID("id");
const string CifString::CIF_DDL_ITEM_CATEGORY_ID("category_id");
const string CifString::CIF_DDL_ITEM_SUB_CATEGORY_ID("sub_category_id");
const string CifString::CIF_DDL_ITEM_METHOD_ID("method_id");
const string CifString::CIF_DDL_ITEM_PARENT_NAME("parent_name");
const string CifString::CIF_DDL_ITEM_CHILD_NAME("child_name");
const string CifString::CIF_DDL_ITEM_CHILD_CATEGORY_ID("child_category_id");
const string CifString::CIF_DDL_ITEM_PARENT_CATEGORY_ID("parent_category_id");
const string CifString::CIF_DDL_ITEM_LINK_GROUP_ID("link_group_id");
const string CifString::CIF_DDL_ITEM_LABEL("label");
const string CifString::CIF_DDL_ITEM_CONTEXT("context");
const string CifString::CIF_DDL_ITEM_CONDITION_ID("condition_id");
const string CifString::CIF_DDL_ITEM_ALIAS_NAME("alias_name");
const string CifString::CIF_DDL_ITEM_DICTIONARY("dictionary");
const string CifString::CIF_DDL_ITEM_TITLE("title");
const string CifString::CIF_DDL_ITEM_VERSION("version");
const string CifString::CIF_DDL_ITEM_NAME("name");
const string CifString::CIF_DDL_ITEM_CODE("code");
const string CifString::CIF_DDL_ITEM_PRIMITIVE_CODE("primitive_code");
const string CifString::CIF_DDL_ITEM_CONSTRUCT("construct");
const string CifString::CIF_DDL_ITEM_ORGANIZATION("organization");
const string CifString::CIF_DDL_ITEM_INDEX("index");
const string CifString::CIF_DDL_ITEM_DIMENSION("dimension");
const string CifString::CIF_DDL_ITEM_DATABLOCK_ID("datablock_id");
const string CifString::CIF_DDL_ITEM_DESCRIPTION("description");
const string CifString::CIF_DDL_ITEM_NDB_DESCRIPTION("ndb_description");
const string CifString::CIF_DDL_ITEM_CASE("case");
const string CifString::CIF_DDL_ITEM_MANDATORY_CODE("mandatory_code");
const string CifString::CIF_DDL_ITEM_DETAIL("detail");
const string CifString::CIF_DDL_ITEM_MAXIMUM("maximum");
const string CifString::CIF_DDL_ITEM_MINIMUM("minimum");
const string CifString::CIF_DDL_ITEM_VALUE("value");
const string CifString::CIF_DDL_ITEM_DEPENDENT_NAME("dependent_name");
const string CifString::CIF_DDL_ITEM_RELATED_NAME("related_name");
const string CifString::CIF_DDL_ITEM_FUNCTION_CODE("function_code");
const string CifString::CIF_DDL_ITEM_OFFSET("file_offset");
const string CifString::CIF_DDL_ITEM_OPERATOR("operator");
const string CifString::CIF_DDL_ITEM_FACTOR("factor");
const string CifString::CIF_DDL_ITEM_FROM_CODE("from_code");
const string CifString::CIF_DDL_ITEM_TO_CODE("to_code");
const string CifString::CIF_DDL_ITEM_UPDATE("update");
const string CifString::CIF_DDL_ITEM_REVISION("revision");
const string CifString::CIF_DDL_ITEM_INLINE("item_inline");
const string CifString::CIF_DDL_ITEM_LANGUAGE("item_language");
const string CifString::CIF_DDL_ITEM_PARENT_ID("parent_id");


const string CifString::UnknownValue("?");
const string CifString::InapplicableValue(".");


void CifString::MakeCifItem(string& cifItem, const string& categoryName,
  const string& attribName)
{
    if (categoryName.empty())
    {
        throw EmptyValueException("Empty category name",
          "CifString::MakeCifItem");
    }

    if (attribName.empty())
    {
        throw EmptyValueException("Empty attribute name",
          "CifString::MakeCifItem");
    }

    cifItem = PREFIX_CHAR + categoryName + JOIN_CHAR + attribName;
}


void CifString::MakeCifItems(vector<string>& cifItems,
  const string& categoryName, const vector<string>& attribsNames)
{
    cifItems.assign(attribsNames.size(), string());

    for (unsigned int attribI = 0; attribI < attribsNames.size(); ++attribI)
    {
        MakeCifItem(cifItems[attribI], categoryName, attribsNames[attribI]);
    }
}


void CifString::GetItemFromCifItem(string& itemName, const string& cifItem)
{
    /* ----------------------------------------------------------------------
     Purpose: CifString::GetItemFromCifItem(itemName, cifItem)

     Get the itemName part of an item name (ie. _<category>.<itemName>)
     Return 1 for success or 0 otherwise.

     *--------------------------------------------------------------------- */

    unsigned int k, ilen;

    if (cifItem.empty() || (cifItem[0] != '_') || (ilen = cifItem.size()) < 4)
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetItemFromCifItem");
    }

    for (k = 1; k < ilen; ++k)
    {
        if (cifItem[k] == JOIN_CHAR)
            break;
    }

    if ((k == 1) || (k == (ilen - 1)))
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetItemFromCifItem");
    }

    itemName.clear();

    for (unsigned int i = k + 1; i < ilen; ++i)
    {
        itemName.push_back(cifItem[i]);
    }
}


void CifString::GetCategoryFromCifItem(string& categoryName,
    const string& cifItem)
/* ----------------------------------------------------------------------
 Purpose: CifString::GetCategoryFromCifItem(categoryName, cifItem)

 Get the category part of an item name (ie. _<category>.<keyword>)
 Return 1 for success or 0 otherwise.

 * ---------------------------------------------------------------------- */
{
    if (cifItem.empty() || (cifItem[0] != '_'))
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetCategoryFromCifItem");
    }

    // Skip the first char and search for a join character
    string::size_type dotIndex = cifItem.find(JOIN_CHAR, 1);

    if (dotIndex == string::npos)
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetCategoryFromCifItem");
#ifdef VLAD_DEL
        // If join char not found, set its index to be the end of the string.
        dotIndex = cifItem.size();
#endif
    }

    categoryName.clear();

    // Copy from the first char to the character prior to the dot char.
    categoryName.assign(cifItem, 1, dotIndex - 1);
}


bool CifString::IsEmptyValue(const string& value)
{
    if (value.empty() || (value == InapplicableValue) ||
      (value == UnknownValue))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool CifString::IsUnknownValue(const string& value)
{
    if (value.empty() || (value == UnknownValue))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool CifString::IsSpecialChar(const char charValue)
{
    switch (charValue)
    {
        case '(':
        case ')':
        case '[':
        case ']':
        case '{':
        case '}':
        {
            return (true);
            break;
        }
        default:
        {
            return (false);
            break;
        }
    }
}


bool CifString::IsSpecialFirstChar(const char charValue)
{
    switch (charValue)
    {
        case '$':
        case '#':
        case '_':
        case ';':
        case '(':
        case ')':
        case '[':
        case ']':
        case '{':
        case '}':
        {
            return (true);
            break;
        }
        default:
        {
            return (false);
            break;
        }
    }
}

