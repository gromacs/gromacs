//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef CIFSTRING_H
#define CIFSTRING_H


#include <string>
#include <vector>


/**
 ** \class CifString
 ** 
 ** \brief Public class that contains CIF string related static methods.
 ** 
 ** This class is not a full abstraction of a CIF string. It only contains
 ** static constants and methods, that are related to a CIF string. A CIF
 ** string is a string, prefixed with an underscore, that consists of a
 ** category name and an item name concatenated by a dot, as specified here:
 **
 ** _categoryName.itemName
 **
 ** The class provides methods for creating a CIF string, extracting category
 ** name and item name from a CIF string.
 */
class CifString
{
public:
    static const char PREFIX_CHAR = '_';
    static const char JOIN_CHAR = '.';

    static const char NULL_CHAR = '?';
    static const char NOT_APPROPRIATE_CHAR = '.';

    static const std::string CIF_DDL_CATEGORY_BLOCK;
    static const std::string CIF_DDL_CATEGORY_DATABLOCK;
    static const std::string CIF_DDL_CATEGORY_DATABLOCK_METHODS;
    static const std::string CIF_DDL_CATEGORY_ITEM;
    static const std::string CIF_DDL_CATEGORY_ITEM_LINKED;
    static const std::string CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP;
    static const std::string CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP_LIST;
    static const std::string CIF_DDL_CATEGORY_CATEGORY;
    static const std::string CIF_DDL_CATEGORY_CATEGORY_EXAMPLES;
    static const std::string CIF_DDL_CATEGORY_NDB_CATEGORY_EXAMPLES;
    static const std::string CIF_DDL_CATEGORY_CATEGORY_KEY;
    static const std::string CIF_DDL_CATEGORY_CATEGORY_GROUP;
    static const std::string CIF_DDL_CATEGORY_CATEGORY_GROUP_LIST;
    static const std::string CIF_DDL_CATEGORY_CATEGORY_METHODS;
    static const std::string CIF_DDL_CATEGORY_SUB_CATEGORY;
    static const std::string CIF_DDL_CATEGORY_SUB_CATEGORY_EXAMPLES;
    static const std::string CIF_DDL_CATEGORY_SUB_CATEGORY_METHODS;
    static const std::string CIF_DDL_CATEGORY_ITEM_SUB_CATEGORY;
    static const std::string CIF_DDL_CATEGORY_ITEM_TYPE;
    static const std::string CIF_DDL_CATEGORY_ITEM_TYPE_CONDITIONS;
    static const std::string CIF_DDL_CATEGORY_ITEM_METHODS;
    static const std::string CIF_DDL_CATEGORY_ITEM_TYPE_LIST;
    static const std::string CIF_DDL_CATEGORY_ITEM_STRUCTURE;
    static const std::string CIF_DDL_CATEGORY_ITEM_STRUCTURE_LIST;
    static const std::string CIF_DDL_CATEGORY_ITEM_DESCRIPTION;
    static const std::string CIF_DDL_CATEGORY_NDB_ITEM_DESCRIPTION;
    static const std::string CIF_DDL_CATEGORY_NDB_CATEGORY_DESCRIPTION;
    static const std::string CIF_DDL_CATEGORY_ITEM_EXAMPLES;
    static const std::string CIF_DDL_CATEGORY_NDB_ITEM_EXAMPLES;
    static const std::string CIF_DDL_CATEGORY_ITEM_DEPENDENT;
    static const std::string CIF_DDL_CATEGORY_ITEM_RELATED;
    static const std::string CIF_DDL_CATEGORY_ITEM_RANGE;
    static const std::string CIF_DDL_CATEGORY_ITEM_ENUMERATION;
    static const std::string CIF_DDL_CATEGORY_NDB_ITEM_ENUMERATION;
    static const std::string CIF_DDL_CATEGORY_ITEM_DEFAULT;
    static const std::string CIF_DDL_CATEGORY_ITEM_ALIASES;
    static const std::string CIF_DDL_CATEGORY_DICTIONARY;
    static const std::string CIF_DDL_CATEGORY_DICTIONARY_HISTORY;
    static const std::string CIF_DDL_CATEGORY_ITEM_UNITS;
    static const std::string CIF_DDL_CATEGORY_ITEM_UNITS_LIST;
    static const std::string CIF_DDL_CATEGORY_ITEM_UNITS_CONVERSION;
    static const std::string CIF_DDL_CATEGORY_METHOD_LIST;

    static const std::string CIF_DDL_ITEM_ID;
    static const std::string CIF_DDL_ITEM_CATEGORY_ID;
    static const std::string CIF_DDL_ITEM_SUB_CATEGORY_ID;
    static const std::string CIF_DDL_ITEM_METHOD_ID;
    static const std::string CIF_DDL_ITEM_PARENT_NAME;
    static const std::string CIF_DDL_ITEM_CHILD_NAME;
    static const std::string CIF_DDL_ITEM_CHILD_CATEGORY_ID;
    static const std::string CIF_DDL_ITEM_PARENT_CATEGORY_ID;
    static const std::string CIF_DDL_ITEM_LINK_GROUP_ID;
    static const std::string CIF_DDL_ITEM_LABEL;
    static const std::string CIF_DDL_ITEM_CONTEXT;
    static const std::string CIF_DDL_ITEM_CONDITION_ID;
    static const std::string CIF_DDL_ITEM_ALIAS_NAME;
    static const std::string CIF_DDL_ITEM_DICTIONARY;
    static const std::string CIF_DDL_ITEM_TITLE;
    static const std::string CIF_DDL_ITEM_VERSION;
    static const std::string CIF_DDL_ITEM_NAME;
    static const std::string CIF_DDL_ITEM_CODE;
    static const std::string CIF_DDL_ITEM_PRIMITIVE_CODE;
    static const std::string CIF_DDL_ITEM_CONSTRUCT;
    static const std::string CIF_DDL_ITEM_ORGANIZATION;
    static const std::string CIF_DDL_ITEM_INDEX;
    static const std::string CIF_DDL_ITEM_DIMENSION;
    static const std::string CIF_DDL_ITEM_DATABLOCK_ID;
    static const std::string CIF_DDL_ITEM_DESCRIPTION;
    static const std::string CIF_DDL_ITEM_NDB_DESCRIPTION;
    static const std::string CIF_DDL_ITEM_CASE;
    static const std::string CIF_DDL_ITEM_MANDATORY_CODE;
    static const std::string CIF_DDL_ITEM_DETAIL;
    static const std::string CIF_DDL_ITEM_MAXIMUM;
    static const std::string CIF_DDL_ITEM_MINIMUM;
    static const std::string CIF_DDL_ITEM_VALUE;
    static const std::string CIF_DDL_ITEM_DEPENDENT_NAME;
    static const std::string CIF_DDL_ITEM_RELATED_NAME;
    static const std::string CIF_DDL_ITEM_FUNCTION_CODE;
    static const std::string CIF_DDL_ITEM_OFFSET;
    static const std::string CIF_DDL_ITEM_OPERATOR;
    static const std::string CIF_DDL_ITEM_FACTOR;
    static const std::string CIF_DDL_ITEM_FROM_CODE;
    static const std::string CIF_DDL_ITEM_TO_CODE;
    static const std::string CIF_DDL_ITEM_UPDATE;
    static const std::string CIF_DDL_ITEM_REVISION;
    static const std::string CIF_DDL_ITEM_INLINE;
    static const std::string CIF_DDL_ITEM_LANGUAGE;
    static const std::string CIF_DDL_ITEM_PARENT_ID;

    static const std::string UnknownValue;
    static const std::string InapplicableValue;

    static void MakeCifItem(std::string& cifItem,
      const std::string& categoryName, const std::string& itemName);
    static void MakeCifItems(std::vector<std::string>& cifItems,
      const std::string& categoryName,
      const std::vector<std::string>& attribsNames);

    static void GetItemFromCifItem(std::string& keyword,
      const std::string& itemName);
    static void GetCategoryFromCifItem(std::string& categoryName,
      const std::string& itemName);

    static bool IsEmptyValue(const std::string& value);
    static bool IsUnknownValue(const std::string& value);

    static bool IsSpecialChar(const char charValue);
    static bool IsSpecialFirstChar(const char charValue);
};

#endif
