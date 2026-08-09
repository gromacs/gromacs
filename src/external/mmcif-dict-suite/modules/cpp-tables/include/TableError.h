//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$

/*!
** \file TableError.h
**
** \brief Utility file that is to be removed soon. 
*/

/* 
  PURPOSE:    Error codes
*/

#ifndef __TABLE_ERROR_H__
#define __TABLE_ERROR_H__

const int NO_TABLE_ERROR = 0;

const int ROW_OUT_OF_BOUNDS = -201;
const int COLUMN_OUT_OF_BOUNDS = -202;
const int NO_TREE_ON_COLUMN = -203;
const int INDEX_NAME_NOT_FOUND = -204;
const int NEW_COLUMN_LENGTH_ZERO = -205;
const int ADD_UPDATE_NULL = -206;
const int COLUMN_NAME_NOT_FOUND = -207;
const int SOME_COLUMN_NAMES_NOT_FOUND = -208;
const int REGEX_COMPILE_FAILED = -209;
const int NO_APPROPRIATE_INDEX = -210;
const int NOT_FOUND = -211;
const int DELETED_ROW = -212;
const int INDEX_CORRUPTED = -213;
const int KEY_ERROR = -214;
const int TABLE_NOT_FOUND = -215;

const int ASSERT_WARNING = -275;

const int ASSERT_NULL_DATA_POINTER = -280;

const int DUPLICATE_ROW = -290;

const int TABLE_WARNING = -350; // anything smaller maybe a big error

const int NULL_COMPARISON = -400;
const int DOUBLE_CONVERSION_ERROR = -401;
const int INTEGER_CONVERSION_ERROR = -402;
const int NULL_SEARCH_LIST = -403;
const int NOT_A_DATATYPE_ERROR = -404;
const int ERROR_NO_FILE_NAVIGATOR = -405;
const int INTERNAL_INCONSISTENCY_ERROR = -406;

#endif
