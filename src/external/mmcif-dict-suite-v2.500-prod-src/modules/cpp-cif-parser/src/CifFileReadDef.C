//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$

/*!
** \file CifFileReadDef.C
**
** \brief Implementation file for CifFileReadDef class.
*/

/* 
  PURPOSE:    Definitions for selective parsing/reading cif file
*/


#include <vector>
#include "GenString.h"
#include "CifFileReadDef.h"


using std::string;
using std::vector;


CifFileReadDef::CifFileReadDef(vector<string> dblist,vector<string> clist, type ctype, type dbtype) {
  
  _datablocklist = dblist;
  _categorylist = clist;
  _datablocklisttype = ctype;
  _categorylisttype = dbtype;
  _numCatsToRead = INVALID_NUM_CATS;
  _numReadCats = 0;
}


void CifFileReadDef::SetDataBlockList(vector<string> dblist,
  type dbtype)
{

  _datablocklist=dblist;
  _datablocklisttype=dbtype;

  SetNumCatsToRead();

}


void CifFileReadDef::SetCategoryList(vector<string>clist,
  type ctype)
{

  _categorylist=clist;
  _categorylisttype=ctype;

  SetNumCatsToRead();

}


void CifFileReadDef::IncreaseNumReadCats()
{

  /*
  ** If number of categories to read is invalid, it makes no
  ** sense to do anything, since the counter becomes invalid.
  */
  if (_numCatsToRead == INVALID_NUM_CATS)
  {
    return;
  }

  _numReadCats++;

}


int CifFileReadDef::AreAllCatsRead()
{

  /*
  ** If number of categories to read is invalid, return false to indicate
  ** that more reading is to be done.
  */
  if (_numCatsToRead == INVALID_NUM_CATS)
  {
    return(false);
  }

  if (_numReadCats == _numCatsToRead)
  {
    return(true);
  }
  else
  {
    return(false);
  }

}


int CifFileReadDef::Category_OK(const string& categoryName){
/*
  Returns 0 - category not OK, which means that named caterogy should be 
  avoided during reading proces. In that case category 
  - is in _categorylist and _categorylisttype is set to D (denied)
  - or is not in _categorylist and _categorylisttype is set to A (accepted)

  Returns 1 - category is OK, which means read that category. In that case
  category
  - is in _categorylist and _categorylisttype is set to A (accepted)
  - or is not in _categorylist and _categorylisttype is set to D (denied)
*/
  int i;
  int len = _categorylist.size();

  if (len == 0) return 1;

    i = 0;
    while(i<len && !String::IsCiEqual(categoryName,_categorylist[i])){
      i++;
    }
    if (((i==len)&&(_categorylisttype==A))||((i!=len)&&(_categorylisttype==D)))
      return (0);
    else
      return (1);
}


int CifFileReadDef::Datablock_OK(const string& datablockName){
/*
  same as Category_OK just chacks _datablocklist and _datablocklisttyepe
*/
  int i;
  int len = _datablocklist.size();

  if (len == 0) return 1;

    i = 0;
    while(i<len && !String::IsCiEqual(datablockName,_datablocklist[i])){
      i++;
    }
    if (((i==len)&&(_datablocklisttype==A))||((i!=len)&&(_datablocklisttype==D)))
      return (0);
    else
      return (1);
}


void CifFileReadDef::SetNumCatsToRead()
{

  if ((_datablocklisttype == D) || (_categorylisttype == D))
  {
    /*
    ** If the read type for blocks or for categories is set to "denied",
    ** the total number of categories can not be calculated. This is because
    ** it is not possible to know, without parsing the whole file, how
    ** many categories will be read. Hence, in this case, parsing of the
    ** whole file is inevitable.
    */
    _numCatsToRead = INVALID_NUM_CATS;
    return;
  }

  /*
  ** When both read types for blocks and categories is set to "accepted",
  ** the valid number is set only if non-ALL is selected for both
  ** categories and data blocks. Non-ALL is indicated by non-zero array size.
  */
  if ((_categorylist.size() == 0) || (_datablocklist.size() == 0))
  {
      _numCatsToRead = INVALID_NUM_CATS;
      return;
  }

  _numCatsToRead = _categorylist.size() * _datablocklist.size();

}
