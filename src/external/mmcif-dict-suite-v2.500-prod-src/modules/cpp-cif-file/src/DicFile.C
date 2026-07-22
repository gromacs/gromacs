//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DicFile.C
**
** \brief Implementation file for DicFile class.
*/

#include <algorithm>

#include "GenString.h"
#include "CifString.h"
#ifdef VLAD_DO_WE_NEED_THIS
#include "CifDataInfo.h"
#include "CifParentChild.h"
#endif
#include "DicFile.h"


using std::ofstream;
using std::ostream;
using std::ios;
using std::endl;


DicFile::DicFile(const eFileMode fileMode, const string& objFileName,
  const bool verbose, const Char::eCompareType caseSense,
  const unsigned int maxLineLength, const string& nullValue) :
  CifFile(fileMode, objFileName, verbose, caseSense, maxLineLength, nullValue),
  _formatP(NULL)
{

    _formatP = new ISTable("ddlformat");

    _formatP->AddColumn("dbName");
    _formatP->AddColumn("type");
    _formatP->AddColumn("catName");

}


DicFile::DicFile(const bool verbose,
  const Char::eCompareType caseSense, const unsigned int maxLineLength,
  const string& nullValue) : CifFile(verbose, caseSense, maxLineLength,
  nullValue), _formatP(NULL)
{

    _formatP = new ISTable("ddlformat");

    _formatP->AddColumn("dbName");
    _formatP->AddColumn("type");
    _formatP->AddColumn("catName");

}


DicFile::~DicFile()
{

    if (_formatP != NULL)
    {
        delete(_formatP);
    }

}


ISTable* DicFile::GetFormatTable()
{

    return(_formatP);

}


int DicFile::WriteFormatted(const string& cifFileName, ISTable* formatP)
{

    int iret = 0;

    ofstream cifo(cifFileName.c_str(), ios::out | ios::trunc);

    if (formatP != NULL)
        iret = WriteFormatted(cifo, formatP);
    else
        iret = WriteFormatted(cifo, _formatP);

    cifo.close();

    return(iret);

}


int DicFile::WriteFormatted(ostream& cifo, ISTable* formatP)
{

    return(WriteFormatted(cifo, this, formatP));

}


int DicFile::WriteFormatted(const string& cifFileName, TableFile* ddl,
  ISTable* formatP)
{

    int iret = 0;

    ofstream cifo(cifFileName.c_str(), ios::out | ios::trunc);

    if (formatP != NULL)
        iret = WriteFormatted(cifo, ddl, formatP);
    else
        iret = WriteFormatted(cifo, ddl, _formatP);

    cifo.close();

    return(iret);

}


int DicFile::WriteFormatted(ostream& cifo, TableFile* ddl, ISTable* formatP)
{

  ISTable *tblP = NULL;
  ISTable *cattbl = NULL;
  ISTable *cattbl2 = NULL;
  ISTable *itemtbl = NULL;
  ISTable *itemtblddl = NULL;
  
  
  unsigned int numColumn;

  int k, ilen;
  unsigned int linePos;
  int i2,i3,i4;
  int cwid, *cwidth=NULL;
  string categoryName;
  string categoryName2;
  string catName;
  string itemName;
  string itemName2;
  string item2;
  vector<string> TableNames;
  string BlockName;
 
  vector<string> list;
  vector<string> list2;
  vector<string> listcat;
  vector<string> listcat2;
  vector<string> listitem;
  vector<string> listitem2;

  vector<unsigned int> listOut;
  vector<unsigned int> listOutcat;
  vector<unsigned int> listOutcat2;
  vector<unsigned int> listOutItem;
  vector<unsigned int> listOutItem2;

  vector<string> target;
  vector<string> target2;
  vector<string> targetcat;
  vector<string> targetcat2;
  vector<string> targetitem;

  listitem2.clear();
  listitem2.push_back("category_id");

  list.push_back("dbName");
  list.push_back("type");

  /* 
     This is loop for all datablocks - we have to assume that dictionary might have more then one datablock
  */
  
  for (unsigned int ib0 = 0; ib0 < _blocks.size(); ++ib0)
  {
    string loopBlockName = _blocks[ib0].GetName();

    Block& block = GetBlock(loopBlockName);

#if DEBUG    
    cerr <<  "Writing data block " << loopBlockName << "["<< ib0 << " of " 
	 << _blocks.size() << "]"<< " Null is " << _nullValue << endl;
#endif

    cifo << "data_" << loopBlockName << endl;

    /*
      Table format has information wich tables have information about datablock, category save frame and item save frame. This table is creatid during reading dictionary, and have to be input in this method.
    */

    target.push_back(loopBlockName);
    target.push_back("data");

    formatP->Search(listOut, target, list);

    if (!listOut.empty())
    {
      for (unsigned int l = 0; l < listOut.size(); ++l)
      {
        //format->GetCell(categoryName, string("catName"), listOut[l]);
        categoryName = (*formatP)(listOut[l], "catName");

        tblP = block.GetTablePtr(categoryName);

        unsigned int numColumn = tblP->GetNumColumns();
        unsigned int numRow    = tblP->GetNumRows();

        cifo << "# " << endl;

        if (numRow <= 0)
        {
	  if (true /*emptyTable*/)
          {
            const vector<string>& colNames = tblP->GetColumnNames();
	    cwid = -1;
	    for (unsigned int i=0; i< numColumn; i++)
            {
	      ilen = colNames[i].size();
	      if (ilen > cwid) cwid=ilen;
	    }
	
	    for (unsigned int i=0; i< numColumn; i++)
            {
              linePos = 0;
	      _PrintItemName(cifo, categoryName, colNames[i], linePos);
              _PrintPostItemSeparator(cifo, linePos);
	      ilen = colNames[i].size();
	      for (k=0; k < 2+cwid-ilen; k++)
                  cifo << " ";
	      linePos = cwid+2;
	  
	      cifo << _nullValue; 
	      if (linePos != 0) cifo << endl;
	    }
	  }
        }
        else if (numRow == 1)
        {
          const vector<string>& colNames = tblP->GetColumnNames();
          const vector<string>& rowValues = tblP->GetRow(0);
          cwid = -1;
          for (unsigned int i=0; i< numColumn; i++)
          {
            ilen = colNames[i].size();
            if (  ilen > cwid) cwid=ilen;
          }
      
          for (unsigned int i=0; i< numColumn; i++)
          {
            linePos = 0;
            _PrintItemName(cifo, categoryName, colNames[i], linePos);
            _PrintPostItemSeparator(cifo, linePos);
            ilen = colNames[i].size();
            for (k=0; k < 2+cwid-ilen; k++)
              cifo << " ";
            linePos = cwid+2;
	  
	    _PrintItemValue(cifo, rowValues[i], linePos);
	    if (linePos != 0)
              cifo << endl;
	  }
        }
        else
        {
	  cifo << "loop_" << endl;
	
          const vector<string>& colNames = tblP->GetColumnNames();
	  for (unsigned int i=0; i< numColumn; i++)
          {
            linePos = 0;
	    _PrintItemName(cifo, categoryName, colNames[i], linePos);
            _PrintPostItemSeparator(cifo, linePos);
	    cifo << endl;
	  }
	
	  // ---
	  cwidth  = new int[numColumn];
      
	  for (unsigned int i=0; i< numColumn; i++) { cwidth[i]=-1+2; }
	  unsigned int j = 0;
	  unsigned int m = 0;
	  while (j < numRow)
          {
            const vector<string>& rowValues = tblP->GetRow(m);
	    if (!rowValues.empty())
            {
	      for (unsigned int i=0; i< numColumn; i++)
              {
	        ilen = rowValues[i].size();
	        if ( _IsQuotableText(rowValues[i]))
                  ilen+=2;
	        if (  ilen > cwidth[i] - 2)
                  cwidth[i]=ilen+2;
	      }
	      j++;
	    }
	    m++;
	  }
	  j = 0;
	  m = 0;
	  while (j < numRow)
          {
	    linePos = 0;
            const vector<string>& rowValues = tblP->GetRow(m);
	    if (!rowValues.empty())
            {
	      for (unsigned int i=0; i< numColumn; i++)
              {
	        ilen=_PrintItemValue(cifo, rowValues[i], linePos); 
	        if (linePos != 0)
                {
		  if (linePos + cwidth[i]-2-ilen < _maxCifLineLength)
		    for (k=0; k < cwidth[i]-2-ilen; k++)
                      cifo << " ";
		  linePos +=    cwidth[i]-2-ilen;
	        }
	      }
	      if (linePos != 0) cifo << endl;
	      j++;
	    }
	    m++;
	  }
	  if (cwidth) delete[] cwidth;
        }
      }
    }

    // category save frames

    // target for format-seaching category for category save frame
    target.clear();
    target.push_back(loopBlockName);
    target.push_back("category");

    // target for format-seaching category for item save frame
    target2.clear();
    target2.push_back(loopBlockName);
    target2.push_back("item");

    // Creates index on table item, in ddl and dictionary, on column named name
    // Assumtion that ddl has table item, and that table hase column name

    itemtbl=block.GetTablePtr("item");
    BlockName = ddl->GetFirstBlockName();
    Block& firstBlock = ddl->GetBlock(BlockName);
    itemtblddl=firstBlock.GetTablePtr("item");
    listitem.clear();
    listitem.push_back("name");

    // Creates index on table category, in ddl and dictionary, on column 
    // named id and index on columns id+implicit_key
    // Assumtion that ddl has table category, and that table hase column id
    // and implicit_key
    cattbl=block.GetTablePtr("category");
    listcat.clear();
    listcat.push_back("id");
    listcat.clear();
    listcat.push_back("implicit_key");

    // creating indices for all tables on column with mandatory_code=implicit
    // assumption: table item has columns category_id and mandatory_code
    block.GetTableNames(TableNames);
    unsigned int num = TableNames.size();
    listitem.clear();
    listitem.push_back("category_id");
    listitem.push_back("mandatory_code");
    for (unsigned int l=0; l<num; l++)
    {
      targetitem.clear();
      targetitem.push_back(TableNames[l]);
      targetitem.push_back("implicit");
      listOutItem.clear();
      // in table item in ddl searchs rows with value <TableName> in column
      // category_id and value "implicit" in column "mandatory_code"
      itemtblddl->Search(listOutItem, targetitem, listitem);
      if (!listOutItem.empty())
      {
	//itemtblddl->GetCell(itemName,string("name"),listOutItem[0]);
        itemName = (*itemtblddl)(listOutItem[0], "name");
	cattbl2=block.GetTablePtr(TableNames[l]);
	CifString::GetItemFromCifItem(item2, itemName);
	listcat2.clear();
	listcat2.push_back(item2);
      }
    }
    TableNames.clear();

    // loop for all categories in current datablock; print save frame
    targetcat.clear();
    targetcat.push_back(loopBlockName);
    listOutcat.clear();
    //searchs all tables that belong to current datablock
    cattbl->Search(listOutcat, targetcat, listcat);
    if (!listOutcat.empty())
    {
      for (i2=0; i2<(int)(listOutcat.size());i2++)
      {
        //cattbl->GetCell(categoryName,string("id"),listOutcat[i2]);
        categoryName = (*cattbl)(listOutcat[i2], "id");
        String::UpperCase(categoryName, catName);
        cifo << "save_" <<catName << endl;
        listOut.clear();
        // searchs in table format all row with <datablock_name> in column 
        // "dbName" and "category" in column "type"
        // that will find category, category_examples....
        // in all this tables are information wich have to be placed in
        // save frame for category
        // the following loop is for all this tables
        formatP->Search(listOut, target, list);
        if (!listOut.empty())
        {
          for (i3=0; i3<(int)(listOut.size());i3++)
          {
	    //format->GetCell(categoryName2,string("catName"),listOut[i3]);
            categoryName2 = (*formatP)(listOut[i3], "catName");
	    cattbl2=block.GetTablePtr(categoryName2);
            // implicit value for "category" is name of datablock for column
	    // "imlicit_key" and name of save frame for column "id", for all
	    // other tables implicit value is name of save frame
	    // assumption: in table item in ddl in column named implicit
	    // is defined wich colum has implicit value
	    // for table category, implicit value "imlicit_key" is name of 
	    // datablock, "id"= name of save frame
	    if (!String::IsCiEqual(categoryName2,"category"))
            {
	      targetitem.clear();
	      targetitem.push_back(categoryName2);
	      targetitem.push_back("implicit");
	      listOutItem.clear();
              // in table item in ddl searchs rows with value <categoryName2> in column
              // category_id and value "implicit" in column "mandatory_code"
	      itemtblddl->Search(listOutItem, targetitem, listitem);
	      if (!listOutItem.empty())
              {
	        //itemtblddl->GetCell(itemName,string("name"),listOutItem[0]);
                itemName = (*itemtblddl)(listOutItem[0], "name");
	        CifString::GetItemFromCifItem(item2, itemName);
	      }
	    }
	    else
	      item2 = "id";
	    listcat2.clear();
	    listcat2.push_back(item2);
	
	    targetcat2.clear();
	    targetcat2.push_back(categoryName);
	    listOutcat2.clear();
	    cattbl2->Search(listOutcat2, targetcat2, listcat2);
	    // writes values for particular category like ATOM_SITE,SYMMETRY,
	    // for all tables that describe category like category, category_key,
	    // category_group
	    if (!listOutcat2.empty())
            {
	      if (listOutcat2.size() == 1)
              {
	        // pair
	        numColumn = cattbl2->GetNumColumns();
                const vector<string>& colNames = cattbl2->GetColumnNames();
                const vector<string>& rowValues =
                  cattbl2->GetRow(listOutcat2[0]);
	        cwid = -1;
	        for (unsigned int i=0; i< numColumn; i++)
                {
	          ilen = colNames[i].size();
	          if (  ilen > cwid)
                    cwid=ilen;
	        }
	    
	        for (unsigned int i=0; i< numColumn; i++)
                {
                  _PrintItemIdent(cifo, linePos);
	          _PrintItemName(cifo, categoryName2, colNames[i], linePos);
                  _PrintPostItemSeparator(cifo, linePos, true);
	          _PrintItemValue(cifo, rowValues[i], linePos);
	          if (linePos != 0) 
                    cifo << endl;
	        }
	      }
	      if (listOutcat2.size() > 1)
              {
	        // loop
	        cifo << "     loop_" << endl;
	        numColumn = cattbl2->GetNumColumns();
                const vector<string>& colNames = cattbl2->GetColumnNames();
	        for (unsigned int i=0; i< numColumn; i++)
                {
                  _PrintItemIdent(cifo, linePos);
	          _PrintItemName(cifo, categoryName2, colNames[i], linePos);
                  _PrintPostItemSeparator(cifo, linePos, true);
	          cifo << endl;
	        }
	
	        // ---
	        cwidth  = new int[numColumn];
	  
	        for (unsigned int i=0; i< numColumn; i++)
                {
                  cwidth[i]=-1+2;
                }
	        for (unsigned int j=0; j < listOutcat2.size(); ++j)
                {
                  const vector<string>& rowValues =
	            cattbl2->GetRow(listOutcat2[j]);
	          for (unsigned int i=0; i< numColumn; i++)
                  {
	              ilen = rowValues[i].size();
	              if ( _IsQuotableText(rowValues[i]))
                        ilen+=2;
	      
	              if (  ilen > cwidth[i])
                        cwidth[i]=ilen;
	          }
	        }
	        for (unsigned int j=0; j < listOutcat2.size(); ++j)
                {
	          linePos = 0;
                  const vector<string>& rowValues =
	            cattbl2->GetRow(listOutcat2[j]);
	          for (unsigned int i=0; i< numColumn; i++)
                  {
	            ilen=_PrintItemValue(cifo, rowValues[i], linePos,eNONE,cwidth[i]); 
	          }
	          if (linePos != 0) cifo << endl;
	        }
	        if (cwidth) delete[] cwidth;
	      }
	    }
          }
        }

        cifo << "     save_" <<endl;
        cifo <<endl;
        // loop for all items in the category
        targetitem.clear();
        targetitem.push_back(categoryName);
        listOutItem.clear();
        // Searchs in the table "item" all items with "category_id"
        // equall to particular category
        itemtbl->Search(listOutItem, targetitem, listitem2);
        if (!listOutItem.empty())
        {
            for (i4=0; i4<(int)(listOutItem.size()); i4++)
            {
              // VLAD - Potential bug ?? should the second argument be "name"
	      //itemtbl->GetCell(itemName,string("name"),listOutItem[i4]);
              itemName = (*itemtbl)(listOutItem[i4], "name");

              cifo << "save_" <<itemName << endl;

	      listOut.clear();
	      formatP->Search(listOut, target2, list);
	
              if (!listOut.empty())
              {
                // writes values for particular item like atom_site.aniso_B[1][1]
                //,symmetry.entry_id
                // for all tables that describe item like item, item_description
                // item_default....
	        for (i3=0; i3<(int)(listOut.size());i3++)
                {
	          //format->GetCell(categoryName2,string("catName"),listOut[i3]);
                  categoryName2 = (*formatP)(listOut[i3], "catName");
	          cattbl2=block.GetTablePtr(categoryName2);
	          targetitem.clear();
	          targetitem.push_back(categoryName2);
	          targetitem.push_back("implicit");
	          listOutItem2.clear();
	          // searchs in table item in ddl all rows with <categoryName2> in
	          // column "category_id" and "implicit" in column "mandatory_code"
	          // assumption : there is only one row with this characteristic
	          // in column name is name of implicit item for this category
	          // later this category will be searched to find
	          // particular item i column name found in this search
	          itemtblddl->Search(listOutItem2, targetitem, listitem);
	          if (!listOutItem2.empty())
                  {
                    //itemtblddl->GetCell(itemName2,string("name"),listOutItem2[0]);
                    itemName2 = (*itemtblddl)(listOutItem2[0], "name");
	            CifString::GetItemFromCifItem(item2, itemName2);
	            listcat2.clear();
	            listcat2.push_back(item2);
	
	            targetcat2.clear();
	            targetcat2.push_back(itemName);
	            listOutcat2.clear();
	            // search one of table related to item
	            cattbl2->Search(listOutcat2, targetcat2, listcat2);
	            // here is part to write value(s) for thiscategory
	            if (!listOutcat2.empty())
                    {
	              if (listOutcat2.size() == 1)
                      {
	                // pair
	                numColumn = cattbl2->GetNumColumns();
                        const vector<string>& colNames =
                          cattbl2->GetColumnNames();
                        const vector<string>& rowValues =
	                  cattbl2->GetRow(listOutcat2[0]);
	                cwid = -1;
	                for (unsigned int i=0; i< numColumn; i++)
                        {
	                  ilen = colNames[i].size();
	                  if (  ilen > cwid) cwid=ilen;
	                }
	    
	                for (unsigned int i=0; i< numColumn; i++)
                        {
	                  linePos=8; //**************
                          _PrintItemIdent(cifo, linePos);
	                  _PrintItemName(cifo, categoryName2, colNames[i], linePos);
                          _PrintPostItemSeparator(cifo, linePos, true);
	                  _PrintItemValue(cifo, rowValues[i], linePos);
	                  if (linePos != 0) cifo << endl;
	                }
	              }
	  
	              if (listOutcat2.size() > 1)
                      {
	                // loop
	                cifo << "     loop_" << endl;
	                numColumn = cattbl2->GetNumColumns();
                        const vector<string>& colNames =
                          cattbl2->GetColumnNames();
	                for (unsigned int i=0; i< numColumn; i++)
                        {
                          _PrintItemIdent(cifo, linePos);
	                  _PrintItemName(cifo, categoryName2, colNames[i], linePos);
                          _PrintPostItemSeparator(cifo, linePos, true);
	                  cifo << endl;
	                }
	    
	                // ---
	                cwidth  = new int[numColumn];
	    
      	                for (unsigned int i=0; i< numColumn; i++) { cwidth[i]=-1+2; }
	                for (unsigned int j=0; j < listOutcat2.size(); ++j)
                        {
                            const vector<string>& rowValues =
	                      cattbl2->GetRow(listOutcat2[j]);
	                    for (unsigned int i=0; i< numColumn; i++)
                            {
		              ilen = rowValues[i].size();
		              if ( _IsQuotableText(rowValues[i])) ilen+=2;
		
		              if (  ilen > cwidth[i]) cwidth[i]=ilen;
	                }
	              }
	              for (unsigned int j=0; j < listOutcat2.size(); ++j)
                      {
	                linePos = 0;
                        const vector<string>& rowValues =
	                  cattbl2->GetRow(listOutcat2[j]);
	                for (unsigned int i=0; i< numColumn; i++)
                        {
		
		          ilen=_PrintItemValue(cifo, rowValues[i], linePos,eNONE,cwidth[i]); 
	                }
	                if (linePos != 0) cifo << endl;
	              }
	    
	              if (cwidth) delete[] cwidth;
	            }
	          }
	        }
	      }
            }
	    cifo << "     save_" <<endl;
	    cifo <<endl;
          }
        }
      }
    }
  }

  return(1); 

}


void DicFile::Compress(CifFile* ddl)
{

    ISTable* catkey;
    ISTable* tbl;
    vector<string> catList;
    vector<string> catList2;
    vector<string> catkeyList;
    vector<unsigned int> listOut2;
    vector<string> catTarget;
    vector<string> catkeyTarget;

    string cell;
    string categoryName;
    string itemName;
    string blockName;

    vector<pair<unsigned int, unsigned int> > duplRows;

    blockName = ddl->GetFirstBlockName();
    Block& firstBlock = ddl->GetBlock(blockName);

    for (unsigned int ib = 0; ib < _blocks.size(); ++ib)
    {
        Block& block = GetBlock(_blocks[ib].GetName());

        catkey = firstBlock.GetTablePtr("category_key");

        catkeyList.clear();
        catkeyList.push_back("id");

        vector<string> tableNames;
        block.GetTableNames(tableNames);

        for (unsigned int jt = 0; jt < tableNames.size(); ++jt)
        {
            tbl = block.GetTablePtr(tableNames[jt]);

            catkeyTarget.clear();
            catkeyTarget.push_back(tableNames[jt]);

            catkey->Search(listOut2, catkeyTarget, catkeyList);
            if (!listOut2.empty())
            {
                catList2.clear();
                for (unsigned int i = 0; i < listOut2.size(); ++i)
                {
                    //catkey->GetCell(cell, string("name"), listOut2[i]);
                    cell = (*catkey)(listOut2[i], "name");
                    CifString::GetItemFromCifItem(itemName, cell);
                    catList2.push_back(itemName);
                }
                // Third argument (false) indicates that the duplicate rows
                // will not be kept, i.e., they will be deleted

#ifndef VLAD_PRINT_DUPL_ROWS
                tbl->FindDuplicateRows(duplRows, catList2, false);
#else
                tbl->FindDuplicateRows(duplRows, catList2, true);

                if (!duplRows.empty())
                {
                    cout << "Duplicate rows for table: \"" << tbl->GetName() <<
                      "\"" << endl;

                    const vector<string>& colNames = tbl->GetColumnNames();

                    for (unsigned int colI = 0; colI < colNames.size(); ++colI)
                    {
                           cout << colNames[colI] << "    ";
                    } 
                    cout << endl;

                    for (unsigned int duplI = 0; duplI < duplRows.size();
                      ++duplI)
                    {
                       for (unsigned int colI = 0; colI < colNames.size();
                         ++colI)
                       {
                           cout << (*tbl)(duplRows[duplI].second,
                             colNames[colI]) << "    ";
                       }
                       cout << endl;
                    }
                }
#endif
            }
            listOut2.clear();
        }
    }
}


CifFile* DicFile::GetRefFile()
{

    CifFile* refFileP = new CifFile();
 
    refFileP->AddBlock("ref_block");
    Block& block = refFileP->GetBlock(refFileP->GetFirstBlockName());

    ISTable* catTableP = new ISTable("category");
    catTableP->AddColumn("id");
    catTableP->AddColumn("mandatory_code");
    catTableP->AddColumn("implicit_key");

    ISTable* tableP = catTableP;

    AddRefRow(*tableP, "datablock", "no", "mmcif_ddl.dic");
    AddRefRow(*tableP, "datablock_methods", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category_key", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category_group", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category_group_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "category_methods", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "sub_category", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "sub_category_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "sub_category_methods", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_aliases", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_default", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_dependent", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_description", "yes", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_enumeration", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_linked", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_methods", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_range", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_related", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_structure", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_structure_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_sub_category", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_type", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_type_conditions", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_type_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_units", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_units_conversion", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "item_units_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "method_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "dictionary", "yes", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "dictionary_history", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_category_description", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_category_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item_description", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item_enumeration", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item_range", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item_type", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "ndb_item", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_category_context", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_context", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_linked_group", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_linked_group_list", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_range", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_type", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_category_description", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_category_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_description", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_enumeration", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_item_examples", "no", "mmcif_ddl.dic"); 
    AddRefRow(*tableP, "pdbx_dictionary_component", "no", "mmcif_ddl.dic");
    AddRefRow(*tableP, "pdbx_dictionary_component_history", "no", "mmcif_ddl.dic");     


    ISTable* itemTableP = new ISTable("item");
    itemTableP->AddColumn("name");
    itemTableP->AddColumn("category_id");
    itemTableP->AddColumn("mandatory_code");

    tableP = itemTableP;

    AddRefRow(*tableP, "_datablock.id","datablock","implicit"); 
    AddRefRow(*tableP, "_datablock.description","datablock","yes");
    AddRefRow(*tableP, "_datablock_methods.datablock_id","datablock_methods","implicit"); 
    AddRefRow(*tableP, "_datablock_methods.method_id","datablock_methods","yes");
    AddRefRow(*tableP, "_category.id","category","yes");
    AddRefRow(*tableP, "_category.description","category","yes");
    AddRefRow(*tableP, "_category.implicit_key","category","implicit"); 
    AddRefRow(*tableP, "_category.mandatory_code","category","yes");
    AddRefRow(*tableP, "_category_examples.id","category_examples","implicit"); 
    AddRefRow(*tableP, "_category_examples.case","category_examples","yes");
    AddRefRow(*tableP, "_category_examples.detail","category_examples","no");
    AddRefRow(*tableP, "_category_key.name","category_key","yes");
    AddRefRow(*tableP, "_category_key.id","category_key","implicit"); 
    AddRefRow(*tableP, "_category_group.id","category_group","yes");
    AddRefRow(*tableP, "_category_group.category_id","category_group"              ,"implicit"); 
    AddRefRow(*tableP, "_category_group_list.id","category_group_list"         ,"yes");
    AddRefRow(*tableP, "_category_group_list.description","category_group_list"         ,"yes");
    AddRefRow(*tableP, "_category_group_list.parent_id","category_group_list"         ,"no");
    AddRefRow(*tableP, "_category_methods.category_id","category_methods"            ,"implicit"); 
    AddRefRow(*tableP, "_category_methods.method_id","category_methods"            ,"yes");
    AddRefRow(*tableP, "_sub_category.id"           ,"sub_category"                ,"yes");
    AddRefRow(*tableP, "_sub_category.description","sub_category"                ,"yes");
    AddRefRow(*tableP, "_sub_category_examples.id","sub_category_examples"       ,"yes");
    AddRefRow(*tableP, "_sub_category_examples.case","sub_category_examples"       ,"yes");
    AddRefRow(*tableP, "_sub_category_examples.detail","sub_category_examples"       ,"no");
    AddRefRow(*tableP, "_sub_category_methods.sub_category_id","sub_category_methods"        ,"yes");
    AddRefRow(*tableP, "_sub_category_methods.method_id"                 ,"sub_category_methods"        ,"yes");
    AddRefRow(*tableP, "_item.name","item"                        ,"implicit"); 
    AddRefRow(*tableP, "_item.mandatory_code","item"                        ,"yes");
    AddRefRow(*tableP, "_item.category_id" ,"item"                        ,"implicit"); 
    AddRefRow(*tableP, "_item_aliases.name","item_aliases"                ,"implicit"); 
    AddRefRow(*tableP, "_item_aliases.alias_name" ,"item_aliases"                ,"yes");
    AddRefRow(*tableP, "_item_aliases.dictionary" ,"item_aliases"                ,"yes");
    AddRefRow(*tableP, "_item_aliases.version"    ,"item_aliases"                ,"yes");
    AddRefRow(*tableP, "_item_default.name"       ,"item_default"                ,"implicit"); 
    AddRefRow(*tableP, "_item_default.value"      ,"item_default"                ,"no");
    AddRefRow(*tableP, "_item_dependent.name"     ,"item_dependent"              ,"implicit"); 
    AddRefRow(*tableP, "_item_dependent.dependent_name" ,"item_dependent"  ,"yes");
    AddRefRow(*tableP, "_item_description.name"         ,"item_description","implicit"); 
    AddRefRow(*tableP, "_item_description.description"  ,"item_description","yes");
    AddRefRow(*tableP, "_item_enumeration.name"         ,"item_enumeration","implicit");
    AddRefRow(*tableP, "_item_enumeration.value"        ,"item_enumeration","yes");
    AddRefRow(*tableP, "_item_enumeration.detail"       ,"item_enumeration","no");
    AddRefRow(*tableP, "_item_enumeration.pdbx_value_display" ,"item_enumeration","no");
    AddRefRow(*tableP, "_item_examples.name"            ,"item_examples"   ,"implicit"); 
    AddRefRow(*tableP, "_item_examples.case"            ,"item_examples"   ,"yes");
    AddRefRow(*tableP, "_item_examples.detail"          ,"item_examples"   ,"no");
    AddRefRow(*tableP, "_item_linked.child_name"        ,"item_linked"     ,"yes");
    AddRefRow(*tableP, "_item_linked.parent_name"       ,"item_linked"     ,"implicit"); 
    AddRefRow(*tableP, "_item_methods.name"             ,"item_methods"    ,"implicit"); 
    AddRefRow(*tableP, "_item_methods.method_id"        ,"item_methods"    ,"yes");
    AddRefRow(*tableP, "_item_range.name"               ,"item_range"      ,"implicit"); 
    AddRefRow(*tableP, "_item_range.minimum"            ,"item_range"      ,"yes");
    AddRefRow(*tableP, "_item_range.maximum"            ,"item_range"      ,"yes");
    AddRefRow(*tableP, "_item_range.ordinal"            ,"item_range"      ,"implicit-ordinal");
    AddRefRow(*tableP, "_item_related.name"             ,"item_related"    ,"implicit"); 
    AddRefRow(*tableP, "_item_related.related_name"     ,"item_related"    ,"yes");
    AddRefRow(*tableP, "_item_related.function_code"    ,"item_related"    ,"yes");
    AddRefRow(*tableP, "_item_structure.name"           ,"item_structure"  ,"implicit"); 
    AddRefRow(*tableP, "_item_structure.code"           ,"item_structure"  ,"yes");
    AddRefRow(*tableP, "_item_structure.organization"   ,"item_structure"  ,"yes");
    AddRefRow(*tableP, "_item_structure_list.code"      ,"item_structure_list" ,"yes");
    AddRefRow(*tableP, "_item_structure_list.index"     ,"item_structure_list" ,"yes");
    AddRefRow(*tableP, "_item_structure_list.dimension" ,"item_structure_list" ,"yes");
    AddRefRow(*tableP, "_item_sub_category.name"        ,"item_sub_category"   ,"implicit"); 
    AddRefRow(*tableP, "_item_sub_category.id"          ,"item_sub_category"   ,"yes");
    AddRefRow(*tableP, "_item_type.name"                ,"item_type"           ,"implicit"); 
    AddRefRow(*tableP, "_item_type.code"                ,"item_type"           ,"yes");
    AddRefRow(*tableP, "_item_type_conditions.name"     ,"item_type_conditions","implicit"); 
    AddRefRow(*tableP, "_item_type_conditions.code"     ,"item_type_conditions","yes");
    AddRefRow(*tableP, "_item_type_list.code"           ,"item_type_list"      ,"yes");
    AddRefRow(*tableP, "_item_type_list.primitive_code" ,"item_type_list"      ,"yes");
    AddRefRow(*tableP, "_item_type_list.construct"      ,"item_type_list"      ,"no");
    AddRefRow(*tableP, "_item_type_list.detail"         ,"item_type_list"      ,"no");
    AddRefRow(*tableP, "_item_units.name"               ,"item_units"          ,"implicit"); 
    AddRefRow(*tableP, "_item_units.code"               ,"item_units"          ,"yes");
    AddRefRow(*tableP, "_item_units_conversion.from_code" ,"item_units_conversion","yes");
    AddRefRow(*tableP, "_item_units_conversion.to_code"   ,"item_units_conversion","yes");
    AddRefRow(*tableP, "_item_units_conversion.operator"  ,"item_units_conversion","yes");
    AddRefRow(*tableP, "_item_units_conversion.factor"    ,"item_units_conversion","yes");
    AddRefRow(*tableP, "_item_units_list.code"            ,"item_units_list"      ,"yes");
    AddRefRow(*tableP, "_item_units_list.detail"          ,"item_units_list"      ,"no");
    AddRefRow(*tableP, "_method_list.id"                  ,"method_list"          ,"yes");
    AddRefRow(*tableP, "_method_list.detail"              ,"method_list"          ,"no");
    AddRefRow(*tableP, "_method_list.inline"              ,"method_list"          ,"yes");
    AddRefRow(*tableP, "_method_list.code"                ,"method_list"          ,"yes");
    AddRefRow(*tableP, "_method_list.language"            ,"method_list"          ,"yes");
    AddRefRow(*tableP, "_dictionary.datablock_id"         ,"dictionary"           ,"implicit"); 
    AddRefRow(*tableP, "_dictionary.title"                ,"dictionary"           ,"yes");
    AddRefRow(*tableP, "_dictionary.version"              ,"dictionary"           ,"yes");
    AddRefRow(*tableP, "_dictionary_history.version"      ,"dictionary_history"   ,"yes");
    AddRefRow(*tableP, "_dictionary_history.update"       ,"dictionary_history"   ,"yes");
    AddRefRow(*tableP, "_dictionary_history.revision"     ,"dictionary_history"   ,"yes");
    AddRefRow(*tableP, "_ndb_category_description.id"     ,"ndb_category_description","implicit"); 
    AddRefRow(*tableP, "_ndb_category_description.description" ,"ndb_category_description"    ,"yes");
    AddRefRow(*tableP, "_ndb_category_examples.id"          ,"ndb_category_examples" ,"implicit"); 
    AddRefRow(*tableP, "_ndb_category_examples.case"        ,"ndb_category_examples" ,"yes");
    AddRefRow(*tableP, "_ndb_category_examples.detail"      ,"ndb_category_examples" ,"no");
    AddRefRow(*tableP, "_ndb_item_description.name"         ,"ndb_item_description"  ,"implicit"); 
    AddRefRow(*tableP, "_ndb_item_description.description"  ,"ndb_item_description"  ,"yes");
    AddRefRow(*tableP, "_ndb_item_enumeration.name"         ,"ndb_item_enumeration"  ,"implicit"); 
    AddRefRow(*tableP, "_ndb_item_enumeration.value"        ,"ndb_item_enumeration"  ,"yes");
    AddRefRow(*tableP, "_ndb_item_enumeration.detail"       ,"ndb_item_enumeration"  ,"no");
    AddRefRow(*tableP, "_ndb_item_examples.case"            ,"ndb_item_examples"     ,"yes");
    AddRefRow(*tableP, "_ndb_item_examples.detail"          ,"ndb_item_examples"     ,"yes");
    AddRefRow(*tableP, "_ndb_item_examples.name"            ,"ndb_item_examples"     ,"implicit"); 
    AddRefRow(*tableP, "_ndb_item_range.ordinal"            ,"ndb_item_range"     ,"implicit-ordinal");
    AddRefRow(*tableP, "_ndb_item_range.name"            ,"ndb_item_range"     ,"implicit");
    AddRefRow(*tableP, "_ndb_item_range.minimum"            ,"ndb_item_range"     ,"yes");
    AddRefRow(*tableP, "_ndb_item_range.maximum"            ,"ndb_item_range"     ,"yes");
    AddRefRow(*tableP, "_ndb_item_type.name"            ,"ndb_item_type"     ,"implicit");
    AddRefRow(*tableP, "_ndb_item_type.code"            ,"ndb_item_type"     ,"yes");
    AddRefRow(*tableP, "_ndb_item.name"            ,"ndb_item"     ,"implicit");
    AddRefRow(*tableP, "_ndb_item.mandatory_code"            ,"ndb_item"     ,"yes");
    AddRefRow(*tableP, "_pdbx_category_context.category_id"            ,"pdbx_category_context"     ,"implicit");
    AddRefRow(*tableP, "_pdbx_category_context.type"            ,"pdbx_category_context"     ,"yes");
    AddRefRow(*tableP, "_pdbx_item_context.item_name"            ,"pdbx_item_context"     ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_context.type"            ,"pdbx_item_context"     ,"yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group.category_id","pdbx_item_linked_group","yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group.link_group_id","pdbx_item_linked_group","yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group.label"        ,"pdbx_item_linked_group","yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group.context"      ,"pdbx_item_linked_group","no");
    AddRefRow(*tableP, "_pdbx_item_linked_group.condition_id" ,"pdbx_item_linked_group","no");
    AddRefRow(*tableP, "_pdbx_item_linked_group_list.child_category_id","pdbx_item_linked_group_list" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group_list.link_group_id" ,"pdbx_item_linked_group_list" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group_list.child_name"    ,"pdbx_item_linked_group_list" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group_list.parent_name"   ,"pdbx_item_linked_group_list" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_linked_group_list.parent_category_id" ,"pdbx_item_linked_group_list" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_range.ordinal" ,"pdbx_item_range" ,"implicit-ordinal");
    AddRefRow(*tableP, "_pdbx_item_range.name" ,"pdbx_item_range" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_range.minimum" ,"pdbx_item_range" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_range.maximum" ,"pdbx_item_range" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_type.name" ,"pdbx_item_type" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_type.code" ,"pdbx_item_type" ,"yes");
    AddRefRow(*tableP, "_pdbx_item.name" ,"pdbx_item" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item.mandatory_code" ,"pdbx_item" ,"yes");
    AddRefRow(*tableP, "_pdbx_category_description.id" ,"pdbx_category_description" ,"implicit");
    AddRefRow(*tableP, "_pdbx_category_description.description" ,"pdbx_category_description" ,"yes");
    AddRefRow(*tableP, "_pdbx_category_examples.id" ,"pdbx_category_examples" ,"implicit");
    AddRefRow(*tableP, "_pdbx_category_examples.case" ,"pdbx_category_examples" ,"yes");
    AddRefRow(*tableP, "_pdbx_category_examples.detail" ,"pdbx_category_examples" ,"no");
    AddRefRow(*tableP, "_pdbx_item_description.name" ,"pdbx_item_description" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_description.description" ,"pdbx_item_description" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_enumeration.name" ,"pdbx_item_enumeration" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_enumeration.value" ,"pdbx_item_enumeration" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_enumeration.detail" ,"pdbx_item_enumeration" ,"no");
    AddRefRow(*tableP, "_pdbx_item_enumeration.value_display" ,"pdbx_item_enumeration" ,"no");    
    AddRefRow(*tableP, "_pdbx_item_examples.name" ,"pdbx_item_examples" ,"implicit");
    AddRefRow(*tableP, "_pdbx_item_examples.case" ,"pdbx_item_examples" ,"yes");
    AddRefRow(*tableP, "_pdbx_item_examples.detail" ,"pdbx_item_examples" ,"no");
    AddRefRow(*tableP, "_pdbx_dictionary_component.datablock_id" ,"pdbx_dictionary_component" ,"implicit");
    AddRefRow(*tableP, "_pdbx_dictionary_component.dictionary_component_id" ,"pdbx_dictionary_component" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component.title" ,"pdbx_dictionary_component" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component.version" ,"pdbx_dictionary_component" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component_history.version" ,"pdbx_dictionary_component_history" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component_history.dictionary_component_id" ,"pdbx_dictionary_component" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component_history.update" ,"pdbx_dictionary_component_history" ,"yes");
    AddRefRow(*tableP, "_pdbx_dictionary_component_history.revision" ,"pdbx_dictionary_component_history" ,"yes");
    
    block.WriteTable(catTableP);
    block.WriteTable(itemTableP);

    return (refFileP);
}


void DicFile::AddRefRow(ISTable& table, const char* first, const char* second,
  const char* third)
{
    vector<string> row;
    row.push_back(first);
    row.push_back(second);
    row.push_back(third);

    table.AddRow(row);
}


void DicFile::WriteItemAliases(const string& cifFileName)
{

    ofstream cifo(cifFileName.c_str(), ios::out | ios::trunc);

    WriteItemAliases(cifo);

    // VLAD - ERROR HANDLING: file will not be closed if WriteItemAliases
    // throws exception, as it can.

    cifo.close();

}


void DicFile::WriteItemAliases(ostream& cifo)
{

    Block& block = GetBlock(GetFirstBlockName());

    ISTable* tbl = block.GetTablePtr("item_aliases");

    if (tbl == NULL)
    {
        return;
    }

    unsigned int numRow = tbl->GetNumRows();

    for (unsigned int rowI = 0; rowI < numRow; ++rowI)
    {
        cifo.setf(ios::left, ios::adjustfield);
        cifo.width(50);
        cifo << (*tbl)(rowI, "name");

        cifo.setf(ios::left, ios::adjustfield);
        cifo.width(30);
        cifo << (*tbl)(rowI, "alias_name");
        cifo << endl;
    }

}


#ifdef VLAD_DO_WE_NEED_THIS
// If we fix types and report in CheckTypesAndRe...
void DicFile::CheckParentChildTypes(const string& diagFileName)
{
    CifDataInfo cifDataInfo(*this);

    Block& block = this->GetBlock(this->GetFirstBlockName());

    CifParentChild cifParentChild(block);


    //
    vector<string> categories = cifDataInfo.GetCatNames();
    sort(categories.begin(), categories.end());

    // For all categories in the dictionary
    for (unsigned int catI = 0; catI < categories.size(); ++catI)
    {
        const string& catName = categories[catI];

        const vector<vector<string> >& origParComboKeys =
          cifParentChild.GetComboKeys(catName);

        for (unsigned int keyI = 0; keyI < origParComboKeys.size(); ++keyI)
        {
            const vector<string>& currOrigParComboKey = origParComboKeys[keyI];

            vector<string> parKeyTypeCodes;

            // Get parent key types
            for (unsigned int parKeyI = 0;
              parKeyI < currOrigParComboKey.size(); ++parKeyI)
            {
                const vector<string>& typeCodes =
                  cifDataInfo.GetItemAttribute(currOrigParComboKey[parKeyI],
                    "item_type", "code");

                if (typeCodes.empty())
                    cout << "NEW - ERROR - Empty type code for item \"" <<
                      currOrigParComboKey[parKeyI] << "\"" << endl;
               
                if (typeCodes.size() > 1)
                    cout <<
                      "NEW - ERROR - More than one type code for item \"" <<
                      currOrigParComboKey[parKeyI] << "\"" << endl;

                parKeyTypeCodes.push_back(typeCodes[0]);

                //cout << "NEW - INFO - Item \"" <<
                //  currOrigParComboKey[parKeyI] << "\" has type code \"" <<
                //  typeCodes[0] << "\"" << endl;
            }

            vector<vector<vector<string> > >& origChildrenKeys =
              cifParentChild.GetChildrenKeys(currOrigParComboKey);

            for (unsigned int childI = 0; childI < origChildrenKeys.size();
              ++childI)
            {
                for (unsigned int childKeyI = 0; childKeyI <
                  origChildrenKeys[childI].size(); ++childKeyI)
                {
                    const vector<string>& currChKey =
                      origChildrenKeys[childI][childKeyI];

                    for (unsigned int chKeyI = 0; chKeyI < currChKey.size();
                      ++chKeyI)
                    {
                        const vector<string>& chTypeCodes =
                          cifDataInfo.GetItemAttribute(currChKey[chKeyI],
                            "item_type", "code");

                        if (chTypeCodes.empty())
                            cout << "NEW - ERROR - Empty type code for item \"" <<
                              currChKey[chKeyI] << "\"" << endl;
               
                        if (chTypeCodes.size() > 1)
                            cout <<
                              "NEW - ERROR - More than one type code for item \"" <<
                              currChKey[chKeyI] << "\"" << endl;
                        if (chTypeCodes[0] != parKeyTypeCodes[chKeyI])
                            cout <<
                              "NEW - ERROR - Key mismatch between child item \"" << currChKey[chKeyI] << "\" with type code \"" << chTypeCodes[0] <<
"\" and parent item \"" << currOrigParComboKey[chKeyI] << "\" with type code \"" << parKeyTypeCodes[chKeyI] << "\"" << endl;
                    }
                }
            } // For every child key compare its types to the parent types
            // eTypeCode iType = _dataInfo._GetDataType(itemName);
        }
    }
    //

}
#endif

