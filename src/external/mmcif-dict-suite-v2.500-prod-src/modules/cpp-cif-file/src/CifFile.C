//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file CifFile.C
**
** \brief Implementation file for CifFile class.
*/


// VLAD - Carefully examine all methods that accept with, for cases where
// with is 0. This may happen if the value is empty !!

#include <stdexcept>
#include <set>
#include <algorithm>

#include "GenString.h"
#include "RcsbFile.h"
#include "CifString.h"
#include "regex.h"
#include "CifParentChild.h"
#include "CifFile.h"


using std::exception;
using std::runtime_error;
using std::string;
using std::make_pair;
using std::vector;
using std::set;
using std::ofstream;
using std::ostream;
using std::ostringstream;
using std::ios;
using std::endl;


CifFile::CifFile(const eFileMode fileMode, const string& objFileName,
  const bool verbose, const Char::eCompareType caseSense,
  const unsigned int maxLineLength, const string& nullValue) :
  TableFile(fileMode, objFileName, caseSense)
{

    Init();

    _maxCifLineLength = maxLineLength;
    _nullValue = nullValue;
    _verbose = verbose;

}


CifFile::CifFile(const bool verbose, const Char::eCompareType caseSense,
  const unsigned int maxLineLength, const string& nullValue) :
  TableFile(caseSense)
{

    Init();

    _maxCifLineLength = maxLineLength;
    _nullValue = nullValue;
    _verbose = verbose;

}


CifFile::CifFile(const bool fake, const bool verbose,
  const unsigned int intCaseSense, const unsigned int maxLineLength,
  const string& nullValue) :
  TableFile((Char::eCompareType)intCaseSense)
{

    Init();

    _maxCifLineLength = maxLineLength;
    _nullValue = nullValue;
    _verbose = verbose;

}


CifFile::~CifFile()
{

}


void CifFile::Init()
{

  _beginDataKeyword = "data_";
  _endDataKeyword = "# ";

  _beginLoopKeyword = "loop_";
  _endLoopKeyword = "";

  _maxCifLineLength = STD_CIF_LINE_LENGTH;
  _nullValue = CifString::UnknownValue;
  _verbose = false;
  _smartPrint = false;
  _quotes = "\'";
  _enumCaseSense = false;

}


void CifFile::SetSrcFileName(const string& srcFileName)
{
    // VLAD - TODO check for file existence and throw exception
    _srcFileName = srcFileName;
}


const string& CifFile::GetSrcFileName()
{
    return (_srcFileName);
}


void CifFile::SetQuoting(eQuoting quoting)
{

    if (quoting == eSINGLE)
    {
        _quotes = "\'";
    }
    else if (quoting == eDOUBLE)
    {
        _quotes = "\"";
    } 
    // VLAD - ERROR HANDLING
}


unsigned int CifFile::GetQuoting()
{

    if (_quotes == "\'")
        return(eSINGLE);
    else if (_quotes == "\"")
        return(eDOUBLE);
    else
        // VLAD - ERROR HANDLING
        return(eDOUBLE);
       
}


void CifFile::SetLooping(const string& catName, bool looping)
{
    _looping[catName] = looping;
}


bool CifFile::GetLooping(const string& catName)
{
    return (_looping[catName]);
}


void CifFile::Write(const string& cifFileName, const bool sortTables,
  const bool writeEmptyTables)
{

    ofstream cifo(cifFileName.c_str(), ios::out | ios::trunc);

    Write(cifo, sortTables, writeEmptyTables);

    cifo.close();

}


void CifFile::Write(ostream& cifo, const bool sortTables,
  const bool writeEmptyTables)
{

    vector<unsigned int> tablesIndices;

    if (sortTables)
    {
        GetSortedTablesIndices(tablesIndices);
    }
    else
    {
        GetTablesIndices(tablesIndices);
    }

    Write(cifo, tablesIndices, writeEmptyTables);

}


void CifFile::Write(const string& cifFileName, const vector<string>& catOrder,
  const bool writeEmptyTables)
{

    ofstream cifo(cifFileName.c_str(), ios::out | ios::trunc);

    Write(cifo, catOrder, writeEmptyTables);

    cifo.close();

}


void CifFile::Write(ostream& cifo, const vector<string>& catOrder,
  const bool writeEmptyTables)
{

    vector<unsigned int> tables;
    unsigned int tableIndex = 0;

    // tables is a list broken into segments of # tables/data block 
    // providing the index of tables to output

    // VLAD : Can this be improved by first searching the category and
    // then doing something about it? Well maybe this is not doable.
    for (unsigned int blockI = 0; blockI < _blocks.size();
      blockI++)
    {
        vector<bool> processedTables(_blocks[blockI]._tables.size());

        for (unsigned int l = 0; l < catOrder.size(); l++)
        {
            tableIndex = _blocks[blockI]._tables.\
              find(catOrder[l]);

            if (tableIndex != _blocks[blockI]._tables.size())
            {
                tables.push_back(tableIndex);
                processedTables[tableIndex] = true;
            }
            else
            {
                processedTables[tableIndex] = false;
            }
        }

	// Now fill in the tables index that were missed -- critical for multi datablock files
	// that tables are 
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            if (!processedTables[tableI])
            {
                tables.push_back(tableI);
            }
        }

    } // for blockI

    Write(cifo, tables, writeEmptyTables);

}


void CifFile::WriteNmrStar(const string& nmrStarFileName,
  const string& globalBlockName, const bool sortTables,
  const bool writeEmptyTables)
{

    string savedBeginDataKeyword = _beginDataKeyword;
    string savedEndDataKeyword = _endDataKeyword;

    string savedBeginLoopKeyword = _beginLoopKeyword;
    string savedEndLoopKeyword = _endLoopKeyword;

    _beginDataKeyword = "save_";
    _endDataKeyword = "save_";

    _beginLoopKeyword = "loop_";
    _endLoopKeyword = "stop_";

    ofstream cifo(nmrStarFileName.c_str(), ios::out | ios::trunc);

    cifo << "data_" << globalBlockName << endl;
    cifo << endl;

    Write(cifo, sortTables, writeEmptyTables);

    cifo.close();

    _beginDataKeyword = savedBeginDataKeyword;
    _endDataKeyword = savedEndDataKeyword;

    _beginLoopKeyword = savedBeginLoopKeyword;
    _endLoopKeyword = savedEndLoopKeyword;

}


int CifFile::DataChecking(CifFile& ref, const string& diagFileName,
  const bool extraDictChecks, const bool extraCifChecks,
  const std::vector<std::string>& skipBlockNames)
{
    _extraDictChecks = extraDictChecks;
    _extraCifChecks = extraCifChecks;

    int ret = 0;

    ofstream log;
    log.open(diagFileName.c_str(), ios::out | ios::app);

    vector<string> refBlockNames;
    ref.GetBlockNames(refBlockNames);

    if (refBlockNames.size() > 1)
    {
        log << "WARNING - dictionary has " << refBlockNames.size() <<
          " datablocks" << endl;
    }

    Block& refBlock = ref.GetBlock(ref.GetFirstBlockName());

    vector<string> BlockNames;
    GetBlockNames(BlockNames);

    for (unsigned int blockI = 0; blockI < BlockNames.size(); ++blockI)
    {
        Block& block = GetBlock(BlockNames[blockI]);

	if (skipBlockNames.size() > 0) {
	  std::string bName = block.GetName();
	  if (std::find(skipBlockNames.begin(), skipBlockNames.end(), bName) != skipBlockNames.end()) {
	      continue;
	  }
	}
	
        ostringstream buf;

        ret = DataChecking(block, refBlock, buf, extraDictChecks,
          extraCifChecks);

        string sBuf = buf.str();

        if (!sBuf.empty())
        {
            log << sBuf << endl;
        }
    }

    if (RcsbFile::IsEmpty(log))
    {
        log.close();
        RcsbFile::Delete(diagFileName);
    }
    else
    {
        log.close();
    }

    return(ret);
}


int CifFile::DataChecking(Block& block, Block& refBlock, ostringstream& buf,
  const bool extraDictChecks, const bool extraCifChecks)
{

    _extraDictChecks = extraDictChecks;
    _extraCifChecks = extraCifChecks;

    int ret1, ret2;

    ret1 = CheckCategories(block, refBlock, buf);
    CheckCategoryKey(block, buf);
    CheckItemsTable(block, buf);
    ret2 = CheckItems(block, refBlock, buf);

    CheckAndRectifyItemTypeCode(block, buf);

    if (ret1 == 0)
        return ret2;
    else
        return ret1;

}


void CifFile::SetEnumCheck(bool caseSense)
{
    _enumCaseSense = caseSense;
}


bool CifFile::GetEnumCheck()
{
    return(_enumCaseSense);
}


const string& CifFile::GetParsingDiags()
{
    return (_parsingDiags);
}


int CifFile::_IsQuotableText(const string& itemValue)
{

  if (itemValue.empty())
      return(0);
    
  if (itemValue[0] == '_')
      return(1);

  for (unsigned int i = 0; i < itemValue.size(); i++)
  {
      if (itemValue[i] == ' ')
          return(1);

      if (itemValue[i] == '\t')
          return(1);
      
      if (itemValue[i] == '\n')
          return(1);

      if (itemValue[i] == '\'' || itemValue[i] == '\"')
          return(1);

      if (i == 0)
      {
          if (CifString::IsSpecialFirstChar(itemValue[i]))
          {
              return(1);
          }
      }
      else
      {
          if (CifString::IsSpecialChar(itemValue[i]))
          {
              return(1);
          }
      }
  }

  if (String::IsCiEqual(itemValue.substr(0, 5), "data_"))
      return(1);
  if (String::IsCiEqual(itemValue.substr(0, 5), "loop_"))
      return(1);
  if (String::IsCiEqual(itemValue.substr(0, 5), "save_"))
      return(1);
  if (String::IsCiEqual(itemValue.substr(0, 5), "stop_"))
      return(1);
  if (String::IsCiEqual(itemValue.substr(0, 7), "global_"))
      return(1);

  return(0);

}


CifFile::eIdentType CifFile::_FindPrintType(const vector<string>& values)
{

    // If at least one value is a number, all should be right justified.
    for (unsigned int valuesI = 0; valuesI < values.size(); ++valuesI)
    {
        if (String::IsNumber(values[valuesI]))
        {
            return(eRIGHT);
        }
    }

    return(eLEFT);

}


void CifFile::_PrintItemIdent(ostream& cifo, unsigned int& linePos)
{

    string ItemNameIdent = "    ";
 
    cifo << ItemNameIdent;

    linePos = ItemNameIdent.size();

}


void CifFile::_PrintItemName(ostream& cifo, const string& category,
  const string& keyword, unsigned int& linePos)
{

    string cifItem;

    CifString::MakeCifItem(cifItem, category, keyword);

    cifo << cifItem;

    linePos += cifItem.size();

}


void CifFile::_PrintPostItemSeparator(ostream& cifo, unsigned int& linePos,
  const bool ident, const unsigned int numSpaces)
{

    const string Space = " ";

    for (unsigned int spaceI = 0; spaceI < numSpaces; spaceI++)
    {
        cifo << Space;
    }

    linePos += numSpaces * Space.size();

    if (!ident)
    {
        return;
    }

    unsigned int start = 36;

    while (linePos > start)
        start += 10;

    for (unsigned int i = 0; i < start - linePos; i++)
        cifo << Space;

    linePos = start;

}


int CifFile::_PrintItemValue(ostream& cifo, const string& itemValue,
  unsigned int& linePos, const eIdentType identType,
  const unsigned int width)
{

    bool multipleLine, multipleWord, embeddedQuotes, embeddedSingleQuotes,
      embeddedDoubleQuotes;
    string Ident;
    unsigned int N = 0;
    // VLAD: with has different meaning for eNONE and others.

    if ((identType == eNONE) && (width != 0))
    {
        Ident = "          ";
    }

    if (linePos == 0)
    {
        cifo << Ident;
        linePos = Ident.size();
    }
  
    if (itemValue.empty())
    {
        if (_maxCifLineLength <= linePos + 2)
        {
            cifo << endl;
            linePos = 0;
        }

        N = _nullValue.size();

        if (identType == eRIGHT)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemValue");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }

        cifo << _nullValue;
        linePos += 1;

        if (identType == eLEFT)
        {
            if (linePos != 0)
            {
                if (width < N)
                {
                    throw runtime_error("CRITICAL ERROR IN: "\
                      "CifFile::_PrintItemValue");
                }

                if (linePos + width - N < _maxCifLineLength)
                {
                    for (unsigned int k = 0; k < width - N; k++)
                        cifo << " ";
                    linePos += width - N;
                }
            }
        }

        cifo << " ";
        linePos += 1;

        return(1);
    }
 
    unsigned int str_len = itemValue.size();

    multipleLine = false;
    multipleWord = false;
    embeddedQuotes = false;
    embeddedSingleQuotes = false;
    embeddedDoubleQuotes = false;
    bool specialChars = false;

    string multipleWordQuotes = _quotes;

    for (unsigned int i = 0; i < str_len; i++)
    {
        if (itemValue[i] == ' ' || itemValue[i] == '\t')
            multipleWord = true;
        else if (itemValue[i] == '\n')
            multipleLine = true;
        else if (itemValue[i] == '\'')
        {
            embeddedSingleQuotes = true;
            embeddedQuotes = true;
            multipleWordQuotes = "\"";
        }
        else if (itemValue[i] == '\"')
        {
            embeddedDoubleQuotes = true;
            embeddedQuotes = true;
            multipleWordQuotes = "\'";
        }
        else if (!specialChars)
        {
            if ((i == 0) && (CifString::IsSpecialFirstChar(itemValue[i])))
                specialChars = true;
            if ((i != 0) && (CifString::IsSpecialChar(itemValue[i])))
                specialChars = true;
        }
    }

    if (itemValue[0] == '_')
        multipleWord = true;
    if (itemValue[0] == ';')
        multipleWord = true;

    if (String::IsCiEqual(itemValue.substr(0, 5), "data_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "loop_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "save_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "stop_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 7), "global_"))
        multipleWord = true;

    if (embeddedQuotes && multipleWord)
        multipleLine = true;
    if (embeddedQuotes)
        multipleWord = true;
    if (specialChars)
        multipleWord = true;

    if (embeddedSingleQuotes && embeddedDoubleQuotes)
        multipleWordQuotes = _quotes;

    if (str_len >= _maxCifLineLength || multipleLine)
    {
        if (linePos != 0)
            cifo << endl;
        cifo << ";" << itemValue;
        cifo << endl;
        cifo << ";" << endl;
        linePos = 0;
    }
    else
    {
        if ((!multipleWord && str_len + 2 + linePos > _maxCifLineLength) ||
          (multipleWord && str_len + 4 + linePos > _maxCifLineLength))
        {
            cifo << endl;
            linePos = 0;
            cifo << Ident;
            linePos += Ident.size();
        }

        string fullItemValue;

        if (multipleWord)
        {
            fullItemValue = multipleWordQuotes + itemValue + multipleWordQuotes;
        }
        else
        {
            fullItemValue = itemValue;
        }

        N = fullItemValue.size();

        if (identType == eRIGHT)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemValue");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }

        cifo << fullItemValue;
        linePos += N;

        if ((identType == eNONE) || (identType == eLEFT))
        {
            cifo << " ";
            linePos++;
        }
    }
 
    if ((identType == eNONE) && (!Ident.empty()))
    {
        if (linePos > Ident.size())
        {
            linePos = linePos + width - N;

            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemValue");
            }

            for (unsigned int i = 0; i < width - N; i++)
                cifo << " ";
        }
    }

    if (identType == eLEFT)
    {
        if (linePos != 0)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemValue");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }
    }

    if (identType == eRIGHT)
    {
        cifo << " ";
        linePos++;
    }

    if (multipleWord || multipleLine || embeddedQuotes) 
        return(str_len + 2);
    else
        return(str_len);
 
}


int CifFile::_PrintItemNameInHeader(ostream& cifo, const string& itemValue,
  unsigned int& linePos, const eIdentType identType,
  const unsigned int width)
{

    bool multipleLine, multipleWord, embeddedQuotes, embeddedSingleQuotes,
      embeddedDoubleQuotes;
    string Ident;
    unsigned int N = 0;
    // VLAD: with has different meaning for eNONE and others.

    if ((identType == eNONE) && (width != 0))
    {
        Ident = "          ";
    }

    if (linePos == 0)
    {
        cifo << Ident;
        linePos = Ident.size();
    }
  
    if (itemValue.empty())
    {
        if (_maxCifLineLength <= linePos + 2)
        {
            cifo << endl << "# ";
            linePos = 2;
        }

        N = _nullValue.size();

        if (identType == eRIGHT)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemNameInHeader");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }

        cifo << _nullValue;
        linePos += 1;

        if (identType == eLEFT)
        {
            if (linePos != 0)
            {
                if (width < N)
                {
                    throw runtime_error("CRITICAL ERROR IN: "\
                      "CifFile::_PrintItemNameInHeader");
                }

                if (linePos + width - N < _maxCifLineLength)
                {
                    for (unsigned int k = 0; k < width - N; k++)
                        cifo << " ";
                    linePos += width - N;
                }
            }
        }

        cifo << " ";
        linePos += 1;

        return(1);
    }
  
    unsigned int str_len = itemValue.size();

    multipleLine = false;
    multipleWord = false;
    embeddedQuotes = false;
    embeddedSingleQuotes = false;
    embeddedDoubleQuotes = false;

    string multipleWordQuotes = _quotes;

    for (unsigned int i = 0; i < str_len; i++)
    {
        if (itemValue[i] == ' ')
            multipleWord = true;
        else if (itemValue[i] == '\n')
            multipleLine = true;
        else if (itemValue[i] == '\'')
        {
            embeddedSingleQuotes = true;
            embeddedQuotes = true;
            multipleWordQuotes = "\"";
        }
        else if (itemValue[i] == '\"')
        {
            embeddedDoubleQuotes = true;
            embeddedQuotes = true;
            multipleWordQuotes = "\'";
        }
    }

    if (itemValue[0] == '_')
        multipleWord = true;
    if (itemValue[0] == ';')
        multipleWord = true;

    if (String::IsCiEqual(itemValue.substr(0, 5), "data_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "loop_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "save_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 5), "stop_"))
        multipleWord = true;
    if (String::IsCiEqual(itemValue.substr(0, 7), "global_"))
        multipleWord = true;

    if (embeddedQuotes && multipleWord)
        multipleLine = true;
    if (embeddedQuotes)
        multipleWord = true;
    if (embeddedSingleQuotes && embeddedDoubleQuotes)
        multipleWordQuotes = _quotes;

    if (str_len >= _maxCifLineLength || multipleLine)
    {
        if (linePos != 0)
            cifo << endl;
        cifo << ";" << itemValue;
        cifo << endl;
        cifo << ";" << endl;
        linePos = 0;
    }
    else
    {
        if ((!multipleWord && str_len + 2 + linePos > _maxCifLineLength) ||
          (multipleWord && str_len + 4 + linePos > _maxCifLineLength))
        {
            cifo << endl << "# ";
            linePos = 2;
            cifo << Ident;
            linePos += Ident.size();
        }

        string fullItemValue;

        if (multipleWord)
        {
            fullItemValue = multipleWordQuotes + itemValue + multipleWordQuotes;
        }
        else
        {
            fullItemValue = itemValue;
        }

        N = fullItemValue.size();

        if (identType == eRIGHT)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemNameInHeader");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }

        cifo << fullItemValue;
        linePos += N;

        if ((identType == eNONE) || (identType == eLEFT))
        {
            cifo << " ";
            linePos++;
        }
    }
 
    if ((identType == eNONE) && (!Ident.empty()))
    {
        if (linePos > Ident.size())
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemNameInHeader");
            }

            linePos = linePos + width - N;

            for (unsigned int i = 0; i < width - N; i++)
                cifo << " ";
        }
    }

    if (identType == eLEFT)
    {
        if (linePos != 0)
        {
            if (width < N)
            {
                throw runtime_error("CRITICAL ERROR IN: "\
                  "CifFile::_PrintItemNameInHeader");
            }

            if (linePos + width - N < _maxCifLineLength)
            {
                for (unsigned int k = 0; k < width - N; k++)
                    cifo << " ";
                linePos += width - N;
            }
        }
    }

    if (identType == eRIGHT)
    {
        cifo << " ";
        linePos++;
    }

    if (multipleWord || multipleLine || embeddedQuotes) 
        return(str_len + 2);
    else
        return(str_len);
 
}



void CifFile::Write(ostream& cifo, vector<unsigned int>& tables,
  const bool writeEmptyTables)
{

    ISTable* tableP = NULL;
 
    unsigned int numColumn, numRow, linePos, cwid, ilen, numPostItemSpaces;

    vector<string> rowValues;
 
    vector<unsigned int> cwidth;
 
    if (_smartPrint)
    {
        numPostItemSpaces = SMART_PRINT_SPACING;
    }
    else
    {
        numPostItemSpaces = STD_PRINT_SPACING;
    }

 
    for (unsigned int blockI = 0, jt = 0;
      blockI < _blocks.size(); blockI++)
    {
        cifo << _beginDataKeyword << _blocks[blockI].GetName() << endl;

        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI, ++jt)
        {
// VLAD - EXAMINE See this same loop in DICFileObjBase or DDLFileObjBase
// VLAD PERFORMANCE. The old BTreeObj code was doing a partial search
// on tids. It searched for "_0_ " in the tree. Even though the exact
// search was never found, it moved the pointer down the tree and reached
// the appropriate "_0_ blah" for example and then processed from that point
// onward the entries that are exact match. This kind of modified proprietory
// search is better than STL. We need to represent this subsearch using
// different STL structures to be able to have better performance than what
// we currently have.

            tableP = _GetTablePtr(blockI, tables[jt]);

            numColumn = tableP->GetNumColumns();
            numRow = tableP->GetNumRows();

            if ((numRow == 0) && (!writeEmptyTables))
                continue;

            cifo << "# " << endl;

            const vector<string>& colNames = tableP->GetColumnNames();

            if ((numRow <= 1) && (!GetLooping(tableP->GetName())))
            {
                // cwid is maxColumnNameSize
                unsigned int longestNameIndex = 0;

                cwid = 0;
                for (unsigned int i = 0; i < numColumn; i++)
                {
                    if (colNames[i].size() > cwid)
                    {
                          cwid = colNames[i].size();
                          longestNameIndex = i;
                    }
                }
 
                string longestCifItem;

                CifString::MakeCifItem(longestCifItem,
                  tableP->GetName(), colNames[longestNameIndex]);

                if (numRow == 0)
                {
                    for (unsigned int i = 0; i < numColumn; i++)
                    {
                        rowValues.push_back(_nullValue);
                    }
                }
                else
                {
                    tableP->GetRow(rowValues, 0);
                }

                for (unsigned int i = 0; i < numColumn; i++)
                {
                    linePos = 0;
                    _PrintItemName(cifo,
                      tableP->GetName(),
                      colNames[i], linePos);
                    _PrintPostItemSeparator(cifo, linePos, false,
                      numPostItemSpaces + colNames[longestNameIndex].size() -
                      colNames[i].size());

                    linePos = longestCifItem.size() + numPostItemSpaces - 1;

                    _PrintItemValue(cifo, rowValues[i], linePos);

                    if (linePos != 0)
                        cifo << endl;
                }
                rowValues.clear();
            }
            else
            {
                cifo << _beginLoopKeyword << endl;

                for (unsigned int i = 0; i < numColumn; i++)
                {
                    linePos = 0;
                    _PrintItemName(cifo,
                      tableP->GetName(),
                      colNames[i], linePos);
                    _PrintPostItemSeparator(cifo, linePos);
                    cifo << endl;
                }

                unsigned int addSpace;

                for (unsigned int i = 0; i < numColumn; i++)
                {
                    if (_smartPrint)
                    {
                        // In smart printing, the column width is determined
                        // based on column values as well as column name
                        if (i == 0)
                        {
                            // "# "
                            addSpace = 2;
                        }
                        else
                        {
                            addSpace = 0;
                        }
                        cwidth.push_back(colNames[i].size() + addSpace);
                    }
                    else
                    { 
                        // Even if values are empty, they should occupy at
                        // least one character.
                        cwidth.push_back(1);
                    }
                }

                for (unsigned int l = 0; l < numRow; l++)
                {
                    tableP->GetRow(rowValues, l);
                    for (unsigned int i = 0; i < numColumn; i++)
                    {
                        ilen = rowValues[i].size();
                        if (_IsQuotableText(rowValues[i]))
                            ilen += 2;
                        if (ilen > cwidth[i])
                        {
                            if ((i == 0) && _smartPrint)
                                cwidth[i] = ilen + 2;
                            else
                                cwidth[i] = ilen;
                        }
                    }
                    rowValues.clear();
                }

                // In smart print mode, determine justification type for
                // each column.
                vector<eIdentType> colPrintType;
                vector<string> col;
                for (unsigned int colI = 0; colI < numColumn; ++colI)
                {
                    if (_smartPrint)
                    {
                        tableP->GetColumn(col, colNames[colI]);

                        colPrintType.push_back(_FindPrintType(col));

                        col.clear();
                    }
                    else
                    {
                        colPrintType.push_back(eLEFT);
                    }
                }

                for (unsigned int l = 0; l < numRow; l++)
                {
                    if (_smartPrint)
                    {
                        if ((l == 0) || ((l % HEADER_SPACING) == 0))
                        {
                            _PrintHeaderedItems(cifo, colNames, cwidth,
                              colPrintType);
                        }
                        // VLAD - PRINT - Add here printing the header
                        // every certain number of lines
                    }
                    linePos = 0;
                    tableP->GetRow(rowValues, l);
                    for (unsigned int i = 0; i < numColumn; i++)
                    {
                        _PrintItemValue(cifo, rowValues[i], linePos,
                          colPrintType[i], cwidth[i]);
                    }

                    if (linePos != 0)
                        cifo << endl;

                    rowValues.clear();
                }

                if (!_endLoopKeyword.empty())
                {
                    cifo << _endLoopKeyword << endl;
                }

                cwidth.clear();
            }
        }

        if (!_endDataKeyword.empty())
        {
            cifo << _endDataKeyword << endl;
        }
    }

}

void CifFile::_PrintHeaderedItems(ostream& cifo,
  const vector<string>& colNames, const vector<unsigned int>& colWidths,
  const vector<eIdentType> colPrintType)
{

    unsigned int linePos = 0;
    unsigned int lessSpace = 0;

    cifo << "# ";

    linePos += 2;

    for (unsigned int colI = 0; colI < colNames.size(); ++colI)
    {
        if (colI == 0)
        {
            lessSpace = 2;
        }
        else
        {
            lessSpace = 0;
        }

        _PrintItemNameInHeader(cifo, colNames[colI], linePos,
          colPrintType[colI], colWidths[colI] - lessSpace);
    }

    cifo << endl;

}


int CifFile::CheckCategories(Block& block, Block& refBlock, ostringstream& log)
{

    int ret = 1;

    ISTable* refCatTableP = refBlock.GetTablePtr("category");

    vector<string> list;
    list.push_back("mandatory_code");

    vector<string> target;
    target.push_back("yes");

    vector<unsigned int> OutList;
    refCatTableP->Search(OutList, target, list);
    if (OutList.empty())
    {
        return(ret);
    }

    for (unsigned int i = 0; i < OutList.size(); i++)
    {
        // For all the mandatory categories
        const string& mandatoryCat = (*refCatTableP)(OutList[i], "id");
        if (!block.IsTablePresent(mandatoryCat))
        {
            log << "ERROR - category \"" << mandatoryCat <<
              "\" is mandatory, but it is not present in datablock \"" <<
              block.GetName() << "\"" << endl;
            ret = 0;
        }
    }

    return(ret);

}


void CifFile::CheckCategoryKey(Block& block, ostringstream& log)
{

    ISTable* catKeyTableP = block.GetTablePtr("category_key");

    ISTable* itemTableP = block.GetTablePtr("item");

    if ((catKeyTableP == NULL) || (itemTableP == NULL))
        return;

    vector<string> searchCols;
    searchCols.push_back("name");

    // For every item in the category key table, check if it is set as
    // mandatory or implicit
    for (unsigned int itemI = 0; itemI < catKeyTableP->GetNumRows(); ++itemI)
    {
        const string& cifItem = (*catKeyTableP)(itemI, "name");

        // Validate key item name
        if (CifString::IsEmptyValue(cifItem))
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", for \"_category_key.id\" == \"" << 
              (*catKeyTableP)(itemI, "id") << "\", invalid value " <<
              "\"_category_key.name\" == \"" << cifItem << "\"" << endl;

            continue;
        }

        vector<string> targets;
        targets.push_back(cifItem);

        unsigned int foundIndex = itemTableP->FindFirst(targets, searchCols);
        if (foundIndex == itemTableP->GetNumRows())
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", item \"" << cifItem << "\" exists in category_key " <<
              " but it is not defined" << endl;
        }
        else
        {
            const string& mandatoryCode = (*itemTableP)(foundIndex,
              "mandatory_code");

            if ((mandatoryCode != "yes") && (mandatoryCode != "implicit") && 
              (mandatoryCode != "implicit-ordinal"))
            {
                log << "ERROR - In block \"" << block.GetName() <<
                  "\", key item \"" << cifItem << "\" is defined as" <<
                  " non-mandatory" << endl;
            }
        }
    }

}


void CifFile::CheckItemsTable(Block& block, ostringstream& log)
{
    ISTable* itemsTableP = block.GetTablePtr("item");

    ISTable* catsTableP = block.GetTablePtr("category");

    if (itemsTableP == NULL)
    {
        // "item" table does not exist. This is the case when checking a CIF
        // file and not dictionary or DDL file.

        return;
    }

#ifdef VLAD_DEBUG
    cout << "Items table: " << *itemsTableP << endl;
#endif

    // Check if "_item.category_id" matches the one in "_item.name"
    for (unsigned int rowI = 0; rowI < itemsTableP->GetNumRows(); ++rowI)
    {
        const string& cifItem = (*itemsTableP)(rowI, "name");

        string catName;
        CifString::GetCategoryFromCifItem(catName, cifItem);

        const string& tableCatName = (*itemsTableP)(rowI, "category_id");

        if (catName != tableCatName)
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", category name mismatch between \"_item.category_id\" == \""
              << tableCatName << "\" and \"_item.name\" == \"" <<
              cifItem << "\" at row# " << rowI + 1 << " in \"item\" table" <<
              endl;
        }

        // See if item's category is defined
        if (!IsCatDefinedInRef(catName, *catsTableP))
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", item \"" << cifItem << "\" does not have its category, \"" <<
              catName << "\", defined." << endl;
        }

    }
}


int CifFile::CheckItems(Block& block, Block& refBlock, ostringstream& log)
{
    int ret = 1;
 
    // Get references to relevant dictionary categories

    ISTable* refCatTableP = refBlock.GetTablePtr("category");

    ISTable* catKeyTableP = refBlock.GetTablePtr("category_key");

    // Set search on "id" column to be case-insensitive, since the "id" item
    // of "category_key" category is implicitly set from the category name of
    // category save frame in the reference file and that category name in
    // category save frame is uppercase.
    catKeyTableP->SetFlags("id", ISTable::DT_STRING | ISTable::CASE_INSENSE);

    ISTable* itemTableP = NULL;
    ISTable* itemTypeTableP = NULL;
    ISTable* itemRangeTableP = NULL;
    ISTable* itemEnumTableP = NULL;

    if (_extraCifChecks)
    {
        itemTableP = refBlock.GetTablePtr("pdbx_item");
        if (itemTableP == NULL)
        {
            itemTableP = refBlock.GetTablePtr("item");
        }

        itemTypeTableP = refBlock.GetTablePtr("pdbx_item_type");
        if (itemTypeTableP == NULL)
        {
            itemTypeTableP = refBlock.GetTablePtr("item_type");
        }

        itemRangeTableP = refBlock.GetTablePtr("pdbx_item_range");
        if (itemRangeTableP == NULL)
        {
            itemRangeTableP = refBlock.GetTablePtr("item_range");
        }

        itemEnumTableP = refBlock.GetTablePtr("pdbx_item_enumeration");
        if (itemEnumTableP == NULL)
        {
            itemEnumTableP = refBlock.GetTablePtr("item_enumeration");
        }
    }
    else
    {
        itemTableP = refBlock.GetTablePtr("item");
        itemTypeTableP = refBlock.GetTablePtr("item_type");
        itemRangeTableP = refBlock.GetTablePtr("item_range");
        itemEnumTableP = refBlock.GetTablePtr("item_enumeration");
    }

    ISTable* itemLinkedTableP = refBlock.GetTablePtr("item_linked");

    ISTable* itemTypeListTableP = refBlock.GetTablePtr("item_type_list");
  
    ISTable* itemAliasesTableP = refBlock.GetTablePtr("item_aliases");

    ISTable* itemDefaultTableP = refBlock.GetTablePtr("item_default");

    CifParentChild parentChild(refBlock);

    vector<string> catNames;
    block.GetTableNames(catNames);

    for (unsigned int catI = 0; catI < catNames.size(); ++catI)
    {
#ifdef JW_DEBUG
        cout << endl;
        cout << "--------------------------------------------------------------" << endl;
        cout << "Checking category: " << catNames[catI] << " " << catI <<
          " of " << catNames.size() << endl;
#endif
#ifdef VLAD_PERF
        cout << "Checking category: " << catNames[catI] << "; " << catI <<
          " of " << catNames.size() << endl;
#endif
        ISTable* catTableP = block.GetTablePtr(catNames[catI]);

#ifdef VLAD_PERF
        log << "Category: \"" << catNames[catI] << "\"; Index: " << catI
          << " of " << catNames.size() << endl;
#endif
        bool skipCatCheck = false;

        // Check if the table is defined in the reference file
        if (!IsCatDefinedInRef(catNames[catI], *refCatTableP))
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", category \"" << catNames[catI] << "\" is not " <<
              "defined in the reference file" << endl;
            skipCatCheck = true;
        }

        const vector<string>& colNames = catTableP->GetColumnNames();

        for (unsigned int itemI = 0; itemI < colNames.size(); ++itemI)
        {
            if (!IsItemDefinedInRef(catNames[catI], colNames[itemI],
              *itemTableP))
            {
                log << "ERROR - In block \"" << block.GetName() <<
                  "\", in category \"" << catNames[catI] << "\", item \"" <<
                  colNames[itemI] << "\" is not defined in the " <<
                  "reference file" << endl;
            }
        }

        if (skipCatCheck)
            continue;

        vector<string> keyAttributes;

        GetKeyAttributes(keyAttributes, catTableP->GetName(), *catKeyTableP);

        CheckKeyItems(block.GetName(), *catTableP, keyAttributes, *itemTableP,
          itemDefaultTableP, log);

        CheckMandatoryItems(block.GetName(), *catTableP, *itemTableP,
          keyAttributes, log);

#ifdef JW_DEBUG

        cout << "Checking parent/child for category: " << catNames[catI] << endl;
#endif
        ret = parentChild.CheckParentChild(block, *catTableP, log);

#ifdef VLAD_PERF
        cout << "Finished parent/child for category: " << catNames[catI] << endl;
#endif

#ifdef JW_DEBUG
        cout << "Finished parent/child for category: " << catNames[catI] << endl;
#endif
        for (unsigned int colI = 0; colI < colNames.size(); ++colI)
        {
#ifdef VLAD_PERF
            cout << "Processing attribute: " << colI + 1 << " of " <<
              colNames.size() << endl;
#endif
            ret = CheckRegExpRangeEnum(block, *catTableP, colNames[colI],
               *itemTypeTableP, *itemTypeListTableP, itemRangeTableP,
               itemEnumTableP, *itemLinkedTableP, itemAliasesTableP, log);
        }

        // Check for null CIF string rows.
        vector<unsigned int> nullRowsIndices;
        FindCifNullRows(nullRowsIndices, *catTableP);
        if (!nullRowsIndices.empty())
        {
            for (unsigned int rowI = 0; rowI < nullRowsIndices.size(); ++rowI)
            {
                log << "WARNING - In block \"" << block.GetName() <<
                  "\", in category \"" << catNames[catI] << "\", row " <<
                  nullRowsIndices[rowI] + 1 << " contains all null values." <<
                  endl;
            }
        }
    }

    return(ret);
}

void CifFile::FindCifNullRows(vector<unsigned int>& nullRowsIndices,
  const ISTable& isTable)
{
    nullRowsIndices.clear();

    const vector<string>& colNames = isTable.GetColumnNames();

    for (unsigned int rowI = 0; rowI < isTable.GetNumRows(); ++rowI)
    {
        bool nullRow = true;
        for (unsigned int colI = 0; colI < colNames.size(); ++colI)
        {
            if (!CifString::IsUnknownValue(isTable(rowI, colNames[colI])))
            {
                nullRow = false;
                break;
            }
        }

        if (nullRow)
        {
            nullRowsIndices.push_back(rowI);
        }
    }
}


// This method converts the string with embedded escape sequences into
// the real string. For example, a two character combination "\t" is
// converted into one character representing horizontal tab.
void CifFile::ConvertEscapedString(const string& inString,
  string& outString)
{

    outString.clear();

    outString.reserve(inString.size());

    // Go through the characters of the input string
    for (unsigned int posI = 0; posI < inString.size(); )
    {
        if (inString[posI] == '\\')
        {
            if (inString[posI + 1] == '\\')
            {
                  outString.push_back('\\');
                  posI += 2;
            }
            else if (inString[posI + 1] == 't')
            {
                  outString.push_back('\t');
                  posI += 2;
            }
            else if (inString[posI + 1] == 'n')
            {
                  outString.push_back('\n');
                  posI += 2;
            }
            else
            {
                  outString.push_back(inString[posI]);
                  posI++;
            }
        } // Current character is '\'
        else
        {
            outString.push_back(inString[posI]);
            posI++;
        } // Curent character is not '\'
    } // for (characters in the input string)

} // CifFile::ConvertEscapedString()


void CifFile::GetAttributeValueIf(string& attribVal,
  const string& blockId, const string& category, const string& attributeA,
  const string& attributeB, const string& valB)
{

    attribVal.clear();

    Block& block = GetBlock(blockId);

    ISTable* t = block.GetTablePtr(category);

    if (t == NULL)
    {
        // VLAD - EXCEPTION - Throw exception here and test all the code
        return;
    }

    unsigned int nRow = t->GetNumRows();

    if ((nRow > 0) && (t->IsColumnPresent(attributeA)))
    {
        vector<string> qCol, qTar;

        qCol.push_back(attributeB);
        qTar.push_back(valB);

        unsigned int findIndex = t->FindFirst(qTar, qCol);

        if (findIndex != nRow)
        {
            attribVal = (*t)(findIndex, attributeA);
        }
    }

}


void CifFile::GetAttributeValuesIf(vector<string>& strings,
  const string& blockId, const string& category, const string& attributeA,
  const string& attributeB, const string& valB)
{

    strings.clear();

    Block& block = GetBlock(blockId);

    ISTable* t = block.GetTablePtr(category);

    if (t == NULL)
    {
        // VLAD - EXCEPTION - Throw exception here and test all the code
        return;
    }

    vector<string> qCol, qTar;
    qCol.push_back(attributeB);
    qTar.push_back(valB);

    vector<unsigned int> iResult;
    t->Search(iResult, qTar, qCol);
    if (!iResult.empty())
    {
        t->GetColumn(strings, attributeA, iResult);
    }

}


void CifFile::SetAttributeValueIf(const string& blockId,
  const string& category, const string& attributeA, const string& valA, 
  const string& attributeB, const string& valB, const bool create)
{

    if (blockId.empty() || category.empty() || valA.empty() ||
      attributeA.empty() || attributeB.empty() || valB.empty())
    {
        return;
    }

    Block& block = GetBlock(blockId);

    ISTable* t = block.GetTablePtr(category);

    if (t == NULL)
    {
        // VLAD - EXCEPTION - Throw exception here and test all the code
        return;
    }

    vector<string> qCol, qTar;
    qTar.push_back(valB);
    qCol.push_back(attributeB);

    vector<unsigned int> iResult;

    t->Search(iResult, qTar, qCol);

    if (!iResult.empty())
    {
        for (unsigned int i = 0; i < iResult.size(); ++i)
        {
            t->UpdateCell(iResult[i], attributeA, valA);
        }
    }
    else
    {
        if (!create)
        {
            return;
        }

        t->AddRow();
        unsigned int irow = t->GetNumRows() - 1;

        t->UpdateCell(irow, attributeA, valA);
        t->UpdateCell(irow, attributeB, valB);
    }

    block.WriteTable(t);

}


void CifFile::SetAttributeValue(const string& blockId,
  const string& category, const string& attribute, const string& value,
  const bool create)
{

    if (blockId.empty() || category.empty() || attribute.empty() ||
      value.empty())
    {
        return;
    }

    Block& block = GetBlock(blockId);

    ISTable* t = block.GetTablePtr(category);

    if (!create && (t == NULL))
    {
        return;
    }

    unsigned int numRows = 0;

    if (t == NULL)
    {
        t = new ISTable(category);
        t->AddColumn(attribute);
        t->AddRow();
    }

    numRows = t->GetNumRows();

    if (numRows == 0)
    {
        return;
    }

    if (t->IsColumnPresent(attribute))
    {
        for (unsigned int rowI = 0; rowI < numRows; ++rowI)
        {
            t->UpdateCell(rowI, attribute, value);
        }

        block.WriteTable(t);
    }
    // else VLAD - EXCEPTION - Throw exception here and test all the code
}

#ifdef VLAD_TO_CIF_FILE_NOT_USED

void del_attribute_value_where(CifFile *fobj, const char *blockId, const char *category, 
const char *attributeB, const char *valB)
{
  ISTable *t=NULL;

  vector<unsigned int>  iResult;
  int i, nRow;
  vector<string> qCol, qTar;

  if ( !fobj  || !blockId || !category || !attributeB || !valB) return;

  Block& block = fobj->GetBlock(blockId);
  t=block.GetTablePtr(category);

  if (t != NULL) {
    nRow = t->GetNumRows();
    ndb_log_message_text(NDB_MSG_DEBUG,"Got %d rows in %s",nRow,category);
    //t->GetColumnIndex(attributeB);

    if (nRow > 0) {
      qTar.clear();  qTar.push_back(valB);
      qCol.clear();  qCol.push_back(attributeB);
      t->Search(iResult, qTar, qCol);

      if (!iResult.empty()) {
	for (i=0; i < (int) iResult.size(); i++) {
	  t->DeleteRow(iResult[i]);
	}
	block.WriteTable(t);
      }
    }
  }

}
#endif // VLAD_TO_CIF_FILE_NOT_USED not defined


void CifFile::GetAttributeValue(string& attribVal, const string& blockId,
  const string& category, const string& attribute)
{
  attribVal.clear();

  Block& block = GetBlock(blockId);
  if (block.IsTablePresent(category))
  {
    ISTable* t = block.GetTablePtr(category);
    unsigned int nRow = t->GetNumRows();
    if ((nRow != 0) && (t->IsColumnPresent(attribute)))
    {
        attribVal = (*t)(0, attribute);

#ifdef VLAD_DELETED
	t->GetColumn(col, attribute);
	if (!col.empty() && (col[0].size() > 0))
        {
	  string = new char[col[0].size() + 1];
	  strcpy(string,col[0].c_str());
	}
#endif
    }
  }

}


void CifFile::GetAttributeValues(vector<string>& strings,
   const string& blockId, const string& category, const string& attribute)
{
  strings.clear();

  Block& block = GetBlock(blockId);
  if (block.IsTablePresent(category))
  {
    ISTable* t = block.GetTablePtr(category);
    unsigned int nRow = t->GetNumRows();

    if ((nRow != 0) && (t->IsColumnPresent(attribute)))
    {
	t->GetColumn(strings, attribute);
        // VLAD - WIERD LOGIC - IT SHOULD BE LEFT TO THE APPLICATION !!!!
	if (strings.empty() || strings[0].empty())
        {
            strings.clear();
	}
    }
  }

}


void CifFile::SetAttributeValues(const string& blockId,
  const string& category, const string& attribute,
  const vector<string>& values)
{

  if (blockId.empty() || category.empty() || attribute.empty() ||
    values.empty())

  {
      return;
  }

  Block& block = GetBlock(blockId);
  ISTable* t = block.GetTablePtr(category);
  
  if (t != NULL)
  {
    unsigned int nRow = t->GetNumRows();
    if ((nRow != 0) && (t->IsColumnPresent(attribute)))
    {
	t->FillColumn(attribute, values);
	block.WriteTable(t);
    }
  }
}


bool CifFile::IsAttributeValueDefined(const string& blockId,
  const string& category, const string& attribute)
{

  if (blockId.empty() || category.empty() || attribute.empty())
      return(false);

  string attribVal;

  GetAttributeValue(attribVal, blockId, category, attribute);

  if (CifString::IsEmptyValue(attribVal))
      return (false);
  else
      return (true);


#ifdef VLAD_DELETED
  ISTable *t=NULL;
 
  vector<string> col;
  int nRow;

  if ((blockId == NULL) || (category == NULL) || (attribute == NULL) ) return(0);
  
  Block& block = fobj->GetBlock(blockId);
  if ((t=block.GetTablePtr(category)) != NULL) {
    nRow = t->GetNumRows();
    ndb_log_message_text(NDB_MSG_DEBUG,"Found %s with %d rows",category,nRow);
    if (nRow > 0) {
      //iCol = -1;
      //ndb_log_message_text(NDB_MSG_DEBUG,"Found %s at index %d",attribute,iCol);
      if (t->IsColumnPresent(attribute)) {
        //iCol = t->GetColumnIndex(attribute);
	t->GetColumn(col, attribute);
	if (col.empty() || CifString::IsEmptyValue(col[0]))
        {
	  return(0);
	} else 
	  return(1);
      }
    } else 
      return(0);
    
  } else 
    return(0);

  return(0);
#endif

}


void CifFile::SetAttributeValueIfNull(const string& blockId, const string& category,
  const string& attribute, const string& value)
{

  if (blockId.empty() || category.empty() || attribute.empty() ||
    value.empty())
  {
      return;
  }

  Block& block = GetBlock(blockId);
  ISTable* t = block.GetTablePtr(category);
  if (t != NULL)
  {
    unsigned int nRow = t->GetNumRows();
    if ((nRow != 0) && (t->IsColumnPresent(attribute)))
    {
        vector<string> col;
	t->GetColumn(col, attribute);
	if (col.empty() || CifString::IsEmptyValue(col[0]))
        {
          vector<string> uv(nRow, value);
	  t->FillColumn(attribute, uv);
	  block.WriteTable(t);
	}
    }
  }
}


bool CifFile::IsCatDefinedInRef(const string& catName, ISTable& catTable)
{

    /*
    ** For a category, method looks into dictionary ("category" table)
    ** to find out whether it exists or not.
    */

    vector<string> keyTarget;
    keyTarget.push_back(catName);

    vector<string> keyList;
    keyList.push_back("id");

    unsigned int catIndex = catTable.FindFirst(keyTarget, keyList);
    if (catIndex == catTable.GetNumRows())
    {
        return(false);
    }

    return(true);
   
}


bool CifFile::IsItemDefinedInRef(const string& catName, const string& itemName,
  ISTable& refItemTable)
{

    /*
    ** For an item, method looks into dictionary ("item" table)
    ** to find out whether the item is defined in the category or not.
    */

    bool found = false;

    vector<string> refItemList;
    refItemList.push_back("category_id");
    refItemList.push_back("name");

    vector<string> refItemTarget;
    refItemTarget.push_back(catName);

    string cifItem;
    CifString::MakeCifItem(cifItem, catName, itemName);

    refItemTarget.push_back(cifItem);

    unsigned int index = refItemTable.FindFirst(refItemTarget, refItemList);

    if (index != refItemTable.GetNumRows())
    {
        found = true;
    }

    return(found);

}


void CifFile::GetKeyAttributes(vector<string>& keyAttributes,
  const string& catTableName, ISTable& catKeyTable)
{
    /*
    ** For a category, method looks into dictionary ("category_key" table)
    ** to find out what its key items are.
    */

    keyAttributes.clear();

    vector<string> keyTarget;
    keyTarget.push_back(catTableName);

    vector<string> keyList;
    keyList.push_back("id");

    vector<unsigned int> OutList;
    catKeyTable.Search(OutList, keyTarget, keyList);

    for (unsigned int i = 0; i < OutList.size(); ++i)
    {
        string keyAttributeName;

        CifString::GetItemFromCifItem(keyAttributeName,
          catKeyTable(OutList[i], "name"));

        keyAttributes.push_back(keyAttributeName);
    }
}


void CifFile::CheckKeyItems(const string& blockName, ISTable& catTable,
  const vector<string>& keyItems, ISTable& itemTable,
  ISTable* itemDefaultTableP, ostringstream& log)
{
    /*
    ** For a category, method checks for existence of key
    ** items and checks if there are duplicate key values.
    */

    if (keyItems.empty())
    {
        log << "ERROR - In block \"" << blockName << "\", no key items " <<
          "found in category \"" << catTable.GetName() << "\"" << endl;
        return;
    }

    bool keyNotFound = false;

    for (unsigned int k = 0; k < keyItems.size(); k++)
    {
        if (!catTable.IsColumnPresent(keyItems[k]))
        {
            // Missing implicit key values are allowed
            if (IsImplicitNatureKey(catTable.GetName(), keyItems[k],
              itemTable))
            {
                continue;
            }

            log << "ERROR - In block \"" << blockName << "\", key item \"" <<
              keyItems[k] << "\" not found in category \"" <<
              catTable.GetName() << "\"" << endl;

            keyNotFound = true;
        }
    }

    if (keyNotFound)
    {
        return;
    }

    if (itemDefaultTableP == NULL)
    {
        return;
    }

    /*
    ** If any of the keys are implicit in nature and no values are present
    ** for them, a copy of a category table has to be created, missing values
    ** filled from default item values and that new category table be
    ** validated.
    */

    ISTable* valCatTableP = &catTable;

    vector<string> implKeyItems;
    GetImplNatureKeysWithMissingValues(implKeyItems, catTable.GetName(),
      catTable, itemTable);

    if (!implKeyItems.empty())
    {
        // Make a copy of category table that is to be fixed
        valCatTableP = new ISTable(catTable);

        FixMissingValuesOfImplNatureKeys(*valCatTableP, implKeyItems,
          *itemDefaultTableP, log);
    }

    // Check the values of key items
    CheckKeyValues(keyItems, *valCatTableP, log);

    vector<pair<unsigned int, unsigned int> > duplRows;
    (*valCatTableP).FindDuplicateRows(duplRows, keyItems, true);
    if (!duplRows.empty())
    {
        for (unsigned int rowI = 0; rowI < duplRows.size(); rowI++)
        {
            bool report = true;
            for (unsigned int keyI = 0; keyI < keyItems.size(); ++keyI)
            {
                const string& cell = (*valCatTableP)(duplRows[rowI].first,
                  keyItems[keyI]);
             
                if (CifString::IsEmptyValue(cell))
                {
                    report = false;
                    break;
                }
            }

            if (!report)
                continue;

            log << "ERROR - In block \"" << blockName << "\", in " <<
              "cateogory \"" << (*valCatTableP).GetName() <<
              "\", values for key item(s):" << endl;

            for (unsigned int keyI = 0; keyI < keyItems.size(); ++keyI)
            {
                log << "  \"" << keyItems[keyI] << "\"," << endl; 
            }

            log << "  in row #" << duplRows[rowI].first + 1 <<
              " are repeated in row #" << duplRows[rowI].second + 1 << endl;
        }
    }
}


void CifFile::CheckMandatoryItems(const string& blockName, ISTable& catTable,
  ISTable& refItemTable, const vector<string>& keyItems, ostringstream& log)
{

    /*
    ** For a category, method looks into dictionary ("item" table)
    ** to find out what its mandatory items are. Then it checks for existence
    ** of those items in the category.
    */

    vector<string> refItemList;
    refItemList.push_back("category_id");
    refItemList.push_back("mandatory_code");

    vector<string> refItemTarget;
    refItemTarget.push_back(catTable.GetName());
    refItemTarget.push_back("yes");

    vector<unsigned int> OutList;
    refItemTable.Search(OutList, refItemTarget, refItemList);
    if (OutList.empty())
    {
        refItemTarget[1] = "implicit";
 
        refItemTable.Search(OutList, refItemTarget, refItemList);
        if (OutList.empty())
        {
            refItemTarget[1] = "implicit-ordinal";
 
            refItemTable.Search(OutList, refItemTarget, refItemList);
            if (OutList.empty())
            {
                log << "ERROR - In block \"" << blockName <<
                  "\", no mandatory items found in category \"" <<
                  catTable.GetName() << "\"" << endl;
            }
        }

    }

    for (unsigned int k = 0; k < OutList.size(); ++k)
    {
        string cell = refItemTable(OutList[k], "name");
        string itemName;
        CifString::GetItemFromCifItem(itemName, cell);
        if (!catTable.IsColumnPresent(itemName))
        {
            log << "ERROR - In block \"" << blockName <<
              "\", mandatory item \"" << itemName <<
              "\" is not in category \"" << catTable.GetName() <<
              "\"" << endl;
            continue;
        }

        // Values for mandatory items must not be null.
        for (unsigned int rowI = 0; rowI < catTable.GetNumRows(); ++rowI)
        {
            if (catTable(rowI, itemName) == CifString::UnknownValue)
            {
#ifdef VLAD_OLD
                log << "ERROR - In block \"" << blockName <<
                  "\", item \"" << cell <<
                  "\", mandatory item has invalid value \"" <<
                  catTable(rowI, itemName) << "\" in row# " << rowI << endl;
#endif
                log << "ERROR - In block \"" << blockName <<
                  "\", mandatory item \"" << cell <<
                  "\" has invalid value \"" << catTable(rowI, itemName) <<
                  "\"";

                bool first = true;
                for (unsigned int i = 0; i < keyItems.size(); ++i)
                {
                    if (!catTable.IsColumnPresent(keyItems[i]))
                    {
                        continue;
                    }

                    if (keyItems[i] != itemName)
                    {
                        if (first)
                            log << " in row: ";
                        else
                            log << ", ";

                        string cifItem;
                        CifString::MakeCifItem(cifItem, catTable.GetName(),
                          keyItems[i]);
                        log << "\"" << cifItem <<
                          "\" == \"" << catTable(rowI, keyItems[i]) <<
                          "\"";
                        first = false;
                    }
                }
                log << endl;
            }
        }
    }

}


int CifFile::CheckRegExpRangeEnum(Block& block, ISTable& catTable,
  const string& attribName, ISTable& itemTypeTable, 
  ISTable& itemTypeListTable, ISTable* itemRangeTableP,
  ISTable* itemEnumTableP, ISTable& parChildTable,
  ISTable* itemAliasesP, ostringstream& log)
{

    /************** REGEX **************/
    /* 
    ** For an item, a check is being made whether that item is in _item_type
    ** table. If it is not, method ends. If it is, its "code" is taken
    ** from _item_type table and regular expression (from "construct" item
    ** in item_type_list category) table for that "code" is taken and verified
    ** if it is a proper regular expression. Then every value for that item
    ** (that is of regular expression kind) is checked against that regular
    ** expression.
    */

    /************** RANGE **************/
    /* 
    ** Cell type is first determinted and then floating range check or
    ** integer range check is performed. For each found range, it is
    ** checked if cell value is withing that range. If one of the limits
    ** is set to unknown or any, the range against that limit is not checked.
    */

    /************** ENUMERATION **************/
    /* 
    ** Cell type is first determinted and then floating enumeration check or
    ** integer enumeration check is performed. Cell value is checked against
    ** all enumerated values. If value is not one of the enumerated values,
    ** check fails.
    */


    /* compile a regex */
    string patternString;
    regex_t preg;
#define NS 10
    regmatch_t pmatch[NS];
    int icmp = 0;
    int ret = 1;
    vector<string> maxlist;
    vector<string> minlist;
    int errCodeRange = -1;
    vector<string> enumlist;
    int errCodeEnumeration = -1;

    int errCodeRegex = -1;

    string cifItemName;
    CifString::MakeCifItem(cifItemName, catTable.GetName(), attribName);

    vector<string> target;
    target.push_back(cifItemName);

    vector<string> nameList;
    nameList.push_back("name");

    string typeCode;
    GetItemTypeCode(typeCode, cifItemName, itemTypeTable);

    if (typeCode.empty())
    {
        return(1);
    }

    /************** REGEX step 1**************/
    errCodeRegex = 0;

    vector<string> ItemTypeLTarget;
    ItemTypeLTarget.push_back(typeCode);

    vector<string> ItemTypeLList;
    ItemTypeLList.push_back("code");

    unsigned int iOut = itemTypeListTable.FindFirst(ItemTypeLTarget,
      ItemTypeLList);

    string primCode = itemTypeListTable(iOut, "primitive_code");

    ConvertEscapedString(itemTypeListTable(iOut, "construct"),
      patternString); 

    icmp = rcsb_regcomp(&preg, patternString.c_str(), REG_EXTENDED);

    rcsb_regfree(&preg);


    /************** RANGE step 1**************/
    vector<unsigned int> OutList;
    if (itemRangeTableP != NULL)
    {
        itemRangeTableP->Search(OutList, target, nameList);
        if (!OutList.empty())
        {
            errCodeRange = 0;
            for (unsigned int l = 0; l < OutList.size(); ++l)
            {
                maxlist.push_back((*itemRangeTableP)(OutList[l], "maximum"));
                minlist.push_back((*itemRangeTableP)(OutList[l], "minimum"));
            }
        }
    }

    /************ ENUMERATION step 1************/
    if (itemEnumTableP != NULL)
    {
        itemEnumTableP->Search(OutList, target, nameList);
        if (!OutList.empty())
        {
            errCodeEnumeration = 0;
            for (unsigned int l = 0; l < OutList.size(); ++l)
            {
                enumlist.push_back((*itemEnumTableP)(OutList[l], "value"));
            }
        }
    }

    string itemName;
    CifString::GetItemFromCifItem(itemName, cifItemName);
    for (unsigned int l = 0; l < catTable.GetNumRows(); l++)
    {
#ifdef VLAD_PERF
                cout << "  Processing regex row: " << l << " of " <<
                  catTable.GetNumRows() << endl;
#endif
        // loop for all rows
        const string& cell = catTable(l, itemName);

        /*************** REGEX step 2**************/
        if (!CifString::IsEmptyValue(cell))
        {
            if ((errCodeRegex >= 0) && (icmp == 0))
            {
                rcsb_regcomp(&preg, patternString.c_str(), REG_EXTENDED);
                ret = rcsb_regexec(&preg, cell.c_str(), NS, pmatch, 0);
                rcsb_regfree(&preg);
                if (ret != 0)
                {
                    log << "ERROR - In block \"" << block.GetName() <<
                       "\", data type pattern of value \"" << cell <<
                       "\" for \"" << cifItemName <<
                       "\" does not match, in row #" << l + 1 << endl;
                    ret = 0;
                }
                else
                {
                    int len = pmatch[0].rm_eo - pmatch[0].rm_so;
                    if ((len != (int)(cell.size())) && (typeCode != "text") && (typeCode != "symop"))
                    {
                        log << "ERROR - In block \"" << block.GetName() <<
                          "\", in item \"" << cifItemName << "\"";

                        if (itemAliasesP != NULL)
                        {
                            // Find out if item has aliases.
                            vector<string> itemAliasesTarget;
                            itemAliasesTarget.push_back(cifItemName);

                            vector<string> itemAliasesColumns;
                            itemAliasesColumns.push_back("name");

                            unsigned int where =
                              itemAliasesP->FindFirst(itemAliasesTarget,
                              itemAliasesColumns);

                            if (where != itemAliasesP->GetNumRows())
                            {
                                // Found that item has an alias
                                log << " (its alias is: \"" <<
                                  (*itemAliasesP)(where, "alias_name") << "\")";
                            }
                        }

                        log << ", in row #" << l + 1 << ", value \"" <<
                          cell << "\" does not "\
                          "correspond to item's type. The expected"\
                          " type is \"" << typeCode << "\"." << endl;
                        ret = 0;
                    }
                }
            }

            if (errCodeRange >= 0)
            {
                ret = CheckCellRange(cell, typeCode, minlist, maxlist);
                if (ret == 0)
                {
                    log << "ERROR - In block \"" << block.GetName() <<
                      "\", value \"" << cell << "\" for \"" << cifItemName <<
                      "\" is out of range, in row #" << l + 1 << endl;
                }
            }

            if (errCodeEnumeration >= 0)
            {
                ret = CheckCellEnum(cell, typeCode, primCode, enumlist);
                if (ret == 0)
                {
                    log << "ERROR - In block \"" << block.GetName() <<
                      "\", for \"" << cifItemName << "\", value \"" << cell <<
                      "\" is not defined in the enumeration list. "\
                      "Acceptable values are: ";
                    for (unsigned int enumI = 0; enumI < enumlist.size();
                      ++enumI)
                    {
                        log << "\"" << enumlist[enumI] << "\"";
                        if (enumI != (enumlist.size() - 1))
                            log << ", ";
                    }
                    log << endl;
                }
            }
        }
    }

    return(ret);
}


int CifFile::CheckCellRange(const string& cell, const string& typeCode, 
  const vector<string>& minlist, const vector<string>& maxlist)
{
    try
    {
        int matched = 1;

        if (typeCode == "float")
        {
            matched = CheckCellFloatRange(cell, minlist, maxlist);
        }
        else
        {
            if (typeCode == "int")
            {
                matched = CheckCellIntRange(cell, minlist, maxlist);
            }
        }

        return(matched);
    }
    catch (exception& exc)
    {
        return (1);
    }
}


int CifFile::CheckCellEnum(const string& cell, const string& typeCode,
  const string& primCode, const vector<string>& enumlist)
{

    int matched = 1;

    if (typeCode == "float")
    {
        matched = CheckCellFloatEnum(cell, enumlist);
    }
    else
    {
        if (typeCode == "int")
        {
            matched = CheckCellIntEnum(cell, enumlist);
        }
        else
        {
            matched = CheckCellOtherEnum(cell, primCode, enumlist);
        }
    }

    return(matched);

}


int CifFile::CheckCellFloatRange(const string& cell,
  const vector<string>& minlist, const vector<string>& maxlist)
{

    int matched = 0;

    double dmin = 0;
    double dmax = 0;

    double dval = String::StringToDouble(cell);

    unsigned int m = 0;
    while (m < minlist.size() && !matched)
    {

        bool checkMin = false;
        if (!CifString::IsEmptyValue(minlist[m]))
        {
            checkMin = true;
            dmin = String::StringToDouble(minlist[m]);
        }

        bool checkMax = false;
        if (!CifString::IsEmptyValue(maxlist[m]))
        {
            checkMax = true;
            dmax = String::StringToDouble(maxlist[m]);
        }

        if ((!checkMin) && (checkMax))
        {
            if (dval < dmax)
                matched = 1;
        }

        if ((!checkMax) && (checkMin))
        {
            if (dval > dmin)
                matched = 1;
        }

        if (checkMax && checkMin)
        {
            if (dmin == dmax)
            {
                if (dval == dmax)
                    matched = 1;
            }
            else
            {
                if ((dval > dmin) && (dval < dmax))
                    matched = 1;
            }
        }

        m++;
    }

    return(matched);

}


int CifFile::CheckCellIntRange(const string& cell,
  const vector<string>& minlist, const vector<string>& maxlist)
{

    int matched = 0;

    int imin = 0;
    int imax = 0;

    int ival = String::StringToInt(cell);

    unsigned int m = 0;
    while (m < minlist.size() && !matched)
    {

        bool checkMin = false;
        if (!CifString::IsEmptyValue(minlist[m]))
        {
            checkMin = true;
            imin = String::StringToInt(minlist[m]);
        }

        bool checkMax = false;
        if (!CifString::IsEmptyValue(maxlist[m]))
        {
            checkMax = true;
            imax = String::StringToInt(maxlist[m]);
        }

        if ((!checkMin) && (checkMax))
        {
            if (ival < imax)
                matched = 1;
        }

        if ((!checkMax) && (checkMin))
        {
            if (ival > imin)
                matched = 1;
        }

        if (checkMax && checkMin)
        {
            if (imin == imax)
            {
                if (ival == imax)
                    matched = 1;
            }
            else
            {
                if ((ival > imin) && (ival < imax))
                    matched = 1;
            }
        }

        m++;
    }

    return(matched);

}


int CifFile::CheckCellFloatEnum(const string& cell,
  const vector<string>& enumlist)
{

    double de;
    double dval;

    int matched = 0;

    dval = String::StringToDouble(cell);

    unsigned int m = 0;
    while (m < enumlist.size() && !matched)
    {
        de = String::StringToDouble(enumlist[m]);
        if (dval == de)
            matched = 1;
        m++;
    }

    return(matched);

}


int CifFile::CheckCellIntEnum(const string& cell,
  const vector<string>& enumlist)
{

    int matched = 0;

    int ival = String::StringToInt(cell);

    unsigned int m = 0;
    while (m < enumlist.size() && !matched)
    {
        int ie = String::StringToInt(enumlist[m]);
        if (ival == ie)
            matched = 1;
        m++;
    }

    return(matched);

}


int CifFile::CheckCellOtherEnum(const string& cell, const string& primCode,
  const vector<string>& enumlist)
{

    int matched = 0;

    unsigned int m = 0;
    while (m < enumlist.size() && !matched)
    {
        if (primCode == "uchar")
        {
            if (_enumCaseSense == true)
            {
                if (String::IsEqual(cell, enumlist[m], Char::eCASE_SENSITIVE))
                    matched = 1;
            }
            else
            {
                if (String::IsCiEqual(cell, enumlist[m]))
                    matched = 1;
            }
        }
        else
        {
            if (cell == enumlist[m])
                matched = 1;
        }
        m++;
    }

    return(matched);

}


void CifFile::CheckAndRectifyItemTypeCode(Block& block, ostringstream& log)
{
    ISTable* itemTypeTableP = block.GetTablePtr("item_type");

    ISTable* itemTableP = block.GetTablePtr("item");

    if ((itemTypeTableP == NULL) || (itemTableP == NULL))
        return;

    // VLAD: Begin: Refactor this with GetCatNames()
    ISTable* categoryTableP = block.GetTablePtr("category");
    if (categoryTableP == NULL)
    {
        log << "CRITICAL: CANNOT FIND DDL CATEGORY: " << "category" << endl;
        return;
    }

    vector<string> categories;
    categoryTableP->GetColumn(categories, "id");
    // VLAD: End: Refactor this with GetCatNames()

    CifParentChild cifParentChild(block);

    for (unsigned int catI = 0; catI < categories.size(); ++catI)
    {
        const string& catName = categories[catI];

        // VLAD: Begin: Refactor this with GetCatItemsNames()
        ISTable* itemCatP = block.GetTablePtr("item");
        if (itemCatP == NULL)
        {
            log << "CRITICAL: CANNOT FIND DDL CATEGORY: " << "item" << endl;
            return;
        }

        vector<string> itemsNames;
        // Get all items of a category
        for (unsigned int itemI = 0; itemI < itemCatP->GetNumRows(); ++itemI)
        {
            const string& itemName = (*itemCatP)(itemI, "name");

            string itemCatName;
            CifString::GetCategoryFromCifItem(itemCatName, itemName);

            if (itemCatName == catName)
            {
                itemsNames.push_back(itemName);
            }
        }
        // VLAD: end: Refactor this with GetCatItemsNames()

        for (unsigned int itemI = 0; itemI < itemsNames.size(); ++itemI)
        {
            string itemTypeCode;
            RectifyItemTypeCode(itemTypeCode, log, block, cifParentChild,
              itemsNames[itemI]);
        } // for (all items)
    } // for (all categories)
}


void CifFile::RectifyItemTypeCode(string& retItemTypeCode,
  std::ostringstream& log, Block& block, CifParentChild& cifParentChild,
  const string& cifItemName)
{
    retItemTypeCode.clear();

    string itemTypeCode;

    vector<string> target;
    target.push_back(cifItemName);

    vector<string> nameList;
    nameList.push_back("name");

    ISTable* itemTypeTableP = block.GetTablePtr("item_type");

    unsigned int iOut = itemTypeTableP->FindFirst(target, nameList);
    if (iOut != itemTypeTableP->GetNumRows())
    {
        itemTypeCode = (*itemTypeTableP)(iOut, "code");
    }

    bool hasItemTypeCode = !itemTypeCode.empty();

    if (!hasItemTypeCode)
    {
        if (_extraDictChecks)
        {
#ifdef EXTENDED_CHECK_MESSAGE_SUPPRESSED
            log << "ERROR - In block \"" << block.GetName() <<
              "\", \"_item_type.code\" not defined for item \"" <<
              cifItemName << "\" and strict checking options prevent "\
              "deducting its type from its parents" << endl;
#endif

            return;
        }
    }

    string catName;
    CifString::GetCategoryFromCifItem(catName, cifItemName);

    vector<vector<string> > parParKeys;
    vector<vector<string> > comboComboKeys;

    cifParentChild.GetParents(parParKeys, comboComboKeys, catName);

    vector<pair<string, string> > parItemsAndTypes;

    // VLAD - LATER - Change this to fomaly go via comboComboKeys iteration
    for (unsigned int allParI = 0; allParI < parParKeys.size(); ++allParI)
    {
        vector<string>& parKeys = parParKeys[allParI];
        vector<string>& comboKeys = comboComboKeys[allParI];

        for (unsigned int keysI = 0; keysI < comboKeys.size(); ++keysI)
        {
            if (cifItemName != comboKeys[keysI])
            {
                continue;
            }

            if (parKeys[keysI] == cifItemName)
            {
                throw runtime_error("ERROR - Parent and child key have the "\
                  "same value \"" + cifItemName + "\"");
            }

            string parItemTypeCode;
            RectifyItemTypeCode(parItemTypeCode, log, block, cifParentChild,
              parKeys[keysI]);

            parItemsAndTypes.push_back(make_pair(parKeys[keysI],
              parItemTypeCode));
        }
    }

#ifdef VLAD_DEBUG
    log << "DEBUG - cifItem:parItemsAndTypes.size() \"" << cifItemName <<
      ":" << parItemsAndTypes.size() << endl;
#endif

    bool hasParents = !parItemsAndTypes.empty();

    // If no parents at all, set parType to null.
    bool hasParTypeCode = false;

    string parTypeCode;

    if (hasParents)
    {
        set<string> parTypes;
        for (unsigned int pairI = 0; pairI < parItemsAndTypes.size(); ++pairI)
        {
            parTypes.insert(parItemsAndTypes[pairI].second);
        }

        // If any type is empty string, report and return

        if (parTypes.size() > 1)
        {
            log << "ERROR - Item \"" << cifItemName << "\" has parents "\
              "with different types:" << endl;

            for (unsigned int pairI = 0; pairI < parItemsAndTypes.size();
              ++pairI)
            {
                log << "  Parent item \"" << parItemsAndTypes[pairI].first <<
                  "\" has item type code \"" <<
                  parItemsAndTypes[pairI].second << "\"" << endl;
            }
        }
        else
        {
            // Take it from the first parent
            hasParTypeCode = true;
            parTypeCode = parItemsAndTypes[0].second;
        }
    }

    if (!hasParents && !hasItemTypeCode)
    {
        // No type code in either the item or the parents
        log << "ERROR - In block \"" << block.GetName() <<
          "\", item \"" << cifItemName << "\"" << 
          " does not have item type code defined and has no parents." << endl;
        return;
    }

    if ((!hasParTypeCode) && (!hasItemTypeCode))
    {
        log << "ERROR - In block \"" << block.GetName() <<
          "\", item \"" << cifItemName << "\"" <<
          " does not have item type code defined and its item type code "\
          " cannot be deducted from its parents." << endl;
        return;
    }

    if (hasItemTypeCode)
    {
        retItemTypeCode = itemTypeCode;

        if (hasParTypeCode && (itemTypeCode != parTypeCode))
        {
            log << "ERROR - In block \"" << block.GetName() <<
              "\", child item \"" << cifItemName <<
              "\" has item type code \"" << itemTypeCode <<
              "\", while its parent(s) have item type code \"" <<
              parTypeCode << "\"" << endl;
        }
    }
    else
    {
#ifdef VLAD_TRACE
        log << "INFO - For item \"" << cifItemName << "\", inserting parent "\
          "item type code \"" << parTypeCode << "\"" << endl;
#endif

        // Add it to the "item_type" table of the item
        itemTypeTableP->AddRow();
        unsigned int lastRowInd = itemTypeTableP->GetNumRows() - 1;

        itemTypeTableP->UpdateCell(lastRowInd, "name",
          cifItemName);
        itemTypeTableP->UpdateCell(lastRowInd, "code",
          parTypeCode);

        retItemTypeCode = parTypeCode;
    }
}


void CifFile::GetItemTypeCode(string& typeCode, const string& cifItemName, 
  ISTable& itemTypeTable)
{

    typeCode.clear();

    vector<string> target;
    target.push_back(cifItemName);

    vector<string> nameList;
    nameList.push_back("name");

    unsigned int iOut = itemTypeTable.FindFirst(target, nameList);
    if (iOut != itemTypeTable.GetNumRows())
    {
        typeCode = itemTypeTable(iOut, "code");
    }

}

void CifFile::CheckKeyValues(const vector<string>& keysAttribs,
  ISTable& catTable, ostringstream& log)
{
    for (unsigned int rowI = 0; rowI < catTable.GetNumRows(); ++rowI)
    {
        for (unsigned int keyI = 0; keyI < keysAttribs.size(); ++keyI)
        {
            const string& value = catTable(rowI, keysAttribs[keyI]);

            if (CifString::IsEmptyValue(value))
            {
                string keyItem;
                CifString::MakeCifItem(keyItem, catTable.GetName(),
                  keysAttribs[keyI]);

#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
                if ((keyItem == "_atom_sites_alt.id") &&
                  (value == CifString::InapplicableValue))
                {
                    continue;
                }
#endif

                log << "ERROR - Key item \"" << keyItem <<
                  "\" has invalid value \"" << value << "\"";

                if (keysAttribs.size() > 1)
                    log << " in row having:";

                log << endl;

                for (unsigned int kI = 0; kI < keysAttribs.size(); ++kI)
                {
                    if (kI == keyI)
                    {
                        // Do not print this item again.
                        continue;
                    }

                    string tmpKeyItem;
                    CifString::MakeCifItem(tmpKeyItem, catTable.GetName(),
                      keysAttribs[kI]);

                    const string& tmpValue = catTable(rowI, keysAttribs[kI]);

                    log << "  item \"" << tmpKeyItem << "\" with value = \"" <<
                      tmpValue << "\"" << endl;
                }
            }
        }
    }
}


bool CifFile::IsImplicitNatureKey(const string& catName,
  const string& attribName, ISTable& itemTable)
{
    vector<string> refItemList;
    refItemList.push_back("name");
    refItemList.push_back("mandatory_code");

    string cifItem;
    CifString::MakeCifItem(cifItem, catName, attribName);

    vector<string> refItemTarget;
    refItemTarget.push_back(cifItem);
    refItemTarget.push_back("implicit");

    vector<unsigned int> OutList;
    itemTable.Search(OutList, refItemTarget, refItemList);
    if (!OutList.empty())
    {
        return true;
    }

    refItemTarget[1] = "implicit-ordinal";
 
    itemTable.Search(OutList, refItemTarget, refItemList);
    if (!OutList.empty())
    {
        return true;
    }

    return false;
}


void CifFile::GetImplNatureKeysWithMissingValues(vector<string>& implKeyItems,
  const string& catName, ISTable& catTable, ISTable& refItemTable)
{
    // Return vector of keys (as CIF items) having implicit nature for which
    // some of values in the category table are unknown
    implKeyItems.clear();

    vector<string> allImplNatKeyItems;
    GetImplNatureKeys(allImplNatKeyItems, catName, refItemTable);

    for (unsigned int keyInd = 0; keyInd < allImplNatKeyItems.size(); ++keyInd)
    {
        string attribName;
        CifString::GetItemFromCifItem(attribName, allImplNatKeyItems[keyInd]);

        if (AreSomeValuesInColumnEmpty(catTable, attribName))
        {
            implKeyItems.push_back(allImplNatKeyItems[keyInd]);
        }
    }
}


void CifFile::FixMissingValuesOfImplNatureKeys(ISTable& catTable,
  const vector<string>& implNatureKeys, ISTable& refItemDefaultTable,
  ostringstream& log)
{
    for (unsigned int keyInd = 0; keyInd < implNatureKeys.size(); ++keyInd)
    {
        string defValue;
        GetItemDefaultValue(defValue, implNatureKeys[keyInd],
          refItemDefaultTable);

        if (CifString::IsEmptyValue(defValue))
        {
            log << "ERROR - mandatory implicit item \"" <<
              implNatureKeys[keyInd] << "\" has no default value defined" <<
              endl;
            continue;
        }

        string attribName;
        CifString::GetItemFromCifItem(attribName, implNatureKeys[keyInd]);
        if (!catTable.IsColumnPresent(attribName))
        {
            catTable.AddColumn(attribName);
        }

        for (unsigned int rowI = 0; rowI < catTable.GetNumRows(); ++rowI)
        {
            const string& value = catTable(rowI, attribName);

            if (CifString::IsEmptyValue(value))
            {
                catTable.UpdateCell(rowI, attribName, defValue);
            }
        }
    }
}


void CifFile::GetImplNatureKeys(vector<string>& implNatureKeys,
  const string& catName, ISTable& refItemTable)
{
    implNatureKeys.clear();

    vector<string> attribsNames;
    attribsNames.push_back("category_id");
    attribsNames.push_back("mandatory_code");

    vector<string> targetsValues;
    targetsValues.push_back(catName);
    targetsValues.push_back("implicit");

    vector<unsigned int> matchRowInds;
    refItemTable.Search(matchRowInds, targetsValues, attribsNames);

    for (unsigned int matchInd = 0; matchInd < matchRowInds.size(); ++matchInd)
    {
        implNatureKeys.push_back(refItemTable(matchRowInds[matchInd], "name"));
    }
    
    // VLAD - TODO - For implicit-ordinal. IMPLEMENT LATER WHEN STABLE
    /*
    matchRowInds.clear();
    targetsValues[1] = "implicit-ordinal";

    refItemTable.Search(matchRowInds, targetsValues, attribsNames);

    for (matchInd = 0; matchInd < matchRowInds.size(); ++matchInd)
    {
        implNatureKeys.push_back(refItemTable(matchRowInds[mathInd], "name"));
    }
    */
}


bool CifFile::AreSomeValuesInColumnEmpty(ISTable& table, const string& colName)
{
    if (!table.IsColumnPresent(colName))
    {
        return true;
    }

    for (unsigned int rowI = 0; rowI < table.GetNumRows(); ++rowI)
    {
        const string& value = table(rowI, colName);

        if (CifString::IsEmptyValue(value))
        {
            return true;
        }
    }

    return false;
}


void CifFile::GetItemDefaultValue(string& defValue, const string& implNatKey,
  ISTable& refItemDefaultTable)
{
    defValue.clear();

    vector<string> attribsNames;
    attribsNames.push_back("name");

    vector<string> targetsValues;
    targetsValues.push_back(implNatKey);

    unsigned int iOut = refItemDefaultTable.FindFirst(targetsValues,
      attribsNames);
    if (iOut != refItemDefaultTable.GetNumRows())
    {
        defValue = refItemDefaultTable(iOut, "value");
    }
}


