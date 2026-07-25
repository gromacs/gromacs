#include <string>

#include "RcsbFile.h"
#include "CifFile.h"
#include "DicFile.h"

#include "CifParserBase.h"
#include "DICParserBase.h"

#include "CifDataInfo.h"
#include "CifCorrector.h"

#include "CifParentChild.h"
#include "CifFileUtil.h"


using std::string;


DicFile* GetDictFile(DicFile* ddlFileP, const string& dictFileName,
  const string& dictSdbFileName, const bool verbose, const eFileMode fileMode)
{
    DicFile* dictFileP = NULL;

    if (!dictFileName.empty())
    {
        // If dictionary text file is specified, fileMode is ignored, as
        // it is always writeable.

        dictFileP = ParseDict(dictFileName, ddlFileP, verbose);

        Block& block = dictFileP->GetBlock(dictFileP->GetFirstBlockName());

        CifParentChild cifParentChild(block);

        cifParentChild.WriteGroupTables(block);

        return(dictFileP);
    }

    if (!dictSdbFileName.empty())
    {
        dictFileP = new DicFile(fileMode, dictSdbFileName, verbose);

        return(dictFileP);
    }

    return(dictFileP);
}


void CheckDict(DicFile* dictFileP, DicFile* ddlFileP,
  const string& dictFileName, const bool extraDictChecks)
{
    string relLogFileName;
    RcsbFile::RelativeFileName(relLogFileName, dictFileName);

    relLogFileName += "-diag.log";

    dictFileP->DataChecking(*ddlFileP, relLogFileName, extraDictChecks);
}


void CheckCif(CifFile* cifFileP, DicFile* dictFileP, const string& cifFileName,
  const bool extraCifChecks,
  const std::vector<std::string>& skipBlockNames)
{
    string relLogFileName;
    RcsbFile::RelativeFileName(relLogFileName, cifFileName);
    
    relLogFileName += "-diag.log";

    cifFileP->DataChecking(*dictFileP, relLogFileName, false, extraCifChecks, skipBlockNames);
}


CifFile* ParseCif(const string& fileName, const bool verbose,
  const Char::eCompareType caseSense, const unsigned int maxLineLength,
  const string& nullValue, const string& parseLogFileName)
{
    CifFile* cifFileP = new CifFile(verbose, caseSense, maxLineLength,
      nullValue);

    cifFileP->SetSrcFileName(fileName);

    CifParser cifParser(cifFileP, verbose);

    cifParser.Parse(fileName, cifFileP->_parsingDiags, parseLogFileName);

    return (cifFileP);
}


CifFile* ParseCifString(const string& cifString, const bool verbose,
  const Char::eCompareType caseSense, const unsigned int maxLineLength,
  const string& nullValue)
{
    CifFile* cifFileP = new CifFile(verbose, caseSense, maxLineLength,
      nullValue);

    CifParser cifParser(cifFileP, verbose);

    cifParser.ParseString(cifString, cifFileP->_parsingDiags);

    return (cifFileP);
}


CifFile* ParseCifSimple(const string& fileName,
  const bool verbose, const unsigned int intCaseSense,
  const unsigned int maxLineLength, const string& nullValue,
  const string& parseLogFileName)
{
    return ParseCif(fileName, verbose, (Char::eCompareType)intCaseSense,
      maxLineLength, nullValue, parseLogFileName);
}


CifFile* ParseCifSelective(const string& fileName,
  const CifFileReadDef& readDef, const bool verbose,
  const unsigned int intCaseSense, const unsigned int maxLineLength,
  const string& nullValue, const string& parseLogFileName)
{
    CifFile* cifFileP = new CifFile(verbose, (Char::eCompareType)intCaseSense,
      maxLineLength, nullValue);

    cifFileP->SetSrcFileName(fileName);

    CifParser cifParser(cifFileP, readDef, verbose);

    cifParser.Parse(fileName, cifFileP->_parsingDiags, parseLogFileName);

    return (cifFileP);
}


DicFile* ParseDict(const string& dictFileName, DicFile* inRefFileP,
  const bool verbose)
{
    CifFile* refFileP = NULL;

    DicFile* dictFileP = new DicFile(verbose);

    dictFileP->SetSrcFileName(dictFileName);

    if (inRefFileP == NULL)
    {
        // DDL file parsing
        refFileP = dictFileP->GetRefFile();
    }
    else
    {
        refFileP = inRefFileP;
    }

    DICParser dicParser(dictFileP, refFileP, verbose);

    dicParser.Parse(dictFileName, dictFileP->_parsingDiags);

    if (inRefFileP != NULL)
    {
        // Compres only if dictionary parsing is done
#ifndef VLAD_TMP_DEL_RESOLVE
        // The existence of Compress() indicates problems. Our code base
        // should never invoke this method !!!!
        dictFileP->Compress(refFileP);
#endif
    }
    else
    {
        delete (refFileP);
    }

    return (dictFileP);
}


void DataCorrection(CifFile& cifFile, DicFile& dicRef)
{
    CifDataInfo cifDataInfo(dicRef);

    CifCorrector::CorrectEnumsSimple(cifFile, cifDataInfo);
}

