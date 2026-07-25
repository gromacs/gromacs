//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef CIFFILEUTIL_H
#define CIFFILEUTIL_H


#include <string>

#include "CifFileReadDef.h"
#include "DicFile.h"
#include "CifFile.h"


DicFile* GetDictFile(DicFile* ddlFileP, const std::string& dictFileName,
  const std::string& dictSdbFileName = std::string(), const bool verbose =
  false, const eFileMode fileMode = READ_MODE);
void CheckDict(DicFile* dictFileP, DicFile* ddlFileP,
  const string& dictFileName, const bool extraDictChecks = false);
void CheckCif(CifFile* cifFileP, DicFile* dictFileP,
  const string& cifFileName, const bool extraCifChecks = false,
  const std::vector<std::string>& skipBlockNames = std::vector<std::string>());

DicFile* ParseDict(const std::string& dictFileName, DicFile* ddlFileP = NULL,
  const bool verbose = false);
CifFile* ParseCif(const std::string& fileName, const bool verbose = false,
  const Char::eCompareType caseSense = Char::eCASE_SENSITIVE,
  const unsigned int maxLineLength = CifFile::STD_CIF_LINE_LENGTH,
  const std::string& nullValue = CifString::UnknownValue,
  const std::string& parseLogFileName = std::string());
CifFile* ParseCifString(const std::string& cifString,
  const bool verbose = false,
  const Char::eCompareType caseSense = Char::eCASE_SENSITIVE,
  const unsigned int maxLineLength = CifFile::STD_CIF_LINE_LENGTH,
  const std::string& nullValue = CifString::UnknownValue);
CifFile* ParseCifSimple(const std::string& fileName,
  const bool verbose = false,
  const unsigned int intCaseSense = 0,
  const unsigned int maxLineLength = CifFile::STD_CIF_LINE_LENGTH,
  const std::string& nullValue = CifString::UnknownValue,
  const std::string& parseLogFileName = std::string());
CifFile* ParseCifSelective(const std::string& fileName,
  const CifFileReadDef& readDef, const bool verbose = false,
  const unsigned int intCaseSense = 0,
  const unsigned int maxLineLength = CifFile::STD_CIF_LINE_LENGTH,
  const std::string& nullValue = CifString::UnknownValue,
  const std::string& parseLogFileName = std::string());

/**
**  Corrects a CIF file with respect to the following:
**    - Sets proper casing of the case-insensitive enumerations
**
**  \param[in] dicRef - reference to a dictionary file. The check is
**    done against the first block in the dictionary file.
**
**  \return None
**
**  \pre None
**
**  \post None
**
**  \exception: None
*/
void DataCorrection(CifFile& cifFile, DicFile& dicRef);

#endif
