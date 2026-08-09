//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/**
** \file CifCorrector.h
**
** CIF corrector class.
*/


#ifndef CIFCORRECTOR_H
#define CIFCORRECTOR_H


#include <string>
#include <vector>

#include "DataInfo.h"
#include "CifFile.h"


class CifCorrector
{
  public:
    static void MakeOutputCifFileName(std::string& outCifFileName,
      const std::string& inCifFileName);
    static CifFile* CreateConfigFile();

    CifCorrector(CifFile& cifFile, DataInfo& dataInfo, DataInfo& pdbxDataInfo,
      CifFile& configFile, const bool verbose = false);
    ~CifCorrector();

    void Correct();
    void CheckAliases();

    void Write(const std::string& outFileName);

    static void CorrectEnumsSimple(CifFile& cifFile, DataInfo& dataInfo,
      const bool verbose = false);

  private:
    CifFile& _cifFile;

    DataInfo& _dataInfo;
    DataInfo& _pdbxDataInfo;

    CifFile& _configFile;

    bool _verbose;

    ISTable* _configTableP;

    void ValidateConfigTable();

    void RemoveItem(const std::string& item);
    void CorrectUpperCase(const std::string& item);
    void RenameItem(const std::string& item, const std::string& refItem);
    void CorrectValues(const std::string& item, const std::string& itemValue,
      const std::string& refItemValue);
    void CorrectEnums();
    void CorrectNumericList(const std::string& item);
    void CorrectNotApplicableValues();
    void CorrectMissingValues(const std::string& item,
      const std::string& refItem);
    void CorrectLabeling(const std::string& item, 
      const std::string& refItem);
    void CorrectBadSequence(const std::string& item,
      const std::string& refCondItem, const std::string& refCondItemValue);

    static unsigned int FindEnumIndex(const std::string& value,
      const std::vector<std::string>& enums);
    void FixNumericList(std::string& outValue, const std::string& inValue);
    void FixNotApplicable(std::string& outValue, const std::string& inValue);
    void FixBadSequence(std::vector<std::string>& outValues,
      std::vector<std::pair<std::vector<std::string>,
      std::vector<std::string> > >& refMap,
      std::vector<std::vector<std::string> >& confRefValues, Block& block,
      const unsigned int numRows);

    void GetSrcValues(vector<string>& srcValues,
      const vector<string>& srcItems, const vector<string>& prevRefItems,
      const vector<string>& confValues, Block& block,
      const unsigned int rowIndex);

    void GetRefValues(std::vector<std::string>& refValues,
      const std::vector<std::string>& srcItems,
      const std::vector<std::string>& srcValues,
      const std::vector<std::string>& refItems, Block& block);
    void GetTableValues(vector<string>& values, const vector<string>& items,
      Block& block, const unsigned int rowIndex);


    void CreateRefMap(std::vector<std::pair<std::vector<std::string>,
      std::vector<std::string> > >& refMap, const std::string& refCondItem);
    void ExtractRefValues(std::vector<std::vector<std::string> >& refValues,
      const std::string& refCondItemValue);

    void GetValuesFromSeparatedString(vector<string>& values,
      const std::string& valueString, const std::string& sep);
    void NumberToLetter(char& letter, const unsigned int number);
};


#endif
