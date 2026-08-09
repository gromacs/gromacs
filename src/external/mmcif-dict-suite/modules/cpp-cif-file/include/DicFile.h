//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


/*!
** \file DicFile.h
**
** \brief Header file for DicFile class.
*/


#ifndef DICFILE_H
#define DICFILE_H


#include <string>
#include <iostream>

#include "GenString.h"
#include "ISTable.h"
#include "CifFile.h"


/**
**  \class DicFile
**
**  \brief Public class that represents a dictionary file, composed of
**    blocks with tables.
**
**  This class represents a dictionary file. In addition to inherited methods
**  from \e CifFile class, this class provides a method for writing the
**  content of "item_aliases" table to a text file.
*/
class DicFile : public CifFile
{
  public:
    using CifFile::Write;

    /**
    **  Constructs a dictionary file.
    **
    **  \param[in] fileMode - dictionary file mode. Possible values are
    **    read-only, create, update and virtual. Detailed description of
    **    file mode is given in \e TableFile documentation.
    **  \param[in] fileName - relative or absolute name of the file
    **    where object persistency is maintained. If \e fileMode specifies
    **    virtual mode, this parameter is ignored.
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    DicFile(const eFileMode fileMode, const std::string& objFileName,
      const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Constructs a dictionary file in virtual mode.
    **
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    DicFile(const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Destructs a dictionary file, by releasing all consumed resources.
    **
    **  \param: Not applicable
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ~DicFile();

    /**
    **  Writes the content of "item_aliases" table to a text file.
    **
    **  \param[in] fileName - relative or absolute name of the text file
    **    to which the content of "item_aliases" table is to be written to.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteItemAliases(const std::string& fileName);


    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    ISTable* GetFormatTable();
 
    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    int WriteFormatted(const std::string& cifFileName, ISTable* formatP = NULL);

    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    int WriteFormatted(const std::string& cifFileName, TableFile* ddl,
      ISTable* formatP = NULL);

    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    void Compress(CifFile* ddl);

    CifFile* GetRefFile();

  protected:
    ISTable* _formatP;

    int WriteFormatted(std::ostream& cifo, ISTable* formatP);
    int WriteFormatted(std::ostream& cifo, TableFile* ddl, ISTable* formatP);

    void WriteItemAliases(std::ostream& cifo);

  private:
    void AddRefRow(ISTable& table, const char* first, const char* second, 
      const char* third);
};

#endif
