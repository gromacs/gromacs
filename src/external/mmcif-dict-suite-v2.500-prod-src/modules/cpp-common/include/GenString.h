//$$FILE$$
//$$VERSION$$
//$$DATE$$
//$$LICENSE$$


#ifndef GENSTRING_H
#define GENSTRING_H


#include <string>
#include <functional>


/**
 ** \class Char
 **
 ** \brief Generic character class that contains character related methods.
 **
 ** This class is a static class that contains generic character related utility
 ** methods.
 */
class Char
{
  public:
    enum eCompareType
    {
        eCASE_SENSITIVE = 0,
        eCASE_INSENSITIVE,
        eWS_INSENSITIVE,  // But case-sensitive
        eAS_INTEGER
    };

    static char ToLower(const char c);
    static char ToUpper(const char c);

    static bool IsCiLess(const char c1, const char c2);

    static bool IsWhiteSpace(const char c);
    static bool IsDigit(const char c);
    static bool IsCarriageReturn(const char c);
    static bool IsPrintable(const char c);

    static void AsciiCodeInHex(const char c, std::string& asciiHexString);
};


/**
 ** \class CharLess
 **
 ** \brief Public class that encapsulates character comparison.
 **
 ** This class encapsulates character comparison. It supports the following
 ** compare types: case-sensitive and case-insensitive.
 */
class CharLess
{
  public:
    CharLess(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    CharLess& operator=(const CharLess& in);

    bool operator()(const char c1, const char c2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class CharEqualTo
 **
 ** \brief Public class that encapsulates generic character equal_to functor.
 **
 ** This class is equal_to functor for generic character. It supports the
 ** following compare types: case-sensitive and case-insensitive.
 */
class CharEqualTo : public std::binary_function<char, char, bool>
{
  public:
    CharEqualTo(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    CharEqualTo& operator=(const CharEqualTo& in);

    bool operator()(const char c1, const char c2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


class WhiteSpace : public std::unary_function<char, bool>
{
  public:
    bool operator()(const char c) const;
    bool operator()(const char c1, const char c2) const;
};


/**
 ** \class StringLess
 **
 ** \brief Public class that encapsulates string comparison.
 **
 ** This class encapsulates string comparison. It supports the following
 ** compare types: case-sensitive, case-insensitive and as-integer.
 */
class StringLess
{
  public:
    StringLess(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    StringLess& operator=(const StringLess& in);

    bool operator()(const std::string& s1, const std::string& s2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class StringEqualTo
 **
 ** \brief Public class that encapsulates generic string equal_to functor.
 **
 ** This class is equal_to functor for generic strings. It supports the
 ** following compare types: case-sensitive, case-insensitive and as-integer.
 */
class StringEqualTo : public std::binary_function<std::string, std::string,
  bool>
{
  public:
    StringEqualTo(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    StringEqualTo& operator=(const StringEqualTo& in);

    bool operator()(const std::string& s1, const std::string& s2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class String
 **
 ** \brief Generic string class that contains string related utility methods.
 **
 ** This class is a static class that contains generic string related utility
 ** methods, such as: converting string to uppercase/lowercase, removing
 ** whitespaces, converting strings to/from integers/real numbers, determining
 ** if string a number, determining whether strings are equal, escaping and
 ** unescaping.
 */
class String
{
  public:
    static void LowerCase(const std::string& inString, std::string& outString);
    static void LowerCase(std::string& inOutString);
    static void UpperCase(const std::string& inString, std::string& outString);
    static void UpperCase(std::string& inOutString);

    static void RemoveWhiteSpace(const std::string& inString,
      std::string& outString);

    static std::string IntToString(int inInteger);
    static std::string DoubleToString(double inDouble);
    static int StringToInt(const std::string& inString);
    static double StringToDouble(const std::string& inString);
    static bool IsScientific(const std::string& number);
    static void ToFixedFormat(std::string& fixedFormat,
      const std::string& number);
    static bool StringToBoolean(const std::string& inString);

    static bool IsNumber(const std::string& inString);

    static bool IsCiEqual(const std::string& firstString,
      const std::string& secondString);
    static bool IsEqual(const std::string& firstString,
      const std::string& secondString,
      const Char::eCompareType compareType);

    static void StripLeadingWs(std::string& resString);
    static void StripTrailingWs(std::string& resString);
    static void StripAndCompressWs(std::string& resString);
    static void rcsb_clean_string(std::string& theString);

    static void UnEscape(std::string& outStr, const std::string& inStr);

    static void Replace(std::string& resString, const std::string& fromStr,
      const std::string& toStr);

  private:
    static std::string::const_iterator GetExpValue(int& expValue,
      const std::string::const_iterator& beg,
      const std::string::const_iterator& end);
    static void GetMantissa(std::string& mantissa, int& addExpValue,
      const std::string::const_iterator& beg,
      const std::string::const_iterator& end);
    static void ScientificNumberToFixed(std::string& fixed,
      const bool isPositive, const std::string& mantissa, const int exponent);
};


inline Char::eCompareType StringLess::GetCompareType()
{
    return (_compareType);
}

inline Char::eCompareType StringEqualTo::GetCompareType()
{
    return (_compareType);
}

#endif
