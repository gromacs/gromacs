/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2011- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \file
 * \brief
 * Declares common string utility and formatting routines.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRINGUTIL_H
#define GMX_UTILITY_STRINGUTIL_H

#include <cstdarg>
#include <cstring>

#include <string>
#include <vector>

namespace gmx
{

//! \addtogroup module_utility
//! \{

/*! \brief
 * Tests whether a string is null or empty.
 *
 * Does not throw.
 */
static inline bool isNullOrEmpty(const char* str)
{
    return str == nullptr || str[0] == '\0';
}

/*! \brief
 * Tests whether a string starts with another string.
 *
 * \param[in] str    String to process.
 * \param[in] prefix Prefix to find.
 * \returns   true if \p str starts with \p prefix.
 *
 * Returns true if \p prefix is empty.
 * Does not throw.
 */
static inline bool startsWith(const std::string& str, const std::string& prefix)
{
    return str.compare(0, prefix.length(), prefix) == 0;
}
//! \copydoc startsWith(const std::string &, const std::string &)
static inline bool startsWith(const char* str, const char* prefix)
{
    return std::strncmp(str, prefix, std::strlen(prefix)) == 0;
}

/*! \brief
 * Tests whether a string ends with another string.
 *
 * \param[in] str    String to process.
 * \param[in] suffix Suffix to find.
 * \returns   true if \p str ends with \p suffix.
 *
 * Returns true if \p suffix is NULL or empty.
 * Does not throw.
 */
bool endsWith(const char* str, const char* suffix);
//! \copydoc endsWith(const char *, const char *)
static inline bool endsWith(const std::string& str, const char* suffix)
{
    return endsWith(str.c_str(), suffix);
}

/*! \brief
 * Tests whether a string contains another as a substring.
 *
 * \param[in] str    String to process.
 * \param[in] substr Substring to find.
 * \returns   true if \p str contains \p substr.
 *
 * Does not throw.
 */
static inline bool contains(const std::string& str, const char* substr)
{
    return str.find(substr) != std::string::npos;
}
//! \copydoc contains(const std::string &str, const char *substr)
static inline bool contains(const std::string& str, const std::string& substr)
{
    return str.find(substr) != std::string::npos;
}

/*!\brief Returns number of space-separated words in zero-terminated char ptr
 *
 * \param s Character pointer to zero-terminated, which will not be changed.
 *
 * \returns number of words in string.
 *
 * \note This routine is mainly meant to support legacy code in GROMACS. For
 *       new source you should try hard to use C++ string objects instead.
 */
std::size_t countWords(const char* s);

/*!\brief Returns the number of space-separated words in a string object
 *
 * \param str Reference to string object, which will not be changed.
 *
 * \returns number of words in string.
 */
std::size_t countWords(const std::string& str);

//! \copydoc endsWith(const std::string &str, const char *suffix)
static inline bool endsWith(const std::string& str, const std::string& suffix)
{
    return endsWith(str, suffix.c_str());
}

/*! \brief
 * Removes a suffix from a string.
 *
 * \param[in] str    String to process.
 * \param[in] suffix Suffix to remove.
 * \returns   \p str with \p suffix removed, or \p str unmodified if it does
 *      not end with \p suffix.
 * \throws    std::bad_alloc if out of memory.
 *
 * Returns \p str if \p suffix is NULL or empty.
 */
std::string stripSuffixIfPresent(const std::string& str, const char* suffix);
/*! \brief
 * Removes leading and trailing whitespace from a string.
 *
 * \param[in] str  String to process.
 * \returns   \p str with leading and trailing whitespaces removed.
 * \throws    std::bad_alloc if out of memory.
 */
std::string stripString(const std::string& str);
#ifdef __GNUC__
#    define gmx_format(archetype, string_index, first_to_check) \
        __attribute__((format(archetype, string_index, first_to_check)))
#else
/*! \brief GCC like function format attribute
 *
 * The format attribute specifies that a function takes printf, scanf, ...
 * style arguments that should be type-checked against a format string.
 * The attribute has to be placed after the function.
 * This attribute is only valid for function declarations and not function
 * definitions (GCC limitation). For member functions the implicit `this`
 * pointer is included in the argument count.
 */
#    define gmx_format(archetype, string_index, first_to_check)
#endif
#ifdef _MSC_VER
#    define gmx_fmtstr _In_ _Printf_format_string_
#else
/*! \brief MSVC like function format attribute
 *
 * Does type checking for printf like format strings in MSVC style.
 * Attribute has to be placed before format string.
 */
#    define gmx_fmtstr
#endif
/*! \brief
 * Formats a string (snprintf() wrapper).
 *
 * \throws  std::bad_alloc if out of memory.
 *
 * This function works like sprintf(), except that it returns an std::string
 * instead of requiring a preallocated buffer.  Arbitrary length output is
 * supported.
 */
std::string formatString(gmx_fmtstr const char* fmt, ...) gmx_format(printf, 1, 2);

/*! \brief
 * Formats a string (vsnprintf() wrapper).
 *
 * \throws  std::bad_alloc if out of memory.
 *
 * This function works like vsprintf(), except that it returns an std::string
 * instead of requiring a preallocated buffer.  Arbitrary length output is
 * supported.
 */
std::string formatStringV(const char* fmt, std::va_list ap);

/*! \brief Function object that wraps a call to formatString() that
 * expects a single conversion argument, for use with algorithms. */
class StringFormatter
{
public:
    /*! \brief Constructor
     *
     * \param[in] format The printf-style format string that will
     *     be applied to convert values of type T to
     *     string. Exactly one argument to the conversion
     *     specification(s) in `format` is supported. */
    explicit StringFormatter(const char* format) : format_(format) {}

    //! Implements the formatting functionality
    template<typename T>
    std::string operator()(const T& value) const
    {
        return formatString(format_, value);
    }

private:
    //! Format string to use
    const char* format_;
};

/*! \brief Function object to implement the same interface as
 * `StringFormatter` to use with strings that should not be formatted
 * further. */
class IdentityFormatter
{
public:
    //! Implements the formatting non-functionality
    std::string operator()(const std::string& value) const { return value; }
};

/*! \brief Formats all the range as strings, and then joins them with
 * a separator in between.
 *
 * \param[in] begin      Iterator the beginning of the range to join.
 * \param[in] end        Iterator the end of the range to join.
 * \param[in] separator  String to put in between the joined strings.
 * \param[in] formatter  Function object to format the objects in
 *     `container` as strings
 * \returns   All objects in the range from `begin` to `end` formatted
 *     as strings and concatenated with `separator` between each pair.
 * \throws    std::bad_alloc if out of memory.
 */
template<typename InputIterator, typename FormatterType>
std::string formatAndJoin(InputIterator begin, InputIterator end, const char* separator, const FormatterType& formatter)
{
    std::string result;
    const char* currentSeparator = "";
    for (InputIterator i = begin; i != end; ++i)
    {
        result.append(currentSeparator);
        result.append(formatter(*i));
        currentSeparator = separator;
    }
    return result;
}

/*! \brief Formats all elements of the container as strings, and then
 * joins them with a separator in between.
 *
 * \param[in] container  Objects to join.
 * \param[in] separator  String to put in between the joined strings.
 * \param[in] formatter  Function object to format the objects in
 *     `container` as strings
 * \returns   All objects from `container` formatted as strings and
 *     concatenated with `separator` between each pair.
 * \throws    std::bad_alloc if out of memory.
 */
template<typename ContainerType, typename FormatterType>
std::string formatAndJoin(const ContainerType& container, const char* separator, const FormatterType& formatter)
{
    return formatAndJoin(container.begin(), container.end(), separator, formatter);
}

/*! \brief
 * Joins strings from a range with a separator in between.
 *
 * \param[in] begin      Iterator the beginning of the range to join.
 * \param[in] end        Iterator the end of the range to join.
 * \param[in] separator  String to put in between the joined strings.
 * \returns   All strings from (`begin`, `end`) concatenated with `separator`
 *     between each pair.
 * \throws    std::bad_alloc if out of memory.
 */
template<typename InputIterator>
std::string joinStrings(InputIterator begin, InputIterator end, const char* separator)
{
    return formatAndJoin(begin, end, separator, IdentityFormatter());
}

/*! \brief
 * Joins strings from a container with a separator in between.
 *
 * \param[in] container  Strings to join.
 * \param[in] separator  String to put in between the joined strings.
 * \returns   All strings from `container` concatenated with `separator`
 *     between each pair.
 * \throws    std::bad_alloc if out of memory.
 */
template<typename ContainerType>
std::string joinStrings(const ContainerType& container, const char* separator)
{
    return joinStrings(container.begin(), container.end(), separator);
}

/*! \brief
 * Joins strings from an array with a separator in between.
 *
 * \param[in] array      Array of strings to join.
 * \param[in] separator  String to put in between the joined strings.
 * \tparam    count      Deduced number of elements in \p array.
 * \returns   All strings from `aray` concatenated with `separator`
 *     between each pair.
 * \throws    std::bad_alloc if out of memory.
 */
template<size_t count>
std::string joinStrings(const char* const (&array)[count], const char* separator)
{
    return joinStrings(array, array + count, separator);
}

/*! \brief
 * Splits a string to whitespace separated tokens.
 *
 * \param[in] str  String to process.
 * \returns   \p str split into tokens at each whitespace sequence.
 * \throws    std::bad_alloc if out of memory.
 *
 * This function works like `split` in Python, i.e., leading and trailing
 * whitespace is ignored, and consecutive whitespaces are treated as a single
 * separator.
 */
std::vector<std::string> splitString(const std::string& str);
/*! \brief
 * Splits a string to tokens separated by a given delimiter.
 *
 * \param[in] str   String to process.
 * \param[in] delim Delimiter to use for splitting.
 * \returns   \p str split into tokens at delimiter.
 * \throws    std::bad_alloc if out of memory.
 *
 * Unlike splitString(), consecutive delimiters will generate empty tokens, as
 * will leading or trailing delimiters.
 * Empty input will return an empty vector.
 */
std::vector<std::string> splitDelimitedString(const std::string& str, char delim);
/*! \brief
 * Splits \c str to tokens separated by delimiter \c delim. Removes
 * leading and trailing whitespace from those strings with std::isspace.
 *
 * \param[in] str   String to process.
 * \param[in] delim Delimiter to use for splitting.
 * \returns   \p str split into tokens at delimiter, with whitespace stripped.
 * \throws    std::bad_alloc if out of memory.
 *
 * Unlike splitString(), consecutive delimiters will generate empty tokens, as
 * will leading or trailing delimiters.
 * Empty input will return an empty vector.
 * Input with only whitespace will return a vector of size 1,
 * that contains an empty token.
 */
std::vector<std::string> splitAndTrimDelimitedString(const std::string& str, char delim);

/*! \brief
 * Replace all occurrences of a string with another string.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \returns   Copy of \p input with all occurrences of \p from replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 *
 * The replacement is greedy and not recursive: starting from the beginning of
 * \p input, each match of \p from is replaced with \p to, and the search for
 * the next match begins after the end of the previous match.
 *
 * Complexity is O(N), where N is length of output.
 *
 * \see replaceAllWords()
 */
std::string replaceAll(const std::string& input, const char* from, const char* to);
//! \copydoc replaceAll(const std::string &, const char *, const char *)
std::string replaceAll(const std::string& input, const std::string& from, const std::string& to);
/*! \brief
 * Replace whole words with others.
 *
 * \param[in] input  Input string.
 * \param[in] from   String to find.
 * \param[in] to     String to use to replace \p from.
 * \returns   Copy of \p input with all \p from words replaced with \p to.
 * \throws    std::bad_alloc if out of memory.
 *
 * Works as replaceAll(), but a match is only considered if it is delimited by
 * non-alphanumeric characters.
 *
 * \see replaceAll()
 */
std::string replaceAllWords(const std::string& input, const char* from, const char* to);
//! \copydoc replaceAllWords(const std::string &, const char *, const char *)
std::string replaceAllWords(const std::string& input, const std::string& from, const std::string& to);

/*! \brief Return whether two strings are equal, ignoring case.
 *
 * Checks if two strings have the same length and if all characters
 * in them match when compared case insensitive.
 * Characters are converted by using std::tolower.
 *
 * \param[in] source Search string to compare against \p target.
 * \param[in] target String to be matched to \p source.
 * \returns True if the strings match.
 */
bool equalCaseInsensitive(const std::string& source, const std::string& target);

/*! \brief Return whether too strings are case sensitive equal, ignoring dashes.
 */
bool equalIgnoreDash(const std::string& source, const std::string& target);

//! Function object for comparisons with \c equalCaseInsensitive
class EqualCaseInsensitive
{
public:
    bool operator()(const std::string& lhs, const std::string& rhs) const
    {
        return gmx::equalCaseInsensitive(lhs, rhs);
    }
};

/*! \brief
 * Checks if at most \p maxLengthOfComparison characters of two strings match case insensitive.
 *
 * The function tests two strings \p source and \p target to see if at most
 * \p maxLengthOfComparison characters match between the two. If fewer characters are present
 * in \p source, only the maximum number of characters in \p source will be compared instead.
 * In this case both \p source and \p target also need to have the same length, or the strings will
 * compare as false, even if \p target matches \p source over the length of \p source.
 *
 * If \p maxLengthOfComparison is 0, the function always returns true.
 * Characters are converted by using std::tolower.
 *
 * \param[in] source Search string to compare against \p target.
 * \param[in] target String to be matched to \p source.
 * \param[in] maxLengthOfComparison The maximum string length to compare.
 * \returns True if the strings match.
 */
bool equalCaseInsensitive(const std::string& source, const std::string& target, size_t maxLengthOfComparison);

/*! \brief
 * Makes the string uppercase.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all characters transformed to uppercase.
 * \throws    std::bad_alloc if out of memory.
 */
std::string toUpperCase(const std::string& text);

/*! \brief
 * Makes the string lowercase.
 *
 * \param[in] text  Input text.
 * \returns   \p text with all characters transformed to lowercase.
 * \throws    std::bad_alloc if out of memory.
 */
std::string toLowerCase(const std::string& text);


class TextLineWrapper;

/*! \brief
 * Stores settings for line wrapping.
 *
 * Methods in this class do not throw.
 *
 * \see TextLineWrapper
 *
 * \inpublicapi
 */
class TextLineWrapperSettings
{
public:
    /*! \brief
     * Initializes default wrapper settings.
     *
     * Default settings are:
     *  - No maximum line width (only explicit line breaks).
     *  - No indentation.
     *  - No continuation characters.
     *  - Do not keep final spaces in input strings.
     */
    TextLineWrapperSettings();

    /*! \brief
     * Sets the maximum length for output lines.
     *
     * \param[in] length  Maximum length for the lines after wrapping.
     *
     * If this method is not called, or is called with zero \p length, the
     * wrapper has no maximum length (only wraps at explicit line breaks).
     */
    void setLineLength(int length) { maxLength_ = length; }
    /*! \brief
     * Sets the indentation for output lines.
     *
     * \param[in] indent  Number of spaces to add for indentation.
     *
     * If this method is not called, the wrapper does not add indentation.
     */
    void setIndent(int indent) { indent_ = indent; }
    /*! \brief
     * Sets the indentation for first output line after a line break.
     *
     * \param[in] indent  Number of spaces to add for indentation.
     *
     * If this method is not called, or called with \p indent equal to -1,
     * the value set with setIndent() is used.
     */
    void setFirstLineIndent(int indent) { firstLineIndent_ = indent; }
    /*! \brief
     * Sets whether final spaces in input should be kept.
     *
     * \param[in] bKeep  Whether to keep spaces at the end of the input.
     *
     * This means that wrapping a string that ends in spaces also keeps
     * those spaces in the output.  This allows using the wrapper for
     * partial lines where the initial part of the line may end in a space.
     * By default, all trailing whitespace is removed.  Note that this
     * option does not affect spaces before an explicit newline: those are
     * always removed.
     */
    void setKeepFinalSpaces(bool bKeep) { bKeepFinalSpaces_ = bKeep; }
    /*! \brief
     * Sets a continuation marker for wrapped lines.
     *
     * \param[in] continuationChar  Character to use to mark continuation
     *      lines.
     *
     * If set to non-zero character code, this character is added at the
     * end of each line where a line break is added by TextLineWrapper
     * (but not after lines produced by explicit line breaks).
     * The default (\c '\0') is to not add continuation markers.
     *
     * Note that currently, the continuation char may cause the output line
     * length to exceed the value set with setLineLength() by at most two
     * characters.
     */
    void setContinuationChar(char continuationChar) { continuationChar_ = continuationChar; }

    //! Returns the maximum length set with setLineLength().
    int lineLength() const { return maxLength_; }
    //! Returns the indentation set with setIndent().
    int indent() const { return indent_; }
    /*! \brief
     * Returns the indentation set with setFirstLineIndent().
     *
     * If setFirstLineIndent() has not been called or has been called with
     * -1, indent() is returned.
     */
    int firstLineIndent() const { return (firstLineIndent_ >= 0 ? firstLineIndent_ : indent_); }

private:
    //! Maximum length of output lines, or <= 0 if no limit.
    int maxLength_;
    //! Number of spaces to indent each output line with.
    int indent_;
    /*! \brief
     * Number of spaces to indent the first line after a newline.
     *
     * If -1, \a indent_ is used.
     */
    int firstLineIndent_;
    //! Whether to keep spaces at end of input.
    bool bKeepFinalSpaces_;
    //! If not \c '\0', mark each wrapping point with this character.
    char continuationChar_;

    //! Needed to access the members.
    friend class TextLineWrapper;
};

/*! \brief
 * Wraps lines to a predefined length.
 *
 * This utility class wraps lines at word breaks to produce lines that are not
 * longer than a predefined length.  Explicit newlines ('\\n') are preserved.
 * Only space is considered a word separator.  If a single word exceeds the
 * maximum line length, it is still printed on a single line.
 * Extra whitespace is stripped from the end of produced lines.
 * Other options on the wrapping, such as the line length or indentation,
 * can be changed using a TextLineWrapperSettings object.
 *
 * Two interfaces to do the wrapping are provided:
 *  -# High-level interface using either wrapToString() (produces a single
 *     string with embedded newlines) or wrapToVector() (produces a vector of
 *     strings with each line as one element).
 *     These methods operate on std::string and wrap the entire input string.
 *  -# Low-level interface using findNextLine() and formatLine().
 *     findNextLine() operates either on a C string or an std::string, and does
 *     not do any memory allocation (so it does not throw).  It finds the next
 *     line to be wrapped, considering the wrapping settings.
 *     formatLine() does whitespace operations on the line found by
 *     findNextLine() and returns an std::string.
 *     These methods allow custom wrapping implementation to either avoid
 *     exceptions or to wrap only a part of the input string.
 *
 * Typical usage:
 * \code
   gmx::TextLineWrapper wrapper;
   wrapper.settings().setLineLength(78);
   printf("%s\n", wrapper.wrapToString(textToWrap).c_str());
   \endcode
 *
 * \inpublicapi
 */
class TextLineWrapper
{
public:
    /*! \brief
     * Constructs a new line wrapper with default settings.
     *
     * Does not throw.
     */
    TextLineWrapper() {}
    /*! \brief
     * Constructs a new line wrapper with given settings.
     *
     * \param[in] settings  Wrapping settings.
     *
     * Does not throw.
     */
    explicit TextLineWrapper(const TextLineWrapperSettings& settings) : settings_(settings) {}

    /*! \brief
     * Provides access to settings of this wrapper.
     *
     * \returns  The settings object for this wrapper.
     *
     * The returned object can be used to modify settings for the wrapper.
     * All subsequent calls to wrapToString() and wrapToVector() use the
     * modified settings.
     *
     * Does not throw.
     */
    TextLineWrapperSettings& settings() { return settings_; }

    //! Returns true if the wrapper would not modify the input string.
    bool isTrivial() const;

    /*! \brief
     * Finds the next line to be wrapped.
     *
     * \param[in] input     String to wrap.
     * \param[in] lineStart Index of first character of the line to find.
     * \returns   Index of first character of the next line.
     *
     * If this is the last line, returns the length of \p input.
     * In determining the length of the returned line, this function
     * considers the maximum line length, leaving space for indentation,
     * and also whitespace stripping behavior.
     * Thus, the line returned may be longer than the maximum line length
     * if it has leading and/or trailing space.
     * When wrapping a line on a space (not on an explicit line break),
     * the returned index is always on a non-whitespace character after the
     * space.
     *
     * To iterate over lines in a string, use the following code:
     * \code
       gmx::TextLineWrapper wrapper;
       // <set desired wrapping settings>
       size_t lineStart = 0;
       size_t length = input.length();
       while (lineStart < length)
       {
           size_t nextLineStart = wrapper.findNextLine(input, lineStart);
           std::string line = wrapper.formatLine(input, lineStart, nextLineStart));
           // <do something with the line>
           lineStart = nextLineStart;
       }
       return result;
       \endcode
     *
     * Does not throw.
     */
    size_t findNextLine(const char* input, size_t lineStart) const;
    //! \copydoc findNextLine(const char *, size_t)const
    size_t findNextLine(const std::string& input, size_t lineStart) const;
    /*! \brief
     * Formats a single line for output according to wrapping settings.
     *
     * \param[in] input     Input string.
     * \param[in] lineStart Index of first character of the line to format.
     * \param[in] lineEnd   Index of first character of the next line.
     * \returns   The line with leading and/or trailing whitespace removed
     *      and indentation applied.
     * \throws    std::bad_alloc if out of memory.
     *
     * Intended to be used on the lines found by findNextLine().
     * When used with the lines returned from findNextLine(), the returned
     * line conforms to the wrapper settings.
     * Trailing whitespace is always stripped (including any newlines,
     * i.e., the return value does not contain a newline).
     */
    std::string formatLine(const std::string& input, size_t lineStart, size_t lineEnd) const;

    /*! \brief
     * Formats a string, producing a single string with all the lines.
     *
     * \param[in] input  String to wrap.
     * \returns   \p input with added newlines such that maximum line
     *      length is not exceeded.
     * \throws    std::bad_alloc if out of memory.
     *
     * Newlines in the input are preserved, including terminal newlines.
     * Note that if the input does not contain a terminal newline, the
     * output does not either.
     */
    std::string wrapToString(const std::string& input) const;
    /*! \brief
     * Formats a string, producing a vector with all the lines.
     *
     * \param[in] input  String to wrap.
     * \returns   \p input split into lines such that maximum line length
     *      is not exceeded.
     * \throws    std::bad_alloc if out of memory.
     *
     * The strings in the returned vector do not contain newlines at the
     * end.
     * Note that a single terminal newline does not affect the output:
     * "line\\n" and "line" both produce the same output (but "line\\n\\n"
     * produces two lines, the second of which is empty).
     */
    std::vector<std::string> wrapToVector(const std::string& input) const;

private:
    TextLineWrapperSettings settings_;
};

//! \}

} // namespace gmx

#endif
