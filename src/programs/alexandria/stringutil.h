#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <string>
#include <vector>

/*! \brief
 * Return the index of the selected term
 * \param[in] opts List of strings to select from
 * \returns the index or -1 if not found
 */
int get_option(const char **opts);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \param[in] elems A vector of substring elements
 * \returns The vector of substring elements (same as input elems)
 */
std::vector<std::string> &split(const std::string        &s,
                                char                      delim,
                                std::vector<std::string> &elems);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \returns A vector of substring elements
 */
std::vector<std::string> split(const std::string &s,
                               char               delim);

/** return new string with f printed. */
extern std::string gmx_ftoa(double f);

/** return new string with f printed. */
extern std::string gmx_itoa(int f);

#endif
