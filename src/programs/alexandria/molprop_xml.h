/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLPROP_XML_H
#define MOLPROP_XML_H

#include "molprop.h"

/*! \brief
 * Write a vector of molprops to an XML file
 *
 * Function that uses the libxml2 library to store an array of molprops
 * to file.
 *
 * \todo Implement using a serialized protocol rather than reading the
 * whole file into memory. Since we have a vector of molprops this should
 * be straightforward.
 * \param[in] fn         The file name to write to
 * \param[in] mpt        The vector of MolProp
 * \param[in] bCompress  Determines whether zlib compression is used when writing
 * \ingroup module_alexandria
 */
void MolPropWrite(const char                            *fn,
                  const std::vector<alexandria::MolProp> mpt,
                  gmx_bool                               bCompress);

/*! \brief
 * Reads a vector of molprops from an XML file
 *
 * Function that uses the libxml2 library to store an array of molprops
 * to file.
 *
 * \todo Implement using a serialized protocol rather than reading the
 * whole file into memory. Since we have a vector of molprops this should
 * be straightforward.
 * \param[in]  fn         The file name to read from
 * \param[out] mpt        The vector of MolProp
 * \ingroup module_alexandria
 */
void MolPropRead(const char                       *fn,
                 std::vector<alexandria::MolProp> &mpt);

#endif
