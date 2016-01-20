/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef POLDATA_XML_H
#define POLDATA_XML_H

#include <string>

#include "gromacs/topology/atomprop.h"

namespace alexandria
{
    class Poldata;

    /*! \brief Store the Poldata class to an XML file
     *
     * \param[in] fileName The filename to save to
     * \param[in] pd       Pointer to a Poldata class instance
     * \param[in] compress Whether or not to write a compressed file
     */
    void writePoldata(const std::string &fileName,
                      const Poldata &pd,
                      bool compress = true);

    /*! \brief Read a Poldata class from an XML file
     *
     * \param[in]  fileName The filename to read from
     * \param[out] pd       The Poldata class instance
     * \param[in]  aps      Atom properties
     */
    void readPoldata(const std::string &fileName,
                     Poldata &pd,
                     const gmx_atomprop_t aps);

} // namespace alexandria

#endif
