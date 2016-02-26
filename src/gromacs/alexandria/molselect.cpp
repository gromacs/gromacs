/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "molselect.h"

#include <algorithm>
#include <vector>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"

#include "stringutil.h"

static const char *ims_names[imsNR] = { "Train", "Test", "Ignore", "Unknown" };

const char *iMolSelectName(iMolSelect ims)
{
    int i = static_cast<int>(ims);

    return ims_names[i];
}

namespace alexandria
{

void MolSelect::read(const char *fn)
{
    gmx::TextReader tr(fn);
    std::string     tmp;
    int             index = 0;

    while(tr.readLine(&tmp)) 
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n') 
        {
            tmp.erase(tmp.length()-1);
        }
        std::vector<std::string> ptr = split(tmp, '|');
        if ((ptr.size() == 2) && (ptr[0].length() > 1)
            && (ptr[1].length() > 1))
        {
            iMolSelect ii;
            int        j;

            for(j = 0; (j < (int)imsNR); j++)
            {
                if (strcasecmp(ims_names[j], ptr[1].c_str()) == 0)
                {
                    break;
                }
            }
            if (j < imsNR)
            {
                ii = static_cast<iMolSelect>(j);
            }
            else
            {
                ii = imsUnknown;
                fprintf(stderr, "Unknown status '%s' for molecule %s on line %d in file %s\n",
                        ptr[1].c_str(), ptr[0].c_str(), index, fn);
            }
            ims_.push_back(IMolSelect(ptr[0], ii, index++));
        }
    }
}

iMolSelect MolSelect::status(const std::string &iupac) const
{
    std::vector<IMolSelect>::const_iterator imi = 
        std::find_if(ims_.begin(), ims_.end(),
                     [iupac](IMolSelect const &i)
                     { return i.iupac().compare(iupac); });

    if (imi != ims_.end())
    {
        return imi->status();
    }

    return imsUnknown;
}

int MolSelect::index(const std::string &iupac) const
{
    std::vector<IMolSelect>::const_iterator imi = 
        std::find_if(ims_.begin(), ims_.end(),
                     [iupac](IMolSelect const &i)
                     { return i.iupac().compare(iupac); });

    if (imi != ims_.end())
    {
        return imi->index();
    }

    return -1;
}

}
