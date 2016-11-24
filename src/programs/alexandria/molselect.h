/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLSELECT_H
#define MOLSELECT_H

#include <string>
#include <vector>

enum iMolSelect {
    imsTrain, imsTest, imsIgnore, imsUnknown, imsNR
};

const char *iMolSelectName(iMolSelect ims);

namespace alexandria
{

class IMolSelect 
{
private:
    std::string iupac_;
    iMolSelect  status_;
    int         index_;
public:
    IMolSelect(const std::string &iupac, iMolSelect status, int index) :
        iupac_(iupac), status_(status), index_(index) {}

    const std::string &iupac() const { return iupac_; }

    iMolSelect status() const { return status_; }

    int index() const { return index_; }
};

class MolSelect 
{
private:
    std::vector<IMolSelect> ims_;

public:
    MolSelect() {};

    void read(const char *filename);

    size_t nMol() const { return ims_.size(); }

    iMolSelect status(const std::string &iupac) const;

    int index(const std::string &iupac) const;

    int count(iMolSelect ims) const
    {
        return std::count_if(ims_.begin(), ims_.end(),
                             [ims](IMolSelect const &i)
                             { return i.status() == ims; });
    }
};

} // namespace

#endif
