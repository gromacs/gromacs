/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLSELECT_H
#define MOLSELECT_H

struct gmx_molselect;

enum iMolSelect {
    imsTrain, imsTest, imsIgnore, imsUnknown, imsNR
};

extern const char *ims_names[imsNR];

gmx_molselect *gmx_molselect_init(const char *fn);

void gmx_molselect_done(gmx_molselect *gms);

iMolSelect gmx_molselect_status(gmx_molselect *gms, const char *iupac);

int gmx_molselect_index(gmx_molselect *gms, const char *iupac);

#endif
