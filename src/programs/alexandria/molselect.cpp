/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "molselect.h"

#include <stdlib.h>
#include <string.h>

#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"

#include "stringutil.h"

const char *ims_names[imsNR] = { "Train", "Test", "Ignore", "Unknown" };

struct t_ims 
{
    char      *iupac;
    iMolSelect status;
    int        index;
};

struct gmx_molselect 
{
    int    nmol;
    t_ims *ims;
};

static int ims_comp(const void *a, const void *b)
{
    t_ims *ia = (t_ims *)a;
    t_ims *ib = (t_ims *)b;

    return strcasecmp(ia->iupac, ib->iupac);
}

gmx_molselect *gmx_molselect_init(const char *fn)
{
    gmx_molselect *gms;
    char         **strings;
    int            i, j;

    snew(gms, 1);
    gms->nmol = get_lines(fn, &strings);
    snew(gms->ims, gms->nmol);
    for (i = 0; (i < gms->nmol); i++)
    {
        std::vector<std::string> ptr = split(strings[i], '|');
        if ((ptr.size() == 2) && (ptr[0].length() > 0))
        {
            gms->ims[i].iupac   = strdup(ptr[0].c_str());
            gms->ims[i].index   = i+1;
            if (ptr[1].length() > 0)
            {
                for (j = 0; (j < (int)imsNR); j++)
                {
                    if (strcasecmp(ims_names[j], ptr[1].c_str()) == 0)
                    {
                        break;
                    }
                }
                if (j < imsNR)
                {
                    gms->ims[i].status = (iMolSelect)j;
                }
                else
                {
                    gms->ims[i].status = imsUnknown;
                    fprintf(stderr, "Unknown status '%s' for molecule %s on line %d in file %s\n",
                            ptr[1].c_str(), ptr[0].c_str(), i, fn);
                }
            }
            else
            {
                gms->ims[i].status = imsUnknown;
                fprintf(stderr, "No status field for molecule %s on line %d in file %s\n",
                        ptr[0].c_str(), i, fn);
            }
        }
        sfree(strings[i]);
    }
    sfree(strings);

    /* Sort the molecules for faster searching down below */
    qsort(gms->ims, gms->nmol, sizeof(gms->ims[0]), ims_comp);

    return gms;
}

void gmx_molselect_done(gmx_molselect *gms)
{
    gmx_molselect *g = (gmx_molselect *)gms;
    int            i;

    for (i = 0; (i < g->nmol); i++)
    {
        sfree(g->ims[i].iupac);
    }
    sfree(g->ims);
    g->ims = NULL;
    sfree(g);
}

iMolSelect gmx_molselect_status(gmx_molselect *gms, const char *iupac)
{
    gmx_molselect *g = (gmx_molselect *)gms;
    t_ims          key;
    t_ims         *ims;

    if (NULL == iupac)
    {
        return imsUnknown;
    }
    key.iupac = strdup(iupac);
    ims       = (t_ims *)bsearch(&key, g->ims, g->nmol, sizeof(g->ims[0]), ims_comp);
    sfree(key.iupac);
    if (NULL != ims)
    {
        return ims->status;
    }
    else
    {
        return imsUnknown;
    }
}

int gmx_molselect_index(gmx_molselect *gms, const char *iupac)
{
    if (NULL == iupac)
    {
        fprintf(stderr, "iupac == NULL in call to gmx_molselect_index.\n");
    }
    else if (0 == strlen(iupac))
    {
        fprintf(stderr, "strlen(iupac) == 0 in call to gmx_molselect_index.\n");
    }
    else
    {
        t_ims  key;
        t_ims *ims;

        key.iupac  = strdup(iupac);
        key.status = imsTest;
        key.index  = 0;
        ims        = (t_ims *)bsearch(&key, gms->ims, gms->nmol, sizeof(gms->ims[0]), ims_comp);
        sfree(key.iupac);
        if (NULL != ims)
        {
            return ims->index;
        }
        else
        {
            fprintf(stderr, "Could not find %s in gmx_molselect_index.\n", iupac);
        }
    }
    return -1;
}
