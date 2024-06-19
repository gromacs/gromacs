/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "nm2type.h"

#include <cstring>

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/fflibutil.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp_impl.h"
#include "gromacs/gmxpreprocess/notset.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/gmxpreprocess/toppush.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

static void rd_nm2type_file(const std::filesystem::path& filename, int* nnm, t_nm2type** nmp)
{
    FILE*      fp;
    bool       bCont;
    char       libfilename[128];
    char       format[128], f1[128];
    char       buf[1024], elem[16], type[16], nbbuf[16], **newbuf;
    int        i, nb, nnnm, line = 1;
    double     qq, mm;
    t_nm2type* nm2t = nullptr;

    fp = fflib_open(filename);
    if (nullptr == fp)
    {
        gmx_fatal(FARGS, "Can not find %s in library directory", filename.string().c_str());
    }

    nnnm = *nnm;
    nm2t = *nmp;
    do
    {
        /* Read a line from the file */
        bCont = (fgets2(buf, 1023, fp) != nullptr);

        if (bCont)
        {
            /* Remove comment */
            strip_comment(buf);
            if (sscanf(buf, "%s%s%lf%lf%d", elem, type, &qq, &mm, &nb) == 5)
            {
                /* If we can read the first four, there probably is more */
                srenew(nm2t, nnnm + 1);
                snew(nm2t[nnnm].blen, nb);
                if (nb > 0)
                {
                    snew(newbuf, nb);
                    strcpy(format, "%*s%*s%*s%*s%*s");
                    for (i = 0; (i < nb); i++)
                    {
                        /* Complicated format statement */
                        strcpy(f1, format);
                        strcat(f1, "%s%lf");
                        if (sscanf(buf, f1, nbbuf, &(nm2t[nnnm].blen[i])) != 2)
                        {
                            gmx_fatal(FARGS, "Error on line %d of %s", line, libfilename);
                        }
                        newbuf[i] = gmx_strdup(nbbuf);
                        strcat(format, "%*s%*s");
                    }
                }
                else
                {
                    newbuf = nullptr;
                }
                nm2t[nnnm].elem   = gmx_strdup(elem);
                nm2t[nnnm].type   = gmx_strdup(type);
                nm2t[nnnm].q      = qq;
                nm2t[nnnm].m      = mm;
                nm2t[nnnm].nbonds = nb;
                nm2t[nnnm].bond   = newbuf;
                nnnm++;
            }
            line++;
        }
    } while (bCont);
    gmx_ffclose(fp);

    *nnm = nnnm;
    *nmp = nm2t;
}

t_nm2type* rd_nm2type(const std::filesystem::path& ffdir, int* nnm)
{
    auto ff       = fflib_search_file_end(ffdir, ".n2t", FALSE);
    *nnm          = 0;
    t_nm2type* nm = nullptr;
    for (const auto& filename : ff)
    {
        rd_nm2type_file(filename, nnm, &nm);
    }
    return nm;
}

void dump_nm2type(FILE* fp, int nnm, t_nm2type nm2t[])
{
    int i, j;

    fprintf(fp, "; nm2type database\n");
    for (i = 0; (i < nnm); i++)
    {
        fprintf(fp, "%-8s %-8s %8.4f %8.4f %-4d", nm2t[i].elem, nm2t[i].type, nm2t[i].q, nm2t[i].m, nm2t[i].nbonds);
        for (j = 0; (j < nm2t[i].nbonds); j++)
        {
            fprintf(fp, " %-5s %6.4f", nm2t[i].bond[j], nm2t[i].blen[j]);
        }
        fprintf(fp, "\n");
    }
}

enum
{
    ematchNone,
    ematchWild,
    ematchElem,
    ematchExact,
    ematchNR
};

static int match_str(const char* atom, const char* template_string)
{
    if (!atom || !template_string)
    { // NOLINT bugprone-branch-clone
        return ematchNone;
    }
    else if (gmx_strcasecmp(atom, template_string) == 0)
    {
        return ematchExact;
    }
    else if (atom[0] == template_string[0])
    {
        return ematchElem;
    }
    else if (strcmp(template_string, "*") == 0)
    {
        return ematchWild;
    }
    else
    {
        return ematchNone;
    }
}

int nm2type(int nnm, t_nm2type nm2t[], t_atoms* atoms, PreprocessingAtomTypes* atype, int* nbonds, InteractionsOfType* bonds)
{
    int cur = 0;
#define prev (1 - cur)
    int   nresolved, nb, maxbond, nqual[2][ematchNR];
    int * bbb, *n_mask, *m_mask, **match;
    char *aname_i, *aname_n;

    maxbond = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        maxbond = std::max(maxbond, nbonds[i]);
    }
    if (debug)
    {
        fprintf(debug, "Max number of bonds per atom is %d\n", maxbond);
    }
    snew(bbb, maxbond);
    snew(n_mask, maxbond);
    snew(m_mask, maxbond);
    snew(match, maxbond);
    for (int i = 0; (i < maxbond); i++)
    {
        snew(match[i], maxbond);
    }

    nresolved = 0;
    for (int i = 0; (i < atoms->nr); i++)
    {
        aname_i = *atoms->atomname[i];
        nb      = 0;
        for (const auto& bond : bonds->interactionTypes)
        {
            int ai = bond.ai();
            int aj = bond.aj();
            if (ai == i)
            {
                bbb[nb++] = aj;
            }
            else if (aj == i)
            {
                bbb[nb++] = ai;
            }
        }
        if (nb != nbonds[i])
        {
            gmx_fatal(FARGS, "Counting number of bonds nb = %d, nbonds[%d] = %d", nb, i, nbonds[i]);
        }
        if (debug)
        {
            fprintf(debug, "%4s has bonds to", aname_i);
            for (int j = 0; (j < nb); j++)
            {
                fprintf(debug, " %4s", *atoms->atomname[bbb[j]]);
            }
            fprintf(debug, "\n");
        }
        int best = -1;
        for (int k = 0; (k < ematchNR); k++)
        {
            nqual[prev][k] = 0;
        }

        /* First check for names */
        for (int k = 0; (k < nnm); k++)
        {
            if (nm2t[k].nbonds == nb)
            {
                int im = match_str(*atoms->atomname[i], nm2t[k].elem);
                if (im > ematchWild)
                {
                    for (int j = 0; (j < ematchNR); j++)
                    {
                        nqual[cur][j] = 0;
                    }

                    /* Fill a matrix with matching quality */
                    for (int m = 0; (m < nb); m++)
                    {
                        const char* aname_m = *atoms->atomname[bbb[m]];
                        for (int n = 0; (n < nb); n++)
                        {
                            aname_n     = nm2t[k].bond[n];
                            match[m][n] = match_str(aname_m, aname_n);
                        }
                    }
                    /* Now pick the best matches */
                    for (int m = 0; (m < nb); m++)
                    {
                        n_mask[m] = 0;
                        m_mask[m] = 0;
                    }
                    for (int j = ematchNR - 1; (j > 0); j--)
                    {
                        for (int m = 0; (m < nb); m++)
                        {
                            for (int n = 0; (n < nb); n++)
                            {
                                if ((n_mask[n] == 0) && (m_mask[m] == 0) && (match[m][n] == j))
                                {
                                    n_mask[n] = 1;
                                    m_mask[m] = 1;
                                    nqual[cur][j]++;
                                }
                            }
                        }
                    }
                    if ((nqual[cur][ematchExact] + nqual[cur][ematchElem] + nqual[cur][ematchWild]) == nb)
                    {
                        if ((nqual[cur][ematchExact] > nqual[prev][ematchExact]) ||

                            ((nqual[cur][ematchExact] == nqual[prev][ematchExact])
                             && (nqual[cur][ematchElem] > nqual[prev][ematchElem]))
                            ||

                            ((nqual[cur][ematchExact] == nqual[prev][ematchExact])
                             && (nqual[cur][ematchElem] == nqual[prev][ematchElem])
                             && (nqual[cur][ematchWild] > nqual[prev][ematchWild])))
                        {
                            best = k;
                            cur  = prev;
                        }
                    }
                }
            }
        }
        if (best != -1)
        {
            int  atomnr = 0;
            real alpha  = 0;

            double      qq   = nm2t[best].q;
            double      mm   = nm2t[best].m;
            const char* type = nm2t[best].type;

            auto atomType = atype->atomTypeFromName(type);
            int  k;
            if (!atomType.has_value())
            {
                atoms->atom[i].qB = alpha;
                atoms->atom[i].m = atoms->atom[i].mB = mm;
                k                                    = atype->addType(
                        atoms->atom[i], type, InteractionOfType({}, {}), atoms->atom[i].type, atomnr);
            }
            else
            {
                k = *atomType;
            }
            atoms->atom[i].type  = k;
            atoms->atom[i].typeB = k;
            atoms->atom[i].q     = qq;
            atoms->atom[i].qB    = qq;
            atoms->atom[i].m     = mm;
            atoms->atom[i].mB    = mm;
            nresolved++;
        }
        else
        {
            fprintf(stderr,
                    "Can not find forcefield for atom %s-%d with %d bonds\n",
                    *atoms->atomname[i],
                    i + 1,
                    nb);
        }
    }
    sfree(bbb);
    sfree(n_mask);
    sfree(m_mask);

    return nresolved;
}
