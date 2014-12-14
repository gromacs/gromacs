/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include "gmxpre.h"

#include "config.h"

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_IO_H
#include <io.h>
#endif

#include "gromacs/fileio/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/gmxfio_int.h"
#include "gromacs/fileio/md5.h"
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"


/* This is the part that reads dummy and ascii files.  */




/* file type functions */
static gmx_bool do_ascread(t_fileio *fio, void *item, int nitem, int eio,
                           const char *desc, const char *srcfile, int line);
static gmx_bool do_ascwrite(t_fileio *fio, const void *item, int nitem, int eio,
                            const char *desc, const char *srcfile, int line);
static gmx_bool do_dummyread(t_fileio *fio, void *item, int nitem, int eio,
                             const char *desc, const char *srcfile, int line);
static gmx_bool do_dummywrite(t_fileio *fio, const void *item, int nitem, int eio,
                              const char *desc, const char *srcfile, int line);


const t_iotype asc_iotype = {do_ascread, do_ascwrite};
const t_iotype dummy_iotype = {do_dummyread, do_dummywrite};






static gmx_bool do_dummyread(t_fileio gmx_unused *fio, void gmx_unused *item, int gmx_unused nitem, int gmx_unused eio,
                             const char gmx_unused *desc, const char gmx_unused *srcfile, int gmx_unused line)
{
    gmx_fatal(FARGS, "File type not set!");
    return FALSE;
}

static gmx_bool do_dummywrite(t_fileio gmx_unused *fio, const void gmx_unused *item, int gmx_unused nitem, int gmx_unused eio,
                              const char gmx_unused *desc, const char gmx_unused *srcfile, int gmx_unused line)
{
    gmx_fatal(FARGS, "File type not set!");
    return FALSE;
}



static void encode_string(int maxlen, char dst[], const char src[])
{
    int i;

    for (i = 0; (i < maxlen - 1) && (src[i] != '\0'); i++)
    {
        if ((src[i] == ' ') || (src[i] == '\t'))
        {
            dst[i] = '_';
        }
        else
        {
            dst[i] = src[i];
        }
    }
    dst[i] = '\0';

    if (i == maxlen)
    {
        fprintf(stderr, "String '%s' truncated to '%s'\n", src, dst);
    }
}

static void decode_string(int maxlen, char dst[], const char src[])
{
    int i;

    for (i = 0; (i < maxlen - 1) && (src[i] != '\0'); i++)
    {
        if (src[i] == '_')
        {
            dst[i] = ' ';
        }
        else
        {
            dst[i] = src[i];
        }
    }
    dst[i] = '\0';

    if (i == maxlen)
    {
        fprintf(stderr, "String '%s' truncated to '%s'\n", src, dst);
    }
}

static gmx_bool do_ascwrite(t_fileio *fio, const void *item, int nitem, int eio,
                            const char *desc, const char *srcfile, int line)
{
    int            i;
    int            res = 0, *iptr;
    real          *ptr;
    char           strbuf[256];
    char           buf[GMX_FIO_BUFLEN];
    unsigned char *ucptr;
    FILE          *fp = fio->fp;

    gmx_fio_check_nitem(eio, nitem, srcfile, line);
    switch (eio)
    {
        case eioREAL:
        case eioFLOAT:
        case eioDOUBLE:
            res = fprintf(fp, "%18.10e%s\n", *((real *) item),
                          gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioINT:
            res = fprintf(fp, "%18d%s\n", *((int *) item), gmx_fio_dbgstr(fio,
                                                                          desc,
                                                                          buf));
            break;
        case eioINT64:
            sprintf(strbuf, "%s%s%s", "%", GMX_PRId64, "\n");
            res = fprintf(fp, strbuf, *((gmx_int64_t *) item),
                          gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioUCHAR:
            res = fprintf(fp, "%4d%s\n", *((unsigned char *) item),
                          gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioNUCHAR:
            ucptr = (unsigned char *) item;
            for (i = 0; (i < nitem); i++)
            {
                res = fprintf(fp, "%4d", (int) ucptr[i]);
            }
            fprintf(fio->fp, "%s\n", gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioUSHORT:
            res = fprintf(fp, "%18d%s\n", *((unsigned short *) item),
                          gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioRVEC:
            ptr = (real *) item;
            res = fprintf(fp, "%18.10e%18.10e%18.10e%s\n", ptr[XX],
                          ptr[YY], ptr[ZZ], gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioNRVEC:
            for (i = 0; (i < nitem); i++)
            {
                ptr = ((rvec *) item)[i];
                res = fprintf(fp, "%18.10e%18.10e%18.10e%s\n", ptr[XX],
                              ptr[YY], ptr[ZZ], gmx_fio_dbgstr(fio, desc, buf));
            }
            break;
        case eioIVEC:
            iptr = (int *) item;
            res  = fprintf(fp, "%18d%18d%18d%s\n", iptr[XX], iptr[YY],
                           iptr[ZZ], gmx_fio_dbgstr(fio, desc, buf));
            break;
        case eioSTRING:
            encode_string(256, strbuf, (char *) item);
            res = fprintf(fp, "%-18s%s\n", strbuf, gmx_fio_dbgstr(fio, desc, buf));
            break;
        default:
            gmx_fio_fe(fio, eio, desc, srcfile, line);
    }
    if ((res <= 0) && fio->bDebug)
    {
        fprintf(stderr,
                "Error writing %s %s to file %s (source %s, line %d)\n",
                eioNames[eio], desc, fio->fn, srcfile, line);
    }

    return (res > 0);
}


static char *next_item(FILE *fp, char *buf, int buflen)
{
    int      rd;
    gmx_bool in_comment = FALSE;
    gmx_bool in_token   = FALSE;
    int      i          = 0;
    /* This routine reads strings from the file fp, strips comment
     * and buffers. For thread-safety reasons, It reads through getc()  */

    rd = getc(fp);
    if (rd == EOF)
    {
        gmx_file("End of file");
    }
    do
    {
        if (in_comment)
        {
            if (rd == '\n')
            {
                in_comment = FALSE;
            }
        }
        else if (in_token)
        {
            if (isspace(rd) || rd == ';')
            {
                break;
            }
            buf[i++] = (char) rd;
        }
        else
        {
            if (!isspace(rd))
            {
                if (rd == ';')
                {
                    in_comment = TRUE;
                }
                else
                {
                    in_token = TRUE;
                    buf[i++] = (char) (rd);
                }
            }
        }
        if (i >= buflen - 2)
        {
            break;
        }
    }
    while ((rd = getc(fp)) != EOF);

    fprintf(stderr, "WARNING, ftpASC file type not tested!\n");

    buf[i] = 0;

    return buf;
}

static gmx_bool do_ascread(t_fileio *fio, void *item, int nitem, int eio,
                           const char *desc, const char *srcfile, int line)
{
    FILE           *fp = fio->fp;
    int             i, m, res = 0, *iptr, ix;
    gmx_int64_t     s;
    double          d, x;
    char            c;
    real           *ptr;
    unsigned char  *ucptr;
    char           *cptr;
#define NEXT_ITEM_BUF_LEN 128
    char            ni_buf[NEXT_ITEM_BUF_LEN];

    gmx_fio_check_nitem(eio, nitem, srcfile, line);
    switch (eio)
    {
        case eioREAL:
        case eioFLOAT:
        case eioDOUBLE:
            res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%lf", &d);
            if (item)
            {
                *((real *) item) = d;
            }
            break;
        case eioINT:
            res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%d", &i);
            if (item)
            {
                *((int *) item) = i;
            }
            break;
        case eioINT64:
            res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN),
                         "%"GMX_SCNd64, &s);
            if (item)
            {
                *((gmx_int64_t *) item) = s;
            }
            break;
        case eioUCHAR:
            res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%c", &c);
            if (item)
            {
                *((unsigned char *) item) = (unsigned char)c;
            }
            break;
        case eioNUCHAR:
            ucptr = (unsigned char *) item;
            for (i = 0; (i < nitem); i++)
            {
                res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%d", &ix);
                if (item)
                {
                    ucptr[i] = ix;
                }
            }
            break;
        case eioUSHORT:
            res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%d", &i);
            if (item)
            {
                *((unsigned short *) item) = i;
            }
            break;
        case eioRVEC:
            ptr = (real *) item;
            for (m = 0; (m < DIM); m++)
            {
                res    = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%lf\n", &x);
                ptr[m] = x;
            }
            break;
        case eioNRVEC:
            assert(item);
            for (i = 0; (i < nitem); i++)
            {
                ptr = ((rvec *) item)[i];
                for (m = 0; (m < DIM); m++)
                {
                    res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%lf\n",
                                 &x);
                    ptr[m] = x;
                }
            }
            break;
        case eioIVEC:
            iptr = (int *) item;
            for (m = 0; (m < DIM); m++)
            {
                res = sscanf(next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN), "%d\n", &ix);
                if (item)
                {
                    iptr[m] = ix;
                }
            }
            break;
        case eioSTRING:
            cptr = next_item(fp, ni_buf, NEXT_ITEM_BUF_LEN);
            if (item)
            {
                decode_string(strlen(cptr) + 1, (char *) item, cptr);
                /* res = sscanf(cptr,"%s",(char *)item);*/
                res = 1;
            }
            break;
        default:
            gmx_fio_fe(fio, eio, desc, srcfile, line);
    }

    if ((res <= 0) && fio->bDebug)
    {
        fprintf(stderr,
                "Error reading %s %s from file %s (source %s, line %d)\n",
                eioNames[eio], desc, fio->fn, srcfile, line);
    }
    return (res > 0);
}
