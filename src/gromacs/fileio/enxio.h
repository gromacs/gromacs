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
#ifndef GMX_FILEIO_ENXIO_H
#define GMX_FILEIO_ENXIO_H

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xdr_datatype.h"
#include "gromacs/legacyheaders/types/energy.h"
#include "gromacs/legacyheaders/types/inputrec.h"
#include "gromacs/legacyheaders/types/state.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_groups_t;

/**************************************************************
 * These are the base datatypes + functions for reading and
 * writing energy files (.edr). They are either called directly
 * (as in the processing tools), or indirectly through mdebin.c
 * during mdrun.
 *
 * The routines in the corresponding c-file enxio.c
 * are based on the lower level routines in gmxfio.c.
 * The file pointer returned from open_enx
 * can also be used with the routines in gmxfio.h
 *
 **************************************************************/

typedef struct {
    char *name;
    char *unit;
} gmx_enxnm_t;

/*
 * Index for the IDs of additional blocks in the energy file.
 * Blocks can be added without sacrificing backward and forward
 * compatibility of the energy files.
 *
 * For backward compatibility, the order of these should not be changed.
 */
enum {
    enxOR,     /* Time and ensemble averaged data for orientation restraints */
    enxORI,    /* Instantaneous data for orientation restraints              */
    enxORT,    /* Order tensor(s) for orientation restraints                 */
    enxDISRE,  /* Distance restraint blocks                                  */

    enxDHCOLL, /* Data about the free energy blocks in this frame.           */
    enxDHHIST, /* BAR histogram                                              */
    enxDH,     /* BAR raw delta H data                                       */

    enxNR      /* Total number of extra blocks in the current code,
                * note that the enxio code can read files written by
                * future code which contain more blocks.
                */
};

/* names for the above enum */
extern const char *enx_block_id_name[];


/* the subblocks that are contained in energy file blocks. Each of these
   has a number of values of a single data type in a .edr file. */
typedef struct
{
    int          nr;        /* number of items in subblock */
    xdr_datatype type;      /* the block type */

    /* the values: pointers for each type */
    float*             fval;
    double*            dval;
    int*               ival;
    gmx_int64_t*       lval;
    unsigned char*     cval;
    char**             sval;

    /* the allocated sizes, defined separately.
       (nonzero sizes can be free()d later): */
    int fval_alloc;
    int dval_alloc;
    int ival_alloc;
    int lval_alloc;
    int cval_alloc;
    int sval_alloc;
} t_enxsubblock;


/* the energy file blocks. Each block contains a number of sub-blocks
   of a single type that contain the actual data. */
typedef struct t_enxblock{
    int            id;         /* block id, from the enx enums above */
    int            nsub;       /* number of subblocks */
    t_enxsubblock *sub;        /* the subblocks */
    int            nsub_alloc; /* number of allocated subblocks */
} t_enxblock;


/* The frames that are read/written */
typedef struct {
    double          t;            /* Timestamp of this frame	                     */
    gmx_int64_t     step;         /* MD step	                             */
    gmx_int64_t     nsteps;       /* The number of steps between frames            */
    double          dt;           /* The MD time step                              */
    int             nsum;         /* The number of terms for the sums in ener      */
    int             nre;          /* Number of energies			     */
    int             e_size;       /* Size (in bytes) of energies		     */
    int             e_alloc;      /* Allocated size (in elements) of ener          */
    t_energy       *ener;         /* The energies                                  */
    int             nblock;       /* Number of following energy blocks             */
    t_enxblock     *block;        /* The blocks                                    */
    int             nblock_alloc; /* The number of blocks allocated                */
} t_enxframe;

/* file handle */
typedef struct ener_file *ener_file_t;

/*
 * An energy file is read like this:
 *
 * ener_file_t fp;
 * t_enxframe *fr;
 *
 * fp = open_enx(...);
 * do_enxnms(fp,...);
 * snew(fr,1);
 * while (do_enx(fp,fr)) {
 * ...
 * }
 * free_enxframe(fr);
 * sfree(fr);
 */
/* New energy reading and writing interface */


/* initialize a pre-allocated frame */
void init_enxframe(t_enxframe *ef);
/* delete a frame's memory (except the ef itself) */
void free_enxframe(t_enxframe *ef);


ener_file_t open_enx(const char *fn, const char *mode);

t_fileio *enx_file_pointer(const ener_file_t ef);

void close_enx(ener_file_t ef);

void do_enxnms(ener_file_t ef, int *nre, gmx_enxnm_t **enms);

void free_enxnms(int n, gmx_enxnm_t *nms);
/* Frees nms and all strings in it */

gmx_bool do_enx(ener_file_t ef, t_enxframe *fr);
/* Reads enx_frames, memory in fr is (re)allocated if necessary */

void get_enx_state(const char *fn, real t,
                   struct gmx_groups_t *groups, t_inputrec *ir,
                   t_state *state);
/*
 * Reads state variables from enx file fn at time t.
 * atoms and ir are required for determining which things must be read.
 * Currently pcoupl and tcoupl state are read from enx.
 */


/* block funtions */

/* allocate n blocks to a frame (if neccesary). Don't touch existing blocks */
void add_blocks_enxframe(t_enxframe *ef, int n);

/* find a block by id number; if prev!=NULL, it searches from
   that block's next block.
   Returns NULL if no block is found with the given id. */
t_enxblock *find_block_id_enxframe(t_enxframe *ef, int id, t_enxblock *prev);


/* allocate n subblocks to a block (if neccesary). Don't touch existing
   subbblocks. */
void add_subblocks_enxblock(t_enxblock *eb, int n);



#ifdef __cplusplus
}
#endif

#endif
