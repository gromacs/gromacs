/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _sema_h
#define _sema_h

static char *SRCID_sema_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) sema.h 1.7 11/23/92"
#endif /* HAVE_IDENT */

#include <sys/types.h>

#ifdef HAVE_SEMA
#include <sys/ipc.h>
#include <sys/sem.h>
#else
/* Dummy declarations... */
#ifndef _sol_
typedef int key_t;
struct sembuf { int i; };
#endif
#endif

typedef struct
{
  int sem_id;
  key_t sem_key;
  char *sem_name;
  struct sembuf sem_ops[1];
} t_sema;

extern void create_sema(t_sema *sema,char *sem_name,key_t key,int count);
     /*
      * Use this routine to create and to initialise the semaphore, before
      * any other process accesses it. The name is a simple debugging tool
      * to identify the semaphore (not needed in the OS routines). The value
      * of key is needed by a process to attach to the specific semaphore.
      */
     
extern void destroy_sema(t_sema *sema);
     /*
      * After a semaphore is not needed anymore, a call to destroy_sema
      * removes it. This routine needs to be called because semaphores 
      * are not automatically disposed of when they are not used anymore.
      */

extern void ini_sema(t_sema *sema,char *sem_name,key_t key);
     /*
      * Use this routine to initialise the semaphore, after it has been 
      * created by a call to create_sema. A call to create_sema implies 
      * the initialisation of the semaphore, so it is not needed to call
      * ini_sema after a call of create_sema. The name is a simple debugging 
      * tool to identify the semaphore (not needed in the OS routines). The 
      * value of key is needed by a process to attach to the specific 
      * semaphore.
      */

extern void p_sema(t_sema *sema);
     /*
      * Applies the P-operation on a semaphore i.e. decrement the value.
      * If the value is zero, suspend execution until it is incremented.
      */

extern void v_sema(t_sema *sema);
     /*
      * Applies the V-operation on a semaphore i.e. increment the value.
      * A process blocked by executing p_sema() may proceed.
      */

#endif	/* _sema_h */
