#include "typedefs.h"
#include "smalloc.h"
#include "assert.h"
#include "symtab.h"

void init_block(t_block *block)
{
  int i;

  block->nr    = 0;
  block->nra   = 0;
  snew(block->index,1);
  block->index[0] = 0;
  block->a     = NULL;
  for(i=0; (i<MAXPROC); i++)
    block->multinr[i]=0;
}

void init_atom (t_atoms *at)
{
  int i;
  
  init_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  at->ngrpname = 0;
  at->atom     = NULL;
  at->resname  = NULL;
  at->atomname = NULL;
  at->grpname  = NULL;
  for(i=0; (i<egcNR); i++) {
    at->grps[i].nr=0;
    at->grps[i].nm_ind=NULL;
  }
}

void init_top (t_topology *top)
{
  int i;
  
  top->name = NULL;
  init_atom (&(top->atoms));
  for (i=0; (i<ebNR); i++)
    init_block(&(top->blocks[i]));
}

void init_inputrec(t_inputrec *ir)
{
  memset(ir,0,(size_t)sizeof(*ir));
}

void done_block(t_block *block)
{
  block->nr    = 0;
  block->nra   = 0;
  sfree(block->index);
  sfree(block->a);
}

void done_atom (t_atoms *at)
{
  done_block(&(at->excl));
  at->nr       = 0;
  at->nres     = 0;
  sfree(at->atom);
  sfree(at->resname);
  sfree(at->atomname);
}

void done_symtab(t_symtab *symtab)
{
  int i;
  t_symbuf *symbuf,*freeptr;
  
  close_symtab(symtab);
  symbuf=symtab->symbuf;
  while (symbuf!=NULL) {
    for (i=0; (i<symbuf->bufsize)&&(i<symtab->nr); i++)
      sfree(symbuf->buf[i]);
    symtab->nr-=i;
    sfree(symbuf->buf);
    freeptr=symbuf;
    symbuf=symbuf->next;
    sfree(freeptr);
  }
  symtab->symbuf=NULL;
  assert(symtab->nr==0);
}

void done_top(t_topology *top)
{
  int i;
  
  done_atom (&(top->atoms));
  done_symtab(&(top->symtab));
  for (i=0; (i<ebNR); i++)
    done_block(&(top->blocks[i]));
}

void done_inputrec(t_inputrec *ir)
{
  int m;
  
  for(m=0; (m<DIM); m++) {
    if (ir->ex[m].a)   sfree(ir->ex[m].a);
    if (ir->ex[m].phi) sfree(ir->ex[m].phi);
    if (ir->et[m].a)   sfree(ir->et[m].a);
    if (ir->et[m].phi) sfree(ir->et[m].phi);
  }
  if (ir->opts.nrdf)    sfree(ir->opts.nrdf);
  if (ir->opts.ref_t)   sfree(ir->opts.ref_t);
  if (ir->opts.tau_t)   sfree(ir->opts.tau_t);
  if (ir->opts.acc)     sfree(ir->opts.acc);
  if (ir->opts.nFreeze) sfree(ir->opts.nFreeze);
}

