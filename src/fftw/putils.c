
/*
 * Copyright (c) 1997,1998 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
 * putils.c -- plan utilities shared by planner.c and rplanner.c
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw-int.h>
#include <stdlib.h>
#include <stdio.h>

int fftw_node_cnt = 0;
int fftw_plan_cnt = 0;

/*
 * These two constants are used for the FFTW_ESTIMATE flag to help
 * create a heuristic plan.  They don't affect FFTW_MEASURE.
 */
#define NOTW_OPTIMAL_SIZE 32
#define TWIDDLE_OPTIMAL_SIZE 12

#define IS_POWER_OF_TWO(n) ((n & (n - 1)) == 0)

/* constructors --- I wish I had ML */
fftw_plan_node *fftw_make_node(void)
{
     fftw_plan_node *p = (fftw_plan_node *)
     fftw_malloc(sizeof(fftw_plan_node));
     p->refcnt = 0;
     fftw_node_cnt++;
     return p;
}

void fftw_use_node(fftw_plan_node *p)
{
     ++p->refcnt;
}

fftw_plan_node *fftw_make_node_notw(int size, const fftw_codelet_desc *config)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = config->type;
     p->nodeu.notw.size = size;
     p->nodeu.notw.codelet = (fftw_notw_codelet *) config->codelet;
     p->nodeu.notw.codelet_desc = config;
     return p;
}

fftw_plan_node *fftw_make_node_real2hc(int size,
				       const fftw_codelet_desc *config)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = config->type;
     p->nodeu.real2hc.size = size;
     p->nodeu.real2hc.codelet = (fftw_real2hc_codelet *) config->codelet;
     p->nodeu.real2hc.codelet_desc = config;
     return p;
}

fftw_plan_node *fftw_make_node_hc2real(int size,
				       const fftw_codelet_desc *config)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = config->type;
     p->nodeu.hc2real.size = size;
     p->nodeu.hc2real.codelet = (fftw_hc2real_codelet *) config->codelet;
     p->nodeu.hc2real.codelet_desc = config;
     return p;
}

fftw_plan_node *fftw_make_node_twiddle(int n,
				       const fftw_codelet_desc *config,
				       fftw_plan_node *recurse,
				       int flags)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = config->type;
     p->nodeu.twiddle.size = config->size;
     p->nodeu.twiddle.codelet = (fftw_twiddle_codelet *) config->codelet;
     p->nodeu.twiddle.recurse = recurse;
     p->nodeu.twiddle.codelet_desc = config;
     fftw_use_node(recurse);
     if (flags & FFTW_MEASURE)
	  p->nodeu.twiddle.tw = fftw_create_twiddle(n, config);
     else
	  p->nodeu.twiddle.tw = 0;
     return p;
}

fftw_plan_node *fftw_make_node_hc2hc(int n, fftw_direction dir,
				     const fftw_codelet_desc *config,
				     fftw_plan_node *recurse,
				     int flags)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = config->type;
     p->nodeu.hc2hc.size = config->size;
     p->nodeu.hc2hc.dir = dir;
     p->nodeu.hc2hc.codelet = (fftw_hc2hc_codelet *) config->codelet;
     p->nodeu.hc2hc.recurse = recurse;
     p->nodeu.hc2hc.codelet_desc = config;
     fftw_use_node(recurse);
     if (flags & FFTW_MEASURE)
	  p->nodeu.hc2hc.tw = fftw_create_twiddle(n, config);
     else
	  p->nodeu.hc2hc.tw = 0;
     return p;
}

fftw_plan_node *fftw_make_node_generic(int n, int size,
				       fftw_generic_codelet *codelet,
				       fftw_plan_node *recurse,
				       int flags)
{
     fftw_plan_node *p = fftw_make_node();

     p->type = FFTW_GENERIC;
     p->nodeu.generic.size = size;
     p->nodeu.generic.codelet = codelet;
     p->nodeu.generic.recurse = recurse;
     fftw_use_node(recurse);

     if (flags & FFTW_MEASURE)
	  p->nodeu.generic.tw = fftw_create_twiddle(n,
					  (const fftw_codelet_desc *) 0);
     else
	  p->nodeu.generic.tw = 0;
     return p;
}

fftw_plan_node *fftw_make_node_rgeneric(int n, int size,
					fftw_direction dir,
					fftw_rgeneric_codelet *codelet,
					fftw_plan_node *recurse,
					int flags)
{
     fftw_plan_node *p = fftw_make_node();

     if (size % 2 == 0 || (n / size) % 2 == 0)
	  fftw_die("invalid size for rgeneric codelet\n");

     p->type = FFTW_RGENERIC;
     p->nodeu.rgeneric.size = size;
     p->nodeu.rgeneric.dir = dir;
     p->nodeu.rgeneric.codelet = codelet;
     p->nodeu.rgeneric.recurse = recurse;
     fftw_use_node(recurse);

     if (flags & FFTW_MEASURE)
	  p->nodeu.rgeneric.tw = fftw_create_twiddle(n,
					  (const fftw_codelet_desc *) 0);
     else
	  p->nodeu.rgeneric.tw = 0;
     return p;
}

/* 
 * Note that these two Rader-related things must go here, rather than
 * in rader.c, in order that putils.c (and rplanner.c) won't depend
 * upon rader.c. 
 */

fftw_rader_data *fftw_rader_top = NULL;

static void fftw_destroy_rader(fftw_rader_data * d)
{
     if (d) {
	  d->refcount--;
	  if (d->refcount <= 0) {
	       fftw_rader_data *cur = fftw_rader_top, *prev = NULL;

	       while (cur && cur != d) {
		    prev = cur;
		    cur = cur->next;
	       }
	       if (!cur)
		    fftw_die("invalid Rader data pointer\n");

	       if (prev)
		    prev->next = d->next;
	       else
		    fftw_rader_top = d->next;

	       fftw_destroy_plan_internal(d->plan);
	       fftw_free(d->omega);
	       fftw_free(d->cdesc);
	       fftw_free(d);
	  }
     }
}

static void destroy_tree(fftw_plan_node *p)
{
     if (p) {
	  --p->refcnt;
	  if (p->refcnt == 0) {
	       switch (p->type) {
		   case FFTW_NOTW:
		   case FFTW_REAL2HC:
		   case FFTW_HC2REAL:
			break;

		   case FFTW_TWIDDLE:
			if (p->nodeu.twiddle.tw)
			     fftw_destroy_twiddle(p->nodeu.twiddle.tw);
			destroy_tree(p->nodeu.twiddle.recurse);
			break;

		   case FFTW_HC2HC:
			if (p->nodeu.hc2hc.tw)
			     fftw_destroy_twiddle(p->nodeu.hc2hc.tw);
			destroy_tree(p->nodeu.hc2hc.recurse);
			break;

		   case FFTW_GENERIC:
			if (p->nodeu.generic.tw)
			     fftw_destroy_twiddle(p->nodeu.generic.tw);
			destroy_tree(p->nodeu.generic.recurse);
			break;

		   case FFTW_RADER:
			if (p->nodeu.rader.tw)
			     fftw_destroy_twiddle(p->nodeu.rader.tw);
			if (p->nodeu.rader.rader_data)
			     fftw_destroy_rader(p->nodeu.rader.rader_data);
			destroy_tree(p->nodeu.rader.recurse);
			break;

		   case FFTW_RGENERIC:
			if (p->nodeu.rgeneric.tw)
			     fftw_destroy_twiddle(p->nodeu.rgeneric.tw);
			destroy_tree(p->nodeu.rgeneric.recurse);
			break;
	       }

	       fftw_free(p);
	       fftw_node_cnt--;
	  }
     }
}

/* create a plan with twiddle factors, and other bells and whistles */
fftw_plan fftw_make_plan(int n, fftw_direction dir,
			 fftw_plan_node *root, int flags,
			 enum fftw_node_type wisdom_type,
			 int wisdom_signature)
{
     fftw_plan p = (fftw_plan) fftw_malloc(sizeof(struct fftw_plan_struct));

     p->n = n;
     p->dir = dir;
     p->flags = flags;
     fftw_use_node(root);
     p->root = root;
     p->cost = 0.0;
     p->wisdom_type = wisdom_type;
     p->wisdom_signature = wisdom_signature;
     p->next = (fftw_plan) 0;
     p->refcnt = 0;
     fftw_plan_cnt++;
     return p;
}

/*
 * complete with twiddle factors (because nodes don't have
 * them when FFTW_ESTIMATE is set)
 */
void fftw_complete_twiddle(fftw_plan_node *p, int n)
{
     int r;
     switch (p->type) {
	 case FFTW_NOTW:
	 case FFTW_REAL2HC:
	 case FFTW_HC2REAL:
	      break;

	 case FFTW_TWIDDLE:
	      r = p->nodeu.twiddle.size;
	      if (!p->nodeu.twiddle.tw)
		   p->nodeu.twiddle.tw =
		       fftw_create_twiddle(n, p->nodeu.twiddle.codelet_desc);
	      fftw_complete_twiddle(p->nodeu.twiddle.recurse, n / r);
	      break;

	 case FFTW_HC2HC:
	      r = p->nodeu.hc2hc.size;
	      if (!p->nodeu.hc2hc.tw)
		   p->nodeu.hc2hc.tw =
		       fftw_create_twiddle(n, p->nodeu.hc2hc.codelet_desc);
	      fftw_complete_twiddle(p->nodeu.hc2hc.recurse, n / r);
	      break;

	 case FFTW_GENERIC:
	      r = p->nodeu.generic.size;
	      if (!p->nodeu.generic.tw)
		   p->nodeu.generic.tw =
		       fftw_create_twiddle(n, (const fftw_codelet_desc *) 0);
	      fftw_complete_twiddle(p->nodeu.generic.recurse, n / r);
	      break;

	 case FFTW_RADER:
	      r = p->nodeu.rader.size;
	      if (!p->nodeu.rader.tw)
		   p->nodeu.rader.tw =
		       fftw_create_twiddle(n, p->nodeu.rader.rader_data->cdesc);
	      fftw_complete_twiddle(p->nodeu.rader.recurse, n / r);
	      break;

	 case FFTW_RGENERIC:
	      r = p->nodeu.rgeneric.size;
	      if (!p->nodeu.rgeneric.tw)
		   p->nodeu.rgeneric.tw =
		       fftw_create_twiddle(n, (const fftw_codelet_desc *) 0);
	      fftw_complete_twiddle(p->nodeu.rgeneric.recurse, n / r);
	      break;

     }
}

void fftw_use_plan(fftw_plan p)
{
     ++p->refcnt;
}

void fftw_destroy_plan_internal(fftw_plan p)
{
     --p->refcnt;

     if (p->refcnt == 0) {
	  destroy_tree(p->root);
	  fftw_plan_cnt--;
	  fftw_free(p);
     }
}

/* end of constructors */

/* management of plan tables */
void fftw_make_empty_table(fftw_plan *table)
{
     *table = (fftw_plan) 0;
}

void fftw_insert(fftw_plan *table, fftw_plan this_plan, int n)
{
     fftw_use_plan(this_plan);
     this_plan->n = n;
     this_plan->next = *table;
     *table = this_plan;
}

fftw_plan fftw_lookup(fftw_plan *table, int n, int flags)
{
     fftw_plan p;

     for (p = *table; p &&
	  ((p->n != n) || (p->flags != flags)); p = p->next);

     return p;
}

void fftw_destroy_table(fftw_plan *table)
{
     fftw_plan p, q;

     for (p = *table; p; p = q) {
	  q = p->next;
	  fftw_destroy_plan_internal(p);
     }
}

double fftw_estimate_node(fftw_plan_node *p)
{
     int k;

     switch (p->type) {
	 case FFTW_NOTW:
	      k = p->nodeu.notw.size;
	      goto common1;

	 case FFTW_REAL2HC:
	      k = p->nodeu.real2hc.size;
	      goto common1;

	 case FFTW_HC2REAL:
	      k = p->nodeu.hc2real.size;
	    common1:
	      return 1.0 + 0.1 * (k - NOTW_OPTIMAL_SIZE) *
		  (k - NOTW_OPTIMAL_SIZE);

	 case FFTW_TWIDDLE:
	      k = p->nodeu.twiddle.size;
	      return 1.0 + 0.1 * (k - TWIDDLE_OPTIMAL_SIZE) *
		  (k - TWIDDLE_OPTIMAL_SIZE)
		  + fftw_estimate_node(p->nodeu.twiddle.recurse);

	 case FFTW_HC2HC:
	      k = p->nodeu.hc2hc.size;
	      return 1.0 + 0.1 * (k - TWIDDLE_OPTIMAL_SIZE) *
		  (k - TWIDDLE_OPTIMAL_SIZE)
		  + fftw_estimate_node(p->nodeu.hc2hc.recurse);

	 case FFTW_GENERIC:
	      k = p->nodeu.generic.size;
	      return 10.0 + k * k
		  + fftw_estimate_node(p->nodeu.generic.recurse);

	 case FFTW_RADER:
	      k = p->nodeu.rader.size;
	      return 10.0 + 10 * k
		  + fftw_estimate_node(p->nodeu.rader.recurse);

	 case FFTW_RGENERIC:
	      k = p->nodeu.rgeneric.size;
	      return 10.0 + k * k
		  + fftw_estimate_node(p->nodeu.rgeneric.recurse);
     }
     return 1.0E20;
}

/* pick the better of two plans and destroy the other one. */
fftw_plan fftw_pick_better(fftw_plan p1, fftw_plan p2)
{
     if (!p1)
	  return p2;

     if (!p2)
	  return p1;

     if (p1->cost > p2->cost) {
	  fftw_destroy_plan_internal(p1);
	  return p2;
     } else {
	  fftw_destroy_plan_internal(p2);
	  return p1;
     }
}

/* find the smallest prime factor of n */
int fftw_factor(int n)
{
     int r;

     /* try 2 */
     if ((n & 1) == 0)
	  return 2;

     /* try odd numbers up to sqrt(n) */
     for (r = 3; r * r <= n; r += 2)
	  if (n % r == 0)
	       return r;

     /* n is prime */
     return n;
}

static void print_node(FILE *f, fftw_plan_node *p, int indent)
{
     if (p) {
	  switch (p->type) {
	      case FFTW_NOTW:
		   fprintf(f, "%*sFFTW_NOTW %d\n", indent, "",
			   p->nodeu.notw.size);
		   break;
	      case FFTW_REAL2HC:
		   fprintf(f, "%*sFFTW_REAL2HC %d\n", indent, "",
			   p->nodeu.real2hc.size);
		   break;
	      case FFTW_HC2REAL:
		   fprintf(f, "%*sFFTW_HC2REAL %d\n", indent, "",
			   p->nodeu.hc2real.size);
		   break;
	      case FFTW_TWIDDLE:
		   fprintf(f, "%*sFFTW_TWIDDLE %d\n", indent, "",
			   p->nodeu.twiddle.size);
		   print_node(f, p->nodeu.twiddle.recurse, indent);
		   break;
	      case FFTW_HC2HC:
		   fprintf(f, "%*sFFTW_HC2HC %d\n", indent, "",
			   p->nodeu.hc2hc.size);
		   print_node(f, p->nodeu.hc2hc.recurse, indent);
		   break;
	      case FFTW_GENERIC:
		   fprintf(f, "%*sFFTW_GENERIC %d\n", indent, "",
			   p->nodeu.generic.size);
		   print_node(f, p->nodeu.generic.recurse, indent);
		   break;
	      case FFTW_RADER:
		   fprintf(f, "%*sFFTW_RADER %d\n", indent, "",
			   p->nodeu.rader.size);

		   fprintf(f, "%*splan for size %d convolution:\n",
			   indent + 4, "", p->nodeu.rader.size - 1);
		   print_node(f, p->nodeu.rader.rader_data->plan->root,
			      indent + 6);

		   print_node(f, p->nodeu.rader.recurse, indent);
		   break;
	      case FFTW_RGENERIC:
		   fprintf(f, "%*sFFTW_RGENERIC %d\n", indent, "",
			   p->nodeu.rgeneric.size);
		   print_node(f, p->nodeu.rgeneric.recurse, indent);
		   break;
	  }
     }
}

void fftw_fprint_plan(FILE *f, fftw_plan p)
{

     fprintf(f, "plan: (cost = %e)\n", p->cost);
     print_node(f, p->root, 0);
}

void fftw_print_plan(fftw_plan p)
{
     fftw_fprint_plan(stdout, p);
}

size_t fftw_sizeof_fftw_real(void)
{
     return(sizeof(fftw_real));
}
