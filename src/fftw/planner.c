/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

/*
 * planner.c -- find the optimal plan
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw.h>
#include <stdlib.h>
#include <stdio.h>

int fftw_node_cnt = 0;
int fftw_plan_cnt = 0;

#define NOTW_OPTIMAL_SIZE 32
#define TWIDDLE_OPTIMAL_SIZE 12

/* wisdom prototypes */
extern int fftw_wisdom_lookup(int n, int flags, fftw_direction dir,
			    enum fftw_node_type *type,
			    int *signature, int replace_p);
extern void fftw_wisdom_add(int n, int flags, fftw_direction dir,
			  enum fftw_node_type type,
			  int signature);

/* constructors --- I wish I had ML */
static fftw_plan_node *make_node(void)
{
     fftw_plan_node *p = (fftw_plan_node *)
     fftw_malloc(sizeof(fftw_plan_node));
     p->refcnt = 0;
     fftw_node_cnt++;
     return p;
}

static void use_node(fftw_plan_node *p)
{
     ++p->refcnt;
}

static fftw_plan_node *make_node_notw(int size, notw_codelet *codelet)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_NOTW;
     p->nodeu.notw.size = size;
     p->nodeu.notw.codelet = codelet;
     return p;
}

static fftw_plan_node *make_node_twiddle(int n, int size, twiddle_codelet *codelet,
					 fftw_plan_node *recurse,
					 int flags)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_TWIDDLE;
     p->nodeu.twiddle.size = size;
     p->nodeu.twiddle.codelet = codelet;
     p->nodeu.twiddle.recurse = recurse;
     use_node(recurse);
     if (flags & FFTW_MEASURE)
	  p->nodeu.twiddle.tw = fftw_create_twiddle(n, size, n / size);
     else
	  p->nodeu.twiddle.tw = 0;
     return p;
}

static fftw_plan_node *make_node_generic(int n, int size,
					 generic_codelet *codelet,
					 fftw_plan_node *recurse,
					 int flags)
{
     fftw_plan_node *p = make_node();

     p->type = FFTW_GENERIC;
     p->nodeu.generic.size = size;
     p->nodeu.generic.codelet = codelet;
     p->nodeu.generic.recurse = recurse;
     use_node(recurse);

     if (flags & FFTW_MEASURE)
	  p->nodeu.generic.tw = fftw_create_twiddle(n, 2, n);
     else
	  p->nodeu.generic.tw = 0;
     return p;
}

static void destroy_tree(fftw_plan_node *p)
{
     if (p) {
	  --p->refcnt;
	  if (p->refcnt == 0) {
	       switch (p->type) {
		   case FFTW_NOTW:
			break;

		   case FFTW_TWIDDLE:
			if (p->nodeu.twiddle.tw)
			     fftw_destroy_twiddle(p->nodeu.twiddle.tw);
			destroy_tree(p->nodeu.twiddle.recurse);
			break;

		   case FFTW_GENERIC:
			if (p->nodeu.generic.tw)
			     fftw_destroy_twiddle(p->nodeu.generic.tw);
			destroy_tree(p->nodeu.generic.recurse);
			break;
	       }

	       fftw_free(p);
	       fftw_node_cnt--;
	  }
     }
}

/* create a plan with twiddle factors, and other bells and whistles */
static fftw_plan make_plan(int n, fftw_direction dir,
			   fftw_plan_node *root, int flags,
			   enum fftw_node_type wisdom_type,
			   int wisdom_signature)
{
     fftw_plan p = (fftw_plan) fftw_malloc(sizeof(struct fftw_plan_struct));

     p->n = n;
     p->dir = dir;
     p->flags = flags;
     use_node(root);
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
static void complete_twiddle(fftw_plan_node *p, int n)
{
     int r;
     switch (p->type) {
	 case FFTW_NOTW:
	      break;

	 case FFTW_TWIDDLE:
	      r = p->nodeu.twiddle.size;
	      if (!p->nodeu.twiddle.tw)
		   p->nodeu.twiddle.tw = fftw_create_twiddle(n, r, n / r);
	      complete_twiddle(p->nodeu.twiddle.recurse, n / r);
	      break;

	 case FFTW_GENERIC:
	      r = p->nodeu.generic.size;
	      if (!p->nodeu.generic.tw)
		   p->nodeu.generic.tw = fftw_create_twiddle(n, 2, n);
	      complete_twiddle(p->nodeu.generic.recurse, n / r);
	      break;
     }
}

static void use_plan(fftw_plan p)
{
     ++p->refcnt;
}

static void destroy_plan(fftw_plan p)
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
static void make_empty_table(fftw_plan *table)
{
     *table = (fftw_plan) 0;
}

static void insert(fftw_plan *table, fftw_plan this_plan, int n)
{
     use_plan(this_plan);
     this_plan->n = n;
     this_plan->next = *table;
     *table = this_plan;
}

static fftw_plan lookup(fftw_plan *table, int n, int flags)
{
     fftw_plan p;

     for (p = *table; p &&
	  ((p->n != n) || (p->flags != flags)); p = p->next);

     return p;
}

static void destroy_table(fftw_plan *table)
{
     fftw_plan p, q;

     for (p = *table; p; p = q) {
	  q = p->next;
	  destroy_plan(p);
     }
}

static double estimate_node(fftw_plan_node *p)
{
     int k;

     switch (p->type) {
	 case FFTW_NOTW:
	      k = p->nodeu.notw.size;
	      return 1.0 + 0.1 * (k - NOTW_OPTIMAL_SIZE) *
		  (k - NOTW_OPTIMAL_SIZE);

	 case FFTW_TWIDDLE:
	      k = p->nodeu.twiddle.size;
	      return 1.0 + 0.1 * (k - TWIDDLE_OPTIMAL_SIZE) *
		  (k - TWIDDLE_OPTIMAL_SIZE)
		  + estimate_node(p->nodeu.twiddle.recurse);

	 case FFTW_GENERIC:
	      k = p->nodeu.generic.size;
	      return 10.0 + k * k
		  + estimate_node(p->nodeu.generic.recurse);
     }
     return 1.0E20;
}

/* auxiliary functions */
static void compute_cost(fftw_plan plan)
{
     if (plan->flags & FFTW_MEASURE)
	  plan->cost = fftw_measure_runtime(plan);
     else {
	  double c;
	  c = plan->n * estimate_node(plan->root);
	  plan->cost = c;
     }
}

/* pick the better of two plans and destroy the other one. */
static fftw_plan pick_better(fftw_plan p1, fftw_plan p2)
{
     if (!p1)
	  return p2;

     if (!p2)
	  return p1;

     if (p1->cost > p2->cost) {
	  destroy_plan(p1);
	  return p2;
     } else {
	  destroy_plan(p2);
	  return p1;
     }
}

/* find the smallest prime factor of n */
static int factor(int n)
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

/* 
 * Some macrology for the planner.  If you have to write
 * the same line of code twice, there must be some bug.
 */
#define NOTW_ITERATOR(p, dir)                                \
      config_notw *p =                                       \
	  p = (dir == FFTW_FORWARD ?                         \
	       fftw_config_notw : fftwi_config_notw)

#define TWIDDLE_ITERATOR(p, dir)                             \
      config_twiddle *p =                                    \
	  p = (dir == FFTW_FORWARD ?                         \
	       fftw_config_twiddle : fftwi_config_twiddle);

#define FORALL_NOTW(p)             \
	 for (; p->size; ++p) 

#define FORALL_TWIDDLE(p)          \
	 for (; p->size; ++p) 

/******************************************
 *      Recursive planner                 *
 ******************************************/
fftw_plan planner(fftw_plan *table, int n, fftw_direction dir, int flags);

/*
 * the planner consists of two parts: one that tries to
 * use accumulated wisdom, and one that does not.
 * A small driver invokes both parts in sequence
 */

/* planner with wisdom: look up the codelet suggested by the wisdom */
fftw_plan planner_wisdom(fftw_plan *table, int n,
			 fftw_direction dir, int flags)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan_node *node;
     int have_wisdom;
     enum fftw_node_type wisdom_type;
     int wisdom_signature;

     /* see if we remember any wisdom for this case */
     have_wisdom = fftw_wisdom_lookup(n, flags, dir, 
				      &wisdom_type, &wisdom_signature, 0);

     if (!have_wisdom)
	  return best;

     if (wisdom_type == FFTW_NOTW) {
	  NOTW_ITERATOR(p, dir);
	       
	  FORALL_NOTW(p) {
	       /* see if wisdom applies */
	       if (wisdom_signature == p->signature &&
		   p->size == n) {
		    node = make_node_notw(n, p->codelet);
		    best = make_plan(n, dir, node, flags,
				     FFTW_NOTW, p->signature);
		    use_plan(best);
		    return best;
	       }
	  }
     }
	  
     if (wisdom_type == FFTW_TWIDDLE) {
	  TWIDDLE_ITERATOR(p, dir);
	       
	  FORALL_TWIDDLE(p) {
	       /* see if wisdom applies */
	       if (wisdom_signature == p->signature &&
		   (n % p->size) == 0) {
		    fftw_plan r = planner(table, n / p->size, dir, flags);
		    node = make_node_twiddle(n, p->size, p->codelet,
					     r->root, flags);
		    best = make_plan(n, dir, node, flags,
				     FFTW_TWIDDLE, p->signature);
		    use_plan(best);
		    destroy_plan(r);
		    return best;
	       }
	  }
     }

     /* 
      * BUG (or: TODO)  Can we have generic wisdom? This is probably
      * an academic question
      */

     return best;
}

/*
 * planner with no wisdom: try all combinations and pick
 * the best
 */
fftw_plan planner_normal(fftw_plan *table, int n, fftw_direction dir,
			   int flags)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan newplan;
     fftw_plan_node *node;

     /* see if we have any codelet that solves the problem */
     {
	  NOTW_ITERATOR(p, dir);
	       
	  FORALL_NOTW(p) {
	       if (p->size == n) {
		    node = make_node_notw(n, p->codelet);
		    newplan = make_plan(n, dir, node, flags,
					FFTW_NOTW, p->signature);
		    use_plan(newplan);
		    compute_cost(newplan);
		    best = pick_better(newplan, best);
	       }
	  }
     }

     /* Then, try all available twiddle codelets */
     {
	  TWIDDLE_ITERATOR(p, dir);
	       
	  FORALL_TWIDDLE(p) {
	       if ((n % p->size) == 0 &&
		   (!best || n != p->size)) {
		    fftw_plan r = planner(table, n / p->size, dir, flags);
		    node = make_node_twiddle(n, p->size, p->codelet,
					     r->root, flags);
		    newplan = make_plan(n, dir, node, flags,
					FFTW_TWIDDLE, p->signature);
		    use_plan(newplan);
		    destroy_plan(r);
		    compute_cost(newplan);
		    best = pick_better(newplan, best);
	       }
	  }
     }

     /* 
      * if no plan has been found so far, resort to generic codelets 
      */
     if (!best) {
	  generic_codelet *codelet = (dir == FFTW_FORWARD ?
			   fftw_twiddle_generic : fftwi_twiddle_generic);
	  int size = factor(n);
	  fftw_plan r = planner(table, n / size, dir, flags);

	  node = make_node_generic(n, size, codelet, r->root, flags);
	  newplan = make_plan(n, dir, node, flags, FFTW_GENERIC, 0);
	  use_plan(newplan);
	  destroy_plan(r);
	  compute_cost(newplan);
	  best = pick_better(newplan, best);
     }

     return best;
}

fftw_plan planner(fftw_plan *table, int n, fftw_direction dir,
			   int flags)
{
     fftw_plan best = (fftw_plan) 0;

     /* see if plan has already been computed */
     best = lookup(table, n, flags);
     if (best) {
	  use_plan(best);
	  return best;
     }

     /* try a wise plan */
     best = planner_wisdom(table, n, dir, flags);

     if (!best) {
	  /* No wisdom.  Plan normally. */
	  best = planner_normal(table, n, dir, flags);
     }

     if (best) {
	  insert(table, best, n);

	  /* remember the wisdom */
	  fftw_wisdom_add(n, flags, dir, best->wisdom_type,
			  best->wisdom_signature);
     }

     return best;
}

fftw_plan fftw_create_plan(int n, fftw_direction dir, int flags)
{
     fftw_plan table;
     fftw_plan p1;

     /* validate parameters */
     if (n <= 0)
	  return (fftw_plan) 0;

     if ((dir != FFTW_FORWARD) && (dir != FFTW_BACKWARD))
	  return (fftw_plan) 0;

     make_empty_table(&table);
     p1 = planner(&table, n, dir, flags);
     destroy_table(&table);

     complete_twiddle(p1->root, n);
     return p1;
}

void fftw_destroy_plan(fftw_plan plan)
{
     destroy_plan(plan);
}

static void print_node(FILE * f, fftw_plan_node *p, int indent)
{
     if (p) {
	  switch (p->type) {
	      case FFTW_NOTW:
		   fprintf(f, "%*sFFTW_NOTW %d\n", indent, "",
			   p->nodeu.notw.size);
		   break;
	      case FFTW_TWIDDLE:
		   fprintf(f, "%*sFFTW_TWIDDLE %d\n", indent, "",
			   p->nodeu.twiddle.size);
		   print_node(f, p->nodeu.twiddle.recurse, indent);
		   break;
	      case FFTW_GENERIC:
		   fprintf(f, "%*sFFTW_GENERIC %d\n", indent, "",
			   p->nodeu.generic.size);
		   print_node(f, p->nodeu.generic.recurse, indent);
		   break;
	  }
     }
}

void fftw_fprint_plan(FILE * f, fftw_plan p)
{
     fprintf(f, "plan: (cost = %e)\n", p->cost);
     print_node(f, p->root, 0);
}

void fftw_print_plan(fftw_plan p)
{
     fftw_fprint_plan(stdout, p);
}
