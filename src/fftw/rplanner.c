/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
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
 * planner.c -- find the optimal plan
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <stdlib.h>
#include <stdio.h>

#include <fftw-int.h>
#include <rfftw.h>

extern fftw_codelet_desc *rfftw_config[];	/* global from rconfig.c */
extern fftw_rgeneric_codelet fftw_hc2hc_forward_generic;
extern fftw_rgeneric_codelet fftw_hc2hc_backward_generic;

/* undocumented debugging hook */
void (*rfftw_plan_hook)(fftw_plan plan) = (void (*)(fftw_plan)) 0;

/* timing rfftw plans: */
static double rfftw_measure_runtime(fftw_plan plan,
				    fftw_real *in, int istride,
				    fftw_real *out, int ostride)
{
     fftw_time begin, end, start;
     double t, tmin;
     int i, iter;
     int n;
     int repeat;

     n = plan->n;

     iter = 1;

     for (;;) {
	  tmin = 1.0E10;
	  for (i = 0; i < n; ++i)
	       in[istride * i] = 0.0;

	  start = fftw_get_time();
	  /* repeat the measurement FFTW_TIME_REPEAT times */
	  for (repeat = 0; repeat < FFTW_TIME_REPEAT; ++repeat) {
	       begin = fftw_get_time();
	       for (i = 0; i < iter; ++i)
		    rfftw(plan, 1, in, istride, 0, out, ostride, 0);
	       end = fftw_get_time();

	       t = fftw_time_to_sec(fftw_time_diff(end, begin));
	       if (t < tmin)
		    tmin = t;

	       /* do not run for too long */
	       t = fftw_time_to_sec(fftw_time_diff(end, start));
	       if (t > FFTW_TIME_LIMIT)
		    break;
	  }

	  if (tmin >= FFTW_TIME_MIN)
	       break;

	  iter *= 2;
     }

     tmin /= (double) iter;
     return tmin;
}

/* auxiliary functions */
static void rcompute_cost(fftw_plan plan,
			  fftw_real *in, int istride,
			  fftw_real *out, int ostride)
{
     if (plan->flags & FFTW_MEASURE)
	  plan->cost = rfftw_measure_runtime(plan, in, istride, out, ostride);
     else {
	  double c;
	  c = plan->n * fftw_estimate_node(plan->root);
	  plan->cost = c;
     }
}

static void run_plan_hooks(fftw_plan p)
{
     if (rfftw_plan_hook && p) {
	  fftw_complete_twiddle(p->root, p->n);
	  rfftw_plan_hook(p);
     }
}

/* macrology */
#define FOR_ALL_RCODELETS(p) \
   fftw_codelet_desc **__q, *p;                         \
   for (__q = &rfftw_config[0]; (p = (*__q)); ++__q)

/******************************************
 *      Recursive planner                 *
 ******************************************/
static fftw_plan rplanner(fftw_plan *table, int n,
			  fftw_direction dir, int flags,
			  fftw_real *, int, fftw_real *, int);

/*
 * the planner consists of two parts: one that tries to
 * use accumulated wisdom, and one that does not.
 * A small driver invokes both parts in sequence
 */

/* planner with wisdom: look up the codelet suggested by the wisdom */
static fftw_plan rplanner_wisdom(fftw_plan *table, int n,
				 fftw_direction dir, int flags,
				 fftw_real *in, int istride,
				 fftw_real *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan_node *node;
     int have_wisdom;
     enum fftw_node_type wisdom_type;
     int wisdom_signature;

     /* see if we remember any wisdom for this case */
     have_wisdom = fftw_wisdom_lookup(n, flags, dir, RFFTW_WISDOM,
				      istride, ostride,
				      &wisdom_type, &wisdom_signature, 0);

     if (!have_wisdom)
	  return best;

     if (wisdom_type == FFTW_REAL2HC || wisdom_type == FFTW_HC2REAL) {
	  FOR_ALL_RCODELETS(p) {
	       if (p->dir == dir && p->type == wisdom_type) {
		    /* see if wisdom applies */
		    if (wisdom_signature == p->signature &&
			p->size == n) {
			 if (wisdom_type == FFTW_REAL2HC)
			      node = fftw_make_node_real2hc(n, p);
			 else
			      node = fftw_make_node_hc2real(n, p);
			 best = fftw_make_plan(n, dir, node, flags,
					       p->type, p->signature);
			 fftw_use_plan(best);
			 run_plan_hooks(best);
			 return best;
		    }
	       }
	  }
     }
     if (wisdom_type == FFTW_HC2HC) {
	  FOR_ALL_RCODELETS(p) {
	       if (p->dir == dir && p->type == wisdom_type) {

		    /* see if wisdom applies */
		    if (wisdom_signature == p->signature &&
			p->size > 1 &&
			(n % p->size) == 0) {
			 fftw_plan r = rplanner(table, n / p->size, dir,
						flags, in, istride, out,
						ostride);
			 if (!r)
			      continue;
			 node = fftw_make_node_hc2hc(n, dir, p,
						     r->root, flags);
			 best = fftw_make_plan(n, dir, node, flags,
					       p->type, p->signature);
			 fftw_use_plan(best);
			 run_plan_hooks(best);
			 fftw_destroy_plan_internal(r);
			 return best;
		    }
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

static fftw_plan rplanner_normal(fftw_plan *table, int n, fftw_direction dir,
				 int flags, fftw_real *in, int istride,
				 fftw_real *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan newplan;
     fftw_plan_node *node;

     /* see if we have any codelet that solves the problem */
     {
	  FOR_ALL_RCODELETS(p) {
	       if (p->dir == dir &&
		   (p->type == FFTW_REAL2HC || p->type == FFTW_HC2REAL)) {
		    if (p->size == n) {
			 if (p->type == FFTW_REAL2HC)
			      node = fftw_make_node_real2hc(n, p);
			 else
			      node = fftw_make_node_hc2real(n, p);
			 newplan = fftw_make_plan(n, dir, node, flags,
						  p->type, p->signature);
			 fftw_use_plan(newplan);
			 run_plan_hooks(newplan);
			 rcompute_cost(newplan, in, istride, out, ostride);
			 best = fftw_pick_better(newplan, best);
		    }
	       }
	  }
     }

     /* Then, try all available twiddle codelets */
     {
	  FOR_ALL_RCODELETS(p) {
	       if (p->dir == dir && p->type == FFTW_HC2HC) {
		    if ((n % p->size) == 0 &&
			p->size > 1 &&
			(!best || n != p->size)) {
			 fftw_plan r = rplanner(table, n / p->size, dir, flags,
					      in, istride, out, ostride);
			 if (!r)
			      continue;
			 node = fftw_make_node_hc2hc(n, dir, p,
						     r->root, flags);
			 newplan = fftw_make_plan(n, dir, node, flags,
						  p->type, p->signature);
			 fftw_use_plan(newplan);
			 run_plan_hooks(newplan);
			 fftw_destroy_plan_internal(r);
			 rcompute_cost(newplan, in, istride, out, ostride);
			 best = fftw_pick_better(newplan, best);
		    }
	       }
	  }
     }

     /* 
      * Resort to generic codelets for unknown factors, but only if
      * n is odd--the rgeneric codelets can't handle even n's.
      */
     if (n % 2 != 0) {
	  fftw_rgeneric_codelet *codelet = (dir == FFTW_FORWARD ?
					    fftw_hc2hc_forward_generic :
					    fftw_hc2hc_backward_generic);
	  int size, prev_size = 0, remaining_factors = n;
	  fftw_plan r;

	  while (remaining_factors > 1) {
	       size = fftw_factor(remaining_factors);
	       remaining_factors /= size;

	       /* don't try the same factor more than once */
	       if (size == prev_size)
		    continue;
	       prev_size = size;

	       /* Look for codelets corresponding to this factor. */
	       {
		    FOR_ALL_RCODELETS(p) {
			 if (p->dir == dir && p->type == FFTW_HC2HC
			     && p->size == size) {
			      size = 0;
			      break;
			 }
		    }
	       }

	       /* 
	        * only try a generic/rader codelet if there were no 
	        * twiddle codelets for this factor 
	        */
	       if (!size)
		    continue;

	       r = rplanner(table, n / size, dir, flags,
			    in, istride, out, ostride);

	       node = fftw_make_node_rgeneric(n, size, dir, codelet,
					      r->root, flags);
	       newplan = fftw_make_plan(n, dir, node, flags, FFTW_RGENERIC, 0);
	       fftw_use_plan(newplan);
	       run_plan_hooks(newplan);
	       fftw_destroy_plan_internal(r);
	       rcompute_cost(newplan, in, istride, out, ostride);
	       best = fftw_pick_better(newplan, best);
	  }
     }
     return best;
}

static fftw_plan rplanner(fftw_plan *table, int n, fftw_direction dir,
			  int flags, fftw_real *in, int istride,
			  fftw_real *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;

     /* see if plan has already been computed */
     best = fftw_lookup(table, n, flags);
     if (best) {
	  fftw_use_plan(best);
	  return best;
     }
     /* try a wise plan */
     best = rplanner_wisdom(table, n, dir, flags,
			    in, istride, out, ostride);

     if (!best) {
	  /* No wisdom.  Plan normally. */
	  best = rplanner_normal(table, n, dir, flags,
				 in, istride, out, ostride);
     }
     if (best) {
	  fftw_insert(table, best, n);

	  /* remember the wisdom */
	  fftw_wisdom_add(n, flags, dir, RFFTW_WISDOM,
			  istride, ostride,
			  best->wisdom_type,
			  best->wisdom_signature);
     }
     return best;
}

fftw_plan rfftw_create_plan_specific(int n, fftw_direction dir, int flags,
				     fftw_real *in, int istride,
				     fftw_real *out, int ostride)
{
     fftw_plan table;
     fftw_plan p1;

     /* validate parameters */
     if (n <= 0)
	  return (fftw_plan) 0;

     if ((dir != FFTW_FORWARD) && (dir != FFTW_BACKWARD))
	  return (fftw_plan) 0;

     fftw_make_empty_table(&table);
     p1 = rplanner(&table, n, dir, flags,
		   in, istride, out, ostride);
     fftw_destroy_table(&table);

     if (p1)
	  fftw_complete_twiddle(p1->root, n);
     return p1;
}

fftw_plan rfftw_create_plan(int n, fftw_direction dir, int flags)
{
     fftw_real *tmp_in;
     fftw_real *tmp_out;
     fftw_plan p;

     if (flags & FFTW_MEASURE) {
	  tmp_in = (fftw_real *) fftw_malloc(n * sizeof(fftw_real));
	  tmp_out = (fftw_real *) fftw_malloc(n * sizeof(fftw_real));
	  if (!tmp_in || !tmp_out)
	       return 0;

	  p = rfftw_create_plan_specific(n, dir, flags,
					 tmp_in, 1, tmp_out, 1);

	  fftw_free(tmp_in);
	  fftw_free(tmp_out);
     } else
	  p = rfftw_create_plan_specific(n, dir, flags,
				 (fftw_real *) 0, 1, (fftw_real *) 0, 1);

     return p;
}

void rfftw_destroy_plan(fftw_plan plan)
{
     fftw_destroy_plan_internal(plan);
}

void rfftw_fprint_plan(FILE *f, fftw_plan p)
{
     fftw_fprint_plan(f, p);
}

void rfftw_print_plan(fftw_plan p)
{
     rfftw_fprint_plan(stdout, p);
}
