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
 * planner.c -- find the optimal plan
 */

/* $Id$ */
#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw-int.h>
#include <stdlib.h>
#include <stdio.h>

extern fftw_generic_codelet fftw_twiddle_generic;
extern fftw_generic_codelet fftwi_twiddle_generic;
extern fftw_codelet_desc *fftw_config[];

static void init_test_array(fftw_complex *arr, int stride, int n)
{
     int j;

     for (j = 0; j < n; ++j) {
	  c_re(arr[stride * j]) = 0.0;
	  c_im(arr[stride * j]) = 0.0;
     }
}

/*
 * The timer keeps doubling the number of iterations
 * until the program runs for more than FFTW_TIME_MIN
 */
static double fftw_measure_runtime(fftw_plan plan,
				   fftw_complex *in, int istride,
				   fftw_complex *out, int ostride)
{
     fftw_time begin, end, start;
     double t, tmax, tmin;
     int i, iter;
     int n;
     int repeat;

     n = plan->n;

     iter = 1;

     for (;;) {
	  tmin = 1.0E10;
	  tmax = -1.0E10;
	  init_test_array(in, istride, n);

	  start = fftw_get_time();
	  /* repeat the measurement FFTW_TIME_REPEAT times */
	  for (repeat = 0; repeat < FFTW_TIME_REPEAT; ++repeat) {
	       begin = fftw_get_time();
	       for (i = 0; i < iter; ++i) {
		    fftw(plan, 1, in, istride, 0, out, ostride, 0);
	       }
	       end = fftw_get_time();

	       t = fftw_time_to_sec(fftw_time_diff(end, begin));
	       if (t < tmin)
		    tmin = t;
	       if (t > tmax)
		    tmax = t;

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
     tmax /= (double) iter;

     return tmin;
}

/* auxiliary functions */
static void compute_cost(fftw_plan plan,
			 fftw_complex *in, int istride,
			 fftw_complex *out, int ostride)
{
     if (plan->flags & FFTW_MEASURE)
	  plan->cost = fftw_measure_runtime(plan, in, istride, out, ostride);
     else {
	  double c;
	  c = plan->n * fftw_estimate_node(plan->root);
	  plan->cost = c;
     }
}

/* macrology */
#define FOR_ALL_CODELETS(p) \
   fftw_codelet_desc **__q, *p;                         \
   for (__q = &fftw_config[0]; (p = (*__q)); ++__q)

/******************************************
 *      Recursive planner                 *
 ******************************************/
static fftw_plan planner(fftw_plan *table, int n, fftw_direction dir, int flags,
			 fftw_complex *, int, fftw_complex *, int);

/*
 * the planner consists of two parts: one that tries to
 * use accumulated wisdom, and one that does not.
 * A small driver invokes both parts in sequence
 */

/* planner with wisdom: look up the codelet suggested by the wisdom */
static fftw_plan planner_wisdom(fftw_plan *table, int n,
				fftw_direction dir, int flags,
				fftw_complex *in, int istride,
				fftw_complex *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan_node *node;
     int have_wisdom;
     enum fftw_node_type wisdom_type;
     int wisdom_signature;

     /* see if we remember any wisdom for this case */
     have_wisdom = fftw_wisdom_lookup(n, flags, dir, FFTW_WISDOM,
				      istride, ostride,
				      &wisdom_type, &wisdom_signature, 0);

     if (!have_wisdom)
	  return best;

     if (wisdom_type == FFTW_NOTW) {
	  FOR_ALL_CODELETS(p) {
	       if (p->dir == dir && p->type == wisdom_type) {
		    /* see if wisdom applies */
		    if (wisdom_signature == p->signature &&
			p->size == n) {
			 node = fftw_make_node_notw(n, p);
			 best = fftw_make_plan(n, dir, node, flags,
					       p->type, p->signature);
			 fftw_use_plan(best);
			 return best;
		    }
	       }
	  }
     }
     if (wisdom_type == FFTW_TWIDDLE) {
	  FOR_ALL_CODELETS(p) {
	       if (p->dir == dir && p->type == wisdom_type) {

		    /* see if wisdom applies */
		    if (wisdom_signature == p->signature &&
			p->size > 1 &&
			(n % p->size) == 0) {
			 fftw_plan r = planner(table, n / p->size, dir, flags,
					       in, istride, out, ostride);
			 node = fftw_make_node_twiddle(n, p,
						       r->root, flags);
			 best = fftw_make_plan(n, dir, node, flags,
					       p->type, p->signature);
			 fftw_use_plan(best);
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
static fftw_plan planner_normal(fftw_plan *table, int n, fftw_direction dir,
				int flags,
				fftw_complex *in, int istride,
				fftw_complex *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;
     fftw_plan newplan;
     fftw_plan_node *node;

     /* see if we have any codelet that solves the problem */
     {
	  FOR_ALL_CODELETS(p) {
	       if (p->dir == dir && p->type == FFTW_NOTW) {
		    if (p->size == n) {
			 node = fftw_make_node_notw(n, p);
			 newplan = fftw_make_plan(n, dir, node, flags,
						  p->type, p->signature);
			 fftw_use_plan(newplan);
			 compute_cost(newplan, in, istride, out, ostride);
			 best = fftw_pick_better(newplan, best);
		    }
	       }
	  }
     }

     /* Then, try all available twiddle codelets */
     {
	  FOR_ALL_CODELETS(p) {
	       if (p->dir == dir && p->type == FFTW_TWIDDLE) {
		    if ((n % p->size) == 0 &&
			p->size > 1 &&
			(!best || n != p->size)) {
			 fftw_plan r = planner(table, n / p->size, dir, flags,
					       in, istride, out, ostride);
			 node = fftw_make_node_twiddle(n, p,
						       r->root, flags);
			 newplan = fftw_make_plan(n, dir, node, flags,
						  p->type, p->signature);
			 fftw_use_plan(newplan);
			 fftw_destroy_plan_internal(r);
			 compute_cost(newplan, in, istride, out, ostride);
			 best = fftw_pick_better(newplan, best);
		    }
	       }
	  }
     }

     /* 
      * resort to generic or rader codelets for unknown factors
      */
     {
	  fftw_generic_codelet *codelet = (dir == FFTW_FORWARD ?
					   fftw_twiddle_generic :
					   fftwi_twiddle_generic);
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
		    FOR_ALL_CODELETS(p) {
			 if (p->dir == dir && p->type == FFTW_TWIDDLE
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

	       r = planner(table, n / size, dir, flags,
			   in, istride, out, ostride);

	       /* Try Rader codelet: */
	       node = fftw_make_node_rader(n, size, dir, r->root, flags);
	       newplan = fftw_make_plan(n, dir, node, flags, FFTW_RADER, 0);
	       fftw_use_plan(newplan);
	       compute_cost(newplan, in, istride, out, ostride);
	       best = fftw_pick_better(newplan, best);

	       if (size < 100) {	/*
					 * only try generic for small 
					 * sizes 
					 */
		    /* Try generic codelet: */
		    node = fftw_make_node_generic(n, size, codelet,
						  r->root, flags);
		    newplan = fftw_make_plan(n, dir, node, flags,
					     FFTW_GENERIC, 0);
		    fftw_use_plan(newplan);
		    compute_cost(newplan, in, istride, out, ostride);
		    best = fftw_pick_better(newplan, best);
	       }
	       fftw_destroy_plan_internal(r);
	  }
     }

     if (!best)
	  fftw_die("bug in planner");

     return best;
}

static fftw_plan planner(fftw_plan *table, int n, fftw_direction dir,
			 int flags,
			 fftw_complex *in, int istride,
			 fftw_complex *out, int ostride)
{
     fftw_plan best = (fftw_plan) 0;

     /* see if plan has already been computed */
     best = fftw_lookup(table, n, flags);
     if (best) {
	  fftw_use_plan(best);
	  return best;
     }
     /* try a wise plan */
     best = planner_wisdom(table, n, dir, flags,
			   in, istride, out, ostride);

     if (!best) {
	  /* No wisdom.  Plan normally. */
	  best = planner_normal(table, n, dir, flags,
				in, istride, out, ostride);
     }
     if (best) {
	  fftw_insert(table, best, n);

	  /* remember the wisdom */
	  fftw_wisdom_add(n, flags, dir, FFTW_WISDOM, istride, ostride,
			  best->wisdom_type,
			  best->wisdom_signature);
     }
     return best;
}

fftw_plan fftw_create_plan_specific(int n, fftw_direction dir, int flags,
				    fftw_complex *in, int istride,
				    fftw_complex *out, int ostride)
{
     fftw_plan table;
     fftw_plan p1;

     /* validate parameters */
     if (n <= 0)
	  return (fftw_plan) 0;

     if ((dir != FFTW_FORWARD) && (dir != FFTW_BACKWARD))
	  return (fftw_plan) 0;

     fftw_make_empty_table(&table);
     p1 = planner(&table, n, dir, flags,
		  in, istride, out, ostride);
     fftw_destroy_table(&table);

     fftw_complete_twiddle(p1->root, n);
     return p1;
}

fftw_plan fftw_create_plan(int n, fftw_direction dir, int flags)
{
     fftw_complex *tmp_in, *tmp_out;
     fftw_plan p;

     if (flags & FFTW_MEASURE) {
	  tmp_in = (fftw_complex *) fftw_malloc(2 * n * sizeof(fftw_complex));
	  if (!tmp_in)
	       return 0;
	  tmp_out = tmp_in + n;

	  p = fftw_create_plan_specific(n, dir, flags,
					tmp_in, 1, tmp_out, 1);

	  fftw_free(tmp_in);
     } else
	  p = fftw_create_plan_specific(n, dir, flags,
			   (fftw_complex *) 0, 1, (fftw_complex *) 0, 1);

     return p;
}

void fftw_destroy_plan(fftw_plan plan)
{
     fftw_destroy_plan_internal(plan);
}
