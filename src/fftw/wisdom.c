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
 * wisdom.c -- manage the wisdom
 */

#include <fftw.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

struct wisdom {
     int n;
     int flags;
     fftw_direction dir;
     enum fftw_node_type type;	/* this is the wisdom */
     int signature;		/* this is the wisdom */
     struct wisdom *next;
};

/* list of wisdom */
static struct wisdom *wisdom_list = (struct wisdom *) 0;

int fftw_wisdom_lookup(int n, int flags, fftw_direction dir,
		     enum fftw_node_type *type,
		     int *signature, int replacep)
{
     struct wisdom *p;

     if (!(flags & FFTW_USE_WISDOM))
	  return 0;		/* simply ignore if wisdom is disabled */

     flags |= FFTW_MEASURE; /* always use (only) wisdom from measurements */

     for (p = wisdom_list; p; p = p->next) {
	  if (p->n == n && p->flags == flags && p->dir == dir) {
	       /* found wisdom */
	       if (replacep) {
		    /* replace old wisdom with new */
		    p->type = *type;
		    p->signature = *signature;
	       } else {
		    *type = p->type;
		    *signature = p->signature;
	       }
	       return 1;
	  }
     }

     return 0;
}

void fftw_wisdom_add(int n, int flags, fftw_direction dir,
		   enum fftw_node_type type,
		   int signature)
{
     struct wisdom *p;

     if (!(flags & FFTW_USE_WISDOM))
	  return;		/* simply ignore if wisdom is disabled */

     if (!(flags & FFTW_MEASURE))
	  return;  /* only measurements produce wisdom */

     if (fftw_wisdom_lookup(n, flags, dir, &type, &signature, 1))
	  return;		/* wisdom overwrote old wisdom */

     p = (struct wisdom *) fftw_malloc(sizeof(struct wisdom));

     p->n = n;
     p->flags = flags;
     p->dir = dir;
     p->type = type;
     p->signature = signature;

     /* remember this wisdom */
     p->next = wisdom_list;
     wisdom_list = p;
}

void fftw_forget_wisdom(void)
{
     while (wisdom_list) {
	  struct wisdom *p;

	  p = wisdom_list;
	  wisdom_list = wisdom_list->next;
	  fftw_free(p);
     }
}

/*
 * user-visible routines, to convert wisdom into strings etc.
 */
#define WISDOM_FORMAT_VERSION "FFTW-1.2"

static void (*emit)(char c, void *data);

static void emit_string(char *s, void *data)
{
     while (*s) 
	  emit(*s++, data);
}

static void emit_int(int n, void *data)
{
     char buf[128];

     sprintf(buf, "%d", n);
     emit_string(buf, data);
}

/* dump wisdom in lisp-like format */
void fftw_export_wisdom(void (*emitter)(char c, void *), void *data)
{
     struct wisdom *p;

     /* install the output handler */
     emit = emitter;

     emit('(',data);
     emit_string(WISDOM_FORMAT_VERSION,data);

     for (p = wisdom_list; p; p = p->next) {
	  emit(' ',data);	/* separator to make the output nicer */
	  emit('(',data);
	  emit_int((int) p->n, data);
	  emit(' ',data);
	  emit_int((int) p->flags, data);
	  emit(' ',data);
	  emit_int((int) p->dir, data);
	  emit(' ',data);
	  emit_int((int) p->type, data);
	  emit(' ',data);
	  emit_int((int) p->signature, data);
	  emit(')',data);
     }
     emit(')',data);
}

/* input part */
static int next_char;
static int (*get_input)(void *data);
static fftw_status input_error;

static void read_char(void *data)
{
     next_char = get_input(data);
     if (next_char == 0 ||
	 next_char == EOF)
	  input_error = FFTW_FAILURE;
}

/* skip blanks, newlines, tabs, etc */
static void eat_blanks(void *data)
{
     while (isspace(next_char))
	  read_char(data);
}

static int read_int(void *data)
{
     int sign = 1;
     int n = 0;

     eat_blanks(data);
     if (next_char == '-') {
	  sign = -1;
	  read_char(data);
	  eat_blanks(data);
     }

     if (!isdigit(next_char)) {
	  /* error, no digit */
	  input_error = FFTW_FAILURE;
	  return 0;
     }

     while (isdigit(next_char)) {
	  n = n * 10 + (next_char - '0');
	  read_char(data);
     }

     return sign * n;
}

#define EXPECT(c)                     \
{				      \
     eat_blanks(data);		      \
     if (input_error == FFTW_FAILURE || \
         next_char != c)	      \
	  return FFTW_FAILURE;	      \
     read_char(data);		      \
}
				      
#define EXPECT_INT(n)                                 \
{				                      \
     n = read_int(data);	                      \
     if (input_error == FFTW_FAILURE)                 \
	  return FFTW_FAILURE;		              \
}				      
				      
#define EXPECT_STRING(s)             \
{                                    \
     char *s1 = s;		     \
     while (*s1) {		     \
	  EXPECT(*s1);		     \
	  ++s1;			     \
     }				     \
}       			      
                                      
fftw_status fftw_import_wisdom(int (*g)(void *), void *data)
{
     int n;
     int flags;
     fftw_direction dir;
     enum fftw_node_type type;
     int signature;

     get_input = g;
     input_error = FFTW_SUCCESS;

     read_char(data);

     eat_blanks(data);
     EXPECT('(');
     eat_blanks(data);
     EXPECT_STRING(WISDOM_FORMAT_VERSION);
     eat_blanks(data);

     while (next_char != ')') {
	  EXPECT('(');
	  EXPECT_INT(n);
	  EXPECT_INT(flags);
	  EXPECT_INT(dir);
	  EXPECT_INT(type);
	  EXPECT_INT(signature);
	  eat_blanks(data);
	  EXPECT(')');

	  /* the wisdom has been read properly. Add it */
	  fftw_wisdom_add(n, flags, dir, type, signature);

	  /* prepare for next morsel of wisdom */
	  eat_blanks(data);
     }

     return FFTW_SUCCESS;
}
