#ifndef _metacode_h
#define _metacode_h
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <types/simple.h>

static char *SRCID_metacode_h = "";

#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

#define FCON "     &" /* continuation line in f77 */
#define max(a,b) (((a) > (b)) ? (a) : (b))


/* Array referencing shortcut */
#define ARRAY(a,idx)          _array(#a,#idx)

typedef struct {
  char typename[150];
  char name[150];
  bool breferenced;
  bool bvector;
  bool bconst;
  char constval[50];
} decl_t;         /* Argument and variable buffer element */

extern int prec;  /* precision (4=single, 8=double) */
extern int IND;   /* current indentation */
extern char *codebuffer; /* buffer to which code is written */
extern char header[10000]; /* buffer for info and loop name */
extern FILE *output;  /* output file */
extern decl_t *decl_list; /* list with args and vars */
extern int ndecl; /* length of above list */
extern int nargs; /* first nargs are arguments, rest is vars */

extern bool bC;   

/* Concatenate a string to a buffer with plus sign between terms. */
void add_to_buffer(char *buffer,char *term);
  
/* initiate output buffers */
void init_metacode(void);

/* write a function to file and empty buffers */
void flush_buffers(void);

/* Return the correct indentation as a string */
char *indent(void);

/* Print a line of code to the output file */
void code(char *fmt, ...);

void newline(void);

/* Add a comment */
void comment(char *s);

/* Define a new floating-point variable */
void declare_real(char *name);
void declare_real_vector(char *name);

void declare_const_real(char *name,double val);

void declare_int(char *name);
void declare_int_vector(char *name);
void declare_int4(char *name);
void declare_int8(char *name);
void declare_intreal(char *name);
void declare_other(char *typename,char *name);

/* Cray vector pragma */
void vector_pragma(void);


char *_array(char *a,char *idx, ...);

void _p_state(char *left,char *right,char *symb);

void file_error(char *fn);

void assign(char *left, char *right, ...);

void increment(char *left,char *right, ...);

void decrement(char *left,char *right, ...);

void add(char *left,char *r1,char *r2, ...);

void subtract(char *left,char *r1,char *r2, ...);

void multiply(char *left,char *r1,char *r2, ...);

void closeit(void);

void usage(int argc,char *argv[]);

int count_lines(char *fn);

void edit_warning(char *fn);

void start_loop(char *lvar,char *from,char *to);

void start_stride_loop(char *lvar,char *from,char *to, char *stride);

void end_loop(void);

void start_if(char *cond);
void do_else(void);
void end_if(void);

void close_func(void);

#endif

