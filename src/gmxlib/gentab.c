/*
 *	@(#) gentab.c 1.9 11/5/92
 *
 *
 *	GROMACS - Groningen Machine for Chemical Simulation
 *	Copyright (c) 1990, 1991, 1992, Groningen University
 *
 */

#include <stdio.h>
#include <math.h>

#define tversion  "\tInvsqrt lookup tables @(#) gentab.c 1.9 11/5/92"
#define copyright "\tCopyright (c) 1992, A. Sijbers, Groningen University"

/*
   Fast sqrt(x) routine (actually 1.0/sqrt(x)) By A. Sijbers (c) 1992
   
   Based on Newton Rhapson iteration, initial value from lookup table

   Floating point representation:

         msb                                   lsb
         3322 2222 2222 1111 1111 1000 0000 0000
   bitnr 1098 7654 3210 9876 5432 1098 7654 3210
         |||| |||| |+-- ---- ---- ---- ---- ---+ fraction
         |+-- ---- + exponent
         +sign

   S = sign
   E = exponent
   F = fraction
                      S  E-127
   IEEE : value = (-1) (2     ) (1.F)     <------ must be float representation

                      S  E-128
   DEC  : value = (-1) (2     ) (0.1F)
*/

#define	EXP_LSB		0x00800000
#define	EXP_SEED_SIZE	256
#define	EXP_MASK	0x7f800000
#define	EXP_SHIFT	23
#define	MAX_FRACT	0x007fffff
#define	FRACT_MASK	0x007fffff
#define	FRACT_SIZE	11              /* significant part of fraction */
#define	FRACT_SHIFT	(EXP_SHIFT-FRACT_SIZE)
#define	FRACT_SEED_SIZE	(1<<(FRACT_SIZE+1))   /* one bit for even/odd */
#define	FRACT_FIRST	(0x3f800000>>FRACT_SHIFT)
#define	NOT_INITIALISED	~0
#define	EXP_ADDR(val)	(((val)&EXP_MASK)>>EXP_SHIFT)
#define	FRACT_ADDR(val)	(((val)&(FRACT_MASK|EXP_LSB))>>FRACT_SHIFT)

typedef unsigned int word;

typedef union 
{
  word bval;
  float fval;
} t_convert;

typedef struct
{
  word exp_seed[EXP_SEED_SIZE];
  word fract_seed[FRACT_SEED_SIZE];
} t_lutab;

t_lutab lookup_table;

float invsqrt(float x)
{
  t_convert result,bit_pattern;
  word  exp,fract;
  float y,lu;
  
  bit_pattern.fval = x;
  exp              = EXP_ADDR(bit_pattern.bval);
  fract            = FRACT_ADDR(bit_pattern.bval);
  result.bval      = (lookup_table.exp_seed[exp] | 
		      lookup_table.fract_seed[fract]);
  lu               = result.fval;
  
  y = (0.5*lu*(3.0-((x*lu)*lu)));  
  return y;
}

void init_table(t_lutab *lookup_table)
{
  word index,addr;
  int new_exp,i,exp;
  t_convert bit_in,bit_out;
  
  for (i=0; i<EXP_SEED_SIZE; i++)
    lookup_table->exp_seed[i]=NOT_INITIALISED;
  for (i=0; i<FRACT_SEED_SIZE; i++) 
    lookup_table->fract_seed[i]=NOT_INITIALISED;
  for (i=0; i<EXP_SEED_SIZE; i++)
    {
      if ((exp=(i-127))<0) new_exp=127-((exp+1)/2); else new_exp=127-(exp/2)-1;
      lookup_table->exp_seed[i]=((new_exp)<<EXP_SHIFT)&EXP_MASK;
    }
  index=FRACT_FIRST; 
  for (i=0; i<FRACT_SEED_SIZE; i++)
    {
      bit_in.bval=(index<<FRACT_SHIFT);
      addr=FRACT_ADDR(bit_in.bval);
      bit_out.fval=(1.0/sqrt(bit_in.fval));
      lookup_table->fract_seed[addr]=(bit_out.bval&FRACT_MASK);
      if (lookup_table->fract_seed[addr]==0) 
        lookup_table->fract_seed[addr]=(MAX_FRACT&FRACT_MASK);
      index++;
    }
}

void dump_asm_list(FILE *fp,word list[],int n)
{
  int i;
  
  for (i=0; i<n; i++)
    {
      if ((i%4)!=0)
        (void) fprintf(fp,",");
      else
        {
          if (i!=0) (void) fprintf(fp,"\n");
          (void) fprintf(fp,"\t.long\t");
        }
      (void) fprintf(fp,"0x%.8x",list[i]);
    }
  (void) fprintf(fp,"\n");
}

static int ilog2(int x)
{
  int i;
  
  for (i=-1; x!=0; i++) x>>=1;
  return i;
}

int write_asm_table(FILE *fp,char *fn,char *table,t_lutab *lookup_table)
{
  int size2log2;
  
  (void) fprintf(fp,"//\n");
  (void) fprintf(fp,"//%s\n",tversion);
  (void) fprintf(fp,"//\n");
  (void) fprintf(fp,"//%s\n",copyright);
  (void) fprintf(fp,"//\n");
  size2log2=ilog2(sizeof(lookup_table->exp_seed[0]));
  (void) fprintf(fp,"exptabs\t=\t\t%d\t// exponent table size\n",
                 sizeof(lookup_table->exp_seed));
  (void) fprintf(fp,"exp_s\t=\t\t%d\t// exponent shift\n",
                 EXP_SHIFT-size2log2);
  (void) fprintf(fp,"exp_m\t=\t\t0x%.4x\t// exponent mask\n",
                 ((unsigned)(EXP_MASK)>>EXP_SHIFT)<<size2log2);
  size2log2=ilog2(sizeof(lookup_table->exp_seed[0]));
  (void) fprintf(fp,"frac_s\t=\t\t%d\t// fraction shift\n",
                 FRACT_SHIFT-size2log2);
  (void) fprintf(fp,"frac_m\t=\t\t0x%.4x\t// fraction mask\n",
                 (FRACT_MASK>>FRACT_SIZE)<<size2log2);
  (void) fprintf(fp,"//\n");
  (void) fprintf(fp,"\t.file\t\t\"%s\"\n",fn);
  (void) fprintf(fp,"//\n");
  (void) fprintf(fp,"\t.data\n");
  (void) fprintf(fp,"\t.align\t\t8\n");
  (void) fprintf(fp,"//\n// start of exponent table\n");
  (void) fprintf(fp,"//\n%s:\n",table);
  dump_asm_list(fp,lookup_table->exp_seed,EXP_SEED_SIZE);
  (void) fprintf(fp,"//\n// start of fraction table\n//\n");
  dump_asm_list(fp,lookup_table->fract_seed,FRACT_SEED_SIZE);
  (void) fprintf(fp,"//\n// end of tables\n//\n");
  return 0;
}

void dump_c_list(FILE *fp,word list[],int n)
{
  int i;
  
  (void) fprintf(fp,"  {");
  for (i=0; i<n; i++)
    {
      if ((i%4)==0) (void) fprintf(fp,"\n    ");
      (void) fprintf(fp,"0x%.8x",list[i]);
      if (i<(n-1)) (void) fprintf(fp,",");
    }
  (void) fprintf(fp,"\n  }");
}

int write_c_table(FILE *fp,char *fn,char *table,t_lutab *lookup_table)
{
  (void) fprintf(fp,"/*\n");
  (void) fprintf(fp,"%s%s\n"," *",tversion);
  (void) fprintf(fp," *\n");
  (void) fprintf(fp," *%s\n",copyright);
  (void) fprintf(fp," */\n");
  (void) fprintf(fp,"t_lutab %s={\n",table);
  dump_c_list(fp,lookup_table->exp_seed,EXP_SEED_SIZE);
  (void) fprintf(fp,",\n");
  dump_c_list(fp,lookup_table->fract_seed,FRACT_SEED_SIZE);
  (void) fprintf(fp,"\n};\n");
  return 0;
}

int write_table(char *fn,char *table,t_lutab *lookup_table,
                int formatter(FILE *fp,char *fn,char *table,
                              t_lutab *lookup_table))
{
  int done;
  FILE *fp;
  
  if ((fp=fopen(fn,"w"))==NULL)
    done=-1;
  else
    {
      if ((done=formatter(fp,fn,table,lookup_table))!=-1)
        done=fclose(fp);
      else
        (void) fclose(fp);
    }
  return done;
}

int table_test(FILE *log)
{
  int i,errcount,exp_mask;
  float x,y0,y1,yd;
  t_convert convert;
  
  errcount=0;
  init_table(&lookup_table);
  convert.fval=1.0;
  exp_mask=convert.bval;
  for (i=1; i<FRACT_SEED_SIZE; i++)
    {
      convert.bval=(i<<FRACT_SHIFT)|(exp_mask);
      x=convert.fval;
      y0=1.0/sqrt(x);
      y1=invsqrt(x);
      yd=fabs(y0-y1);
      
      if (yd!=0.0) {
	errcount++;
	fprintf(log,"x=%e, y0=%e, y1=%e, yd=%e\n",x,y0,y1,yd);
      }
    }
  return errcount;
}

int main(int argc)
{
  int done;
  
  init_table(&lookup_table);
  if (argc>1)
    done=table_test(stdout);
  /*done=write_table("lutab.i","lutab",&lookup_table,write_asm_table);*/
  done=write_table("sqrtab.c","lookup_table",&lookup_table,write_c_table);
  return 0;
}  
