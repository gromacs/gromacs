/***************************************************************
**                                                            **
**        D I F F E R E N T I A L     E V O L U T I O N       **
**                                                            **
** Program: de.c                                              **
** Version: 3.6                                               **
**                                                            **
** Authors: Dr. Rainer Storn                                  **
**          c/o ICSI, 1947 Center Street, Suite 600           **
**          Berkeley, CA 94707                                **
**          Tel.:   510-642-4274 (extension 192)              **
**          Fax.:   510-643-7684                              **
**          E-mail: storn@icsi.berkeley.edu                   **
**          WWW: http://http.icsi.berkeley.edu/~storn/        **
**          on leave from                                     **
**          Siemens AG, ZFE T SN 2, Otto-Hahn Ring 6          **
**          D-81739 Muenchen, Germany                         **
**          Tel:    636-40502                                 **
**          Fax:    636-44577                                 **
**          E-mail: rainer.storn@zfe.siemens.de               **
**                                                            **
**          Kenneth Price                                     **
**          836 Owl Circle                                    **
**          Vacaville, CA 95687                               **
**          E-mail: kprice@solano.community.net               ** 
**                                                            **
** This program implements some variants of Differential      **
** Evolution (DE) as described in part in the techreport      **
** tr-95-012.ps of ICSI. You can get this report either via   **
** ftp.icsi.berkeley.edu/pub/techreports/1995/tr-95-012.ps.Z  **
** or via WWW: http://http.icsi.berkeley.edu/~storn/litera.html*
** A more extended version of tr-95-012.ps is submitted for   **
** publication in the Journal Evolutionary Computation.       ** 
**                                                            **
** You may use this program for any purpose, give it to any   **
** person or change it according to your needs as long as you **
** are referring to Rainer Storn and Ken Price as the origi-  **
** nators of the the DE idea.                                 **
** If you have questions concerning DE feel free to contact   **
** us. We also will be happy to know about your experiences   **
** with DE and your suggestions of improvement.               **
**                                                            **
***************************************************************/
/**H*O*C**************************************************************
**                                                                  **
** No.!Version! Date ! Request !    Modification           ! Author **
** ---+-------+------+---------+---------------------------+------- **
**  1 + 3.1  +5/18/95+   -     + strategy DE/rand-to-best/1+  Storn **
**    +      +       +         + included                  +        **
**  1 + 3.2  +6/06/95+C.Fleiner+ change loops into memcpy  +  Storn **
**  2 + 3.2  +6/06/95+   -     + update comments           +  Storn **
**  1 + 3.3  +6/15/95+ K.Price + strategy DE/best/2 incl.  +  Storn **
**  2 + 3.3  +6/16/95+   -     + comments and beautifying  +  Storn **
**  3 + 3.3  +7/13/95+   -     + upper and lower bound for +  Storn **
**    +      +       +         + initialization            +        **
**  1 + 3.4  +2/12/96+   -     + increased printout prec.  +  Storn **
**  1 + 3.5  +5/28/96+   -     + strategies revisited      +  Storn **
**  2 + 3.5  +5/28/96+   -     + strategy DE/rand/2 incl.  +  Storn **
**  1 + 3.6  +8/06/96+ K.Price + Binomial Crossover added  +  Storn **
**  2 + 3.6  +9/30/96+ K.Price + cost variance output      +  Storn **
**  3 + 3.6  +9/30/96+   -     + alternative to ASSIGND    +  Storn **
**  4 + 3.6  +10/1/96+   -    + variable checking inserted +  Storn **
**  5 + 3.6  +10/1/96+   -     + strategy indic. improved  +  Storn **
**                                                                  **
***H*O*C*E***********************************************************/

/* Adopted for use in GROMACS by David van der Spoel, Oct 2001 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "typedefs.h"
#include "smalloc.h"
#include "futil.h"
#include "genalg.h"
#include "fatal.h"
#include "random.h"
#include "vec.h"
#include "main.h"

static char  *strat[] = {
  "DE/best/1/exp",          "DE/rand/1/exp",
  "DE/rand-to-best/1/exp",  "DE/best/2/exp",
  "DE/rand/2/exp",          "DE/best/1/bin",
  "DE/rand/1/bin",          "DE/rand-to-best/1/bin",
  "DE/best/2/bin",          "DE/rand/2/bin"
};

/*------------------------Macros----------------------------------------*/
#define assignd(D,a,b) memcpy((a),(b),sizeof(a[0])*D)   
/* quick copy by Claudio, works only for small arrays, but is faster */

static real **make2d(int n,int m)
{
  int  i;
  real **r;
  
  snew(r,n);
  for(i=0;(i<n); i++)
    snew(r[i],m);
  return r;
}

static void copy2range(int D,real c[],t_range r[])
{
  int i;
  
  for(i=0; (i<D); i++) {
    /* Range check */
    while (c[i] < r[i].rmin)
      c[i] += r[i].dr;
    while (c[i] > r[i].rmax)
      c[i] -= r[i].dr;
    /*    if (c[i] < r[i].rmin)
	  c[i] = r[i].rmin;
	  if (c[i] > r[i].rmax)
	  c[i] = r[i].rmax;
    */
    r[i].rval = c[i];
  }
}

t_genalg *init_ga(char *infile,int D,t_range range[])
{
  FILE     *fpin_ptr;
  t_genalg *ga;
  double   ff,cr;
  int      i,j;
  
  /*------Initializations----------------------------*/
  snew(ga,1);
  /*-----Read input data------------------------------------------------*/
  fpin_ptr   = ffopen(infile,"r");
  fscanf(fpin_ptr,"%d",&ga->NP);             /*---choice of strategy-----------------*/
  fscanf(fpin_ptr,"%d",&ga->strategy);       /*---choice of strategy-----------------*/
  fscanf(fpin_ptr,"%lf",&ff);            /*---weight factor----------------------*/
  fscanf(fpin_ptr,"%lf",&cr);            /*---crossing over factor---------------*/
  fscanf(fpin_ptr,"%d",&ga->seed);           /*---random seed------------------------*/
  fclose(fpin_ptr);
  
  ga->FF   = ff;
  ga->CR   = cr;
  ga->D    = D;
  ga->ipop = 0;
  ga->gen  = 0;
  
  /* Allocate memory */
  ga->pold = make2d(ga->NP,ga->D);
  ga->pnew = make2d(ga->NP,ga->D);
  snew(ga->tmp,ga->D);
  snew(ga->best,ga->D);
  snew(ga->bestit,ga->D);
  snew(ga->cost,ga->NP);
  snew(ga->rmsf,ga->NP);
  snew(ga->energy,ga->NP);
    
  /*-----Checking input variables for proper range--------------*/
  if ((ga->CR < 0) || (ga->CR > 1.0)) 
    fatal_error(0,"CR=%f, should be ex [0,1]",ga->CR);
  if (ga->seed <= 0) 
    fatal_error(0,"seed=%d, should be > 0",ga->seed);
  if ((ga->strategy < 0) || (ga->strategy > 10)) 
    fatal_error(0,"strategy=%d, should be ex {1-10}",ga->strategy);

  /* spread initial population members */
  for (i=0; (i<ga->NP); i++) {
    for (j=0; (j<ga->D); j++) {
      ga->pold[i][j] = value_rand(&(range[j]),&ga->seed);
    }
  }

  fprintf(stdlog,"-----------------------------------------------\n");  
  fprintf(stdlog,"Genetic algorithm parameters\n");
  fprintf(stdlog,"-----------------------------------------------\n");  
  fprintf(stdlog,"Population size:       %d\n",ga->NP);   
  fprintf(stdlog,"Strategy:              %s\n",strat[ga->strategy]);
  fprintf(stdlog,"Weight factor:         %g\n",ga->FF);
  fprintf(stdlog,"Crossing over factor:  %g\n",ga->CR);
  fprintf(stdlog,"Random seed:           %d\n",ga->seed);
  fprintf(stdlog,"-----------------------------------------------\n");  
  
  return ga;
}

void update_ga(FILE *fpout_ptr,t_range range[],t_genalg *ga)
{
  static int  i_init=0;   /* Initialisation related stuff       */
  int   i, j, L, n;      /* counting variables                 */
  int   r1,r2,r3,r4,r5;  /* placeholders for random indexes    */

  if (i_init < ga->NP) {
    /* Copy data for first force evaluation to range array  */
    copy2range(ga->D,ga->pold[i_init],range);
    
    i_init++;
    return;
  }
  else {
    /* Now starts real genetic stuff, a new trial set is made */
    if (ga->ipop == ga->NP) {
      ga->gen++;
      i=ga->ipop=0;
    }
    else 
      i=ga->ipop;

    do {                        /* Pick a random population member */
      /* Endless loop for ga->NP < 2 !!!     */
      r1 = (int)(rando(&ga->seed)*ga->NP);
    } while(r1==i);            
      
    do {                        /* Pick a random population member */
      /* Endless loop for ga->NP < 3 !!!     */
      r2 = (int)(rando(&ga->seed)*ga->NP);
      } while((r2==i) || (r2==r1));
    
    do {  
      /* Pick a random population member */
      /* Endless loop for ga->NP < 4 !!!     */
      r3 = (int)(rando(&ga->seed)*ga->NP);
    } while((r3==i) || (r3==r1) || (r3==r2));
    
    do {
      /* Pick a random population member */
      /* Endless loop for ga->NP < 5 !!!     */
      r4 = (int)(rando(&ga->seed)*ga->NP);
    } while((r4==i) || (r4==r1) || (r4==r2) || (r4==r3));
    
    do {
      /* Pick a random population member */
      /* Endless loop for ga->NP < 6 !!!     */
      r5 = (int)(rando(&ga->seed)*ga->NP);
    } while((r5==i) || (r5==r1) || (r5==r2) || (r5==r3) || (r5==r4));
    
    
    /* Choice of strategy
     * We have tried to come up with a sensible naming-convention: DE/x/y/z
     * DE :  stands for Differential Evolution
     * x  :  a string which denotes the vector to be perturbed
     * y  :  number of difference vectors taken for perturbation of x
     * z  :  crossover method (exp = exponential, bin = binomial)
     *
     * There are some simple rules which are worth following:
     * 1)  ga->FF is usually between 0.5 and 1 (in rare cases > 1)
     * 2)  ga->CR is between 0 and 1 with 0., 0.3, 0.7 and 1. being worth to 
     *     be tried first
     * 3)  To start off ga->NP = 10*ga->D is a reasonable choice. Increase ga->NP if 
     *     misconvergence happens.
     * 4)  If you increase ga->NP, ga->FF usually has to be decreased
     * 5)  When the DE/ga->best... schemes fail DE/rand... usually works and 
     *     vice versa
     * EXPONENTIAL ga->CROSSOVER
     *-------DE/ga->best/1/exp-------
     *-------Our oldest strategy but still not bad. However, we have found several
     *-------optimization problems where misconvergence occurs.
     */
    assignd(ga->D,ga->tmp,ga->pold[i]);
    n = (int)(rando(&ga->seed)*ga->D);
    L = 0;
    
    switch (ga->strategy) {
    case 1:
      /* strategy DE0 (not in our paper) */
      do {                       
	ga->tmp[n] = ga->bestit[n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
	n = (n+1)%ga->D;
	L++;
      } while((rando(&ga->seed) < ga->CR) && (L < ga->D));
      break;
      
      /* DE/rand/1/exp
       * This is one of my favourite strategies. It works especially 
       * well when the ga->bestit[]"-schemes experience misconvergence. 
       * Try e.g. ga->FF=0.7 and ga->CR=0.5 * as a first guess.
       */
    case 2:
      /* strategy DE1 in the techreport */
      do {                       
	ga->tmp[n] = ga->pold[r1][n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
	n = (n+1)%ga->D;
	L++;
      } while((rando(&ga->seed) < ga->CR) && (L < ga->D));
      break;
      
      /* DE/rand-to-ga->best/1/exp 
       * This strategy seems to be one of the ga->best strategies. 
       * Try ga->FF=0.85 and ga->CR=1.
       * If you get misconvergence try to increase ga->NP. 
       * If this doesn't help you should play around with all three 
       * control variables.
       */
    case 3:
      /* similar to DE2 but generally better */
      do {                       
	ga->tmp[n] = ga->tmp[n] + ga->FF*(ga->bestit[n] - ga->tmp[n]) + 
	  ga->FF*(ga->pold[r1][n]-ga->pold[r2][n]);
	n = (n+1)%ga->D;
	L++;
      } while((rando(&ga->seed) < ga->CR) && (L < ga->D));
      break;
      
      /* DE/ga->best/2/exp is another powerful strategy worth trying */
    case 4:
      do {                           
	ga->tmp[n] = ga->bestit[n] + 
	  (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
	n = (n+1)%ga->D;
	L++;
      } while((rando(&ga->seed) < ga->CR) && (L < ga->D));
      break;
      
      /*----DE/rand/2/exp seems to be a robust optimizer for many functions-----*/
    case 5:
      do {                           
	ga->tmp[n] = ga->pold[r5][n] + 
	  (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
	n = (n+1)%ga->D;
	L++;
      } while((rando(&ga->seed) < ga->CR) && (L < ga->D));
      break;
      
      /*===Essentially same strategies but BINOMIAL ga->CROSSOVER===*/
      
      /*-------DE/ga->best/1/bin------*/
    case 6:
      for (L=0; L<ga->D; L++) {
	/* perform D binomial trials */
	if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1))) {
	  /* change at least one parameter */
	    ga->tmp[n] = ga->bestit[n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
	}
	n = (n+1)%ga->D;
      }
      break;
      
      /*-------DE/rand/1/bin------*/
    case 7:
      for (L=0; L<ga->D; L++) {
	/* perform D binomial trials */
	if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1))) {
	  /* change at least one parameter */
	  ga->tmp[n] = ga->pold[r1][n] + ga->FF*(ga->pold[r2][n]-ga->pold[r3][n]);
	}
	n = (n+1)%ga->D;
      }
      break;
      
      /*-------DE/rand-to-ga->best/1/bin------*/
    case 8:
      for (L=0; L<ga->D; L++) {
	/* perform ga->D binomial trials */
	if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1))) {
	  /* change at least one parameter */
	  ga->tmp[n] = ga->tmp[n] + ga->FF*(ga->bestit[n] - ga->tmp[n]) + 
	    ga->FF*(ga->pold[r1][n]-ga->pold[r2][n]);
	}
	n = (n+1)%ga->D;
      }
      break;
      
      /*-------DE/ga->best/2/bin------*/
    case 9:
      for (L=0; L<ga->D; L++) {
	/* perform ga->D binomial trials */
	if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1))) {
	  /* change at least one parameter */
	  ga->tmp[n] = ga->bestit[n] + 
	    (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
	}
	n = (n+1)%ga->D;
      }
      break;
      
      /*-------DE/rand/2/bin-------*/
    default:
      for (L=0; L<ga->D; L++) {
	/* perform ga->D binomial trials */
	if ((rando(&ga->seed) < ga->CR) || (L == (ga->D-1))) {
	  /* change at least one parameter */
	  ga->tmp[n] = ga->pold[r5][n] + 
	    (ga->pold[r1][n]+ga->pold[r2][n]-ga->pold[r3][n]-ga->pold[r4][n])*ga->FF;
	}
	n = (n+1)%ga->D;
      }
      break;
    }
    
    /*===Trial mutation now in ga->tmp[]. Test how good this choice really was.==*/
    copy2range(ga->D,ga->tmp,range);
  }     
}

bool print_ga(FILE *fp,t_genalg *ga,real rmsf,real energy,t_range range[],
	      real tol)
{
  static int nfeval=0;          /* number of function evaluations     */
  static bool bImproved;
  real trial_cost;
  real cvar;            /* computes the cost variance         */
  real cmean;           /* mean cost                          */
  int  i,j;
  real **pswap;
  
  /* When we get here we have done an initial evaluation for all
   * animals in the population
   */
  trial_cost = cost(rmsf,energy);
  if (nfeval < ga->NP) {
    ga->cost[nfeval]   = trial_cost;
    ga->rmsf[nfeval]   = rmsf;
    ga->energy[nfeval] = energy;
    nfeval++;
    return FALSE;
  }
  if (ga->ipop == 0)
    bImproved = FALSE;
    
  /* First iteration after first round of trials */
  if (nfeval == ga->NP) {
    /* Evaluate who is ga->best */  
    ga->imin = 0;
    for (j=1; (j<ga->NP); j++) {
      if (ga->cost[j] < ga->cost[ga->imin]) 
	ga->imin = j;
    }
    assignd(ga->D,ga->best,ga->pold[ga->imin]);   /* save best member ever          */
    assignd(ga->D,ga->bestit,ga->pold[ga->imin]); /* save best member of generation */
  }

  if (trial_cost <= ga->cost[ga->ipop]) {
    /* improved objective function value ? */
    ga->cost[ga->ipop]   = trial_cost;         
    ga->rmsf[ga->ipop]   = rmsf;         
    ga->energy[ga->ipop] = energy;         
    assignd(ga->D,ga->pnew[ga->ipop],ga->tmp);
    if (trial_cost < ga->cost[ga->imin]) {
      /* Was this a new minimum? */
      /* if so, reset cmin to new low...*/
      ga->imin = ga->ipop;
      assignd(ga->D,ga->best,ga->tmp);
      bImproved = TRUE;
    }                           
  }                            
  else {
    assignd(ga->D,ga->pnew[ga->ipop],ga->pold[ga->ipop]); 
    /* replace target with old value */
  }
  /* Increase population member count */
  ga->ipop++;
  
  if (ga->ipop == ga->NP) {
    /* End mutation loop through pop. */
    assignd(ga->D,ga->bestit,ga->best);  
    /* Save ga->best population member of current iteration */
    
    /* swap population arrays. New generation becomes old one */
    pswap = ga->pold;
    ga->pold  = ga->pnew;
    ga->pnew  = pswap;
    
    /*----Compute the energy variance (just for monitoring purposes)-----------*/
    /* compute the mean value first */
    cmean = 0.;          
    for (j=0; (j<ga->NP); j++) 
      cmean += ga->cost[j];
    cmean = cmean/ga->NP;
    
    /* now the variance              */
    cvar = 0.;   
    for (j=0; (j<ga->NP); j++) 
      cvar += sqr(ga->cost[j] - cmean);
    cvar = cvar/(ga->NP-1);
    
    /*----Output part----------------------------------------------------------*/
    if (1 || bImproved || (nfeval == ga->NP)) {
      fprintf(fp,"\nGen: %6d Cost:%8.3f  Ener: %8.3f RMSF: %8.3f  <Cost>: %8.3f\n",
	      ga->gen,ga->cost[ga->imin],
	      ga->energy[ga->imin],ga->rmsf[ga->imin],cmean);
      
      for (j=0; (j<ga->D); j++)
	fprintf(fp,"\tbest[%d]=%-15.10g\n",j,ga->best[j]);

      if (ga->cost[ga->imin] < tol) {
	for (i=0; (i<ga->NP); i++) {
	  fprintf(fp,"Animal: %3d Cost:%8.3f  Ener: %8.3f RMSF: %8.3f%s\n",
		  i,ga->cost[i],ga->energy[i],ga->rmsf[i],
		  (i == ga->imin) ? " ***" : "");
	  for (j=0; (j<ga->D); j++)
	    fprintf(fp,"\tParam[%d][%d]=%15.10g\n",i,j,ga->pold[i][j]);
	}
	return TRUE;
      }
      fflush(fp);
    }
  }
  nfeval++;
  
  return FALSE;
}

