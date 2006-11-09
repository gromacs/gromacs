#include "gmx_wallcycle.h"
#include "gmx_cyclecounter.h"
#include "smalloc.h"

#ifdef GMX_MPI
#include <mpi.h>
#endif

typedef struct gmx_wallcycle {
  int          n;
  gmx_cycles_t c;
  gmx_cycles_t start;
  gmx_cycles_t last;
} gmx_wallcycle_t_t;

static char *wcn[ewcNR] =
  { "Run", "Domain decomp.", "Neighbor search", "Force", "PME mesh", "PME mesh", "Update" };

bool wallcycle_have_counter(void)
{
  return gmx_cycles_have_counter();
}

gmx_wallcycle_t wallcycle_init(void)
{
  gmx_wallcycle_t_t *wc;

  if (wallcycle_have_counter()) {
    snew(wc,ewcNR);
  } else {
    wc = NULL;
  }

  return wc;
}

void wallcycle_start(gmx_wallcycle_t wc, int ewc)
{
  if (wc) {
    wc[ewc].start = gmx_cycles_read();
  }
}

void wallcycle_stop(gmx_wallcycle_t wc, int ewc)
{
  if (wc) {
    wc[ewc].last = gmx_cycles_read() - wc[ewc].start;
    wc[ewc].c += wc[ewc].last;
    wc[ewc].n++;
  }
}

double wallcycle_lastcycle(gmx_wallcycle_t wc, int ewc)
{
  double c;
  
  if (wc) {
    c = (double)wc[ewc].last;
    if (ewc == ewcFORCE)
      /* Remove the PME mesh part from the force count */
      c -= (double)wc[ewcPMEMESH].last;
  } else {
    c = 0;
  }

  return c;
}

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc,double cycles[])
{
  double buf[ewcNR];
  int    i;

  if (wc) {
    for(i=0; i<ewcNR; i++)
      cycles[i] = (double)wc[i].c;

    /* Remove the PME mesh part from the force count */
    cycles[ewcFORCE] -= cycles[ewcPMEMESH];

    /* Correct the PME mesh only call count */
    wc[ewcPMEMESH_SEP].n = wc[ewcFORCE].n;

#ifdef GMX_MPI    
    if (cr->nnodes > 1) {
      MPI_Allreduce(cycles,buf,ewcNR,MPI_DOUBLE,MPI_SUM,cr->mpi_comm_mysim);
      for(i=0; i<ewcNR; i++)
	cycles[i] = buf[i];
    }
#endif
  }
}

static void print_cycles(FILE *fplog, double c2t, char *name, int nnodes,
			 int n, gmx_cycles_t c, gmx_cycles_t tot)
{
  char num[11];
  
  if (c > 0) {
    if (n > 0)
      sprintf(num,"%10d",n);
    else
      sprintf(num,"          ");
    fprintf(fplog," %-19s %4d %10s %12.3f %9.1f   %5.1f\n",
	    name,nnodes,num,c*1e-9,c*c2t,100*(double)c/(double)tot);
  }
}

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
		     gmx_wallcycle_t wc, double cycles[])
{
  double c2t,tot,sum;
  int    i,npp;
  char   *myline = "-----------------------------------------------------------------------";
  
  if (wc) {
    npp = nnodes - npme;
    tot = cycles[ewcRUN];
    /* Conversion factor from cycles to seconds */
    if (tot > 0)
      c2t = nnodes*realtime/tot;
    else
      c2t = 0;

    fprintf(fplog,"     R E A L   C Y C L E   A N D   T I M E   A C C O U N T I N G\n\n");

    fprintf(fplog," Computing:         Nodes     Number     G-Cycles   Seconds     %c\n",'%');
    fprintf(fplog,"%s\n",myline);
    sum = 0;
    for(i=ewcRUN+1; i<ewcNR; i++) {
      print_cycles(fplog,c2t,wcn[i],i==ewcPMEMESH_SEP ? npme : npp,
		   wc[i].n,cycles[i],tot);
      sum += cycles[i];
    }
    print_cycles(fplog,c2t,"Rest",nnodes,0,tot-sum,tot);
    fprintf(fplog,"%s\n",myline);
    print_cycles(fplog,c2t,"Total",nnodes,0,tot,tot);
    fprintf(fplog,"%s\n",myline);
  }
}
