#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
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
  { "Run", "Domain decomp.", "Vsite constr.", "Send X to PME", "Comm. coord.", "Neighbor search", "Force", "Wait + Comm. F", "PME mesh", "PME mesh", "Wait + Comm. X/F", "Wait + Recv. PME F", "Vsite spread", "Write traj.", "Update", "Constraints", "Comm. energies", "Test" };

/* variables for testing/debugging */
static bool              wc_barrier=FALSE;
static gmx_wallcycle_t_t *wc_all=NULL;
static int               wc_depth=0;
static int               ewc_prev=-1;
static gmx_cycles_t      cycle_prev;
#ifdef GMX_MPI
static MPI_Comm          wc_mpi_comm_mygroup;
#endif

bool wallcycle_have_counter(void)
{
  return gmx_cycles_have_counter();
}

gmx_wallcycle_t wallcycle_init(FILE *fplog,t_commrec *cr)
{
  gmx_wallcycle_t_t *wc;

#ifdef GMX_MPI
    if (PAR(cr) && getenv("GMX_CYCLE_BARRIER") != NULL) {
      fprintf(fplog,"\nWill call MPI_Barrier before each cycle start/stop call\n\n");
      wc_barrier = TRUE;
      wc_mpi_comm_mygroup = cr->mpi_comm_mygroup;
    }
#endif

  if (wallcycle_have_counter()) {
    snew(wc,ewcNR);
    if (getenv("GMX_CYCLE_ALL") != NULL) {
      fprintf(fplog,"\nWill time all the code during the run\n\n");
      snew(wc_all,ewcNR*ewcNR);
    }
  } else {
    wc = NULL;
  }

  return wc;
}

static void wallcycle_all_start(int ewc,gmx_cycles_t cycle)
{
  ewc_prev = ewc;
  cycle_prev = cycle;
}

static void wallcycle_all_stop(int ewc,gmx_cycles_t cycle)
{
  wc_all[ewc_prev*ewcNR+ewc].n += 1;
  wc_all[ewc_prev*ewcNR+ewc].c += cycle - cycle_prev;
}

void wallcycle_start(gmx_wallcycle_t wc, int ewc)
{
  gmx_cycles_t cycle;

#ifdef GMX_MPI
  if (wc_barrier)
    MPI_Barrier(wc_mpi_comm_mygroup);
#endif

  if (wc) {
    cycle = gmx_cycles_read();
    wc[ewc].start = cycle;
    if (wc_all) {
      wc_depth++;
      if (ewc == ewcRUN)
	wallcycle_all_start(ewc,cycle);
      else if (wc_depth == 2) {
	wallcycle_all_stop(ewc,cycle);
      }
    }
  }
}

void wallcycle_stop(gmx_wallcycle_t wc, int ewc)
{
  gmx_cycles_t cycle;

#ifdef GMX_MPI
    if (wc_barrier)
      MPI_Barrier(wc_mpi_comm_mygroup);
#endif

  if (wc) {
    cycle = gmx_cycles_read();
    wc[ewc].last = cycle - wc[ewc].start;
    wc[ewc].c += wc[ewc].last;
    wc[ewc].n++;
    if (wc_all) {
      wc_depth--;
      if (ewc == ewcRUN)
	wallcycle_all_stop(ewc,cycle);
      else if (wc_depth == 1)
	wallcycle_all_start(ewc,cycle);
    }
  }
}

double wallcycle_lastcycle(gmx_wallcycle_t wc, int ewc)
{
  double c;
  
  if (wc) {
    c = (double)wc[ewc].last;
    if (ewc == ewcFORCE) {
      /* Remove the PME mesh part from the force count */
      c -= (double)wc[ewcPMEMESH].last;
    }
    if (ewc == ewcUPDATE) {
      /* Remove the constraint part from the force count */
      c -= (double)wc[ewcCONSTR].last;
    }
  } else {
    c = 0;
  }

  return c;
}

void wallcycle_sum(t_commrec *cr, gmx_wallcycle_t wc,double cycles[])
{
  double buf[ewcNR],*cyc_all,*buf_all;
  int    i;

  if (wc) {
    if (wc[ewcPMEMESH_SEP].n > 0) {
      /* This must be a PME only node, calculate the Wait + Comm. time */
      wc[ewcPMEWAITCOMM].c = wc[ewcRUN].c - wc[ewcPMEMESH_SEP].c;
    } else {
      /* Correct the PME mesh only call count */
      wc[ewcPMEMESH_SEP].n = wc[ewcFORCE].n;
      wc[ewcPMEWAITCOMM].n = wc[ewcFORCE].n;
    }

    /* Store the cycles in a double buffer for summing */
    for(i=0; i<ewcNR; i++)
      cycles[i] = (double)wc[i].c;

    /* Remove the PME mesh part from the force count */
    cycles[ewcFORCE] -= cycles[ewcPMEMESH];

    /* Remove the constraint part from the update count */
    cycles[ewcUPDATE] -= cycles[ewcCONSTR];

#ifdef GMX_MPI    
    if (cr->nnodes > 1) {
      MPI_Allreduce(cycles,buf,ewcNR,MPI_DOUBLE,MPI_SUM,cr->mpi_comm_mysim);
      for(i=0; i<ewcNR; i++)
	cycles[i] = buf[i];
      if (wc_all) {
	snew(cyc_all,ewcNR*ewcNR);
	snew(buf_all,ewcNR*ewcNR);
	for(i=0; i<ewcNR*ewcNR; i++)
	  cyc_all[i] = wc_all[i].c;
	MPI_Allreduce(cyc_all,buf_all,ewcNR*ewcNR,MPI_DOUBLE,MPI_SUM,
		      cr->mpi_comm_mysim);
	for(i=0; i<ewcNR*ewcNR; i++)
	  wc_all[i].c = buf_all[i];
	sfree(buf_all);
	sfree(cyc_all);
      }
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
    fprintf(fplog," %-19s %4d %10s %12.3f %10.1f   %5.1f\n",
	    name,nnodes,num,c*1e-9,c*c2t,100*(double)c/(double)tot);
  }
}

void wallcycle_print(FILE *fplog, int nnodes, int npme, double realtime,
		     gmx_wallcycle_t wc, double cycles[])
{
  double c2t,tot,sum;
  int    i,j,npp;
  char   buf[STRLEN];
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

    fprintf(fplog," Computing:         Nodes     Number     G-Cycles    Seconds     %c\n",'%');
    fprintf(fplog,"%s\n",myline);
    sum = 0;
    for(i=ewcRUN+1; i<ewcNR; i++) {
      print_cycles(fplog,c2t,wcn[i],
		   (i==ewcPMEMESH_SEP || i==ewcPMEWAITCOMM) ? npme : npp,
		   wc[i].n,cycles[i],tot);
      sum += cycles[i];
    }
    if (wc_all) {
      for(i=0; i<ewcNR; i++) {
	for(j=0; j<ewcNR; j++) {
	  sprintf(buf,"%-9s",wcn[i]);
	  buf[9] = ' ';
	  sprintf(buf+10,"%-9s",wcn[j]);
	  buf[19] = '\0';
	  print_cycles(fplog,c2t,buf,
		       (i==ewcPMEMESH_SEP || i==ewcPMEWAITCOMM) ? npme : npp,
		       wc_all[i*ewcNR+j].n,wc_all[i*ewcNR+j].c,tot);
	  sum += wc_all[i*ewcNR+j].c;
	}
      }
    }
    print_cycles(fplog,c2t,"Rest",npp,0,tot-sum,tot);
    fprintf(fplog,"%s\n",myline);
    print_cycles(fplog,c2t,"Total",nnodes,0,tot,tot);
    fprintf(fplog,"%s\n",myline);
  }
}
