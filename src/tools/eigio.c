#include "smalloc.h"
#include "vec.h"
#include "eigio.h"
#include "trnio.h"
#include "tpxio.h"

void read_eigenvectors(char *file,int *natoms,bool *bFit,
		       rvec **xref,bool *bDMR,
		       rvec **xav,bool *bDMA,
		       int *nvec, int **eignr, rvec ***eigvec)
{
  t_trnheader head;
  int    status,i,snew_size;
  rvec   *x;
  matrix box;
  bool   bOK;

  *bDMR=FALSE;

  /* read (reference (t=-1) and) average (t=0) structure */
  status=open_trn(file,"r");
  fread_trnheader(status,&head,&bOK);
  *natoms=head.natoms;
  snew(*xav,*natoms);
  fread_htrn(status,&head,box,*xav,NULL,NULL);
  if ((head.t>=-1.1) && (head.t<=-0.9)) {
    snew(*xref,*natoms);
    for(i=0; i<*natoms; i++)
      copy_rvec((*xav)[i],(*xref)[i]);
    *bDMR = (head.lambda > 0.5);
    *bFit = (head.lambda > -0.5);
    if (*bFit)
      fprintf(stderr,"Read %smass weighted reference structure with %d atoms from %s\n", *bDMR ? "" : "non ",*natoms,file);
    else {
       fprintf(stderr,"Eigenvectors in %s were determined without fitting\n",
	       file);
       sfree(*xref);
       *xref=NULL;
    }
    fread_trnheader(status,&head,&bOK);
    fread_htrn(status,&head,box,*xav,NULL,NULL);
  }
  else {
    *bFit=TRUE;
    *xref=NULL;
  }
  *bDMA = (head.lambda > 0.5);
  if ((head.t<=-0.01) || (head.t>=0.01))
    fprintf(stderr,"WARNING: %s does not start with t=0, which should be the "
	    "average structure. This might not be a eigenvector file. "
	    "Some things might go wrong.\n",
	    file);
  else
    fprintf(stderr,
	    "Read %smass weighted average/minimum structure with %d atoms from %s\n",
	    *bDMA ? "" : "non ",*natoms,file);
  
  snew(x,*natoms);
  snew_size=0;
  *nvec=0;
  while (fread_trnheader(status,&head,&bOK)) {
    fread_htrn(status,&head,box,x,NULL,NULL);
    if (*nvec >= snew_size) {
      snew_size+=10;
      srenew(*eignr,snew_size);
      srenew(*eigvec,snew_size);
    }
    i=(int)(head.t+0.01);
    if ((head.t-i<=-0.01) || (head.t-i>=0.01))
      fatal_error(0,"%s contains a frame with non-integer time (%f), this "
		  "time should be an eigenvector index. "
		  "This might not be a eigenvector file.",file,head.t);
    (*eignr)[*nvec]=i-1;
    snew((*eigvec)[*nvec],*natoms);
    for(i=0; i<*natoms; i++)
      copy_rvec(x[i],(*eigvec)[*nvec][i]);
    (*nvec)++;
  }
  sfree(x);
  fprintf(stderr,"Read %d eigenvectors (dim=%d)\n\n",*nvec,*natoms*DIM);
}

void write_eigenvectors(char *trnname,int natoms,real mat[],
			bool bReverse,int begin,int end,
			int WriteXref,rvec *xref,bool bDMR,
			rvec xav[],bool bDMA)
{
  int    trnout;
  int    ndim,i,j,d,vec;
  matrix zerobox;
  rvec   *x;
  
  ndim = natoms*DIM;
  clear_mat(zerobox);
  snew(x,natoms);

  fprintf (stderr,
	   "\nWriting %saverage structure\nand eigenvectors 1 to %d to %s\n",
	   (WriteXref==eWXR_YES) ? "reference and " : "",
	   end,trnname);
  
  trnout = open_tpx(trnname,"w");
  if (WriteXref == eWXR_YES)
    /* misuse lambda: 0/1 mass weighted fit no/yes */  
    fwrite_trn(trnout,-1,-1,bDMR ? 1.0 : 0.0,zerobox,natoms,xref,NULL,NULL);
  else if (WriteXref == eWXR_NOFIT)
    /* misuse lambda: -1 no fit */  
    fwrite_trn(trnout,-1,-1,-1.0,zerobox,natoms,x,NULL,NULL);

  /* misuse lambda: 0/1 mass weighted analysis no/yes */ 
  fwrite_trn(trnout,0,0,bDMA ? 1.0 : 0.0,zerobox,natoms,xav,NULL,NULL);

  for(i=begin; i<=end; i++) {
    if (!bReverse)
      vec = i-1;
    else
      vec = ndim-i;
    for (j=0; j<natoms; j++)
      for(d=0; d<DIM; d++)
	x[j][d]=mat[vec*ndim+DIM*j+d];
    fwrite_trn(trnout,i,(real)i,0,zerobox,natoms,x,NULL,NULL);
  }
  close_trn(trnout);

  sfree(x);
}
