#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This program generates inner-loop source code for the GROMACS
 * Molecular Dynamics package. The innerloops are nearly optimal, i.e.
 * there are no superfluous declarations or statements in the code.
 * The loops are generated in either C or Fortran 77
 * The invsqrt (1/sqrt) routine can be inlined if the defines in the
 * Makefile.$GMXCPU are set accordingly.
 * There are two command line options:
 * First  (0/1) determines C or F77
 * Second (0/1) determines inline invsqrt or not.
 * Currently the program is called from the Makefile with the inline option
 * on always.
 *
 * There are a few more defines that influence the behaviour of mkinl:
 * if compiled with -DUSEVECTOR code which can be vectorized will be generated
 * if compiled with -DHAVE_PRAGMA a Cray compatible pragma will be inserted
 *    before all inner loops that instruct the compiler to ignore dependencies
 *    between variables in the inner loop.
 *    I am curious to hear whether these options make a difference on a 
 *    real vector computer
 *
 * (c) David van der Spoel, May 1999
 */

typedef int bool;

#define FALSE 0
#define TRUE  1

#ifdef CINVSQRT
#define DECREASE_LATENCY
#endif

#ifdef FINVSQRT
#define DECREASE_LATENCY
#endif

#ifndef USEVECTOR
#define PREFETCH 1
#endif

#ifdef _SGI_
static bool bDelayInvsqrt = TRUE;
#else
static bool bDelayInvsqrt = FALSE;
#endif

int CIND = 2;
int FIND = 6;
#define FCON "     &"
#define max(a,b) ((a) > (b)) ? (a) : (b)

/* Global variables */
static FILE *fp    = NULL;
static bool bC     = TRUE;
static bool bLJC   = FALSE;
static bool bEquiv = TRUE; /* Don't make false: generates type errors */
static bool bCoul,bBHAM,bRF,bTab,bWater,bFREE,bInline,bEwald;
static int  prec = 4;
static char buf[256];
static char buf2[256];

char *_arr(char *a,char *n,int one)
{
  char buf[128];
  
  if (bC)
    sprintf(buf,"%s[%s]",a,n);
  else if (one)
    sprintf(buf,"%s(%s)+1",a,n);
  else
    sprintf(buf,"%s(%s)",a,n);
  
  /* Giant memory leak */
  return strdup(buf);
}

#define ARR(a,n)  _arr(#a,#n,0)
#define IARR(a,n) _arr(#a,#n,1)

char *indent(void)
{
  static char indbuf[80];
  int i,n;

  n = max(0,(bC ? CIND : FIND));
  for(i=0; (i<n); i++)
    indbuf[i] = ' ';
  indbuf[i] = '\0';
  
  return indbuf;
}

void p_line(char *s)
{
  fprintf(fp,"%s%s\n",indent(),s);
}

void newline(void)
{
  fprintf(fp,"\n");
}

void comment(char *s)
{
  if (bC)
    fprintf(fp,"\n%s/* %s */\n",indent(),s);
  else {
    FIND--;
    fprintf(fp,"\nC%s%s\n",indent(),s);
    FIND++;
  }
}

void p_real(char *name)
{
  if (bC)
    fprintf(fp,"%sreal %s;\n",indent(),name);
  else 
    fprintf(fp,"%sreal*%1d  %s\n",indent(),prec,name);
}


void p_real4(char *name)
/* lookup seed for invsqrt must be 4 bytes even in double prec. */
{
  if(bC)
    fprintf(fp,"%sfloat %s;\n",indent(),name);
  else
    fprintf(fp,"%sreal*4  %s\n",indent(),name);
}

void p_creal(char *name,double val)
{
  if (bC)
    fprintf(fp,"%sconst real %s = %f;\n",indent(),name,val);
  else 
    fprintf(fp,"%sreal*%1d  %s\n",indent(),prec,name);
}

void p_int(char *name)
{
  if (bC)
    fprintf(fp,"%sint  %s;\n",indent(),name);
  else {
    fprintf(fp,"%sinteger*%d %s\n",indent(),(int)sizeof(int),name);
  }
}

void p_ivdep(void)
{
#ifdef HAVE_PRAGMA
  if (bC)
    fprintf(fp,"#pragma ivdep\n");
  else
    fprintf(fp,"cdir$ivdep\n");
#endif
}

void _p_state(char *left,char *right,char *symb)
{
  if (bC) {
    if (CIND+16+3+strlen(right) > 78) {
      fprintf(fp,"%s%-16s %2s \n",indent(),left,symb);
      CIND+=2;
      fprintf(fp,"%s%s;\n",indent(),right);
      CIND-=2;
    }
    else
      fprintf(fp,"%s%-16s %2s %s;\n",indent(),left,symb,right);
  }
  else {
    if (FIND+16+3+strlen(right) > 72) {
      fprintf(fp,"%s%-16s = \n%s",indent(),left,FCON);
      FIND-=6-3;
      fprintf(fp,"%s%s\n",indent(),right);
      FIND+=6-3;
    }
    else
      fprintf(fp,"%s%-16s = %s\n",indent(),left,right); 
  }
}

#define p_state(l,r) _p_state(l,r,"=")

void p_incr(char *left,char *right)
{
  char sum[256];
  
  sprintf(sum,"%s + %s",left,right);
  p_state(left,sum);
}

void p_decr(char *left,char *right)
{
  char sum[256];
  
  sprintf(sum,"%s - %s",left,right);
  p_state(left,sum);
}

void p_add(char *left,char *r1,char *r2)
{
  sprintf(buf,"%s + %s",r1,r2);
  p_state(left,buf);
}

void p_sub(char *left,char *r1,char *r2)
{
  sprintf(buf,"%s - %s",r1,r2);
  p_state(left,buf);
}

void fseed(void)
{
  comment("Invsqrt variables");
  /*p_line("include 'seed.inc'");*/
  p_int("nexp,nfract,maxfract");
  p_int("fractshift,fractmask,fractf");
  p_int("expshift,expmask,explsb");
  p_int("expseed,fracseed");
  if (!bEquiv || !bInline) {
    p_int("fl2i");
    p_real4("i2fl");
  }
  p_line("parameter (nexp=256,nfract=4096,maxfract=8388607)");
  p_line("parameter (fractshift=12,fractmask=8388607,fractf=1065353216)");
  p_line("parameter (expshift=23,expmask=2139095040,explsb=8388608)");
  p_line("common /seed/ expseed(nexp),fracseed(nfract)");
  p_line("save /seed/");
}

void p_finvsqrt(void)
{
  p_line("function expaddr(val)");
  fseed();
  p_line("integer*4 val");
  p_state("expaddr","rshift(and(val,expmask),expshift)");
  p_line("return");
  p_line("end");
  newline();
      
  p_line("function fractaddr(val)");
  fseed();
  p_line("integer*4 val");
  p_state("fractaddr","rshift(and(val,or(fractmask,explsb)),fractshift)");
  p_line("return");
  p_line("end");
  newline();
      
#ifdef DOUBLE        
  p_line("function invsqrt(x)");
  fseed();
  p_line("real*8    invsqrt,x,y");
  p_line("real*4    lu,xin");
  p_line("integer*4 exp,addr,bval,result");
  newline();
  
  if (bEquiv) {
    p_line("equivalence(xin,bval)");
    p_line("equivalence(lu,result)");
    p_state("xin","x");
  }
  else
    p_state("bval","fl2i(xin)");
  p_state("exp","expaddr(bval)");
  p_state("addr","fractaddr(bval)");
#ifdef DECREASE_LATENCY
  p_state("result","or(expseed(exp+1),fracseed(addr+1))");
  if (!bEquiv)
    p_state("lu","i2fl(result)");
  newline();
  
  p_state("y","(0.5*lu*(3.0-((x*lu)*lu)))");
  p_state("invsqrt","(0.5*y*(3.0-((x*y)*y)))");
  newline();
  
#else
  p_state("result","or(expseed(exp+1),fracseed(addr+1))");
  if (!bEquiv)
    p_state("lu","i2fl(result)");
  newline();
  
  p_state("y","(0.5*lu*(3.0-((x*lu)*lu)))");
  p_state("invsqrt","(0.5*y*(3.0-((x*y)*y)))");
#endif /* DECREASE_LATENCY */
  newline();

  p_line("return");
  p_line("end");
#else /* ! DOUBLE */
  p_line("function invsqrt(x)");
  fseed();
  p_real4("invsqrt,x,y");
#ifdef DECREASE_LATENCY
  if (bWater) {
    p_real4("luO,luH1,luH2,rsqluO,rsqluH1,rsqluH2,tmpO,tmpH1,tmpH2,xin");
    p_line("integer*4 exp,addr,bval,resultrsqO,resultrsqH1,resultrsqH2");
    newline();
    
    if (bEquiv) {
      p_line("equivalence(xin,bval)");
      p_line("equivalence(luO,resultrsqO)");
      p_line("equivalence(luH1,resultrsqH1)");
      p_line("equivalence(luH2,resultrsqH2)");
      p_state("xin","x");
    }
    else
      p_state("bval","fl2i(x)");
  }
  else {
    p_line("real*4    lu,xin");
    p_line("integer*4 exp,addr,bval,result");
    newline();
    
    if (bEquiv) {
      p_line("equivalence(xin,bval)");
      p_line("equivalence(lu,result)");
      p_state("xin","x");
    }
    else
      p_state("bval","fl2i(x)");
  }
#else
  p_line("real*4    lu,xin");
  p_line("integer*4 exp,addr,bval,result");
  newline();
  
  if (bEquiv) {
    p_line("equivalence(xin,bval)");
    p_line("equivalence(lu,result)");
    p_state("xin","x");
  }
  else
    p_state("bval","fl2i(x)");
#endif /* DECREASE_LATENCY */

    
  p_state("exp","expaddr(bval)");
  p_state("addr","fractaddr(bval)");
#ifdef DECREASE_LATENCY
  if (bWater) {
    p_state("resultrsqO","or(expseed(exp+1),fracseed(addr+1))");
    if (!bEquiv)
      p_state("luO","i2fl(resultrsqO)");
    p_state("invsqrt","(0.5*luO*(3.0-((x*luO)*luO)))");

    p_state("resultH1","or(expseed(exp+1),fracseed(addr+1))");
    if (!bEquiv)
      p_state("luH1","i2fl(resultH1)");
    p_state("invsqrt","(0.5*luH1*(3.0-((x*luH1)*luH1)))");

    p_state("resultH2","or(expseed(exp+1),fracseed(addr+1))");
    if (!bEquiv)
      p_state("luH2","i2fl(resultH2)");
    p_state("invsqrt","(0.5*luH2*(3.0-((x*luH2)*luH2)))");
  }
  else {
    p_state("result","or(expseed(exp+1),fracseed(addr+1))");
    if (!bEquiv)
      p_state("lu","i2fl(result)");
    
    p_state("invsqrt","(0.5*lu*(3.0-((x*lu)*lu)))");
  }
#else 
  p_state("result","or(expseed(exp+1),fracseed(addr+1))");
  if (!bEquiv)
    p_state("lu","i2fl(result)");

  p_state("invsqrt","(0.5*lu*(3.0-((x*lu)*lu)))");
#endif /* DECREASE_LATENCY */

  p_line("return");
  p_line("end");
#endif /* DOUBLE */
  if (!bEquiv) {
    p_line("function fl2i(x)");
    comment("convert float to integer (when a float is passed to this routine)");
    p_line("integer*4 x,fl2i");
    p_state("fl2i","x");
    p_line("return");
    p_line("end");
    newline();
    
    p_line("function i2fl(ix)");
    comment("convert integer to float (when an integer is passed to this routine)");
    p_line("real*4 ix,i2fl");
    p_state("i2fl","ix");
    p_line("return");
    p_line("end");
    newline();
  }
}

void p_invsqrt1(char *tmp,char *right)
{
  char rstlbuf[16];

  comment("Doing fast invsqrt");
  if (bC) {
    p_state("bit_pattern.fval",right);
    p_state("exp_addr","EXP_ADDR(bit_pattern.bval)");
    p_state("fract","FRACT_ADDR(bit_pattern.bval)");
    p_state("result.bval","lookup_table.exp_seed[exp_addr] | lookup_table.fract_seed[fract]");
    p_state(tmp,"result.fval");
  }
  else {
    p_state("xin",right);
    p_state("iexp","rshift(and(bval,expmask),expshift)");
    p_state("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
#ifdef DECREASE_LATENCY
    if (bWater) {
      sprintf(rstlbuf,"result%s",right);
      p_state(rstlbuf,"or(expseed(iexp+1),fracseed(addr+1))");
    }
    else {
      p_state("result","or(expseed(iexp+1),fracseed(addr+1))");
    }
#else
    p_state("result","or(expseed(iexp+1),fracseed(addr+1))");
#endif /* DECREASE_LATENCY */
  }
}

void p_invsqrt2(char *left,char *tmp,char *right)
{
  /* Language independent! */
  sprintf(buf,"(half*%s*(three-((%s*lu)*lu)))",tmp,right);
#ifdef DOUBLE
  p_state("y1",buf);
  sprintf(buf,"y2=(half*y1*(three-((%s*y1)*y1)))",right);
  p_state(left,buf);
#else
  p_state(left,buf);
#endif /* DOUBLE */
}

int p_invsqrt(char *left,char *right)
{
  bool bInvsqrt = FALSE;
  
#ifdef CINVSQRT
  if (bC)
    bInvsqrt = TRUE;
#endif
#ifdef FINVSQRT
  if (!bC)
    bInvsqrt = TRUE;
#endif

  if (bInvsqrt) {
    if (bInline) {
      p_invsqrt1("lu",right);
      p_invsqrt2(left,"lu",right);
    }
    else {
      sprintf(buf,"invsqrt(%s)",right);
      p_state(left,buf);
    }
  }
  else {
    sprintf(buf,"1.0/sqrt(%s)",right);
#ifdef HAVE_SQRTF
#ifndef DOUBLE
    if (bC)
      sprintf(buf,"1.0/sqrtf(%s)",right);
#endif
#endif
    p_state(left,buf);
  }
#ifdef DOUBLE
  return 10;
#else
  return 5;
#endif
}

int p_waterinvsqrt(void)
{
  char *w[] = { "O", "H1", "H2" };
  char vbuf[16],rbuf[32];
  int  m;
  bool bInvsqrt = FALSE;
  
  if (bC)
#ifdef CINVSQRT
    bInvsqrt = TRUE
#endif
      ;
  else
#ifdef FINVSQRT
    bInvsqrt = TRUE
#endif
      ;
  
  if (bInvsqrt) {
    if (bInline) {
      comment("Water invsqrt optimization");
      for(m=0; (m<3); m++) {
	sprintf(vbuf,"lu%s",w[m]);
	sprintf(rbuf,"rsq%s",w[m]);
	p_invsqrt1(vbuf,rbuf);
      }
      for(m=0; (m<3); m++) {
	sprintf(vbuf,"rsqlu%s",w[m]);
	sprintf(rbuf,"rsq%s*lu%s",w[m],w[m]);
	p_state(vbuf,rbuf);
      }
      for(m=0; (m<3); m++) {
	sprintf(vbuf,"tmp%s",w[m]);
	sprintf(rbuf,"rsqlu%s*lu%s",w[m],w[m]);
	p_state(vbuf,rbuf);
      }
      for(m=0; (m<3); m++) {
	sprintf(vbuf,"rinv1%s",w[m]);
	sprintf(rbuf,"half*lu%s*(three-tmp%s)",w[m],w[m]);
	p_state(vbuf,rbuf);
      }
    }
    else {
      for(m=0; (m<3); m++) {
	sprintf(vbuf,"rinv1%s",w[m]);
	sprintf(rbuf,"invsqrt(rsq%s)",w[m]);
	p_state(vbuf,rbuf);
      }
    }
  } 
  else {
    for(m=0; (m<3); m++) {
      sprintf(vbuf,"rinv1%s",w[m]);
      sprintf(rbuf,"1.0/sqrt(rsq%s)",w[m]);
#ifdef HAVE_SQRTF
#ifndef DOUBLE
      if (bC)
	sprintf(rbuf,"1.0/sqrtf(rsq%s)",w[m]);
#endif
#endif
      p_state(vbuf,rbuf);
    }
  }
#ifdef DOUBLE
  return 10;
#else
  return 5;
#endif
}

void fheader(void)
{
  int oprec;
  
  if (bC) 
    fprintf(fp,"#include <stdio.h>\n"
	    "#include <math.h>\n"
	    "#include \"typedefs.h\"\n"
	    "#include \"lutab.h\"\n" 
	    /* "#include \"fatal.h\"\n" */
	    "#include \"inner.h\"\n\n");
  else {
#ifdef FINVSQRT
    fprintf(fp,"%ssubroutine ffillbuf\n",indent());
    comment("This subroutine initiates the buffer for invsqrt computation");
    fseed();
    p_int("i,exp,newexp,bval,addr,indx,findex");
    oprec = prec;
    prec  = 4;
    p_real4("fval,bflt");
    prec  = oprec;
    if (bEquiv) {
      p_line("equivalence(fval,findex)");
      p_line("equivalence(bflt,bval)");
    }
    p_line("do i=1,nexp");
    FIND += 2;
    p_state("expseed(i)","-1");
    FIND -= 2;
    p_line("end do");
    p_line("do i=1,nfract");
    FIND += 2;
    p_state("fracseed(i)","-1");
    FIND-=2;
    p_line("end do");
    p_line("do i=0,nexp-1");
    FIND += 2;
    p_state("exp","(i-127)");
    p_line("if (exp .lt. 0) then");
    FIND += 2;
    p_line("newexp=127-((exp+1)/2)");
    FIND -= 2;
    p_line("else");
    FIND += 2;
    p_state("newexp","127-(exp/2)-1");
    FIND -= 2;
    p_line("endif");
    p_state("expseed(i+1)","and(lshift(newexp,expshift),expmask)");
    FIND -= 2;
    p_line("end do");

    p_state("indx","rshift(fractf,fractshift)");
    p_line("do i=0,nfract-1");
    FIND += 2;
    p_state("bval","lshift(indx,fractshift)");
    p_state("addr","rshift(and(bval,or(fractmask,explsb)),fractshift)");
    if (bEquiv) {
      p_state("fval","(1.0/sqrt(bflt))");
      p_state("fracseed(addr+1)","and(findex,fractmask)");
    }
    else {
      p_state("fval","(1.0/sqrt(i2fl(bval)))");
      p_state("fracseed(addr+1)","and(fl2i(fval),fractmask)");
    }
    p_line("if (fracseed(addr+1) .eq. 0) then");
    FIND += 2;
    p_state("fracseed(addr+1)","and(maxfract,fractmask)");
    FIND -= 2;
    p_line("endif");
    p_state("indx","indx+1");
    FIND -= 2;
    p_line("end do");
      
    p_line("return");
    p_line("end");
    fprintf(fp,"\n");
    
    if (!bEquiv || !bInline)
      p_finvsqrt();
#else
    ;
#endif 
  }
}

void fname(char *name)
{
  if (bC)
    fprintf(fp,"void %s(\n",name);
  else
    fprintf(fp,"%ssubroutine %s(\n",indent(),name);
}

void fwarning(void)
{
  fprintf(fp,bC ?
	  "  /**********************************************************\n" 
	  "   * This function was generated automatically by mkinl.\n"
	  "   *                     DO NOT EDIT\n"
	  "   *           Copyright David van der Spoel 1999\n"
	  "   **********************************************************/"
	  "\n" :
	  "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n"
	  "C     This function was generated automatically by mkinl.\n" 
	  "C                        DO NOT EDIT\n"
	  "C              Copyright David van der Spoel 1999\n"
	  "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n");
}

void fcall_args(void)
{
  if (bBHAM || !bCoul)
    fprintf(fp,"%s",bC ? 
	    "\tntype,type,nbfp,Vnb,\n\t" :
	    FCON "ntype,type,nbfp,vnb,\n");
  if (bRF)
    fprintf(fp,"%s",bC ? "krf,\n\t" : FCON "   krf,\n");
  if (bEwald)
    fprintf(fp,"%s",bC ? "ewc,\n\t" : FCON "   ewc,\n");
  if (bTab)
    fprintf(fp,"%s",bC ? 
	    "tabscale,VFtab,\n\t" :
	    FCON "  tabscale,vftab,\n");
  if (bBHAM && bTab)
    fprintf(fp,"%s",bC ? 
	    "exptabscale,\n\t" : 
	    FCON "  exptabscale,\n");
  if (bFREE) 
    fprintf(fp,"%s",bC ?
	    "chargeB,typeB,lambda,dvdlambda,\n\t" :
	    FCON "  chargeB,typeB,lambda,dvdlambda,\n");
#ifdef USEVECTOR
  fprintf(fp,"%s",bC ? "fbuf,\n\t" : FCON "  fbuf,\n");
#endif	
  /* Petty formatting... */
  if (bC && bCoul && !bTab && !bRF)
    fprintf(fp,"\t");
    
  fprintf(fp,"%s",bC ?
	  "nri,iinr,shift,gid,jindex,jjnr,pos,fshift,\n\t"
	  "facel,charge,faction,Vc,shiftvec);\n" :
	  FCON "  nri,iinr,shift,gid,jindex,jjnr,pos,fshift,\n"
	  FCON "  facel,charge,faction,Vc,shiftvec)\n");
}

void fargs(void)
{
  if (bBHAM || !bCoul)
    fprintf(fp,"%s",bC ? 
	    "\tint  ntype,"
	    "\n\tint  type[],\n\treal nbfp[],\n\treal Vnb[],\n\t" :
	    FCON "  ntype,\n"
	    FCON "  type,\n" FCON "  nbfp,\n" FCON "  vnb,\n");
  if (bRF)
    fprintf(fp,"%s",bC ? "real krf,\n\t" : FCON "   krf,\n");
  if (bEwald)
    fprintf(fp,"%s",bC ? "real ewc,\n\t" : FCON "   ewc,\n");    

  if (bTab)
    fprintf(fp,"%s",bC ? 
	    "real tabscale,\n\treal VFtab[],\n\t" :
	    FCON "  tabscale,\n" FCON "  vftab,\n");
  if (bBHAM && bTab)
    fprintf(fp,"%s",bC ? 
	    "real exptabscale,\n\t" : 
	    FCON "  exptabscale,\n");
  if (bFREE) 
    fprintf(fp,"%s",bC ?
	    "real chargeB[],\n\tint  typeB[],\n\t"
	    "real lambda,\n\treal *dvdlambda,\n\t" :
	    FCON "  chargeB,\n" FCON "  typeB,\n"
	    FCON "  lambda,\n"  FCON "  dvdlambda,\n");
#ifdef USEVECTOR
  fprintf(fp,"%s",bC ? "real fbuf[],\n\t" : FCON "  fbuf,\n");
#endif	
  /* Petty formatting... */
  if (bC && bCoul && !bTab && !bRF)
    fprintf(fp,"\t");
    
  fprintf(fp,"%s",bC ?
	  "int  nri,\n\tint  iinr[],\n\tint  shift[],\n\t"
	  "int  gid[],\n\tint  jindex[],\n\tint  jjnr[],\n\t"
	  "real pos[],\n\treal fshift[],\n\t"
	  "real facel,\n\treal charge[],\n\t"
	  "real faction[],\n\treal Vc[],\n\treal shiftvec[])\n{\n" :
	  FCON "  nri,\n"     FCON "  iinr,\n"   FCON "  shift,\n"    
	  FCON "  gid,\n"     FCON "  jindex,\n" FCON "  jjnr,\n"    
	  FCON "  pos,\n"     FCON "  fshift,\n"
	  FCON "  facel,\n"   FCON "  charge,\n"
	  FCON "  faction,\n" FCON "  Vc,\n"    
	  FCON "  shiftvec)\n");
  if (!bC) {
    fprintf(fp,"\n%simplicit none\n",indent());
    comment("Parameters passed to this function");
    if (bBHAM || !bCoul) {
      comment("Van der Waals parameters");
      p_int("ntype,type(*)");
      p_real("nbfp(*),vnb(*)");
    }
    if (bRF) {
      comment("Reaction field parameters");
      p_real("krf");
    }
    if (bEwald) {
      comment("Ewald coefficient");
      p_real("ewc");
    }
    if (bTab) {
      comment("Table lookup parameters");
      p_real("tabscale,vftab(*)");
      if (bBHAM)
	p_real("exptabscale");
    }
    if (bFREE) {
      comment("Free energy parameters");
      p_real("chargeB(*),lambda,dvdlambda");
      p_int("typeB(*)");
    }
#ifdef USEVECTOR
    comment("Buffer for vectorization");
    p_real("fbuf(*)");
#endif
    comment("General parameters");
    p_int("nri,iinr(*),jindex(*),jjnr(*),shift(*),gid(*)");
    p_real("pos(*),facel,charge(*),faction(*),vc(*)");
    p_real("shiftvec(*),fshift(*)");
    fprintf(fp,"\n");
  }
  fwarning();
}

void flocal(void)
{
  fprintf(fp,"\n");
  comment("Local variables");
  p_creal("nul",0.0);
  if (bLJC) {
    comment("Lennard Jones");
    p_int("tjA,ntiA");
    p_real("vnb6,vnb12,vnbtot");
  }
  if (bRF) {
    comment("Reaction field");
    p_creal("two",2.0);
    p_real("krsqO");
#ifndef SIMPLEWATER
    if (bWater)
      p_real("krsqH1,krsqH2");
#endif
  }
  if (bEwald) {
      comment("Ewald correction");
      p_creal("isp",0.564189583547756);
      if (bLJC || bBHAM)
	p_real("ewtmp"); /* needed to break long lines into two */
  }
  if (bTab) {
    comment("Table stuff");
    p_creal("one",1.0);
    p_creal("two",2.0);
    p_real("r1,r1t,fijC,eps,eps2,Y,F,Fp,Geps,Heps2,VV,FF");
    if (bBHAM || bLJC) {
      p_real("c6,fijD,fijR");
    }
    if (bBHAM) {
      comment("Buckingham");
      p_int("tjA,ntiA");
      p_real("a,b,vnbexp,vnb6,vnbtot");
    }
    if (bLJC)
      p_real("c12");
    p_int("n0,n1,nnn");
  }
  else {
    if (bLJC) {
      comment("More Lennard Jones");
      p_creal("twelve",12.0);
      p_creal("six",6.0);
      p_real("rinv6");
    }
    if (bBHAM) {
      comment("Buckingham");
      p_creal("six",6.0);
      p_int("tjA,ntiA");
      p_real("vnbexp,vnb6,vnbtot,r1,br,rinv6");
    }
  }
  if (bFREE) {
    comment("Free energy stuff");
    p_int("tjB,ntiB");
    p_real("L1,dvdl,qB,c6a,c6b,qqA,qqB,qqq");
    p_real("lam2,lam3,lam4,L12,L13,L14");
    if (bLJC)
      p_real("c12a,c12b");
  }
#ifdef CINVSQRT
  if (bC && bInline) {
    comment("Fast invsqrt stuff");
    p_creal("half",0.5);
    p_creal("three",3.0);
    fprintf(fp,"%st_convert  result,bit_pattern;\n",indent());
    fprintf(fp,"%sword       exp_addr,fract;\n",indent());
#ifdef DECREASE_LATENCY
    if (bWater) {
      p_real("rsqluO,rsqluH1,rsqluH2,tmpO,tmpH1,tmpH2");
      p_real4("luO,luH1,luH2");
    }
    else
      p_real4("lu");
#else
    p_real4("lu");
#endif /* DECREASE_LATENCY */
#ifdef DOUBLE
    p_real("y1,y2");
#endif /* DOUBLE */
  }
#endif /* CINVSQRT */

  comment("General and coulomb stuff");
  p_int("ii,k,n,jnr,ii3,nj0,nj1,is3,j3,ggid");
#ifdef USEVECTOR
  p_int("k3");
#endif  
  p_real("fxJ,fyJ,fzJ,vctot,fxO,fyO,fzO");
  p_real("ixO,iyO,izO,dxO,dyO,dzO");
  p_real("txO,tyO,tzO,vcO,fsO,qO,rsqO,rinv1O,rinv2O");
  if (!bFREE)
    p_real("qqO,qj");
  p_real("jx,jy,jz,shX,shY,shZ");
#ifdef FINVSQRT
  if (!bC) {
    comment("Fast invsqrt function");
    if (bInline) {
      fseed();
      p_int("fractaddr,expaddr");
#ifdef DECREASE_LATENCY
      if (bWater) {
	p_real("xin,luO,luH1,luH2,half,three");
      }	
      else {
	p_real("xin,lu,half,three");
      }
#else
      p_real("xin,lu,half,three");
#endif /* DECREASE_LATENCY */
#ifdef DOUBLE
      p_real("y1,y2");
#endif /* DOUBLE */
#ifdef DECREASE_LATENCY
      if (bWater) {
	p_int("iexp,addr,bval,resultrsqO,resultrsqH1,resultrsqH2");
	if (bEquiv) {
	  fprintf(fp,"%sequivalence(xin,bval)\n",indent());
	  fprintf(fp,"%sequivalence(luO,resultrsqO)\n",indent());
	  fprintf(fp,"%sequivalence(luH1,resultrsqH1)\n",indent());
	  fprintf(fp,"%sequivalence(luH2,resultrsqH2)\n",indent());
	}
      }
      else {
	p_int("iexp,addr,bval,result");
	if (bEquiv) {
	  fprintf(fp,"%sequivalence(xin,bval)\n",indent());
	  fprintf(fp,"%sequivalence(lu,result)\n",indent());
	}
      }
#else
      p_int("iexp,addr,bval,result");
      if (bEquiv) {
	fprintf(fp,"%sequivalence(xin,bval)\n",indent());
	fprintf(fp,"%sequivalence(lu,result)\n",indent());
      }
#endif /* DECREASE_LATENCY */
    }
    else
      p_real("invsqrt");
  }
#endif /* FINVSQRT */
#ifdef DECREASE_LATENCY
  if (bWater) {
    comment("Water stuff");
    p_real("fxH1,fyH1,fzH1,fxH2,fyH2,fzH2");
    p_real("ixH1,iyH1,izH1,ixH2,iyH2,izH2");
    p_real("dxH1,dyH1,dzH1,dxH2,dyH2,dzH2");
    p_real("txH1,tyH1,tzH1,txH2,tyH2,tzH2");
    p_real("rsqH1,rsqH2,rinv1H1,rinv1H2");
    if (!bC)
      p_real("rsqluO,rsqluH1,rsqluH2,tmpO,tmpH1,tmpH2");
    if (!bTab)
      p_real("rinv2H1,rinv2H2");
    p_real("vcH1,vcH2,fsH1,fsH2,qH,qqH");
  }
#else
  if (bWater) {
    comment("Water stuff");
#ifndef SIMPLEWATER
    p_real("fxH1,fyH1,fzH1,fxH2,fyH2,fzH2");
    p_real("ixH1,iyH1,izH1,ixH2,iyH2,izH2");
    p_real("dxH1,dyH1,dzH1,dxH2,dyH2,dzH2");
    p_real("txH1,tyH1,tzH1,txH2,tyH2,tzH2");
    p_real("rsqH1,rsqH2,rinv1H1,rinv1H2");
    if (!bTab)
      p_real("rinv2H1,rinv2H2");
    p_real("vcH1,vcH2,fsH1,fsH2,qH,qqH");
#else
    p_real("qH");
#endif
  }
#endif /* DECREASE_LATENCY */
  if(!bC && bEwald && !bTab)
    p_real("cerfc");
}

void flocal_init(void)
{
  fprintf(fp,"\n");
  if (!bC) {
    p_state("nul","0.0");
    if (bRF) 
      p_state("two","2.0");
    if(bEwald)
	p_state("isp","0.564189583547756");
    if (bTab) {
      p_state("one","1.0");
      p_state("two","2.0");
    }
    else {
      if (bLJC) {
	p_state("twelve","12.0");
	p_state("six","6.0");
      }
      if (bBHAM) 
	p_state("six","6.0");
    }
#ifdef FINVSQRT
    if (bInline) {
      p_state("half","0.5");
      p_state("three","3.0");
    }
#endif
  }
  if (bFREE) {
    p_state("L1","one - lambda");
    p_state("dvdl","nul");
    p_state("lam2","lambda*lambda");
    p_state("L12","L1*L1");
    p_state("lam3","lambda*lam2");
    p_state("L13","L1*L1*L1");
    p_state("lam4","lam2*lam2");
    p_state("L14","L12*L12");
  }
  
  fprintf(fp,"\n");
}

int fextract(char *s)
{
  p_state("nnn",s);
  p_state("Y",ARR(VFtab,nnn));
  p_state("F",ARR(VFtab,nnn+1)); 
  sprintf(buf,"eps*%s",ARR(VFtab,nnn+2));
  p_state("Geps",buf);
  sprintf(buf,"eps2*%s",ARR(VFtab,nnn+3));
  p_state("Heps2",buf);
  p_state("Fp","F+Geps+Heps2"); 
  p_state("VV","Y+eps*Fp"); 
  p_state("FF","Fp+Geps+two*Heps2");
  
  return 9;
}

int p_sqr(char *atom)
{
  char rsq[12],lm[12],rm[32],lm2[32],rrsq[64];
  int  m;
  
  for(m='x'; (m<='z'); m++) {  
    sprintf(lm,"d%c%s",m,atom);
    sprintf(lm2,"%s*%s",lm,lm);
    sprintf(rm,"i%c%s - j%c",m,atom,m);
    p_state(lm,rm);
    if (m == 'x')
      strcpy(rrsq,lm2);
    else {
      strcat(rrsq," + ");
      strcat(rrsq,lm2);
    }
  }
  sprintf(rsq,"rsq%s",atom);
  p_state(rsq,rrsq);
  
  return 8;
}

int fetch_j(void)
{
  int nflop = 0;
  
  comment("Unpack neighbourlist");
  p_state("jnr",IARR(jjnr,k));
  if (bC)
    p_state("j3","3*jnr");
  else
    p_state("j3","3*jnr-2");
#ifdef USEVECTOR
  if (bC)
    p_state("k3","3*(k-nj0)");
  else
    p_state("k3","3*(k-nj0)+1");
#endif

  return nflop;
}

void start_loop(char *lvar,char *from,char *to)
{
  if (bC) {
    fprintf(fp,"%sfor(%s=%s; (%s<%s); %s++) {\n",
	    indent(),lvar,from,lvar,to,lvar);
    CIND += 2;
  }
  else {
    fprintf(fp,"%sdo %s=%s,%s\n",indent(),lvar,from,to);
    FIND += 2;
  }
}

void end_loop(void)
{
  if (bC) {
    CIND -= 2;
    fprintf(fp,"%s}\n",indent());
  }
  else {
    FIND -= 2;
    fprintf(fp,"%send do\n",indent());
  }
}

int unpack_f(void)
{
  int nflop = 0;

  comment("Unpack forces for vectorization");
  p_ivdep();
  start_loop("k","nj0","nj1");
  nflop += fetch_j();
  
  p_state(ARR(fbuf,k3),ARR(faction,j3));
  p_state(ARR(fbuf,k3+1),ARR(faction,j3+1));
  p_state(ARR(fbuf,k3+2),ARR(faction,j3+2));
  
  end_loop();
  
  return nflop;
}

void finnerloop(char *loopname)
{
  int nflop=0;

#ifdef USEVECTOR  
  nflop += unpack_f();
#endif
  comment("Inner loop (over j-particles) starts right here");
  p_ivdep();
  start_loop("k","nj0","nj1");
  nflop += fetch_j();
  if (!bFREE) 
    p_state("qj",ARR(charge,jnr));
  
  p_state("jx",ARR(pos,j3));
  p_state("jy",ARR(pos,j3+1));
  p_state("jz",ARR(pos,j3+2));  
  
#ifdef DECREASE_LATENCY
  bDelayInvsqrt = FALSE;
  if (bWater) {
    comment("Ramones forever");
    comment("Determine square for Oxygen");
    nflop += p_sqr("O");
    comment("Determine square for first Hydrogen");
    nflop += p_sqr("H1");
    comment("Determine square for second Hydrogen");
    nflop += p_sqr("H2");
    nflop += p_waterinvsqrt();
  }
  else {
    comment("First one is for oxygen, with LJ");
    nflop += p_sqr("O");
    nflop += p_invsqrt("rinv1O","rsqO");
  }
#else
  comment("First one is for oxygen, with LJ");
  nflop += p_sqr("O");
  nflop += p_invsqrt("rinv1O","rsqO");

  if (!bDelayInvsqrt) {  
    if (bWater) {
      comment("Now for first hydrogen");
      nflop += p_sqr("H1");
      nflop += p_invsqrt("rinv1H1","rsqH1");
      
      comment("Now for second hydrogen");
      nflop += p_sqr("H2");
      nflop += p_invsqrt("rinv1H2","rsqH2");
    }
  }
#endif /* DECREASE_LATENCY */ 

#ifdef PREFETCH
  comment("Prefetch the forces");
  p_state("fxJ",ARR(faction,j3));
  p_state("fyJ",ARR(faction,j3+1));
  p_state("fzJ",ARR(faction,j3+2));
  nflop += 3;
#endif

#ifndef DECREASE_LATENCY
  comment("O block");
#endif
  if (!bFREE) {
    p_state("qqO","qO*qj");
    nflop ++;
  }
  else {
    sprintf(buf,"qO*%s",ARR(charge,jnr));
    p_state("qqA",buf);
    sprintf(buf,"qB*%s",ARR(chargeB,jnr));
    p_state("qqB",buf);
    p_state("qqq","(L12*qqA + lam2*qqB)");
    nflop += 5;
  }

#ifdef DECREASE_LATENCY
  if (bWater) {
    comment("inverse squares (for forces)");
    p_state("rinv2O","rinv1O*rinv1O");
    if (!bTab) {
      p_state("rinv2H1","rinv1H1*rinv1H1");
      p_state("rinv2H2","rinv1H2*rinv1H2");
      nflop += 2;
    }
    nflop++;
  }
  else {
    p_state("rinv2O","rinv1O*rinv1O");
    nflop ++;
  }
#else
  p_state("rinv2O","rinv1O*rinv1O");
  nflop ++;
#endif /* DECREASE_LATENCY */

  if (bTab) {
    p_state("r1","one/rinv1O");
    p_state("r1t","r1*tabscale");
    p_state("n0","r1t");
    if (bC)
      p_state("n1","12*n0");
    else 
      p_state("n1","12*n0+1");
    p_state("eps","r1t-n0");
    p_state("eps2","eps*eps");

    nflop += 8;
    
    nflop += fextract("n1");
  
    if (bFREE) {
      p_state("vcO","qqq*VV");
      p_state("fijC","qqq*FF");
      p_incr("vctot","vcO");
      p_incr("dvdl","two*(lambda*qqB - L1*qqA)*VV");
      nflop += 9;    
    }
    else {
      p_state("vcO","qqO*VV");
      p_state("fijC","qqO*FF");
      p_incr("vctot","vcO");
      nflop += 3;
    }
    
    if (bLJC || bBHAM) {
      comment("Dispersion table");
      sprintf(buf,"ntiA+%d*%s",bLJC ? 2 : 3,IARR(type,jnr));
      p_state("tjA",buf);
      if (bFREE) {
	sprintf(buf,"ntiB+%d*%s",bLJC ? 2 : 3,IARR(typeB,jnr));
	p_state("tjB",buf);
	p_state("c6a",ARR(nbfp,tjA));
	p_state("c6b",ARR(nbfp,tjB));
	p_state("c6","L14*c6a + L13*c6b");
	nflop += 3;
      }
      else 
	p_state("c6",ARR(nbfp,tjA));
    
      nflop += fextract("n1+4");
      p_state("vnb6","c6*VV");
      p_state("fijD","c6*FF");
      nflop += 2;
      if (bFREE) {
	p_incr("dvdl","3.0*(lam2*c6b  - L12*c6a)*VV");
	nflop += 6;
      }
            
      comment("Repulsion table");
      if (bLJC) {
	comment("Lennard Jones");
	sprintf(buf,"%s",ARR(nbfp,tjA+1));
	if (bFREE) {
	  p_state("c12a",buf);
	  p_state("c12b",ARR(nbfp,tjB+1));
	  p_state("c12","L14*c12a + lam4*c12b");
	  nflop += 3;
	}
	else 
	  p_state("c12",buf);
	
	nflop += fextract("n1+8");
	
	p_state("vnb12","c12*VV");
	p_state("fijR","c12*FF");
	p_state("vnbtot","vnbtot+vnb12+vnb6");
      
	comment("Total force");
	p_state("fsO","-(fijD+fijR+fijC)*tabscale*rinv1O");
	nflop += 9;
	
	if (bFREE) {
	  p_incr("dvdl","(4.0*(lam3*c12b - L13*c12a)*VV)");
	  nflop += 6;
	}
      }
      else {
	comment("Buckingham");
	p_state("b",ARR(nbfp,tjA+1));
	p_state("r1t","b*r1*exptabscale");
	p_state("n0","r1t");
	if (bC)
	  p_state("n1","12*n0");
	else
	  p_state("n1","12*n0+1");
	p_state("eps","r1t-n0");
	p_state("eps2","eps*eps");
	p_state("a",ARR(nbfp,tjA));
	nflop += 4;
	
	nflop += fextract("n1+8");
	p_state("vnbexp","a*VV");
	p_state("fijR","a*b*FF");
	p_state("vnbtot","vnbtot + vnbexp + vnb6");
	comment("Total force");
	p_state("fsO","-((fijD + fijC)*tabscale + fijR*exptabscale)*rinv1O");
	
	nflop += 10;
      }
    }
    else {
      comment("Coulomb on O only");
      p_state("fsO","-fijC*tabscale*rinv1O");
      nflop += 3;
    }
  }
  else {
    comment("No tables...");
    if (bLJC) {
      comment("Lennard Jones");
      p_state("rinv6","rinv2O*rinv2O*rinv2O");
      sprintf(buf,"ntiA+2*%s",IARR(type,jnr));
      p_state("tjA",buf);
      sprintf(buf,"rinv6*%s",ARR(nbfp,tjA));
      p_state("vnb6",buf);
      sprintf(buf,"rinv6*rinv6*%s",ARR(nbfp,tjA+1));
      p_state("vnb12",buf);
      p_state("vnbtot","vnbtot+vnb12-vnb6");
      nflop += 7;
    }
    else if (bBHAM) {
      comment("Buckingham");
      p_state("r1","rsqO*rinv1O");
      p_state("rinv6","rinv2O*rinv2O*rinv2O");
      sprintf(buf,"ntiA+3*%s",IARR(type,jnr));
      p_state("tjA",buf);
      sprintf(buf,"r1*%s",ARR(nbfp,tjA+1));
      p_state("br",buf);
      sprintf(buf,"exp(-br)*%s",ARR(nbfp,tjA));
      p_state("vnbexp",buf);
      sprintf(buf,"rinv6*%s",ARR(nbfp,tjA+2));
      p_state("vnb6",buf);
      p_state("vnbtot","vnbtot+vnbexp-vnb6");
      nflop += 13;
    }
    if (bRF) {
      comment("Reaction field stuff");
      p_state("krsqO","krf*rsqO");
      p_state("vcO","qqO*(rinv1O+krsqO)");
      nflop += 3;
      if (bLJC) {
	p_state("fsO",
		"(twelve*vnb12-six*vnb6+qqO*(rinv1O-two*krsqO))*rinv2O");
	nflop += 8;
      }
      else if (bBHAM) {
	p_state("fsO","(br*vnbexp-six*vnb6+qqO*(rinv1O-two*krsqO))*rinv2O");
	nflop += 8;
      }
      else {
	p_state("fsO","qqO*(rinv1O-two*krsqO)*rinv2O");
	nflop += 4;
      }
    }
    else if (bEwald) {
      comment("Ewald corrected terms");
      if(bC)
	p_state("vcO","qqO*rinv1O*erfc(ewc/rinv1O)");
      else
	p_state("vcO","qqO*rinv1O*cerfc(ewc/rinv1O)");	    
      nflop += 20; /* wild guess for erfc? */
      if (bLJC) {
	p_state("ewtmp","(vcO+2*qqO*ewc*exp(-ewc*ewc*rsqO)*isp)");
	p_state("fsO","(twelve*vnb12-six*vnb6+ewtmp)*rinv2O");
	nflop += 5;
      } 
      else if (bBHAM) {
	p_state("ewtmp","(vcO+2*qqO*ewc*exp(-ewc*ewc*rsqO)*isp)");
	p_state("fsO","(br*vnbexp-six*vnb6+ewtmp)*rinv2O"); 
	nflop += 5;
      } 
      else {
	p_state("fsO","rinv2O*(vcO+2*qqO*ewc*exp(-ewc*ewc*rsqO)*isp)");
	nflop += 1;
      } 
    } 
    else {
      comment("No reaction field or Ewald correction");
      p_state("vcO","qqO*rinv1O");
      nflop += 1;
      if (bLJC) {
	p_state("fsO","(twelve*vnb12-six*vnb6+vcO)*rinv2O");
	nflop += 5;
      }
      else if (bBHAM) {
	p_state("fsO","(br*vnbexp-six*vnb6+vcO)*rinv2O");
	nflop += 5;
      }
      else {
	p_state("fsO","vcO*rinv2O");
	nflop += 1;
      }
    }
    p_incr("vctot","vcO");
    nflop += 1;
  }

#ifdef DECREASE_LATENCY
  if (bWater) {
    p_state("qqH","qH*qj");
    nflop += 1;
  }
#endif /* DECREASE_LATENCY */
  
  comment("Sum the forces");   
  p_state("txO","dxO * fsO");
  p_state("tyO","dyO * fsO");
  p_state("tzO","dzO * fsO");
  p_incr("fxO","txO");
  p_incr("fyO","tyO");
  p_incr("fzO","tzO");
#ifndef PREFETCH
  if (!bWater) {
#ifdef USEVECTOR
    p_sub(ARR(faction,j3),ARR(fbuf,k3),"txO");
    p_sub(ARR(faction,j3+1),ARR(fbuf,k3+1),"tyO");
    p_sub(ARR(faction,j3+2),ARR(fbuf,k3+2),"tzO");
#else
    p_sub(ARR(faction,j3),ARR(faction,j3),"txO");
    p_sub(ARR(faction,j3+1),ARR(faction,j3+1),"tyO");
    p_sub(ARR(faction,j3+2),ARR(faction,j3+2),"tzO");
#endif
  }
  else {
#ifdef USEVECTOR
    p_sub("fxJ",ARR(fbuf,k3),"txO");
    p_sub("fyJ",ARR(fbuf,k3+1),"tyO");
    p_sub("fzJ",ARR(fbuf,k3+2),"tzO");
#else
    p_sub("fxJ",ARR(faction,j3),"txO");
    p_sub("fyJ",ARR(faction,j3+1),"tyO");
    p_sub("fzJ",ARR(faction,j3+2),"tzO");
#endif
  }
  nflop += 9;
#else

  if (!bWater) {
    p_sub(ARR(faction,j3),"fxJ","txO");
    p_sub(ARR(faction,j3+1),"fyJ","tyO");
    p_sub(ARR(faction,j3+2),"fzJ","tzO");
  }
  else {
    p_state("fxJ","fxJ - txO");
    p_state("fyJ","fyJ - tyO");
    p_state("fzJ","fzJ - tzO");
  }
  nflop += 6;
#endif
  
  if (bWater) {
#ifndef DECREASE_LATENCY
    p_state("qqH","qH*qj");
    nflop += 1;
#endif
    comment("H1 block");
    if (bDelayInvsqrt) {
      nflop += p_sqr("H1");
      nflop += p_invsqrt("rinv1H1","rsqH1");
    }
    if (bTab) {
      p_state("r1","one/rinv1H1");
      p_state("r1t","r1*tabscale");
      p_state("n0","r1t");
      if (bC)
	p_state("n1","12*n0");
      else
	p_state("n1","12*n0+1");
      p_state("eps","r1t-n0");
      p_state("eps2","eps*eps");
      nflop += 8;
	
      comment("Coulomb");
      nflop += fextract("n1");
      p_state("vcH1","qqH*VV");
      p_state("fijC","qqH*FF");
      p_state("fsH1","-fijC*tabscale*rinv1H1");
      nflop += 5;
    }
    else {
#ifndef DECREASE_LATENCY
      p_state("rinv2H1","rinv1H1*rinv1H1");
      nflop += 1;
#endif
      if (bRF) {
	comment("Reaction field");
	p_state("krsqH1","krf*rsqH1");
	p_state("vcH1","qqH*(rinv1H1+krsqH1)");
	p_state("fsH1","qqH*(rinv1H1-two*krsqH1)*rinv2H1");
	nflop += 7;
      }
      else if (bEwald) {
	comment("Ewald corrected terms");
	if(bC)
	  p_state("vcH1","qqH*rinv1H1*erfc(ewc/rinv1H1)");
	else
	  p_state("vcH1","qqH*rinv1H1*cerfc(ewc/rinv1H1)");	      
	p_state("fsH1","vcH1*rinv2H1+2*ewc*exp(-ewc*ewc*rsqH1)*rinv2H1*isp");        
      } else {
	comment("No reaction field or Ewald correction");
	p_state("vcH1","qqH*rinv1H1");
	p_state("fsH1","vcH1*rinv2H1");
	nflop += 2;
      }
    }
    p_incr("vctot","vcH1");
    
    comment("Forces on H1");
    p_state("txH1","dxH1 * fsH1");
    p_state("fxH1","fxH1 + txH1");
    p_decr("fxJ", "txH1");
    p_state("tyH1","dyH1 * fsH1");
    p_state("fyH1","fyH1 + tyH1");
    p_decr("fyJ", "tyH1");
    p_state("tzH1","dzH1 * fsH1");
    p_state("fzH1","fzH1 + tzH1");
    p_decr("fzJ", "tzH1");
    nflop += 10;
    
    comment("H2 block");
    if (bDelayInvsqrt) {
      nflop += p_sqr("H2");
      nflop += p_invsqrt("rinv1H2","rsqH2");
    }
    if (bTab) {
      p_state("r1","one/rinv1H2");
      p_state("r1t","r1*tabscale");
      p_state("n0","r1t");
      if (bC)
	p_state("n1","12*n0");
      else
	p_state("n1","12*n0+1");
      p_state("eps","r1t-n0");
      p_state("eps2","eps*eps");
      nflop += 8;
      
      comment("Coulomb");
      nflop += fextract("n1");
      p_state("vcH2","qqH*VV");
      p_state("fijC","qqH*FF");
      p_state("fsH2","-fijC*tabscale*rinv1H2");
      nflop += 5;
    }
    else {
#ifndef DECREASE_LATENCY
      p_state("rinv2H2","rinv1H2*rinv1H2");
      nflop += 1;
#endif
      if (bRF) {
	comment("Reaction field");
	p_state("krsqH2","krf*rsqH2");
	p_state("vcH2","qqH*(rinv1H2+krsqH2)");
	p_state("fsH2","qqH*(rinv1H2-two*krsqH2)*rinv2H2");
	nflop += 7;
      }
      else if (bEwald) {
	comment("Ewald corrected terms");
	if(bC)
	  p_state("vcH2","qqH*rinv1H2*erfc(ewc/rinv1H2)");
	else
	  p_state("vcH2","qqH*rinv1H2*cerfc(ewc/rinv1H2)");	      
	p_state("fsH2","vcH2*rinv2H2+2*ewc*exp(-ewc*ewc*rsqH2)*rinv2H2*isp");        
      }
      else {
	comment("No reaction field");
	p_state("vcH2","qqH*rinv1H2");
	p_state("fsH2","vcH2*rinv2H2");
	nflop += 2;
      }
    }
    p_incr("vctot","vcH2");
    
    comment("Forces on H2");
    p_state("txH2","dxH2*fsH2");
    p_incr("fxH2","txH2");
    p_state(ARR(faction,j3),"fxJ - txH2");
    p_state("tyH2","dyH2*fsH2");
    p_incr("fyH2","tyH2");
    p_state(ARR(faction,j3+1),"fyJ - tyH2");
    p_state("tzH2","dzH2*fsH2");
    p_incr("fzH2","tzH2");
    p_state(ARR(faction,j3+2),"fzJ - tzH2");
    nflop += 10;
  }

  comment("One more ultra fast GROMACS innerloop");
  sprintf(buf,"Innerloop %s costs %d flops",loopname,nflop);
  comment(buf);
  end_loop();
}

int fouter1(void)
{
  int nflop = 0;
  
  comment("Outer loop (over i particles) starts here");
#ifdef USE_OMP
  if (!bC) {
    fprintf(fp,"!$OMP DO\n");
  }
#endif
  start_loop("n",bC ? "0" : "1","nri");

  comment("Unpack shift vector");
  sprintf(buf,"3*%s",IARR(shift,n));     p_state("is3",buf);
  
  sprintf(buf,"%s",ARR(shiftvec,is3));   p_state("shX",buf);
  sprintf(buf,"%s",ARR(shiftvec,is3+1)); p_state("shY",buf);
  sprintf(buf,"%s",ARR(shiftvec,is3+2)); p_state("shZ",buf);
  
  comment("Unpack I particle");
  sprintf(buf,"%s",IARR(iinr,n));
  p_state("ii",buf);
  if (bC)
    p_state("ii3","3*ii");
  else
    p_state("ii3","3*ii-2");
    
  comment("Local variables for energy");  
  p_state("vctot","nul");
  if (bLJC || bBHAM) {
    p_state("vnbtot","nul");
    comment("Atom types for Van der Waals interactions");
    sprintf(buf,"%d*ntype*%s",bLJC ? 2 : 3,ARR(type,ii));
    p_state("ntiA",buf);
    if (bFREE) {
      sprintf(buf,"%d*ntype*%s",bLJC ? 2 : 3,ARR(typeB,ii));
      p_state("ntiB",buf);
    }
  }
  comment("Charge of i particle(s) divided by 4 pi eps0");
  sprintf(buf,"facel*%s",ARR(charge,ii));
  p_state("qO",buf);
  if (bFREE) {
    sprintf(buf,"facel*%s",ARR(chargeB,ii));
    p_state("qB",buf);
  }
  nflop++;
  
  if (bWater) {
    sprintf(buf,"facel*%s",ARR(charge,ii+1));
    p_state("qH",buf);
    nflop++;
  }
  
  comment("Bounds for the innerloop");
  p_state("nj0",IARR(jindex,n));
  p_state("nj1",ARR(jindex,n+1));
  
  comment("Compute shifted i position");
  p_add("ixO","shX",ARR(pos,ii3));
  p_add("iyO","shY",ARR(pos,ii3+1));
  p_add("izO","shZ",ARR(pos,ii3+2));
  
  p_state("fxO","nul");
  p_state("fyO","nul");
  p_state("fzO","nul");
  
  nflop += 3;
  
  return nflop;
}

void fouterloop(char *loopname)
{
  int nflop;
  
  nflop = fouter1();
  
  if (bWater) {
    p_add("ixH1","shX",ARR(pos,ii3+3));
    p_add("iyH1","shY",ARR(pos,ii3+4));
    p_add("izH1","shZ",ARR(pos,ii3+5));
    p_add("ixH2","shX",ARR(pos,ii3+6));
    p_add("iyH2","shY",ARR(pos,ii3+7));
    p_add("izH2","shZ",ARR(pos,ii3+8));
    nflop += 6;
    
    p_state("fxH1","nul");
    p_state("fyH1","nul");
    p_state("fzH1","nul");
    p_state("fxH2","nul");
    p_state("fyH2","nul");
    p_state("fzH2","nul");
  }
  
  finnerloop(loopname);

  sprintf(buf,"%s",ARR(faction,ii3));
  p_incr(buf,"fxO");
  sprintf(buf,"%s",ARR(faction,ii3+1));
  p_incr(buf,"fyO");
  sprintf(buf,"%s",ARR(faction,ii3+2));
  p_incr(buf,"fzO");
  nflop += 3;
  if (bWater) {
    p_incr(ARR(fshift,is3),  "fxO+fxH1+fxH2");
    p_incr(ARR(fshift,is3+1),"fyO+fyH1+fyH2");
    p_incr(ARR(fshift,is3+2),"fzO+fzH1+fzH2");
    p_incr(ARR(faction,ii3+3),"fxH1");
    p_incr(ARR(faction,ii3+4),"fyH1");
    p_incr(ARR(faction,ii3+5),"fzH1");
    p_incr(ARR(faction,ii3+6),"fxH2");
    p_incr(ARR(faction,ii3+7),"fyH2");
    p_incr(ARR(faction,ii3+8),"fzH2");
    nflop += 15;
  }
  else {
    p_incr(ARR(fshift,is3),"fxO");
    p_incr(ARR(fshift,is3+1),"fyO");
    p_incr(ARR(fshift,is3+2),"fzO");
    nflop += 3;
  }  
  
  comment("Update energies");
  sprintf(buf,"%s",IARR(gid,n));
  p_state("ggid",buf);
  sprintf(buf,"%s",ARR(Vc,ggid));
  p_incr(buf,"vctot");
  nflop++;
  
  if (bLJC || bBHAM) {
    sprintf(buf,"%s",ARR(Vnb,ggid));
    p_incr(buf,"vnbtot");
    nflop++;
  }
  
  sprintf(buf,"Outer loop costs %d flops",nflop);
  comment(buf);
  end_loop();
#ifdef USE_OMP
  if (!bC)
    fprintf(fp,"!$OMP END DO\n");
#endif
  if (bFREE) {
    if (bC)
      p_incr("*dvdlambda","dvdl");
    else
      p_incr("dvdlambda","dvdl");
  }
}

void fouterloop2(char *loopname)
{
  int  nflop;
  bool bBSave,bCSave,bLJCSave;
  
  if (!bWater) {
    fprintf(stderr,"fouterloop2 called without water\n");
    exit(1);
  }

  comment("Silly expanded water loop because of broken SGI compiler");    
  nflop = fouter1();

  /* Oxygen loop */  
  bWater = FALSE;
  finnerloop(loopname);
  comment("Update forces on oxygen");
  p_incr(ARR(faction,ii3),"fxO");
  p_incr(ARR(faction,ii3+1),"fyO");
  p_incr(ARR(faction,ii3+2),"fzO");
  p_incr(ARR(fshift,is3),"fxO");
  p_incr(ARR(fshift,is3+1),"fyO");
  p_incr(ARR(fshift,is3+2),"fzO");
  nflop += 6;

  comment("Now starts first hydrogen");
  p_add("ixO","shX",ARR(pos,ii3+3));
  p_add("iyO","shY",ARR(pos,ii3+4));
  p_add("izO","shZ",ARR(pos,ii3+5));
  nflop += 3;
  
  comment("Copy the hydrogen charge to the oxygen");
  p_state("qO","qH");
  bBSave   = bBHAM;
  bLJCSave = bLJC;
  bCSave   = bCoul;
  
  /* Change options for the hydrogen */
  bCoul = TRUE;
  bBHAM = FALSE;
  bLJC  = FALSE;
  
  p_state("fxO","nul");
  p_state("fyO","nul");
  p_state("fzO","nul");
  finnerloop(loopname);
  comment("Update forces on hydrogen 1");
  p_incr(ARR(faction,ii3+3),"fxO");
  p_incr(ARR(faction,ii3+4),"fyO");
  p_incr(ARR(faction,ii3+5),"fzO");
  p_incr(ARR(fshift,is3),"fxO");
  p_incr(ARR(fshift,is3+1),"fyO");
  p_incr(ARR(fshift,is3+2),"fzO");
  nflop += 6;
  
  comment("Now starts second hydrogen");
  p_add("ixO","shX",ARR(pos,ii3+6));
  p_add("iyO","shY",ARR(pos,ii3+7));
  p_add("izO","shZ",ARR(pos,ii3+8));
  nflop += 3;
  
  p_state("fxO","nul");
  p_state("fyO","nul");
  p_state("fzO","nul");
  finnerloop(loopname);
  comment("Update forces on hydrogen 2");
  p_incr(ARR(faction,ii3+6),"fxO");
  p_incr(ARR(faction,ii3+7),"fyO");
  p_incr(ARR(faction,ii3+8),"fzO");
  p_incr(ARR(fshift,is3),"fxO");
  p_incr(ARR(fshift,is3+1),"fyO");
  p_incr(ARR(fshift,is3+2),"fzO");
  nflop += 6;
  
  comment("Update energies");
  sprintf(buf,"%s",IARR(gid,n));
  p_state("ggid",buf);
  sprintf(buf,"%s",ARR(Vc,ggid));
  p_incr(buf,"vctot");
  nflop++;
  
  bBHAM = bBSave;
  bLJC  = bLJCSave;
  bCoul = bCSave;
  
  if (bLJC || bBHAM) {
    sprintf(buf,"%s",ARR(Vnb,ggid));
    p_incr(buf,"vnbtot");
    nflop++;
  }
  
  sprintf(buf,"Outer loop costs %d flops",nflop);
  comment(buf);
  end_loop();

}

void closeit(void)
{
  fwarning();
  if (bC)
    fprintf(fp,"}\n\n");
  else
    fprintf(fp,"%sreturn\n%send\n\n",indent(),indent());
}

void make_funcnm(char *buf)
{
  sprintf(buf,"%s",bC ? "c_" : "f77");
  if (bBHAM)
    strcat(buf,"bham");
  else if (bCoul)
    strcat(buf,"coul");
  else
    strcat(buf,"ljc");
  if (bTab)
    strcat(buf,"tab"); 
  else if (bRF)
    strcat(buf,"rf");
}

void doit(void)
{
  char loopname[256];
  
  bLJC = (!bCoul && !bBHAM);
  
  if (bFREE && (!bTab || bCoul || bRF || bWater || bEwald)) {
    fprintf(stderr,"Unsupported combination of options to doit\n");
    return;
  }
  if(bRF && bEwald) {
    fprintf(stderr,"Unsupported combination of options to doit\n");
    return;
  }
  make_funcnm(loopname);
  if (bWater)
    strcat(loopname,"water");
  if (bFREE)
    strcat(loopname,"free");
  if(bEwald)
    strcat(loopname,"ewald");
  fname(loopname);
  fargs();
  flocal();
  flocal_init();
  
#ifdef SIMPLEWATER
  if (bWater)
    fouterloop2(loopname);
  else
    fouterloop(loopname);
#else
  fouterloop(loopname);
#endif

  closeit();
}

void error(int argc,char *argv[])
{
  fprintf(stderr,"Usage: %s bC bInline\n",argv[0]);
  fprintf(stderr,"\tbC      = 0: F77 code,            1: C code\n");
  fprintf(stderr,"\tbInline = 0: Don't inline invsqrt 1: Do inline\n");
  exit(-1);
}

int count_lines(char *fn)
{
  FILE *fp;
  int nl=0;
  
  if ((fp = fopen(fn,"r")) == NULL) {
    perror(fn);
    exit(1);
  }

  while (fgets(buf,255,fp) != NULL)
    nl++;
  fclose(fp);
  
  return nl;
}

int main(int argc,char *argv[])
{
  bool bbb[12][5] = {
    { 1, 0, 0, 0, 0 },
    { 1, 0, 0, 0, 1 },
    { 0, 0, 0, 0, 0 },
    { 0, 0, 0, 0, 1 },
    { 0, 1, 0, 0, 0 },
    { 0, 1, 0, 0, 1 },
    { 1, 0, 1, 0, 0 },
    { 0, 0, 1, 0, 0 },
    { 0, 1, 1, 0, 0 },
    { 1, 0, 0, 1, 0 },
    { 0, 0, 0, 1, 0 },
    { 0, 1, 0, 1, 0 },
  };
  int i,j,nfunc,nlines;
  char fn[256];
  
  if (argc != 3) 
    error(argc,argv);
  
  bC      = atoi(argv[1]);
  bInline = FALSE;

#ifdef CINVSQRT
  if (bC) 
    bInline = atoi(argv[2]);
#endif
#ifdef FINVSQRT
  if (!bC)
    bInline = atoi(argv[2]);
#endif

#ifdef DOUBLE
  prec   = 8;
#else
  prec   = 4;
#endif
  sprintf(fn,"%s",bC ? "innerc.c" : "innerf.f");
  fprintf(stderr,">>> %s is the GROMACS innerloop code generator\n",argv[0]);
  fprintf(stderr,">>> It will generate %s precision %s code in file %s\n",
	  (prec == 8) ? "double" : "single",
	  bC ? "C" : "Fortran",fn);
  if (bInline)
    fprintf(stderr,">>> The software 1.0/sqrt routine will be inlined\n");
#ifdef SIMPLEWATER
  fprintf(stderr,">>> Will expand the water loops to have 3 inner loops\n");
#endif
#ifdef USEVECTOR
  fprintf(stderr,">>> Will generate vectorizable code (hopefully)\n");
#endif
  if ((fp = fopen(fn,"w")) == NULL) {
    perror(fn);
    exit(1);
  }
  fheader();

  bFREE  = TRUE;  
  bCoul  = FALSE;
  bRF    = FALSE;
  bEwald = FALSE; /* no free ewald yet */
  bTab   = TRUE;
  bWater = FALSE;
  bBHAM  = FALSE;
  doit();
  bBHAM  = TRUE;
  doit();
  bFREE  = FALSE;
  nfunc  = 2;
   
  for(i=0; (i<2); i++) {
    for(j=0; (j<12); j++) {
      bCoul  = bbb[j][0];
      bBHAM  = bbb[j][1];
      bRF    = bbb[j][2];
      bTab   = bbb[j][3];
      bEwald = bbb[j][4];
      bWater = i;
      doit();
      nfunc++;
    }
  }
  fclose(fp);
  nlines = count_lines(fn);
  fprintf(stderr,">>> Generated %d lines of code in %d functions\n",
	  nlines,nfunc);
    
  return 0;
}
