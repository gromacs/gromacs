#ifndef _inner_h
#define _inner_h
#include "typedefs.h"

define(`LJC_ARGS',`int SCALARG(ntype),int type[],real nbfp[],real Vnb[]')
define(`BHAM_ARGS',`LJC_ARGS')
define(`FEP_ARGS',`real chargeB[],int typeB[],real SCALARG(lambda),real *dvdlambda')
define(`VEC_ARGS',`real fbuf[]')
define(`RF_ARGS',`real SCALARG(krf)')
define(`EWALD_ARGS',`real SCALARG(ewaldcoeff)')
define(`TAB_ARGS',`real SCALARG(tabscale),real VFtab[]')
define(`BHTAB_ARGS',`real SCALARG(tabscale_exp)')
ifdef(`USEVECTOR',`define(`ALL_ARGS',`real fbuf[],int SCALARG(nri),int iinr[],int `shift'[],int gid[],int jindex[],int jjnr[],real pos[],real fshift[],real SCALARG(facel),real charge[],real faction[],real Vc[],real shiftvec[]')',`define(`ALL_ARGS',`int SCALARG(nri),int iinr[],int `shift'[],int gid[],int jindex[],int jjnr[],real pos[],real fshift[],real SCALARG(facel),real charge[],real faction[],real Vc[],real shiftvec[]')')
		      
extern void FUNC(ljctabfree)(LJC_ARGS,TAB_ARGS,FEP_ARGS,ALL_ARGS);

extern void FUNC(bhamtabfree)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,FEP_ARGS,ALL_ARGS);

extern void FUNC(ljcrf)(LJC_ARGS,RF_ARGS,ALL_ARGS);

extern void FUNC(bhamrf)(BHAM_ARGS,RF_ARGS,ALL_ARGS);

extern void FUNC(coulrf)(RF_ARGS,ALL_ARGS);

extern void FUNC(ljctab)(LJC_ARGS,TAB_ARGS,ALL_ARGS);

extern void FUNC(bhamtab)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,ALL_ARGS);

extern void FUNC(coultab)(TAB_ARGS,ALL_ARGS);

extern void FUNC(ljc)(LJC_ARGS,ALL_ARGS);

extern void FUNC(ljcewald)(LJC_ARGS,EWALD_ARGS,ALL_ARGS);

extern void FUNC(bham)(BHAM_ARGS,ALL_ARGS);

extern void FUNC(bhamewald)(BHAM_ARGS,EWALD_ARGS,ALL_ARGS);

extern void FUNC(coul)(ALL_ARGS);

extern void FUNC(coulewald)(EWALD_ARGS,ALL_ARGS);

extern void FUNC(ljcrfwater)(LJC_ARGS,RF_ARGS,ALL_ARGS);

extern void FUNC(bhamrfwater)(BHAM_ARGS,RF_ARGS,ALL_ARGS);

extern void FUNC(coulrfwater)(RF_ARGS,ALL_ARGS);

extern void FUNC(ljctabwater)(LJC_ARGS,TAB_ARGS,ALL_ARGS);

extern void FUNC(bhamtabwater)(BHAM_ARGS,TAB_ARGS,BHTAB_ARGS,ALL_ARGS);

extern void FUNC(coultabwater)(TAB_ARGS,ALL_ARGS);

extern void FUNC(ljcwater)(LJC_ARGS,ALL_ARGS);

extern void FUNC(ljcwaterewald)(LJC_ARGS,EWALD_ARGS,ALL_ARGS);

extern void FUNC(bhamwater)(BHAM_ARGS,ALL_ARGS);

extern void FUNC(bhamwaterewald)(BHAM_ARGS,EWALD_ARGS,ALL_ARGS);

extern void FUNC(coulwater)(ALL_ARGS);

extern void FUNC(coulwaterewald)(EWALD_ARGS,ALL_ARGS);

#endif
