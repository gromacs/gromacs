;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).
;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text



.globl nb_kernel301_ia32_sse
.globl _nb_kernel301_ia32_sse
nb_kernel301_ia32_sse:	
_nb_kernel301_ia32_sse:	
.equiv          nb301_p_nri,            8
.equiv          nb301_iinr,             12
.equiv          nb301_jindex,           16
.equiv          nb301_jjnr,             20
.equiv          nb301_shift,            24
.equiv          nb301_shiftvec,         28
.equiv          nb301_fshift,           32
.equiv          nb301_gid,              36
.equiv          nb301_pos,              40
.equiv          nb301_faction,          44
.equiv          nb301_charge,           48
.equiv          nb301_p_facel,          52
.equiv          nb301_argkrf,           56
.equiv          nb301_argcrf,           60
.equiv          nb301_Vc,               64
.equiv          nb301_type,             68
.equiv          nb301_p_ntype,          72
.equiv          nb301_vdwparam,         76
.equiv          nb301_Vvdw,             80
.equiv          nb301_p_tabscale,       84
.equiv          nb301_VFtab,            88
.equiv          nb301_invsqrta,         92
.equiv          nb301_dvda,             96
.equiv          nb301_p_gbtabscale,     100
.equiv          nb301_GBtab,            104
.equiv          nb301_p_nthreads,       108
.equiv          nb301_count,            112
.equiv          nb301_mtx,              116
.equiv          nb301_outeriter,        120
.equiv          nb301_inneriter,        124
.equiv          nb301_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb301_ixO,              0
.equiv          nb301_iyO,              16
.equiv          nb301_izO,              32
.equiv          nb301_ixH1,             48
.equiv          nb301_iyH1,             64
.equiv          nb301_izH1,             80
.equiv          nb301_ixH2,             96
.equiv          nb301_iyH2,             112
.equiv          nb301_izH2,             128
.equiv          nb301_iqO,              144
.equiv          nb301_iqH,              160
.equiv          nb301_dxO,              176
.equiv          nb301_dyO,              192
.equiv          nb301_dzO,              208
.equiv          nb301_dxH1,             224
.equiv          nb301_dyH1,             240
.equiv          nb301_dzH1,             256
.equiv          nb301_dxH2,             272
.equiv          nb301_dyH2,             288
.equiv          nb301_dzH2,             304
.equiv          nb301_qqO,              320
.equiv          nb301_qqH,              336
.equiv          nb301_rinvO,            352
.equiv          nb301_rinvH1,           368
.equiv          nb301_rinvH2,           384
.equiv          nb301_rO,               400
.equiv          nb301_rH1,              416
.equiv          nb301_rH2,              432
.equiv          nb301_tsc,              448
.equiv          nb301_two,              464
.equiv          nb301_vctot,            480
.equiv          nb301_fixO,             496
.equiv          nb301_fiyO,             512
.equiv          nb301_fizO,             528
.equiv          nb301_fixH1,            544
.equiv          nb301_fiyH1,            560
.equiv          nb301_fizH1,            576
.equiv          nb301_fixH2,            592
.equiv          nb301_fiyH2,            608
.equiv          nb301_fizH2,            624
.equiv          nb301_fjx,              640
.equiv          nb301_fjy,              656
.equiv          nb301_fjz,              672
.equiv          nb301_half,             688
.equiv          nb301_three,            704
.equiv          nb301_is3,              720
.equiv          nb301_ii3,              724
.equiv          nb301_innerjjnr,        728
.equiv          nb301_innerk,           732
.equiv          nb301_n,                736
.equiv          nb301_nn1,              740
.equiv          nb301_nri,              744
.equiv          nb301_nouter,           748
.equiv          nb301_ninner,           752
.equiv          nb301_salign,           756
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 760		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb301_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb301_p_nri]
	mov ecx, [ecx]
	mov [esp + nb301_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301_nouter], eax
	mov [esp + nb301_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb301_half], eax
	movss xmm1, [esp + nb301_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb301_half],  xmm1
	movaps [esp + nb301_two],  xmm2
	movaps [esp + nb301_three],  xmm3
	mov eax, [ebp + nb301_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0 
	movaps [esp + nb301_tsc], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb301_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb301_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb301_iqO], xmm3
	movaps [esp + nb301_iqH], xmm4
	
.nb301_threadloop:
        mov   esi, [ebp + nb301_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb301_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb301_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb301_n], eax
        mov [esp + nb301_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301_outerstart
        jmp .nb301_end
	
.nb301_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb301_nouter]
	mov [esp + nb301_nouter], ebx

.nb301_outer:
	mov   eax, [ebp + nb301_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb301_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb301_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb301_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301_pos]    ;# eax = base of pos[]  
	mov   [esp + nb301_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb301_ixO], xmm3
	movaps [esp + nb301_iyO], xmm4
	movaps [esp + nb301_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb301_ixH1], xmm0
	movaps [esp + nb301_iyH1], xmm1
	movaps [esp + nb301_izH1], xmm2
	movaps [esp + nb301_ixH2], xmm3
	movaps [esp + nb301_iyH2], xmm4
	movaps [esp + nb301_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb301_vctot], xmm4
	movaps [esp + nb301_fixO], xmm4
	movaps [esp + nb301_fiyO], xmm4
	movaps [esp + nb301_fizO], xmm4
	movaps [esp + nb301_fixH1], xmm4
	movaps [esp + nb301_fiyH1], xmm4
	movaps [esp + nb301_fizH1], xmm4
	movaps [esp + nb301_fixH2], xmm4
	movaps [esp + nb301_fiyH2], xmm4
	movaps [esp + nb301_fizH2], xmm4
	
	mov   eax, [ebp + nb301_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb301_pos]
	mov   edi, [ebp + nb301_faction]	
	mov   eax, [ebp + nb301_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb301_ninner]
	mov   [esp + nb301_ninner], ecx
	add   edx, 0
	mov   [esp + nb301_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301_unroll_loop
	jmp   .nb301_odd_inner
.nb301_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb301_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb301_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb301_iqO]
	mulps  xmm4, [esp + nb301_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb301_qqO], xmm3
	movaps  [esp + nb301_qqH], xmm4	

	mov esi, [ebp + nb301_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb301_ixO]
	movaps xmm5, [esp + nb301_iyO]
	movaps xmm6, [esp + nb301_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb301_dxO], xmm4
	movaps [esp + nb301_dyO], xmm5
	movaps [esp + nb301_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb301_ixH1]
	movaps xmm5, [esp + nb301_iyH1]
	movaps xmm6, [esp + nb301_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb301_dxH1], xmm4
	movaps [esp + nb301_dyH1], xmm5
	movaps [esp + nb301_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb301_ixH2]
	movaps xmm4, [esp + nb301_iyH2]
	movaps xmm5, [esp + nb301_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb301_dxH2], xmm3
	movaps [esp + nb301_dyH2], xmm4
	movaps [esp + nb301_dzH2], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301_half]
	movaps  [esp + nb301_rinvO], xmm4	;# rinvO in xmm4 
	mulps   xmm7, xmm4
	movaps  [esp + nb301_rO], xmm7	

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301_half]
	movaps  [esp + nb301_rinvH1], xmm4	;# rinvH1 in xmm4 
	mulps   xmm6, xmm4
	movaps  [esp + nb301_rH1], xmm6

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301_half]
	movaps  [esp + nb301_rinvH2], xmm4	;# rinvH2 in xmm4 
	mulps   xmm5, xmm4
	movaps  [esp + nb301_rH2], xmm5

	;# do O interactions 
	;# rO is still in xmm7 
	mulps   xmm7, [esp + nb301_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3

	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd mm0, eax   
    	movd mm1, ebx
    	movd mm2, ecx
    	movd mm3, edx

    	mov  esi, [ebp + nb301_VFtab]
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
    	    
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	mulps  xmm7, [esp + nb301_two]       ;# two*Heps2 
    	movaps xmm0, [esp + nb301_qqO]
    	addps  xmm7, xmm6
    	addps  xmm7, xmm5 ;# xmm7=FF 
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	mulps  xmm0, xmm7 ;# fijC=FF*qq 

    	;# at this point mm5 contains vcoul and xmm0 fijC 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [esp + nb301_vctot]
    	movaps [esp + nb301_vctot], xmm5 
	xorps  xmm4, xmm4

	mulps  xmm0, [esp + nb301_tsc]
	mulps  xmm0, [esp + nb301_rinvO]	
	subps  xmm4, xmm0

	movaps xmm0, [esp + nb301_dxO]
	movaps xmm1, [esp + nb301_dyO]
	movaps xmm2, [esp + nb301_dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4	;# tx in xmm0-xmm2 

	;# update O forces 
	movaps xmm3, [esp + nb301_fixO]
	movaps xmm4, [esp + nb301_fiyO]
	movaps xmm7, [esp + nb301_fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb301_fixO], xmm3
	movaps [esp + nb301_fiyO], xmm4
	movaps [esp + nb301_fizO], xmm7
	;# update j forces with water O 
	movaps [esp + nb301_fjx], xmm0
	movaps [esp + nb301_fjy], xmm1
	movaps [esp + nb301_fjz], xmm2

	;# Done with O interactions - now H1! 
	movaps xmm7, [esp + nb301_rH1]
	mulps   xmm7, [esp + nb301_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7
	
    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	mulps  xmm7, [esp + nb301_two]       ;# two*Heps2 
    	movaps xmm0, [esp + nb301_qqH]
    	addps  xmm7, xmm6
    	addps  xmm7, xmm5 ;# xmm7=FF 
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	mulps  xmm7, xmm0 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm7 fijC 
    	;# increment vcoul 
	xorps  xmm4, xmm4
    	addps  xmm5, [esp + nb301_vctot]
	mulps  xmm7, [esp + nb301_rinvH1]
    	movaps [esp + nb301_vctot], xmm5 
	mulps  xmm7, [esp + nb301_tsc]
	subps xmm4, xmm7

	movaps xmm0, [esp + nb301_dxH1]
	movaps xmm1, [esp + nb301_dyH1]
	movaps xmm2, [esp + nb301_dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb301_fixH1]
	movaps xmm4, [esp + nb301_fiyH1]
	movaps xmm7, [esp + nb301_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb301_fixH1], xmm3
	movaps [esp + nb301_fiyH1], xmm4
	movaps [esp + nb301_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [esp + nb301_fjx]
	addps  xmm1, [esp + nb301_fjy]
	addps  xmm2, [esp + nb301_fjz]
	movaps [esp + nb301_fjx], xmm0
	movaps [esp + nb301_fjy], xmm1
	movaps [esp + nb301_fjz], xmm2

	;# Done with H1, finally we do H2 interactions 
	movaps xmm7, [esp + nb301_rH2]
	mulps   xmm7, [esp + nb301_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 
	
    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	mulps  xmm7, [esp + nb301_two]       ;# two*Heps2 
    	movaps xmm0, [esp + nb301_qqH]
    	addps  xmm7, xmm6
    	addps  xmm7, xmm5 ;# xmm7=FF 
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	mulps  xmm7, xmm0 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm0 fijC 
    	;# increment vcoul 
	xorps  xmm4, xmm4
    	addps  xmm5, [esp + nb301_vctot]
	mulps  xmm7, [esp + nb301_rinvH2]
    	movaps [esp + nb301_vctot], xmm5 
	mulps  xmm7, [esp + nb301_tsc]
	subps  xmm4, xmm7

	movaps xmm0, [esp + nb301_dxH2]
	movaps xmm1, [esp + nb301_dyH2]
	movaps xmm2, [esp + nb301_dzH2]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

    	movd eax, mm0   
    	movd ebx, mm1
    	movd ecx, mm2
    	movd edx, mm3
	
	;# update H2 forces 
	movaps xmm3, [esp + nb301_fixH2]
	movaps xmm4, [esp + nb301_fiyH2]
	movaps xmm7, [esp + nb301_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb301_fixH2], xmm3
	movaps [esp + nb301_fiyH2], xmm4
	movaps [esp + nb301_fizH2], xmm7

	mov edi, [ebp + nb301_faction]
	;# update j forces 
	addps xmm0, [esp + nb301_fjx]
	addps xmm1, [esp + nb301_fjy]
	addps xmm2, [esp + nb301_fjz]

	movlps xmm4, [edi + eax*4]
	movlps xmm7, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm7, [edi + edx*4]
	
	movaps xmm3, xmm4
	shufps xmm3, xmm7, 136  ;# constant 10001000
	shufps xmm4, xmm7, 221  ;# constant 11011101			      
	;# xmm3 has fjx, xmm4 has fjy 
	subps xmm3, xmm0
	subps xmm4, xmm1
	;# unpack them back for storing 
	movaps xmm7, xmm3
	unpcklps xmm7, xmm4
	unpckhps xmm3, xmm4	
	movlps [edi + eax*4], xmm7
	movlps [edi + ecx*4], xmm3
	movhps [edi + ebx*4], xmm7
	movhps [edi + edx*4], xmm3
	;# finally z forces 
	movss  xmm0, [edi + eax*4 + 8]
	movss  xmm1, [edi + ebx*4 + 8]
	movss  xmm3, [edi + ecx*4 + 8]
	movss  xmm4, [edi + edx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm4, xmm2
	movss  [edi + eax*4 + 8], xmm0
	movss  [edi + ebx*4 + 8], xmm1
	movss  [edi + ecx*4 + 8], xmm3
	movss  [edi + edx*4 + 8], xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb301_innerk],  4
	jl    .nb301_odd_inner
	jmp   .nb301_unroll_loop
.nb301_odd_inner:	
	add dword ptr [esp + nb301_innerk],  4
	jnz   .nb301_odd_loop
	jmp   .nb301_updateouterdata
.nb301_odd_loop:
	mov   edx, [esp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb301_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb301_iqO]
	mov esi, [ebp + nb301_charge] 
	movhps xmm4, [esp + nb301_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb301_qqO], xmm3	;# use oxygen qq for storage 

	mov esi, [ebp + nb301_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb301_ixO]
	movss xmm4, [esp + nb301_iyO]
	movss xmm5, [esp + nb301_izO]
		
	movlps xmm6, [esp + nb301_ixH1]
	movlps xmm7, [esp + nb301_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb301_iyH1]
	movlps xmm7, [esp + nb301_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb301_izH1]
	movlps xmm7, [esp + nb301_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + nb301_dxO], xmm3
	movaps [esp + nb301_dyO], xmm4
	movaps [esp + nb301_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb301_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb301_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000	
	
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [esp + nb301_rinvO], xmm0
	
	mulps xmm4, [esp + nb301_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ;# mm6/mm7 contain lu indices 
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm7, mm7
    	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
    	movd mm0, eax   
    	movd mm1, ecx
    	movd mm2, edx

    	mov  esi, [ebp + nb301_VFtab]
    	movd eax, mm6
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	mulps  xmm7, [esp + nb301_two]       ;# two*Heps2 
    	movaps xmm0, [esp + nb301_qqO]
    	addps  xmm7, xmm6
    	addps  xmm7, xmm5 ;# xmm7=FF 
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	mulps  xmm0, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm0 fijC 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [esp + nb301_vctot]
    	movaps [esp + nb301_vctot], xmm5

	xorps xmm4, xmm4
	mulps  xmm0, [esp + nb301_tsc]
	mulps  xmm0, [esp + nb301_rinvO]	
	subps  xmm4, xmm0
		
    	movd eax, mm0   
    	movd ecx, mm1
    	movd edx, mm2	
		
	movaps xmm0, [esp + nb301_dxO]
	movaps xmm1, [esp + nb301_dyO]
	movaps xmm2, [esp + nb301_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force) 
	movss  xmm3, [esp + nb301_fixO]	
	movss  xmm4, [esp + nb301_fiyO]	
	movss  xmm5, [esp + nb301_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb301_fixO], xmm3	
	movss  [esp + nb301_fiyO], xmm4	
	movss  [esp + nb301_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb301_fixH1]
	addss  xmm4, [esp + nb301_fiyH1]
	addss  xmm5, [esp + nb301_fizH1]
	movss  [esp + nb301_fixH1], xmm3	
	movss  [esp + nb301_fiyH1], xmm4	
	movss  [esp + nb301_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb301_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb301_fixH2]
	addss  xmm4, [esp + nb301_fiyH2]
	addss  xmm5, [esp + nb301_fizH2]
	movss  [esp + nb301_fixH2], xmm3	
	movss  [esp + nb301_fiyH2], xmm4	
	movss  [esp + nb301_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [edi + eax*4]
	movss  xmm7, [edi + eax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1    ;# z sum in xmm2 
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4],     xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec dword ptr [esp + nb301_innerk]
	jz    .nb301_updateouterdata
	jmp   .nb301_odd_loop
.nb301_updateouterdata:
	mov   ecx, [esp + nb301_ii3]
	mov   edi, [ebp + nb301_faction]
	mov   esi, [ebp + nb301_fshift]
	mov   edx, [esp + nb301_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb301_fixO]
	movaps xmm1, [esp + nb301_fiyO]
	movaps xmm2, [esp + nb301_fizO]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb301_fixH1]
	movaps xmm1, [esp + nb301_fiyH1]
	movaps xmm2, [esp + nb301_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 12]
	movss  xmm4, [edi + ecx*4 + 16]
	movss  xmm5, [edi + ecx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 12], xmm3
	movss  [edi + ecx*4 + 16], xmm4
	movss  [edi + ecx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb301_fixH2]
	movaps xmm1, [esp + nb301_fiyH2]
	movaps xmm2, [esp + nb301_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4 + 24]
	movss  xmm4, [edi + ecx*4 + 28]
	movss  xmm5, [edi + ecx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 24], xmm3
	movss  [edi + ecx*4 + 28], xmm4
	movss  [edi + ecx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [esi + edx*4],    xmm3
	movss  [esi + edx*4 + 8], xmm4

	;# get n from stack
	mov esi, [esp + nb301_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb301_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb301_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 

        ;# finish if last 
        mov ecx, [esp + nb301_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301_n], esi
        jmp .nb301_outer
.nb301_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301_end
        ;# non-zero, do one more workunit
        jmp   .nb301_threadloop
.nb301_end:
	emms

	mov eax, [esp + nb301_nouter]
	mov ebx, [esp + nb301_ninner]
	mov ecx, [ebp + nb301_outeriter]
	mov edx, [ebp + nb301_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb301_salign]
	add esp, eax
	add esp, 760
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
	




.globl nb_kernel301nf_ia32_sse
.globl _nb_kernel301nf_ia32_sse
nb_kernel301nf_ia32_sse:	
_nb_kernel301nf_ia32_sse:	
.equiv          nb301nf_p_nri,          8
.equiv          nb301nf_iinr,           12
.equiv          nb301nf_jindex,         16
.equiv          nb301nf_jjnr,           20
.equiv          nb301nf_shift,          24
.equiv          nb301nf_shiftvec,       28
.equiv          nb301nf_fshift,         32
.equiv          nb301nf_gid,            36
.equiv          nb301nf_pos,            40
.equiv          nb301nf_faction,        44
.equiv          nb301nf_charge,         48
.equiv          nb301nf_p_facel,        52
.equiv          nb301nf_argkrf,         56
.equiv          nb301nf_argcrf,         60
.equiv          nb301nf_Vc,             64
.equiv          nb301nf_type,           68
.equiv          nb301nf_p_ntype,        72
.equiv          nb301nf_vdwparam,       76
.equiv          nb301nf_Vvdw,           80
.equiv          nb301nf_p_tabscale,     84
.equiv          nb301nf_VFtab,          88
.equiv          nb301nf_invsqrta,       92
.equiv          nb301nf_dvda,           96
.equiv          nb301nf_p_gbtabscale,   100
.equiv          nb301nf_GBtab,          104
.equiv          nb301nf_p_nthreads,     108
.equiv          nb301nf_count,          112
.equiv          nb301nf_mtx,            116
.equiv          nb301nf_outeriter,      120
.equiv          nb301nf_inneriter,      124
.equiv          nb301nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb301nf_ixO,            0
.equiv          nb301nf_iyO,            16
.equiv          nb301nf_izO,            32
.equiv          nb301nf_ixH1,           48
.equiv          nb301nf_iyH1,           64
.equiv          nb301nf_izH1,           80
.equiv          nb301nf_ixH2,           96
.equiv          nb301nf_iyH2,           112
.equiv          nb301nf_izH2,           128
.equiv          nb301nf_iqO,            144
.equiv          nb301nf_iqH,            160
.equiv          nb301nf_qqO,            176
.equiv          nb301nf_qqH,            192
.equiv          nb301nf_rinvO,          208
.equiv          nb301nf_rinvH1,         224
.equiv          nb301nf_rinvH2,         240
.equiv          nb301nf_rO,             256
.equiv          nb301nf_rH1,            272
.equiv          nb301nf_rH2,            288
.equiv          nb301nf_tsc,            304
.equiv          nb301nf_vctot,          320
.equiv          nb301nf_half,           336
.equiv          nb301nf_three,          352
.equiv          nb301nf_is3,            368
.equiv          nb301nf_ii3,            372
.equiv          nb301nf_innerjjnr,      376
.equiv          nb301nf_innerk,         380
.equiv          nb301nf_n,              384
.equiv          nb301nf_nn1,            388
.equiv          nb301nf_nri,            392
.equiv          nb301nf_nouter,         396
.equiv          nb301nf_ninner,         400
.equiv          nb301nf_salign,         404
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 408		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb301nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb301nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb301nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb301nf_nouter], eax
	mov [esp + nb301nf_ninner], eax


	mov eax, [ebp + nb301nf_p_tabscale]
	movss xmm3, [eax]
	shufps xmm3, xmm3, 0 
	movaps [esp + nb301nf_tsc], xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb301nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb301nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb301nf_iqO], xmm3
	movaps [esp + nb301nf_iqH], xmm4
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb301nf_half], eax
	movss xmm1, [esp + nb301nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb301nf_half],  xmm1
	movaps [esp + nb301nf_three],  xmm3

.nb301nf_threadloop:
        mov   esi, [ebp + nb301nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb301nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb301nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb301nf_n], eax
        mov [esp + nb301nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301nf_outerstart
        jmp .nb301nf_end
			
.nb301nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb301nf_nouter]
	mov [esp + nb301nf_nouter], ebx

.nb301nf_outer:
	mov   eax, [ebp + nb301nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb301nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb301nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb301nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb301nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb301nf_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb301nf_ixO], xmm3
	movaps [esp + nb301nf_iyO], xmm4
	movaps [esp + nb301nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 12]
	addss xmm1, [eax + ebx*4 + 16]
	addss xmm2, [eax + ebx*4 + 20]		
	addss xmm3, [eax + ebx*4 + 24]
	addss xmm4, [eax + ebx*4 + 28]
	addss xmm5, [eax + ebx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb301nf_ixH1], xmm0
	movaps [esp + nb301nf_iyH1], xmm1
	movaps [esp + nb301nf_izH1], xmm2
	movaps [esp + nb301nf_ixH2], xmm3
	movaps [esp + nb301nf_iyH2], xmm4
	movaps [esp + nb301nf_izH2], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb301nf_vctot], xmm4
	
	mov   eax, [ebp + nb301nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb301nf_pos]
	mov   eax, [ebp + nb301nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb301nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb301nf_ninner]
	mov   [esp + nb301nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb301nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301nf_unroll_loop
	jmp   .nb301nf_odd_inner
.nb301nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb301nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb301nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb301nf_iqO]
	mulps  xmm4, [esp + nb301nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [esp + nb301nf_qqO], xmm3
	movaps  [esp + nb301nf_qqH], xmm4	

	mov esi, [ebp + nb301nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [esp + nb301nf_ixO]
	movaps xmm5, [esp + nb301nf_iyO]
	movaps xmm6, [esp + nb301nf_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2
	
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb301nf_ixH1]
	movaps xmm5, [esp + nb301nf_iyH1]
	movaps xmm6, [esp + nb301nf_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb301nf_ixH2]
	movaps xmm4, [esp + nb301nf_iyH2]
	movaps xmm5, [esp + nb301nf_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301nf_half]
	movaps  [esp + nb301nf_rinvO], xmm4	;# rinvO in xmm4 
	mulps   xmm7, xmm4
	movaps  [esp + nb301nf_rO], xmm7	

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301nf_half]
	movaps  [esp + nb301nf_rinvH1], xmm4	;# rinvH1 in xmm4 
	mulps   xmm6, xmm4
	movaps  [esp + nb301nf_rH1], xmm6

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb301nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb301nf_half]
	movaps  [esp + nb301nf_rinvH2], xmm4	;# rinvH2 in xmm4 
	mulps   xmm5, xmm4
	movaps  [esp + nb301nf_rH2], xmm5

	;# do O interactions 
	;# rO is still in xmm7 
	mulps   xmm7, [esp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3

	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd mm0, eax   
    	movd mm1, ebx
    	movd mm2, ecx
    	movd mm3, edx

    	mov  esi, [ebp + nb301nf_VFtab]
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp    
    	movaps xmm0, [esp + nb301nf_qqO]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  

    	;# at this point mm5 contains vcoul 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [esp + nb301nf_vctot]
    	movaps [esp + nb301nf_vctot], xmm5
	
	;# Done with O interactions - now H1! 
	movaps xmm7, [esp + nb301nf_rH1]
	mulps   xmm7, [esp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [esp + nb301nf_qqH]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul 
    	addps  xmm5, [esp + nb301nf_vctot]
    	movaps [esp + nb301nf_vctot], xmm5 

	;# Done with H1, finally we do H2 interactions 
	movaps xmm7, [esp + nb301nf_rH2]
	mulps   xmm7, [esp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
		
    	movd eax, mm6
    	psrlq mm6, 32
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd ebx, mm6
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm5, [esi + ebx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm7, [esi + ebx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [esp + nb301nf_qqH]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul 
    	addps  xmm5, [esp + nb301nf_vctot]
    	movaps [esp + nb301nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb301nf_innerk],  4
	jl    .nb301nf_odd_inner
	jmp   .nb301nf_unroll_loop
.nb301nf_odd_inner:	
	add dword ptr [esp + nb301nf_innerk],  4
	jnz   .nb301nf_odd_loop
	jmp   .nb301nf_updateouterdata
.nb301nf_odd_loop:
	mov   edx, [esp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb301nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb301nf_iqO]
	mov esi, [ebp + nb301nf_charge] 
	movhps xmm4, [esp + nb301nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb301nf_qqO], xmm3	;# use oxygen qq for storage 

	mov esi, [ebp + nb301nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb301nf_ixO]
	movss xmm4, [esp + nb301nf_iyO]
	movss xmm5, [esp + nb301nf_izO]
		
	movlps xmm6, [esp + nb301nf_ixH1]
	movlps xmm7, [esp + nb301nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb301nf_iyH1]
	movlps xmm7, [esp + nb301nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb301nf_izH1]
	movlps xmm7, [esp + nb301nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb301nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb301nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000	
	
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [esp + nb301nf_rinvO], xmm0
	
	mulps xmm4, [esp + nb301nf_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ;# mm6/mm7 contain lu indices 
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm7, mm7
    	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
    	movd mm0, eax   
    	movd mm1, ecx
    	movd mm2, edx

    	mov  esi, [ebp + nb301nf_VFtab]
    	movd eax, mm6
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd edx, mm7

    	movlps xmm5, [esi + eax*4]
    	movlps xmm7, [esi + ecx*4]
    	movhps xmm7, [esi + edx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# constant 10001000
    	shufps xmm5, xmm7, 221  ;# constant 11011101

    	movlps xmm7, [esi + eax*4 + 8]
    	movlps xmm3, [esi + ecx*4 + 8]
    	movhps xmm3, [esi + edx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# constant 10001000
    	shufps xmm7, xmm3, 221  ;# constant 11011101
    	;# coulomb table ready, in xmm4-xmm7      
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [esp + nb301nf_qqO]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [esp + nb301nf_vctot]
    	movaps [esp + nb301nf_vctot], xmm5

	dec dword ptr [esp + nb301nf_innerk]
	jz    .nb301nf_updateouterdata
	jmp   .nb301nf_odd_loop
.nb301nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb301nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb301nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb301nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb301nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 

        ;# finish if last 
        mov ecx, [esp + nb301nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb301nf_n], esi
        jmp .nb301nf_outer
.nb301nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb301nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301nf_end
        ;# non-zero, do one more workunit
        jmp   .nb301nf_threadloop
.nb301nf_end:
	emms

	mov eax, [esp + nb301nf_nouter]
	mov ebx, [esp + nb301nf_ninner]
	mov ecx, [ebp + nb301nf_outeriter]
	mov edx, [ebp + nb301nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb301nf_salign]
	add esp, eax
	add esp, 408
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
	
