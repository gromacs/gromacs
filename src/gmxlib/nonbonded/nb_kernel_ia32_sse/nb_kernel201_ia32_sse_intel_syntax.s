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




	

	

.globl nb_kernel201_ia32_sse
.globl _nb_kernel201_ia32_sse
nb_kernel201_ia32_sse:	
_nb_kernel201_ia32_sse:	
.equiv          nb201_p_nri,            8
.equiv          nb201_iinr,             12
.equiv          nb201_jindex,           16
.equiv          nb201_jjnr,             20
.equiv          nb201_shift,            24
.equiv          nb201_shiftvec,         28
.equiv          nb201_fshift,           32
.equiv          nb201_gid,              36
.equiv          nb201_pos,              40
.equiv          nb201_faction,          44
.equiv          nb201_charge,           48
.equiv          nb201_p_facel,          52
.equiv          nb201_argkrf,           56
.equiv          nb201_argcrf,           60
.equiv          nb201_Vc,               64
.equiv          nb201_type,             68
.equiv          nb201_p_ntype,          72
.equiv          nb201_vdwparam,         76
.equiv          nb201_Vvdw,             80
.equiv          nb201_p_tabscale,       84
.equiv          nb201_VFtab,            88
.equiv          nb201_invsqrta,         92
.equiv          nb201_dvda,             96
.equiv          nb201_p_gbtabscale,     100
.equiv          nb201_GBtab,            104
.equiv          nb201_p_nthreads,       108
.equiv          nb201_count,            112
.equiv          nb201_mtx,              116
.equiv          nb201_outeriter,        120
.equiv          nb201_inneriter,        124
.equiv          nb201_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb201_ixO,              0
.equiv          nb201_iyO,              16
.equiv          nb201_izO,              32
.equiv          nb201_ixH1,             48
.equiv          nb201_iyH1,             64
.equiv          nb201_izH1,             80
.equiv          nb201_ixH2,             96
.equiv          nb201_iyH2,             112
.equiv          nb201_izH2,             128
.equiv          nb201_iqO,              144
.equiv          nb201_iqH,              160
.equiv          nb201_dxO,              176
.equiv          nb201_dyO,              192
.equiv          nb201_dzO,              208
.equiv          nb201_dxH1,             224
.equiv          nb201_dyH1,             240
.equiv          nb201_dzH1,             256
.equiv          nb201_dxH2,             272
.equiv          nb201_dyH2,             288
.equiv          nb201_dzH2,             304
.equiv          nb201_qqO,              320
.equiv          nb201_qqH,              336
.equiv          nb201_vctot,            352
.equiv          nb201_fixO,             384
.equiv          nb201_fiyO,             400
.equiv          nb201_fizO,             416
.equiv          nb201_fixH1,            432
.equiv          nb201_fiyH1,            448
.equiv          nb201_fizH1,            464
.equiv          nb201_fixH2,            480
.equiv          nb201_fiyH2,            496
.equiv          nb201_fizH2,            512
.equiv          nb201_fjx,              528
.equiv          nb201_fjy,              544
.equiv          nb201_fjz,              560
.equiv          nb201_half,             576
.equiv          nb201_three,            592
.equiv          nb201_two,              608
.equiv          nb201_krf,              624
.equiv          nb201_crf,              640
.equiv          nb201_krsqO,            656
.equiv          nb201_krsqH1,           672
.equiv          nb201_krsqH2,           688
.equiv          nb201_is3,              704
.equiv          nb201_ii3,              708
.equiv          nb201_innerjjnr,        712
.equiv          nb201_innerk,           716
.equiv          nb201_n,                720
.equiv          nb201_nn1,              724
.equiv          nb201_nri,              728
.equiv          nb201_nouter,           732
.equiv          nb201_ninner,           736
.equiv          nb201_salign,           740
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 744		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb201_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb201_p_nri]
	mov ecx, [ecx]
	mov [esp + nb201_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb201_nouter], eax
	mov [esp + nb201_ninner], eax


	mov esi, [ebp + nb201_argkrf]
	mov edi, [ebp + nb201_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb201_krf], xmm5
	movaps [esp + nb201_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb201_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb201_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb201_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb201_iqO], xmm3
	movaps [esp + nb201_iqH], xmm4
			
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb201_half], eax
	movss xmm1, [esp + nb201_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb201_half],  xmm1
	movaps [esp + nb201_two],  xmm2
	movaps [esp + nb201_three],  xmm3


.nb201_threadloop:
        mov   esi, [ebp + nb201_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb201_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb201_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb201_n], eax
        mov [esp + nb201_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201_outerstart
        jmp .nb201_end
	
.nb201_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb201_nouter]
	mov [esp + nb201_nouter], ebx

.nb201_outer:
	mov   eax, [ebp + nb201_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb201_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb201_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb201_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb201_pos]    ;# eax = base of pos[]  
	mov   [esp + nb201_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb201_ixO], xmm3
	movaps [esp + nb201_iyO], xmm4
	movaps [esp + nb201_izO], xmm5

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
	movaps [esp + nb201_ixH1], xmm0
	movaps [esp + nb201_iyH1], xmm1
	movaps [esp + nb201_izH1], xmm2
	movaps [esp + nb201_ixH2], xmm3
	movaps [esp + nb201_iyH2], xmm4
	movaps [esp + nb201_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb201_vctot], xmm4
	movaps [esp + nb201_fixO], xmm4
	movaps [esp + nb201_fiyO], xmm4
	movaps [esp + nb201_fizO], xmm4
	movaps [esp + nb201_fixH1], xmm4
	movaps [esp + nb201_fiyH1], xmm4
	movaps [esp + nb201_fizH1], xmm4
	movaps [esp + nb201_fixH2], xmm4
	movaps [esp + nb201_fiyH2], xmm4
	movaps [esp + nb201_fizH2], xmm4
	
	mov   eax, [ebp + nb201_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb201_pos]
	mov   edi, [ebp + nb201_faction]	
	mov   eax, [ebp + nb201_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb201_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb201_ninner]
	mov   [esp + nb201_ninner], ecx
	add   edx, 0
	mov   [esp + nb201_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201_unroll_loop
	jmp   .nb201_odd_inner
.nb201_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb201_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb201_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb201_iqO]
	mulps  xmm4, [esp + nb201_iqH]

	movaps  [esp + nb201_qqO], xmm3
	movaps  [esp + nb201_qqH], xmm4

	mov esi, [ebp + nb201_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb201_ixO]
	movaps xmm5, [esp + nb201_iyO]
	movaps xmm6, [esp + nb201_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb201_dxO], xmm4
	movaps [esp + nb201_dyO], xmm5
	movaps [esp + nb201_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb201_ixH1]
	movaps xmm5, [esp + nb201_iyH1]
	movaps xmm6, [esp + nb201_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb201_dxH1], xmm4
	movaps [esp + nb201_dyH1], xmm5
	movaps [esp + nb201_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [esp + nb201_ixH2]
	movaps xmm4, [esp + nb201_iyH2]
	movaps xmm5, [esp + nb201_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb201_dxH2], xmm3
	movaps [esp + nb201_dyH2], xmm4
	movaps [esp + nb201_dzH2], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb201_krf]	
	mulps  xmm1, [esp + nb201_krf]	
	mulps  xmm2, [esp + nb201_krf]	

	movaps [esp + nb201_krsqH2], xmm0
	movaps [esp + nb201_krsqH1], xmm1
	movaps [esp + nb201_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 

	movaps xmm0, xmm7
	movaps xmm1, [esp + nb201_krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb201_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm1, [esp + nb201_two]
	subps  xmm7, xmm1
	mulps  xmm0, [esp + nb201_qqO]
	mulps  xmm7, [esp + nb201_qqO]

	mulps  xmm4, xmm7	;# total fsO in xmm4 

	addps  xmm0, [esp + nb201_vctot]
	movaps [esp + nb201_vctot], xmm0

	movaps xmm0, [esp + nb201_dxO]
	movaps xmm1, [esp + nb201_dyO]
	movaps xmm2, [esp + nb201_dzO]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update O forces 
	movaps xmm3, [esp + nb201_fixO]
	movaps xmm4, [esp + nb201_fiyO]
	movaps xmm7, [esp + nb201_fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb201_fixO], xmm3
	movaps [esp + nb201_fiyO], xmm4
	movaps [esp + nb201_fizO], xmm7
	;# update j forces with water O 
	movaps [esp + nb201_fjx], xmm0
	movaps [esp + nb201_fjy], xmm1
	movaps [esp + nb201_fjz], xmm2

	;# H1 interactions 
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + nb201_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb201_crf] ;# xmm6=rinv+ krsq-crf 
	mulps   xmm0, [esp + nb201_two]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm6, [esp + nb201_qqH] ;# vcoul 
	mulps   xmm7, [esp + nb201_qqH]
	mulps  xmm4, xmm7		;# total fsH1 in xmm4 
	
	addps  xmm6, [esp + nb201_vctot]

	movaps xmm0, [esp + nb201_dxH1]
	movaps xmm1, [esp + nb201_dyH1]
	movaps xmm2, [esp + nb201_dzH1]
	movaps [esp + nb201_vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb201_fixH1]
	movaps xmm4, [esp + nb201_fiyH1]
	movaps xmm7, [esp + nb201_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb201_fixH1], xmm3
	movaps [esp + nb201_fiyH1], xmm4
	movaps [esp + nb201_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [esp + nb201_fjx]
	addps  xmm1, [esp + nb201_fjy]
	addps  xmm2, [esp + nb201_fjz]
	movaps [esp + nb201_fjx], xmm0
	movaps [esp + nb201_fjy], xmm1
	movaps [esp + nb201_fjz], xmm2

	;# H2 interactions 
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + nb201_krsqH2]
	addps   xmm5, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm5, [esp + nb201_crf] ;# xmm5=rinv+ krsq-crf 
	mulps   xmm0, [esp + nb201_two]
	subps   xmm7, xmm0	;# xmm7=rinv-2*krsq 
	mulps   xmm5, [esp + nb201_qqH] ;# vcoul 
	mulps   xmm7, [esp + nb201_qqH]
	mulps  xmm4, xmm7		;# total fsH2 in xmm4 
	
	addps  xmm5, [esp + nb201_vctot]

	movaps xmm0, [esp + nb201_dxH2]
	movaps xmm1, [esp + nb201_dyH2]
	movaps xmm2, [esp + nb201_dzH2]
	movaps [esp + nb201_vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H2 forces 
	movaps xmm3, [esp + nb201_fixH2]
	movaps xmm4, [esp + nb201_fiyH2]
	movaps xmm7, [esp + nb201_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb201_fixH2], xmm3
	movaps [esp + nb201_fiyH2], xmm4
	movaps [esp + nb201_fizH2], xmm7

	mov edi, [ebp + nb201_faction]
	;# update j forces 
	addps xmm0, [esp + nb201_fjx]
	addps xmm1, [esp + nb201_fjy]
	addps xmm2, [esp + nb201_fjz]

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
	sub dword ptr [esp + nb201_innerk],  4
	jl    .nb201_odd_inner
	jmp   .nb201_unroll_loop
.nb201_odd_inner:	
	add dword ptr [esp + nb201_innerk],  4
	jnz   .nb201_odd_loop
	jmp   .nb201_updateouterdata
.nb201_odd_loop:
	mov   edx, [esp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb201_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb201_iqO]
	mov esi, [ebp + nb201_charge] 
	movhps xmm4, [esp + nb201_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb201_qqO], xmm3	;# use oxygen qq for storage 

	mov esi, [ebp + nb201_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb201_ixO]
	movss xmm4, [esp + nb201_iyO]
	movss xmm5, [esp + nb201_izO]
		
	movlps xmm6, [esp + nb201_ixH1]
	movlps xmm7, [esp + nb201_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb201_iyH1]
	movlps xmm7, [esp + nb201_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb201_izH1]
	movlps xmm7, [esp + nb201_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [esp + nb201_dxO], xmm3
	movaps [esp + nb201_dyO], xmm4
	movaps [esp + nb201_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb201_krf]
	movaps [esp + nb201_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb201_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb201_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000

	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [esp + nb201_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [esp + nb201_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [esp + nb201_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [esp + nb201_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [esp + nb201_qqO] 	;# xmm1=coul part of fs 

	
	mulps  xmm4, xmm1	;# xmm4=total fscal 
	addps  xmm0, [esp + nb201_vctot]
	movaps [esp + nb201_vctot], xmm0
	
	movaps xmm0, [esp + nb201_dxO]
	movaps xmm1, [esp + nb201_dyO]
	movaps xmm2, [esp + nb201_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [esp + nb201_fixO]	
	movss  xmm4, [esp + nb201_fiyO]	
	movss  xmm5, [esp + nb201_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb201_fixO], xmm3	
	movss  [esp + nb201_fiyO], xmm4	
	movss  [esp + nb201_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb201_fixH1]
	addss  xmm4, [esp + nb201_fiyH1]
	addss  xmm5, [esp + nb201_fizH1]
	movss  [esp + nb201_fixH1], xmm3	
	movss  [esp + nb201_fiyH1], xmm4	
	movss  [esp + nb201_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb201_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb201_fixH2]
	addss  xmm4, [esp + nb201_fiyH2]
	addss  xmm5, [esp + nb201_fizH2]
	movss  [esp + nb201_fixH2], xmm3	
	movss  [esp + nb201_fiyH2], xmm4	
	movss  [esp + nb201_fizH2], xmm5	;# updated the H2 force 

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

	dec dword ptr [esp + nb201_innerk]
	jz    .nb201_updateouterdata
	jmp   .nb201_odd_loop
.nb201_updateouterdata:
	mov   ecx, [esp + nb201_ii3]
	mov   edi, [ebp + nb201_faction]
	mov   esi, [ebp + nb201_fshift]
	mov   edx, [esp + nb201_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb201_fixO]
	movaps xmm1, [esp + nb201_fiyO]
	movaps xmm2, [esp + nb201_fizO]

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
	movaps xmm0, [esp + nb201_fixH1]
	movaps xmm1, [esp + nb201_fiyH1]
	movaps xmm2, [esp + nb201_fizH1]

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
	movaps xmm0, [esp + nb201_fixH2]
	movaps xmm1, [esp + nb201_fiyH2]
	movaps xmm2, [esp + nb201_fizH2]

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
	mov esi, [esp + nb201_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb201_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb201_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb201_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb201_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb201_n], esi
        jmp .nb201_outer
.nb201_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb201_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201_end
        ;# non-zero, do one more workunit
        jmp   .nb201_threadloop
.nb201_end:
	emms

	mov eax, [esp + nb201_nouter]
	mov ebx, [esp + nb201_ninner]
	mov ecx, [ebp + nb201_outeriter]
	mov edx, [ebp + nb201_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb201_salign]
	add esp, eax
	add esp, 744
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


	

.globl nb_kernel201nf_ia32_sse
.globl _nb_kernel201nf_ia32_sse
nb_kernel201nf_ia32_sse:	
_nb_kernel201nf_ia32_sse:	
.equiv          nb201nf_p_nri,          8
.equiv          nb201nf_iinr,           12
.equiv          nb201nf_jindex,         16
.equiv          nb201nf_jjnr,           20
.equiv          nb201nf_shift,          24
.equiv          nb201nf_shiftvec,       28
.equiv          nb201nf_fshift,         32
.equiv          nb201nf_gid,            36
.equiv          nb201nf_pos,            40
.equiv          nb201nf_faction,        44
.equiv          nb201nf_charge,         48
.equiv          nb201nf_p_facel,        52
.equiv          nb201nf_argkrf,         56
.equiv          nb201nf_argcrf,         60
.equiv          nb201nf_Vc,             64
.equiv          nb201nf_type,           68
.equiv          nb201nf_p_ntype,        72
.equiv          nb201nf_vdwparam,       76
.equiv          nb201nf_Vvdw,           80
.equiv          nb201nf_p_tabscale,     84
.equiv          nb201nf_VFtab,          88
.equiv          nb201nf_invsqrta,       92
.equiv          nb201nf_dvda,           96
.equiv          nb201nf_p_gbtabscale,   100
.equiv          nb201nf_GBtab,          104
.equiv          nb201nf_p_nthreads,     108
.equiv          nb201nf_count,          112
.equiv          nb201nf_mtx,            116
.equiv          nb201nf_outeriter,      120
.equiv          nb201nf_inneriter,      124
.equiv          nb201nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb201nf_ixO,            0
.equiv          nb201nf_iyO,            16
.equiv          nb201nf_izO,            32
.equiv          nb201nf_ixH1,           48
.equiv          nb201nf_iyH1,           64
.equiv          nb201nf_izH1,           80
.equiv          nb201nf_ixH2,           96
.equiv          nb201nf_iyH2,           112
.equiv          nb201nf_izH2,           128
.equiv          nb201nf_iqO,            144
.equiv          nb201nf_iqH,            160
.equiv          nb201nf_qqO,            176
.equiv          nb201nf_qqH,            192
.equiv          nb201nf_vctot,          208
.equiv          nb201nf_half,           224
.equiv          nb201nf_three,          240
.equiv          nb201nf_krf,            256
.equiv          nb201nf_crf,            272
.equiv          nb201nf_krsqO,          288
.equiv          nb201nf_krsqH1,         304
.equiv          nb201nf_krsqH2,         320
.equiv          nb201nf_is3,            336
.equiv          nb201nf_ii3,            340
.equiv          nb201nf_innerjjnr,      344
.equiv          nb201nf_innerk,         348
.equiv          nb201nf_n,              352
.equiv          nb201nf_nn1,            356
.equiv          nb201nf_nri,            360
.equiv          nb201nf_nouter,         364
.equiv          nb201nf_ninner,         368
.equiv          nb201nf_salign,         372
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 376		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb201nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb201nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb201nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb201nf_nouter], eax
	mov [esp + nb201nf_ninner], eax

	mov esi, [ebp + nb201nf_argkrf]
	mov edi, [ebp + nb201nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	shufps xmm6, xmm6, 0
	movaps [esp + nb201nf_krf], xmm5
	movaps [esp + nb201nf_crf], xmm6
	
	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb201nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + nb201nf_charge]
	movss xmm3, [edx + ebx*4]	
	movss xmm4, [edx + ebx*4 + 4]	
	mov esi, [ebp + nb201nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb201nf_iqO], xmm3
	movaps [esp + nb201nf_iqH], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb201nf_half], eax
	movss xmm1, [esp + nb201nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb201nf_half],  xmm1
	movaps [esp + nb201nf_three],  xmm3

.nb201nf_threadloop:
        mov   esi, [ebp + nb201nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb201nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb201nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb201nf_n], eax
        mov [esp + nb201nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201nf_outerstart
        jmp .nb201nf_end
			
.nb201nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb201nf_nouter]
	mov [esp + nb201nf_nouter], ebx

.nb201nf_outer:
	mov   eax, [ebp + nb201nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb201nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb201nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb201nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb201nf_pos]    ;# eax = base of pos[]  
	mov   [esp + nb201nf_ii3], ebx

	addss xmm3, [eax + ebx*4]
	addss xmm4, [eax + ebx*4 + 4]
	addss xmm5, [eax + ebx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb201nf_ixO], xmm3
	movaps [esp + nb201nf_iyO], xmm4
	movaps [esp + nb201nf_izO], xmm5

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
	movaps [esp + nb201nf_ixH1], xmm0
	movaps [esp + nb201nf_iyH1], xmm1
	movaps [esp + nb201nf_izH1], xmm2
	movaps [esp + nb201nf_ixH2], xmm3
	movaps [esp + nb201nf_iyH2], xmm4
	movaps [esp + nb201nf_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb201nf_vctot], xmm4
	
	mov   eax, [ebp + nb201nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb201nf_pos]
	mov   eax, [ebp + nb201nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb201nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb201nf_ninner]
	mov   [esp + nb201nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb201nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201nf_unroll_loop
	jmp   .nb201nf_odd_inner
.nb201nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 

	add dword ptr [esp + nb201nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb201nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [esp + nb201nf_iqO]
	mulps  xmm4, [esp + nb201nf_iqH]

	movaps  [esp + nb201nf_qqO], xmm3
	movaps  [esp + nb201nf_qqH], xmm4

	mov esi, [ebp + nb201nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb201nf_ixO]
	movaps xmm5, [esp + nb201nf_iyO]
	movaps xmm6, [esp + nb201nf_izO]

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
	movaps xmm4, [esp + nb201nf_ixH1]
	movaps xmm5, [esp + nb201nf_iyH1]
	movaps xmm6, [esp + nb201nf_izH1]

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
	movaps xmm3, [esp + nb201nf_ixH2]
	movaps xmm4, [esp + nb201nf_iyH2]
	movaps xmm5, [esp + nb201nf_izH2]

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

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [esp + nb201nf_krf]	
	mulps  xmm1, [esp + nb201nf_krf]	
	mulps  xmm2, [esp + nb201nf_krf]	

	movaps [esp + nb201nf_krsqH2], xmm0
	movaps [esp + nb201nf_krsqH1], xmm1
	movaps [esp + nb201nf_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201nf_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201nf_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb201nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb201nf_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 

	movaps xmm0, xmm7
	movaps xmm1, [esp + nb201nf_krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [esp + nb201nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb201nf_qqO]

	addps  xmm0, [esp + nb201nf_vctot]
	movaps [esp + nb201nf_vctot], xmm0

	;# H1 interactions 
	movaps  xmm7, xmm6
	movaps  xmm0, [esp + nb201nf_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [esp + nb201nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps   xmm6, [esp + nb201nf_qqH] ;# vcoul 
	addps   xmm6, [esp + nb201nf_vctot]
	
	;# H2 interactions 
	movaps  xmm7, xmm5
	movaps  xmm0, [esp + nb201nf_krsqH2]
	addps   xmm5, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm5, [esp + nb201nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps   xmm5, [esp + nb201nf_qqH] ;# vcoul 
	addps  xmm6, xmm5
	movaps [esp + nb201nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [esp + nb201nf_innerk],  4
	jl    .nb201nf_odd_inner
	jmp   .nb201nf_unroll_loop
.nb201nf_odd_inner:	
	add dword ptr [esp + nb201nf_innerk],  4
	jnz   .nb201nf_odd_loop
	jmp   .nb201nf_updateouterdata
.nb201nf_odd_loop:
	mov   edx, [esp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb201nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb201nf_iqO]
	mov esi, [ebp + nb201nf_charge] 
	movhps xmm4, [esp + nb201nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb201nf_qqO], xmm3	;# use oxygen qq for storage 

	mov esi, [ebp + nb201nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb201nf_ixO]
	movss xmm4, [esp + nb201nf_iyO]
	movss xmm5, [esp + nb201nf_izO]
		
	movlps xmm6, [esp + nb201nf_ixH1]
	movlps xmm7, [esp + nb201nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb201nf_iyH1]
	movlps xmm7, [esp + nb201nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb201nf_izH1]
	movlps xmm7, [esp + nb201nf_izH2]
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

	movaps xmm0, xmm4
	mulps xmm0, [esp + nb201nf_krf]
	movaps [esp + nb201nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb201nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb201nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000	
	
	movaps xmm3, [esp + nb201nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [esp + nb201nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [esp + nb201nf_qqO]	;# xmm0=vcoul 
	addps  xmm0, [esp + nb201nf_vctot]
	movaps [esp + nb201nf_vctot], xmm0

	dec dword ptr [esp + nb201nf_innerk]
	jz    .nb201nf_updateouterdata
	jmp   .nb201nf_odd_loop
.nb201nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb201nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb201nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb201nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb201nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb201nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb201nf_n], esi
        jmp .nb201nf_outer
.nb201nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb201nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201nf_end
        ;# non-zero, do one more workunit
        jmp   .nb201nf_threadloop
.nb201nf_end:
	emms

	mov eax, [esp + nb201nf_nouter]
	mov ebx, [esp + nb201nf_ninner]
	mov ecx, [ebp + nb201nf_outeriter]
	mov edx, [ebp + nb201nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb201nf_salign]
	add esp, eax
	add esp, 376
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
