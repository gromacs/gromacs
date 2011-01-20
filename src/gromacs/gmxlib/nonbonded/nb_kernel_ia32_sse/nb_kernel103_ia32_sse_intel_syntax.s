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


.globl nb_kernel103_ia32_sse
.globl _nb_kernel103_ia32_sse
nb_kernel103_ia32_sse:	
_nb_kernel103_ia32_sse:	
.equiv          nb103_p_nri,            8
.equiv          nb103_iinr,             12
.equiv          nb103_jindex,           16
.equiv          nb103_jjnr,             20
.equiv          nb103_shift,            24
.equiv          nb103_shiftvec,         28
.equiv          nb103_fshift,           32
.equiv          nb103_gid,              36
.equiv          nb103_pos,              40
.equiv          nb103_faction,          44
.equiv          nb103_charge,           48
.equiv          nb103_p_facel,          52
.equiv          nb103_p_krf,            56
.equiv          nb103_p_crf,            60
.equiv          nb103_Vc,               64
.equiv          nb103_type,             68
.equiv          nb103_p_ntype,          72
.equiv          nb103_vdwparam,         76
.equiv          nb103_Vvdw,             80
.equiv          nb103_p_tabscale,       84
.equiv          nb103_VFtab,            88
.equiv          nb103_invsqrta,         92
.equiv          nb103_dvda,             96
.equiv          nb103_p_gbtabscale,     100
.equiv          nb103_GBtab,            104
.equiv          nb103_p_nthreads,       108
.equiv          nb103_count,            112
.equiv          nb103_mtx,              116
.equiv          nb103_outeriter,        120
.equiv          nb103_inneriter,        124
.equiv          nb103_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb103_ixH1,             0
.equiv          nb103_iyH1,             16
.equiv          nb103_izH1,             32
.equiv          nb103_ixH2,             48
.equiv          nb103_iyH2,             64
.equiv          nb103_izH2,             80
.equiv          nb103_ixM,              96
.equiv          nb103_iyM,              112
.equiv          nb103_izM,              128
.equiv          nb103_iqH,              144
.equiv          nb103_iqM,              160
.equiv          nb103_dxH1,             176
.equiv          nb103_dyH1,             192
.equiv          nb103_dzH1,             208
.equiv          nb103_dxH2,             224
.equiv          nb103_dyH2,             240
.equiv          nb103_dzH2,             256
.equiv          nb103_dxM,              272
.equiv          nb103_dyM,              288
.equiv          nb103_dzM,              304
.equiv          nb103_qqH,              320
.equiv          nb103_qqM,              336
.equiv          nb103_vctot,            352
.equiv          nb103_fixH1,            368
.equiv          nb103_fiyH1,            384
.equiv          nb103_fizH1,            400
.equiv          nb103_fixH2,            416
.equiv          nb103_fiyH2,            432
.equiv          nb103_fizH2,            448
.equiv          nb103_fixM,             464
.equiv          nb103_fiyM,             480
.equiv          nb103_fizM,             496
.equiv          nb103_fjx,              512
.equiv          nb103_fjy,              528
.equiv          nb103_fjz,              544
.equiv          nb103_half,             560
.equiv          nb103_three,            576
.equiv          nb103_is3,              592
.equiv          nb103_ii3,              596
.equiv          nb103_innerjjnr,        600
.equiv          nb103_innerk,           604
.equiv          nb103_n,                608
.equiv          nb103_nn1,              612
.equiv          nb103_nri,              616
.equiv          nb103_nouter,           620
.equiv          nb103_ninner,           624
.equiv          nb103_salign,           628
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 632		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb103_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb103_p_nri]
	mov ecx, [ecx]
	mov [esp + nb103_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb103_nouter], eax
	mov [esp + nb103_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb103_half], eax
	movss xmm1, [esp + nb103_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb103_half],  xmm1
	movaps [esp + nb103_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb103_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb103_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb103_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb103_iqH], xmm3
	movaps [esp + nb103_iqM], xmm4
	
.nb103_threadloop:
        mov   esi, [ebp + nb103_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb103_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                            ;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb103_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb103_n], eax
        mov [esp + nb103_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb103_outerstart
        jmp .nb103_end
	
.nb103_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb103_nouter]
	mov [esp + nb103_nouter], ebx

.nb103_outer:
	mov   eax, [ebp + nb103_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb103_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb103_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb103_iinr]   	;# ecx = pointer into iinr[] 
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb103_pos]	;# eax = base of pos[]  
	mov   [esp + nb103_ii3], ebx

	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]	
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb103_ixH1], xmm3
	movaps [esp + nb103_iyH1], xmm4
	movaps [esp + nb103_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb103_ixH2], xmm0
	movaps [esp + nb103_iyH2], xmm1
	movaps [esp + nb103_izH2], xmm2
	movaps [esp + nb103_ixM], xmm3
	movaps [esp + nb103_iyM], xmm4
	movaps [esp + nb103_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb103_vctot], xmm4
	movaps [esp + nb103_fixH1], xmm4
	movaps [esp + nb103_fiyH1], xmm4
	movaps [esp + nb103_fizH1], xmm4
	movaps [esp + nb103_fixH2], xmm4
	movaps [esp + nb103_fiyH2], xmm4
	movaps [esp + nb103_fizH2], xmm4
	movaps [esp + nb103_fixM], xmm4
	movaps [esp + nb103_fiyM], xmm4
	movaps [esp + nb103_fizM], xmm4
	
	mov   eax, [ebp + nb103_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb103_pos]
	mov   edi, [ebp + nb103_faction]	
	mov   eax, [ebp + nb103_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb103_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb103_ninner]
	mov   [esp + nb103_ninner], ecx
	add   edx, 0
	mov   [esp + nb103_innerk], edx	;# number of innerloop atoms 
	jge   .nb103_unroll_loop
	jmp   .nb103_odd_inner
.nb103_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb103_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb103_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb103_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb103_iqH]
	mulps  xmm4, [esp + nb103_iqM]

	movaps  [esp + nb103_qqH], xmm3
	movaps  [esp + nb103_qqM], xmm4	

	mov esi, [ebp + nb103_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb103_ixH1]
	movaps xmm5, [esp + nb103_iyH1]
	movaps xmm6, [esp + nb103_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb103_dxH1], xmm4
	movaps [esp + nb103_dyH1], xmm5
	movaps [esp + nb103_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [esp + nb103_ixH2]
	movaps xmm5, [esp + nb103_iyH2]
	movaps xmm6, [esp + nb103_izH2]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb103_dxH2], xmm4
	movaps [esp + nb103_dyH2], xmm5
	movaps [esp + nb103_dzH2], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [esp + nb103_ixM]
	movaps xmm4, [esp + nb103_iyM]
	movaps xmm5, [esp + nb103_izM]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [esp + nb103_dxM], xmm3
	movaps [esp + nb103_dyM], xmm4
	movaps [esp + nb103_dzM], xmm5
	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	mulps  xmm7, [esp + nb103_qqH]	;# xmm7=vcoul 
	
	mulps  xmm4, xmm7	;# total fsH1 in xmm4 

	addps  xmm7, [esp + nb103_vctot]
	
	movaps [esp + nb103_vctot], xmm7

	movaps xmm0, [esp + nb103_dxH1]
	movaps xmm1, [esp + nb103_dyH1]
	movaps xmm2, [esp + nb103_dzH1]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H1 forces 
	movaps xmm3, [esp + nb103_fixH1]
	movaps xmm4, [esp + nb103_fiyH1]
	movaps xmm7, [esp + nb103_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb103_fixH1], xmm3
	movaps [esp + nb103_fiyH1], xmm4
	movaps [esp + nb103_fizH1], xmm7
	;# update j forces with water O 
	movaps [esp + nb103_fjx], xmm0
	movaps [esp + nb103_fjy], xmm1
	movaps [esp + nb103_fjz], xmm2

	;# H2 interactions 
	movaps  xmm4, xmm6	
	mulps   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulps  xmm6, [esp + nb103_qqH]	;# xmm6=vcoul 
	mulps  xmm4, xmm6		;# total fsH2 in xmm4 
	
	addps  xmm6, [esp + nb103_vctot]

	movaps xmm0, [esp + nb103_dxH2]
	movaps xmm1, [esp + nb103_dyH2]
	movaps xmm2, [esp + nb103_dzH2]
	movaps [esp + nb103_vctot], xmm6
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update H2 forces 
	movaps xmm3, [esp + nb103_fixH2]
	movaps xmm4, [esp + nb103_fiyH2]
	movaps xmm7, [esp + nb103_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb103_fixH2], xmm3
	movaps [esp + nb103_fiyH2], xmm4
	movaps [esp + nb103_fizH2], xmm7
	;# update j forces with water H2 
	addps  xmm0, [esp + nb103_fjx]
	addps  xmm1, [esp + nb103_fjy]
	addps  xmm2, [esp + nb103_fjz]
	movaps [esp + nb103_fjx], xmm0
	movaps [esp + nb103_fjy], xmm1
	movaps [esp + nb103_fjz], xmm2

	;# M interactions 
	movaps  xmm4, xmm5	
	mulps   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulps  xmm5, [esp + nb103_qqM]	;# xmm5=vcoul 
	mulps  xmm4, xmm5		;# total fsM in xmm4 
	
	addps  xmm5, [esp + nb103_vctot]

	movaps xmm0, [esp + nb103_dxM]
	movaps xmm1, [esp + nb103_dyM]
	movaps xmm2, [esp + nb103_dzM]
	movaps [esp + nb103_vctot], xmm5
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4

	;# update M forces 
	movaps xmm3, [esp + nb103_fixM]
	movaps xmm4, [esp + nb103_fiyM]
	movaps xmm7, [esp + nb103_fizM]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [esp + nb103_fixM], xmm3
	movaps [esp + nb103_fiyM], xmm4
	movaps [esp + nb103_fizM], xmm7

	mov edi, [ebp + nb103_faction]
	;# update j forces 
	addps xmm0, [esp + nb103_fjx]
	addps xmm1, [esp + nb103_fjy]
	addps xmm2, [esp + nb103_fjz]

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
	sub dword ptr [esp + nb103_innerk],  4
	jl    .nb103_odd_inner
	jmp   .nb103_unroll_loop
.nb103_odd_inner:	
	add dword ptr [esp + nb103_innerk],  4
	jnz   .nb103_odd_loop
	jmp   .nb103_updateouterdata
.nb103_odd_loop:
	mov   edx, [esp + nb103_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb103_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb103_iqM]
	mov esi, [ebp + nb103_charge] 
	movhps xmm4, [esp + nb103_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb103_qqM], xmm3	;# use dummy qq for storage 

	mov esi, [ebp + nb103_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb103_ixM]
	movss xmm4, [esp + nb103_iyM]
	movss xmm5, [esp + nb103_izM]
		
	movlps xmm6, [esp + nb103_ixH1]
	movlps xmm7, [esp + nb103_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb103_iyH1]
	movlps xmm7, [esp + nb103_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb103_izH1]
	movlps xmm7, [esp + nb103_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# use dummy dx for storage
	movaps [esp + nb103_dxM], xmm3
	movaps [esp + nb103_dyM], xmm4
	movaps [esp + nb103_dzM], xmm5

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
	movaps xmm1, [esp + nb103_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb103_half]
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
	movaps xmm3, [esp + nb103_qqM]

	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm4, xmm3	;# xmm4=total fscal 
	addps  xmm3, [esp + nb103_vctot]

	movaps xmm0, [esp + nb103_dxM]
	movaps xmm1, [esp + nb103_dyM]
	movaps xmm2, [esp + nb103_dzM]

	movaps [esp + nb103_vctot], xmm3

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [esp + nb103_fixM]	
	movss  xmm4, [esp + nb103_fiyM]	
	movss  xmm5, [esp + nb103_fizM]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esp + nb103_fixM], xmm3	
	movss  [esp + nb103_fiyM], xmm4	
	movss  [esp + nb103_fizM], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [esp + nb103_fixH1]
	addss  xmm4, [esp + nb103_fiyH1]
	addss  xmm5, [esp + nb103_fizH1]
	movss  [esp + nb103_fixH1], xmm3	
	movss  [esp + nb103_fiyH1], xmm4	
	movss  [esp + nb103_fizH1], xmm5	;# updated the H1 force 

	mov edi, [ebp + nb103_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [esp + nb103_fixH2]
	addss  xmm4, [esp + nb103_fiyH2]
	addss  xmm5, [esp + nb103_fizH2]
	movss  [esp + nb103_fixH2], xmm3	
	movss  [esp + nb103_fiyH2], xmm4	
	movss  [esp + nb103_fizH2], xmm5	;# updated the H2 force 

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
	addss    xmm2, xmm1	;# z sum in xmm2 
	subps    xmm6, xmm0
	subss    xmm7, xmm2
	
	movlps [edi + eax*4], xmm6
	movss  [edi + eax*4 + 8], xmm7

	dec   dword ptr [esp + nb103_innerk]
	jz    .nb103_updateouterdata
	jmp   .nb103_odd_loop
.nb103_updateouterdata:
	mov   ecx, [esp + nb103_ii3]
	mov   edi, [ebp + nb103_faction]
	mov   esi, [ebp + nb103_fshift]
	mov   edx, [esp + nb103_is3]

	;# accumulate H1 forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb103_fixH1]
	movaps xmm1, [esp + nb103_fiyH1]
	movaps xmm2, [esp + nb103_fizH1]

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
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H2 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb103_fixH2]
	movaps xmm1, [esp + nb103_fiyH2]
	movaps xmm2, [esp + nb103_fizH2]

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

	;# accumulate M i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb103_fixM]
	movaps xmm1, [esp + nb103_fiyM]
	movaps xmm2, [esp + nb103_fizM]

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
	movss  xmm3, [edi + ecx*4 + 36]
	movss  xmm4, [edi + ecx*4 + 40]
	movss  xmm5, [edi + ecx*4 + 44]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4 + 36], xmm3
	movss  [edi + ecx*4 + 40], xmm4
	movss  [edi + ecx*4 + 44], xmm5

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
	mov esi, [esp + nb103_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb103_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb103_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb103_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb103_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb103_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb103_n], esi
        jmp .nb103_outer
.nb103_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb103_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103_end
        ;# non-zero, do one more workunit
        jmp   .nb103_threadloop
.nb103_end:
	emms

	mov eax, [esp + nb103_nouter]
	mov ebx, [esp + nb103_ninner]
	mov ecx, [ebp + nb103_outeriter]
	mov edx, [ebp + nb103_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb103_salign]
	add esp, eax
	add esp, 632
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret



.globl nb_kernel103nf_ia32_sse
.globl _nb_kernel103nf_ia32_sse
nb_kernel103nf_ia32_sse:	
_nb_kernel103nf_ia32_sse:	
.equiv          nb103nf_p_nri,          8
.equiv          nb103nf_iinr,           12
.equiv          nb103nf_jindex,         16
.equiv          nb103nf_jjnr,           20
.equiv          nb103nf_shift,          24
.equiv          nb103nf_shiftvec,       28
.equiv          nb103nf_fshift,         32
.equiv          nb103nf_gid,            36
.equiv          nb103nf_pos,            40
.equiv          nb103nf_faction,        44
.equiv          nb103nf_charge,         48
.equiv          nb103nf_p_facel,        52
.equiv          nb103nf_p_krf,          56
.equiv          nb103nf_p_crf,          60
.equiv          nb103nf_Vc,             64
.equiv          nb103nf_type,           68
.equiv          nb103nf_p_ntype,        72
.equiv          nb103nf_vdwparam,       76
.equiv          nb103nf_Vvdw,           80
.equiv          nb103nf_p_tabscale,     84
.equiv          nb103nf_VFtab,          88
.equiv          nb103nf_invsqrta,       92
.equiv          nb103nf_dvda,           96
.equiv          nb103nf_p_gbtabscale,   100
.equiv          nb103nf_GBtab,          104
.equiv          nb103nf_p_nthreads,     108
.equiv          nb103nf_count,          112
.equiv          nb103nf_mtx,            116
.equiv          nb103nf_outeriter,      120
.equiv          nb103nf_inneriter,      124
.equiv          nb103nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb103nf_ixH1,           0
.equiv          nb103nf_iyH1,           16
.equiv          nb103nf_izH1,           32
.equiv          nb103nf_ixH2,           48
.equiv          nb103nf_iyH2,           64
.equiv          nb103nf_izH2,           80
.equiv          nb103nf_ixM,            96
.equiv          nb103nf_iyM,            112
.equiv          nb103nf_izM,            128
.equiv          nb103nf_iqH,            144
.equiv          nb103nf_iqM,            160
.equiv          nb103nf_vctot,          176
.equiv          nb103nf_half,           192
.equiv          nb103nf_three,          208
.equiv          nb103nf_qqH,            224
.equiv          nb103nf_qqM,            240
.equiv          nb103nf_is3,            256
.equiv          nb103nf_ii3,            260
.equiv          nb103nf_innerjjnr,      264
.equiv          nb103nf_innerk,         268
.equiv          nb103nf_n,              272
.equiv          nb103nf_nn1,            276
.equiv          nb103nf_nri,            280
.equiv          nb103nf_nouter,         284
.equiv          nb103nf_ninner,         288
.equiv          nb103nf_salign,         292
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 296		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb103nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb103nf_p_nri]
	mov ecx, [ecx]
	mov [esp + nb103nf_nri], ecx

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb103nf_nouter], eax
	mov [esp + nb103nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb103nf_half], eax
	movss xmm1, [esp + nb103nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb103nf_half],  xmm1
	movaps [esp + nb103nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   ecx, [ebp + nb103nf_iinr]   	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]		;# ebx =ii 

	mov   edx, [ebp + nb103nf_charge]
	movss xmm3, [edx + ebx*4 + 4]	
	movss xmm4, [edx + ebx*4 + 12]	
	mov esi, [ebp + nb103nf_p_facel]
	movss xmm5, [esi]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [esp + nb103nf_iqH], xmm3
	movaps [esp + nb103nf_iqM], xmm4
	
.nb103nf_threadloop:
        mov   esi, [ebp + nb103nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb103nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                            ;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb103nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb103nf_n], eax
        mov [esp + nb103nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb103nf_outerstart
        jmp .nb103nf_end
	
.nb103nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb103nf_nouter]
	mov [esp + nb103nf_nouter], ebx

.nb103nf_outer:
	mov   eax, [ebp + nb103nf_shift]  	;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]	;# ebx=3*is 
	mov   [esp + nb103nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb103nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb103nf_iinr]   	;# ecx = pointer into iinr[] 
	mov   ebx, [ecx + esi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb103nf_pos]	;# eax = base of pos[]  
	mov   [esp + nb103nf_ii3], ebx

	addss xmm3, [eax + ebx*4 + 12]
	addss xmm4, [eax + ebx*4 + 16]
	addss xmm5, [eax + ebx*4 + 20]	
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb103nf_ixH1], xmm3
	movaps [esp + nb103nf_iyH1], xmm4
	movaps [esp + nb103nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [eax + ebx*4 + 24]
	addss xmm1, [eax + ebx*4 + 28]
	addss xmm2, [eax + ebx*4 + 32]		
	addss xmm3, [eax + ebx*4 + 36]
	addss xmm4, [eax + ebx*4 + 40]
	addss xmm5, [eax + ebx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [esp + nb103nf_ixH2], xmm0
	movaps [esp + nb103nf_iyH2], xmm1
	movaps [esp + nb103nf_izH2], xmm2
	movaps [esp + nb103nf_ixM], xmm3
	movaps [esp + nb103nf_iyM], xmm4
	movaps [esp + nb103nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [esp + nb103nf_vctot], xmm4

	mov   eax, [ebp + nb103nf_jindex]
	mov   ecx, [eax + esi*4]	 	;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   esi, [ebp + nb103nf_pos]
	mov   edi, [ebp + nb103nf_faction]	
	mov   eax, [ebp + nb103nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb103nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb103nf_ninner]
	mov   [esp + nb103nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb103nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb103nf_unroll_loop
	jmp   .nb103nf_odd_inner
.nb103nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb103nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]     	;# eax-edx=jnr1-4 

	add dword ptr [esp + nb103nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb103nf_charge]	;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [esp + nb103nf_iqH]
	mulps  xmm4, [esp + nb103nf_iqM]

	movaps  [esp + nb103nf_qqH], xmm3
	movaps  [esp + nb103nf_qqM], xmm4	

	mov esi, [ebp + nb103nf_pos]   	;# base of pos[] 

	lea   eax, [eax + eax*2] 	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2] 	;# replace jnr with j3 
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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [esp + nb103nf_ixH1]
	movaps xmm5, [esp + nb103nf_iyH1]
	movaps xmm6, [esp + nb103nf_izH1]

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

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [esp + nb103nf_ixH2]
	movaps xmm5, [esp + nb103nf_iyH2]
	movaps xmm6, [esp + nb103nf_izH2]

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

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [esp + nb103nf_ixM]
	movaps xmm4, [esp + nb103nf_iyM]
	movaps xmm5, [esp + nb103nf_izM]

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
	;# rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103nf_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103nf_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [esp + nb103nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [esp + nb103nf_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions - xmm7=rinv
	mulps  xmm7, [esp + nb103nf_qqH]	;# xmm7=vcoul 
	addps  xmm7, [esp + nb103nf_vctot]
	movaps [esp + nb103nf_vctot], xmm7

	;# H2 interactions - xmm6=rinv
	mulps  xmm6, [esp + nb103nf_qqH]	;# xmm6=vcoul 
	addps  xmm6, xmm7
	movaps [esp + nb103nf_vctot], xmm6

	;# M interactions  - xmm5=rinv
	mulps  xmm5, [esp + nb103nf_qqM]	;# xmm5=vcoul 
	addps  xmm5, xmm6
	movaps [esp + nb103nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb103nf_innerk],  4
	jl    .nb103nf_odd_inner
	jmp   .nb103nf_unroll_loop
.nb103nf_odd_inner:	
	add dword ptr [esp + nb103nf_innerk],  4
	jnz   .nb103nf_odd_loop
	jmp   .nb103nf_updateouterdata
.nb103nf_odd_loop:
	mov   edx, [esp + nb103nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	add dword ptr [esp + nb103nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [esp + nb103nf_iqM]
	mov esi, [ebp + nb103nf_charge] 
	movhps xmm4, [esp + nb103nf_iqH]     
	movss xmm3, [esi + eax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [esp + nb103nf_qqM], xmm3	;# use dummy qq for storage 

	mov esi, [ebp + nb103nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [esi + eax*4]
	movss xmm1, [esi + eax*4 + 4]
	movss xmm2, [esi + eax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [esp + nb103nf_ixM]
	movss xmm4, [esp + nb103nf_iyM]
	movss xmm5, [esp + nb103nf_izM]
		
	movlps xmm6, [esp + nb103nf_ixH1]
	movlps xmm7, [esp + nb103nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [esp + nb103nf_iyH1]
	movlps xmm7, [esp + nb103nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [esp + nb103nf_izH1]
	movlps xmm7, [esp + nb103nf_izH2]
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
	movaps xmm1, [esp + nb103nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb103nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000	- xmm0=rinv
	movaps xmm3, [esp + nb103nf_qqM]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm3, [esp + nb103nf_vctot]
	movaps [esp + nb103nf_vctot], xmm3

	dec   dword ptr [esp + nb103nf_innerk]
	jz    .nb103nf_updateouterdata
	jmp   .nb103nf_odd_loop
.nb103nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb103nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb103nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb103nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   eax, [ebp + nb103nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [esp + nb103nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb103nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb103nf_n], esi
        jmp .nb103nf_outer
.nb103nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb103nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103nf_end
        ;# non-zero, do one more workunit
        jmp   .nb103nf_threadloop
.nb103nf_end:
	emms

	mov eax, [esp + nb103nf_nouter]
	mov ebx, [esp + nb103nf_ninner]
	mov ecx, [ebp + nb103nf_outeriter]
	mov edx, [ebp + nb103nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb103nf_salign]
	add esp, eax
	add esp, 296
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret

