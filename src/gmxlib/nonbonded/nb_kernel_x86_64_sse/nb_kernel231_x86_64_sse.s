;#
;# $Id$
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
.equiv          .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as




.globl nb_kernel231_x86_64_sse
.globl _nb_kernel231_x86_64_sse
nb_kernel231_x86_64_sse:	
_nb_kernel231_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb231_fshift,           16
.equiv          nb231_gid,              24
.equiv          nb231_pos,              32
.equiv          nb231_faction,          40
.equiv          nb231_charge,           48
.equiv          nb231_p_facel,          56
.equiv          nb231_argkrf,           64
.equiv          nb231_argcrf,           72
.equiv          nb231_Vc,               80
.equiv          nb231_type,             88
.equiv          nb231_p_ntype,          96
.equiv          nb231_vdwparam,         104
.equiv          nb231_Vvdw,             112
.equiv          nb231_p_tabscale,       120
.equiv          nb231_VFtab,            128
.equiv          nb231_invsqrta,         136
.equiv          nb231_dvda,             144
.equiv          nb231_p_gbtabscale,     152
.equiv          nb231_GBtab,            160
.equiv          nb231_p_nthreads,       168
.equiv          nb231_count,            176
.equiv          nb231_mtx,              184
.equiv          nb231_outeriter,        192
.equiv          nb231_inneriter,        200
.equiv          nb231_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb231_ixO,              0
.equiv          nb231_iyO,              16
.equiv          nb231_izO,              32
.equiv          nb231_ixH1,             48
.equiv          nb231_iyH1,             64
.equiv          nb231_izH1,             80
.equiv          nb231_ixH2,             96
.equiv          nb231_iyH2,             112
.equiv          nb231_izH2,             128
.equiv          nb231_iqO,              144
.equiv          nb231_iqH,              160
.equiv          nb231_dxO,              176
.equiv          nb231_dyO,              192
.equiv          nb231_dzO,              208
.equiv          nb231_dxH1,             224
.equiv          nb231_dyH1,             240
.equiv          nb231_dzH1,             256
.equiv          nb231_dxH2,             272
.equiv          nb231_dyH2,             288
.equiv          nb231_dzH2,             304
.equiv          nb231_qqO,              320
.equiv          nb231_qqH,              336
.equiv          nb231_c6,               352
.equiv          nb231_c12,              368
.equiv          nb231_tsc,              384
.equiv          nb231_fstmp ,           400
.equiv          nb231_vctot,            416
.equiv          nb231_Vvdwtot,          432
.equiv          nb231_fixO,             448
.equiv          nb231_fiyO,             464
.equiv          nb231_fizO,             480
.equiv          nb231_fixH1,            496
.equiv          nb231_fiyH1,            512
.equiv          nb231_fizH1,            528
.equiv          nb231_fixH2,            544
.equiv          nb231_fiyH2,            560
.equiv          nb231_fizH2,            576
.equiv          nb231_fjx,              592
.equiv          nb231_fjy,              608
.equiv          nb231_fjz,              624
.equiv          nb231_half,             640
.equiv          nb231_three,            656
.equiv          nb231_two,              672
.equiv          nb231_krf,              688
.equiv          nb231_crf,              704
.equiv          nb231_krsqO,            720
.equiv          nb231_krsqH1,           736
.equiv          nb231_krsqH2,           752
.equiv          nb231_rinvO,            768
.equiv          nb231_rinvH1,           784
.equiv          nb231_rinvH2,           800
.equiv          nb231_facel,            816
.equiv          nb231_iinr,             824
.equiv          nb231_jindex,           832
.equiv          nb231_jjnr,             840
.equiv          nb231_shift,            848
.equiv          nb231_shiftvec,         856
.equiv          nb231_innerjjnr,        864
.equiv          nb231_nri,              872
.equiv          nb231_is3,              876
.equiv          nb231_ii3,              880
.equiv          nb231_ntia,             884
.equiv          nb231_innerk,           888
.equiv          nb231_n,                892
.equiv          nb231_nn1,              896
.equiv          nb231_nouter,           900
.equiv          nb231_ninner,           904

        push rbp
        mov  rbp, rsp
        push rbx

        femms
        sub rsp, 920	; # local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb231_nouter], eax
	mov [rsp + nb231_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb231_nri], edi
	mov [rsp + nb231_iinr], rsi
	mov [rsp + nb231_jindex], rdx
	mov [rsp + nb231_jjnr], rcx
	mov [rsp + nb231_shift], r8
	mov [rsp + nb231_shiftvec], r9
	mov rsi, [rbp + nb231_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb231_facel], xmm0

	mov rax, [rbp + nb231_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb231_tsc], xmm3

	mov rsi, [rbp + nb231_argkrf]
	mov rdi, [rbp + nb231_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb231_krf], xmm1
	movaps [rsp + nb231_crf], xmm2


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb231_half], eax
	movss xmm1, [rsp + nb231_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb231_half],  xmm1
	movaps [rsp + nb231_two],  xmm2
	movaps [rsp + nb231_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb231_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb231_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb231_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb231_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb231_iqO], xmm3
	movaps [rsp + nb231_iqH], xmm4
	
	mov   rdx, [rbp + nb231_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb231_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb231_ntia], ecx		

.nb231_threadloop:
        mov   rsi, [rbp + nb231_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb231_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb231_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb231_n], eax
        mov [rsp + nb231_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231_outerstart
        jmp .nb231_end

.nb231_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb231_nouter]
	mov [rsp + nb231_nouter], ebx

.nb231_outer:
	mov   rax, [rsp + nb231_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax +rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb231_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb231_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb231_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb231_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb231_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb231_ixO], xmm3
	movaps [rsp + nb231_iyO], xmm4
	movaps [rsp + nb231_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb231_ixH1], xmm0
	movaps [rsp + nb231_iyH1], xmm1
	movaps [rsp + nb231_izH1], xmm2
	movaps [rsp + nb231_ixH2], xmm3
	movaps [rsp + nb231_iyH2], xmm4
	movaps [rsp + nb231_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb231_vctot], xmm4
	movaps [rsp + nb231_Vvdwtot], xmm4
	movaps [rsp + nb231_fixO], xmm4
	movaps [rsp + nb231_fiyO], xmm4
	movaps [rsp + nb231_fizO], xmm4
	movaps [rsp + nb231_fixH1], xmm4
	movaps [rsp + nb231_fiyH1], xmm4
	movaps [rsp + nb231_fizH1], xmm4
	movaps [rsp + nb231_fixH2], xmm4
	movaps [rsp + nb231_fiyH2], xmm4
	movaps [rsp + nb231_fizH2], xmm4
	
	mov   rax, [rsp + nb231_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb231_pos]
	mov   rdi, [rbp + nb231_faction]	
	mov   rax, [rsp + nb231_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb231_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb231_ninner]
	mov   [rsp + nb231_ninner], ecx
	add   edx, 0
	mov   [rsp + nb231_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231_unroll_loop
	jmp   .nb231_odd_inner
.nb231_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb231_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb231_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb231_iqO]
	mulps  xmm4, [rsp + nb231_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb231_qqO], xmm3
	movaps  [rsp + nb231_qqH], xmm4
	
	mov rsi, [rbp + nb231_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb231_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb231_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb231_c6], xmm4
	movaps [rsp + nb231_c12], xmm6

	mov rsi, [rbp + nb231_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb231_ixO]
	movaps xmm5, [rsp + nb231_iyO]
	movaps xmm6, [rsp + nb231_izO]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb231_dxO], xmm4
	movaps [rsp + nb231_dyO], xmm5
	movaps [rsp + nb231_dzO], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	movaps xmm7, xmm4
	;# rsqO in xmm7 
	
	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb231_ixH1]
	movaps xmm5, [rsp + nb231_iyH1]
	movaps xmm6, [rsp + nb231_izH1]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb231_dxH1], xmm4
	movaps [rsp + nb231_dyH1], xmm5
	movaps [rsp + nb231_dzH1], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb231_ixH2]
	movaps xmm4, [rsp + nb231_iyH2]
	movaps xmm5, [rsp + nb231_izH2]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# store dr 
	movaps [rsp + nb231_dxH2], xmm3
	movaps [rsp + nb231_dyH2], xmm4
	movaps [rsp + nb231_dzH2], xmm5
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
	
	mulps  xmm0, [rsp + nb231_krf]	
	mulps  xmm1, [rsp + nb231_krf]	
	mulps  xmm2, [rsp + nb231_krf]	

	movaps [rsp + nb231_krsqH2], xmm0
	movaps [rsp + nb231_krsqH1], xmm1
	movaps [rsp + nb231_krsqO], xmm2
	
	;# start with rsqO (still in xmm7) - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231_half]
	movaps  [rsp + nb231_rinvO], xmm4	;# rinvO in xmm0, rsqO in xmm7
	
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231_half]
	movaps  [rsp + nb231_rinvH1], xmm4	;# rinvH1 in xmm6 

	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231_half]
	movaps  [rsp + nb231_rinvH2], xmm4	;# rinvH2 in xmm5 
	
	
	;# do O table interactions - rsqO in xmm7.
	mulps xmm7, [rsp + nb231_rinvO]
	mulps xmm7, [rsp + nb231_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx	

	mov  rsi, [rbp + nb231_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb231_c6]
	mulps  xmm7, xmm4	 ;# fijD 
	mulps  xmm5, xmm4	 ;# Vvdw6 
	mulps  xmm7, [rsp + nb231_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addps  xmm5, [rsp + nb231_Vvdwtot]
	movaps [rsp + nb231_fstmp], xmm7
	movaps [rsp + nb231_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb231_c12]
	mulps  xmm7, xmm4 ;# fijR 
	mulps  xmm5, xmm4 ;# Vvdw12 
	mulps  xmm7, [rsp + nb231_tsc]
	addps  xmm7, [rsp + nb231_fstmp]
	movaps [rsp + nb231_fstmp], xmm7

	addps  xmm5, [rsp + nb231_Vvdwtot]
	movaps [rsp + nb231_Vvdwtot], xmm5
	
	movaps xmm2, [rsp + nb231_rinvO]
	movaps xmm1, [rsp + nb231_krsqO]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm2, [rsp + nb231_qqO]
	mulps  xmm3, [rsp + nb231_qqO]

	mulps  xmm3, xmm0
	subps  xmm3, [rsp + nb231_fstmp]
	mulps  xmm3, xmm0

	addps  xmm2, [rsp + nb231_vctot]
	movaps [rsp + nb231_vctot], xmm2

	movaps xmm0, [rsp + nb231_dxO]
	movaps xmm1, [rsp + nb231_dyO]
	movaps xmm2, [rsp + nb231_dzO]
	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3	

	;# update O forces 
	movaps xmm3, [rsp + nb231_fixO]
	movaps xmm4, [rsp + nb231_fiyO]
	movaps xmm7, [rsp + nb231_fizO]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [rsp + nb231_fixO], xmm3
	movaps [rsp + nb231_fiyO], xmm4
	movaps [rsp + nb231_fizO], xmm7
	;# update j forces with water O 
	movaps [rsp + nb231_fjx], xmm0
	movaps [rsp + nb231_fjy], xmm1
	movaps [rsp + nb231_fjz], xmm2

	;# H1 interactions 
	movaps xmm2, [rsp + nb231_rinvH1]
	movaps xmm1, [rsp + nb231_krsqH1]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm0, xmm0
	mulps  xmm2, [rsp + nb231_qqH]
	mulps  xmm3, [rsp + nb231_qqH]	
	
	mulps  xmm3, xmm0
	addps  xmm2, [rsp + nb231_vctot]
	movaps [rsp + nb231_vctot], xmm2

	movaps xmm0, [rsp + nb231_dxH1]
	movaps xmm1, [rsp + nb231_dyH1]
	movaps xmm2, [rsp + nb231_dzH1]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3

	;# update H1 forces 
	movaps xmm3, [rsp + nb231_fixH1]
	movaps xmm4, [rsp + nb231_fiyH1]
	movaps xmm7, [rsp + nb231_fizH1]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [rsp + nb231_fixH1], xmm3
	movaps [rsp + nb231_fiyH1], xmm4
	movaps [rsp + nb231_fizH1], xmm7
	;# update j forces with water H1 
	addps  xmm0, [rsp + nb231_fjx]
	addps  xmm1, [rsp + nb231_fjy]
	addps  xmm2, [rsp + nb231_fjz]
	movaps [rsp + nb231_fjx], xmm0
	movaps [rsp + nb231_fjy], xmm1
	movaps [rsp + nb231_fjz], xmm2

	;# H2 interactions 
	movaps xmm2, [rsp + nb231_rinvH2]
	movaps xmm1, [rsp + nb231_krsqH2]
	movaps xmm3, xmm2
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231_crf] ;# rinv+krsq-crf
	movaps xmm0, xmm3
	subps  xmm3, xmm1
	subps  xmm3, xmm1 ;# rinv-2*krsq
	mulps  xmm0, xmm0
	mulps  xmm2, [rsp + nb231_qqH]
	mulps  xmm3, [rsp + nb231_qqH]	
	
	mulps  xmm3, xmm0
	addps  xmm2, [rsp + nb231_vctot]
	movaps [rsp + nb231_vctot], xmm2

	movaps xmm0, [rsp + nb231_dxH2]
	movaps xmm1, [rsp + nb231_dyH2]
	movaps xmm2, [rsp + nb231_dzH2]

	mulps  xmm0, xmm3
	mulps  xmm1, xmm3
	mulps  xmm2, xmm3
	;# update H2 forces 
	movaps xmm3, [rsp + nb231_fixH2]
	movaps xmm4, [rsp + nb231_fiyH2]
	movaps xmm7, [rsp + nb231_fizH2]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm7, xmm2
	movaps [rsp + nb231_fixH2], xmm3
	movaps [rsp + nb231_fiyH2], xmm4
	movaps [rsp + nb231_fizH2], xmm7

	mov rdi, [rbp + nb231_faction]
	;# update j forces 
	addps xmm0, [rsp + nb231_fjx]
	addps xmm1, [rsp + nb231_fjy]
	addps xmm2, [rsp + nb231_fjz]

	movlps xmm4, [rdi + rax*4]
	movlps xmm7, [rdi + rcx*4]
	movhps xmm4, [rdi + rbx*4]
	movhps xmm7, [rdi + rdx*4]
	
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
	movlps [rdi + rax*4], xmm7
	movlps [rdi + rcx*4], xmm3
	movhps [rdi + rbx*4], xmm7
	movhps [rdi + rdx*4], xmm3
	;# finally z forces 
	movss  xmm0, [rdi + rax*4 + 8]
	movss  xmm1, [rdi + rbx*4 + 8]
	movss  xmm3, [rdi + rcx*4 + 8]
	movss  xmm4, [rdi + rdx*4 + 8]
	subss  xmm0, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm1, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm3, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm4, xmm2
	movss  [rdi + rax*4 + 8], xmm0
	movss  [rdi + rbx*4 + 8], xmm1
	movss  [rdi + rcx*4 + 8], xmm3
	movss  [rdi + rdx*4 + 8], xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb231_innerk],  4
	jl    .nb231_odd_inner
	jmp   .nb231_unroll_loop
.nb231_odd_inner:	
	add dword ptr [rsp + nb231_innerk],  4
	jnz   .nb231_odd_loop
	jmp   .nb231_updateouterdata
.nb231_odd_loop:
	mov   rdx, [rsp + nb231_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb231_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb231_iqO]
	mov rsi, [rbp + nb231_charge] 
	movhps xmm4, [rsp + nb231_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb231_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb231_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb231_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb231_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb231_c6], xmm6
	movaps [rsp + nb231_c12], xmm7

	mov rsi, [rbp + nb231_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb231_ixO]
	movss xmm4, [rsp + nb231_iyO]
	movss xmm5, [rsp + nb231_izO]
		
	movlps xmm6, [rsp + nb231_ixH1]
	movlps xmm7, [rsp + nb231_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb231_iyH1]
	movlps xmm7, [rsp + nb231_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb231_izH1]
	movlps xmm7, [rsp + nb231_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb231_dxO], xmm3
	movaps [rsp + nb231_dyO], xmm4
	movaps [rsp + nb231_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb231_krf]
	movaps [rsp + nb231_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb231_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb231_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [rsp + nb231_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  rsi, [rbp + nb231_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb231_c6]
	mulss  xmm7, xmm4	 ;# fijD 
	mulss  xmm5, xmm4	 ;# Vvdw6 
	mulss  xmm7, [rsp + nb231_tsc]
	;# put scalar force on stack Update Vvdwtot directly 
	addss  xmm5, [rsp + nb231_Vvdwtot]
	movss [rsp + nb231_fstmp], xmm7
	movss [rsp + nb231_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb231_two]	;# two*Heps2 
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb231_c12]
	mulss  xmm7, xmm4 ;# fijR 
	mulss  xmm5, xmm4 ;# Vvdw12 
	mulss  xmm7, [rsp + nb231_tsc]
	addss  xmm7, [rsp + nb231_fstmp]
	movss [rsp + nb231_fstmp], xmm7

	addss  xmm5, [rsp + nb231_Vvdwtot]
	movss [rsp + nb231_Vvdwtot], xmm5

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm4, xmm0
	movaps xmm3, [rsp + nb231_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	mulps  xmm3, [rsp + nb231_two]
	subps  xmm0, [rsp + nb231_crf] ;# xmm0=rinv+ krsq-crf 
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [rsp + nb231_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [rsp + nb231_qqO] 	;# xmm1=coul part of fs 

	mulps  xmm1, xmm4
	subss  xmm1, [rsp + nb231_fstmp]
	mulps  xmm4, xmm1
	
	addps  xmm0, [rsp + nb231_vctot]
	movaps [rsp + nb231_vctot], xmm0
	
	movaps xmm0, [rsp + nb231_dxO]
	movaps xmm1, [rsp + nb231_dyO]
	movaps xmm2, [rsp + nb231_dzO]

	movd eax, mm0	

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb231_fixO]	
	movss  xmm4, [rsp + nb231_fiyO]	
	movss  xmm5, [rsp + nb231_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb231_fixO], xmm3	
	movss  [rsp + nb231_fiyO], xmm4	
	movss  [rsp + nb231_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# constant 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# constant 11100110
	shufps xmm5, xmm5, 230 ;# constant 11100110
	addss  xmm3, [rsp + nb231_fixH1]
	addss  xmm4, [rsp + nb231_fiyH1]
	addss  xmm5, [rsp + nb231_fizH1]
	movss  [rsp + nb231_fixH1], xmm3	
	movss  [rsp + nb231_fiyH1], xmm4	
	movss  [rsp + nb231_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb231_faction]
	shufps xmm3, xmm3, 231 ;# constant 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# constant 11100111
	shufps xmm5, xmm5, 231 ;# constant 11100111
	addss  xmm3, [rsp + nb231_fixH2]
	addss  xmm4, [rsp + nb231_fiyH2]
	addss  xmm5, [rsp + nb231_fizH2]
	movss  [rsp + nb231_fixH2], xmm3	
	movss  [rsp + nb231_fiyH2], xmm4	
	movss  [rsp + nb231_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
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
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb231_innerk]
	jz    .nb231_updateouterdata
	jmp   .nb231_odd_loop
.nb231_updateouterdata:
	mov   ecx, [rsp + nb231_ii3]
	mov   rdi, [rbp + nb231_faction]
	mov   rsi, [rbp + nb231_fshift]
	mov   edx, [rsp + nb231_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb231_fixO]
	movaps xmm1, [rsp + nb231_fiyO]
	movaps xmm2, [rsp + nb231_fizO]

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
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# constant 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb231_fixH1]
	movaps xmm1, [rsp + nb231_fiyH1]
	movaps xmm2, [rsp + nb231_fizH1]

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
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb231_fixH2]
	movaps xmm1, [rsp + nb231_fiyH2]
	movaps xmm2, [rsp + nb231_fizH2]

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
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# constant 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	addps  xmm3, xmm6
	addss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb231_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb231_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb231_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb231_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb231_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb231_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb231_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb231_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb231_n], esi
        jmp .nb231_outer
.nb231_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb231_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb231_end
        ;# non-zero, do one more workunit
        jmp   .nb231_threadloop
.nb231_end:
	mov eax, [rsp + nb231_nouter]
	mov ebx, [rsp + nb231_ninner]
	mov rcx, [rbp + nb231_outeriter]
	mov rdx, [rbp + nb231_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 920
	femms

	pop rbx
	pop	rbp
	ret




.globl nb_kernel231nf_x86_64_sse
.globl _nb_kernel231nf_x86_64_sse
nb_kernel231nf_x86_64_sse:	
_nb_kernel231nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb231nf_fshift,         16
.equiv          nb231nf_gid,            24
.equiv          nb231nf_pos,            32
.equiv          nb231nf_faction,        40
.equiv          nb231nf_charge,         48
.equiv          nb231nf_p_facel,        56
.equiv          nb231nf_argkrf,         64
.equiv          nb231nf_argcrf,         72
.equiv          nb231nf_Vc,             80
.equiv          nb231nf_type,           88
.equiv          nb231nf_p_ntype,        96
.equiv          nb231nf_vdwparam,       104
.equiv          nb231nf_Vvdw,           112
.equiv          nb231nf_p_tabscale,     120
.equiv          nb231nf_VFtab,          128
.equiv          nb231nf_invsqrta,       136
.equiv          nb231nf_dvda,           144
.equiv          nb231nf_p_gbtabscale,   152
.equiv          nb231nf_GBtab,          160
.equiv          nb231nf_p_nthreads,     168
.equiv          nb231nf_count,          176
.equiv          nb231nf_mtx,            184
.equiv          nb231nf_outeriter,      192
.equiv          nb231nf_inneriter,      200
.equiv          nb231nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb231nf_ixO,            0
.equiv          nb231nf_iyO,            16
.equiv          nb231nf_izO,            32
.equiv          nb231nf_ixH1,           48
.equiv          nb231nf_iyH1,           64
.equiv          nb231nf_izH1,           80
.equiv          nb231nf_ixH2,           96
.equiv          nb231nf_iyH2,           112
.equiv          nb231nf_izH2,           128
.equiv          nb231nf_iqO,            144
.equiv          nb231nf_iqH,            160
.equiv          nb231nf_qqO,            176
.equiv          nb231nf_qqH,            192
.equiv          nb231nf_c6,             208
.equiv          nb231nf_c12,            224
.equiv          nb231nf_vctot,          240
.equiv          nb231nf_Vvdwtot,        256
.equiv          nb231nf_half,           272
.equiv          nb231nf_three,          288
.equiv          nb231nf_krf,            304
.equiv          nb231nf_crf,            320
.equiv          nb231nf_rinvO,          336
.equiv          nb231nf_rinvH1,         352
.equiv          nb231nf_rinvH2,         368
.equiv          nb231nf_krsqO,          384
.equiv          nb231nf_krsqH1,         400
.equiv          nb231nf_krsqH2,         416
.equiv          nb231nf_tsc,            432
.equiv          nb231nf_facel,          448
.equiv          nb231nf_iinr,           456
.equiv          nb231nf_jindex,         464
.equiv          nb231nf_jjnr,           472
.equiv          nb231nf_shift,          480
.equiv          nb231nf_shiftvec,       488
.equiv          nb231nf_innerjjnr,      496
.equiv          nb231nf_nri,            504
.equiv          nb231nf_ntia,           508
.equiv          nb231nf_is3,            512
.equiv          nb231nf_ii3,            516
.equiv          nb231nf_innerk,         520
.equiv          nb231nf_n,              524
.equiv          nb231nf_nn1,            528
.equiv          nb231nf_nouter,         532
.equiv          nb231nf_ninner,         536
	
        push rbp
        mov  rbp, rsp
        push rbx
	
        femms
        sub rsp, 552         ; # local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb231nf_nouter], eax
	mov [rsp + nb231nf_ninner], eax
	
	mov edi, [rdi]
	mov [rsp + nb231nf_nri], edi
	mov [rsp + nb231nf_iinr], rsi
	mov [rsp + nb231nf_jindex], rdx
	mov [rsp + nb231nf_jjnr], rcx
	mov [rsp + nb231nf_shift], r8
	mov [rsp + nb231nf_shiftvec], r9
	mov rsi, [rbp + nb231nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb231nf_facel], xmm0

	mov rax, [rbp + nb231nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb231nf_tsc], xmm3

	mov rsi, [rbp + nb231nf_argkrf]
	mov rdi, [rbp + nb231nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb231nf_krf], xmm1
	movaps [rsp + nb231nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb231nf_half], eax
	movss xmm1, [rsp + nb231nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb231nf_half],  xmm1
	movaps [rsp + nb231nf_three],  xmm3	

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb231nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb231nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb231nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb231nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb231nf_iqO], xmm3
	movaps [rsp + nb231nf_iqH], xmm4
	
	mov   rdx, [rbp + nb231nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb231nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb231nf_ntia], ecx		

.nb231nf_threadloop:
        mov   rsi, [rbp + nb231nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb231nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb231nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb231nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb231nf_n], eax
        mov [rsp + nb231nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb231nf_outerstart
        jmp .nb231nf_end

.nb231nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb231nf_nouter]
	mov [rsp + nb231nf_nouter], ebx

.nb231nf_outer:
	mov   rax, [rsp + nb231nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax +rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb231nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb231nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb231nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb231nf_pos]    ;# eax = base of pos[]  
	mov   [rsp + nb231nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb231nf_ixO], xmm3
	movaps [rsp + nb231nf_iyO], xmm4
	movaps [rsp + nb231nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb231nf_ixH1], xmm0
	movaps [rsp + nb231nf_iyH1], xmm1
	movaps [rsp + nb231nf_izH1], xmm2
	movaps [rsp + nb231nf_ixH2], xmm3
	movaps [rsp + nb231nf_iyH2], xmm4
	movaps [rsp + nb231nf_izH2], xmm5
	
	;# clear vctot
	xorps xmm4, xmm4
	movaps [rsp + nb231nf_vctot], xmm4
	movaps [rsp + nb231nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb231nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb231nf_pos]
	mov   rax, [rsp + nb231nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb231nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb231nf_ninner]
	mov   [rsp + nb231nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb231nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb231nf_unroll_loop
	jmp   .nb231nf_odd_inner
.nb231nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb231nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb231nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb231nf_iqO]
	mulps  xmm4, [rsp + nb231nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb231nf_qqO], xmm3
	movaps  [rsp + nb231nf_qqH], xmm4
	
	mov rsi, [rbp + nb231nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb231nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb231nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb231nf_c6], xmm4
	movaps [rsp + nb231nf_c12], xmm6

	mov rsi, [rbp + nb231nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb231nf_ixO]
	movaps xmm5, [rsp + nb231nf_iyO]
	movaps xmm6, [rsp + nb231nf_izO]

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
	movaps xmm4, [rsp + nb231nf_ixH1]
	movaps xmm5, [rsp + nb231nf_iyH1]
	movaps xmm6, [rsp + nb231nf_izH1]

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
	movaps xmm3, [rsp + nb231nf_ixH2]
	movaps xmm4, [rsp + nb231nf_iyH2]
	movaps xmm5, [rsp + nb231nf_izH2]

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
	
	mulps  xmm0, [rsp + nb231nf_krf]	
	mulps  xmm1, [rsp + nb231nf_krf]	
	mulps  xmm2, [rsp + nb231nf_krf]	

	movaps [rsp + nb231nf_krsqH2], xmm0
	movaps [rsp + nb231nf_krsqH1], xmm1
	movaps [rsp + nb231nf_krsqO], xmm2
	
	;# start with rsqO (still in xmm7) - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231nf_half]
	movaps  [rsp + nb231nf_rinvO], xmm4	;# rinvO in xmm0, rsqO in xmm7
	
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231nf_half]
	movaps  [rsp + nb231nf_rinvH1], xmm4	;# rinvH1 in xmm6 

	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb231nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# constant 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb231nf_half]
	movaps  [rsp + nb231nf_rinvH2], xmm4	;# rinvH2 in xmm5 
	
	
	;# do O table interactions - rsqO in xmm7.
	mulps xmm7, [rsp + nb231nf_rinvO]
	mulps xmm7, [rsp + nb231nf_tsc] ;# rtab
	
	movhlps xmm5, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps  xmm7, xmm6	
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3

	mov  rsi, [rbp + nb231nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb231nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	addps  xmm5, [rsp + nb231nf_Vvdwtot]
	movaps [rsp + nb231nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb231nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 

	addps  xmm5, [rsp + nb231nf_Vvdwtot]
	movaps [rsp + nb231nf_Vvdwtot], xmm5
	
	movaps xmm2, [rsp + nb231nf_rinvO]
	movaps xmm1, [rsp + nb231nf_krsqO]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [rsp + nb231nf_qqO]

	addps  xmm2, [rsp + nb231nf_vctot]
	movaps [rsp + nb231nf_vctot], xmm2

	;# H1 interactions 
	movaps xmm2, [rsp + nb231nf_rinvH1]
	movaps xmm1, [rsp + nb231nf_krsqH1]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [rsp + nb231nf_qqH]
	
	addps  xmm2, [rsp + nb231nf_vctot]
	movaps [rsp + nb231nf_vctot], xmm2

	;# H2 interactions 
	movaps xmm2, [rsp + nb231nf_rinvH2]
	movaps xmm1, [rsp + nb231nf_krsqH2]
	addps  xmm2, xmm1 ;# rinv+krsq
	subps  xmm2, [rsp + nb231nf_crf] ;# rinv+krsq-crf
	mulps  xmm2, [rsp + nb231nf_qqH]
	addps  xmm2, [rsp + nb231nf_vctot]
	movaps [rsp + nb231nf_vctot], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb231nf_innerk],  4
	jl    .nb231nf_odd_inner
	jmp   .nb231nf_unroll_loop
.nb231nf_odd_inner:	
	add dword ptr [rsp + nb231nf_innerk],  4
	jnz   .nb231nf_odd_loop
	jmp   .nb231nf_updateouterdata
.nb231nf_odd_loop:
	mov   rdx, [rsp + nb231nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb231nf_innerjjnr],  4

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb231nf_iqO]
	mov rsi, [rbp + nb231nf_charge] 
	movhps xmm4, [rsp + nb231nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb231nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb231nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb231nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb231nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# constant 11111100
	shufps xmm7, xmm7, 253  ;# constant 11111101
	movaps [rsp + nb231nf_c6], xmm6
	movaps [rsp + nb231nf_c12], xmm7

	mov rsi, [rbp + nb231nf_pos]
	lea   eax, [eax + eax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb231nf_ixO]
	movss xmm4, [rsp + nb231nf_iyO]
	movss xmm5, [rsp + nb231nf_izO]
		
	movlps xmm6, [rsp + nb231nf_ixH1]
	movlps xmm7, [rsp + nb231nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb231nf_iyH1]
	movlps xmm7, [rsp + nb231nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb231nf_izH1]
	movlps xmm7, [rsp + nb231nf_izH2]
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
	mulps xmm0, [rsp + nb231nf_krf]
	movaps [rsp + nb231nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb231nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb231nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# constant 11100000
		
	;# rsq still in xmm4, rinv in xmm0.
	mulps  xmm4, xmm0
	mulps  xmm4, [rsp + nb231nf_tsc] ;# rtab
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss  xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  rsi, [rbp + nb231nf_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss xmm4, [rsp + nb231nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 

	addss  xmm5, [rsp + nb231nf_Vvdwtot]
	movss [rsp + nb231nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss xmm4, [rsp + nb231nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	addss  xmm5, [rsp + nb231nf_Vvdwtot]
	movss [rsp + nb231nf_Vvdwtot], xmm5

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm4, xmm0
	movaps xmm3, [rsp + nb231nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 

	subps  xmm0, [rsp + nb231nf_crf] ;# xmm0=rinv+ krsq-crf 

	mulps  xmm0, [rsp + nb231nf_qqO]	;# xmm0=vcoul 
	
	addps  xmm0, [rsp + nb231nf_vctot]
	movaps [rsp + nb231nf_vctot], xmm0
	
	dec dword ptr [rsp + nb231nf_innerk]
	jz    .nb231nf_updateouterdata
	jmp   .nb231nf_odd_loop
.nb231nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb231nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb231nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb231nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb231nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb231nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb231nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb231nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb231nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb231nf_n], esi
        jmp .nb231nf_outer
.nb231nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb231nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb231nf_end
        ;# non-zero, do one more workunit
        jmp   .nb231nf_threadloop
.nb231nf_end:
	mov eax, [rsp + nb231nf_nouter]
	mov ebx, [rsp + nb231nf_ninner]
	mov rcx, [rbp + nb231nf_outeriter]
	mov rdx, [rbp + nb231nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 552
	femms

	pop rbx
	pop	rbp
	ret


