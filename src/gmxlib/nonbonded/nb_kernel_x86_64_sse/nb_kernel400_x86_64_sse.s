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


	

.globl nb_kernel400_x86_64_sse
.globl _nb_kernel400_x86_64_sse
nb_kernel400_x86_64_sse:	
_nb_kernel400_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb400_fshift,           16
.equiv          nb400_gid,              24
.equiv          nb400_pos,              32
.equiv          nb400_faction,          40
.equiv          nb400_charge,           48
.equiv          nb400_p_facel,          56
.equiv          nb400_argkrf,           64
.equiv          nb400_argcrf,           72
.equiv          nb400_Vc,               80
.equiv          nb400_type,             88
.equiv          nb400_p_ntype,          96
.equiv          nb400_vdwparam,         104
.equiv          nb400_Vvdw,             112
.equiv          nb400_p_tabscale,       120
.equiv          nb400_VFtab,            128
.equiv          nb400_invsqrta,         136
.equiv          nb400_dvda,             144
.equiv          nb400_p_gbtabscale,     152
.equiv          nb400_GBtab,            160
.equiv          nb400_p_nthreads,       168
.equiv          nb400_count,            176
.equiv          nb400_mtx,              184
.equiv          nb400_outeriter,        192
.equiv          nb400_inneriter,        200
.equiv          nb400_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb400_ix,               0
.equiv          nb400_iy,               16
.equiv          nb400_iz,               32
.equiv          nb400_iq,               48
.equiv          nb400_dx,               64
.equiv          nb400_dy,               80
.equiv          nb400_dz,               96
.equiv          nb400_two,              112
.equiv          nb400_gbtsc,            128
.equiv          nb400_qq,               144
.equiv          nb400_r,                160
.equiv          nb400_vctot,            176
.equiv          nb400_fix,              192
.equiv          nb400_fiy,              208
.equiv          nb400_fiz,              224
.equiv          nb400_half,             240
.equiv          nb400_three,            256
.equiv          nb400_isai,             272
.equiv          nb400_isaprod,          288
.equiv          nb400_dvdasum,          304
.equiv          nb400_gbscale,          320
.equiv          nb400_nri,              336
.equiv          nb400_iinr,             344
.equiv          nb400_jindex,           352
.equiv          nb400_jjnr,             360
.equiv          nb400_shift,            368
.equiv          nb400_shiftvec,         376
.equiv          nb400_facel,            384
.equiv          nb400_innerjjnr,        392
.equiv          nb400_is3,              400
.equiv          nb400_ii3,              404
.equiv          nb400_ii,               408
.equiv          nb400_innerk,           412
.equiv          nb400_n,                416
.equiv          nb400_nn1,              420
.equiv          nb400_nouter,           424
.equiv          nb400_ninner,           428
.equiv          nb400_jnra,             432
.equiv          nb400_jnrb,             436
.equiv          nb400_jnrc,             440
.equiv          nb400_jnrd,             444

	push rbp
	mov  rbp, rsp
	push rbx

	
	femms
	sub rsp, 456		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb400_nouter], eax
	mov [rsp + nb400_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb400_nri], edi
	mov [rsp + nb400_iinr], rsi
	mov [rsp + nb400_jindex], rdx
	mov [rsp + nb400_jjnr], rcx
	mov [rsp + nb400_shift], r8
	mov [rsp + nb400_shiftvec], r9
	mov rsi, [rbp + nb400_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb400_facel], xmm0

	mov rbx, [rbp + nb400_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb400_gbtsc], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb400_half], eax
	movss xmm1, [rsp + nb400_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb400_half],  xmm1
	movaps [rsp + nb400_two],  xmm2
	movaps [rsp + nb400_three],  xmm3

.nb400_threadloop:
        mov   rsi, [rbp + nb400_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb400_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb400_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb400_n], eax
        mov [rsp + nb400_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400_outerstart
        jmp .nb400_end

.nb400_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb400_nouter]
	mov [rsp + nb400_nouter], ebx

.nb400_outer:
	mov   rax, [rsp + nb400_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb400_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb400_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb400_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 
	mov   [rsp + nb400_ii], ebx
	
	mov   rdx, [rbp + nb400_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb400_facel]
	shufps xmm3, xmm3, 0


	mov   rdx, [rbp + nb400_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb400_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb400_iq], xmm3
	movaps [rsp + nb400_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb400_ix], xmm0
	movaps [rsp + nb400_iy], xmm1
	movaps [rsp + nb400_iz], xmm2

	mov   [rsp + nb400_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb400_vctot], xmm4
	movaps [rsp + nb400_dvdasum], xmm4
	movaps [rsp + nb400_fix], xmm4
	movaps [rsp + nb400_fiy], xmm4
	movaps [rsp + nb400_fiz], xmm4
	
	mov   rax, [rsp + nb400_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb400_pos]
	mov   rdi, [rbp + nb400_faction]	
	mov   rax, [rsp + nb400_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb400_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb400_ninner]
	mov   [rsp + nb400_ninner], ecx
	add   edx, 0
	mov   [rsp + nb400_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400_unroll_loop
	jmp   .nb400_finish_inner
.nb400_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb400_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb400_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isaj
	mov rsi, [rbp + nb400_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb400_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all isaj in xmm3 
	mulps  xmm2, xmm3
	
	movaps [rsp + nb400_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb400_gbtsc]
	movaps [rsp + nb400_gbscale], xmm1
	
	mov rsi, [rbp + nb400_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb400_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb400_qq], xmm3	

	
	mov rsi, [rbp + nb400_pos]       ;# base of pos[] 

	mov [rsp + nb400_jnra], eax
	mov [rsp + nb400_jnrb], ebx
	mov [rsp + nb400_jnrc], ecx
	mov [rsp + nb400_jnrd], edx
	
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	
	lea   ecx, [ecx + ecx*2]    
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

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [rsp + nb400_ix]
	movaps xmm5, [rsp + nb400_iy]
	movaps xmm6, [rsp + nb400_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb400_dx], xmm4
	movaps [rsp + nb400_dy], xmm5
	movaps [rsp + nb400_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb400_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb400_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb400_r], xmm4
	mulps xmm4, [rsp + nb400_gbscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  rsi, [rbp + nb400_GBtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	;# load coulomb table
	movaps xmm4, [rsi + rax*4]
	movaps xmm5, [rsi + rbx*4]
	movaps xmm6, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm3, xmm6
	shufps xmm3, xmm7, 0xEE 
	shufps xmm6, xmm7, 0x44
	movaps xmm7, xmm4
	shufps xmm7, xmm5, 0xEE
	shufps xmm4, xmm5, 0x44
	movaps xmm5, xmm4
	shufps xmm5, xmm6, 0xDD
	shufps xmm4, xmm6, 0x88
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x88
	shufps xmm7, xmm3, 0xDD
	;# coulomb table ready, in xmm4-xmm7  	
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [rsp + nb400_two]	;# two*Heps2 
	movaps xmm3, [rsp + nb400_qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC

	;# get jnr from stack
	mov eax, [rsp + nb400_jnra]
	mov ebx, [rsp + nb400_jnrb]
	mov ecx, [rsp + nb400_jnrc]
	mov edx, [rsp + nb400_jnrd]
	
	mov rsi, [rbp + nb400_dvda]
	
	;# Calculate dVda
	xorps  xmm7, xmm7
	mulps xmm3, [rsp + nb400_gbscale]
	movaps xmm6, xmm3
	mulps  xmm6, [rsp + nb400_r]
	addps  xmm6, xmm5
	addps  xmm5, [rsp + nb400_vctot]
	movaps [rsp + nb400_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum
	addps  xmm7, [rsp + nb400_dvdasum]
	movaps [rsp + nb400_dvdasum], xmm7 

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	movaps  xmm5, xmm6
	movaps  xmm4, xmm7
	shufps  xmm5, xmm5, 0x1
	shufps  xmm4, xmm4, 0x1
	;# xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
	addss  xmm6, [rsi + rax*4]
	addss  xmm5, [rsi + rbx*4]
	addss  xmm7, [rsi + rcx*4]
	addss  xmm4, [rsi + rdx*4]
	movss  [rsi + rax*4], xmm6
	movss  [rsi + rbx*4], xmm5
	movss  [rsi + rcx*4], xmm7
	movss  [rsi + rdx*4], xmm4
	
	xorps  xmm4, xmm4	
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [rsp + nb400_dx]
	movaps xmm1, [rsp + nb400_dy]
	movaps xmm2, [rsp + nb400_dz]

	movd eax, mm0	
	movd ebx, mm1
	movd ecx, mm2
	movd edx, mm3

	mov    rdi, [rbp + nb400_faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [rsp + nb400_fix]
	movaps xmm4, [rsp + nb400_fiy]
	movaps xmm5, [rsp + nb400_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [rsp + nb400_fix], xmm3
	movaps [rsp + nb400_fiy], xmm4
	movaps [rsp + nb400_fiz], xmm5
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm4, [rdi + rax*4]
	movlps xmm6, [rdi + rcx*4]
	movhps xmm4, [rdi + rbx*4]
	movhps xmm6, [rdi + rdx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 136  ;# 10001000
	shufps xmm4, xmm6, 221  ;# 11011101			      

	;# now xmm3-xmm5 contains fjx, fjy, fjz 
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	;# unpack them back so we can store them - first x & y in xmm3/xmm4 

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	;# xmm6(l)=x & y for j1, (h) for j2 
	;# xmm3(l)=x & y for j3, (h) for j4 
	movlps [rdi + rax*4], xmm6
	movlps [rdi + rcx*4], xmm3
	
	movhps [rdi + rbx*4], xmm6
	movhps [rdi + rdx*4], xmm3

	;# and the z forces 
	movss  xmm4, [rdi + rax*4 + 8]
	movss  xmm5, [rdi + rbx*4 + 8]
	movss  xmm6, [rdi + rcx*4 + 8]
	movss  xmm7, [rdi + rdx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 229  ;# 11100101
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 234  ;# 11101010
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 255  ;# 11111111
	subss  xmm7, xmm2
	movss  [rdi + rax*4 + 8], xmm4
	movss  [rdi + rbx*4 + 8], xmm5
	movss  [rdi + rcx*4 + 8], xmm6
	movss  [rdi + rdx*4 + 8], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb400_innerk],  4
	jl    .nb400_finish_inner
	jmp   .nb400_unroll_loop
.nb400_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb400_innerk],  4
	mov   edx, [rsp + nb400_innerk]
	and   edx, 2
	jnz   .nb400_dopair
	jmp   .nb400_checksingle
.nb400_dopair:	
	mov   rcx, [rsp + nb400_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb400_innerjjnr],  8

	xorps xmm2, xmm2
	movaps xmm6, xmm2
	
	;# load isaj
	mov rsi, [rbp + nb400_invsqrta]
	movss xmm2, [rsi + rax*4]
	movss xmm3, [rsi + rbx*4]
	unpcklps xmm2, xmm3	;# isaj in xmm2(0,1)
	mulps  xmm2, [rsp + nb400_isai]
	movaps [rsp + nb400_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb400_gbtsc]
	movaps [rsp + nb400_gbscale], xmm1	
	
	mov rsi, [rbp + nb400_charge]    ;# base of charge[] 	
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	unpcklps xmm3, xmm6 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm2, [rsp + nb400_iq]
	mulps  xmm3, xmm2
	movaps [rsp + nb400_qq], xmm3

	mov rdi, [rbp + nb400_pos]	

	movd  mm0, eax		;# copy jnr to mm0/mm1
	movd  mm1, ebx
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	mov    rdi, [rbp + nb400_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb400_ix]
	movaps xmm5, [rsp + nb400_iy]
	movaps xmm6, [rsp + nb400_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb400_dx], xmm4
	movaps [rsp + nb400_dy], xmm5
	movaps [rsp + nb400_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb400_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb400_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [rsp + nb400_r], xmm4
	mulps xmm4, [rsp + nb400_gbscale]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb400_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	;# load coulomb table
	movaps xmm4, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm6, xmm4
	unpcklps xmm4, xmm7  	;# Y1 Y2 F1 F2 
	unpckhps xmm6, xmm7     ;# G1 G2 H1 H2
	movhlps  xmm5, xmm4    	;# F1 F2 
	movhlps  xmm7, xmm6     ;# H1 H2
	;# coulomb table ready, in xmm4-xmm7  	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm7, [rsp + nb400_two]	;# two*Heps2 
	movaps xmm3, [rsp + nb400_qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC
	
	;# get jnr from mm0/mm1
	movd ecx, mm0
	movd edx, mm1

	mov rsi, [rbp + nb400_dvda]
	
	;# Calculate dVda
	xorps  xmm7, xmm7
	mulps xmm3, [rsp + nb400_gbscale]
	movaps xmm6, xmm3
	mulps  xmm6, [rsp + nb400_r]
	addps  xmm6, xmm5
	addps  xmm5, [rsp + nb400_vctot]
	movaps [rsp + nb400_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum
	addps  xmm7, [rsp + nb400_dvdasum]
	movaps [rsp + nb400_dvdasum], xmm7 

	;# update j atoms dvdaj
	movaps xmm7, xmm6
	shufps xmm7, xmm7, 0x1
	addss  xmm6, [rsi + rcx*4]
	addss  xmm7, [rsi + rdx*4]
	movss  [rsi + rcx*4], xmm6
	movss  [rsi + rdx*4], xmm7
	
	xorps  xmm4, xmm4	
	mulps xmm3, xmm0
	subps  xmm4, xmm3

	movaps xmm0, [rsp + nb400_dx]
	movaps xmm1, [rsp + nb400_dy]
	movaps xmm2, [rsp + nb400_dz]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [rsp + nb400_fix]
	movaps xmm4, [rsp + nb400_fiy]
	movaps xmm5, [rsp + nb400_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [rsp + nb400_fix], xmm3
	movaps [rsp + nb400_fiy], xmm4
	movaps [rsp + nb400_fiz], xmm5
	;# update the fj's 
	movss   xmm3, [rdi + rax*4]
	movss   xmm4, [rdi + rax*4 + 4]
	movss   xmm5, [rdi + rax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [rdi + rax*4], xmm3
	movss   [rdi + rax*4 + 4], xmm4
	movss   [rdi + rax*4 + 8], xmm5	

	shufps  xmm0, xmm0, 225  ;# 11100001
	shufps  xmm1, xmm1, 225  ;# 11100001
	shufps  xmm2, xmm2, 225  ;# 11100001

	movss   xmm3, [rdi + rbx*4]
	movss   xmm4, [rdi + rbx*4 + 4]
	movss   xmm5, [rdi + rbx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [rdi + rbx*4], xmm3
	movss   [rdi + rbx*4 + 4], xmm4
	movss   [rdi + rbx*4 + 8], xmm5	

.nb400_checksingle:				
	mov   edx, [rsp + nb400_innerk]
	and   edx, 1
	jnz    .nb400_dosingle
	jmp    .nb400_updateouterdata
.nb400_dosingle:
	mov rsi, [rbp + nb400_charge]
	mov rdx, [rbp + nb400_invsqrta]
	mov rdi, [rbp + nb400_pos]
	mov   rcx, [rsp + nb400_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm2, xmm2
	movaps xmm6, xmm2
	movss xmm2, [rdx + rax*4]	;# isaj
	mulss xmm2, [rsp + nb400_isai]
	movss [rsp + nb400_isaprod], xmm2	
	movss xmm1, xmm2
	mulss xmm1, [rsp + nb400_gbtsc]
	movss [rsp + nb400_gbscale], xmm1	
	
	mulss  xmm2, [rsp + nb400_iq]
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulss  xmm6, xmm2
	movss [rsp + nb400_qq], xmm6

	movd  mm0, eax
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movss xmm4, [rsp + nb400_ix]
	movss xmm5, [rsp + nb400_iy]
	movss xmm6, [rsp + nb400_iz]

	;# calc dr 
	subss xmm4, xmm0
	subss xmm5, xmm1
	subss xmm6, xmm2

	;# store dr 
	movss [rsp + nb400_dx], xmm4
	movss [rsp + nb400_dy], xmm5
	movss [rsp + nb400_dz], xmm6
	;# square it 
	mulss xmm4,xmm4
	mulss xmm5,xmm5
	mulss xmm6,xmm6
	addss xmm4, xmm5
	addss xmm4, xmm6
	;# rsq in xmm4 

	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movss xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [rsp + nb400_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [rsp + nb400_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 

	mulss xmm4, xmm0	;# xmm4=r 
	movss [rsp + nb400_r], xmm4
	mulss xmm4, [rsp + nb400_gbscale]

	cvttss2si ebx, xmm4     ;# mm6 contain lu indices 
	cvtsi2ss xmm6, ebx
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 

	shl  ebx, 2

	mov  rsi, [rbp + nb400_GBtab]

	movaps xmm4, [rsi + rbx*4]	
	movhlps xmm6, xmm4
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm7, [rsp + nb400_two]	;# two*Heps2 
	movss xmm3, [rsp + nb400_qq]
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, xmm3 ;# vcoul=qq*VV  
	mulss  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 

	movd ebx, mm0
	mov rsi, [rbp + nb400_dvda]
	
	;# Calculate dVda
	xorps xmm7, xmm7
	mulss xmm3, [rsp + nb400_gbscale]
	movaps xmm6, xmm3
	mulss  xmm6, [rsp + nb400_r]
	addss  xmm6, xmm5
	addss  xmm5, [rsp + nb400_vctot]
	movss [rsp + nb400_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum
	addps  xmm7, [rsp + nb400_dvdasum]
	movaps [rsp + nb400_dvdasum], xmm7 

	;# update j atoms dvdaj
	addss  xmm6, [rsi + rbx*4]
	movss  [rsi + rbx*4], xmm6
	
	xorps  xmm4, xmm4	
	mulss xmm3, xmm0
	subss  xmm4, xmm3

	mov    rdi, [rbp + nb400_faction]

	movss xmm0, [rsp + nb400_dx]
	movss xmm1, [rsp + nb400_dy]
	movss xmm2, [rsp + nb400_dz]

	mulss  xmm0, xmm4
	mulss  xmm1, xmm4
	mulss  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movss xmm3, [rsp + nb400_fix]
	movss xmm4, [rsp + nb400_fiy]
	movss xmm5, [rsp + nb400_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss [rsp + nb400_fix], xmm3
	movss [rsp + nb400_fiy], xmm4
	movss [rsp + nb400_fiz], xmm5
	;# update fj 
	
	movss   xmm3, [rdi + rax*4]
	movss   xmm4, [rdi + rax*4 + 4]
	movss   xmm5, [rdi + rax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2	
	movss   [rdi + rax*4], xmm3
	movss   [rdi + rax*4 + 4], xmm4
	movss   [rdi + rax*4 + 8], xmm5	
.nb400_updateouterdata:
	mov   ecx, [rsp + nb400_ii3]
	mov   rdi, [rbp + nb400_faction]
	mov   rsi, [rbp + nb400_fshift]
	mov   edx, [rsp + nb400_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb400_fix]
	movaps xmm1, [rsp + nb400_fiy]
	movaps xmm2, [rsp + nb400_fiz]

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

	;# increment fshift force  
	movss  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 4]
	movss  xmm5, [rsi + rdx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsi + rdx*4],     xmm3
	movss  [rsi + rdx*4 + 4], xmm4
	movss  [rsi + rdx*4 + 8], xmm5

	;# get n from stack
	mov esi, [rsp + nb400_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb400_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb400_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb400_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate dVda and update it 
	movaps xmm7, [rsp + nb400_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
	
	mov edx, [rsp + nb400_ii]
	mov rax, [rbp + nb400_dvda]
	addss xmm7, [rax + rdx*4]
	movss [rax + rdx*4], xmm7
	
        ;# finish if last 
        mov ecx, [rsp + nb400_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb400_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb400_n], esi
        jmp .nb400_outer
.nb400_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb400_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb400_end
        ;# non-zero, do one more workunit
        jmp   .nb400_threadloop
.nb400_end:

	mov eax, [rsp + nb400_nouter]
	mov ebx, [rsp + nb400_ninner]
	mov rcx, [rbp + nb400_outeriter]
	mov rdx, [rbp + nb400_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 456
	femms

	pop rbx
	pop	rbp
	ret


	

.globl nb_kernel400nf_x86_64_sse
.globl _nb_kernel400nf_x86_64_sse
nb_kernel400nf_x86_64_sse:	
_nb_kernel400nf_x86_64_sse:	
.equiv          nb400nf_fshift,         16
.equiv          nb400nf_gid,            24
.equiv          nb400nf_pos,            32
.equiv          nb400nf_faction,        40
.equiv          nb400nf_charge,         48
.equiv          nb400nf_p_facel,        56
.equiv          nb400nf_argkrf,         64
.equiv          nb400nf_argcrf,         72
.equiv          nb400nf_Vc,             80
.equiv          nb400nf_type,           88
.equiv          nb400nf_p_ntype,        96
.equiv          nb400nf_vdwparam,       104
.equiv          nb400nf_Vvdw,           112
.equiv          nb400nf_p_tabscale,     120
.equiv          nb400nf_VFtab,          128
.equiv          nb400nf_invsqrta,       136
.equiv          nb400nf_dvda,           144
.equiv          nb400nf_p_gbtabscale,   152
.equiv          nb400nf_GBtab,          160
.equiv          nb400nf_p_nthreads,     168
.equiv          nb400nf_count,          176
.equiv          nb400nf_mtx,            184
.equiv          nb400nf_outeriter,      192
.equiv          nb400nf_inneriter,      200
.equiv          nb400nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb400nf_ix,             0
.equiv          nb400nf_iy,             16
.equiv          nb400nf_iz,             32
.equiv          nb400nf_iq,             48
.equiv          nb400nf_gbtsc,          64
.equiv          nb400nf_qq,             80
.equiv          nb400nf_vctot,          96
.equiv          nb400nf_half,           112
.equiv          nb400nf_three,          128
.equiv          nb400nf_isai,           144
.equiv          nb400nf_isaprod,        160
.equiv          nb400nf_gbscale,        176
.equiv          nb400nf_nri,            192
.equiv          nb400nf_iinr,           200
.equiv          nb400nf_jindex,         208
.equiv          nb400nf_jjnr,           216
.equiv          nb400nf_shift,          224
.equiv          nb400nf_shiftvec,       232
.equiv          nb400nf_facel,          240
.equiv          nb400nf_innerjjnr,      248
.equiv          nb400nf_is3,            256
.equiv          nb400nf_ii3,            260
.equiv          nb400nf_innerk,         264
.equiv          nb400nf_n,              268
.equiv          nb400nf_nn1,            272
.equiv          nb400nf_nouter,         276
.equiv          nb400nf_ninner,         280


	push rbp
	mov  rbp, rsp
	push rbx

	
	femms
	sub rsp, 296		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb400nf_nouter], eax
	mov [rsp + nb400nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb400nf_nri], edi
	mov [rsp + nb400nf_iinr], rsi
	mov [rsp + nb400nf_jindex], rdx
	mov [rsp + nb400nf_jjnr], rcx
	mov [rsp + nb400nf_shift], r8
	mov [rsp + nb400nf_shiftvec], r9
	mov rsi, [rbp + nb400nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb400nf_facel], xmm0

	mov rbx, [rbp + nb400nf_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb400nf_gbtsc],  xmm4



	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb400nf_half], eax
	movss xmm1, [rsp + nb400nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb400nf_half],  xmm1
	movaps [rsp + nb400nf_three],  xmm3

.nb400nf_threadloop:
        mov   rsi, [rbp + nb400nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb400nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb400nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb400nf_n], eax
        mov [rsp + nb400nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400nf_outerstart
        jmp .nb400nf_end

.nb400nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb400nf_nouter]
	mov [rsp + nb400nf_nouter], ebx

.nb400nf_outer:
	mov   rax, [rsp + nb400nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb400nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb400nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb400nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 
	
	mov   rdx, [rbp + nb400nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb400nf_facel]
	shufps xmm3, xmm3, 0

	mov   rdx, [rbp + nb400nf_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb400nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb400nf_iq], xmm3
	movaps [rsp + nb400nf_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb400nf_ix], xmm0
	movaps [rsp + nb400nf_iy], xmm1
	movaps [rsp + nb400nf_iz], xmm2

	mov   [rsp + nb400nf_ii3], ebx
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb400nf_vctot], xmm4
	
	mov   rax, [rsp + nb400nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb400nf_pos]
	mov   rdi, [rbp + nb400nf_faction]	
	mov   rax, [rsp + nb400nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb400nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb400nf_ninner]
	mov   [rsp + nb400nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb400nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400nf_unroll_loop
	jmp   .nb400nf_finish_inner
.nb400nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb400nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb400nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isa2
	mov rsi, [rbp + nb400nf_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb400nf_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm2, xmm3
	
	movaps [rsp + nb400nf_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb400nf_gbtsc]
	movaps [rsp + nb400nf_gbscale], xmm1
	
	mov rsi, [rbp + nb400nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb400nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb400nf_qq], xmm3	

	
	mov rsi, [rbp + nb400nf_pos]       ;# base of pos[] 

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

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [rsp + nb400nf_ix]
	movaps xmm5, [rsp + nb400nf_iy]
	movaps xmm6, [rsp + nb400nf_iz]

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
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb400nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb400nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	mulps xmm4, [rsp + nb400nf_gbscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  rsi, [rbp + nb400nf_GBtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	;# load coulomb table
	movaps xmm4, [rsi + rax*4]
	movaps xmm5, [rsi + rbx*4]
	movaps xmm6, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm3, xmm6
	shufps xmm3, xmm7, 0xEE 
	shufps xmm6, xmm7, 0x44
	movaps xmm7, xmm4
	shufps xmm7, xmm5, 0xEE
	shufps xmm4, xmm5, 0x44
	movaps xmm5, xmm4
	shufps xmm5, xmm6, 0xDD
	shufps xmm4, xmm6, 0x88
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x88
	shufps xmm7, xmm3, 0xDD
	;# coulomb table ready, in xmm4-xmm7  	
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 
	movaps xmm3, [rsp + nb400nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	addps  xmm5, [rsp + nb400nf_vctot]
	movaps [rsp + nb400nf_vctot], xmm5 
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb400nf_innerk],  4
	jl    .nb400nf_finish_inner
	jmp   .nb400nf_unroll_loop
.nb400nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb400nf_innerk],  4
	mov   edx, [rsp + nb400nf_innerk]
	and   edx, 2
	jnz   .nb400nf_dopair
	jmp   .nb400nf_checksingle
.nb400nf_dopair:	
	mov   rcx, [rsp + nb400nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb400nf_innerjjnr],  8

	xorps xmm2, xmm2
	movaps xmm6, xmm2
	
	;# load isa2
	mov rsi, [rbp + nb400nf_invsqrta]
	movss xmm2, [rsi + rax*4]
	movss xmm3, [rsi + rbx*4]
	unpcklps xmm2, xmm3	;# isa2 in xmm3(0,1)
	mulps  xmm2, [rsp + nb400nf_isai]
	movaps [rsp + nb400nf_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb400nf_gbtsc]
	movaps [rsp + nb400nf_gbscale], xmm1	
	
	mov rsi, [rbp + nb400nf_charge]    ;# base of charge[] 	
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	unpcklps xmm3, xmm6 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm2, [rsp + nb400nf_iq]
	mulps  xmm3, xmm2
	movaps [rsp + nb400nf_qq], xmm3

	mov rdi, [rbp + nb400nf_pos]	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	mov    rdi, [rbp + nb400nf_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb400nf_ix]
	movaps xmm5, [rsp + nb400nf_iy]
	movaps xmm6, [rsp + nb400nf_iz]

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
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb400nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb400nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb400nf_gbscale]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb400nf_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	;# load coulomb table
	movaps xmm4, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm6, xmm4
	unpcklps xmm4, xmm7  	;# Y1 Y2 F1 F2 
	unpckhps xmm6, xmm7     ;# G1 G2 H1 H2
	movhlps  xmm5, xmm4    	;# F1 F2 
	movhlps  xmm7, xmm6     ;# H1 H2
	;# coulomb table ready, in xmm4-xmm7  	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb400nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	addps  xmm5, [rsp + nb400nf_vctot]
	movaps [rsp + nb400nf_vctot], xmm5 

.nb400nf_checksingle:				
	mov   edx, [rsp + nb400nf_innerk]
	and   edx, 1
	jnz    .nb400nf_dosingle
	jmp    .nb400nf_updateouterdata
.nb400nf_dosingle:
	mov rsi, [rbp + nb400nf_charge]
	mov rdx, [rbp + nb400nf_invsqrta]
	mov rdi, [rbp + nb400nf_pos]
	mov   rcx, [rsp + nb400nf_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm2, xmm2
	movaps xmm6, xmm2
	movss xmm2, [rdx + rax*4]	;# isa2
	mulss xmm2, [rsp + nb400nf_isai]
	movss [rsp + nb400nf_isaprod], xmm2	
	movss xmm1, xmm2
	mulss xmm1, [rsp + nb400nf_gbtsc]
	movss [rsp + nb400nf_gbscale], xmm1	
	
	mulss  xmm2, [rsp + nb400nf_iq]
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulss  xmm6, xmm2
	movss [rsp + nb400nf_qq], xmm6
	
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movss xmm4, [rsp + nb400nf_ix]
	movss xmm5, [rsp + nb400nf_iy]
	movss xmm6, [rsp + nb400nf_iz]

	;# calc dr 
	subss xmm4, xmm0
	subss xmm5, xmm1
	subss xmm6, xmm2

	;# square it 
	mulss xmm4,xmm4
	mulss xmm5,xmm5
	mulss xmm6,xmm6
	addss xmm4, xmm5
	addss xmm4, xmm6
	;# rsq in xmm4 

	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movss xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [rsp + nb400nf_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [rsp + nb400nf_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 

	mulss xmm4, xmm0	;# xmm4=r 
	mulss xmm4, [rsp + nb400nf_gbscale]

	cvttss2si ebx, xmm4     ;# mm6 contain lu indices 
	cvtsi2ss xmm6, ebx
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 

	shl  ebx, 2

	mov  rsi, [rbp + nb400nf_GBtab]

	movaps xmm4, [rsi + rbx*4]	
	movhlps xmm6, xmm4
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	movss xmm3, [rsp + nb400nf_qq]
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, xmm3 ;# vcoul=qq*VV  
	addss  xmm5, [rsp + nb400nf_vctot]
	movss [rsp + nb400nf_vctot], xmm5 
.nb400nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb400nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb400nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb400nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb400nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb400nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb400nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb400nf_n], esi
        jmp .nb400nf_outer
.nb400nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb400nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb400nf_end
        ;# non-zero, do one more workunit
        jmp   .nb400nf_threadloop
.nb400nf_end:

	mov eax, [rsp + nb400nf_nouter]
	mov ebx, [rsp + nb400nf_ninner]
	mov rcx, [rbp + nb400nf_outeriter]
	mov rdx, [rbp + nb400nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 296
	femms

	pop rbx
	pop	rbp
	ret

