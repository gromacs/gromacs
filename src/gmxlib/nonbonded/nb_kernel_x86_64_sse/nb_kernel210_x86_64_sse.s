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




.globl nb_kernel210_x86_64_sse
.globl _nb_kernel210_x86_64_sse
nb_kernel210_x86_64_sse:	
_nb_kernel210_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb210_fshift,           16
.equiv          nb210_gid,              24
.equiv          nb210_pos,              32
.equiv          nb210_faction,          40
.equiv          nb210_charge,           48
.equiv          nb210_p_facel,          56
.equiv          nb210_argkrf,           64
.equiv          nb210_argcrf,           72
.equiv          nb210_Vc,               80
.equiv          nb210_type,             88
.equiv          nb210_p_ntype,          96
.equiv          nb210_vdwparam,         104
.equiv          nb210_Vvdw,             112
.equiv          nb210_p_tabscale,       120
.equiv          nb210_VFtab,            128
.equiv          nb210_invsqrta,         136
.equiv          nb210_dvda,             144
.equiv          nb210_p_gbtabscale,     152
.equiv          nb210_GBtab,            160
.equiv          nb210_p_nthreads,       168
.equiv          nb210_count,            176
.equiv          nb210_mtx,              184
.equiv          nb210_outeriter,        192
.equiv          nb210_inneriter,        200
.equiv          nb210_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb210_ix,               0
.equiv          nb210_iy,               16
.equiv          nb210_iz,               32
.equiv          nb210_iq,               48
.equiv          nb210_dx,               64
.equiv          nb210_dy,               80
.equiv          nb210_dz,               96
.equiv          nb210_c6,               112
.equiv          nb210_c12,              128
.equiv          nb210_six,              144
.equiv          nb210_twelve,           160
.equiv          nb210_vctot,            176
.equiv          nb210_Vvdwtot,          192
.equiv          nb210_fix,              208
.equiv          nb210_fiy,              224
.equiv          nb210_fiz,              240
.equiv          nb210_half,             256
.equiv          nb210_three,            272
.equiv          nb210_two,              288
.equiv          nb210_krf,              304
.equiv          nb210_crf,              320
.equiv          nb210_nri,              336
.equiv          nb210_iinr,             344
.equiv          nb210_jindex,           352
.equiv          nb210_jjnr,             360
.equiv          nb210_shift,            368
.equiv          nb210_shiftvec,         376
.equiv          nb210_facel,            384
.equiv          nb210_innerjjnr,        392
.equiv          nb210_is3,              400
.equiv          nb210_ii3,              404
.equiv          nb210_ntia,             408
.equiv          nb210_innerk,           412
.equiv          nb210_n,                416
.equiv          nb210_nn1,              420
.equiv          nb210_ntype,            424
.equiv          nb210_nouter,           428
.equiv          nb210_ninner,           432


	push rbp
	mov  rbp, rsp
	push rbx

	
	femms
	sub rsp, 440		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb210_nouter], eax
	mov [rsp + nb210_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb210_nri], edi
	mov [rsp + nb210_iinr], rsi
	mov [rsp + nb210_jindex], rdx
	mov [rsp + nb210_jjnr], rcx
	mov [rsp + nb210_shift], r8
	mov [rsp + nb210_shiftvec], r9
	mov rdi, [rbp + nb210_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb210_ntype], edi
	mov rsi, [rbp + nb210_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb210_facel], xmm0


	mov rsi, [rbp + nb210_argkrf]
	mov rdi, [rbp + nb210_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb210_krf], xmm1
	movaps [rsp + nb210_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb210_half], eax
	movss xmm1, [rsp + nb210_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + nb210_half],  xmm1
	movaps [rsp + nb210_two],  xmm2
	movaps [rsp + nb210_three],  xmm3
	movaps [rsp + nb210_six],  xmm4
	movaps [rsp + nb210_twelve],  xmm5

.nb210_threadloop:
        mov   rsi, [rbp + nb210_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb210_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb210_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb210_n], eax
        mov [rsp + nb210_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210_outerstart
        jmp .nb210_end

.nb210_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb210_nouter]
	mov [rsp + nb210_nouter], ebx

.nb210_outer:
	mov   rax, [rsp + nb210_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb210_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb210_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb210_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb210_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb210_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb210_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb210_ntype]
    	shl   edx, 1
    	mov   [rsp + nb210_ntia], edx
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb210_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb210_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb210_ix], xmm0
	movaps [rsp + nb210_iy], xmm1
	movaps [rsp + nb210_iz], xmm2

	mov   [rsp + nb210_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb210_vctot], xmm4
	movaps [rsp + nb210_Vvdwtot], xmm4
	movaps [rsp + nb210_fix], xmm4
	movaps [rsp + nb210_fiy], xmm4
	movaps [rsp + nb210_fiz], xmm4
	
	mov   rax, [rsp + nb210_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb210_pos]
	mov   rdi, [rbp + nb210_faction]	
	mov   rax, [rsp + nb210_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb210_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb210_ninner]
	mov   [rsp + nb210_ninner], ecx
	add   edx, 0
	mov   [rsp + nb210_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210_unroll_loop
	jmp   .nb210_finish_inner
.nb210_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb210_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb210_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb210_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb210_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov rsi, [rbp + nb210_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb210_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb210_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb210_c6], xmm4
	movaps [rsp + nb210_c12], xmm6
	
	mov rsi, [rbp + nb210_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
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
	movaps xmm4, [rsp + nb210_ix]
	movaps xmm5, [rsp + nb210_iy]
	movaps xmm6, [rsp + nb210_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb210_dx], xmm4
	movaps [rsp + nb210_dy], xmm5
	movaps [rsp + nb210_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 
	
	movaps xmm7, [rsp + nb210_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210_crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulps  xmm7, [rsp + nb210_two]
	mulps  xmm1, [rsp + nb210_c6]
	mulps  xmm2, [rsp + nb210_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb210_Vvdwtot]
	mulps  xmm1, [rsp + nb210_six]
	mulps  xmm2, [rsp + nb210_twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm6, [rsp + nb210_vctot]

	movaps xmm0, [rsp + nb210_dx]
	movaps xmm1, [rsp + nb210_dy]
	movaps xmm2, [rsp + nb210_dz]

	movaps [rsp + nb210_vctot], xmm6
	movaps [rsp + nb210_Vvdwtot], xmm5

	mov    rdi, [rbp + nb210_faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [rsp + nb210_fix]
	movaps xmm4, [rsp + nb210_fiy]
	movaps xmm5, [rsp + nb210_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [rsp + nb210_fix], xmm3
	movaps [rsp + nb210_fiy], xmm4
	movaps [rsp + nb210_fiz], xmm5
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
	sub dword ptr [rsp + nb210_innerk],  4
	jl    .nb210_finish_inner
	jmp   .nb210_unroll_loop
.nb210_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb210_innerk],  4
	mov   edx, [rsp + nb210_innerk]
	and   edx, 2
	jnz   .nb210_dopair
	jmp   .nb210_checksingle
.nb210_dopair:	
	mov rsi, [rbp + nb210_charge]

    mov   rcx, [rsp + nb210_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb210_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 12 ;# 00001100 
	shufps xmm3, xmm3, 88 ;# 01011000 ;# xmm3(0,1) has the charges 

	mov rsi, [rbp + nb210_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb210_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb210_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb210_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb210_c6], xmm4
	movaps [rsp + nb210_c12], xmm6	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb210_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	mov    rdi, [rbp + nb210_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb210_ix]
	movaps xmm5, [rsp + nb210_iy]
	movaps xmm6, [rsp + nb210_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb210_dx], xmm4
	movaps [rsp + nb210_dy], xmm5
	movaps [rsp + nb210_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [rsp + nb210_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210_crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulps  xmm7, [rsp + nb210_two]	
	mulps  xmm1, [rsp + nb210_c6]
	mulps  xmm2, [rsp + nb210_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb210_Vvdwtot]
	mulps  xmm1, [rsp + nb210_six]
	mulps  xmm2, [rsp + nb210_twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0	
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm6, [rsp + nb210_vctot]

	movaps xmm0, [rsp + nb210_dx]
	movaps xmm1, [rsp + nb210_dy]
	movaps xmm2, [rsp + nb210_dz]

	movaps [rsp + nb210_vctot], xmm6
	movaps [rsp + nb210_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [rsp + nb210_fix]
	movaps xmm4, [rsp + nb210_fiy]
	movaps xmm5, [rsp + nb210_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [rsp + nb210_fix], xmm3
	movaps [rsp + nb210_fiy], xmm4
	movaps [rsp + nb210_fiz], xmm5
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

.nb210_checksingle:				
	mov   edx, [rsp + nb210_innerk]
	and   edx, 1
	jnz    .nb210_dosingle
	jmp    .nb210_updateouterdata
.nb210_dosingle:			
	mov rsi, [rbp + nb210_charge]
	mov rdi, [rbp + nb210_pos]
	mov   rcx, [rsp + nb210_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [rcx]
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 	

	mov rsi, [rbp + nb210_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb210_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb210_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb210_c6], xmm4
	movaps [rsp + nb210_c12], xmm6	
	
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb210_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb210_ix]
	movaps xmm5, [rsp + nb210_iy]
	movaps xmm6, [rsp + nb210_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [rsp + nb210_dx], xmm4
	movaps [rsp + nb210_dy], xmm5
	movaps [rsp + nb210_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [rsp + nb210_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210_crf]	
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul 
	mulps  xmm7, [rsp + nb210_two]
	mulps  xmm1, [rsp + nb210_c6]
	mulps  xmm2, [rsp + nb210_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [rsp + nb210_Vvdwtot]
	mulps  xmm1, [rsp + nb210_six]
	mulps  xmm2, [rsp + nb210_twelve]
	subps  xmm2, xmm1
	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addss  xmm6, [rsp + nb210_vctot]
	
	mov    rdi, [rbp + nb210_faction]

	movaps xmm0, [rsp + nb210_dx]
	movaps xmm1, [rsp + nb210_dy]
	movaps xmm2, [rsp + nb210_dz]

	movss [rsp + nb210_vctot], xmm6
	movss [rsp + nb210_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [rsp + nb210_fix]
	movaps xmm4, [rsp + nb210_fiy]
	movaps xmm5, [rsp + nb210_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [rsp + nb210_fix], xmm3
	movaps [rsp + nb210_fiy], xmm4
	movaps [rsp + nb210_fiz], xmm5
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
.nb210_updateouterdata:
	mov   ecx, [rsp + nb210_ii3]
	mov   rdi, [rbp + nb210_faction]
	mov   rsi, [rbp + nb210_fshift]
	mov   edx, [rsp + nb210_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb210_fix]
	movaps xmm1, [rsp + nb210_fiy]
	movaps xmm2, [rsp + nb210_fiz]

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
	mov esi, [rsp + nb210_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb210_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb210_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb210_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb210_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb210_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb210_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb210_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210_n], esi
        jmp .nb210_outer
.nb210_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb210_end
        ;# non-zero, do one more workunit
        jmp   .nb210_threadloop
.nb210_end:
	mov eax, [rsp + nb210_nouter]
	mov ebx, [rsp + nb210_ninner]
	mov rcx, [rbp + nb210_outeriter]
	mov rdx, [rbp + nb210_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 440
	femms

	pop rbx
	pop	rbp
	ret







.globl nb_kernel210nf_x86_64_sse
.globl _nb_kernel210nf_x86_64_sse
nb_kernel210nf_x86_64_sse:	
_nb_kernel210nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb210nf_fshift,         16
.equiv          nb210nf_gid,            24
.equiv          nb210nf_pos,            32
.equiv          nb210nf_faction,        40
.equiv          nb210nf_charge,         48
.equiv          nb210nf_p_facel,        56
.equiv          nb210nf_argkrf,         64
.equiv          nb210nf_argcrf,         72
.equiv          nb210nf_Vc,             80
.equiv          nb210nf_type,           88
.equiv          nb210nf_p_ntype,        96
.equiv          nb210nf_vdwparam,       104
.equiv          nb210nf_Vvdw,           112
.equiv          nb210nf_p_tabscale,     120
.equiv          nb210nf_VFtab,          128
.equiv          nb210nf_invsqrta,       136
.equiv          nb210nf_dvda,           144
.equiv          nb210nf_p_gbtabscale,   152
.equiv          nb210nf_GBtab,          160
.equiv          nb210nf_p_nthreads,     168
.equiv          nb210nf_count,          176
.equiv          nb210nf_mtx,            184
.equiv          nb210nf_outeriter,      192
.equiv          nb210nf_inneriter,      200
.equiv          nb210nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb210nf_ix,             0
.equiv          nb210nf_iy,             16
.equiv          nb210nf_iz,             32
.equiv          nb210nf_iq,             48
.equiv          nb210nf_c6,             64
.equiv          nb210nf_c12,            80
.equiv          nb210nf_vctot,          96
.equiv          nb210nf_Vvdwtot,        112
.equiv          nb210nf_half,           128
.equiv          nb210nf_three,          144
.equiv          nb210nf_krf,            160
.equiv          nb210nf_crf,            176
.equiv          nb210nf_nri,            192
.equiv          nb210nf_iinr,           200
.equiv          nb210nf_jindex,         208
.equiv          nb210nf_jjnr,           216
.equiv          nb210nf_shift,          224
.equiv          nb210nf_shiftvec,       232
.equiv          nb210nf_facel,          240
.equiv          nb210nf_innerjjnr,      248
.equiv          nb210nf_is3,            256
.equiv          nb210nf_ii3,            260
.equiv          nb210nf_ntia,           264
.equiv          nb210nf_innerk,         268
.equiv          nb210nf_n,              272
.equiv          nb210nf_nn1,            276
.equiv          nb210nf_ntype,          280
.equiv          nb210nf_nouter,         284
.equiv          nb210nf_ninner,         288

	push rbp
	mov  rbp, rsp
	push rbx
	
	femms
	sub rsp, 312		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb210nf_nouter], eax
	mov [rsp + nb210nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb210nf_nri], edi
	mov [rsp + nb210nf_iinr], rsi
	mov [rsp + nb210nf_jindex], rdx
	mov [rsp + nb210nf_jjnr], rcx
	mov [rsp + nb210nf_shift], r8
	mov [rsp + nb210nf_shiftvec], r9
	mov rdi, [rbp + nb210nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb210nf_ntype], edi
	mov rsi, [rbp + nb210nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb210nf_facel], xmm0


	mov rsi, [rbp + nb210nf_argkrf]
	mov rdi, [rbp + nb210nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb210nf_krf], xmm1
	movaps [rsp + nb210nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb210nf_half], eax
	movss xmm1, [rsp + nb210nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb210nf_half],  xmm1
	movaps [rsp + nb210nf_three],  xmm3


.nb210nf_threadloop:
        mov   rsi, [rbp + nb210nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb210nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [rsi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb210nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb210nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb210nf_n], eax
        mov [rsp + nb210nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb210nf_outerstart
        jmp .nb210nf_end

.nb210nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb210nf_nouter]
	mov [rsp + nb210nf_nouter], ebx

.nb210nf_outer:
	mov   rax, [rsp + nb210nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [rsp + nb210nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb210nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb210nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb210nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb210nf_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb210nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb210nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb210nf_ntia], edx
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   rax, [rbp + nb210nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb210nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb210nf_ix], xmm0
	movaps [rsp + nb210nf_iy], xmm1
	movaps [rsp + nb210nf_iz], xmm2

	mov   [rsp + nb210nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb210nf_vctot], xmm4
	movaps [rsp + nb210nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb210nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb210nf_pos]
	mov   rax, [rsp + nb210nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb210nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb210nf_ninner]
	mov   [rsp + nb210nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb210nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb210nf_unroll_loop
	jmp   .nb210nf_finish_inner
.nb210nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb210nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb210nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb210nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb210nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov rsi, [rbp + nb210nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb210nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb210nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb210nf_c6], xmm4
	movaps [rsp + nb210nf_c12], xmm6
	
	mov rsi, [rbp + nb210nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
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
	movaps xmm4, [rsp + nb210nf_ix]
	movaps xmm5, [rsp + nb210nf_iy]
	movaps xmm6, [rsp + nb210nf_iz]

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
	
	movaps xmm7, [rsp + nb210nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210nf_crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulps  xmm1, [rsp + nb210nf_c6]
	mulps  xmm2, [rsp + nb210nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb210nf_Vvdwtot]
	addps  xmm6, [rsp + nb210nf_vctot]
	movaps [rsp + nb210nf_vctot], xmm6
	movaps [rsp + nb210nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb210nf_innerk],  4
	jl    .nb210nf_finish_inner
	jmp   .nb210nf_unroll_loop
.nb210nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb210nf_innerk],  4
	mov   edx, [rsp + nb210nf_innerk]
	and   edx, 2
	jnz   .nb210nf_dopair
	jmp   .nb210nf_checksingle
.nb210nf_dopair:	
	mov rsi, [rbp + nb210nf_charge]

    	mov   rcx, [rsp + nb210nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb210nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 12 ;# 00001100 
	shufps xmm3, xmm3, 88 ;# 01011000 ;# xmm3(0,1) has the charges 

	mov rsi, [rbp + nb210nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb210nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb210nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb210nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb210nf_c6], xmm4
	movaps [rsp + nb210nf_c12], xmm6	
	
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb210nf_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb210nf_ix]
	movaps xmm5, [rsp + nb210nf_iy]
	movaps xmm6, [rsp + nb210nf_iz]

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

	movaps xmm7, [rsp + nb210nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210nf_crf]
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulps  xmm1, [rsp + nb210nf_c6]
	mulps  xmm2, [rsp + nb210nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb210nf_Vvdwtot]
	addps  xmm6, [rsp + nb210nf_vctot]
	movaps [rsp + nb210nf_vctot], xmm6
	movaps [rsp + nb210nf_Vvdwtot], xmm5

.nb210nf_checksingle:				
	mov   edx, [rsp + nb210nf_innerk]
	and   edx, 1
	jnz    .nb210nf_dosingle
	jmp    .nb210nf_updateouterdata
.nb210nf_dosingle:			
	mov rsi, [rbp + nb210nf_charge]
	mov rdi, [rbp + nb210nf_pos]
	mov   rcx, [rsp + nb210nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [rcx]
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 	

	mov rsi, [rbp + nb210nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb210nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb210nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb210nf_c6], xmm4
	movaps [rsp + nb210nf_c12], xmm6	
	
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb210nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb210nf_ix]
	movaps xmm5, [rsp + nb210nf_iy]
	movaps xmm6, [rsp + nb210nf_iz]

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

	movaps xmm7, [rsp + nb210nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb210nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb210nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	movaps xmm1, xmm4
	subps  xmm6, [rsp + nb210nf_crf]	
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm6, xmm3	;# xmm6=vcoul 
	mulps  xmm1, [rsp + nb210nf_c6]
	mulps  xmm2, [rsp + nb210nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [rsp + nb210nf_Vvdwtot]
	addss  xmm6, [rsp + nb210nf_vctot]
	movss [rsp + nb210nf_vctot], xmm6
	movss [rsp + nb210nf_Vvdwtot], xmm5

.nb210nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb210nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb210nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb210nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb210nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb210nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb210nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb210nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb210nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210nf_n], esi
        jmp .nb210nf_outer
.nb210nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb210nf_end
        ;# non-zero, do one more workunit
        jmp   .nb210nf_threadloop
.nb210nf_end:

	mov eax, [rsp + nb210nf_nouter]
	mov ebx, [rsp + nb210nf_ninner]
	mov rcx, [rbp + nb210nf_outeriter]
	mov rdx, [rbp + nb210nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 312
	femms

	pop rbx
	pop	rbp
	ret
