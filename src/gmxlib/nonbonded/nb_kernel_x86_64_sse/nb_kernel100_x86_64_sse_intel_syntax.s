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


 

.globl nb_kernel100_x86_64_sse
.globl _nb_kernel100_x86_64_sse
nb_kernel100_x86_64_sse:	
_nb_kernel100_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb100_fshift,           16
.equiv          nb100_gid,              24
.equiv          nb100_pos,              32
.equiv          nb100_faction,          40
.equiv          nb100_charge,           48
.equiv          nb100_p_facel,          56
.equiv          nb100_argkrf,           64
.equiv          nb100_argcrf,           72
.equiv          nb100_Vc,               80
.equiv          nb100_type,             88
.equiv          nb100_p_ntype,          96
.equiv          nb100_vdwparam,         104
.equiv          nb100_Vvdw,             112
.equiv          nb100_p_tabscale,       120
.equiv          nb100_VFtab,            128
.equiv          nb100_invsqrta,         136
.equiv          nb100_dvda,             144
.equiv          nb100_p_gbtabscale,     152
.equiv          nb100_GBtab,            160
.equiv          nb100_p_nthreads,       168
.equiv          nb100_count,            176
.equiv          nb100_mtx,              184
.equiv          nb100_outeriter,        192
.equiv          nb100_inneriter,        200
.equiv          nb100_work,             208
	;# facel, krf,crf, tabscale, gbtabscale passed in xmm regs.
 	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100_ix,               0
.equiv          nb100_iy,               16
.equiv          nb100_iz,               32
.equiv          nb100_iq,               48
.equiv          nb100_dx,               64
.equiv          nb100_dy,               80
.equiv          nb100_dz,               96
.equiv          nb100_vctot,            112
.equiv          nb100_fix,              128
.equiv          nb100_fiy,              144
.equiv          nb100_fiz,              160
.equiv          nb100_half,             176
.equiv          nb100_three,            192
.equiv          nb100_innerjjnr,        208
.equiv          nb100_iinr,             216
.equiv          nb100_jindex,           224
.equiv          nb100_jjnr,             232
.equiv          nb100_shift,            240
.equiv          nb100_shiftvec,         248
.equiv          nb100_facel,            256
.equiv          nb100_is3,              264
.equiv          nb100_ii3,              268
.equiv          nb100_innerk,           272
.equiv          nb100_n,                276
.equiv          nb100_nn1,              280
.equiv          nb100_nouter,           284
.equiv          nb100_ninner,           288
.equiv          nb100_nri,              292

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms
    sub rsp, 304	; # local variable stack space (n*16+8)                                                         
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb100_nouter], eax
	mov [rsp + nb100_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb100_nri], edi
	mov [rsp + nb100_iinr], rsi
	mov [rsp + nb100_jindex], rdx
	mov [rsp + nb100_jjnr], rcx
	mov [rsp + nb100_shift], r8
	mov [rsp + nb100_shiftvec], r9
	mov rsi, [rbp + nb100_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb100_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb100_half], eax
	movss xmm1, [rsp + nb100_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb100_half],  xmm1
	movaps [rsp + nb100_three],  xmm3

.nb100_threadloop:
    mov   rsi, [rbp + nb100_count]          ;# pointer to sync counter
    mov   eax, [rsi]
.nb100_spinlock:
    mov   ebx, eax                          ;# ebx=*count=nn0
    add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
    pause                                   ;# -> better p4 performance
    jnz .nb100_spinlock

    ;# if(nn1>nri) nn1=nri
    mov  ecx, [rsp + nb100_nri]
    mov  edx, ecx
    sub ecx, ebx
    cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
    ;# Cleared the spinlock if we got here.
    ;# eax contains nn0, ebx contains nn1.
    mov  [rsp + nb100_n], eax
    mov  [rsp + nb100_nn1], ebx
    sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax							;# copy n to esi
    jg  .nb100_outerstart
    jmp .nb100_end		

.nb100_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb100_nouter]
	mov [rsp + nb100_nouter], ebx

.nb100_outer:
	mov   rax, [rsp + nb100_shift]      	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]				;# ebx=shift[n] 

	lea	rbx, [rbx + rbx*2]    		;# rbx=3*is 
	mov    [rsp + nb100_is3],ebx    		;# store is3 

	mov   	rax, [rsp + nb100_shiftvec]   	;# eax = base of shiftvec[] 

	movss xmm10, [rax + rbx*4]
	movss xmm11, [rax + rbx*4 + 4]
	movss xmm12, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb100_iinr]       	;# rcx = pointer into iinr[]
	mov   ebx, [rcx + rsi*4]	 		   ;# ebx =ii 

	mov   rdx, [rbp + nb100_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb100_facel]
	shufps xmm3, xmm3, 0

	lea   rbx, [rbx + rbx*2]		;# rbx = 3*ii=ii3 

    
    mov   rax, [rbp + nb100_pos]    ;# rax = base of pos[]  
	addss xmm10, [rax + rbx*4]
	addss xmm11, [rax + rbx*4 + 4]
	addss xmm12, [rax + rbx*4 + 8]

	movaps [rsp + nb100_iq], xmm3
	
	shufps xmm10, xmm10, 0
	shufps xmm11, xmm11, 0
	shufps xmm12, xmm12, 0

    movaps [rsp + nb100_ix], xmm10
    movaps [rsp + nb100_iy], xmm11
    movaps [rsp + nb100_iz], xmm12
    
	mov   [rsp + nb100_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorps xmm12, xmm12
	movaps xmm13, xmm12
	movaps xmm14, xmm12
	movaps xmm15, xmm12
	
	mov   rax, [rsp + nb100_jindex]
	mov   ecx, [rax + rsi*4]		     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx    			         ;# number of innerloop atoms 

	mov   rax, [rsp + nb100_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb100_innerjjnr], rax   ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb100_ninner]
	mov   [rsp + nb100_ninner], ecx
	add   edx, 0 ;# to check sign
	mov   [rsp + nb100_innerk], edx    ;# number of innerloop atoms 

	jge   .nb100_unroll_loop
	jmp   .nb100_finish_inner

.nb100_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
    
	add qword ptr [rsp + nb100_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   r8, [rax + rax*2]     ;# j3
	lea   r9, [rbx + rbx*2]	

	lea   r10, [rcx + rcx*2]    
	lea   r11, [rdx + rdx*2]	

	mov rdi, [rbp + nb100_pos]
	;# load coordinates    
	movlps xmm1, [rdi + r8*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + r9*4]	;# x2 y2 - - 
	movlps xmm3, [rdi + r10*4]	;# x3 y3 - -
	movlps xmm4, [rdi + r11*4]	;# x4 y4 - -

	movss xmm5, [rdi + r8*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + r9*4 + 8]	;# z2 - - - 
	movss xmm7, [rdi + r10*4 + 8]	;# z3 - - - 
	movss xmm8, [rdi + r11*4 + 8]	;# z4 - - - 

    unpcklps xmm1, xmm3 ;# x1 x3 y1 y3
    unpcklps xmm2, xmm4 ;# x2 x4 y2 y4
    unpcklps xmm5, xmm7 ;# z1 z3 -  -
    unpcklps xmm6, xmm8 ;# z2 z4 -  -

    movaps xmm3, xmm1

	mov rsi, [rbp + nb100_charge]
    unpcklps xmm1, xmm2 ;# x1 x2 x3 x4
    unpckhps xmm3, xmm2 ;# y1 y2 y3 y4
    unpcklps xmm5, xmm6 ;# z1 z2 z3 z4

	movss xmm0, [rsi + rax*4]
	movss xmm2, [rsi + rcx*4]
	movss xmm7, [rsi + rbx*4]
	movss xmm8, [rsi + rdx*4]

	;# calc dr  
	subps xmm1, [rsp + nb100_ix]
	subps xmm3, [rsp + nb100_iy]
	subps xmm5, [rsp + nb100_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm3
    movaps xmm11, xmm5
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm3,xmm3
	mulps xmm5,xmm5
	addps xmm1, xmm3
	addps xmm1, xmm5
	;# rsq in xmm1

    unpcklps xmm0, xmm2  ;# jqa jqc - -
    unpcklps xmm7, xmm8  ;# jqb jqd - -

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm2, xmm5
	mulps xmm5, xmm5
    unpcklps xmm0, xmm7  ;# jqa jqb jqc jqd
	movaps xmm4, [rsp + nb100_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
    mulps xmm0, [rsp + nb100_iq]
	mulps xmm4, [rsp + nb100_half]	
	movaps xmm1, xmm4
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 
    
    ;# calculate coulomb interaction, xmm0=qq
	mulps  xmm0, xmm1	;# xmm0=vcoul 
	mulps  xmm4, xmm0	;# xmm4=fscal 

    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm0

	mov rsi, [rbp + nb100_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + r8*4] ;# x1 y1 - -
	movlps xmm1, [rsi + r10*4] ;# x3 y3 - -
	movhps xmm0, [rsi + r9*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + r11*4] ;# x3 y3 x4 y4

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

    movaps xmm8, xmm9
    unpcklps xmm9, xmm10 ;# x1 y1 x2 y2
    unpckhps xmm8, xmm10 ;# x3 y3 x4 y4
    
    ;# update fjx and fjy
	addps  xmm0, xmm9
	addps  xmm1, xmm8
	
	movlps [rsi + r8*4], xmm0
	movlps [rsi + r10*4], xmm1
	movhps [rsi + r9*4], xmm0
	movhps [rsi + r11*4], xmm1
    
    ;# xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd  xmm5, xmm11, 1  ;# fjz2 - - -
    movhlps xmm4,  xmm11     ;# fjz3 - - -
    pshufd  xmm3,  xmm11, 3  ;# fjz4 - - -
    
	addss  xmm11, [rsi + r8*4 + 8]
	addss  xmm5,  [rsi + r9*4 + 8]
	addss  xmm4,  [rsi + r10*4 + 8]
	addss  xmm3,  [rsi + r11*4 + 8]    
	movss  [rsi + r8*4 + 8], xmm11
	movss  [rsi + r9*4 + 8], xmm5
	movss  [rsi + r10*4 + 8], xmm4
	movss  [rsi + r11*4 + 8], xmm3
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb100_innerk],  4
	jl    .nb100_finish_inner
	jmp   .nb100_unroll_loop
.nb100_finish_inner:
    ;# check if at least two particles remain 
    add dword ptr [rsp + nb100_innerk],  4
    mov   edx, [rsp + nb100_innerk]
    and   edx, 2
    jnz   .nb100_dopair
    jmp   .nb100_checksingle
.nb100_dopair:  
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
    
	add qword ptr [rsp + nb100_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb100_charge]
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rbx*4]

    unpcklps xmm0, xmm1  ;# jqa jqb - -
	mulps xmm0, [rsp + nb100_iq]    ;#qq

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	mov rdi, [rbp + nb100_pos]
	;# load coordinates    
	movlps xmm4, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm5, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm6, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm4, xmm5 ;# x1 x2 y1 y2
    movhlps  xmm5, xmm4 ;# y1 y2 -  -
    unpcklps xmm6, xmm7 ;# z1 z2 -  -

	;# calc dr  
	subps xmm4, [rsp + nb100_ix]
	subps xmm5, [rsp + nb100_iy]
	subps xmm6, [rsp + nb100_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm4
    movaps xmm10, xmm5
    movaps xmm11, xmm6
    
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm4
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb100_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 	
    subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm1, [rsp + nb100_half]	
	movaps xmm4, xmm1
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

    xorps xmm6, xmm6
    
    ;# calculate coulomb interaction, xmm0=qq
	mulps  xmm0, xmm1	;# xmm0=vcoul 
	mulps  xmm4, xmm0	;# xmm4=fscal 

    movlhps xmm0, xmm6
    
    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm0
    
    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

    movlhps xmm9, xmm6
    movlhps xmm10, xmm6
    movlhps xmm11, xmm6
    
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb100_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2

    unpcklps xmm9, xmm10  ;# x1 y1 x2 y2
    addps    xmm0, xmm9

	movlps [rsi + rax*4], xmm0
	movhps [rsi + rbx*4], xmm0
    
    ;# z forces
    pshufd xmm8, xmm11, 1
    addss  xmm11, [rsi + rax*4 + 8] 
    addss  xmm8,  [rsi + rbx*4 + 8]
    movss  [rsi + rax*4 + 8], xmm11
    movss  [rsi + rbx*4 + 8], xmm8

.nb100_checksingle:                             
    mov   edx, [rsp + nb100_innerk]
    and   edx, 1
    jnz    .nb100_dosingle
    jmp    .nb100_updateouterdata

.nb100_dosingle:	
    mov rcx, [rsp + nb100_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb100_charge]
	movss xmm0, [rsi + rax*4]       ;# jq
	mulss xmm0, [rsp + nb100_iq]    ;# qq

	lea   rax, [rax + rax*2]        ;# replace jnr with j3 

	mov rdi, [rbp + nb100_pos]
	movss xmm4, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm5, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm6, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm4, [rsp + nb100_ix]
	subss xmm5, [rsp + nb100_iy]
	subss xmm6, [rsp + nb100_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm4
    movaps xmm10, xmm5
    movaps xmm11, xmm6
    
	;# square it 
	mulss xmm4,xmm4
	mulss xmm5,xmm5
	mulss xmm6,xmm6
	addss xmm4, xmm5
	addss xmm4, xmm6
	;# rsq in xmm4 

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm4
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm1, [rsp + nb100_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 	
    subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm1, [rsp + nb100_half]	
	movaps xmm4, xmm1
	mulss  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

    ;# calculate coulomb interaction, xmm0=qq
	mulss  xmm0, xmm1	;# xmm0=vcoul 
	mulss  xmm4, xmm0	;# xmm4=fscal 

    ;# add potential to vctot (sum in xmm12)
	addss  xmm12, xmm0

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm4
	mulss  xmm10, xmm4
	mulss  xmm11, xmm4

	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11
	mov rsi, [rbp + nb100_faction]

    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb100_updateouterdata:

	mov   ecx, [rsp + nb100_ii3]
	mov   rsi, [rbp + nb100_fshift]
	mov   edx, [rsp + nb100_is3]
	mov   rdi, [rbp + nb100_faction]

	;# accumulate i forces in xmm13, xmm14, xmm15
	movhlps xmm0, xmm13
	movhlps xmm1, xmm14
	movhlps xmm2, xmm15
	addps  xmm0, xmm13
	addps  xmm1, xmm14
	addps  xmm2, xmm15 
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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# increment fshift force  
	movss  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 4]
	movss  xmm5, [rsi + rdx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rsi + rdx*4],     xmm3
	movss  [rsi + rdx*4 + 4], xmm4
	movss  [rsi + rdx*4 + 8], xmm5

	;# get n from stack
	mov  esi, [rsp + nb100_n]
    ;# get group index for i particle 
    mov   rdx, [rbp + nb100_gid]      	;# base of gid[]
    mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	;# accumulate 
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb100_Vc]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	

        ;# finish if last 
        mov  ecx, [rsp + nb100_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100_outerend

        ;# not last, iterate outer loop once more!  
        mov  [rsp + nb100_n], esi
        jmp .nb100_outer
.nb100_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb100_nri]
	;# n is already loaded in esi
        sub   ecx, esi 
        jz .nb100_end
        ;# non-zero, do one more workunit
        jmp   .nb100_threadloop
.nb100_end:
	
	mov eax, [rsp + nb100_nouter]
	mov ebx, [rsp + nb100_ninner]
	mov rcx, [rbp + nb100_outeriter]
	mov rdx, [rbp + nb100_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 304
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret







.globl nb_kernel100nf_x86_64_sse
.globl _nb_kernel100nf_x86_64_sse
nb_kernel100nf_x86_64_sse:	
_nb_kernel100nf_x86_64_sse:	
.equiv          nb100nf_fshift,         16
.equiv          nb100nf_gid,            24
.equiv          nb100nf_pos,            32
.equiv          nb100nf_faction,        40
.equiv          nb100nf_charge,         48
.equiv          nb100nf_p_facel,        56
.equiv          nb100nf_argkrf,         64
.equiv          nb100nf_argcrf,         72
.equiv          nb100nf_Vc,             80
.equiv          nb100nf_type,           88
.equiv          nb100nf_p_ntype,        96
.equiv          nb100nf_vdwparam,       104
.equiv          nb100nf_Vvdw,           112
.equiv          nb100nf_p_tabscale,     120
.equiv          nb100nf_VFtab,          128
.equiv          nb100nf_invsqrta,       136
.equiv          nb100nf_dvda,           144
.equiv          nb100nf_p_gbtabscale,   152
.equiv          nb100nf_GBtab,          160
.equiv          nb100nf_p_nthreads,     168
.equiv          nb100nf_count,          176
.equiv          nb100nf_mtx,            184
.equiv          nb100nf_outeriter,      192
.equiv          nb100nf_inneriter,      200
.equiv          nb100nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100nf_ix,             0
.equiv          nb100nf_iy,             16
.equiv          nb100nf_iz,             32
.equiv          nb100nf_iq,             48
.equiv          nb100nf_vctot,          64
.equiv          nb100nf_half,           80
.equiv          nb100nf_three,          96
.equiv          nb100nf_nri,            112
.equiv          nb100nf_iinr,           120
.equiv          nb100nf_jindex,         128
.equiv          nb100nf_jjnr,           136
.equiv          nb100nf_shift,          144
.equiv          nb100nf_shiftvec,       152
.equiv          nb100nf_facel,          160
.equiv          nb100nf_innerjjnr,      168
.equiv          nb100nf_is3,            176
.equiv          nb100nf_ii3,            180
.equiv          nb100nf_innerk,         184
.equiv          nb100nf_n,              188
.equiv          nb100nf_nn1,            192
.equiv          nb100nf_nouter,         196
.equiv          nb100nf_ninner,         200

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms
    sub rsp, 208		; # local variable stack space (n*16+8)
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block


	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb100nf_nouter], eax
	mov [rsp + nb100nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb100nf_nri], edi
	mov [rsp + nb100nf_iinr], rsi
	mov [rsp + nb100nf_jindex], rdx
	mov [rsp + nb100nf_jjnr], rcx
	mov [rsp + nb100nf_shift], r8
	mov [rsp + nb100nf_shiftvec], r9
	mov rsi, [rbp + nb100nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb100nf_facel], xmm0




	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb100nf_half], eax
	movss xmm1, [rsp + nb100nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb100nf_half],  xmm1
	movaps [rsp + nb100nf_three],  xmm3

.nb100nf_threadloop:
        mov   rsi, [rbp + nb100nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb100nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb100nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb100nf_n], eax
        mov [rsp + nb100nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100nf_outerstart
        jmp .nb100nf_end

.nb100nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb100nf_nouter]
	mov [rsp + nb100nf_nouter], ebx

.nb100nf_outer:
	mov   rax, [rsp + nb100nf_shift]      	;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]			;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb100nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb100nf_shiftvec]   	;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb100nf_iinr]       	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    		;# ebx =ii 

	mov   rdx, [rbp + nb100nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb100nf_facel]
	shufps xmm3, xmm3, 0
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb100nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb100nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb100nf_ix], xmm0
	movaps [rsp + nb100nf_iy], xmm1
	movaps [rsp + nb100nf_iz], xmm2

	mov   [rsp + nb100nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb100nf_vctot], xmm4

	mov   rax, [rsp + nb100nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb100nf_pos]

	mov   rax, [rsp + nb100nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb100nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb100nf_ninner]
	mov   [rsp + nb100nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb100nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100nf_unroll_loop
	jmp   .nb100nf_finish_inner
.nb100nf_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb100nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb100nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm5, [rsp + nb100nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# 10001000	      
	mov rsi, [rbp + nb100nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	mulps xmm3, xmm5
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [rsi + rax*4]	;# x1 y1 - - 
	movlps xmm5, [rsi + rcx*4]	;# x3 y3 - - 
	movss xmm2, [rsi + rax*4 + 8]	;# z1 -  - - 
	movss xmm6, [rsi + rcx*4 + 8]   ;# z3 -  - - 

	movhps xmm4, [rsi + rbx*4]	;# x1 y1 x2 y2 
	movhps xmm5, [rsi + rdx*4]	;# x3 y3 x4 y4 

	movss xmm0, [rsi + rbx*4 + 8]	;# z2 - - - 
	movss xmm1, [rsi + rdx*4 + 8]	;# z4 - - - 

	shufps xmm2, xmm0, 0		;# z1 z1 z2 z2 
	shufps xmm6, xmm1, 0		;# z3 z3 z4 z4 
	
	movaps xmm0, xmm4		;# x1 y1 x2 y2 	
	movaps xmm1, xmm4		;# x1 y1 x2 y2 

	shufps xmm2, xmm6, 136  ;# 10001000	;# z1 z2 z3 z4 
	
	shufps xmm0, xmm5, 136  ;# 10001000	;# x1 x2 x3 x4 
	shufps xmm1, xmm5, 221  ;# 11011101	;# y1 y2 y3 y4 		

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movaps xmm4, [rsp + nb100nf_ix]
	movaps xmm5, [rsp + nb100nf_iy]
	movaps xmm6, [rsp + nb100nf_iz]

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
	movaps xmm1, [rsp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb100nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	movaps xmm5, [rsp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm5, xmm3
	movaps [rsp + nb100nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb100nf_innerk],  4
	jl    .nb100nf_finish_inner
	jmp   .nb100nf_unroll_loop
.nb100nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb100nf_innerk],  4
	mov   edx, [rsp + nb100nf_innerk]
	and   edx, 2
	jnz   .nb100nf_dopair
	jmp   .nb100nf_checksingle
.nb100nf_dopair:	
	mov rsi, [rbp + nb100nf_charge]
	mov rdi, [rbp + nb100nf_pos]
    	mov   rcx, [rsp + nb100nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb100nf_innerjjnr],  8

	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# 00001000 ;# xmm3(0,1) has the charges 

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb100nf_iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	;# move nb100nf_ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb100nf_ix]
	movaps xmm5, [rsp + nb100nf_iy]
	movaps xmm6, [rsp + nb100nf_iz]

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
	movaps xmm1, [rsp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb100nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm5, [rsp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm5, xmm3
	movaps [rsp + nb100nf_vctot], xmm5
.nb100nf_checksingle:				
	mov   edx, [rsp + nb100nf_innerk]
	and   edx, 1
	jnz    .nb100nf_dosingle
	jmp    .nb100nf_updateouterdata
.nb100nf_dosingle:			
	mov rsi, [rbp + nb100nf_charge]
	mov rdi, [rbp + nb100nf_pos]
	mov   rcx, [rsp + nb100nf_innerjjnr]
	mov   eax, [rcx]	
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb100nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb100nf_ix]
	movaps xmm5, [rsp + nb100nf_iy]
	movaps xmm6, [rsp + nb100nf_iz]

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
	movaps xmm1, [rsp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb100nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm5, [rsp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addss  xmm5, xmm3
	movaps [rsp + nb100nf_vctot], xmm5

.nb100nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb100nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb100nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb100nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb100nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb100nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb100nf_n], esi
        jmp .nb100nf_outer
.nb100nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb100nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100nf_end
        ;# non-zero, do one more workunit
        jmp   .nb100nf_threadloop
.nb100nf_end:
	
	mov eax, [rsp + nb100nf_nouter]
	mov ebx, [rsp + nb100nf_ninner]
	mov rcx, [rbp + nb100nf_outeriter]
	mov rdx, [rbp + nb100nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 208
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
