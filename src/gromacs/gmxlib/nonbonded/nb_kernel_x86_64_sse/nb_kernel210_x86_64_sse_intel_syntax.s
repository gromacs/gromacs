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
.equiv          nb210_c6,               64
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
	sub rsp, 432		;# local variable stack space (n*16+8)
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
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
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
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
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
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
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
	xorps xmm12, xmm12
	movaps xmm13, xmm12
	movaps xmm14, xmm12
	movaps xmm15, xmm12
	movaps [rsp + nb210_vctot], xmm12
	movaps [rsp + nb210_Vvdwtot], xmm12
	
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
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb210_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 
    
	add qword ptr [rsp + nb210_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]    
	lea   rdx, [r11 + r11*2]	

	mov rdi, [rbp + nb210_pos]
	;# load coordinates
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rcx*4]	;# x3 y3 - - 
	movhps xmm1, [rdi + rbx*4]	;# x2 y2 - -
	movhps xmm2, [rdi + rdx*4]	;# x4 y4 - -

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rcx*4 + 8]	;# z2 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z3 - - - 
	movss xmm8, [rdi + rdx*4 + 8]	;# z4 - - - 
    movlhps xmm5, xmm7 ;# jzOa  -  jzOb  -
    movlhps xmm6, xmm8 ;# jzOc  -  jzOd -

	mov rsi, [rbp + nb210_charge]

    movaps xmm4, xmm1
    unpcklps xmm1, xmm2  ;# jxa jxc jya jyc        
    unpckhps xmm4, xmm2  ;# jxb jxd jyb jyd
    movaps xmm2, xmm1
    unpcklps xmm1, xmm4 ;# x
    unpckhps xmm2, xmm4 ;# y
    shufps   xmm5, xmm6,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

	;# calc dr  
	subps xmm1, [rsp + nb210_ix]
	subps xmm2, [rsp + nb210_iy]
	subps xmm5, [rsp + nb210_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm2
    movaps xmm11, xmm5
    
    ;# load charges
	movss xmm12, [rsi + r8*4]
	movss xmm3,  [rsi + r10*4]
	movss xmm4,  [rsi + r9*4]
	movss xmm6,  [rsi + r11*4]
    unpcklps xmm12, xmm3  ;# jqa jqc - -
    unpcklps xmm4, xmm6  ;# jqb jqd - -
    unpcklps xmm12, xmm4  ;# jqa jqb jqc jqd
	mulps xmm12, [rsp + nb210_iq]  ;# qq
    

	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5
	;# rsq in xmm1
    
 	mov rsi, [rbp + nb210_type]
	mov r12d, [rsi + r8*4]
	mov r13d, [rsi + r9*4]
	mov r14d, [rsi + r10*4]
	mov r15d, [rsi + r11*4]

    movaps xmm7, [rsp + nb210_krf] 
    mulps  xmm7, xmm1    ;# krsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm6, xmm5
	mulps xmm5, xmm5

	shl r12d, 1	
	shl r13d, 1	
	shl r14d, 1	
	shl r15d, 1	

	movaps xmm4, [rsp + nb210_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mov rsi, [rbp + nb210_vdwparam]

    mov edi, [rsp + nb210_ntia]
	add r12d, edi
	add r13d, edi
	add r14d, edi
	add r15d, edi

	mulps xmm4, xmm6	
	mulps xmm4, [rsp + nb210_half]	
	movaps xmm1, xmm4
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

    ;# load c6/c12
	movlps xmm0, [rsi + r12*4]
	movlps xmm2, [rsi + r14*4]
	movhps xmm0, [rsi + r13*4]
	movhps xmm2, [rsi + r15*4]

    movaps xmm3, xmm1
    
    subps  xmm1, xmm7
    subps  xmm1, xmm7   ;# rinv-2*krsq
    addps  xmm3, xmm7   ;# rinv+krsq
    
    subps  xmm3, [rsp + nb210_crf] ;# rinv+krsq-crf
    mulps  xmm3, xmm12  ;# vcoul=qq*(rinv+krsq-crf)
    mulps  xmm1, xmm12  ;# fijC
        
	movaps xmm5, xmm4
	mulps  xmm5, xmm4   ;# rinv4
	mulps  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulps  xmm5, xmm5	;# xmm5=rinv12
    
	movaps xmm8, xmm0
	shufps xmm0, xmm2, 136  ;# 10001000
	shufps xmm8, xmm2, 221  ;# 11011101	
    
    ;# add to vctot
    addps  xmm3, [rsp + nb210_vctot]
    movaps [rsp + nb210_vctot], xmm3
    
	mulps  xmm6, xmm0   ;# vvdw6=c6*rinv6
	mulps  xmm5, xmm8   ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subps  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulps  xmm6, [rsp + nb210_six]
	mulps  xmm7, [rsp + nb210_twelve]
	subps  xmm7, xmm6
    addps  xmm1, xmm7
	mulps  xmm4, xmm1	;# xmm4=total fscal 
    
	mov rsi, [rbp + nb210_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -

    ;# add potential to Vvdwtot 
	addps  xmm5, [rsp + nb210_Vvdwtot]
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

    movaps [rsp + nb210_Vvdwtot], xmm5

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
	
	movlps [rsi + rax*4], xmm0
	movlps [rsi + rcx*4], xmm1
	movhps [rsi + rbx*4], xmm0
	movhps [rsi + rdx*4], xmm1
    
    ;# xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd  xmm10, xmm11, 1  ;# fjz2 - - -
    movhlps xmm9,  xmm11     ;# fjz3 - - -
    pshufd  xmm8,  xmm11, 3  ;# fjz4 - - -
    
	addss  xmm11, [rsi + rax*4 + 8]
	addss  xmm10, [rsi + rbx*4 + 8]
	addss  xmm9,  [rsi + rcx*4 + 8]
	addss  xmm8,  [rsi + rdx*4 + 8]    
	movss  [rsi + rax*4 + 8], xmm11
	movss  [rsi + rbx*4 + 8], xmm10
	movss  [rsi + rcx*4 + 8], xmm9
	movss  [rsi + rdx*4 + 8], xmm8
	
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
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb210_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
    
	add qword ptr [rsp + nb210_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb210_charge]
	movss xmm12, [rsi + rax*4]
	movss xmm2, [rsi + rbx*4]

    unpcklps xmm12, xmm2  ;# jqa jqb - -
	mulps xmm12, [rsp + nb210_iq]   ;# qq

	mov rsi, [rbp + nb210_type]
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb210_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb210_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps  xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101
	    
    movaps [rsp + nb210_c6], xmm0
    movaps [rsp + nb210_c12], xmm3

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb210_pos]
    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm1, xmm2 ;# x1 x2 y1 y2
    movhlps  xmm2, xmm1 ;# y1 y2 -  -
    unpcklps xmm5, xmm6 ;# z1 z2 -  -

	;# calc dr  
	subps xmm1, [rsp + nb210_ix]
	subps xmm2, [rsp + nb210_iy]
	subps xmm5, [rsp + nb210_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm2
    movaps xmm11, xmm5
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5
	;# rsq in xmm1
    
    movaps xmm7, [rsp + nb210_krf] 
    mulps  xmm7, xmm1    ;# krsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm6, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb210_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm6	
	mulps xmm4, [rsp + nb210_half]	
	movaps xmm1, xmm4
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

    movaps xmm3, xmm1
    
    subps  xmm1, xmm7
    subps  xmm1, xmm7   ;# rinv-2*krsq
    addps  xmm3, xmm7   ;# rinv+krsq
    
    subps  xmm3, [rsp + nb210_crf] ;# rinv+krsq-crf
    mulps  xmm3, xmm12  ;# vcoul=qq*(rinv+krsq-crf)
    mulps  xmm1, xmm12  ;# fijC
        
	movaps xmm5, xmm4
	mulps  xmm5, xmm4   ;# rinv4
	mulps  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulps  xmm5, xmm5	;# xmm5=rinv12
    

    ;# add to vctot
    addps  xmm3, [rsp + nb210_vctot]
    movlps [rsp + nb210_vctot], xmm3
    
	mulps  xmm6, [rsp + nb210_c6]   ;# vvdw6=c6*rinv6
	mulps  xmm5, [rsp + nb210_c12]  ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subps  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulps  xmm6, [rsp + nb210_six]
	mulps  xmm7, [rsp + nb210_twelve]
	subps  xmm7, xmm6
    addps  xmm1, xmm7
	mulps  xmm4, xmm1	;# xmm4=total fscal 
    
    xorps  xmm7, xmm7
    movlhps xmm5, xmm7
    
    ;# add potential to Vvdwtot 
	addps  xmm5, [rsp + nb210_Vvdwtot]
    movaps [rsp + nb210_Vvdwtot], xmm5
    
    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

    movlhps xmm9, xmm7
    movlhps xmm10, xmm7
    movlhps xmm11, xmm7
    
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb210_faction]
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

.nb210_checksingle:                             
    mov   edx, [rsp + nb210_innerk]
    and   edx, 1
    jnz    .nb210_dosingle
    jmp    .nb210_updateouterdata

.nb210_dosingle:	
    mov rcx, [rsp + nb210_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb210_charge]
	movss xmm12, [rsi + rax*4]   
    mulss xmm12, [rsp + nb210_iq] ;# qq

	mov rsi, [rbp + nb210_type]
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb210_ntia]
	add r12d, edi

	mov rsi, [rbp + nb210_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]

    movaps [rsp + nb210_c6], xmm0
    movaps [rsp + nb210_c12], xmm3

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	mov rdi, [rbp + nb210_pos]
	;# load coordinates
	movss xmm1, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm2, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm5, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm1, [rsp + nb210_ix]
	subss xmm2, [rsp + nb210_iy]
	subss xmm5, [rsp + nb210_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm2
    movaps xmm11, xmm5
    
	;# square it 
	mulss xmm1,xmm1
	mulss xmm2,xmm2
	mulss xmm5,xmm5
	addss xmm1, xmm2
	addss xmm1, xmm5
	;# rsq in xmm1
    
    movaps xmm7, [rsp + nb210_krf] 
    mulss  xmm7, xmm1    ;# krsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm1
	movaps xmm6, xmm5
	mulss xmm5, xmm5
	movaps xmm4, [rsp + nb210_three]
	mulss xmm5, xmm1	;# rsq*lu*lu 	
    subss xmm4, xmm5	;# 30-rsq*lu*lu 
	mulss xmm4, xmm6	
	mulss xmm4, [rsp + nb210_half]	
	movaps xmm1, xmm4
	mulss  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

    movaps xmm3, xmm1
    
    subss  xmm1, xmm7
    subss  xmm1, xmm7   ;# rinv-2*krsq
    addss  xmm3, xmm7   ;# rinv+krsq
    
    subss  xmm3, [rsp + nb210_crf] ;# rinv+krsq-crf
    mulss  xmm3, xmm12  ;# vcoul=qq*(rinv+krsq-crf)
    mulss  xmm1, xmm12  ;# fijC
        
	movaps xmm5, xmm4
	mulss  xmm5, xmm4   ;# rinv4
	mulss  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulss  xmm5, xmm5	;# xmm5=rinv12
    

    ;# add to vctot
    addss  xmm3, [rsp + nb210_vctot]
    movss [rsp + nb210_vctot], xmm3
    
	mulss  xmm6, [rsp + nb210_c6]   ;# vvdw6=c6*rinv6
	mulss  xmm5, [rsp + nb210_c12]   ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subss  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulss  xmm6, [rsp + nb210_six]
	mulss  xmm7, [rsp + nb210_twelve]
	subss  xmm7, xmm6
    addss  xmm1, xmm7
	mulss  xmm4, xmm1	;# xmm4=total fscal 

    ;# add potential to Vvdwtot 
	addss  xmm5, [rsp + nb210_Vvdwtot]
    movss  [rsp + nb210_Vvdwtot], xmm5
    
    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm4
	mulss  xmm10, xmm4
	mulss  xmm11, xmm4
    
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb210_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb210_updateouterdata:
	mov   ecx, [rsp + nb210_ii3]
	mov   rdi, [rbp + nb210_faction]
	mov   rsi, [rbp + nb210_fshift]
	mov   edx, [rsp + nb210_is3]

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
	
	;# accumulate total potential energy and update it 
	movaps xmm12, [rsp + nb210_Vvdwtot]
    ;# accumulate
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb210_Vvdw]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb210_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb210_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210_n], esi
        jmp .nb210_outer
.nb210_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210_end
        ;# non-zero, do one more workunit
        jmp   .nb210_threadloop
.nb210_end:
	mov eax, [rsp + nb210_nouter]
	mov ebx, [rsp + nb210_ninner]
	mov rcx, [rbp + nb210_outeriter]
	mov rdx, [rbp + nb210_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 432
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
	sub rsp, 304		;# local variable stack space (n*16+8)
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
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
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
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
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
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
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

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	mulps xmm3, xmm2
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

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
	
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
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
	
	lea   rax, [rax + rax*2]
	
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
        jz .nb210nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb210nf_n], esi
        jmp .nb210nf_outer
.nb210nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb210nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb210nf_end
        ;# non-zero, do one more workunit
        jmp   .nb210nf_threadloop
.nb210nf_end:

	mov eax, [rsp + nb210nf_nouter]
	mov ebx, [rsp + nb210nf_ninner]
	mov rcx, [rbp + nb210nf_outeriter]
	mov rdx, [rbp + nb210nf_inneriter]
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
