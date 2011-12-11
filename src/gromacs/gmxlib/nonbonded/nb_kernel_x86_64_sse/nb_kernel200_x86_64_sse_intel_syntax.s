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



.globl nb_kernel200_x86_64_sse
.globl _nb_kernel200_x86_64_sse
nb_kernel200_x86_64_sse:	
_nb_kernel200_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb200_fshift,           16
.equiv          nb200_gid,              24
.equiv          nb200_pos,              32
.equiv          nb200_faction,          40
.equiv          nb200_charge,           48
.equiv          nb200_p_facel,          56
.equiv          nb200_argkrf,           64
.equiv          nb200_argcrf,           72
.equiv          nb200_Vc,               80
.equiv          nb200_type,             88
.equiv          nb200_p_ntype,          96
.equiv          nb200_vdwparam,         104
.equiv          nb200_Vvdw,             112
.equiv          nb200_p_tabscale,       120
.equiv          nb200_VFtab,            128
.equiv          nb200_invsqrta,         136
.equiv          nb200_dvda,             144
.equiv          nb200_p_gbtabscale,     152
.equiv          nb200_GBtab,            160
.equiv          nb200_p_nthreads,       168
.equiv          nb200_count,            176
.equiv          nb200_mtx,              184
.equiv          nb200_outeriter,        192
.equiv          nb200_inneriter,        200
.equiv          nb208_work,             200
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200_ix,               0
.equiv          nb200_iy,               16
.equiv          nb200_iz,               32
.equiv          nb200_iq,               48
.equiv          nb200_dx,               64
.equiv          nb200_dy,               80
.equiv          nb200_dz,               96
.equiv          nb200_vctot,            112
.equiv          nb200_fix,              128
.equiv          nb200_fiy,              144
.equiv          nb200_fiz,              160
.equiv          nb200_half,             176
.equiv          nb200_three,            192
.equiv          nb200_two,              208
.equiv          nb200_krf,              224
.equiv          nb200_crf,              240
.equiv          nb200_innerjjnr,        256
.equiv          nb200_nri,              264
.equiv          nb200_iinr,             272
.equiv          nb200_jindex,           280
.equiv          nb200_jjnr,             288
.equiv          nb200_shift,            296
.equiv          nb200_shiftvec,         304
.equiv          nb200_facel,            312
.equiv          nb200_is3,              320
.equiv          nb200_ii3,              324
.equiv          nb200_innerk,           328
.equiv          nb200_n,                332
.equiv          nb200_nn1,              336
.equiv          nb200_nouter,           340
.equiv          nb200_ninner,           344

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
	sub rsp, 352		;# local variable stack space (n*16+8)
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
	mov [rsp + nb200_nouter], eax
	mov [rsp + nb200_ninner], eax

	


	mov edi, [rdi]
	mov [rsp + nb200_nri], edi
	mov [rsp + nb200_iinr], rsi
	mov [rsp + nb200_jindex], rdx
	mov [rsp + nb200_jjnr], rcx
	mov [rsp + nb200_shift], r8
	mov [rsp + nb200_shiftvec], r9
	mov rsi, [rbp + nb200_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb200_facel], xmm0


	mov rsi, [rbp + nb200_argkrf]
	mov rdi, [rbp + nb200_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb200_krf], xmm1
	movaps [rsp + nb200_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb200_half], eax
	movss xmm1, [rsp + nb200_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb200_half],  xmm1
	movaps [rsp + nb200_two],  xmm2
	movaps [rsp + nb200_three],  xmm3

.nb200_threadloop:
        mov   rsi, [rbp + nb200_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb200_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb200_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb200_n], eax
        mov [rsp + nb200_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200_outerstart
        jmp .nb200_end

	;# assume we have at least one i particle - start directly 	
.nb200_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb200_nouter]
	mov [rsp + nb200_nouter], ebx

.nb200_outer:
	mov   rax, [rsp + nb200_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb200_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb200_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb200_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb200_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb200_facel]
	shufps xmm3, xmm3, 0
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb200_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb200_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb200_ix], xmm0
	movaps [rsp + nb200_iy], xmm1
	movaps [rsp + nb200_iz], xmm2

	mov   [rsp + nb200_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorps xmm12, xmm12
	movaps xmm13, xmm12
	movaps xmm14, xmm12
	movaps xmm15, xmm12
	
	mov   rax, [rsp + nb200_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rax, [rsp + nb200_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb200_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb200_ninner]
	mov   [rsp + nb200_ninner], ecx
	add   edx, 0
	mov   [rsp + nb200_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200_unroll_loop
	jmp   .nb200_finish_inner
.nb200_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb200_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb200_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   rax, [r8 + r8*2]     ;# j3
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]     
	lea   rdx, [r11 + r11*2]	

	mov rdi, [rbp + nb200_pos]
	;# load coordinates
    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 
	movlps xmm3, [rdi + rcx*4]	;# x3 y3 - -
	movlps xmm4, [rdi + rdx*4]	;# x4 y4 - -

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 
	movss xmm7, [rdi + rcx*4 + 8]	;# z3 - - - 
	movss xmm8, [rdi + rdx*4 + 8]	;# z4 - - - 

    unpcklps xmm1, xmm3 ;# x1 x3 y1 y3
    unpcklps xmm2, xmm4 ;# x2 x4 y2 y4
    unpcklps xmm5, xmm7 ;# z1 z3 -  -
    unpcklps xmm6, xmm8 ;# z2 z4 -  -

    movaps xmm3, xmm1

    unpcklps xmm1, xmm2 ;# x1 x2 x3 x4
    unpckhps xmm3, xmm2 ;# y1 y2 y3 y4
    unpcklps xmm5, xmm6 ;# z1 z2 z3 z4
	mov rsi, [rbp + nb200_charge]

	;# calc dr  
	subps xmm1, [rsp + nb200_ix]
	subps xmm3, [rsp + nb200_iy]
	subps xmm5, [rsp + nb200_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm3
    movaps xmm11, xmm5

	movss xmm0, [rsi + r8*4]
	movss xmm2, [rsi + r10*4]
	movss xmm6, [rsi + r9*4]
	movss xmm8, [rsi + r11*4]
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm3,xmm3
	mulps xmm5,xmm5
	addps xmm1, xmm3
	addps xmm1, xmm5
	;# rsq in xmm1

	movaps xmm7, [rsp + nb200_krf]

    unpcklps xmm0, xmm2  ;# jqa jqc - -
    unpcklps xmm6, xmm8  ;# jqb jqd - -

	rsqrtps xmm5, xmm1
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb200_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 			
	movaps xmm3, [rsp + nb200_half]
	mulps  xmm7, xmm1	;# xmm7=krsq 
	subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
    
    unpcklps xmm0, xmm6  ;# jqa jqb jqc jqd
	mulps xmm0, [rsp + nb200_iq]    ;#qq

	mulps xmm3, xmm4	;# xmm3=rinv 
	movaps xmm1, xmm3
	mulps  xmm1, xmm1	;# xmm1=rinvsq 
	movaps xmm6, xmm3
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subps  xmm6, [rsp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulps  xmm6, xmm0	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addps  xmm7, xmm7

	subps  xmm3, xmm7
	mulps  xmm0, xmm3	
	mulps  xmm1, xmm0	;# xmm1=total fscal 

    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm6

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm1
	mulps  xmm10, xmm1
	mulps  xmm11, xmm1

	mov rsi, [rbp + nb200_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

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
	sub dword ptr [rsp + nb200_innerk],  4
	jl    .nb200_finish_inner
	jmp   .nb200_unroll_loop
.nb200_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb200_innerk],  4
	mov   edx, [rsp + nb200_innerk]
	and   edx, 2
	jnz   .nb200_dopair
	jmp   .nb200_checksingle
.nb200_dopair:	
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb200_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
    
	add qword ptr [rsp + nb200_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb200_charge]
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rbx*4]

    unpcklps xmm0, xmm1  ;# jqa jqb - -
	mulps xmm0, [rsp + nb200_iq]    ;#qq

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb200_pos]
    
	movlps xmm4, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm5, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm6, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm4, xmm5 ;# x1 x2 y1 y2
    movhlps  xmm5, xmm4 ;# y1 y2 -  -
    unpcklps xmm6, xmm7 ;# z1 z2 -  -

	;# calc dr  
	subps xmm4, [rsp + nb200_ix]
	subps xmm5, [rsp + nb200_iy]
	subps xmm6, [rsp + nb200_iz]

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

	movaps xmm7, [rsp + nb200_krf]

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb200_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm3, [rsp + nb200_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm3, xmm1	;# xmm3=rinv 
	movaps xmm4, xmm3
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm3
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subps  xmm6, [rsp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulps  xmm6, xmm0	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addps  xmm7, xmm7

	subps  xmm3, xmm7
	mulps  xmm0, xmm3	
	mulps  xmm4, xmm0	;# xmm4=total fscal 

    xorps  xmm5, xmm5
    movlhps xmm6, xmm5

    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm6

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

    movlhps xmm9, xmm5
    movlhps xmm10, xmm5
    movlhps xmm11, xmm5
    
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb200_faction]
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

.nb200_checksingle:				
	mov   edx, [rsp + nb200_innerk]
	and   edx, 1
	jnz    .nb200_dosingle
	jmp    .nb200_updateouterdata
.nb200_dosingle:			
    mov rcx, [rsp + nb200_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb200_charge]
	movss xmm0, [rsi + rax*4]       ;# jq
	mulss xmm0, [rsp + nb200_iq]    ;# qq

	lea   rax, [rax + rax*2]        ;# replace jnr with j3 

	mov rdi, [rbp + nb200_pos]
	movss xmm4, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm5, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm6, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm4, [rsp + nb200_ix]
	subss xmm5, [rsp + nb200_iy]
	subss xmm6, [rsp + nb200_iz]

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

	movaps xmm7, [rsp + nb200_krf]

	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm1, [rsp + nb200_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm3, [rsp + nb200_half]
	mulss  xmm7, xmm4	;# xmm7=krsq 
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm3, xmm1	;# xmm3=rinv 
	movaps xmm4, xmm3
	mulss  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm3
	addss  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subss  xmm6, [rsp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulss  xmm6, xmm0	;# xmm6=vcoul=qq*(rinv+krsq-crf) 
	addss  xmm7, xmm7

	subss  xmm3, xmm7
	mulss  xmm0, xmm3	
	mulss  xmm4, xmm0	;# xmm4=total fscal 

    ;# add potential to vctot (sum in xmm12)
	addss  xmm12, xmm6

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm4
	mulss  xmm10, xmm4
	mulss  xmm11, xmm4

	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb200_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb200_updateouterdata:
	mov   ecx, [rsp + nb200_ii3]
	mov   rdi, [rbp + nb200_faction]
	mov   rsi, [rbp + nb200_fshift]
	mov   edx, [rsp + nb200_is3]

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
	mov esi, [rsp + nb200_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb200_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	;# accumulate 
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb200_Vc]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb200_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb200_n], esi
        jmp .nb200_outer
.nb200_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb200_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200_end
        ;# non-zero, do one more workunit
        jmp   .nb200_threadloop
.nb200_end:

	mov eax, [rsp + nb200_nouter]
	mov ebx, [rsp + nb200_ninner]
	mov rcx, [rbp + nb200_outeriter]
	mov rdx, [rbp + nb200_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 352
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



.globl nb_kernel200nf_x86_64_sse
.globl _nb_kernel200nf_x86_64_sse
nb_kernel200nf_x86_64_sse:	
_nb_kernel200nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb200nf_fshift,         16
.equiv          nb200nf_gid,            24
.equiv          nb200nf_pos,            32
.equiv          nb200nf_faction,        40
.equiv          nb200nf_charge,         48
.equiv          nb200nf_p_facel,        56
.equiv          nb200nf_argkrf,         64
.equiv          nb200nf_argcrf,         72
.equiv          nb200nf_Vc,             80
.equiv          nb200nf_type,           88
.equiv          nb200nf_p_ntype,        96
.equiv          nb200nf_vdwparam,       104
.equiv          nb200nf_Vvdw,           112
.equiv          nb200nf_p_tabscale,     120
.equiv          nb200nf_VFtab,          128
.equiv          nb200nf_invsqrta,       136
.equiv          nb200nf_dvda,           144
.equiv          nb200nf_p_gbtabscale,   152
.equiv          nb200nf_GBtab,          160
.equiv          nb200nf_p_nthreads,     168
.equiv          nb200nf_count,          176
.equiv          nb200nf_mtx,            184
.equiv          nb200nf_outeriter,      192
.equiv          nb200nf_inneriter,      200
.equiv          nb208nf_work,           200
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200nf_ix,             0
.equiv          nb200nf_iy,             16
.equiv          nb200nf_iz,             32
.equiv          nb200nf_iq,             48
.equiv          nb200nf_vctot,          64
.equiv          nb200nf_half,           80
.equiv          nb200nf_three,          96
.equiv          nb200nf_krf,            112
.equiv          nb200nf_crf,            128
.equiv          nb200nf_is3,            144
.equiv          nb200nf_ii3,            148
.equiv          nb200nf_innerjjnr,      152
.equiv          nb200nf_nri,            160
.equiv          nb200nf_iinr,           168
.equiv          nb200nf_jindex,         176
.equiv          nb200nf_jjnr,           184
.equiv          nb200nf_shift,          192
.equiv          nb200nf_shiftvec,       200
.equiv          nb200nf_facel,          208
.equiv          nb200nf_innerk,         216
.equiv          nb200nf_n,              220
.equiv          nb200nf_nn1,            224
.equiv          nb200nf_nouter,         228
.equiv          nb200nf_ninner,         232

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
	sub rsp, 240		;# local variable stack space (n*16+8)
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
	mov [rsp + nb200nf_nouter], eax
	mov [rsp + nb200nf_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb200nf_nri], edi
	mov [rsp + nb200nf_iinr], rsi
	mov [rsp + nb200nf_jindex], rdx
	mov [rsp + nb200nf_jjnr], rcx
	mov [rsp + nb200nf_shift], r8
	mov [rsp + nb200nf_shiftvec], r9
	mov rsi, [rbp + nb200nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb200nf_facel], xmm0

	mov rsi, [rbp + nb200nf_argkrf]
	mov rdi, [rbp + nb200nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb200nf_krf], xmm1
	movaps [rsp + nb200nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb200nf_half], eax
	movss xmm1, [rsp + nb200nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb200nf_half],  xmm1
	movaps [rsp + nb200nf_three],  xmm3

.nb200nf_threadloop:
        mov   rsi, [rbp + nb200nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb200nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb200nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb200nf_n], eax
        mov [rsp + nb200nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200nf_outerstart
        jmp .nb200nf_end

.nb200nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb200nf_nouter]
	mov [rsp + nb200nf_nouter], ebx

.nb200nf_outer:
	mov   rax, [rsp + nb200nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb200nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb200nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb200nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb200nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb200nf_facel]
	shufps xmm3, xmm3, 0
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb200nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb200nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb200nf_ix], xmm0
	movaps [rsp + nb200nf_iy], xmm1
	movaps [rsp + nb200nf_iz], xmm2

	mov   [rsp + nb200nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb200nf_vctot], xmm4
	
	mov   rax, [rsp + nb200nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb200nf_pos]
	mov   rax, [rsp + nb200nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb200nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb200nf_ninner]
	mov   [rsp + nb200nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb200nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200nf_unroll_loop
	jmp   .nb200nf_finish_inner
.nb200nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb200nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb200nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb200nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb200nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  

	mov rsi, [rbp + nb200nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb200nf_ix]
	movaps xmm5, [rsp + nb200nf_iy]
	movaps xmm6, [rsp + nb200nf_iz]

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
	
	movaps xmm7, [rsp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addps  xmm6, [rsp + nb200nf_vctot]
	movaps [rsp + nb200nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb200nf_innerk],  4
	jl    .nb200nf_finish_inner
	jmp   .nb200nf_unroll_loop
.nb200nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb200nf_innerk],  4
	mov   edx, [rsp + nb200nf_innerk]
	and   edx, 2
	jnz   .nb200nf_dopair
	jmp   .nb200nf_checksingle
.nb200nf_dopair:	
	mov rsi, [rbp + nb200nf_charge]

    	mov   rcx, [rsp + nb200nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb200nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 12 ;# 00001100 
	shufps xmm3, xmm3, 88 ;# 01011000 ;# xmm3(0,1) has the charges 	

	mov rdi, [rbp + nb200nf_pos]	
	
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb200nf_iq]

	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	;# move ix-iz to xmm4-xmm6 
	
	movaps xmm4, [rsp + nb200nf_ix]
	movaps xmm5, [rsp + nb200nf_iy]
	movaps xmm6, [rsp + nb200nf_iz]

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

	movaps xmm7, [rsp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 

	addps  xmm6, [rsp + nb200nf_vctot]
	movaps [rsp + nb200nf_vctot], xmm6

.nb200nf_checksingle:				
	mov   edx, [rsp + nb200nf_innerk]
	and   edx, 1
	jnz    .nb200nf_dosingle
	jmp    .nb200nf_updateouterdata
.nb200nf_dosingle:			
	mov rsi, [rbp + nb200nf_charge]
	mov rdi, [rbp + nb200nf_pos]
	mov   rcx, [rsp + nb200nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [rcx]
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb200nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb200nf_ix]
	movaps xmm5, [rsp + nb200nf_iy]
	movaps xmm6, [rsp + nb200nf_iz]

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

	movaps xmm7, [rsp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul 
	addss  xmm6, [rsp + nb200nf_vctot]
	movss [rsp + nb200nf_vctot], xmm6

.nb200nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb200nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb200nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb200nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb200nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb200nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb200nf_n], esi
        jmp .nb200nf_outer
.nb200nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb200nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200nf_end
        ;# non-zero, do one more workunit
        jmp   .nb200nf_threadloop
.nb200nf_end:

	mov eax, [rsp + nb200nf_nouter]
	mov ebx, [rsp + nb200nf_ninner]
	mov rcx, [rbp + nb200nf_outeriter]
	mov rdx, [rbp + nb200nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 240
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
