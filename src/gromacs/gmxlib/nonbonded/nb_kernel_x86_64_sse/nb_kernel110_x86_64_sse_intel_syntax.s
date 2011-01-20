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




.globl nb_kernel110_x86_64_sse
.globl _nb_kernel110_x86_64_sse
nb_kernel110_x86_64_sse:	
_nb_kernel110_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb110_fshift,           16
.equiv          nb110_gid,              24
.equiv          nb110_pos,              32
.equiv          nb110_faction,          40
.equiv          nb110_charge,           48
.equiv          nb110_p_facel,          56
.equiv          nb110_argkrf,           64
.equiv          nb110_argcrf,           72
.equiv          nb110_Vc,               80
.equiv          nb110_type,             88
.equiv          nb110_p_ntype,          96
.equiv          nb110_vdwparam,         104
.equiv          nb110_Vvdw,             112
.equiv          nb110_p_tabscale,       120
.equiv          nb110_VFtab,            128
.equiv          nb110_invsqrta,         136
.equiv          nb110_dvda,             144
.equiv          nb110_p_gbtabscale,     152
.equiv          nb110_GBtab,            160
.equiv          nb110_p_nthreads,       168
.equiv          nb110_count,            176
.equiv          nb110_mtx,              184
.equiv          nb110_outeriter,        192
.equiv          nb110_inneriter,        200
.equiv          nb110_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb110_ix,               0
.equiv          nb110_iy,               16
.equiv          nb110_iz,               32
.equiv          nb110_iq,               48
.equiv          nb110_qq,               64
.equiv          nb110_c6,               112
.equiv          nb110_c12,              128
.equiv          nb110_six,              144
.equiv          nb110_twelve,           160
.equiv          nb110_vctot,            176
.equiv          nb110_Vvdwtot,          192
.equiv          nb110_fix,              208
.equiv          nb110_fiy,              224
.equiv          nb110_fiz,              240
.equiv          nb110_half,             256
.equiv          nb110_three,            272
.equiv          nb110_nri,              288
.equiv          nb110_iinr,             296
.equiv          nb110_jindex,           304
.equiv          nb110_jjnr,             312
.equiv          nb110_shift,            320
.equiv          nb110_shiftvec,         328
.equiv          nb110_facel,            336
.equiv          nb110_innerjjnr,        344
.equiv          nb110_is3,              352
.equiv          nb110_ii3,              356
.equiv          nb110_ntia,             360
.equiv          nb110_innerk,           364
.equiv          nb110_n,                368
.equiv          nb110_nn1,              372
.equiv          nb110_ntype,            376
.equiv          nb110_nouter,           380
.equiv          nb110_ninner,           384

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
	sub rsp, 400
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
	mov [rsp + nb110_nouter], eax
	mov [rsp + nb110_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb110_nri], edi
	mov [rsp + nb110_iinr], rsi
	mov [rsp + nb110_jindex], rdx
	mov [rsp + nb110_jjnr], rcx
	mov [rsp + nb110_shift], r8
	mov [rsp + nb110_shiftvec], r9
	mov rdi, [rbp + nb110_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb110_ntype], edi
	mov rsi, [rbp + nb110_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb110_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb110_half], eax
	movss xmm1, [rsp + nb110_half]
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
	movaps [rsp + nb110_half],  xmm1
	movaps [rsp + nb110_three],  xmm3
	movaps [rsp + nb110_six],  xmm4
	movaps [rsp + nb110_twelve],  xmm5

.nb110_threadloop:
        mov   rsi, [rbp + nb110_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb110_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb110_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb110_n], eax
        mov [rsp + nb110_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110_outerstart
        jmp .nb110_end

.nb110_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb110_nouter]
	mov [rsp + nb110_nouter], ebx

.nb110_outer:
	mov   rax, [rsp + nb110_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb110_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb110_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb110_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx +rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb110_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb110_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb110_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb110_ntype]
    	shl   edx, 1
    	mov   [rsp + nb110_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb110_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb110_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb110_ix], xmm0
	movaps [rsp + nb110_iy], xmm1
	movaps [rsp + nb110_iz], xmm2

	mov   [rsp + nb110_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm12, xmm12
	movaps [rsp + nb110_vctot], xmm12
	movaps xmm13, xmm12
	movaps xmm14, xmm12
	movaps xmm15, xmm12
	
	mov   rax, [rsp + nb110_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rax, [rsp + nb110_jjnr]

	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb110_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb110_ninner]
	mov   [rsp + nb110_ninner], ecx
	add   edx, 0
	mov   [rsp + nb110_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110_unroll_loop
	jmp   .nb110_finish_inner
.nb110_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 
    
	add qword ptr [rsp + nb110_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]    
	lea   rdx, [r11 + r11*2]	

	mov rdi, [rbp + nb110_pos]
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

	mov rsi, [rbp + nb110_charge]
	mov rdi, [rbp + nb110_type]

    movaps xmm4, xmm1
    unpcklps xmm1, xmm2  ;# jxa jxc jya jyc        
    unpckhps xmm4, xmm2  ;# jxb jxd jyb jyd

	movss xmm7, [rsi + r8*4]
	movss xmm8, [rsi + r9*4]
	movss xmm10, [rsi + r10*4]
	movss xmm11, [rsi + r11*4]

    movaps xmm2, xmm1
    unpcklps xmm1, xmm4 ;# x
    unpckhps xmm2, xmm4 ;# y
    shufps   xmm5, xmm6,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

	mov r8d, [rdi + r8*4]
	mov r9d, [rdi + r9*4]
	mov r10d, [rdi + r10*4]
	mov r11d, [rdi + r11*4]

    unpcklps xmm7, xmm10
    unpcklps xmm8, xmm11
    
	;# calc dr  
	subps xmm1, [rsp + nb110_ix]
	subps xmm2, [rsp + nb110_iy]
	subps xmm5, [rsp + nb110_iz]


    unpcklps xmm7, xmm8 ;# jq
    
	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm2
    movaps xmm11, xmm5

    mulps xmm7, [rsp + nb110_iq]
    
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	
    mov edi, [rsp + nb110_ntia]

	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5
	;# rsq in xmm1
	mov rsi, [rbp + nb110_vdwparam]
    
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm6, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb110_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm6	
	mulps xmm4, [rsp + nb110_half]	
	movaps xmm1, xmm4
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

	movlps xmm3, [rsi + r8*4]
	movlps xmm2, [rsi + r10*4]
	movhps xmm3, [rsi + r9*4]
	movhps xmm2, [rsi + r11*4]

	movaps xmm5, xmm4
	mulps  xmm5, xmm4   ;# rinv4
	mulps  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulps  xmm5, xmm5	;# xmm5=rinv12
    
    ;# coulomb stuff
    mulps  xmm1, xmm7 ;# vcoul=rinv*qq
    movaps xmm8, xmm1  ;# fijC

	movaps xmm0, xmm3
	shufps xmm0, xmm2, 136  ;# 10001000
	shufps xmm3, xmm2, 221  ;# 11011101

    ;# add to vctot
    addps  xmm1, [rsp + nb110_vctot]
    movaps [rsp + nb110_vctot], xmm1
    
	mulps  xmm6, xmm0   ;# vvdw6=c6*rinv6
	mulps  xmm5, xmm3   ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subps  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulps  xmm6, [rsp + nb110_six]
	mulps  xmm7, [rsp + nb110_twelve]
	subps  xmm7, xmm6
    addps  xmm8, xmm7
	mulps  xmm4, xmm8	;# xmm4=total fscal 
    
	mov rsi, [rbp + nb110_faction]
	;# the fj's - start by combining x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

    ;# add potential to Vvdwtot (sum in xmm12)
	addps  xmm12, xmm5

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

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
	sub dword ptr [rsp + nb110_innerk],  4
	jl    .nb110_finish_inner
	jmp   .nb110_unroll_loop
.nb110_finish_inner:
    ;# check if at least two particles remain 
    add dword ptr [rsp + nb110_innerk],  4
    mov   edx, [rsp + nb110_innerk]
    and   edx, 2
    jnz   .nb110_dopair
    jmp   .nb110_checksingle
.nb110_dopair:  
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
    
	add qword ptr [rsp + nb110_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb110_charge]
	movss xmm0, [rsi + r8*4]
	movss xmm2, [rsi + r9*4]

    unpcklps xmm0, xmm2  ;# jqa jqb - -
	mulps xmm0, [rsp + nb110_iq] 
    movaps [rsp + nb110_qq], xmm0

	mov rsi, [rbp + nb110_type]
	mov r12d, [rsi + r8*4]
	mov r13d, [rsi + r9*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb110_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb110_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps  xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101
	    
    ;# xmm0=c6
    ;# xmm3=c12

	lea   rax, [r8 + r8*2]     ;# j3 
	lea   rbx, [r9 + r9*2]	

	;# load coordinates
	mov rdi, [rbp + nb110_pos]    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm1, xmm2 ;# x1 x2 y1 y2
    movhlps  xmm2, xmm1 ;# y1 y2 -  -
    unpcklps xmm5, xmm6 ;# z1 z2 -  -

	;# calc dr  
	subps xmm1, [rsp + nb110_ix]
	subps xmm2, [rsp + nb110_iy]
	subps xmm5, [rsp + nb110_iz]

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
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm6, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb110_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm6	
	mulps xmm4, [rsp + nb110_half]	
	movaps xmm1, xmm4
	mulps  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

	movaps xmm5, xmm4
	mulps  xmm5, xmm4   ;# rinv4
	mulps  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulps  xmm5, xmm5	;# xmm5=rinv12
    
    ;# coulomb stuff
    mulps  xmm1, [rsp + nb110_qq]  ;# vcoul=rinv*qq
    movaps xmm8, xmm1  ;# fijC

    ;# add to vctot
    addps  xmm1, [rsp + nb110_vctot]
    movlps [rsp + nb110_vctot], xmm1
    
	mulps  xmm6, xmm0   ;# vvdw6=c6*rinv6
	mulps  xmm5, xmm3   ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subps  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulps  xmm6, [rsp + nb110_six]
	mulps  xmm7, [rsp + nb110_twelve]
	subps  xmm7, xmm6
    addps  xmm8, xmm7
	mulps  xmm4, xmm8	;# xmm4=total fscal 
    
    xorps  xmm7, xmm7
    movlhps xmm5, xmm7
    
    ;# add potential to Vvdwtot (sum in xmm12)
	addps  xmm12, xmm5

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

    movlhps xmm9, xmm7
    movlhps xmm10, xmm7
    movlhps xmm11, xmm7
    
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb110_faction]
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

.nb110_checksingle:                             
    mov   edx, [rsp + nb110_innerk]
    and   edx, 1
    jnz    .nb110_dosingle
    jmp    .nb110_updateouterdata

.nb110_dosingle:	
    mov rcx, [rsp + nb110_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb110_charge]
	movss xmm0, [rsi + rax*4]
    mulss xmm0, [rsp + nb110_iq]
    movaps [rsp + nb110_qq], xmm0

	mov rsi, [rbp + nb110_type]
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb110_ntia]
	add r12d, edi

	mov rsi, [rbp + nb110_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]

    ;# xmm0=c6
    ;# xmm3=c12

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	mov rdi, [rbp + nb110_pos]
	;# load coordinates
	movss xmm1, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm2, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm5, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm1, [rsp + nb110_ix]
	subss xmm2, [rsp + nb110_iy]
	subss xmm5, [rsp + nb110_iz]

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
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm1
	movaps xmm6, xmm5
	mulss xmm5, xmm5
	movaps xmm4, [rsp + nb110_three]
	mulss xmm5, xmm1	;# rsq*lu*lu 	
    subss xmm4, xmm5	;# 30-rsq*lu*lu 
	mulss xmm4, xmm6	
	mulss xmm4, [rsp + nb110_half]	
	movaps xmm1, xmm4
	mulss  xmm4, xmm4	
    ;# xmm1=rinv
    ;# xmm4=rinvsq 

	movaps xmm5, xmm4
	mulss  xmm5, xmm4   ;# rinv4
	mulss  xmm5, xmm4	;# rinv6
	movaps xmm6, xmm5
	mulss  xmm5, xmm5	;# xmm5=rinv12
    
    ;# coulomb stuff
    mulss  xmm1, [rsp + nb110_qq]  ;# vcoul=rinv*qq
    movaps xmm8, xmm1  ;# fijC

    ;# add to vctot
    addss  xmm1, [rsp + nb110_vctot]
    movss [rsp + nb110_vctot], xmm1
    
	mulss  xmm6, xmm0   ;# vvdw6=c6*rinv6
	mulss  xmm5, xmm3   ;# vvdw12=c12*rinv12     
	movaps xmm7, xmm5
	subss  xmm5, xmm6	;# Vvdw=Vvdw12-Vvdw6
     
	mulss  xmm6, [rsp + nb110_six]
	mulss  xmm7, [rsp + nb110_twelve]
	subss  xmm7, xmm6
    addss  xmm8, xmm7
	mulss  xmm4, xmm8	;# xmm4=total fscal 
    
    ;# add potential to Vvdwtot (sum in xmm12)
	addss  xmm12, xmm5

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm4
	mulss  xmm10, xmm4
	mulss  xmm11, xmm4
    
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb110_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb110_updateouterdata:
	mov   ecx, [rsp + nb110_ii3]
	mov   rdi, [rbp + nb110_faction]
	mov   rsi, [rbp + nb110_fshift]
	mov   edx, [rsp + nb110_is3]

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
	mov esi, [rsp + nb110_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb110_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb110_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb110_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total potential energy and update it 
	;# accumulate 
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb110_Vvdw]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb110_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb110_n], esi
        jmp .nb110_outer
.nb110_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb110_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110_end
        ;# non-zero, do one more workunit
        jmp   .nb110_threadloop
.nb110_end:


	mov eax, [rsp + nb110_nouter]
	mov ebx, [rsp + nb110_ninner]
	mov rcx, [rbp + nb110_outeriter]
	mov rdx, [rbp + nb110_inneriter]
	mov [rcx], eax
	mov [rdx], ebx
	add rsp, 400
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


.globl nb_kernel110nf_x86_64_sse
.globl _nb_kernel110nf_x86_64_sse
nb_kernel110nf_x86_64_sse:	
_nb_kernel110nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb110nf_fshift,         16
.equiv          nb110nf_gid,            24
.equiv          nb110nf_pos,            32
.equiv          nb110nf_faction,        40
.equiv          nb110nf_charge,         48
.equiv          nb110nf_p_facel,        56
.equiv          nb110nf_argkrf,         64
.equiv          nb110nf_argcrf,         72
.equiv          nb110nf_Vc,             80
.equiv          nb110nf_type,           88
.equiv          nb110nf_p_ntype,        96
.equiv          nb110nf_vdwparam,       104
.equiv          nb110nf_Vvdw,           112
.equiv          nb110nf_p_tabscale,     120
.equiv          nb110nf_VFtab,          128
.equiv          nb110nf_invsqrta,       136
.equiv          nb110nf_dvda,           144
.equiv          nb110nf_p_gbtabscale,   152
.equiv          nb110nf_GBtab,          160
.equiv          nb110nf_p_nthreads,     168
.equiv          nb110nf_count,          176
.equiv          nb110nf_mtx,            184
.equiv          nb110nf_outeriter,      192
.equiv          nb110nf_inneriter,      200
.equiv          nb110nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb110nf_ix,             0
.equiv          nb110nf_iy,             16
.equiv          nb110nf_iz,             32
.equiv          nb110nf_iq,             48
.equiv          nb110nf_c6,             64
.equiv          nb110nf_c12,            80
.equiv          nb110nf_vctot,          96
.equiv          nb110nf_Vvdwtot,        112
.equiv          nb110nf_half,           128
.equiv          nb110nf_three,          144
.equiv          nb110nf_nri,            160
.equiv          nb110nf_iinr,           168
.equiv          nb110nf_jindex,         176
.equiv          nb110nf_jjnr,           184
.equiv          nb110nf_shift,          192
.equiv          nb110nf_shiftvec,       200
.equiv          nb110nf_facel,          208
.equiv          nb110nf_innerjjnr,      216
.equiv          nb110nf_is3,            224
.equiv          nb110nf_ii3,            228
.equiv          nb110nf_ntia,           232
.equiv          nb110nf_innerk,         236
.equiv          nb110nf_n,              240
.equiv          nb110nf_nn1,            244
.equiv          nb110nf_ntype,          248
.equiv          nb110nf_nouter,         252
.equiv          nb110nf_ninner,         256

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
	sub rsp, 272
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
	mov [rsp + nb110nf_nouter], eax
	mov [rsp + nb110nf_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb110nf_nri], edi
	mov [rsp + nb110nf_iinr], rsi
	mov [rsp + nb110nf_jindex], rdx
	mov [rsp + nb110nf_jjnr], rcx
	mov [rsp + nb110nf_shift], r8
	mov [rsp + nb110nf_shiftvec], r9
	mov rdi, [rbp + nb110nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb110nf_ntype], edi
	mov rsi, [rbp + nb110nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb110nf_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb110nf_half], eax
	movss xmm1, [rsp + nb110nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb110nf_half],  xmm1
	movaps [rsp + nb110nf_three],  xmm3

.nb110nf_threadloop:
        mov   rsi, [rbp + nb110nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb110nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb110nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb110nf_n], eax
        mov [rsp + nb110nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110nf_outerstart
        jmp .nb110nf_end

.nb110nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb110nf_nouter]
	mov [rsp + nb110nf_nouter], ebx

.nb110nf_outer:
	mov   rax, [rsp + nb110nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb110nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb110nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb110nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx +rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb110nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb110nf_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb110nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb110nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb110nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb110nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb110nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb110nf_ix], xmm0
	movaps [rsp + nb110nf_iy], xmm1
	movaps [rsp + nb110nf_iz], xmm2

	mov   [rsp + nb110nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb110nf_vctot], xmm4
	movaps [rsp + nb110nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb110nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb110nf_pos]
	mov   rax, [rsp + nb110nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb110nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb110nf_ninner]
	mov   [rsp + nb110nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb110nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110nf_unroll_loop
	jmp   .nb110nf_finish_inner
.nb110nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb110nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb110nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb110nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb110nf_iq]
	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov rsi, [rbp + nb110nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb110nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb110nf_ntia]
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

	movaps [rsp + nb110nf_c6], xmm4
	movaps [rsp + nb110nf_c12], xmm6
	
	mov rsi, [rbp + nb110nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb110nf_ix]
	movaps xmm5, [rsp + nb110nf_iy]
	movaps xmm6, [rsp + nb110nf_iz]

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
	movaps xmm1, [rsp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb110nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [rsp + nb110nf_c6]
	mulps  xmm2, [rsp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb110nf_Vvdwtot]
	addps  xmm3, [rsp + nb110nf_vctot]
	movaps [rsp + nb110nf_vctot], xmm3
	movaps [rsp + nb110nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb110nf_innerk],  4
	jl    .nb110nf_finish_inner
	jmp   .nb110nf_unroll_loop
.nb110nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb110nf_innerk],  4
	mov   edx, [rsp + nb110nf_innerk]
	and   edx, 2
	jnz   .nb110nf_dopair
	jmp   .nb110nf_checksingle
.nb110nf_dopair:	
	mov rsi, [rbp + nb110nf_charge]

    	mov   rcx, [rsp + nb110nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb110nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 12 ;# 00001100 
	shufps xmm3, xmm3, 88 ;# 01011000 ;# xmm3(0,1) has the charges 

	mov rsi, [rbp + nb110nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb110nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb110nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb110nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb110nf_c6], xmm4
	movaps [rsp + nb110nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb110nf_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb110nf_ix]
	movaps xmm5, [rsp + nb110nf_iy]
	movaps xmm6, [rsp + nb110nf_iz]

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
	movaps xmm1, [rsp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb110nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [rsp + nb110nf_c6]
	mulps  xmm2, [rsp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb110nf_Vvdwtot]
	addps  xmm3, [rsp + nb110nf_vctot]
	movaps [rsp + nb110nf_vctot], xmm3
	movaps [rsp + nb110nf_Vvdwtot], xmm5

.nb110nf_checksingle:				
	mov   edx, [rsp + nb110nf_innerk]
	and   edx, 1
	jnz    .nb110nf_dosingle
	jmp    .nb110nf_updateouterdata
.nb110nf_dosingle:			
	mov rsi, [rbp + nb110nf_charge]
	mov rdi, [rbp + nb110nf_pos]
	mov   rcx, [rsp + nb110nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [rcx]
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 	

	mov rsi, [rbp + nb110nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb110nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb110nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb110nf_c6], xmm4
	movaps [rsp + nb110nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb110nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb110nf_ix]
	movaps xmm5, [rsp + nb110nf_iy]
	movaps xmm6, [rsp + nb110nf_iz]

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
	movaps xmm1, [rsp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb110nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [rsp + nb110nf_c6]
	mulps  xmm2, [rsp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [rsp + nb110nf_Vvdwtot]
	addss  xmm3, [rsp + nb110nf_vctot]
	movss [rsp + nb110nf_vctot], xmm3
	movss [rsp + nb110nf_Vvdwtot], xmm5

.nb110nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb110nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb110nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb110nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb110nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb110nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb110nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb110nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb110nf_n], esi
        jmp .nb110nf_outer
.nb110nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb110nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110nf_end
        ;# non-zero, do one more workunit
        jmp   .nb110nf_threadloop
.nb110nf_end:


	mov eax, [rsp + nb110nf_nouter]
	mov ebx, [rsp + nb110nf_ninner]
	mov rcx, [rbp + nb110nf_outeriter]
	mov rdx, [rbp + nb110nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 272
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
