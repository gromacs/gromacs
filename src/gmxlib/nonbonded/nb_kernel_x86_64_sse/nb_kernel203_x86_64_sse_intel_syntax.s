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



.globl nb_kernel203_x86_64_sse
.globl _nb_kernel203_x86_64_sse
nb_kernel203_x86_64_sse:	
_nb_kernel203_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb203_fshift,           16
.equiv          nb203_gid,              24
.equiv          nb203_pos,              32
.equiv          nb203_faction,          40
.equiv          nb203_charge,           48
.equiv          nb203_p_facel,          56
.equiv          nb203_argkrf,           64
.equiv          nb203_argcrf,           72
.equiv          nb203_Vc,               80
.equiv          nb203_type,             88
.equiv          nb203_p_ntype,          96
.equiv          nb203_vdwparam,         104
.equiv          nb203_Vvdw,             112
.equiv          nb203_p_tabscale,       120
.equiv          nb203_VFtab,            128
.equiv          nb203_invsqrta,         136
.equiv          nb203_dvda,             144
.equiv          nb203_p_gbtabscale,     152
.equiv          nb203_GBtab,            160
.equiv          nb203_p_nthreads,       168
.equiv          nb203_count,            176
.equiv          nb203_mtx,              184
.equiv          nb203_outeriter,        192
.equiv          nb203_inneriter,        200
.equiv          nb203_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb203_ixH1,             0
.equiv          nb203_iyH1,             16
.equiv          nb203_izH1,             32
.equiv          nb203_ixH2,             48
.equiv          nb203_iyH2,             64
.equiv          nb203_izH2,             80
.equiv          nb203_ixM,              96
.equiv          nb203_iyM,              112
.equiv          nb203_izM,              128
.equiv          nb203_iqH,              144
.equiv          nb203_iqM,              160
.equiv          nb203_dxH1,             176
.equiv          nb203_dyH1,             192
.equiv          nb203_dzH1,             208
.equiv          nb203_dxH2,             224
.equiv          nb203_dyH2,             240
.equiv          nb203_dzH2,             256
.equiv          nb203_dxM,              272
.equiv          nb203_dyM,              288
.equiv          nb203_dzM,              304
.equiv          nb203_qqH,              320
.equiv          nb203_qqM,              336
.equiv          nb203_vctot,            352
.equiv          nb203_fixH1,            384
.equiv          nb203_fiyH1,            400
.equiv          nb203_fizH1,            416
.equiv          nb203_fixH2,            432
.equiv          nb203_fiyH2,            448
.equiv          nb203_fizH2,            464
.equiv          nb203_fixM,             480
.equiv          nb203_fiyM,             496
.equiv          nb203_fizM,             512
.equiv          nb203_fjx,              528
.equiv          nb203_fjy,              544
.equiv          nb203_fjz,              560
.equiv          nb203_half,             576
.equiv          nb203_three,            592
.equiv          nb203_two,              608
.equiv          nb203_krf,              624
.equiv          nb203_crf,              640
.equiv          nb203_krsqH1,           656
.equiv          nb203_krsqH2,           672
.equiv          nb203_krsqM,            688
.equiv          nb203_is3,              704
.equiv          nb203_ii3,              708
.equiv          nb203_innerjjnr,        712
.equiv          nb203_nri,              720
.equiv          nb203_iinr,             728
.equiv          nb203_jindex,           736
.equiv          nb203_jjnr,             744
.equiv          nb203_shift,            752
.equiv          nb203_shiftvec,         760
.equiv          nb203_facel,            768
.equiv          nb203_innerk,           776
.equiv          nb203_n,                780
.equiv          nb203_nn1,              784
.equiv          nb203_nouter,           788
.equiv          nb203_ninner,           792

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

	sub rsp, 800		;# local variable stack space (n*16+8)
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
	mov [rsp + nb203_nouter], eax
	mov [rsp + nb203_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb203_nri], edi
	mov [rsp + nb203_iinr], rsi
	mov [rsp + nb203_jindex], rdx
	mov [rsp + nb203_jjnr], rcx
	mov [rsp + nb203_shift], r8
	mov [rsp + nb203_shiftvec], r9
	mov rsi, [rbp + nb203_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb203_facel], xmm0


	mov rsi, [rbp + nb203_argkrf]
	mov rdi, [rbp + nb203_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb203_krf], xmm1
	movaps [rsp + nb203_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb203_half], eax
	movss xmm1, [rsp + nb203_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb203_half],  xmm1
	movaps [rsp + nb203_two],  xmm2
	movaps [rsp + nb203_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb203_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb203_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb203_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb203_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb203_iqH], xmm3
	movaps [rsp + nb203_iqM], xmm4
	
.nb203_threadloop:
        mov   rsi, [rbp + nb203_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb203_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb203_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb203_n], eax
        mov [rsp + nb203_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203_outerstart
        jmp .nb203_end
	
.nb203_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb203_nouter]
	mov [rsp + nb203_nouter], ebx

.nb203_outer:
	mov   rax, [rsp + nb203_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb203_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb203_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb203_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb203_pos]	;# rax = base of pos[]  
	mov   [rsp + nb203_ii3], ebx

	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb203_ixH1], xmm3
	movaps [rsp + nb203_iyH1], xmm4
	movaps [rsp + nb203_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb203_ixH2], xmm0
	movaps [rsp + nb203_iyH2], xmm1
	movaps [rsp + nb203_izH2], xmm2
	movaps [rsp + nb203_ixM], xmm3
	movaps [rsp + nb203_iyM], xmm4
	movaps [rsp + nb203_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb203_vctot], xmm4
	movaps [rsp + nb203_fixH1], xmm4
	movaps [rsp + nb203_fiyH1], xmm4
	movaps [rsp + nb203_fizH1], xmm4
	movaps [rsp + nb203_fixH2], xmm4
	movaps [rsp + nb203_fiyH2], xmm4
	movaps [rsp + nb203_fizH2], xmm4
	movaps [rsp + nb203_fixM], xmm4
	movaps [rsp + nb203_fiyM], xmm4
	movaps [rsp + nb203_fizM], xmm4
	
	mov   rax, [rsp + nb203_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb203_pos]
	mov   rdi, [rbp + nb203_faction]	
	mov   rax, [rsp + nb203_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb203_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb203_ninner]
	mov   [rsp + nb203_ninner], ecx
	add   edx, 0
	mov   [rsp + nb203_innerk], edx	;# number of innerloop atoms 
	jge   .nb203_unroll_loop
	jmp   .nb203_odd_inner
.nb203_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb203_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb203_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb203_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb203_iqH]
	mulps  xmm4, [rsp + nb203_iqM]

	movaps  [rsp + nb203_qqH], xmm3
	movaps  [rsp + nb203_qqM], xmm4

	mov rsi, [rbp + nb203_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
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

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb203_ixH1]
    subps xmm1, [rsp + nb203_iyH1]
    subps xmm2, [rsp + nb203_izH1]
    subps xmm3, [rsp + nb203_ixH2]
    subps xmm4, [rsp + nb203_iyH2]
    subps xmm5, [rsp + nb203_izH2]
    subps xmm6, [rsp + nb203_ixM]
    subps xmm7, [rsp + nb203_iyM]
    subps xmm8, [rsp + nb203_izM]
    
	movaps [rsp + nb203_dxH1], xmm0
	movaps [rsp + nb203_dyH1], xmm1
	movaps [rsp + nb203_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb203_dxH2], xmm3
	movaps [rsp + nb203_dyH2], xmm4
	movaps [rsp + nb203_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb203_dxM], xmm6
	movaps [rsp + nb203_dyM], xmm7
	movaps [rsp + nb203_dzM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for j atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb203_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm4, [rsp + nb203_half]
	mulps   xmm9, xmm4  ;# rinvH1 
	mulps   xmm10, xmm4 ;# rinvH2
    mulps   xmm11, xmm4 ;# rinvM
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb203_krf]    
    mulps  xmm9, xmm9   ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, xmm2  ;# k*rsq
    mulps  xmm3, xmm2
    mulps  xmm6, xmm2
    movaps xmm2, xmm0 ;# copy of k*rsq
    movaps xmm5, xmm3
    movaps xmm8, xmm6
    addps  xmm2, xmm1  ;# rinv+krsq
    addps  xmm5, xmm4
    addps  xmm8, xmm7
    movaps xmm14, [rsp + nb203_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb203_qqH]
    movaps xmm13, [rsp + nb203_qqM]    
    mulps  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm12
    mulps  xmm7, xmm13
    addps  xmm2, [rsp + nb203_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps [rsp + nb203_vctot], xmm2
    
    mulps  xmm1, xmm9   ;# fscal
    mulps  xmm4, xmm10
    mulps  xmm7, xmm11

	;# move j forces to local temp variables 
    movlps xmm9, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm10, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm9, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm10, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm11, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm12, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm6,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm8,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm11, xmm6 ;# jza  -  jzb  -
    movlhps xmm12, xmm8 ;# jzc  -  jzd -
    
    shufps xmm11, xmm12,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm9: jxa jya jxb jyb 
    ;# xmm10: jxc jyc jxd jyd
    ;# xmm11: jza jzb jzc jzd

    movaps xmm0, xmm1
    movaps xmm2, xmm1
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm7
    movaps xmm8, xmm7

	mulps xmm0, [rsp + nb203_dxH1]
	mulps xmm1, [rsp + nb203_dyH1]
	mulps xmm2, [rsp + nb203_dzH1]
	mulps xmm3, [rsp + nb203_dxH2]
	mulps xmm4, [rsp + nb203_dyH2]
	mulps xmm5, [rsp + nb203_dzH2]
	mulps xmm6, [rsp + nb203_dxM]
	mulps xmm7, [rsp + nb203_dyM]
	mulps xmm8, [rsp + nb203_dzM]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb203_fixH1]
    addps xmm1, [rsp + nb203_fiyH1]
    addps xmm2, [rsp + nb203_fizH1]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb203_fixH2]
    addps xmm4, [rsp + nb203_fiyH2]
    addps xmm5, [rsp + nb203_fizH2]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb203_fixM]
    addps xmm7, [rsp + nb203_fiyM]
    addps xmm8, [rsp + nb203_fizM]

    movaps [rsp + nb203_fixH1], xmm0
    movaps [rsp + nb203_fiyH1], xmm1
    movaps [rsp + nb203_fizH1], xmm2
    movaps [rsp + nb203_fixH2], xmm3
    movaps [rsp + nb203_fiyH2], xmm4
    movaps [rsp + nb203_fizH2], xmm5
    movaps [rsp + nb203_fixM], xmm6
    movaps [rsp + nb203_fiyM], xmm7
    movaps [rsp + nb203_fizM], xmm8
    
    ;# xmm9 = fjx
    ;# xmm10 = fjy
    ;# xmm11 = fjz
    movaps xmm15, xmm13
    unpcklps xmm13, xmm14
    unpckhps xmm15, xmm14
    
    addps xmm9, xmm13
    addps xmm10, xmm15

    movhlps  xmm12, xmm11 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm9
    movhps [rdi + rbx*4], xmm9
    movlps [rdi + rcx*4], xmm10
    movhps [rdi + rdx*4], xmm10
    movss  [rdi + rax*4 + 8], xmm11
    movss  [rdi + rcx*4 + 8], xmm12
    shufps xmm11, xmm11, 1
    shufps xmm12, xmm12, 1
    movss  [rdi + rbx*4 + 8], xmm11
    movss  [rdi + rdx*4 + 8], xmm12

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb203_innerk],  4
	jl    .nb203_odd_inner
	jmp   .nb203_unroll_loop
.nb203_odd_inner:	
	add dword ptr [rsp + nb203_innerk],  4
	jnz   .nb203_odd_loop
	jmp   .nb203_updateouterdata
.nb203_odd_loop:
	mov   rdx, [rsp + nb203_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb203_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb203_iqM]
	mov rsi, [rbp + nb203_charge] 
	movhps xmm4, [rsp + nb203_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb203_qqM], xmm3	;# use dummy qq for storage 

	mov rsi, [rbp + nb203_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb203_ixM]
	movss xmm1, [rsp + nb203_iyM]
	movss xmm2, [rsp + nb203_izM]
	
	movlps xmm6, [rsp + nb203_ixH1]
	movlps xmm7, [rsp + nb203_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb203_iyH1]
	movlps xmm7, [rsp + nb203_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb203_izH1]
	movlps xmm7, [rsp + nb203_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	;# use dummy dx for storage
	movaps [rsp + nb203_dxM], xmm3
	movaps [rsp + nb203_dyM], xmm4
	movaps [rsp + nb203_dzM], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb203_krf]
	movaps [rsp + nb203_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb203_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb203_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000

	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm1, xmm0	;# xmm1=rinv 
	movaps xmm3, [rsp + nb203_krsqM]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [rsp + nb203_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [rsp + nb203_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [rsp + nb203_qqM]	;# xmm0=vcoul 
	mulps  xmm1, [rsp + nb203_qqM] 	;# xmm1=coul part of fs 

	
	mulps  xmm4, xmm1	;# xmm4=total fscal 
	addps  xmm0, [rsp + nb203_vctot]
	movaps [rsp + nb203_vctot], xmm0
	
	movaps xmm0, [rsp + nb203_dxM]
	movaps xmm1, [rsp + nb203_dyM]
	movaps xmm2, [rsp + nb203_dzM]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb203_fixM]	
	movss  xmm4, [rsp + nb203_fiyM]	
	movss  xmm5, [rsp + nb203_fizM]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb203_fixM], xmm3	
	movss  [rsp + nb203_fiyM], xmm4	
	movss  [rsp + nb203_fizM], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb203_fixH1]
	addss  xmm4, [rsp + nb203_fiyH1]
	addss  xmm5, [rsp + nb203_fizH1]
	movss  [rsp + nb203_fixH1], xmm3	
	movss  [rsp + nb203_fiyH1], xmm4	
	movss  [rsp + nb203_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb203_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb203_fixH2]
	addss  xmm4, [rsp + nb203_fiyH2]
	addss  xmm5, [rsp + nb203_fizH2]
	movss  [rsp + nb203_fixH2], xmm3	
	movss  [rsp + nb203_fiyH2], xmm4	
	movss  [rsp + nb203_fizH2], xmm5	;# updated the H2 force 

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
	addss    xmm2, xmm1	;# z sum in xmm2 
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb203_innerk]
	jz    .nb203_updateouterdata
	jmp   .nb203_odd_loop
.nb203_updateouterdata:
	mov   ecx, [rsp + nb203_ii3]
	mov   rdi, [rbp + nb203_faction]
	mov   rsi, [rbp + nb203_fshift]
	mov   edx, [rsp + nb203_is3]

	;# accumulate  H1 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb203_fixH1]
	movaps xmm1, [rsp + nb203_fiyH1]
	movaps xmm2, [rsp + nb203_fizH1]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H2 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb203_fixH2]
	movaps xmm1, [rsp + nb203_fiyH2]
	movaps xmm2, [rsp + nb203_fizH2]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate m i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb203_fixM]
	movaps xmm1, [rsp + nb203_fiyM]
	movaps xmm2, [rsp + nb203_fizM]

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
	movss  xmm3, [rdi + rcx*4 + 36]
	movss  xmm4, [rdi + rcx*4 + 40]
	movss  xmm5, [rdi + rcx*4 + 44]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 36], xmm3
	movss  [rdi + rcx*4 + 40], xmm4
	movss  [rdi + rcx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb203_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb203_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb203_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb203_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb203_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb203_n], esi
        jmp .nb203_outer
.nb203_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb203_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203_end
        ;# non-zero, do one more workunit
        jmp   .nb203_threadloop
.nb203_end:

	mov eax, [rsp + nb203_nouter]
	mov ebx, [rsp + nb203_ninner]
	mov rcx, [rbp + nb203_outeriter]
	mov rdx, [rbp + nb203_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 800
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


	




.globl nb_kernel203nf_x86_64_sse
.globl _nb_kernel203nf_x86_64_sse
nb_kernel203nf_x86_64_sse:	
_nb_kernel203nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb203nf_fshift,         16
.equiv          nb203nf_gid,            24
.equiv          nb203nf_pos,            32
.equiv          nb203nf_faction,        40
.equiv          nb203nf_charge,         48
.equiv          nb203nf_p_facel,        56
.equiv          nb203nf_argkrf,         64
.equiv          nb203nf_argcrf,         72
.equiv          nb203nf_Vc,             80
.equiv          nb203nf_type,           88
.equiv          nb203nf_p_ntype,        96
.equiv          nb203nf_vdwparam,       104
.equiv          nb203nf_Vvdw,           112
.equiv          nb203nf_p_tabscale,     120
.equiv          nb203nf_VFtab,          128
.equiv          nb203nf_invsqrta,       136
.equiv          nb203nf_dvda,           144
.equiv          nb203nf_p_gbtabscale,   152
.equiv          nb203nf_GBtab,          160
.equiv          nb203nf_p_nthreads,     168
.equiv          nb203nf_count,          176
.equiv          nb203nf_mtx,            184
.equiv          nb203nf_outeriter,      192
.equiv          nb203nf_inneriter,      200
.equiv          nb203nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb203nf_ixH1,           0
.equiv          nb203nf_iyH1,           16
.equiv          nb203nf_izH1,           32
.equiv          nb203nf_ixH2,           48
.equiv          nb203nf_iyH2,           64
.equiv          nb203nf_izH2,           80
.equiv          nb203nf_ixM,            96
.equiv          nb203nf_iyM,            112
.equiv          nb203nf_izM,            128
.equiv          nb203nf_iqH,            144
.equiv          nb203nf_iqM,            160
.equiv          nb203nf_qqH,            176
.equiv          nb203nf_qqM,            192
.equiv          nb203nf_vctot,          208
.equiv          nb203nf_half,           224
.equiv          nb203nf_three,          240
.equiv          nb203nf_krf,            256
.equiv          nb203nf_crf,            272
.equiv          nb203nf_krsqH1,         288
.equiv          nb203nf_krsqH2,         304
.equiv          nb203nf_krsqM,          320
.equiv          nb203nf_is3,            336
.equiv          nb203nf_ii3,            340
.equiv          nb203nf_innerjjnr,      344
.equiv          nb203nf_nri,            352
.equiv          nb203nf_iinr,           360
.equiv          nb203nf_jindex,         368
.equiv          nb203nf_jjnr,           376
.equiv          nb203nf_shift,          384
.equiv          nb203nf_shiftvec,       392
.equiv          nb203nf_facel,          400
.equiv          nb203nf_innerk,         408
.equiv          nb203nf_n,              412
.equiv          nb203nf_nn1,            416
.equiv          nb203nf_nouter,         420
.equiv          nb203nf_ninner,         424

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
	mov [rsp + nb203nf_nouter], eax
	mov [rsp + nb203nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb203nf_nri], edi
	mov [rsp + nb203nf_iinr], rsi
	mov [rsp + nb203nf_jindex], rdx
	mov [rsp + nb203nf_jjnr], rcx
	mov [rsp + nb203nf_shift], r8
	mov [rsp + nb203nf_shiftvec], r9
	mov rsi, [rbp + nb203nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb203nf_facel], xmm0

	mov rsi, [rbp + nb203nf_argkrf]
	mov rdi, [rbp + nb203nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb203nf_krf], xmm1
	movaps [rsp + nb203nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb203nf_half], eax
	movss xmm1, [rsp + nb203nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb203nf_half],  xmm1
	movaps [rsp + nb203nf_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb203nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb203nf_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb203nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb203nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb203nf_iqH], xmm3
	movaps [rsp + nb203nf_iqM], xmm4
	
.nb203nf_threadloop:
        mov   rsi, [rbp + nb203nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb203nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   rbx, 1                           ;# rbx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb203nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb203nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb203nf_n], eax
        mov [rsp + nb203nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb203nf_outerstart
        jmp .nb203nf_end
	
.nb203nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb203nf_nouter]
	mov [rsp + nb203nf_nouter], ebx

.nb203nf_outer:
	mov   rax, [rsp + nb203nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb203nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb203nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb203nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb203nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb203nf_ii3], ebx

	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb203nf_ixH1], xmm3
	movaps [rsp + nb203nf_iyH1], xmm4
	movaps [rsp + nb203nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb203nf_ixH2], xmm0
	movaps [rsp + nb203nf_iyH2], xmm1
	movaps [rsp + nb203nf_izH2], xmm2
	movaps [rsp + nb203nf_ixM], xmm3
	movaps [rsp + nb203nf_iyM], xmm4
	movaps [rsp + nb203nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb203nf_vctot], xmm4
	
	mov   rax, [rsp + nb203nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb203nf_pos]
	mov   rdi, [rbp + nb203nf_faction]	
	mov   rax, [rsp + nb203nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb203nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb203nf_ninner]
	mov   [rsp + nb203nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb203nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb203nf_unroll_loop
	jmp   .nb203nf_odd_inner
.nb203nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb203nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb203nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb203nf_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb203nf_iqH]
	mulps  xmm4, [rsp + nb203nf_iqM]

	movaps  [rsp + nb203nf_qqH], xmm3
	movaps  [rsp + nb203nf_qqM], xmm4

	mov rsi, [rbp + nb203nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb203nf_ixH1]
	movaps xmm5, [rsp + nb203nf_iyH1]
	movaps xmm6, [rsp + nb203nf_izH1]

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
	;# rsqH1 in xmm7 

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [rsp + nb203nf_ixH2]
	movaps xmm5, [rsp + nb203nf_iyH2]
	movaps xmm6, [rsp + nb203nf_izH2]

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
	;# rsqH2 in xmm6 

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [rsp + nb203nf_ixM]
	movaps xmm4, [rsp + nb203nf_iyM]
	movaps xmm5, [rsp + nb203nf_izM]

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

	movaps xmm0, xmm5
	movaps xmm1, xmm6
	movaps xmm2, xmm7

	mulps  xmm0, [rsp + nb203nf_krf]	
	mulps  xmm1, [rsp + nb203nf_krf]	
	mulps  xmm2, [rsp + nb203nf_krf]	

	movaps [rsp + nb203nf_krsqM], xmm0
	movaps [rsp + nb203nf_krsqH2], xmm1
	movaps [rsp + nb203nf_krsqH1], xmm2
	
	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb203nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb203nf_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb203nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb203nf_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb203nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb203nf_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions - xmm7=rinv
	addps xmm7, [rsp + nb203nf_krsqH1]
	subps xmm7, [rsp + nb203nf_crf] ;# xmm7=rinv+ krsq-crf 
	mulps xmm7, [rsp + nb203nf_qqH]
	addps xmm7, [rsp + nb203nf_vctot]

	;# H2 interactions - xmm6=rinv
	addps xmm6, [rsp + nb203nf_krsqH2]
	subps xmm6, [rsp + nb203nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps xmm6, [rsp + nb203nf_qqH]
	addps xmm6, xmm7

	;# M interactions - xmm5=rinv
	addps xmm5, [rsp + nb203nf_krsqM]
	subps xmm5, [rsp + nb203nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps xmm5, [rsp + nb203nf_qqM]
	addps xmm5, xmm6
	movaps [rsp + nb203nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb203nf_innerk],  4
	jl    .nb203nf_odd_inner
	jmp   .nb203nf_unroll_loop
.nb203nf_odd_inner:	
	add dword ptr [rsp + nb203nf_innerk],  4
	jnz   .nb203nf_odd_loop
	jmp   .nb203nf_updateouterdata
.nb203nf_odd_loop:
	mov   rdx, [rsp + nb203nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb203nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb203nf_iqM]
	mov rsi, [rbp + nb203nf_charge] 
	movhps xmm4, [rsp + nb203nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb203nf_qqM], xmm3	;# use dummy qq for storage 

	mov rsi, [rbp + nb203nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb203nf_ixM]
	movss xmm4, [rsp + nb203nf_iyM]
	movss xmm5, [rsp + nb203nf_izM]
	
	movlps xmm6, [rsp + nb203nf_ixH1]
	movlps xmm7, [rsp + nb203nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb203nf_iyH1]
	movlps xmm7, [rsp + nb203nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb203nf_izH1]
	movlps xmm7, [rsp + nb203nf_izH2]
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
	mulps xmm0, [rsp + nb203nf_krf]
	movaps [rsp + nb203nf_krsqM], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb203nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb203nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000

	;# xmm0=rinv 
	addps  xmm0, [rsp + nb203nf_krsqM]
	subps  xmm0, [rsp + nb203nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [rsp + nb203nf_qqM]	;# xmm0=vcoul 
	addps  xmm0, [rsp + nb203nf_vctot]
	movaps [rsp + nb203nf_vctot], xmm0
	
	dec dword ptr [rsp + nb203nf_innerk]
	jz    .nb203nf_updateouterdata
	jmp   .nb203nf_odd_loop
.nb203nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb203nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb203nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb203nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb203nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb203nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb203nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb203nf_n], esi
        jmp .nb203nf_outer
.nb203nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb203nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb203nf_end
        ;# non-zero, do one more workunit
        jmp   .nb203nf_threadloop
.nb203nf_end:

	mov eax, [rsp + nb203nf_nouter]
	mov ebx, [rsp + nb203nf_ninner]
	mov rcx, [rbp + nb203nf_outeriter]
	mov rdx, [rbp + nb203nf_inneriter]
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
