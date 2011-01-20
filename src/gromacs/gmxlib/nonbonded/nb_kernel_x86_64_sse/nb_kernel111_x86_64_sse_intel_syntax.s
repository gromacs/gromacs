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


.globl nb_kernel111_x86_64_sse
.globl _nb_kernel111_x86_64_sse
nb_kernel111_x86_64_sse:	
_nb_kernel111_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb111_fshift,           16
.equiv          nb111_gid,              24
.equiv          nb111_pos,              32
.equiv          nb111_faction,          40
.equiv          nb111_charge,           48
.equiv          nb111_p_facel,          56
.equiv          nb111_argkrf,           64
.equiv          nb111_argcrf,           72
.equiv          nb111_Vc,               80
.equiv          nb111_type,             88
.equiv          nb111_p_ntype,          96
.equiv          nb111_vdwparam,         104
.equiv          nb111_Vvdw,             112
.equiv          nb111_p_tabscale,       120
.equiv          nb111_VFtab,            128
.equiv          nb111_invsqrta,         136
.equiv          nb111_dvda,             144
.equiv          nb111_p_gbtabscale,     152
.equiv          nb111_GBtab,            160
.equiv          nb111_p_nthreads,       168
.equiv          nb111_count,            176
.equiv          nb111_mtx,              184
.equiv          nb111_outeriter,        192
.equiv          nb111_inneriter,        200
.equiv          nb111_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb111_ixO,              0
.equiv          nb111_iyO,              16
.equiv          nb111_izO,              32
.equiv          nb111_ixH1,             48
.equiv          nb111_iyH1,             64
.equiv          nb111_izH1,             80
.equiv          nb111_ixH2,             96
.equiv          nb111_iyH2,             112
.equiv          nb111_izH2,             128
.equiv          nb111_iqO,              144
.equiv          nb111_iqH,              160
.equiv          nb111_dxO,              176
.equiv          nb111_dyO,              192
.equiv          nb111_dzO,              208
.equiv          nb111_dxH1,             224
.equiv          nb111_dyH1,             240
.equiv          nb111_dzH1,             256
.equiv          nb111_dxH2,             272
.equiv          nb111_dyH2,             288
.equiv          nb111_dzH2,             304
.equiv          nb111_qqO,              320
.equiv          nb111_qqH,              336
.equiv          nb111_c6,               352
.equiv          nb111_c12,              368
.equiv          nb111_six,              384
.equiv          nb111_twelve,           400
.equiv          nb111_vctot,            416
.equiv          nb111_Vvdwtot,          432
.equiv          nb111_fixO,             448
.equiv          nb111_fiyO,             464
.equiv          nb111_fizO,             480
.equiv          nb111_fixH1,            496
.equiv          nb111_fiyH1,            512
.equiv          nb111_fizH1,            528
.equiv          nb111_fixH2,            544
.equiv          nb111_fiyH2,            560
.equiv          nb111_fizH2,            576
.equiv          nb111_fjx,              592
.equiv          nb111_fjy,              608
.equiv          nb111_fjz,              624
.equiv          nb111_half,             640
.equiv          nb111_three,            656
.equiv          nb111_is3,              672
.equiv          nb111_ii3,              676
.equiv          nb111_nri,              692
.equiv          nb111_iinr,             700
.equiv          nb111_jindex,           708
.equiv          nb111_jjnr,             716
.equiv          nb111_shift,            724
.equiv          nb111_shiftvec,         732
.equiv          nb111_facel,            740
.equiv          nb111_innerjjnr,        748
.equiv          nb111_ntia,             756
.equiv          nb111_innerk,           760
.equiv          nb111_n,                764
.equiv          nb111_nn1,              768
.equiv          nb111_nouter,           772
.equiv          nb111_ninner,           776

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
	sub rsp, 784
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
	mov [rsp + nb111_nouter], eax
	mov [rsp + nb111_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb111_nri], edi
	mov [rsp + nb111_iinr], rsi
	mov [rsp + nb111_jindex], rdx
	mov [rsp + nb111_jjnr], rcx
	mov [rsp + nb111_shift], r8
	mov [rsp + nb111_shiftvec], r9
	mov rsi, [rbp + nb111_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb111_facel], xmm0


	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb111_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb111_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb111_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb111_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb111_iqO], xmm3
	movaps [rsp + nb111_iqH], xmm4
	
	mov   rdx, [rbp + nb111_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb111_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb111_ntia], ecx		

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb111_half], eax
	movss xmm1, [rsp + nb111_half]
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
	movaps [rsp + nb111_half],  xmm1
	movaps [rsp + nb111_three],  xmm3
	movaps [rsp + nb111_six],  xmm4
	movaps [rsp + nb111_twelve],  xmm5

.nb111_threadloop:
        mov   rsi, [rbp + nb111_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb111_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb111_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb111_n], eax
        mov [rsp + nb111_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111_outerstart
        jmp .nb111_end

.nb111_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb111_nouter]
	mov [rsp + nb111_nouter], ebx

.nb111_outer:
	mov   rax, [rsp + nb111_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb111_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb111_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb111_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb111_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb111_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb111_ixO], xmm3
	movaps [rsp + nb111_iyO], xmm4
	movaps [rsp + nb111_izO], xmm5

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
	movaps [rsp + nb111_ixH1], xmm0
	movaps [rsp + nb111_iyH1], xmm1
	movaps [rsp + nb111_izH1], xmm2
	movaps [rsp + nb111_ixH2], xmm3
	movaps [rsp + nb111_iyH2], xmm4
	movaps [rsp + nb111_izH2], xmm5
	
	;# clear potentials and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb111_vctot], xmm4
	movaps [rsp + nb111_Vvdwtot], xmm4
	movaps [rsp + nb111_fixO], xmm4
	movaps [rsp + nb111_fiyO], xmm4
	movaps [rsp + nb111_fizO], xmm4
	movaps [rsp + nb111_fixH1], xmm4
	movaps [rsp + nb111_fiyH1], xmm4
	movaps [rsp + nb111_fizH1], xmm4
	movaps [rsp + nb111_fixH2], xmm4
	movaps [rsp + nb111_fiyH2], xmm4
	movaps [rsp + nb111_fizH2], xmm4
	
	mov   rax, [rsp + nb111_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb111_pos]
	mov   rdi, [rbp + nb111_faction]	
	mov   rax, [rsp + nb111_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb111_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb111_ninner]
	mov   [rsp + nb111_ninner], ecx
	add   edx, 0
	mov   [rsp + nb111_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111_unroll_loop
	jmp   .nb111_odd_inner
.nb111_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb111_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb111_pos]       ;# base of pos[] 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]     ;# replace jnr with j3 
	lea   rdx, [r11 + r11*2]	

	;# move four j coordinates to xmm0-xmm2 	
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
	mov rsi, [rbp + nb111_charge]    ;# base of charge[] 
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2

    subps xmm0, [rsp + nb111_ixO]
    subps xmm1, [rsp + nb111_iyO]
    subps xmm2, [rsp + nb111_izO]
    subps xmm3, [rsp + nb111_ixH1]
    subps xmm4, [rsp + nb111_iyH1]
    subps xmm5, [rsp + nb111_izH1]
    subps xmm6, [rsp + nb111_ixH2]
    subps xmm7, [rsp + nb111_iyH2]
    subps xmm8, [rsp + nb111_izH2]
    
	movaps [rsp + nb111_dxO], xmm0
	movaps [rsp + nb111_dyO], xmm1
	movaps [rsp + nb111_dzO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb111_dxH1], xmm3
	movaps [rsp + nb111_dyH1], xmm4
	movaps [rsp + nb111_dzH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb111_dxH2], xmm6
	movaps [rsp + nb111_dyH2], xmm7
	movaps [rsp + nb111_dzH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	movss xmm9, [rsi + r8*4]
	movss xmm10, [rsi + r10*4]
	movss xmm11, [rsi + r9*4]
	movss xmm12, [rsi + r11*4]

    unpcklps xmm9, xmm10
    unpcklps xmm11, xmm12
    unpcklps xmm9, xmm11
	movaps xmm7, xmm9
	mulps  xmm9, [rsp + nb111_iqO]
	mulps  xmm7, [rsp + nb111_iqH]

	movaps  [rsp + nb111_qqO], xmm9
	movaps  [rsp + nb111_qqH], xmm7

	mov rsi, [rbp + nb111_type]
	mov r8d, [rsi + r8*4]
	mov r9d, [rsi + r9*4]
	mov r10d, [rsi + r10*4]
	mov r11d, [rsi + r11*4]

	;# start doing invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	mov rsi, [rbp + nb111_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	

	movaps  xmm9, [rsp + nb111_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mov edi, [rsp + nb111_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm0, [rsp + nb111_half]
	mulps   xmm9, xmm0  ;# rinvO 
	mulps   xmm10, xmm0 ;# rinvH1
    mulps   xmm11, xmm0 ;# rinvH2
	
    movlps xmm3, [rsi + r8*4]
	movlps xmm2, [rsi + r10*4]
	movhps xmm3, [rsi + r9*4]
	movhps xmm2, [rsi + r11*4]

	movaps xmm5, xmm3
	shufps xmm3, xmm2, 136  ;# 10001000
	shufps xmm5, xmm2, 221  ;# 11011101

	;# interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9    ;# rinvsq
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    movaps xmm12, xmm9
    mulps  xmm12, xmm12 ;# rinv4
    mulps  xmm12, xmm9  ;# rinv6
    mulps  xmm0, [rsp + nb111_qqO] 
    mulps  xmm1, [rsp + nb111_qqH] 
    mulps  xmm2, [rsp + nb111_qqH] 
    movaps xmm13, xmm12 ;# rinv6
    mulps xmm12, xmm12 ;# rinv12
	mulps  xmm13, xmm3
	mulps  xmm12, xmm5
    movaps xmm14, xmm12
    subps  xmm14, xmm13
    
	addps  xmm14, [rsp + nb111_Vvdwtot]
	mulps  xmm13, [rsp + nb111_six]
	mulps  xmm12, [rsp + nb111_twelve]
	movaps [rsp + nb111_Vvdwtot], xmm14
    subps  xmm12, xmm13 ;# LJ fscal        

    addps  xmm12, xmm0
    
    mulps  xmm9, xmm12
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb111_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb111_vctot], xmm0
    
	mov   rdi, [rbp + nb111_faction]	
	;# move j  forces to local temp variables 
    movlps xmm0, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm1, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm0, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm1, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm2, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm3, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm4, [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5, [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm2, xmm4 ;# jza  -  jzb  -
    movlhps xmm3, xmm5 ;# jzc  -  jzd -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm0: jxa jya jxb jyb 
    ;# xmm1: jxc jyc jxd jyd
    ;# xmm2: jza jzb jzc jzd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb111_dxO]
	mulps xmm8, [rsp + nb111_dyO]
	mulps xmm9, [rsp + nb111_dzO]
	mulps xmm10, [rsp + nb111_dxH1]
	mulps xmm11, [rsp + nb111_dyH1]
	mulps xmm12, [rsp + nb111_dzH1]
	mulps xmm13, [rsp + nb111_dxH2]
	mulps xmm14, [rsp + nb111_dyH2]
	mulps xmm15, [rsp + nb111_dzH2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb111_fixO]
    addps xmm8, [rsp + nb111_fiyO]
    addps xmm9, [rsp + nb111_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb111_fixH1]
    addps xmm11, [rsp + nb111_fiyH1]
    addps xmm12, [rsp + nb111_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb111_fixH2]
    addps xmm14, [rsp + nb111_fiyH2]
    addps xmm15, [rsp + nb111_fizH2]

    movaps [rsp + nb111_fixO], xmm7
    movaps [rsp + nb111_fiyO], xmm8
    movaps [rsp + nb111_fizO], xmm9
    movaps [rsp + nb111_fixH1], xmm10
    movaps [rsp + nb111_fiyH1], xmm11
    movaps [rsp + nb111_fizH1], xmm12
    movaps [rsp + nb111_fixH2], xmm13
    movaps [rsp + nb111_fiyH2], xmm14
    movaps [rsp + nb111_fizH2], xmm15

    ;# xmm3 = fjx , xmm4 = fjy, xmm5=fjz
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2
    
    movlps [rdi + rax*4], xmm0
    movhps [rdi + rbx*4], xmm0
    movlps [rdi + rcx*4], xmm1
    movhps [rdi + rdx*4], xmm1
    movss  [rdi + rax*4 + 8], xmm2
    movss  [rdi + rcx*4 + 8], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 8], xmm2
    movss  [rdi + rdx*4 + 8], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb111_innerk],  4
	jl    .nb111_odd_inner
	jmp   .nb111_unroll_loop
.nb111_odd_inner:	
	add dword ptr [rsp + nb111_innerk],  4
	jnz   .nb111_odd_loop
	jmp   .nb111_updateouterdata
.nb111_odd_loop:
	mov   rdx, [rsp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb111_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb111_iqO]
	mov rsi, [rbp + nb111_charge] 
	movhps xmm4, [rsp + nb111_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb111_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb111_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb111_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb111_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb111_c6], xmm6
	movaps [rsp + nb111_c12], xmm7

	mov rsi, [rbp + nb111_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb111_ixO]
	movss xmm1, [rsp + nb111_iyO]
	movss xmm2, [rsp + nb111_izO]
	
	movlps xmm6, [rsp + nb111_ixH1]
	movlps xmm7, [rsp + nb111_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb111_iyH1]
	movlps xmm7, [rsp + nb111_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb111_izH1]
	movlps xmm7, [rsp + nb111_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb111_dxO], xmm3
	movaps [rsp + nb111_dyO], xmm4
	movaps [rsp + nb111_dzO], xmm5

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
	movaps xmm1, [rsp + nb111_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb111_half]
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
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	movaps xmm3, [rsp + nb111_qqO]
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [rsp + nb111_c6]
	mulps  xmm2, [rsp + nb111_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb111_Vvdwtot]
	mulss  xmm1, [rsp + nb111_six]
	mulss  xmm2, [rsp + nb111_twelve]
	subss  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm3, [rsp + nb111_vctot]

	movaps xmm0, [rsp + nb111_dxO]
	movaps xmm1, [rsp + nb111_dyO]
	movaps xmm2, [rsp + nb111_dzO]

	movaps [rsp + nb111_vctot], xmm3
	movaps [rsp + nb111_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb111_fixO]	
	movss  xmm4, [rsp + nb111_fiyO]	
	movss  xmm5, [rsp + nb111_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb111_fixO], xmm3	
	movss  [rsp + nb111_fiyO], xmm4	
	movss  [rsp + nb111_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb111_fixH1]
	addss  xmm4, [rsp + nb111_fiyH1]
	addss  xmm5, [rsp + nb111_fizH1]
	movss  [rsp + nb111_fixH1], xmm3	
	movss  [rsp + nb111_fiyH1], xmm4	
	movss  [rsp + nb111_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb111_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb111_fixH2]
	addss  xmm4, [rsp + nb111_fiyH2]
	addss  xmm5, [rsp + nb111_fizH2]
	movss  [rsp + nb111_fixH2], xmm3	
	movss  [rsp + nb111_fiyH2], xmm4	
	movss  [rsp + nb111_fizH2], xmm5	;# updated the H2 force 

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
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb111_innerk]
	jz    .nb111_updateouterdata
	jmp   .nb111_odd_loop
.nb111_updateouterdata:
	mov   ecx, [rsp + nb111_ii3]
	mov   rdi, [rbp + nb111_faction]
	mov   rsi, [rbp + nb111_fshift]
	mov   edx, [rsp + nb111_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb111_fixO]
	movaps xmm1, [rsp + nb111_fiyO]
	movaps xmm2, [rsp + nb111_fizO]

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
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb111_fixH1]
	movaps xmm1, [rsp + nb111_fiyH1]
	movaps xmm2, [rsp + nb111_fizH1]

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
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb111_fixH2]
	movaps xmm1, [rsp + nb111_fiyH2]
	movaps xmm2, [rsp + nb111_fizH2]

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

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb111_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb111_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb111_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb111_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb111_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb111_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb111_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb111_n], esi
        jmp .nb111_outer
.nb111_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb111_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111_end
        ;# non-zero, do one more workunit
        jmp   .nb111_threadloop
.nb111_end:


	mov eax, [rsp + nb111_nouter]
	mov ebx, [rsp + nb111_ninner]
	mov rcx, [rbp + nb111_outeriter]
	mov rdx, [rbp + nb111_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 784
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



.globl nb_kernel111nf_x86_64_sse
.globl _nb_kernel111nf_x86_64_sse
nb_kernel111nf_x86_64_sse:	
_nb_kernel111nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb111nf_fshift,         16
.equiv          nb111nf_gid,            24
.equiv          nb111nf_pos,            32
.equiv          nb111nf_faction,        40
.equiv          nb111nf_charge,         48
.equiv          nb111nf_p_facel,        56
.equiv          nb111nf_argkrf,         64
.equiv          nb111nf_argcrf,         72
.equiv          nb111nf_Vc,             80
.equiv          nb111nf_type,           88
.equiv          nb111nf_p_ntype,        96
.equiv          nb111nf_vdwparam,       104
.equiv          nb111nf_Vvdw,           112
.equiv          nb111nf_p_tabscale,     120
.equiv          nb111nf_VFtab,          128
.equiv          nb111nf_invsqrta,       136
.equiv          nb111nf_dvda,           144
.equiv          nb111nf_p_gbtabscale,   152
.equiv          nb111nf_GBtab,          160
.equiv          nb111nf_p_nthreads,     168
.equiv          nb111nf_count,          176
.equiv          nb111nf_mtx,            184
.equiv          nb111nf_outeriter,      192
.equiv          nb111nf_inneriter,      200
.equiv          nb111nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb111nf_ixO,            0
.equiv          nb111nf_iyO,            16
.equiv          nb111nf_izO,            32
.equiv          nb111nf_ixH1,           48
.equiv          nb111nf_iyH1,           64
.equiv          nb111nf_izH1,           80
.equiv          nb111nf_ixH2,           96
.equiv          nb111nf_iyH2,           112
.equiv          nb111nf_izH2,           128
.equiv          nb111nf_iqO,            144
.equiv          nb111nf_iqH,            160
.equiv          nb111nf_qqO,            176
.equiv          nb111nf_qqH,            192
.equiv          nb111nf_c6,             208
.equiv          nb111nf_c12,            224
.equiv          nb111nf_vctot,          240
.equiv          nb111nf_Vvdwtot,        256
.equiv          nb111nf_half,           272
.equiv          nb111nf_three,          288
.equiv          nb111nf_is3,            304
.equiv          nb111nf_ii3,            308
.equiv          nb111nf_nri,            324
.equiv          nb111nf_iinr,           332
.equiv          nb111nf_jindex,         340
.equiv          nb111nf_jjnr,           348
.equiv          nb111nf_shift,          356
.equiv          nb111nf_shiftvec,       364
.equiv          nb111nf_facel,          372
.equiv          nb111nf_innerjjnr,      380
.equiv          nb111nf_ntia,           388
.equiv          nb111nf_innerk,         392
.equiv          nb111nf_n,              396
.equiv          nb111nf_nn1,            400
.equiv          nb111nf_nouter,         404
.equiv          nb111nf_ninner,         408

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
	sub rsp, 416
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
	mov [rsp + nb111nf_nouter], eax
	mov [rsp + nb111nf_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb111nf_nri], edi
	mov [rsp + nb111nf_iinr], rsi
	mov [rsp + nb111nf_jindex], rdx
	mov [rsp + nb111nf_jjnr], rcx
	mov [rsp + nb111nf_shift], r8
	mov [rsp + nb111nf_shiftvec], r9
	mov rsi, [rbp + nb111nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb111nf_facel], xmm0

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb111nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb111nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb111nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb111nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb111nf_iqO], xmm3
	movaps [rsp + nb111nf_iqH], xmm4
	
	mov   rdx, [rbp + nb111nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb111nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb111nf_ntia], ecx		

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb111nf_half], eax
	movss xmm1, [rsp + nb111nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb111nf_half],  xmm1
	movaps [rsp + nb111nf_three],  xmm3

.nb111nf_threadloop:
        mov   rsi, [rbp + nb111nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb111nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb111nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb111nf_n], eax
        mov [rsp + nb111nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111nf_outerstart
        jmp .nb111nf_end

.nb111nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb111nf_nouter]
	mov [rsp + nb111nf_nouter], ebx

.nb111nf_outer:
	mov   rax, [rsp + nb111nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb111nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb111nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb111nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb111nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb111nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb111nf_ixO], xmm3
	movaps [rsp + nb111nf_iyO], xmm4
	movaps [rsp + nb111nf_izO], xmm5

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
	movaps [rsp + nb111nf_ixH1], xmm0
	movaps [rsp + nb111nf_iyH1], xmm1
	movaps [rsp + nb111nf_izH1], xmm2
	movaps [rsp + nb111nf_ixH2], xmm3
	movaps [rsp + nb111nf_iyH2], xmm4
	movaps [rsp + nb111nf_izH2], xmm5
	
	;# clear vctot and Vvdwtot
	xorps xmm4, xmm4
	movaps [rsp + nb111nf_vctot], xmm4
	movaps [rsp + nb111nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb111nf_jindex]
	mov   ecx, [rax + rsi*4]             ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb111nf_pos]
	mov   rax, [rsp + nb111nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb111nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb111nf_ninner]
	mov   [rsp + nb111nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb111nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111nf_unroll_loop
	jmp   .nb111nf_odd_inner
.nb111nf_unroll_loop:
        
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb111nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb111nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb111nf_iqO]
	mulps  xmm4, [rsp + nb111nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb111nf_qqO], xmm3
	movaps  [rsp + nb111nf_qqH], xmm4
	
	mov rsi, [rbp + nb111nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb111nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb111nf_ntia]
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

	movaps [rsp + nb111nf_c6], xmm4
	movaps [rsp + nb111nf_c12], xmm6

	mov rsi, [rbp + nb111nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
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

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb111nf_ixO]
	movaps xmm5, [rsp + nb111nf_iyO]
	movaps xmm6, [rsp + nb111nf_izO]

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
	movaps xmm4, [rsp + nb111nf_ixH1]
	movaps xmm5, [rsp + nb111nf_iyH1]
	movaps xmm6, [rsp + nb111nf_izH1]

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
	movaps xmm3, [rsp + nb111nf_ixH2]
	movaps xmm4, [rsp + nb111nf_iyH2]
	movaps xmm5, [rsp + nb111nf_izH2]

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

	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb111nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb111nf_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb111nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb111nf_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb111nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb111nf_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movaps  xmm4, xmm7	
	mulps   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm7, [rsp + nb111nf_qqO]	;# xmm7=vcoul 
	
	mulps  xmm1, [rsp + nb111nf_c6]
	mulps  xmm2, [rsp + nb111nf_c12]
	movaps xmm3, xmm2
	subps  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addps  xmm3, [rsp + nb111nf_Vvdwtot]
	addps  xmm7, [rsp + nb111nf_vctot]	
	movaps [rsp + nb111nf_Vvdwtot], xmm3
	movaps [rsp + nb111nf_vctot], xmm7

	;# H1 & H2 interactions 
	addps  xmm6, xmm5           ;# add H2 rinv 
	mulps  xmm6, [rsp + nb111nf_qqH]	;# xmm6=vcoul 
	addps  xmm6, [rsp + nb111nf_vctot]
	movaps [rsp + nb111nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb111nf_innerk],  4
	jl    .nb111nf_odd_inner
	jmp   .nb111nf_unroll_loop
.nb111nf_odd_inner:	
	add dword ptr [rsp + nb111nf_innerk],  4
	jnz   .nb111nf_odd_loop
	jmp   .nb111nf_updateouterdata
.nb111nf_odd_loop:
	mov   rdx, [rsp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb111nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb111nf_iqO]
	mov rsi, [rbp + nb111nf_charge] 
	movhps xmm4, [rsp + nb111nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb111nf_qqO], xmm3	;# use oxygen qq for storage 

	xorps xmm6, xmm6
	mov rsi, [rbp + nb111nf_type]
	mov ebx, [rsi + rax*4]
	mov rsi, [rbp + nb111nf_vdwparam]
	shl ebx, 1	
	add ebx, [rsp + nb111nf_ntia]
	movlps xmm6, [rsi + rbx*4]
	movaps xmm7, xmm6
	shufps xmm6, xmm6, 252  ;# 11111100
	shufps xmm7, xmm7, 253  ;# 11111101
	movaps [rsp + nb111nf_c6], xmm6
	movaps [rsp + nb111nf_c12], xmm7

	mov rsi, [rbp + nb111nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb111nf_ixO]
	movss xmm4, [rsp + nb111nf_iyO]
	movss xmm5, [rsp + nb111nf_izO]
	
	movlps xmm6, [rsp + nb111nf_ixH1]
	movlps xmm7, [rsp + nb111nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb111nf_iyH1]
	movlps xmm7, [rsp + nb111nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb111nf_izH1]
	movlps xmm7, [rsp + nb111nf_izH2]
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
	movaps xmm1, [rsp + nb111nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb111nf_half]
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
	movaps xmm1, xmm4
	mulss  xmm1, xmm4
	movaps xmm3, [rsp + nb111nf_qqO]
	mulss  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulss  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [rsp + nb111nf_c6]
	mulps  xmm2, [rsp + nb111nf_c12]
	movaps xmm5, xmm2
	subss  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [rsp + nb111nf_Vvdwtot]
	addps  xmm3, [rsp + nb111nf_vctot]
	movaps [rsp + nb111nf_vctot], xmm3
	movaps [rsp + nb111nf_Vvdwtot], xmm5

	dec dword ptr [rsp + nb111nf_innerk]
	jz    .nb111nf_updateouterdata
	jmp   .nb111nf_odd_loop
.nb111nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb111nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb111nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	movaps xmm7, [rsp + nb111nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb111nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb111nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb111nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb111nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb111nf_n], esi
        jmp .nb111nf_outer
.nb111nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb111nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111nf_end
        ;# non-zero, do one more workunit
        jmp   .nb111nf_threadloop
.nb111nf_end:


	mov eax, [rsp + nb111nf_nouter]
	mov ebx, [rsp + nb111nf_ninner]
	mov rcx, [rbp + nb111nf_outeriter]
	mov rdx, [rbp + nb111nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 416
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
