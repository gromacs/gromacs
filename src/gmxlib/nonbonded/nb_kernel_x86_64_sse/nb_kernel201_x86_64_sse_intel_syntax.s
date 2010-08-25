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




	

	

.globl nb_kernel201_x86_64_sse
.globl _nb_kernel201_x86_64_sse
nb_kernel201_x86_64_sse:	
_nb_kernel201_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb201_fshift,           16
.equiv          nb201_gid,              24
.equiv          nb201_pos,              32
.equiv          nb201_faction,          40
.equiv          nb201_charge,           48
.equiv          nb201_p_facel,          56
.equiv          nb201_argkrf,           64
.equiv          nb201_argcrf,           72
.equiv          nb201_Vc,               80
.equiv          nb201_type,             88
.equiv          nb201_p_ntype,          96
.equiv          nb201_vdwparam,         104
.equiv          nb201_Vvdw,             112
.equiv          nb201_p_tabscale,       120
.equiv          nb201_VFtab,            128
.equiv          nb201_invsqrta,         136
.equiv          nb201_dvda,             144
.equiv          nb201_p_gbtabscale,     152
.equiv          nb201_GBtab,            160
.equiv          nb201_p_nthreads,       168
.equiv          nb201_count,            176
.equiv          nb201_mtx,              184
.equiv          nb201_outeriter,        192
.equiv          nb201_inneriter,        200
.equiv          nb201_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb201_ixO,              0
.equiv          nb201_iyO,              16
.equiv          nb201_izO,              32
.equiv          nb201_ixH1,             48
.equiv          nb201_iyH1,             64
.equiv          nb201_izH1,             80
.equiv          nb201_ixH2,             96
.equiv          nb201_iyH2,             112
.equiv          nb201_izH2,             128
.equiv          nb201_iqO,              144
.equiv          nb201_iqH,              160
.equiv          nb201_dxO,              176
.equiv          nb201_dyO,              192
.equiv          nb201_dzO,              208
.equiv          nb201_dxH1,             224
.equiv          nb201_dyH1,             240
.equiv          nb201_dzH1,             256
.equiv          nb201_dxH2,             272
.equiv          nb201_dyH2,             288
.equiv          nb201_dzH2,             304
.equiv          nb201_qqO,              320
.equiv          nb201_qqH,              336
.equiv          nb201_vctot,            352
.equiv          nb201_fixO,             384
.equiv          nb201_fiyO,             400
.equiv          nb201_fizO,             416
.equiv          nb201_fixH1,            432
.equiv          nb201_fiyH1,            448
.equiv          nb201_fizH1,            464
.equiv          nb201_fixH2,            480
.equiv          nb201_fiyH2,            496
.equiv          nb201_fizH2,            512
.equiv          nb201_fjx,              528
.equiv          nb201_fjy,              544
.equiv          nb201_fjz,              560
.equiv          nb201_half,             576
.equiv          nb201_three,            592
.equiv          nb201_two,              608
.equiv          nb201_krf,              624
.equiv          nb201_crf,              640
.equiv          nb201_krsqO,            656
.equiv          nb201_krsqH1,           672
.equiv          nb201_krsqH2,           688
.equiv          nb201_nri,              704
.equiv          nb201_iinr,             712
.equiv          nb201_jindex,           720
.equiv          nb201_jjnr,             728
.equiv          nb201_shift,            736
.equiv          nb201_shiftvec,         744
.equiv          nb201_facel,            752
.equiv          nb201_innerjjnr,        760
.equiv          nb201_is3,              768
.equiv          nb201_ii3,              772
.equiv          nb201_innerk,           776
.equiv          nb201_n,                780
.equiv          nb201_nn1,              784
.equiv          nb201_nouter,           788
.equiv          nb201_ninner,           792

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
	mov [rsp + nb201_nouter], eax
	mov [rsp + nb201_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb201_nri], edi
	mov [rsp + nb201_iinr], rsi
	mov [rsp + nb201_jindex], rdx
	mov [rsp + nb201_jjnr], rcx
	mov [rsp + nb201_shift], r8
	mov [rsp + nb201_shiftvec], r9
	mov rsi, [rbp + nb201_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb201_facel], xmm0


	mov rsi, [rbp + nb201_argkrf]
	mov rdi, [rbp + nb201_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb201_krf], xmm1
	movaps [rsp + nb201_crf], xmm2

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb201_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb201_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb201_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb201_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb201_iqO], xmm3
	movaps [rsp + nb201_iqH], xmm4
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb201_half], eax
	movss xmm1, [rsp + nb201_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb201_half],  xmm1
	movaps [rsp + nb201_two],  xmm2
	movaps [rsp + nb201_three],  xmm3


.nb201_threadloop:
        mov   rsi, [rbp + nb201_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb201_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb201_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb201_n], eax
        mov [rsp + nb201_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201_outerstart
        jmp .nb201_end
	
.nb201_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb201_nouter]
	mov [rsp + nb201_nouter], ebx

.nb201_outer:
	mov   rax, [rsp + nb201_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb201_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb201_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb201_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb201_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb201_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb201_ixO], xmm3
	movaps [rsp + nb201_iyO], xmm4
	movaps [rsp + nb201_izO], xmm5

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
	movaps [rsp + nb201_ixH1], xmm0
	movaps [rsp + nb201_iyH1], xmm1
	movaps [rsp + nb201_izH1], xmm2
	movaps [rsp + nb201_ixH2], xmm3
	movaps [rsp + nb201_iyH2], xmm4
	movaps [rsp + nb201_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb201_vctot], xmm4
	movaps [rsp + nb201_fixO], xmm4
	movaps [rsp + nb201_fiyO], xmm4
	movaps [rsp + nb201_fizO], xmm4
	movaps [rsp + nb201_fixH1], xmm4
	movaps [rsp + nb201_fiyH1], xmm4
	movaps [rsp + nb201_fizH1], xmm4
	movaps [rsp + nb201_fixH2], xmm4
	movaps [rsp + nb201_fiyH2], xmm4
	movaps [rsp + nb201_fizH2], xmm4
	
	mov   rax, [rsp + nb201_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb201_pos]
	mov   rdi, [rbp + nb201_faction]	
	mov   rax, [rsp + nb201_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb201_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb201_ninner]
	mov   [rsp + nb201_ninner], ecx
	add   edx, 0
	mov   [rsp + nb201_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201_unroll_loop
	jmp   .nb201_odd_inner
.nb201_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb201_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb201_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb201_iqO]
	mulps  xmm4, [rsp + nb201_iqH]

	movaps  [rsp + nb201_qqO], xmm3
	movaps  [rsp + nb201_qqH], xmm4

	mov rsi, [rbp + nb201_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

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
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb201_ixO]
    subps xmm1, [rsp + nb201_iyO]
    subps xmm2, [rsp + nb201_izO]
    subps xmm3, [rsp + nb201_ixH1]
    subps xmm4, [rsp + nb201_iyH1]
    subps xmm5, [rsp + nb201_izH1]
    subps xmm6, [rsp + nb201_ixH2]
    subps xmm7, [rsp + nb201_iyH2]
    subps xmm8, [rsp + nb201_izH2]
    
	movaps [rsp + nb201_dxO], xmm0
	movaps [rsp + nb201_dyO], xmm1
	movaps [rsp + nb201_dzO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb201_dxH1], xmm3
	movaps [rsp + nb201_dyH1], xmm4
	movaps [rsp + nb201_dzH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb201_dxH2], xmm6
	movaps [rsp + nb201_dyH2], xmm7
	movaps [rsp + nb201_dzH2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

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
		
	movaps  xmm9, [rsp + nb201_three]
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

	movaps  xmm4, [rsp + nb201_half]
	mulps   xmm9, xmm4  ;# rinvO
	mulps   xmm10, xmm4 ;# rinvH1
    mulps   xmm11, xmm4 ;# rinvH2
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, xmm9 ;# copy of rinv
    movaps xmm4, xmm10
    movaps xmm7, xmm11
    movaps xmm2, [rsp + nb201_krf]    
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
    movaps xmm14, [rsp + nb201_crf]
    subps  xmm2, xmm14   ;# rinv+krsq-crf
    subps  xmm5, xmm14
    subps  xmm8, xmm14
    movaps xmm12, [rsp + nb201_qqO]
    movaps xmm13, [rsp + nb201_qqH]    
    mulps  xmm2, xmm12 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm5, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    mulps  xmm8, xmm13 ;# voul=qq*(rinv+ krsq-crf)
    addps  xmm0, xmm0 ;# 2*krsq
    addps  xmm3, xmm3 
    addps  xmm6, xmm6 
    subps  xmm1, xmm0 ;# rinv-2*krsq
    subps  xmm4, xmm3
    subps  xmm7, xmm6
    mulps  xmm1, xmm12   ;# (rinv-2*krsq)*qq
    mulps  xmm4, xmm13
    mulps  xmm7, xmm13
    addps  xmm2, [rsp + nb201_vctot]
    addps  xmm5, xmm8
    addps  xmm2, xmm5
    movaps [rsp + nb201_vctot], xmm2
    
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

	mulps xmm0, [rsp + nb201_dxO]
	mulps xmm1, [rsp + nb201_dyO]
	mulps xmm2, [rsp + nb201_dzO]
	mulps xmm3, [rsp + nb201_dxH1]
	mulps xmm4, [rsp + nb201_dyH1]
	mulps xmm5, [rsp + nb201_dzH1]
	mulps xmm6, [rsp + nb201_dxH2]
	mulps xmm7, [rsp + nb201_dyH2]
	mulps xmm8, [rsp + nb201_dzH2]

    movaps xmm13, xmm0
    movaps xmm14, xmm1
    addps xmm11, xmm2
    addps xmm0, [rsp + nb201_fixO]
    addps xmm1, [rsp + nb201_fiyO]
    addps xmm2, [rsp + nb201_fizO]

    addps xmm13,  xmm3
    addps xmm14, xmm4
    addps xmm11, xmm5
    addps xmm3, [rsp + nb201_fixH1]
    addps xmm4, [rsp + nb201_fiyH1]
    addps xmm5, [rsp + nb201_fizH1]

    addps xmm13, xmm6
    addps xmm14, xmm7
    addps xmm11, xmm8
    addps xmm6, [rsp + nb201_fixH2]
    addps xmm7, [rsp + nb201_fiyH2]
    addps xmm8, [rsp + nb201_fizH2]

    movaps [rsp + nb201_fixO], xmm0
    movaps [rsp + nb201_fiyO], xmm1
    movaps [rsp + nb201_fizO], xmm2
    movaps [rsp + nb201_fixH1], xmm3
    movaps [rsp + nb201_fiyH1], xmm4
    movaps [rsp + nb201_fizH1], xmm5
    movaps [rsp + nb201_fixH2], xmm6
    movaps [rsp + nb201_fiyH2], xmm7
    movaps [rsp + nb201_fizH2], xmm8
    
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
	sub dword ptr [rsp + nb201_innerk],  4
	jl    .nb201_odd_inner
	jmp   .nb201_unroll_loop
.nb201_odd_inner:	
	add dword ptr [rsp + nb201_innerk],  4
	jnz   .nb201_odd_loop
	jmp   .nb201_updateouterdata
.nb201_odd_loop:
	mov   rdx, [rsp + nb201_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb201_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb201_iqO]
	mov rsi, [rbp + nb201_charge] 
	movhps xmm4, [rsp + nb201_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb201_qqO], xmm3	;# use oxygen qq for storage 

	mov rsi, [rbp + nb201_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb201_ixO]
	movss xmm1, [rsp + nb201_iyO]
	movss xmm2, [rsp + nb201_izO]
	
	movlps xmm6, [rsp + nb201_ixH1]
	movlps xmm7, [rsp + nb201_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb201_iyH1]
	movlps xmm7, [rsp + nb201_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb201_izH1]
	movlps xmm7, [rsp + nb201_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb201_dxO], xmm3
	movaps [rsp + nb201_dyO], xmm4
	movaps [rsp + nb201_dzO], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	movaps xmm0, xmm4
	mulps xmm0, [rsp + nb201_krf]
	movaps [rsp + nb201_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb201_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb201_half]
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
	movaps xmm3, [rsp + nb201_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [rsp + nb201_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm3, [rsp + nb201_two]
	subps  xmm1, xmm3	;# xmm1=rinv-2*krsq 
	mulps  xmm0, [rsp + nb201_qqO]	;# xmm0=vcoul 
	mulps  xmm1, [rsp + nb201_qqO] 	;# xmm1=coul part of fs 

	
	mulps  xmm4, xmm1	;# xmm4=total fscal 
	addps  xmm0, [rsp + nb201_vctot]
	movaps [rsp + nb201_vctot], xmm0
	
	movaps xmm0, [rsp + nb201_dxO]
	movaps xmm1, [rsp + nb201_dyO]
	movaps xmm2, [rsp + nb201_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb201_fixO]	
	movss  xmm4, [rsp + nb201_fiyO]	
	movss  xmm5, [rsp + nb201_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb201_fixO], xmm3	
	movss  [rsp + nb201_fiyO], xmm4	
	movss  [rsp + nb201_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb201_fixH1]
	addss  xmm4, [rsp + nb201_fiyH1]
	addss  xmm5, [rsp + nb201_fizH1]
	movss  [rsp + nb201_fixH1], xmm3	
	movss  [rsp + nb201_fiyH1], xmm4	
	movss  [rsp + nb201_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb201_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb201_fixH2]
	addss  xmm4, [rsp + nb201_fiyH2]
	addss  xmm5, [rsp + nb201_fizH2]
	movss  [rsp + nb201_fixH2], xmm3	
	movss  [rsp + nb201_fiyH2], xmm4	
	movss  [rsp + nb201_fizH2], xmm5	;# updated the H2 force 

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

	dec dword ptr [rsp + nb201_innerk]
	jz    .nb201_updateouterdata
	jmp   .nb201_odd_loop
.nb201_updateouterdata:
	mov   ecx, [rsp + nb201_ii3]
	mov   rdi, [rbp + nb201_faction]
	mov   rsi, [rbp + nb201_fshift]
	mov   edx, [rsp + nb201_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb201_fixO]
	movaps xmm1, [rsp + nb201_fiyO]
	movaps xmm2, [rsp + nb201_fizO]

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
	movaps xmm0, [rsp + nb201_fixH1]
	movaps xmm1, [rsp + nb201_fiyH1]
	movaps xmm2, [rsp + nb201_fizH1]

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
	movaps xmm0, [rsp + nb201_fixH2]
	movaps xmm1, [rsp + nb201_fiyH2]
	movaps xmm2, [rsp + nb201_fizH2]

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
	mov esi, [rsp + nb201_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb201_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb201_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb201_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb201_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb201_n], esi
        jmp .nb201_outer
.nb201_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb201_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201_end
        ;# non-zero, do one more workunit
        jmp   .nb201_threadloop
.nb201_end:

	mov eax, [rsp + nb201_nouter]
	mov ebx, [rsp + nb201_ninner]
	mov rcx, [rbp + nb201_outeriter]
	mov rdx, [rbp + nb201_inneriter]
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

	

.globl nb_kernel201nf_x86_64_sse
.globl _nb_kernel201nf_x86_64_sse
nb_kernel201nf_x86_64_sse:	
_nb_kernel201nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb201nf_fshift,         16
.equiv          nb201nf_gid,            24
.equiv          nb201nf_pos,            32
.equiv          nb201nf_faction,        40
.equiv          nb201nf_charge,         48
.equiv          nb201nf_p_facel,        56
.equiv          nb201nf_argkrf,         64
.equiv          nb201nf_argcrf,         72
.equiv          nb201nf_Vc,             80
.equiv          nb201nf_type,           88
.equiv          nb201nf_p_ntype,        96
.equiv          nb201nf_vdwparam,       104
.equiv          nb201nf_Vvdw,           112
.equiv          nb201nf_p_tabscale,     120
.equiv          nb201nf_VFtab,          128
.equiv          nb201nf_invsqrta,       136
.equiv          nb201nf_dvda,           144
.equiv          nb201nf_p_gbtabscale,   152
.equiv          nb201nf_GBtab,          160
.equiv          nb201nf_p_nthreads,     168
.equiv          nb201nf_count,          176
.equiv          nb201nf_mtx,            184
.equiv          nb201nf_outeriter,      192
.equiv          nb201nf_inneriter,      200
.equiv          nb201nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb201nf_ixO,            0
.equiv          nb201nf_iyO,            16
.equiv          nb201nf_izO,            32
.equiv          nb201nf_ixH1,           48
.equiv          nb201nf_iyH1,           64
.equiv          nb201nf_izH1,           80
.equiv          nb201nf_ixH2,           96
.equiv          nb201nf_iyH2,           112
.equiv          nb201nf_izH2,           128
.equiv          nb201nf_iqO,            144
.equiv          nb201nf_iqH,            160
.equiv          nb201nf_qqO,            176
.equiv          nb201nf_qqH,            192
.equiv          nb201nf_vctot,          208
.equiv          nb201nf_half,           224
.equiv          nb201nf_three,          240
.equiv          nb201nf_krf,            256
.equiv          nb201nf_crf,            272
.equiv          nb201nf_krsqO,          288
.equiv          nb201nf_krsqH1,         304
.equiv          nb201nf_krsqH2,         320
.equiv          nb201nf_is3,            336
.equiv          nb201nf_ii3,            340
.equiv          nb201nf_innerjjnr,      344
.equiv          nb201nf_nri,            352
.equiv          nb201nf_iinr,           360
.equiv          nb201nf_jindex,         368
.equiv          nb201nf_jjnr,           376
.equiv          nb201nf_shift,          384
.equiv          nb201nf_shiftvec,       392
.equiv          nb201nf_facel,          400
.equiv          nb201nf_innerk,         408
.equiv          nb201nf_n,              412
.equiv          nb201nf_nn1,            416
.equiv          nb201nf_nouter,         420
.equiv          nb201nf_ninner,         424

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
	mov [rsp + nb201nf_nouter], eax
	mov [rsp + nb201nf_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb201nf_nri], edi
	mov [rsp + nb201nf_iinr], rsi
	mov [rsp + nb201nf_jindex], rdx
	mov [rsp + nb201nf_jjnr], rcx
	mov [rsp + nb201nf_shift], r8
	mov [rsp + nb201nf_shiftvec], r9
	mov rsi, [rbp + nb201nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb201nf_facel], xmm0

	mov rsi, [rbp + nb201nf_argkrf]
	mov rdi, [rbp + nb201nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb201nf_krf], xmm1
	movaps [rsp + nb201nf_crf], xmm2

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb201nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb201nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb201nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb201nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb201nf_iqO], xmm3
	movaps [rsp + nb201nf_iqH], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb201nf_half], eax
	movss xmm1, [rsp + nb201nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb201nf_half],  xmm1
	movaps [rsp + nb201nf_three],  xmm3

.nb201nf_threadloop:
        mov   rsi, [rbp + nb201nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb201nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb201nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb201nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb201nf_n], eax
        mov [rsp + nb201nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb201nf_outerstart
        jmp .nb201nf_end
	
.nb201nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb201nf_nouter]
	mov [rsp + nb201nf_nouter], ebx

.nb201nf_outer:
	mov   rax, [rsp + nb201nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb201nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb201nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb201nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb201nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb201nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb201nf_ixO], xmm3
	movaps [rsp + nb201nf_iyO], xmm4
	movaps [rsp + nb201nf_izO], xmm5

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
	movaps [rsp + nb201nf_ixH1], xmm0
	movaps [rsp + nb201nf_iyH1], xmm1
	movaps [rsp + nb201nf_izH1], xmm2
	movaps [rsp + nb201nf_ixH2], xmm3
	movaps [rsp + nb201nf_iyH2], xmm4
	movaps [rsp + nb201nf_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb201nf_vctot], xmm4
	
	mov   rax, [rsp + nb201nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb201nf_pos]
	mov   rax, [rsp + nb201nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb201nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb201nf_ninner]
	mov   [rsp + nb201nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb201nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb201nf_unroll_loop
	jmp   .nb201nf_odd_inner
.nb201nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb201nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb201nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb201nf_iqO]
	mulps  xmm4, [rsp + nb201nf_iqH]

	movaps  [rsp + nb201nf_qqO], xmm3
	movaps  [rsp + nb201nf_qqH], xmm4

	mov rsi, [rbp + nb201nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb201nf_ixO]
	movaps xmm5, [rsp + nb201nf_iyO]
	movaps xmm6, [rsp + nb201nf_izO]

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
	movaps xmm4, [rsp + nb201nf_ixH1]
	movaps xmm5, [rsp + nb201nf_iyH1]
	movaps xmm6, [rsp + nb201nf_izH1]

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
	movaps xmm3, [rsp + nb201nf_ixH2]
	movaps xmm4, [rsp + nb201nf_iyH2]
	movaps xmm5, [rsp + nb201nf_izH2]

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

	mulps  xmm0, [rsp + nb201nf_krf]	
	mulps  xmm1, [rsp + nb201nf_krf]	
	mulps  xmm2, [rsp + nb201nf_krf]	

	movaps [rsp + nb201nf_krsqH2], xmm0
	movaps [rsp + nb201nf_krsqH1], xmm1
	movaps [rsp + nb201nf_krsqO], xmm2
	
	;# start with rsqO - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb201nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb201nf_half]
	movaps  xmm7, xmm4	;# rinvO in xmm7 
	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb201nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb201nf_half]
	movaps  xmm6, xmm4	;# rinvH1 in xmm6 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb201nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb201nf_half]
	movaps  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 

	movaps xmm0, xmm7
	movaps xmm1, [rsp + nb201nf_krsqO]
	addps  xmm0, xmm1
	subps  xmm0, [rsp + nb201nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [rsp + nb201nf_qqO]

	addps  xmm0, [rsp + nb201nf_vctot]
	movaps [rsp + nb201nf_vctot], xmm0

	;# H1 interactions 
	movaps  xmm7, xmm6
	movaps  xmm0, [rsp + nb201nf_krsqH1]
	addps   xmm6, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm6, [rsp + nb201nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps   xmm6, [rsp + nb201nf_qqH] ;# vcoul 
	addps   xmm6, [rsp + nb201nf_vctot]
	
	;# H2 interactions 
	movaps  xmm7, xmm5
	movaps  xmm0, [rsp + nb201nf_krsqH2]
	addps   xmm5, xmm0	;# xmm6=rinv+ krsq 
	subps   xmm5, [rsp + nb201nf_crf] ;# xmm5=rinv+ krsq-crf 
	mulps   xmm5, [rsp + nb201nf_qqH] ;# vcoul 
	addps  xmm6, xmm5
	movaps [rsp + nb201nf_vctot], xmm6

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb201nf_innerk],  4
	jl    .nb201nf_odd_inner
	jmp   .nb201nf_unroll_loop
.nb201nf_odd_inner:	
	add dword ptr [rsp + nb201nf_innerk],  4
	jnz   .nb201nf_odd_loop
	jmp   .nb201nf_updateouterdata
.nb201nf_odd_loop:
	mov   rdx, [rsp + nb201nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb201nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb201nf_iqO]
	mov rsi, [rbp + nb201nf_charge] 
	movhps xmm4, [rsp + nb201nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb201nf_qqO], xmm3	;# use oxygen qq for storage 

	mov rsi, [rbp + nb201nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb201nf_ixO]
	movss xmm4, [rsp + nb201nf_iyO]
	movss xmm5, [rsp + nb201nf_izO]
	
	movlps xmm6, [rsp + nb201nf_ixH1]
	movlps xmm7, [rsp + nb201nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb201nf_iyH1]
	movlps xmm7, [rsp + nb201nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb201nf_izH1]
	movlps xmm7, [rsp + nb201nf_izH2]
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
	mulps xmm0, [rsp + nb201nf_krf]
	movaps [rsp + nb201nf_krsqO], xmm0
	
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb201nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb201nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	
	
	movaps xmm3, [rsp + nb201nf_krsqO]
	addps  xmm0, xmm3	;# xmm0=rinv+ krsq 
	subps  xmm0, [rsp + nb201nf_crf] ;# xmm0=rinv+ krsq-crf 
	mulps  xmm0, [rsp + nb201nf_qqO]	;# xmm0=vcoul 
	addps  xmm0, [rsp + nb201nf_vctot]
	movaps [rsp + nb201nf_vctot], xmm0

	dec dword ptr [rsp + nb201nf_innerk]
	jz    .nb201nf_updateouterdata
	jmp   .nb201nf_odd_loop
.nb201nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb201nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb201nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb201nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb201nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb201nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb201nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb201nf_n], esi
        jmp .nb201nf_outer
.nb201nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb201nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb201nf_end
        ;# non-zero, do one more workunit
        jmp   .nb201nf_threadloop
.nb201nf_end:

	mov eax, [rsp + nb201nf_nouter]
	mov ebx, [rsp + nb201nf_ninner]
	mov rcx, [rbp + nb201nf_outeriter]
	mov rdx, [rbp + nb201nf_inneriter]
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
