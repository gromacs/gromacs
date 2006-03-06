##
## $Id$
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##






.globl nb_kernel410_x86_64_sse
.globl _nb_kernel410_x86_64_sse
nb_kernel410_x86_64_sse:        
_nb_kernel410_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb410_fshift, 16
.set nb410_gid, 24
.set nb410_pos, 32
.set nb410_faction, 40
.set nb410_charge, 48
.set nb410_p_facel, 56
.set nb410_argkrf, 64
.set nb410_argcrf, 72
.set nb410_Vc, 80
.set nb410_type, 88
.set nb410_p_ntype, 96
.set nb410_vdwparam, 104
.set nb410_Vvdw, 112
.set nb410_p_tabscale, 120
.set nb410_VFtab, 128
.set nb410_invsqrta, 136
.set nb410_dvda, 144
.set nb410_p_gbtabscale, 152
.set nb410_GBtab, 160
.set nb410_p_nthreads, 168
.set nb410_count, 176
.set nb410_mtx, 184
.set nb410_outeriter, 192
.set nb410_inneriter, 200
.set nb410_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb410_ix, 0
.set nb410_iy, 16
.set nb410_iz, 32
.set nb410_iq, 48
.set nb410_dx, 64
.set nb410_dy, 80
.set nb410_dz, 96
.set nb410_two, 112
.set nb410_six, 128
.set nb410_twelve, 144
.set nb410_gbtsc, 160
.set nb410_qq, 176
.set nb410_c6, 192
.set nb410_c12, 208
.set nb410_fscal, 224
.set nb410_vctot, 240
.set nb410_Vvdwtot, 256
.set nb410_fix, 272
.set nb410_fiy, 288
.set nb410_fiz, 304
.set nb410_half, 320
.set nb410_three, 336
.set nb410_r, 352
.set nb410_isai, 368
.set nb410_isaprod, 384
.set nb410_dvdasum, 400
.set nb410_gbscale, 416
.set nb410_nri, 432
.set nb410_iinr, 440
.set nb410_jindex, 448
.set nb410_jjnr, 456
.set nb410_shift, 464
.set nb410_shiftvec, 472
.set nb410_facel, 480
.set nb410_innerjjnr, 488
.set nb410_is3, 496
.set nb410_ii3, 500
.set nb410_ii, 504
.set nb410_ntia, 508
.set nb410_innerk, 512
.set nb410_n, 516
.set nb410_nn1, 520
.set nb410_ntype, 524
.set nb410_nouter, 528
.set nb410_ninner, 532
.set nb410_jnra, 536
.set nb410_jnrb, 540
.set nb410_jnrc, 544
.set nb410_jnrd, 548

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $568,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb410_nouter(%rsp)
        movl %eax,nb410_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb410_nri(%rsp)
        movq %rsi,nb410_iinr(%rsp)
        movq %rdx,nb410_jindex(%rsp)
        movq %rcx,nb410_jjnr(%rsp)
        movq %r8,nb410_shift(%rsp)
        movq %r9,nb410_shiftvec(%rsp)
        movq nb410_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb410_ntype(%rsp)
        movq nb410_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb410_facel(%rsp)

        movq nb410_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb410_gbtsc(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb410_half(%rsp)
        movss nb410_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb410_half(%rsp)
        movaps %xmm2,nb410_two(%rsp)
        movaps %xmm3,nb410_three(%rsp)
        movaps %xmm4,nb410_six(%rsp)
        movaps %xmm5,nb410_twelve(%rsp)

_nb_kernel410_x86_64_sse.nb410_threadloop: 
        movq  nb410_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel410_x86_64_sse.nb410_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410_x86_64_sse.nb410_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb410_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb410_n(%rsp)
        movl %ebx,nb410_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel410_x86_64_sse.nb410_outerstart
        jmp _nb_kernel410_x86_64_sse.nb410_end

_nb_kernel410_x86_64_sse.nb410_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410_nouter(%rsp),%ebx
        movl %ebx,nb410_nouter(%rsp)

_nb_kernel410_x86_64_sse.nb410_outer: 
        movq  nb410_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb410_is3(%rsp)      ## store is3 

        movq  nb410_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb410_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 
        movl  %ebx,nb410_ii(%rsp)

        movq  nb410_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb410_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb410_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movq  nb410_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb410_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb410_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb410_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb410_iq(%rsp)
        movaps %xmm4,nb410_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb410_ix(%rsp)
        movaps %xmm1,nb410_iy(%rsp)
        movaps %xmm2,nb410_iz(%rsp)

        movl  %ebx,nb410_ii3(%rsp)

        ## clear vctot and i forces 
        xorps %xmm13,%xmm13
        movaps %xmm13,%xmm12
        movaps %xmm13,nb410_Vvdwtot(%rsp)
        movaps %xmm13,nb410_dvdasum(%rsp)
        movaps %xmm13,%xmm14
        movaps %xmm13,%xmm15

        movq  nb410_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb410_pos(%rbp),%rsi
        movq  nb410_faction(%rbp),%rdi
        movq  nb410_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb410_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb410_ninner(%rsp),%ecx
        movl  %ecx,nb410_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb410_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel410_x86_64_sse.nb410_unroll_loop
        jmp   _nb_kernel410_x86_64_sse.nb410_finish_inner
_nb_kernel410_x86_64_sse.nb410_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb410_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb410_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        ## load isaj
        movq nb410_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7
        movaps nb410_isai(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all isaj in xmm3 
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb410_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb410_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb410_gbscale(%rsp)

        movq nb410_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        mulps nb410_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb410_qq(%rsp)

    ## vdw parameters
        movq nb410_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        movl (%rsi,%rcx,4),%r14d
        movl (%rsi,%rdx,4),%r15d
        shll %r12d
        shll %r13d
        shll %r14d
        shll %r15d
    movl nb410_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d
        addl %edi,%r14d
        addl %edi,%r15d

        movq nb410_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm3
        movhps (%rsi,%r15,4),%xmm7

        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb410_c6(%rsp)
    movaps %xmm3,nb410_c12(%rsp)

        movq nb410_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## jnr 
        lea  (%rbx,%rbx,2),%r9
        lea  (%rcx,%rcx,2),%r10
        lea  (%rdx,%rdx,2),%r11

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm4
        movlps (%rsi,%r10,4),%xmm5
        movss 8(%rsi,%r8,4),%xmm2
        movss 8(%rsi,%r10,4),%xmm6

        movhps (%rsi,%r9,4),%xmm4
        movhps (%rsi,%r11,4),%xmm5

        movss 8(%rsi,%r9,4),%xmm0
        movss 8(%rsi,%r11,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## calc dr 
        subps nb410_ix(%rsp),%xmm0
        subps nb410_iy(%rsp),%xmm1
        subps nb410_iz(%rsp),%xmm2

        ## store dr 
        movaps %xmm0,nb410_dx(%rsp)
        movaps %xmm1,nb410_dy(%rsp)
        movaps %xmm2,nb410_dz(%rsp)

        ## square it 
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
    movaps %xmm0,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb410_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb410_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb410_r(%rsp)
        mulps nb410_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6

    ## multiply by 4
    pslld   $2,%xmm5

    ## move to integer registers
    movhlps %xmm5,%xmm7
    movd    %xmm5,%r12d
    movd    %xmm7,%r14d
    pshufd $1,%xmm5,%xmm5
    pshufd $1,%xmm7,%xmm7
    movd    %xmm5,%r13d
    movd    %xmm7,%r15d

    ## calculate eps
    subps     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb410_GBtab(%rbp),%rsi

    movaps %xmm0,%xmm9 ## rinv
    mulps  %xmm9,%xmm9 ## rinvsq
    movaps %xmm9,%xmm10 ## rinvsq
    mulps  %xmm10,%xmm10 ## rinv4
    mulps  %xmm9,%xmm10 ## rinv6
    movaps %xmm10,%xmm11
    mulps  %xmm11,%xmm11 ## rinv12

    ## load table data
        movlps (%rsi,%r12,4),%xmm5
        movlps (%rsi,%r14,4),%xmm7
        movhps (%rsi,%r13,4),%xmm5
        movhps (%rsi,%r15,4),%xmm7

    movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

    mulps  nb410_c6(%rsp),%xmm10      ## vvdw6=c6*rinv6
        mulps  nb410_c12(%rsp),%xmm11     ## vvdw12=c12*rinv12     

        movaps %xmm11,%xmm9
        subps  %xmm10,%xmm11    ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addps  nb410_Vvdwtot(%rsp),%xmm11
    movaps %xmm11,nb410_Vvdwtot(%rsp)

        movlps 8(%rsi,%r12,4),%xmm7
        movlps 8(%rsi,%r14,4),%xmm8
        movhps 8(%rsi,%r13,4),%xmm7
        movhps 8(%rsi,%r15,4),%xmm8

    movaps %xmm7,%xmm6

        shufps $136,%xmm8,%xmm6 ## 10001000
        shufps $221,%xmm8,%xmm7 ## 11011101
    ## table data ready in xmm4-xmm7

    mulps  %xmm1,%xmm7  ## Heps
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm1,%xmm7      ## Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        addps  %xmm7,%xmm7      ## two*Heps2 
        movaps nb410_qq(%rsp),%xmm3
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

    ## LJ forces
    mulps  nb410_six(%rsp),%xmm10
    mulps  nb410_twelve(%rsp),%xmm9
    subps  %xmm10,%xmm9
    mulps  %xmm0,%xmm9 ## (12*vnb12-6*vnb6)*rinv

        movq nb410_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb410_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb410_r(%rsp),%xmm6
        addps  %xmm5,%xmm6

    ## increment vctot (sum in xmm12)
        addps  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## update dvdasum
    addps  nb410_dvdasum(%rsp),%xmm7
    movaps %xmm7,nb410_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        movaps  %xmm6,%xmm5
        movaps  %xmm7,%xmm4
        shufps $0x1,%xmm5,%xmm5
        shufps $0x1,%xmm4,%xmm4

        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        addss  (%rsi,%rcx,4),%xmm7
        addss  (%rsi,%rdx,4),%xmm4
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)
        movss  %xmm7,(%rsi,%rcx,4)
        movss  %xmm4,(%rsi,%rdx,4)

    subps  %xmm3,%xmm9
    mulps  %xmm0,%xmm9 ## fscal

    movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

    mulps   nb410_dx(%rsp),%xmm9
    mulps   nb410_dy(%rsp),%xmm10
    mulps   nb410_dz(%rsp),%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb410_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movlps (%rsi,%r10,4),%xmm1 ## x3 y3 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2
        movhps (%rsi,%r11,4),%xmm1 ## x3 y3 x4 y4

    movaps %xmm9,%xmm8
    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    unpckhps %xmm10,%xmm8 ## x3 y3 x4 y4

    ## update fjx and fjy
        addps  %xmm9,%xmm0
        addps  %xmm8,%xmm1

        movlps %xmm0,(%rsi,%r8,4)
        movlps %xmm1,(%rsi,%r10,4)
        movhps %xmm0,(%rsi,%r9,4)
        movhps %xmm1,(%rsi,%r11,4)

    ## xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd $1,%xmm11,%xmm10 ## fjz2 - - -
    movhlps %xmm11,%xmm9     ## fjz3 - - -
    pshufd $3,%xmm11,%xmm8  ## fjz4 - - -

        addss  8(%rsi,%r8,4),%xmm11
        addss  8(%rsi,%r9,4),%xmm10
        addss  8(%rsi,%r10,4),%xmm9
        addss  8(%rsi,%r11,4),%xmm8
        movss  %xmm11,8(%rsi,%r8,4)
        movss  %xmm10,8(%rsi,%r9,4)
        movss  %xmm9,8(%rsi,%r10,4)
        movss  %xmm8,8(%rsi,%r11,4)

        ## should we do one more iteration? 
        subl $4,nb410_innerk(%rsp)
        jl    _nb_kernel410_x86_64_sse.nb410_finish_inner
        jmp   _nb_kernel410_x86_64_sse.nb410_unroll_loop
_nb_kernel410_x86_64_sse.nb410_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb410_innerk(%rsp)
        movl  nb410_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel410_x86_64_sse.nb410_dopair
        jmp   _nb_kernel410_x86_64_sse.nb410_checksingle
_nb_kernel410_x86_64_sse.nb410_dopair: 
        movq  nb410_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb410_innerjjnr(%rsp)

        ## load isaj
        movq nb410_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm2
        movss (%rsi,%rbx,4),%xmm6
    unpcklps %xmm6,%xmm2

        mulps  nb410_isai(%rsp),%xmm2

        movaps %xmm2,nb410_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb410_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb410_gbscale(%rsp)

    mulps nb410_iq(%rsp),%xmm2
        movq nb410_charge(%rbp),%rsi     ## base of charge[] 
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
    unpcklps %xmm6,%xmm3


        mulps %xmm2,%xmm3
        movaps %xmm3,nb410_qq(%rsp)

     ## vdw parameters
        movq nb410_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        movl (%rsi,%rbx,4),%r13d
        shll %r12d
        shll %r13d
    movl nb410_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb410_vdwparam(%rbp),%rsi
        movlps (%rsi,%r12,4),%xmm3
        movhps (%rsi,%r13,4),%xmm3

    xorps %xmm7,%xmm7
        movaps %xmm3,%xmm0
        shufps $136,%xmm7,%xmm0 ## 10001000
        shufps $221,%xmm7,%xmm3 ## 11011101

    movaps %xmm0,nb410_c6(%rsp)
    movaps %xmm3,nb410_c12(%rsp)

        movq nb410_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## j3
        lea  (%rbx,%rbx,2),%r9

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%r8,4),%xmm4       ## x1 y1 - - 
        movlps (%rsi,%r9,4),%xmm5       ## x2 y2 - - 

        movss 8(%rsi,%r8,4),%xmm6       ## z1 - - - 
        movss 8(%rsi,%r9,4),%xmm7       ## z2 - - - 

    unpcklps %xmm5,%xmm4 ## x1 x2 y1 y2
    movhlps  %xmm4,%xmm5 ## y1 y2 -  -
    unpcklps %xmm7,%xmm6 ## z1 z2 -  -

        ## calc dr 
        subps nb410_ix(%rsp),%xmm4
        subps nb410_iy(%rsp),%xmm5
        subps nb410_iz(%rsp),%xmm6

        ## store dr 
        movaps %xmm4,nb410_dx(%rsp)
        movaps %xmm5,nb410_dy(%rsp)
        movaps %xmm6,nb410_dz(%rsp)

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb410_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb410_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb410_r(%rsp)
        mulps nb410_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttps2dq %xmm4,%xmm5

    ## convert back to float
    cvtdq2ps  %xmm5,%xmm6

    ## multiply by 4
    pslld   $2,%xmm5

    ## move to integer registers
    movd    %xmm5,%r12d
    pshufd $1,%xmm5,%xmm5
    movd    %xmm5,%r13d

    ## calculate eps
    subps     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb410_GBtab(%rbp),%rsi

    movaps %xmm0,%xmm9 ## rinv
    mulps  %xmm9,%xmm9 ## rinvsq
    movaps %xmm9,%xmm10 ## rinvsq
    mulps  %xmm10,%xmm10 ## rinv4
    mulps  %xmm9,%xmm10 ## rinv6
    movaps %xmm10,%xmm11
    mulps  %xmm11,%xmm11 ## rinv12

    ## load table data
        movlps (%rsi,%r12,4),%xmm4  ## Y1 F1
        movlps (%rsi,%r13,4),%xmm5  ## Y2 F2
    unpcklps %xmm5,%xmm4        ## Y1 Y2 F1 F2
    movhlps  %xmm4,%xmm5        ## F1 F2

    mulps  nb410_c6(%rsp),%xmm10      ## vvdw6=c6*rinv6
        mulps  nb410_c12(%rsp),%xmm11     ## vvdw12=c12*rinv12     

        movaps %xmm11,%xmm9
        subps  %xmm10,%xmm11    ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addps  nb410_Vvdwtot(%rsp),%xmm11
    movlps %xmm11,nb410_Vvdwtot(%rsp)

        movlps 8(%rsi,%r12,4),%xmm6      ## G1 H1
        movlps 8(%rsi,%r13,4),%xmm7      ## G2 H2
    unpcklps %xmm7,%xmm6             ## G1 G2
    movhlps  %xmm6,%xmm7             ## H1 H2
    ## table data ready in xmm4-xmm7

    mulps  %xmm1,%xmm7  ## Heps
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm1,%xmm7      ## Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        addps  %xmm7,%xmm7      ## two*Heps2 
        movaps nb410_qq(%rsp),%xmm3

        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

    ## LJ forces
    mulps  nb410_six(%rsp),%xmm10
    mulps  nb410_twelve(%rsp),%xmm9
    subps  %xmm10,%xmm9
    mulps  %xmm0,%xmm9 ## (12*vnb12-6*vnb6)*rinv

    ## zero upper part of vcoul 
    xorps %xmm2,%xmm2
    movlhps %xmm2,%xmm5

        movq nb410_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulps nb410_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulps  nb410_r(%rsp),%xmm6
        addps  %xmm5,%xmm6

    xorps  %xmm4,%xmm4
    ## increment vctot (sum in xmm12)
        addps  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subps  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## zero upper half of dvda
    movlhps %xmm4,%xmm7

    ## update dvdasum
    addps  nb410_dvdasum(%rsp),%xmm7
    movaps %xmm7,nb410_dvdasum(%rsp)

        ## update j atoms dvdaj
        movaps  %xmm6,%xmm5
        shufps $0x1,%xmm5,%xmm5

        ## xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
        addss  (%rsi,%rax,4),%xmm6
        addss  (%rsi,%rbx,4),%xmm5
        movss  %xmm6,(%rsi,%rax,4)
        movss  %xmm5,(%rsi,%rbx,4)

    xorps %xmm7,%xmm7

    subps  %xmm3,%xmm9
    mulps  %xmm0,%xmm9 ## fscal

    movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

    mulps   nb410_dx(%rsp),%xmm9
    mulps   nb410_dy(%rsp),%xmm10
    mulps   nb410_dz(%rsp),%xmm11

    movlhps  %xmm7,%xmm9
    movlhps  %xmm7,%xmm10
    movlhps  %xmm7,%xmm11

        ## accumulate i forces
    addps %xmm9,%xmm13
    addps %xmm10,%xmm14
    addps %xmm11,%xmm15

        movq nb410_faction(%rbp),%rsi
        ## the fj's - start by accumulating x & y forces from memory 
        movlps (%rsi,%r8,4),%xmm0 ## x1 y1 - -
        movhps (%rsi,%r9,4),%xmm0 ## x1 y1 x2 y2

    unpcklps %xmm10,%xmm9 ## x1 y1 x2 y2
    addps    %xmm9,%xmm0

        movlps %xmm0,(%rsi,%r8,4)
        movhps %xmm0,(%rsi,%r9,4)

    ## z forces
    pshufd $1,%xmm11,%xmm8
    addss  8(%rsi,%r8,4),%xmm11
    addss  8(%rsi,%r9,4),%xmm8
    movss  %xmm11,8(%rsi,%r8,4)
    movss  %xmm8,8(%rsi,%r9,4)

_nb_kernel410_x86_64_sse.nb410_checksingle:     
        movl  nb410_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410_x86_64_sse.nb410_dosingle
        jmp    _nb_kernel410_x86_64_sse.nb410_updateouterdata
_nb_kernel410_x86_64_sse.nb410_dosingle: 
        movq nb410_charge(%rbp),%rsi
        movq nb410_invsqrta(%rbp),%rdx
        movq nb410_pos(%rbp),%rdi
        movq  nb410_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb410_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movaps nb410_isai(%rsp),%xmm2
        mulss  %xmm3,%xmm2

        movss %xmm2,nb410_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulss nb410_gbtsc(%rsp),%xmm1
        movss %xmm1,nb410_gbscale(%rsp)

    mulss nb410_iq(%rsp),%xmm2
        movq nb410_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        mulss %xmm2,%xmm3
        movss %xmm3,nb410_qq(%rsp)

    ## vdw parameters
        movq nb410_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r12d
        shll %r12d
    movl nb410_ntia(%rsp),%edi
        addl %edi,%r12d

        movq nb410_vdwparam(%rbp),%rsi
        movss (%rsi,%r12,4),%xmm0
        movss 4(%rsi,%r12,4),%xmm3
    movaps %xmm0,nb410_c6(%rsp)
    movaps %xmm3,nb410_c12(%rsp)

        movq nb410_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r8     ## jnr 

        ## move four coordinates to xmm0-xmm2   
        movss (%rsi,%r8,4),%xmm4
        movss 4(%rsi,%r8,4),%xmm5
        movss 8(%rsi,%r8,4),%xmm6

        ## calc dr 
        subss nb410_ix(%rsp),%xmm4
        subss nb410_iy(%rsp),%xmm5
        subss nb410_iz(%rsp),%xmm6

        ## store dr 
        movaps %xmm4,nb410_dx(%rsp)
        movaps %xmm5,nb410_dy(%rsp)
        movaps %xmm6,nb410_dz(%rsp)

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movaps nb410_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb410_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 
        mulss %xmm0,%xmm4       ## xmm4=r
        movaps %xmm4,nb410_r(%rsp)
        mulss nb410_gbscale(%rsp),%xmm4

    ## truncate and convert to integers
    cvttss2si %xmm4,%r12d

    ## convert back to float
    cvtsi2ss  %r12d,%xmm6

    ## multiply by 4
    shll  $2,%r12d

    ## calculate eps
    subss     %xmm6,%xmm4
    movaps    %xmm4,%xmm1 ##eps

        movq nb410_GBtab(%rbp),%rsi

    movaps %xmm0,%xmm9 ## rinv
    mulss  %xmm9,%xmm9 ## rinvsq
    movaps %xmm9,%xmm10 ## rinvsq
    mulss  %xmm10,%xmm10 ## rinv4
    mulss  %xmm9,%xmm10 ## rinv6
    movaps %xmm10,%xmm11
    mulss  %xmm11,%xmm11 ## rinv12

    ## load table data
        movss (%rsi,%r12,4),%xmm4
        movss 4(%rsi,%r12,4),%xmm5
        movss 8(%rsi,%r12,4),%xmm6
        movss 12(%rsi,%r12,4),%xmm7
    ## table data ready in xmm4-xmm7

    mulss  nb410_c6(%rsp),%xmm10      ## vvdw6=c6*rinv6
        mulss  nb410_c12(%rsp),%xmm11     ## vvdw12=c12*rinv12     

        movaps %xmm11,%xmm9
        subss  %xmm10,%xmm11    ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addss  nb410_Vvdwtot(%rsp),%xmm11
    movss %xmm11,nb410_Vvdwtot(%rsp)

    mulss  %xmm1,%xmm7  ## Heps
        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm1,%xmm7      ## Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        addss  %xmm7,%xmm7      ## two*Heps2 
        movss  nb410_qq(%rsp),%xmm3
        addss  %xmm6,%xmm7
        addss  %xmm5,%xmm7 ## xmm7=FF 
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulss  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point xmm5 contains vcoul and xmm3 fijC

    ## LJ forces
    mulss  nb410_six(%rsp),%xmm10
    mulss  nb410_twelve(%rsp),%xmm9
    subss  %xmm10,%xmm9
    mulss  %xmm0,%xmm9 ## (12*vnb12-6*vnb6)*rinv

        movq nb410_dvda(%rbp),%rsi

        ## Calculate dVda
        xorps  %xmm7,%xmm7
        mulss nb410_gbscale(%rsp),%xmm3
        movaps %xmm3,%xmm6
        mulss  nb410_r(%rsp),%xmm6
        addss  %xmm5,%xmm6

    ## increment vctot (sum in xmm12)
        addss  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subss  %xmm6,%xmm7
        movaps %xmm7,%xmm6

    ## update dvdasum
    addss  nb410_dvdasum(%rsp),%xmm7
    movss %xmm7,nb410_dvdasum(%rsp)

        ## update j atoms dvdaj
        addss  (%rsi,%rax,4),%xmm6
        movss  %xmm6,(%rsi,%rax,4)

    subss  %xmm3,%xmm9
    mulss  %xmm0,%xmm9 ## fscal

    movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

    mulss   nb410_dx(%rsp),%xmm9
    mulss   nb410_dy(%rsp),%xmm10
    mulss   nb410_dz(%rsp),%xmm11

        ## accumulate i forces
    addss %xmm9,%xmm13
    addss %xmm10,%xmm14
    addss %xmm11,%xmm15

        movq nb410_faction(%rbp),%rsi
    ## add to j forces
    addss  (%rsi,%r8,4),%xmm9
    addss  4(%rsi,%r8,4),%xmm10
    addss  8(%rsi,%r8,4),%xmm11
    movss  %xmm9,(%rsi,%r8,4)
    movss  %xmm10,4(%rsi,%r8,4)
    movss  %xmm11,8(%rsi,%r8,4)

_nb_kernel410_x86_64_sse.nb410_updateouterdata: 
        movl  nb410_ii3(%rsp),%ecx
        movq  nb410_faction(%rbp),%rdi
        movq  nb410_fshift(%rbp),%rsi
        movl  nb410_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm0
        movhlps %xmm14,%xmm1
        movhlps %xmm15,%xmm2
        addps  %xmm13,%xmm0
        addps  %xmm14,%xmm1
        addps  %xmm15,%xmm2
    movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 


        ## increment i force 
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## increment fshift force  
        movss  (%rsi,%rdx,4),%xmm3
        movss  4(%rsi,%rdx,4),%xmm4
        movss  8(%rsi,%rdx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,4(%rsi,%rdx,4)
        movss  %xmm5,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb410_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb410_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        ## accumulate 
        movhlps %xmm12,%xmm6
        addps  %xmm6,%xmm12     ## pos 0-1 in xmm12 have the sum now 
        movaps %xmm12,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm12

        ## add earlier value from mem 
        movq  nb410_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm12
        ## move back to mem 
        movss %xmm12,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb410_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb410_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate dVda and update it 
        movaps nb410_dvdasum(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        movl nb410_ii(%rsp),%edx
        movq nb410_dvda(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb410_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410_x86_64_sse.nb410_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410_n(%rsp)
        jmp _nb_kernel410_x86_64_sse.nb410_outer
_nb_kernel410_x86_64_sse.nb410_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410_x86_64_sse.nb410_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410_x86_64_sse.nb410_threadloop
_nb_kernel410_x86_64_sse.nb410_end: 

        movl nb410_nouter(%rsp),%eax
        movl nb410_ninner(%rsp),%ebx
        movq nb410_outeriter(%rbp),%rcx
        movq nb410_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $568,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



.globl nb_kernel410nf_x86_64_sse
.globl _nb_kernel410nf_x86_64_sse
nb_kernel410nf_x86_64_sse:      
_nb_kernel410nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb410nf_fshift, 16
.set nb410nf_gid, 24
.set nb410nf_pos, 32
.set nb410nf_faction, 40
.set nb410nf_charge, 48
.set nb410nf_p_facel, 56
.set nb410nf_argkrf, 64
.set nb410nf_argcrf, 72
.set nb410nf_Vc, 80
.set nb410nf_type, 88
.set nb410nf_p_ntype, 96
.set nb410nf_vdwparam, 104
.set nb410nf_Vvdw, 112
.set nb410nf_p_tabscale, 120
.set nb410nf_VFtab, 128
.set nb410nf_invsqrta, 136
.set nb410nf_dvda, 144
.set nb410nf_p_gbtabscale, 152
.set nb410nf_GBtab, 160
.set nb410nf_p_nthreads, 168
.set nb410nf_count, 176
.set nb410nf_mtx, 184
.set nb410nf_outeriter, 192
.set nb410nf_inneriter, 200
.set nb410nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb410nf_ix, 0
.set nb410nf_iy, 16
.set nb410nf_iz, 32
.set nb410nf_iq, 48
.set nb410nf_gbtsc, 64
.set nb410nf_qq, 80
.set nb410nf_c6, 96
.set nb410nf_c12, 112
.set nb410nf_vctot, 128
.set nb410nf_Vvdwtot, 144
.set nb410nf_half, 160
.set nb410nf_three, 176
.set nb410nf_isai, 192
.set nb410nf_isaprod, 208
.set nb410nf_gbscale, 224
.set nb410nf_nri, 240
.set nb410nf_iinr, 248
.set nb410nf_jindex, 256
.set nb410nf_jjnr, 264
.set nb410nf_shift, 272
.set nb410nf_shiftvec, 280
.set nb410nf_facel, 288
.set nb410nf_innerjjnr, 296
.set nb410nf_is3, 304
.set nb410nf_ii3, 308
.set nb410nf_ntia, 312
.set nb410nf_innerk, 316
.set nb410nf_n, 320
.set nb410nf_nn1, 324
.set nb410nf_ntype, 328
.set nb410nf_nouter, 332
.set nb410nf_ninner, 336

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $360,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb410nf_nouter(%rsp)
        movl %eax,nb410nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb410nf_nri(%rsp)
        movq %rsi,nb410nf_iinr(%rsp)
        movq %rdx,nb410nf_jindex(%rsp)
        movq %rcx,nb410nf_jjnr(%rsp)
        movq %r8,nb410nf_shift(%rsp)
        movq %r9,nb410nf_shiftvec(%rsp)
        movq nb410nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb410nf_ntype(%rsp)
        movq nb410nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb410nf_facel(%rsp)

        movq nb410nf_p_gbtabscale(%rbp),%rbx
        movss (%rbx),%xmm4
        shufps $0,%xmm4,%xmm4
        movaps %xmm4,nb410nf_gbtsc(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb410nf_half(%rsp)
        movss nb410nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb410nf_half(%rsp)
        movaps %xmm3,nb410nf_three(%rsp)

_nb_kernel410nf_x86_64_sse.nb410nf_threadloop: 
        movq  nb410nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel410nf_x86_64_sse.nb410nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410nf_x86_64_sse.nb410nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb410nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb410nf_n(%rsp)
        movl %ebx,nb410nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel410nf_x86_64_sse.nb410nf_outerstart
        jmp _nb_kernel410nf_x86_64_sse.nb410nf_end

_nb_kernel410nf_x86_64_sse.nb410nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410nf_nouter(%rsp),%ebx
        movl %ebx,nb410nf_nouter(%rsp)

_nb_kernel410nf_x86_64_sse.nb410nf_outer: 
        movq  nb410nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## ebx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb410nf_is3(%rsp)            ## store is3 

        movq  nb410nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb410nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii

        movq  nb410nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        mulss nb410nf_facel(%rsp),%xmm3
        shufps $0,%xmm3,%xmm3

        movq  nb410nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movss (%rdx,%rbx,4),%xmm4
        shufps $0,%xmm4,%xmm4

        movq  nb410nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb410nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb410nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb410nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addss (%rax,%rbx,4),%xmm0
        addss 4(%rax,%rbx,4),%xmm1
        addss 8(%rax,%rbx,4),%xmm2

        movaps %xmm3,nb410nf_iq(%rsp)
        movaps %xmm4,nb410nf_isai(%rsp)

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movaps %xmm0,nb410nf_ix(%rsp)
        movaps %xmm1,nb410nf_iy(%rsp)
        movaps %xmm2,nb410nf_iz(%rsp)

        movl  %ebx,nb410nf_ii3(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb410nf_vctot(%rsp)
        movaps %xmm4,nb410nf_Vvdwtot(%rsp)

        movq  nb410nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb410nf_pos(%rbp),%rsi
        movq  nb410nf_faction(%rbp),%rdi
        movq  nb410nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb410nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb410nf_ninner(%rsp),%ecx
        movl  %ecx,nb410nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb410nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel410nf_x86_64_sse.nb410nf_unroll_loop
        jmp   _nb_kernel410nf_x86_64_sse.nb410nf_finish_inner
_nb_kernel410nf_x86_64_sse.nb410nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb410nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 
        addq $16,nb410nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        ## load isa2
        movq nb410nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7
        movaps nb410nf_isai(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm3,%xmm2

        movaps %xmm2,nb410nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb410nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb410nf_gbscale(%rsp)

        movq nb410nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        mulps nb410nf_iq(%rsp),%xmm2
        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb410nf_qq(%rsp)

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb410nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb410nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb410nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx
        addl %edi,%ecx
        addl %edi,%edx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb410nf_c6(%rsp)
        movaps %xmm6,nb410nf_c12(%rsp)

        movq nb410nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   

        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ix-iz to xmm4-xmm6 
        movaps nb410nf_ix(%rsp),%xmm4
        movaps nb410nf_iy(%rsp),%xmm5
        movaps nb410nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb410nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb410nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb410nf_gbscale(%rsp),%xmm4

        movhlps %xmm4,%xmm5
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm5,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        cvtpi2ps %mm7,%xmm5
        movlhps %xmm5,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movq nb410nf_GBtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        ## load coulomb table
        movaps (%rsi,%rax,4),%xmm4
        movaps (%rsi,%rbx,4),%xmm5
        movaps (%rsi,%rcx,4),%xmm6
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm6,%xmm3
        shufps $0xEE,%xmm7,%xmm3
        shufps $0x44,%xmm7,%xmm6
        movaps %xmm4,%xmm7
        shufps $0xEE,%xmm5,%xmm7
        shufps $0x44,%xmm5,%xmm4
        movaps %xmm4,%xmm5
        shufps $0xDD,%xmm6,%xmm5
        shufps $0x88,%xmm6,%xmm4
        movaps %xmm7,%xmm6
        shufps $0x88,%xmm3,%xmm6
        shufps $0xDD,%xmm3,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 

        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb410nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## update vctot
        addps  nb410nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb410nf_vctot(%rsp)

        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb410nf_c6(%rsp),%xmm6
        mulps  nb410nf_c12(%rsp),%xmm4
        movaps nb410nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb410nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb410nf_innerk(%rsp)
        jl    _nb_kernel410nf_x86_64_sse.nb410nf_finish_inner
        jmp   _nb_kernel410nf_x86_64_sse.nb410nf_unroll_loop
_nb_kernel410nf_x86_64_sse.nb410nf_finish_inner: 
        ## check if at least two particles remain 
        addl $4,nb410nf_innerk(%rsp)
        movl  nb410nf_innerk(%rsp),%edx
        andl  $2,%edx
        jnz   _nb_kernel410nf_x86_64_sse.nb410nf_dopair
        jmp   _nb_kernel410nf_x86_64_sse.nb410nf_checksingle
_nb_kernel410nf_x86_64_sse.nb410nf_dopair: 
        movq  nb410nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        movl  4(%rcx),%ebx
        addq $8,nb410nf_innerjjnr(%rsp)

        xorps %xmm2,%xmm2
        movaps %xmm2,%xmm6

        ## load isa2
        movq nb410nf_invsqrta(%rbp),%rsi
        movss (%rsi,%rax,4),%xmm2
        movss (%rsi,%rbx,4),%xmm3
        unpcklps %xmm3,%xmm2    ## isa2 in xmm3(0,1)
        mulps  nb410nf_isai(%rsp),%xmm2
        movaps %xmm2,nb410nf_isaprod(%rsp)
        movaps %xmm2,%xmm1
        mulps nb410nf_gbtsc(%rsp),%xmm1
        movaps %xmm1,nb410nf_gbscale(%rsp)

        movq nb410nf_charge(%rbp),%rsi     ## base of charge[]  
        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rbx,4),%xmm6
        unpcklps %xmm6,%xmm3 ## 00001000 ;# xmm3(0,1) has the charges 

        mulps  nb410nf_iq(%rsp),%xmm2
        mulps  %xmm2,%xmm3
        movaps %xmm3,nb410nf_qq(%rsp)

        movq nb410nf_type(%rbp),%rsi
        movl  %eax,%ecx
        movl  %ebx,%edx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb410nf_vdwparam(%rbp),%rsi
        shll %ecx
        shll %edx
        movl nb410nf_ntia(%rsp),%edi
        addl %edi,%ecx
        addl %edi,%edx
        movlps (%rsi,%rcx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm6
        movq nb410nf_pos(%rbp),%rdi

        movaps %xmm6,%xmm4
        shufps $8,%xmm4,%xmm4 ## 00001000        
        shufps $13,%xmm6,%xmm6 ## 00001101
        movlhps %xmm7,%xmm4
        movlhps %xmm7,%xmm6

        movaps %xmm4,nb410nf_c6(%rsp)
        movaps %xmm6,nb410nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax
        lea  (%rbx,%rbx,2),%rbx
        ## move coordinates to xmm0-xmm2 
        movlps (%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2
        movhps (%rdi,%rbx,4),%xmm1
        movss 8(%rdi,%rbx,4),%xmm0

        movlhps %xmm7,%xmm3

        shufps $0,%xmm0,%xmm2

        movaps %xmm1,%xmm0

        shufps $136,%xmm2,%xmm2 ## 10001000

        shufps $136,%xmm0,%xmm0 ## 10001000
        shufps $221,%xmm1,%xmm1 ## 11011101

        movq   nb410nf_faction(%rbp),%rdi
        ## move ix-iz to xmm4-xmm6 
        xorps   %xmm7,%xmm7

        movaps nb410nf_ix(%rsp),%xmm4
        movaps nb410nf_iy(%rsp),%xmm5
        movaps nb410nf_iz(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb410nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb410nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        mulps %xmm0,%xmm4       ## xmm4=r 
        mulps nb410nf_gbscale(%rsp),%xmm4

        cvttps2pi %xmm4,%mm6    ## mm6 contain lu indices 
        cvtpi2ps %mm6,%xmm6
        subps %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6

        movq nb410nf_GBtab(%rbp),%rsi
        movd %mm6,%ecx
        psrlq $32,%mm6
        movd %mm6,%edx

        ## load coulomb table
        movaps (%rsi,%rcx,4),%xmm4
        movaps (%rsi,%rdx,4),%xmm7
        ## transpose, using xmm3 for scratch
        movaps %xmm4,%xmm6
        unpcklps %xmm7,%xmm4    ## Y1 Y2 F1 F2 
        unpckhps %xmm7,%xmm6    ## G1 G2 H1 H2
        movhlps  %xmm4,%xmm5    ## F1 F2 
        movhlps  %xmm6,%xmm7    ## H1 H2
        ## coulomb table ready, in xmm4-xmm7    

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp      
        movaps nb410nf_qq(%rsp),%xmm3
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm3,%xmm5 ## vcoul=qq*VV  

        addps  nb410nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb410nf_vctot(%rsp)

        ## L-J 
        movaps %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## at this point mm5 contains vcoul and mm3 fijC 
        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 

        movaps %xmm4,%xmm6
        mulps  %xmm4,%xmm6

        mulps  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulps  nb410nf_c6(%rsp),%xmm6
        mulps  nb410nf_c12(%rsp),%xmm4
        movaps nb410nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movaps %xmm7,nb410nf_Vvdwtot(%rsp)

_nb_kernel410nf_x86_64_sse.nb410nf_checksingle: 
        movl  nb410nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410nf_x86_64_sse.nb410nf_dosingle
        jmp    _nb_kernel410nf_x86_64_sse.nb410nf_updateouterdata
_nb_kernel410nf_x86_64_sse.nb410nf_dosingle: 
        movq nb410nf_charge(%rbp),%rsi
        movq nb410nf_invsqrta(%rbp),%rdx
        movq nb410nf_pos(%rbp),%rdi
        movq  nb410nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax
        xorps  %xmm2,%xmm2
        movaps %xmm2,%xmm6
        movss (%rdx,%rax,4),%xmm2       ## isa2
        mulss nb410nf_isai(%rsp),%xmm2
        movss %xmm2,nb410nf_isaprod(%rsp)
        movss %xmm2,%xmm1
        mulss nb410nf_gbtsc(%rsp),%xmm1
        movss %xmm1,nb410nf_gbscale(%rsp)

        mulss  nb410nf_iq(%rsp),%xmm2
        movss (%rsi,%rax,4),%xmm6       ## xmm6(0) has the charge       
        mulss  %xmm2,%xmm6
        movss %xmm6,nb410nf_qq(%rsp)

        movq nb410nf_type(%rbp),%rsi
        movl %eax,%ecx
        movl (%rsi,%rcx,4),%ecx
        movq nb410nf_vdwparam(%rbp),%rsi
        shll %ecx
        addl nb410nf_ntia(%rsp),%ecx
        movlps (%rsi,%rcx,4),%xmm6
        movaps %xmm6,%xmm4
        shufps $252,%xmm4,%xmm4 ## 11111100     
        shufps $253,%xmm6,%xmm6 ## 11111101     

        movaps %xmm4,nb410nf_c6(%rsp)
        movaps %xmm6,nb410nf_c12(%rsp)

        lea  (%rax,%rax,2),%rax

        ## move coordinates to xmm0-xmm2 
        movss (%rdi,%rax,4),%xmm0
        movss 4(%rdi,%rax,4),%xmm1
        movss 8(%rdi,%rax,4),%xmm2

        movaps nb410nf_ix(%rsp),%xmm4
        movaps nb410nf_iy(%rsp),%xmm5
        movaps nb410nf_iz(%rsp),%xmm6

        ## calc dr 
        subss %xmm0,%xmm4
        subss %xmm1,%xmm5
        subss %xmm2,%xmm6

        ## square it 
        mulss %xmm4,%xmm4
        mulss %xmm5,%xmm5
        mulss %xmm6,%xmm6
        addss %xmm5,%xmm4
        addss %xmm6,%xmm4
        ## rsq in xmm4 

        rsqrtss %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulss %xmm5,%xmm5
        movss nb410nf_three(%rsp),%xmm1
        mulss %xmm4,%xmm5       ## rsq*lu*lu                    
        movss nb410nf_half(%rsp),%xmm0
        subss %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulss %xmm2,%xmm1
        mulss %xmm1,%xmm0       ## xmm0=rinv 

        mulss %xmm0,%xmm4       ## xmm4=r 
        mulss nb410nf_gbscale(%rsp),%xmm4

        cvttss2si %xmm4,%ebx    ## mm6 contain lu indices 
        cvtsi2ss %ebx,%xmm6
        subss %xmm6,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%ebx
        movq nb410nf_GBtab(%rbp),%rsi

        movaps (%rsi,%rbx,4),%xmm4
        movhlps %xmm4,%xmm6
        movaps %xmm4,%xmm5
        movaps %xmm6,%xmm7
        shufps $1,%xmm5,%xmm5
        shufps $1,%xmm7,%xmm7
        ## table ready in xmm4-xmm7 

        mulss  %xmm1,%xmm6      ## xmm6=Geps 
        mulss  %xmm2,%xmm7      ## xmm7=Heps2 
        addss  %xmm6,%xmm5
        addss  %xmm7,%xmm5      ## xmm5=Fp      
        movss nb410nf_qq(%rsp),%xmm3
        mulss  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addss  %xmm4,%xmm5 ## xmm5=VV 
        mulss  %xmm3,%xmm5 ## vcoul=qq*VV  
        addss  nb410nf_vctot(%rsp),%xmm5
        movss %xmm5,nb410nf_vctot(%rsp)

        ## L-J 
        movaps %xmm0,%xmm4
        mulss  %xmm0,%xmm4      ## xmm4=rinvsq 

        movaps %xmm4,%xmm6
        mulss  %xmm4,%xmm6

        mulss  %xmm4,%xmm6      ## xmm6=rinvsix 
        movaps %xmm6,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb410nf_c6(%rsp),%xmm6
        mulss  nb410nf_c12(%rsp),%xmm4
        movss nb410nf_Vvdwtot(%rsp),%xmm7
        addps  %xmm4,%xmm7
        subps  %xmm6,%xmm7
        movss %xmm7,nb410nf_Vvdwtot(%rsp)

_nb_kernel410nf_x86_64_sse.nb410nf_updateouterdata: 
        ## get n from stack
        movl nb410nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb410nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb410nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb410nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb410nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb410nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb410nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410nf_x86_64_sse.nb410nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410nf_n(%rsp)
        jmp _nb_kernel410nf_x86_64_sse.nb410nf_outer
_nb_kernel410nf_x86_64_sse.nb410nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410nf_x86_64_sse.nb410nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410nf_x86_64_sse.nb410nf_threadloop
_nb_kernel410nf_x86_64_sse.nb410nf_end: 

        movl nb410nf_nouter(%rsp),%eax
        movl nb410nf_ninner(%rsp),%ebx
        movq nb410nf_outeriter(%rbp),%rcx
        movq nb410nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $360,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




