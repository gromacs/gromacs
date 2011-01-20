##
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











.globl nb_kernel201_x86_64_sse
.globl _nb_kernel201_x86_64_sse
nb_kernel201_x86_64_sse:        
_nb_kernel201_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb201_fshift, 16
.set nb201_gid, 24
.set nb201_pos, 32
.set nb201_faction, 40
.set nb201_charge, 48
.set nb201_p_facel, 56
.set nb201_argkrf, 64
.set nb201_argcrf, 72
.set nb201_Vc, 80
.set nb201_type, 88
.set nb201_p_ntype, 96
.set nb201_vdwparam, 104
.set nb201_Vvdw, 112
.set nb201_p_tabscale, 120
.set nb201_VFtab, 128
.set nb201_invsqrta, 136
.set nb201_dvda, 144
.set nb201_p_gbtabscale, 152
.set nb201_GBtab, 160
.set nb201_p_nthreads, 168
.set nb201_count, 176
.set nb201_mtx, 184
.set nb201_outeriter, 192
.set nb201_inneriter, 200
.set nb201_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201_ixO, 0
.set nb201_iyO, 16
.set nb201_izO, 32
.set nb201_ixH1, 48
.set nb201_iyH1, 64
.set nb201_izH1, 80
.set nb201_ixH2, 96
.set nb201_iyH2, 112
.set nb201_izH2, 128
.set nb201_iqO, 144
.set nb201_iqH, 160
.set nb201_dxO, 176
.set nb201_dyO, 192
.set nb201_dzO, 208
.set nb201_dxH1, 224
.set nb201_dyH1, 240
.set nb201_dzH1, 256
.set nb201_dxH2, 272
.set nb201_dyH2, 288
.set nb201_dzH2, 304
.set nb201_qqO, 320
.set nb201_qqH, 336
.set nb201_vctot, 352
.set nb201_fixO, 384
.set nb201_fiyO, 400
.set nb201_fizO, 416
.set nb201_fixH1, 432
.set nb201_fiyH1, 448
.set nb201_fizH1, 464
.set nb201_fixH2, 480
.set nb201_fiyH2, 496
.set nb201_fizH2, 512
.set nb201_fjx, 528
.set nb201_fjy, 544
.set nb201_fjz, 560
.set nb201_half, 576
.set nb201_three, 592
.set nb201_two, 608
.set nb201_krf, 624
.set nb201_crf, 640
.set nb201_krsqO, 656
.set nb201_krsqH1, 672
.set nb201_krsqH2, 688
.set nb201_nri, 704
.set nb201_iinr, 712
.set nb201_jindex, 720
.set nb201_jjnr, 728
.set nb201_shift, 736
.set nb201_shiftvec, 744
.set nb201_facel, 752
.set nb201_innerjjnr, 760
.set nb201_is3, 768
.set nb201_ii3, 772
.set nb201_innerk, 776
.set nb201_n, 780
.set nb201_nn1, 784
.set nb201_nouter, 788
.set nb201_ninner, 792

        push %rbp
        movq %rsp,%rbp
        push %rbx


        push %r12
        push %r13
        push %r14
        push %r15

        subq $808,%rsp          ## local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb201_nouter(%rsp)
        movl %eax,nb201_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb201_nri(%rsp)
        movq %rsi,nb201_iinr(%rsp)
        movq %rdx,nb201_jindex(%rsp)
        movq %rcx,nb201_jjnr(%rsp)
        movq %r8,nb201_shift(%rsp)
        movq %r9,nb201_shiftvec(%rsp)
        movq nb201_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb201_facel(%rsp)


        movq nb201_argkrf(%rbp),%rsi
        movq nb201_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb201_krf(%rsp)
        movaps %xmm2,nb201_crf(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb201_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb201_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb201_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb201_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb201_iqO(%rsp)
        movaps %xmm4,nb201_iqH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb201_half(%rsp)
        movss nb201_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb201_half(%rsp)
        movaps %xmm2,nb201_two(%rsp)
        movaps %xmm3,nb201_three(%rsp)


_nb_kernel201_x86_64_sse.nb201_threadloop: 
        movq  nb201_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel201_x86_64_sse.nb201_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201_x86_64_sse.nb201_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201_n(%rsp)
        movl %ebx,nb201_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201_x86_64_sse.nb201_outerstart
        jmp _nb_kernel201_x86_64_sse.nb201_end

_nb_kernel201_x86_64_sse.nb201_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201_nouter(%rsp),%ebx
        movl %ebx,nb201_nouter(%rsp)

_nb_kernel201_x86_64_sse.nb201_outer: 
        movq  nb201_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb201_is3(%rsp)      ## store is3 

        movq  nb201_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb201_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb201_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb201_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb201_ixO(%rsp)
        movaps %xmm4,nb201_iyO(%rsp)
        movaps %xmm5,nb201_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb201_ixH1(%rsp)
        movaps %xmm1,nb201_iyH1(%rsp)
        movaps %xmm2,nb201_izH1(%rsp)
        movaps %xmm3,nb201_ixH2(%rsp)
        movaps %xmm4,nb201_iyH2(%rsp)
        movaps %xmm5,nb201_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb201_vctot(%rsp)
        movaps %xmm4,nb201_fixO(%rsp)
        movaps %xmm4,nb201_fiyO(%rsp)
        movaps %xmm4,nb201_fizO(%rsp)
        movaps %xmm4,nb201_fixH1(%rsp)
        movaps %xmm4,nb201_fiyH1(%rsp)
        movaps %xmm4,nb201_fizH1(%rsp)
        movaps %xmm4,nb201_fixH2(%rsp)
        movaps %xmm4,nb201_fiyH2(%rsp)
        movaps %xmm4,nb201_fizH2(%rsp)

        movq  nb201_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb201_pos(%rbp),%rsi
        movq  nb201_faction(%rbp),%rdi
        movq  nb201_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb201_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb201_ninner(%rsp),%ecx
        movl  %ecx,nb201_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb201_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel201_x86_64_sse.nb201_unroll_loop
        jmp   _nb_kernel201_x86_64_sse.nb201_odd_inner
_nb_kernel201_x86_64_sse.nb201_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb201_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb201_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb201_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb201_iqO(%rsp),%xmm3
        mulps  nb201_iqH(%rsp),%xmm4

        movaps  %xmm3,nb201_qqO(%rsp)
        movaps  %xmm4,nb201_qqH(%rsp)

        movq nb201_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx     ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four j coordinates to xmm0-xmm2         
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

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb201_ixO(%rsp),%xmm0
    subps nb201_iyO(%rsp),%xmm1
    subps nb201_izO(%rsp),%xmm2
    subps nb201_ixH1(%rsp),%xmm3
    subps nb201_iyH1(%rsp),%xmm4
    subps nb201_izH1(%rsp),%xmm5
    subps nb201_ixH2(%rsp),%xmm6
    subps nb201_iyH2(%rsp),%xmm7
    subps nb201_izH2(%rsp),%xmm8

        movaps %xmm0,nb201_dxO(%rsp)
        movaps %xmm1,nb201_dyO(%rsp)
        movaps %xmm2,nb201_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb201_dxH1(%rsp)
        movaps %xmm4,nb201_dyH1(%rsp)
        movaps %xmm5,nb201_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb201_dxH2(%rsp)
        movaps %xmm7,nb201_dyH2(%rsp)
        movaps %xmm8,nb201_dzH2(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb201_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  nb201_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvO
        mulps   %xmm4,%xmm10 ## rinvH1
    mulps   %xmm4,%xmm11 ## rinvH2

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb201_krf(%rsp),%xmm2
    mulps  %xmm9,%xmm9  ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  %xmm2,%xmm0 ## k*rsq
    mulps  %xmm2,%xmm3
    mulps  %xmm2,%xmm6
    movaps %xmm0,%xmm2 ## copy of k*rsq
    movaps %xmm3,%xmm5
    movaps %xmm6,%xmm8
    addps  %xmm1,%xmm2 ## rinv+krsq
    addps  %xmm4,%xmm5
    addps  %xmm7,%xmm8
    movaps nb201_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb201_qqO(%rsp),%xmm12
    movaps nb201_qqH(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm13,%xmm4
    mulps  %xmm13,%xmm7
    addps  nb201_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,nb201_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm9 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm10 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm9 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm10 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm11    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm12    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm6     ## jzb
    movss  8(%rdi,%rdx,4),%xmm8     ## jzd
    movlhps %xmm6,%xmm11 ## jza  -  jzb  -
    movlhps %xmm8,%xmm12 ## jzc  -  jzd -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jza jzb jzc jzd

    ## xmm9: jxa jya jxb jyb 
    ## xmm10: jxc jyc jxd jyd
    ## xmm11: jza jzb jzc jzd

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb201_dxO(%rsp),%xmm0
        mulps nb201_dyO(%rsp),%xmm1
        mulps nb201_dzO(%rsp),%xmm2
        mulps nb201_dxH1(%rsp),%xmm3
        mulps nb201_dyH1(%rsp),%xmm4
        mulps nb201_dzH1(%rsp),%xmm5
        mulps nb201_dxH2(%rsp),%xmm6
        mulps nb201_dyH2(%rsp),%xmm7
        mulps nb201_dzH2(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb201_fixO(%rsp),%xmm0
    addps nb201_fiyO(%rsp),%xmm1
    addps nb201_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb201_fixH1(%rsp),%xmm3
    addps nb201_fiyH1(%rsp),%xmm4
    addps nb201_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb201_fixH2(%rsp),%xmm6
    addps nb201_fiyH2(%rsp),%xmm7
    addps nb201_fizH2(%rsp),%xmm8

    movaps %xmm0,nb201_fixO(%rsp)
    movaps %xmm1,nb201_fiyO(%rsp)
    movaps %xmm2,nb201_fizO(%rsp)
    movaps %xmm3,nb201_fixH1(%rsp)
    movaps %xmm4,nb201_fiyH1(%rsp)
    movaps %xmm5,nb201_fizH1(%rsp)
    movaps %xmm6,nb201_fixH2(%rsp)
    movaps %xmm7,nb201_fiyH2(%rsp)
    movaps %xmm8,nb201_fizH2(%rsp)

    ## xmm9 = fjx
    ## xmm10 = fjy
    ## xmm11 = fjz
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12 ## fjzc fjzd

    movlps %xmm9,(%rdi,%rax,4)
    movhps %xmm9,(%rdi,%rbx,4)
    movlps %xmm10,(%rdi,%rcx,4)
    movhps %xmm10,(%rdi,%rdx,4)
    movss  %xmm11,8(%rdi,%rax,4)
    movss  %xmm12,8(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,8(%rdi,%rbx,4)
    movss  %xmm12,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb201_innerk(%rsp)
        jl    _nb_kernel201_x86_64_sse.nb201_odd_inner
        jmp   _nb_kernel201_x86_64_sse.nb201_unroll_loop
_nb_kernel201_x86_64_sse.nb201_odd_inner: 
        addl $4,nb201_innerk(%rsp)
        jnz   _nb_kernel201_x86_64_sse.nb201_odd_loop
        jmp   _nb_kernel201_x86_64_sse.nb201_updateouterdata
_nb_kernel201_x86_64_sse.nb201_odd_loop: 
        movq  nb201_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb201_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb201_iqO(%rsp),%xmm4
        movq nb201_charge(%rbp),%rsi
        movhps nb201_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb201_qqO(%rsp)    ## use oxygen qq for storage 

        movq nb201_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb201_ixO(%rsp),%xmm0
        movss nb201_iyO(%rsp),%xmm1
        movss nb201_izO(%rsp),%xmm2

        movlps nb201_ixH1(%rsp),%xmm6
        movlps nb201_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb201_iyH1(%rsp),%xmm6
        movlps nb201_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb201_izH1(%rsp),%xmm6
        movlps nb201_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb201_dxO(%rsp)
        movaps %xmm4,nb201_dyO(%rsp)
        movaps %xmm5,nb201_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb201_krf(%rsp),%xmm0
        movaps %xmm0,nb201_krsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb201_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb201_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000

        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb201_krsqO(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb201_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201_two(%rsp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb201_qqO(%rsp),%xmm0    ## xmm0=vcoul 
        mulps  nb201_qqO(%rsp),%xmm1    ## xmm1=coul part of fs 


        mulps  %xmm1,%xmm4      ## xmm4=total fscal 
        addps  nb201_vctot(%rsp),%xmm0
        movaps %xmm0,nb201_vctot(%rsp)

        movaps nb201_dxO(%rsp),%xmm0
        movaps nb201_dyO(%rsp),%xmm1
        movaps nb201_dzO(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb201_fixO(%rsp),%xmm3
        movss  nb201_fiyO(%rsp),%xmm4
        movss  nb201_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb201_fixO(%rsp)
        movss  %xmm4,nb201_fiyO(%rsp)
        movss  %xmm5,nb201_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb201_fixH1(%rsp),%xmm3
        addss  nb201_fiyH1(%rsp),%xmm4
        addss  nb201_fizH1(%rsp),%xmm5
        movss  %xmm3,nb201_fixH1(%rsp)
        movss  %xmm4,nb201_fiyH1(%rsp)
        movss  %xmm5,nb201_fizH1(%rsp)          ## updated the H1 force 

        movq nb201_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb201_fixH2(%rsp),%xmm3
        addss  nb201_fiyH2(%rsp),%xmm4
        addss  nb201_fizH2(%rsp),%xmm5
        movss  %xmm3,nb201_fixH2(%rsp)
        movss  %xmm4,nb201_fiyH2(%rsp)
        movss  %xmm5,nb201_fizH2(%rsp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2   ## z sum in xmm2 
        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb201_innerk(%rsp)
        jz    _nb_kernel201_x86_64_sse.nb201_updateouterdata
        jmp   _nb_kernel201_x86_64_sse.nb201_odd_loop
_nb_kernel201_x86_64_sse.nb201_updateouterdata: 
        movl  nb201_ii3(%rsp),%ecx
        movq  nb201_faction(%rbp),%rdi
        movq  nb201_fshift(%rbp),%rsi
        movl  nb201_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb201_fixO(%rsp),%xmm0
        movaps nb201_fiyO(%rsp),%xmm1
        movaps nb201_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

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

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb201_fixH1(%rsp),%xmm0
        movaps nb201_fiyH1(%rsp),%xmm1
        movaps nb201_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

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
        movss  12(%rdi,%rcx,4),%xmm3
        movss  16(%rdi,%rcx,4),%xmm4
        movss  20(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,12(%rdi,%rcx,4)
        movss  %xmm4,16(%rdi,%rcx,4)
        movss  %xmm5,20(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb201_fixH2(%rsp),%xmm0
        movaps nb201_fiyH2(%rsp),%xmm1
        movaps nb201_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

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
        movss  24(%rdi,%rcx,4),%xmm3
        movss  28(%rdi,%rcx,4),%xmm4
        movss  32(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,24(%rdi,%rcx,4)
        movss  %xmm4,28(%rdi,%rcx,4)
        movss  %xmm5,32(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4
        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb201_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb201_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb201_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb201_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb201_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201_x86_64_sse.nb201_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201_n(%rsp)
        jmp _nb_kernel201_x86_64_sse.nb201_outer
_nb_kernel201_x86_64_sse.nb201_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201_x86_64_sse.nb201_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201_x86_64_sse.nb201_threadloop
_nb_kernel201_x86_64_sse.nb201_end: 

        movl nb201_nouter(%rsp),%eax
        movl nb201_ninner(%rsp),%ebx
        movq nb201_outeriter(%rbp),%rcx
        movq nb201_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $808,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel201nf_x86_64_sse
.globl _nb_kernel201nf_x86_64_sse
nb_kernel201nf_x86_64_sse:      
_nb_kernel201nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb201nf_fshift, 16
.set nb201nf_gid, 24
.set nb201nf_pos, 32
.set nb201nf_faction, 40
.set nb201nf_charge, 48
.set nb201nf_p_facel, 56
.set nb201nf_argkrf, 64
.set nb201nf_argcrf, 72
.set nb201nf_Vc, 80
.set nb201nf_type, 88
.set nb201nf_p_ntype, 96
.set nb201nf_vdwparam, 104
.set nb201nf_Vvdw, 112
.set nb201nf_p_tabscale, 120
.set nb201nf_VFtab, 128
.set nb201nf_invsqrta, 136
.set nb201nf_dvda, 144
.set nb201nf_p_gbtabscale, 152
.set nb201nf_GBtab, 160
.set nb201nf_p_nthreads, 168
.set nb201nf_count, 176
.set nb201nf_mtx, 184
.set nb201nf_outeriter, 192
.set nb201nf_inneriter, 200
.set nb201nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb201nf_ixO, 0
.set nb201nf_iyO, 16
.set nb201nf_izO, 32
.set nb201nf_ixH1, 48
.set nb201nf_iyH1, 64
.set nb201nf_izH1, 80
.set nb201nf_ixH2, 96
.set nb201nf_iyH2, 112
.set nb201nf_izH2, 128
.set nb201nf_iqO, 144
.set nb201nf_iqH, 160
.set nb201nf_qqO, 176
.set nb201nf_qqH, 192
.set nb201nf_vctot, 208
.set nb201nf_half, 224
.set nb201nf_three, 240
.set nb201nf_krf, 256
.set nb201nf_crf, 272
.set nb201nf_krsqO, 288
.set nb201nf_krsqH1, 304
.set nb201nf_krsqH2, 320
.set nb201nf_is3, 336
.set nb201nf_ii3, 340
.set nb201nf_innerjjnr, 344
.set nb201nf_nri, 352
.set nb201nf_iinr, 360
.set nb201nf_jindex, 368
.set nb201nf_jjnr, 376
.set nb201nf_shift, 384
.set nb201nf_shiftvec, 392
.set nb201nf_facel, 400
.set nb201nf_innerk, 408
.set nb201nf_n, 412
.set nb201nf_nn1, 416
.set nb201nf_nouter, 420
.set nb201nf_ninner, 424

        push %rbp
        movq %rsp,%rbp
        push %rbx


        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)
        emms

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb201nf_nouter(%rsp)
        movl %eax,nb201nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb201nf_nri(%rsp)
        movq %rsi,nb201nf_iinr(%rsp)
        movq %rdx,nb201nf_jindex(%rsp)
        movq %rcx,nb201nf_jjnr(%rsp)
        movq %r8,nb201nf_shift(%rsp)
        movq %r9,nb201nf_shiftvec(%rsp)
        movq nb201nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb201nf_facel(%rsp)

        movq nb201nf_argkrf(%rbp),%rsi
        movq nb201nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb201nf_krf(%rsp)
        movaps %xmm2,nb201nf_crf(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb201nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb201nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb201nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb201nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb201nf_iqO(%rsp)
        movaps %xmm4,nb201nf_iqH(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb201nf_half(%rsp)
        movss nb201nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb201nf_half(%rsp)
        movaps %xmm3,nb201nf_three(%rsp)

_nb_kernel201nf_x86_64_sse.nb201nf_threadloop: 
        movq  nb201nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel201nf_x86_64_sse.nb201nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel201nf_x86_64_sse.nb201nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb201nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb201nf_n(%rsp)
        movl %ebx,nb201nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel201nf_x86_64_sse.nb201nf_outerstart
        jmp _nb_kernel201nf_x86_64_sse.nb201nf_end

_nb_kernel201nf_x86_64_sse.nb201nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb201nf_nouter(%rsp),%ebx
        movl %ebx,nb201nf_nouter(%rsp)

_nb_kernel201nf_x86_64_sse.nb201nf_outer: 
        movq  nb201nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb201nf_is3(%rsp)            ## store is3 

        movq  nb201nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb201nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb201nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb201nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb201nf_ixO(%rsp)
        movaps %xmm4,nb201nf_iyO(%rsp)
        movaps %xmm5,nb201nf_izO(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%rax,%rbx,4),%xmm0
        addss 16(%rax,%rbx,4),%xmm1
        addss 20(%rax,%rbx,4),%xmm2
        addss 24(%rax,%rbx,4),%xmm3
        addss 28(%rax,%rbx,4),%xmm4
        addss 32(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb201nf_ixH1(%rsp)
        movaps %xmm1,nb201nf_iyH1(%rsp)
        movaps %xmm2,nb201nf_izH1(%rsp)
        movaps %xmm3,nb201nf_ixH2(%rsp)
        movaps %xmm4,nb201nf_iyH2(%rsp)
        movaps %xmm5,nb201nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb201nf_vctot(%rsp)

        movq  nb201nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb201nf_pos(%rbp),%rsi
        movq  nb201nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb201nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb201nf_ninner(%rsp),%ecx
        movl  %ecx,nb201nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb201nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel201nf_x86_64_sse.nb201nf_unroll_loop
        jmp   _nb_kernel201nf_x86_64_sse.nb201nf_odd_inner
_nb_kernel201nf_x86_64_sse.nb201nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb201nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb201nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb201nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb201nf_iqO(%rsp),%xmm3
        mulps  nb201nf_iqH(%rsp),%xmm4

        movaps  %xmm3,nb201nf_qqO(%rsp)
        movaps  %xmm4,nb201nf_qqH(%rsp)

        movq nb201nf_pos(%rbp),%rsi        ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb201nf_ixO(%rsp),%xmm4
        movaps nb201nf_iyO(%rsp),%xmm5
        movaps nb201nf_izO(%rsp),%xmm6

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
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb201nf_ixH1(%rsp),%xmm4
        movaps nb201nf_iyH1(%rsp),%xmm5
        movaps nb201nf_izH1(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb201nf_ixH2(%rsp),%xmm3
        movaps nb201nf_iyH2(%rsp),%xmm4
        movaps nb201nf_izH2(%rsp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb201nf_krf(%rsp),%xmm0
        mulps  nb201nf_krf(%rsp),%xmm1
        mulps  nb201nf_krf(%rsp),%xmm2

        movaps %xmm0,nb201nf_krsqH2(%rsp)
        movaps %xmm1,nb201nf_krsqH1(%rsp)
        movaps %xmm2,nb201nf_krsqO(%rsp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb201nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb201nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 

        movaps %xmm7,%xmm0
        movaps nb201nf_krsqO(%rsp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb201nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201nf_qqO(%rsp),%xmm0

        addps  nb201nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb201nf_vctot(%rsp)

        ## H1 interactions 
        movaps  %xmm6,%xmm7
        movaps  nb201nf_krsqH1(%rsp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb201nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps   nb201nf_qqH(%rsp),%xmm6   ## vcoul 
        addps   nb201nf_vctot(%rsp),%xmm6

        ## H2 interactions 
        movaps  %xmm5,%xmm7
        movaps  nb201nf_krsqH2(%rsp),%xmm0
        addps   %xmm0,%xmm5     ## xmm6=rinv+ krsq 
        subps   nb201nf_crf(%rsp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps   nb201nf_qqH(%rsp),%xmm5   ## vcoul 
        addps  %xmm5,%xmm6
        movaps %xmm6,nb201nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb201nf_innerk(%rsp)
        jl    _nb_kernel201nf_x86_64_sse.nb201nf_odd_inner
        jmp   _nb_kernel201nf_x86_64_sse.nb201nf_unroll_loop
_nb_kernel201nf_x86_64_sse.nb201nf_odd_inner: 
        addl $4,nb201nf_innerk(%rsp)
        jnz   _nb_kernel201nf_x86_64_sse.nb201nf_odd_loop
        jmp   _nb_kernel201nf_x86_64_sse.nb201nf_updateouterdata
_nb_kernel201nf_x86_64_sse.nb201nf_odd_loop: 
        movq  nb201nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb201nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb201nf_iqO(%rsp),%xmm4
        movq nb201nf_charge(%rbp),%rsi
        movhps nb201nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb201nf_qqO(%rsp)          ## use oxygen qq for storage 

        movq nb201nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb201nf_ixO(%rsp),%xmm3
        movss nb201nf_iyO(%rsp),%xmm4
        movss nb201nf_izO(%rsp),%xmm5

        movlps nb201nf_ixH1(%rsp),%xmm6
        movlps nb201nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb201nf_iyH1(%rsp),%xmm6
        movlps nb201nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb201nf_izH1(%rsp),%xmm6
        movlps nb201nf_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb201nf_krf(%rsp),%xmm0
        movaps %xmm0,nb201nf_krsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb201nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb201nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 

        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      

        movaps nb201nf_krsqO(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb201nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb201nf_qqO(%rsp),%xmm0          ## xmm0=vcoul 
        addps  nb201nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb201nf_vctot(%rsp)

        decl nb201nf_innerk(%rsp)
        jz    _nb_kernel201nf_x86_64_sse.nb201nf_updateouterdata
        jmp   _nb_kernel201nf_x86_64_sse.nb201nf_odd_loop
_nb_kernel201nf_x86_64_sse.nb201nf_updateouterdata: 
        ## get n from stack
        movl nb201nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb201nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb201nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb201nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb201nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel201nf_x86_64_sse.nb201nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb201nf_n(%rsp)
        jmp _nb_kernel201nf_x86_64_sse.nb201nf_outer
_nb_kernel201nf_x86_64_sse.nb201nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb201nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel201nf_x86_64_sse.nb201nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel201nf_x86_64_sse.nb201nf_threadloop
_nb_kernel201nf_x86_64_sse.nb201nf_end: 

        movl nb201nf_nouter(%rsp),%eax
        movl nb201nf_ninner(%rsp),%ebx
        movq nb201nf_outeriter(%rbp),%rcx
        movq nb201nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

