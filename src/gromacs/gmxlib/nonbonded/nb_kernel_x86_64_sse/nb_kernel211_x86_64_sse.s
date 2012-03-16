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






.globl nb_kernel211_x86_64_sse
.globl _nb_kernel211_x86_64_sse
nb_kernel211_x86_64_sse:        
_nb_kernel211_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb211_fshift, 16
.set nb211_gid, 24
.set nb211_pos, 32
.set nb211_faction, 40
.set nb211_charge, 48
.set nb211_p_facel, 56
.set nb211_argkrf, 64
.set nb211_argcrf, 72
.set nb211_Vc, 80
.set nb211_type, 88
.set nb211_p_ntype, 96
.set nb211_vdwparam, 104
.set nb211_Vvdw, 112
.set nb211_p_tabscale, 120
.set nb211_VFtab, 128
.set nb211_invsqrta, 136
.set nb211_dvda, 144
.set nb211_p_gbtabscale, 152
.set nb211_GBtab, 160
.set nb211_p_nthreads, 168
.set nb211_count, 176
.set nb211_mtx, 184
.set nb211_outeriter, 192
.set nb211_inneriter, 200
.set nb211_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb211_ixO, 0
.set nb211_iyO, 16
.set nb211_izO, 32
.set nb211_ixH1, 48
.set nb211_iyH1, 64
.set nb211_izH1, 80
.set nb211_ixH2, 96
.set nb211_iyH2, 112
.set nb211_izH2, 128
.set nb211_iqO, 144
.set nb211_iqH, 160
.set nb211_dxO, 176
.set nb211_dyO, 192
.set nb211_dzO, 208
.set nb211_dxH1, 224
.set nb211_dyH1, 240
.set nb211_dzH1, 256
.set nb211_dxH2, 272
.set nb211_dyH2, 288
.set nb211_dzH2, 304
.set nb211_qqO, 320
.set nb211_qqH, 336
.set nb211_c6, 352
.set nb211_c12, 368
.set nb211_six, 384
.set nb211_twelve, 400
.set nb211_vctot, 416
.set nb211_Vvdwtot, 432
.set nb211_fixO, 448
.set nb211_fiyO, 464
.set nb211_fizO, 480
.set nb211_fixH1, 496
.set nb211_fiyH1, 512
.set nb211_fizH1, 528
.set nb211_fixH2, 544
.set nb211_fiyH2, 560
.set nb211_fizH2, 576
.set nb211_fjx, 592
.set nb211_fjy, 608
.set nb211_fjz, 624
.set nb211_half, 640
.set nb211_three, 656
.set nb211_two, 672
.set nb211_krf, 688
.set nb211_crf, 704
.set nb211_krsqO, 720
.set nb211_krsqH1, 736
.set nb211_krsqH2, 752
.set nb211_nri, 768
.set nb211_iinr, 776
.set nb211_jindex, 784
.set nb211_jjnr, 792
.set nb211_shift, 800
.set nb211_shiftvec, 808
.set nb211_facel, 816
.set nb211_innerjjnr, 824
.set nb211_is3, 832
.set nb211_ii3, 836
.set nb211_ntia, 840
.set nb211_innerk, 844
.set nb211_n, 848
.set nb211_nn1, 852
.set nb211_nouter, 856
.set nb211_ninner, 860

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $872,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb211_nouter(%rsp)
        movl %eax,nb211_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb211_nri(%rsp)
        movq %rsi,nb211_iinr(%rsp)
        movq %rdx,nb211_jindex(%rsp)
        movq %rcx,nb211_jjnr(%rsp)
        movq %r8,nb211_shift(%rsp)
        movq %r9,nb211_shiftvec(%rsp)
        movq nb211_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb211_facel(%rsp)


        movq nb211_argkrf(%rbp),%rsi
        movq nb211_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb211_krf(%rsp)
        movaps %xmm2,nb211_crf(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb211_half(%rsp)
        movss nb211_half(%rsp),%xmm1
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
        movaps %xmm1,nb211_half(%rsp)
        movaps %xmm2,nb211_two(%rsp)
        movaps %xmm3,nb211_three(%rsp)
        movaps %xmm4,nb211_six(%rsp)
        movaps %xmm5,nb211_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb211_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb211_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb211_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb211_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb211_iqO(%rsp)
        movaps %xmm4,nb211_iqH(%rsp)

        movq  nb211_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb211_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211_ntia(%rsp)

_nb_kernel211_x86_64_sse.nb211_threadloop: 
        movq  nb211_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel211_x86_64_sse.nb211_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211_x86_64_sse.nb211_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211_n(%rsp)
        movl %ebx,nb211_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211_x86_64_sse.nb211_outerstart
        jmp _nb_kernel211_x86_64_sse.nb211_end

_nb_kernel211_x86_64_sse.nb211_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211_nouter(%rsp),%ebx
        movl %ebx,nb211_nouter(%rsp)

_nb_kernel211_x86_64_sse.nb211_outer: 
        movq  nb211_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb211_is3(%rsp)      ## store is3 

        movq  nb211_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb211_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb211_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb211_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb211_ixO(%rsp)
        movaps %xmm4,nb211_iyO(%rsp)
        movaps %xmm5,nb211_izO(%rsp)

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
        movaps %xmm0,nb211_ixH1(%rsp)
        movaps %xmm1,nb211_iyH1(%rsp)
        movaps %xmm2,nb211_izH1(%rsp)
        movaps %xmm3,nb211_ixH2(%rsp)
        movaps %xmm4,nb211_iyH2(%rsp)
        movaps %xmm5,nb211_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb211_vctot(%rsp)
        movaps %xmm4,nb211_Vvdwtot(%rsp)
        movaps %xmm4,nb211_fixO(%rsp)
        movaps %xmm4,nb211_fiyO(%rsp)
        movaps %xmm4,nb211_fizO(%rsp)
        movaps %xmm4,nb211_fixH1(%rsp)
        movaps %xmm4,nb211_fiyH1(%rsp)
        movaps %xmm4,nb211_fizH1(%rsp)
        movaps %xmm4,nb211_fixH2(%rsp)
        movaps %xmm4,nb211_fiyH2(%rsp)
        movaps %xmm4,nb211_fizH2(%rsp)

        movq  nb211_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb211_pos(%rbp),%rsi
        movq  nb211_faction(%rbp),%rdi
        movq  nb211_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb211_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb211_ninner(%rsp),%ecx
        movl  %ecx,nb211_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb211_innerk(%rsp)      ## number of innerloop atoms 

        jge   _nb_kernel211_x86_64_sse.nb211_unroll_loop
        jmp   _nb_kernel211_x86_64_sse.nb211_odd_inner
_nb_kernel211_x86_64_sse.nb211_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb211_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb211_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb211_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb211_iqO(%rsp),%xmm3
        mulps  nb211_iqH(%rsp),%xmm4

        movaps  %xmm3,nb211_qqO(%rsp)
        movaps  %xmm4,nb211_qqH(%rsp)

        movq nb211_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb211_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb211_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        movlps (%rsi,%r8,4),%xmm6
        movlps (%rsi,%r10,4),%xmm7
        movhps (%rsi,%r9,4),%xmm6
        movhps (%rsi,%r11,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movaps %xmm4,nb211_c6(%rsp)
        movaps %xmm6,nb211_c12(%rsp)

        movq nb211_pos(%rbp),%rsi        ## base of pos[] 

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

    subps nb211_ixO(%rsp),%xmm0
    subps nb211_iyO(%rsp),%xmm1
    subps nb211_izO(%rsp),%xmm2
    subps nb211_ixH1(%rsp),%xmm3
    subps nb211_iyH1(%rsp),%xmm4
    subps nb211_izH1(%rsp),%xmm5
    subps nb211_ixH2(%rsp),%xmm6
    subps nb211_iyH2(%rsp),%xmm7
    subps nb211_izH2(%rsp),%xmm8

        movaps %xmm0,nb211_dxO(%rsp)
        movaps %xmm1,nb211_dyO(%rsp)
        movaps %xmm2,nb211_dzO(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb211_dxH1(%rsp)
        movaps %xmm4,nb211_dyH1(%rsp)
        movaps %xmm5,nb211_dzH1(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb211_dxH2(%rsp)
        movaps %xmm7,nb211_dyH2(%rsp)
        movaps %xmm8,nb211_dzH2(%rsp)
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

        movaps  nb211_three(%rsp),%xmm9
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

        movaps  nb211_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvO
        mulps   %xmm4,%xmm10 ## rinvH1
    mulps   %xmm4,%xmm11 ## rinvH2

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb211_krf(%rsp),%xmm2
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
    movaps %xmm9,%xmm12
    mulps  %xmm12,%xmm12 ## rinv4
    mulps  %xmm9,%xmm12 ## rinv6
    subps  nb211_crf(%rsp),%xmm2     ## rinv+krsq-crf
    subps  nb211_crf(%rsp),%xmm5
    subps  nb211_crf(%rsp),%xmm8
    mulps  nb211_qqO(%rsp),%xmm2   ## voul=qq*(rinv+ krsq-crf)
    mulps  nb211_qqH(%rsp),%xmm5   ## voul=qq*(rinv+ krsq-crf)
    mulps  nb211_qqH(%rsp),%xmm8   ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    movaps %xmm12,%xmm13 ## rinv6
    mulps %xmm12,%xmm12 ## rinv12
        mulps  nb211_c6(%rsp),%xmm13
        mulps  nb211_c12(%rsp),%xmm12
    movaps %xmm12,%xmm14
    subps  %xmm13,%xmm14
    mulps  nb211_qqO(%rsp),%xmm1     ## (rinv-2*krsq)*qq
    mulps  nb211_qqH(%rsp),%xmm4
    mulps  nb211_qqH(%rsp),%xmm7
    addps  nb211_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,nb211_vctot(%rsp)

        addps  nb211_Vvdwtot(%rsp),%xmm14
        mulps  nb211_six(%rsp),%xmm13
        mulps  nb211_twelve(%rsp),%xmm12
        movaps %xmm14,nb211_Vvdwtot(%rsp)
    subps  %xmm13,%xmm12 ## LJ fscal        

    addps %xmm12,%xmm1

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        movq  nb211_faction(%rbp),%rdi
        ## move j  forces to local temp variables 
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

        mulps nb211_dxO(%rsp),%xmm0
        mulps nb211_dyO(%rsp),%xmm1
        mulps nb211_dzO(%rsp),%xmm2
        mulps nb211_dxH1(%rsp),%xmm3
        mulps nb211_dyH1(%rsp),%xmm4
        mulps nb211_dzH1(%rsp),%xmm5
        mulps nb211_dxH2(%rsp),%xmm6
        mulps nb211_dyH2(%rsp),%xmm7
        mulps nb211_dzH2(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb211_fixO(%rsp),%xmm0
    addps nb211_fiyO(%rsp),%xmm1
    addps nb211_fizO(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb211_fixH1(%rsp),%xmm3
    addps nb211_fiyH1(%rsp),%xmm4
    addps nb211_fizH1(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb211_fixH2(%rsp),%xmm6
    addps nb211_fiyH2(%rsp),%xmm7
    addps nb211_fizH2(%rsp),%xmm8

    movaps %xmm0,nb211_fixO(%rsp)
    movaps %xmm1,nb211_fiyO(%rsp)
    movaps %xmm2,nb211_fizO(%rsp)
    movaps %xmm3,nb211_fixH1(%rsp)
    movaps %xmm4,nb211_fiyH1(%rsp)
    movaps %xmm5,nb211_fizH1(%rsp)
    movaps %xmm6,nb211_fixH2(%rsp)
    movaps %xmm7,nb211_fiyH2(%rsp)
    movaps %xmm8,nb211_fizH2(%rsp)

    ## xmm9 = fjx
    ## xmm10 = fjy
    ## xmm11 = fjz
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12

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
        subl $4,nb211_innerk(%rsp)
        jl    _nb_kernel211_x86_64_sse.nb211_odd_inner
        jmp   _nb_kernel211_x86_64_sse.nb211_unroll_loop
_nb_kernel211_x86_64_sse.nb211_odd_inner: 
        addl $4,nb211_innerk(%rsp)
        jnz   _nb_kernel211_x86_64_sse.nb211_odd_loop
        jmp   _nb_kernel211_x86_64_sse.nb211_updateouterdata
_nb_kernel211_x86_64_sse.nb211_odd_loop: 
        movq  nb211_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb211_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb211_iqO(%rsp),%xmm4
        movq nb211_charge(%rbp),%rsi
        movhps nb211_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb211_qqO(%rsp)    ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb211_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb211_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb211_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb211_c6(%rsp)
        movaps %xmm7,nb211_c12(%rsp)
        movq nb211_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb211_ixO(%rsp),%xmm0
        movss nb211_iyO(%rsp),%xmm1
        movss nb211_izO(%rsp),%xmm2

        movlps nb211_ixH1(%rsp),%xmm6
        movlps nb211_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb211_iyH1(%rsp),%xmm6
        movlps nb211_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb211_izH1(%rsp),%xmm6
        movlps nb211_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        movaps %xmm3,nb211_dxO(%rsp)
        movaps %xmm4,nb211_dyO(%rsp)
        movaps %xmm5,nb211_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb211_krf(%rsp),%xmm0
        movaps %xmm0,nb211_krsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb211_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb211_half(%rsp),%xmm0
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
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211_c6(%rsp),%xmm1
        mulps  nb211_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb211_Vvdwtot(%rsp),%xmm5
        mulss  nb211_six(%rsp),%xmm1
        mulss  nb211_twelve(%rsp),%xmm2
        subss  %xmm1,%xmm2

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb211_krsqO(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        mulps  nb211_two(%rsp),%xmm3
        subps  nb211_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb211_qqO(%rsp),%xmm0    ## xmm0=vcoul 
        mulps  nb211_qqO(%rsp),%xmm1    ## xmm1=coul part of fs 

        addps %xmm1,%xmm2       ## total fs 

        mulps  %xmm2,%xmm4      ## xmm4=total fscal 
        addps  nb211_vctot(%rsp),%xmm0
        movaps %xmm0,nb211_vctot(%rsp)

        movaps nb211_dxO(%rsp),%xmm0
        movaps nb211_dyO(%rsp),%xmm1
        movaps nb211_dzO(%rsp),%xmm2

        movaps %xmm5,nb211_Vvdwtot(%rsp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb211_fixO(%rsp),%xmm3
        movss  nb211_fiyO(%rsp),%xmm4
        movss  nb211_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb211_fixO(%rsp)
        movss  %xmm4,nb211_fiyO(%rsp)
        movss  %xmm5,nb211_fizO(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb211_fixH1(%rsp),%xmm3
        addss  nb211_fiyH1(%rsp),%xmm4
        addss  nb211_fizH1(%rsp),%xmm5
        movss  %xmm3,nb211_fixH1(%rsp)
        movss  %xmm4,nb211_fiyH1(%rsp)
        movss  %xmm5,nb211_fizH1(%rsp)          ## updated the H1 force 

        movq nb211_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb211_fixH2(%rsp),%xmm3
        addss  nb211_fiyH2(%rsp),%xmm4
        addss  nb211_fizH2(%rsp),%xmm5
        movss  %xmm3,nb211_fixH2(%rsp)
        movss  %xmm4,nb211_fiyH2(%rsp)
        movss  %xmm5,nb211_fizH2(%rsp)          ## updated the H2 force 

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

        decl nb211_innerk(%rsp)
        jz    _nb_kernel211_x86_64_sse.nb211_updateouterdata
        jmp   _nb_kernel211_x86_64_sse.nb211_odd_loop
_nb_kernel211_x86_64_sse.nb211_updateouterdata: 
        movl  nb211_ii3(%rsp),%ecx
        movq  nb211_faction(%rbp),%rdi
        movq  nb211_fshift(%rbp),%rsi
        movl  nb211_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb211_fixO(%rsp),%xmm0
        movaps nb211_fiyO(%rsp),%xmm1
        movaps nb211_fizO(%rsp),%xmm2

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
        movaps nb211_fixH1(%rsp),%xmm0
        movaps nb211_fiyH1(%rsp),%xmm1
        movaps nb211_fizH1(%rsp),%xmm2

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
        movaps nb211_fixH2(%rsp),%xmm0
        movaps nb211_fiyH2(%rsp),%xmm1
        movaps nb211_fizH2(%rsp),%xmm2

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
        movl nb211_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb211_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb211_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb211_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb211_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb211_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb211_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211_x86_64_sse.nb211_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211_n(%rsp)
        jmp _nb_kernel211_x86_64_sse.nb211_outer
_nb_kernel211_x86_64_sse.nb211_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211_x86_64_sse.nb211_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211_x86_64_sse.nb211_threadloop
_nb_kernel211_x86_64_sse.nb211_end: 
        movl nb211_nouter(%rsp),%eax
        movl nb211_ninner(%rsp),%ebx
        movq nb211_outeriter(%rbp),%rcx
        movq nb211_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $872,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel211nf_x86_64_sse
.globl _nb_kernel211nf_x86_64_sse
nb_kernel211nf_x86_64_sse:      
_nb_kernel211nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb211nf_fshift, 16
.set nb211nf_gid, 24
.set nb211nf_pos, 32
.set nb211nf_faction, 40
.set nb211nf_charge, 48
.set nb211nf_p_facel, 56
.set nb211nf_argkrf, 64
.set nb211nf_argcrf, 72
.set nb211nf_Vc, 80
.set nb211nf_type, 88
.set nb211nf_p_ntype, 96
.set nb211nf_vdwparam, 104
.set nb211nf_Vvdw, 112
.set nb211nf_p_tabscale, 120
.set nb211nf_VFtab, 128
.set nb211nf_invsqrta, 136
.set nb211nf_dvda, 144
.set nb211nf_p_gbtabscale, 152
.set nb211nf_GBtab, 160
.set nb211nf_p_nthreads, 168
.set nb211nf_count, 176
.set nb211nf_mtx, 184
.set nb211nf_outeriter, 192
.set nb211nf_inneriter, 200
.set nb211nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb211nf_ixO, 0
.set nb211nf_iyO, 16
.set nb211nf_izO, 32
.set nb211nf_ixH1, 48
.set nb211nf_iyH1, 64
.set nb211nf_izH1, 80
.set nb211nf_ixH2, 96
.set nb211nf_iyH2, 112
.set nb211nf_izH2, 128
.set nb211nf_iqO, 144
.set nb211nf_iqH, 160
.set nb211nf_qqO, 176
.set nb211nf_qqH, 192
.set nb211nf_c6, 208
.set nb211nf_c12, 224
.set nb211nf_vctot, 240
.set nb211nf_Vvdwtot, 256
.set nb211nf_half, 272
.set nb211nf_three, 288
.set nb211nf_krf, 304
.set nb211nf_crf, 320
.set nb211nf_krsqO, 336
.set nb211nf_krsqH1, 352
.set nb211nf_krsqH2, 368
.set nb211nf_nri, 384
.set nb211nf_iinr, 392
.set nb211nf_jindex, 400
.set nb211nf_jjnr, 408
.set nb211nf_shift, 416
.set nb211nf_shiftvec, 424
.set nb211nf_facel, 432
.set nb211nf_innerjjnr, 440
.set nb211nf_ntia, 448
.set nb211nf_is3, 452
.set nb211nf_ii3, 456
.set nb211nf_innerk, 460
.set nb211nf_n, 464
.set nb211nf_nn1, 468
.set nb211nf_nouter, 472
.set nb211nf_ninner, 476

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $488,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb211nf_nouter(%rsp)
        movl %eax,nb211nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb211nf_nri(%rsp)
        movq %rsi,nb211nf_iinr(%rsp)
        movq %rdx,nb211nf_jindex(%rsp)
        movq %rcx,nb211nf_jjnr(%rsp)
        movq %r8,nb211nf_shift(%rsp)
        movq %r9,nb211nf_shiftvec(%rsp)
        movq nb211nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb211nf_facel(%rsp)


        movq nb211nf_argkrf(%rbp),%rsi
        movq nb211nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb211nf_krf(%rsp)
        movaps %xmm2,nb211nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb211nf_half(%rsp)
        movss nb211nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb211nf_half(%rsp)
        movaps %xmm3,nb211nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb211nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb211nf_charge(%rbp),%rdx
        movss (%rdx,%rbx,4),%xmm3
        movss 4(%rdx,%rbx,4),%xmm4
        movq nb211nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb211nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb211nf_iqO(%rsp)
        movaps %xmm4,nb211nf_iqH(%rsp)

        movq  nb211nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb211nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb211nf_ntia(%rsp)

_nb_kernel211nf_x86_64_sse.nb211nf_threadloop: 
        movq  nb211nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel211nf_x86_64_sse.nb211nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel211nf_x86_64_sse.nb211nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb211nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb211nf_n(%rsp)
        movl %ebx,nb211nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel211nf_x86_64_sse.nb211nf_outerstart
        jmp _nb_kernel211nf_x86_64_sse.nb211nf_end

_nb_kernel211nf_x86_64_sse.nb211nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb211nf_nouter(%rsp),%ebx
        movl %ebx,nb211nf_nouter(%rsp)

_nb_kernel211nf_x86_64_sse.nb211nf_outer: 
        movq  nb211nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb211nf_is3(%rsp)            ## store is3 

        movq  nb211nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb211nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb211nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb211nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3
        addss 4(%rax,%rbx,4),%xmm4
        addss 8(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb211nf_ixO(%rsp)
        movaps %xmm4,nb211nf_iyO(%rsp)
        movaps %xmm5,nb211nf_izO(%rsp)

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
        movaps %xmm0,nb211nf_ixH1(%rsp)
        movaps %xmm1,nb211nf_iyH1(%rsp)
        movaps %xmm2,nb211nf_izH1(%rsp)
        movaps %xmm3,nb211nf_ixH2(%rsp)
        movaps %xmm4,nb211nf_iyH2(%rsp)
        movaps %xmm5,nb211nf_izH2(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb211nf_vctot(%rsp)
        movaps %xmm4,nb211nf_Vvdwtot(%rsp)

        movq  nb211nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb211nf_pos(%rbp),%rsi
        movq  nb211nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb211nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb211nf_ninner(%rsp),%ecx
        movl  %ecx,nb211nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb211nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel211nf_x86_64_sse.nb211nf_unroll_loop
        jmp   _nb_kernel211nf_x86_64_sse.nb211nf_odd_inner
_nb_kernel211nf_x86_64_sse.nb211nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb211nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx           ## eax-edx=jnr1-4 

        addq $16,nb211nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb211nf_charge(%rbp),%rsi     ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4           ## and in xmm4 
        mulps  nb211nf_iqO(%rsp),%xmm3
        mulps  nb211nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb211nf_qqO(%rsp)
        movaps  %xmm4,nb211nf_qqH(%rsp)

        movq nb211nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb211nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb211nf_ntia(%rsp),%edi
        addq %rdi,%rax
        addq %rdi,%rbx
        addq %rdi,%rcx
        addq %rdi,%rdx

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

        movaps %xmm4,nb211nf_c6(%rsp)
        movaps %xmm6,nb211nf_c12(%rsp)

        movq nb211nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movaps nb211nf_ixO(%rsp),%xmm4
        movaps nb211nf_iyO(%rsp),%xmm5
        movaps nb211nf_izO(%rsp),%xmm6

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
        movaps nb211nf_ixH1(%rsp),%xmm4
        movaps nb211nf_iyH1(%rsp),%xmm5
        movaps nb211nf_izH1(%rsp),%xmm6

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
        movaps nb211nf_ixH2(%rsp),%xmm3
        movaps nb211nf_iyH2(%rsp),%xmm4
        movaps nb211nf_izH2(%rsp),%xmm5

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

        mulps  nb211nf_krf(%rsp),%xmm0
        mulps  nb211nf_krf(%rsp),%xmm1
        mulps  nb211nf_krf(%rsp),%xmm2

        movaps %xmm0,nb211nf_krsqH2(%rsp)
        movaps %xmm1,nb211nf_krsqH1(%rsp)
        movaps %xmm2,nb211nf_krsqO(%rsp)

        ## start with rsqO - seed in xmm2       
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvO in xmm7 
        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH1 in xmm6 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb211nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb211nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do O interactions 
        movaps  %xmm7,%xmm4
        mulps   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        movaps %xmm4,%xmm1
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211nf_c6(%rsp),%xmm1
        mulps  nb211nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm3
        subps  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addps  nb211nf_Vvdwtot(%rsp),%xmm3

        movaps %xmm7,%xmm0
        movaps nb211nf_krsqO(%rsp),%xmm1
        addps  %xmm1,%xmm0
        subps  nb211nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        subps  %xmm1,%xmm7
        mulps  nb211nf_qqO(%rsp),%xmm0
        addps  nb211nf_vctot(%rsp),%xmm0
        movaps %xmm3,nb211nf_Vvdwtot(%rsp)
        movaps %xmm0,nb211nf_vctot(%rsp)

        ## H1 interactions 
        movaps  nb211nf_krsqH1(%rsp),%xmm0
        addps   %xmm0,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb211nf_crf(%rsp),%xmm6
        mulps   nb211nf_qqH(%rsp),%xmm6   ## vcoul 
        addps  nb211nf_vctot(%rsp),%xmm6
        movaps %xmm6,nb211nf_vctot(%rsp)

        ## H2 interactions 
        movaps  %xmm5,%xmm7 ## rinv 
        movaps  nb211nf_krsqH2(%rsp),%xmm0
        addps   %xmm0,%xmm5     ## xmm5=rinv+ krsq 
        subps   nb211nf_crf(%rsp),%xmm5
        mulps   nb211nf_qqH(%rsp),%xmm5   ## vcoul 
        addps   nb211nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb211nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb211nf_innerk(%rsp)
        jl    _nb_kernel211nf_x86_64_sse.nb211nf_odd_inner
        jmp   _nb_kernel211nf_x86_64_sse.nb211nf_unroll_loop
_nb_kernel211nf_x86_64_sse.nb211nf_odd_inner: 
        addl $4,nb211nf_innerk(%rsp)
        jnz   _nb_kernel211nf_x86_64_sse.nb211nf_odd_loop
        jmp   _nb_kernel211nf_x86_64_sse.nb211nf_updateouterdata
_nb_kernel211nf_x86_64_sse.nb211nf_odd_loop: 
        movq  nb211nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb211nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb211nf_iqO(%rsp),%xmm4
        movq nb211nf_charge(%rbp),%rsi
        movhps nb211nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb211nf_qqO(%rsp)          ## use oxygen qq for storage 

        xorps %xmm6,%xmm6
        movq nb211nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb211nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb211nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb211nf_c6(%rsp)
        movaps %xmm7,nb211nf_c12(%rsp)

        movq nb211nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb211nf_ixO(%rsp),%xmm3
        movss nb211nf_iyO(%rsp),%xmm4
        movss nb211nf_izO(%rsp),%xmm5

        movlps nb211nf_ixH1(%rsp),%xmm6
        movlps nb211nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb211nf_iyH1(%rsp),%xmm6
        movlps nb211nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb211nf_izH1(%rsp),%xmm6
        movlps nb211nf_izH2(%rsp),%xmm7
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
        mulps nb211nf_krf(%rsp),%xmm0
        movaps %xmm0,nb211nf_krsqO(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb211nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb211nf_half(%rsp),%xmm0
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
        movaps %xmm4,%xmm1
        mulss  %xmm4,%xmm1
        mulss  %xmm4,%xmm1      ## xmm1=rinvsix 
        movaps %xmm1,%xmm2
        mulss  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulps  nb211nf_c6(%rsp),%xmm1
        mulps  nb211nf_c12(%rsp),%xmm2
        movaps %xmm2,%xmm5
        subss  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addps  nb211nf_Vvdwtot(%rsp),%xmm5
        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb211nf_krsqO(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb211nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb211nf_qqO(%rsp),%xmm0          ## xmm0=vcoul 
        addps  nb211nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb211nf_vctot(%rsp)
        movaps %xmm5,nb211nf_Vvdwtot(%rsp)

        decl nb211nf_innerk(%rsp)
        jz    _nb_kernel211nf_x86_64_sse.nb211nf_updateouterdata
        jmp   _nb_kernel211nf_x86_64_sse.nb211nf_odd_loop
_nb_kernel211nf_x86_64_sse.nb211nf_updateouterdata: 
        ## get n from stack
        movl nb211nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb211nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb211nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb211nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb211nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb211nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb211nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel211nf_x86_64_sse.nb211nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb211nf_n(%rsp)
        jmp _nb_kernel211nf_x86_64_sse.nb211nf_outer
_nb_kernel211nf_x86_64_sse.nb211nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb211nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel211nf_x86_64_sse.nb211nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel211nf_x86_64_sse.nb211nf_threadloop
_nb_kernel211nf_x86_64_sse.nb211nf_end: 
        movl nb211nf_nouter(%rsp),%eax
        movl nb211nf_ninner(%rsp),%ebx
        movq nb211nf_outeriter(%rbp),%rcx
        movq nb211nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $488,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



