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





.globl nb_kernel313_x86_64_sse
.globl _nb_kernel313_x86_64_sse
nb_kernel313_x86_64_sse:        
_nb_kernel313_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb313_fshift, 16
.set nb313_gid, 24
.set nb313_pos, 32
.set nb313_faction, 40
.set nb313_charge, 48
.set nb313_p_facel, 56
.set nb313_argkrf, 64
.set nb313_argcrf, 72
.set nb313_Vc, 80
.set nb313_type, 88
.set nb313_p_ntype, 96
.set nb313_vdwparam, 104
.set nb313_Vvdw, 112
.set nb313_p_tabscale, 120
.set nb313_VFtab, 128
.set nb313_invsqrta, 136
.set nb313_dvda, 144
.set nb313_p_gbtabscale, 152
.set nb313_GBtab, 160
.set nb313_p_nthreads, 168
.set nb313_count, 176
.set nb313_mtx, 184
.set nb313_outeriter, 192
.set nb313_inneriter, 200
.set nb313_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb313_ixO, 0
.set nb313_iyO, 16
.set nb313_izO, 32
.set nb313_ixH1, 48
.set nb313_iyH1, 64
.set nb313_izH1, 80
.set nb313_ixH2, 96
.set nb313_iyH2, 112
.set nb313_izH2, 128
.set nb313_ixM, 144
.set nb313_iyM, 160
.set nb313_izM, 176
.set nb313_iqM, 192
.set nb313_iqH, 208
.set nb313_epsH1, 224
.set nb313_epsH2, 240
.set nb313_epsM, 256
.set nb313_dxH1, 272
.set nb313_dyH1, 288
.set nb313_dzH1, 304
.set nb313_dxH2, 320
.set nb313_dyH2, 336
.set nb313_dzH2, 352
.set nb313_dxM, 368
.set nb313_dyM, 384
.set nb313_dzM, 400
.set nb313_qqM, 416
.set nb313_qqH, 432
.set nb313_rinvH1, 448
.set nb313_rinvH2, 464
.set nb313_rinvM, 480
.set nb313_rH1, 496
.set nb313_rH2, 512
.set nb313_rM, 528
.set nb313_tsc, 544
.set nb313_two, 560
.set nb313_c6, 576
.set nb313_c12, 592
.set nb313_six, 608
.set nb313_twelve, 624
.set nb313_vctot, 640
.set nb313_Vvdwtot, 656
.set nb313_fixO, 672
.set nb313_fiyO, 688
.set nb313_fizO, 704
.set nb313_fixH1, 720
.set nb313_fiyH1, 736
.set nb313_fizH1, 752
.set nb313_fixH2, 768
.set nb313_fiyH2, 784
.set nb313_fizH2, 800
.set nb313_fixM, 816
.set nb313_fiyM, 832
.set nb313_fizM, 848
.set nb313_fjx, 864
.set nb313_fjy, 880
.set nb313_fjz, 896
.set nb313_half, 912
.set nb313_three, 928
.set nb313_is3, 944
.set nb313_ii3, 948
.set nb313_nri, 952
.set nb313_iinr, 960
.set nb313_jindex, 968
.set nb313_jjnr, 976
.set nb313_shift, 984
.set nb313_shiftvec, 992
.set nb313_facel, 1000
.set nb313_innerjjnr, 1008
.set nb313_ntia, 1016
.set nb313_innerk, 1020
.set nb313_n, 1024
.set nb313_nn1, 1028
.set nb313_nouter, 1032
.set nb313_ninner, 1036
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1048,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb313_nouter(%rsp)
        movl %eax,nb313_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb313_nri(%rsp)
        movq %rsi,nb313_iinr(%rsp)
        movq %rdx,nb313_jindex(%rsp)
        movq %rcx,nb313_jjnr(%rsp)
        movq %r8,nb313_shift(%rsp)
        movq %r9,nb313_shiftvec(%rsp)
        movq nb313_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb313_facel(%rsp)

        movq nb313_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb313_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb313_half(%rsp)
        movss nb313_half(%rsp),%xmm1
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
        movaps %xmm1,nb313_half(%rsp)
        movaps %xmm2,nb313_two(%rsp)
        movaps %xmm3,nb313_three(%rsp)
        movaps %xmm4,nb313_six(%rsp)
        movaps %xmm5,nb313_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb313_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb313_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb313_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb313_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb313_iqM(%rsp)
        movaps %xmm4,nb313_iqH(%rsp)

        movq  nb313_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb313_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313_ntia(%rsp)
_nb_kernel313_x86_64_sse.nb313_threadloop: 
        movq  nb313_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel313_x86_64_sse.nb313_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313_x86_64_sse.nb313_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313_n(%rsp)
        movl %ebx,nb313_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313_x86_64_sse.nb313_outerstart
        jmp _nb_kernel313_x86_64_sse.nb313_end

_nb_kernel313_x86_64_sse.nb313_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313_nouter(%rsp),%ebx
        movl %ebx,nb313_nouter(%rsp)

_nb_kernel313_x86_64_sse.nb313_outer: 
        movq  nb313_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb313_is3(%rsp)      ## store is3 

        movq  nb313_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb313_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb313_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb313_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb313_ixO(%rsp)
        movaps %xmm4,nb313_iyO(%rsp)
        movaps %xmm5,nb313_izO(%rsp)
        movaps %xmm6,nb313_ixH1(%rsp)
        movaps %xmm7,nb313_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb313_izH1(%rsp)
        movaps %xmm0,nb313_ixH2(%rsp)
        movaps %xmm1,nb313_iyH2(%rsp)
        movaps %xmm2,nb313_izH2(%rsp)
        movaps %xmm3,nb313_ixM(%rsp)
        movaps %xmm4,nb313_iyM(%rsp)
        movaps %xmm5,nb313_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb313_vctot(%rsp)
        movaps %xmm4,nb313_Vvdwtot(%rsp)
        movaps %xmm4,nb313_fixO(%rsp)
        movaps %xmm4,nb313_fiyO(%rsp)
        movaps %xmm4,nb313_fizO(%rsp)
        movaps %xmm4,nb313_fixH1(%rsp)
        movaps %xmm4,nb313_fiyH1(%rsp)
        movaps %xmm4,nb313_fizH1(%rsp)
        movaps %xmm4,nb313_fixH2(%rsp)
        movaps %xmm4,nb313_fiyH2(%rsp)
        movaps %xmm4,nb313_fizH2(%rsp)
        movaps %xmm4,nb313_fixM(%rsp)
        movaps %xmm4,nb313_fiyM(%rsp)
        movaps %xmm4,nb313_fizM(%rsp)

        movq  nb313_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb313_pos(%rbp),%rsi
        movq  nb313_faction(%rbp),%rdi
        movq  nb313_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb313_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb313_ninner(%rsp),%ecx
        movl  %ecx,nb313_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb313_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel313_x86_64_sse.nb313_unroll_loop
        jmp   _nb_kernel313_x86_64_sse.nb313_odd_inner
_nb_kernel313_x86_64_sse.nb313_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb313_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb313_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb313_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb313_iqM(%rsp),%xmm3
        mulps  nb313_iqH(%rsp),%xmm4

        movaps  %xmm3,nb313_qqM(%rsp)
        movaps  %xmm4,nb313_qqH(%rsp)

        movq nb313_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb313_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb313_ntia(%rsp),%edi
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

        movaps %xmm4,nb313_c6(%rsp)
        movaps %xmm6,nb313_c12(%rsp)

        movq nb313_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
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

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    ## O interaction
    ## copy to xmm3-xmm5
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5

    subps nb313_ixO(%rsp),%xmm3
    subps nb313_iyO(%rsp),%xmm4
    subps nb313_izO(%rsp),%xmm5

    movaps %xmm3,%xmm13
    movaps %xmm4,%xmm14
    movaps %xmm5,%xmm15

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3

    ## calc 1/rsq
    rcpps %xmm3,%xmm5
    movaps nb313_two(%rsp),%xmm4
    mulps %xmm5,%xmm3
    subps %xmm3,%xmm4
    mulps %xmm5,%xmm4       ## xmm4=rinvsq

    movaps %xmm4,%xmm3      ## rinvsq
    mulps  %xmm4,%xmm4      ## rinv4
    mulps  %xmm3,%xmm4      ## rinv6
    movaps %xmm4,%xmm5
    mulps  %xmm5,%xmm5      ## rinv12
    mulps  nb313_c6(%rsp),%xmm4
    mulps  nb313_c12(%rsp),%xmm5
    movaps %xmm5,%xmm6
    subps  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulps  nb313_six(%rsp),%xmm4
    mulps  nb313_twelve(%rsp),%xmm5
    subps  %xmm4,%xmm5
    mulps  %xmm5,%xmm3  ## fscal

    addps  nb313_Vvdwtot(%rsp),%xmm6
    movaps %xmm6,nb313_Vvdwtot(%rsp)

    mulps  %xmm3,%xmm13 ## fx
    mulps  %xmm3,%xmm14 ## fy
    mulps  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movaps %xmm13,nb313_fjx(%rsp)
    movaps %xmm14,nb313_fjy(%rsp)
    movaps %xmm15,nb313_fjz(%rsp)

    ## increment i O force
    addps nb313_fixO(%rsp),%xmm13
    addps nb313_fiyO(%rsp),%xmm14
    addps nb313_fizO(%rsp),%xmm15
    movaps %xmm13,nb313_fixO(%rsp)
    movaps %xmm14,nb313_fiyO(%rsp)
    movaps %xmm15,nb313_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.        
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb313_ixH1(%rsp),%xmm0
    subps nb313_iyH1(%rsp),%xmm1
    subps nb313_izH1(%rsp),%xmm2
    subps nb313_ixH2(%rsp),%xmm3
    subps nb313_iyH2(%rsp),%xmm4
    subps nb313_izH2(%rsp),%xmm5
    subps nb313_ixM(%rsp),%xmm6
    subps nb313_iyM(%rsp),%xmm7
    subps nb313_izM(%rsp),%xmm8

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps %xmm0,nb313_dxH1(%rsp)
        movaps %xmm1,nb313_dyH1(%rsp)
        movaps %xmm2,nb313_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb313_dxH2(%rsp)
        movaps %xmm4,nb313_dyH2(%rsp)
        movaps %xmm5,nb313_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb313_dxM(%rsp)
        movaps %xmm7,nb313_dyM(%rsp)
        movaps %xmm8,nb313_dzM(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb313_three(%rsp),%xmm9
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

        movaps  nb313_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1 
        mulps   %xmm4,%xmm10 ## rinvH2
    mulps   %xmm4,%xmm11 ## rinvM

        movaps  %xmm9,nb313_rinvH1(%rsp)
        movaps  %xmm10,nb313_rinvH2(%rsp)
        movaps  %xmm11,nb313_rinvM(%rsp)

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps nb313_tsc(%rsp),%xmm1
    mulps  %xmm9,%xmm0 ## r
    mulps  %xmm10,%xmm3
    mulps  %xmm11,%xmm6
    mulps  %xmm1,%xmm0 ## rtab
    mulps  %xmm1,%xmm3
    mulps  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttps2dq %xmm0,%xmm1
    cvttps2dq %xmm3,%xmm4
    cvttps2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2ps  %xmm1,%xmm2
    cvtdq2ps  %xmm4,%xmm5
    cvtdq2ps  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    movhlps %xmm1,%xmm13
    movhlps %xmm4,%xmm14
    movhlps %xmm7,%xmm15
    movd    %xmm1,%eax
    movd    %xmm4,%r8d
    movd    %xmm7,%r12d
    movd    %xmm13,%ecx
    movd    %xmm14,%r10d
    movd    %xmm15,%r14d
    pshufd $1,%xmm1,%xmm1
    pshufd $1,%xmm4,%xmm4
    pshufd $1,%xmm7,%xmm7
    pshufd $1,%xmm13,%xmm13
    pshufd $1,%xmm14,%xmm14
    pshufd $1,%xmm15,%xmm15
    movd    %xmm1,%ebx
    movd    %xmm4,%r9d
    movd    %xmm7,%r13d
    movd    %xmm13,%edx
    movd    %xmm14,%r11d
    movd    %xmm15,%r15d

    movq nb313_VFtab(%rbp),%rsi

    ## calculate eps
    subps     %xmm2,%xmm0
    subps     %xmm5,%xmm3
    subps     %xmm8,%xmm6

    movaps    %xmm0,nb313_epsH1(%rsp)
    movaps    %xmm3,nb313_epsH2(%rsp)
    movaps    %xmm6,nb313_epsM(%rsp)

    ## Load LOTS of table data
        movlps (%rsi,%rax,4),%xmm1
        movlps (%rsi,%r8,4),%xmm5
        movlps (%rsi,%r12,4),%xmm9

        movlps (%rsi,%rcx,4),%xmm3
        movlps (%rsi,%r10,4),%xmm7
        movlps (%rsi,%r14,4),%xmm11

        movhps (%rsi,%rbx,4),%xmm1
        movhps (%rsi,%r9,4),%xmm5
        movhps (%rsi,%r13,4),%xmm9

        movhps (%rsi,%rdx,4),%xmm3
        movhps (%rsi,%r11,4),%xmm7
        movhps (%rsi,%r15,4),%xmm11

    movaps %xmm1,%xmm0
    movaps %xmm5,%xmm4
    movaps %xmm9,%xmm8
        shufps $136,%xmm3,%xmm0 ## 10001000
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $136,%xmm11,%xmm8 ## 10001000
        shufps $221,%xmm3,%xmm1 ## 11011101
        shufps $221,%xmm7,%xmm5 ## 11011101
        shufps $221,%xmm11,%xmm9 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm3
        movlps 8(%rsi,%r8,4),%xmm7
        movlps 8(%rsi,%r12,4),%xmm11

        movlps 8(%rsi,%rcx,4),%xmm12
        movlps 8(%rsi,%r10,4),%xmm13
        movlps 8(%rsi,%r14,4),%xmm14

        movhps 8(%rsi,%rbx,4),%xmm3
        movhps 8(%rsi,%r9,4),%xmm7
        movhps 8(%rsi,%r13,4),%xmm11

        movhps 8(%rsi,%rdx,4),%xmm12
        movhps 8(%rsi,%r11,4),%xmm13
        movhps 8(%rsi,%r15,4),%xmm14

    movaps %xmm3,%xmm2
    movaps %xmm7,%xmm6
    movaps %xmm11,%xmm10

        shufps $136,%xmm12,%xmm2 ## 10001000
        shufps $136,%xmm13,%xmm6 ## 10001000
        shufps $136,%xmm14,%xmm10 ## 10001000
        shufps $221,%xmm12,%xmm3 ## 11011101
        shufps $221,%xmm13,%xmm7 ## 11011101
        shufps $221,%xmm14,%xmm11 ## 11011101
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movaps nb313_epsH1(%rsp),%xmm12
    movaps nb313_epsH2(%rsp),%xmm13
    movaps nb313_epsM(%rsp),%xmm14

    mulps  %xmm12,%xmm3  ## Heps
    mulps  %xmm13,%xmm7
    mulps  %xmm14,%xmm11
    mulps  %xmm12,%xmm2  ## Geps
    mulps  %xmm13,%xmm6
    mulps  %xmm14,%xmm10
    mulps  %xmm12,%xmm3  ## Heps2
    mulps  %xmm13,%xmm7
    mulps  %xmm14,%xmm11

    addps  %xmm2,%xmm1  ## F+Geps
    addps  %xmm6,%xmm5
    addps  %xmm10,%xmm9
    addps  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addps  %xmm7,%xmm5
    addps  %xmm11,%xmm9
    addps  %xmm3,%xmm3   ## 2*Heps2
    addps  %xmm7,%xmm7
    addps  %xmm11,%xmm11
    addps  %xmm2,%xmm3   ## 2*Heps2+Geps
    addps  %xmm6,%xmm7
    addps  %xmm10,%xmm11
    addps  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addps  %xmm5,%xmm7
    addps  %xmm9,%xmm11
    mulps  %xmm12,%xmm1  ## eps*Fp
    mulps  %xmm13,%xmm5
    mulps  %xmm14,%xmm9
    movaps nb313_qqH(%rsp),%xmm12
    movaps nb313_qqM(%rsp),%xmm13
    addps  %xmm0,%xmm1    ## VV
    addps  %xmm4,%xmm5
    addps  %xmm8,%xmm9
    mulps  %xmm12,%xmm1  ## VV*qq = vcoul
    mulps  %xmm12,%xmm5
    mulps  %xmm13,%xmm9
    mulps  %xmm12,%xmm3   ## FF*qq = fij
    mulps  %xmm12,%xmm7
    mulps  %xmm13,%xmm11

    ## accumulate vctot
    addps  nb313_vctot(%rsp),%xmm1
    addps  %xmm9,%xmm5
    addps  %xmm5,%xmm1
    movaps %xmm1,nb313_vctot(%rsp)

    movaps nb313_tsc(%rsp),%xmm10
    mulps  %xmm10,%xmm3 ## fscal
    mulps  %xmm10,%xmm7
    mulps  %xmm11,%xmm10

    movd %mm0,%eax
    movd %mm1,%ebx
    movd %mm2,%ecx
    movd %mm3,%edx

        ## move j forces to local temp variables 
    movq nb313_faction(%rbp),%rdi
    movlps (%rdi,%rax,4),%xmm11 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm12 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm11 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm12 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm13    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm14    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm2     ## jzb
    movss  8(%rdi,%rdx,4),%xmm5     ## jzd
    movlhps %xmm2,%xmm13 ## jza  -  jzb  -
    movlhps %xmm5,%xmm14 ## jzc  -  jzd -

    shufps $136,%xmm14,%xmm13 ## 10001000 => jza jzb jzc jzd

    ## xmm11: jxa jya jxb jyb 
    ## xmm12: jxc jyc jxd jyd
    ## xmm13: jza jzb jzc jzd

    xorps  %xmm0,%xmm0
    xorps  %xmm4,%xmm4
    xorps  %xmm8,%xmm8

    mulps  nb313_rinvH1(%rsp),%xmm3
    mulps  nb313_rinvH2(%rsp),%xmm7
    mulps  nb313_rinvM(%rsp),%xmm10

    subps  %xmm3,%xmm0
    subps  %xmm7,%xmm4
    subps  %xmm10,%xmm8

    movaps %xmm0,%xmm1
    movaps %xmm0,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm8,%xmm6
    movaps %xmm8,%xmm7

        mulps nb313_dxH1(%rsp),%xmm0
        mulps nb313_dyH1(%rsp),%xmm1
        mulps nb313_dzH1(%rsp),%xmm2
        mulps nb313_dxH2(%rsp),%xmm3
        mulps nb313_dyH2(%rsp),%xmm4
        mulps nb313_dzH2(%rsp),%xmm5
        mulps nb313_dxM(%rsp),%xmm6
        mulps nb313_dyM(%rsp),%xmm7
        mulps nb313_dzM(%rsp),%xmm8

    ## fetch forces from O interaction
    movaps nb313_fjx(%rsp),%xmm14
    movaps nb313_fjy(%rsp),%xmm15
    addps  nb313_fjz(%rsp),%xmm13

    addps %xmm0,%xmm14
    addps %xmm1,%xmm15
    addps %xmm2,%xmm13
    addps nb313_fixH1(%rsp),%xmm0
    addps nb313_fiyH1(%rsp),%xmm1
    addps nb313_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm14
    addps %xmm4,%xmm15
    addps %xmm5,%xmm13
    addps nb313_fixH2(%rsp),%xmm3
    addps nb313_fiyH2(%rsp),%xmm4
    addps nb313_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm14
    addps %xmm7,%xmm15
    addps %xmm8,%xmm13
    addps nb313_fixM(%rsp),%xmm6
    addps nb313_fiyM(%rsp),%xmm7
    addps nb313_fizM(%rsp),%xmm8

    movaps %xmm0,nb313_fixH1(%rsp)
    movaps %xmm1,nb313_fiyH1(%rsp)
    movaps %xmm2,nb313_fizH1(%rsp)
    movaps %xmm3,nb313_fixH2(%rsp)
    movaps %xmm4,nb313_fiyH2(%rsp)
    movaps %xmm5,nb313_fizH2(%rsp)
    movaps %xmm6,nb313_fixM(%rsp)
    movaps %xmm7,nb313_fiyM(%rsp)
    movaps %xmm8,nb313_fizM(%rsp)

    ## xmm14 = fjx
    ## xmm15 = fjy
    ## xmm13 = fjz
    movaps %xmm14,%xmm0
    unpcklps %xmm15,%xmm14
    unpckhps %xmm15,%xmm0

    addps  %xmm14,%xmm11
    addps  %xmm0,%xmm12

    movhlps  %xmm13,%xmm14 ## fjzc fjzd

    movlps %xmm11,(%rdi,%rax,4)
    movhps %xmm11,(%rdi,%rbx,4)
    movlps %xmm12,(%rdi,%rcx,4)
    movhps %xmm12,(%rdi,%rdx,4)
    movss  %xmm13,8(%rdi,%rax,4)
    movss  %xmm14,8(%rdi,%rcx,4)
    shufps $1,%xmm13,%xmm13
    shufps $1,%xmm14,%xmm14
    movss  %xmm13,8(%rdi,%rbx,4)
    movss  %xmm14,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb313_innerk(%rsp)
        jl    _nb_kernel313_x86_64_sse.nb313_odd_inner
        jmp   _nb_kernel313_x86_64_sse.nb313_unroll_loop
_nb_kernel313_x86_64_sse.nb313_odd_inner: 
        addl $4,nb313_innerk(%rsp)
        jnz   _nb_kernel313_x86_64_sse.nb313_odd_loop
        jmp   _nb_kernel313_x86_64_sse.nb313_updateouterdata
_nb_kernel313_x86_64_sse.nb313_odd_loop: 
        movq  nb313_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb313_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb313_iqM(%rsp),%xmm4
        movq nb313_charge(%rbp),%rsi
        movhps nb313_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb313_qqM(%rsp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb313_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb313_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb313_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb313_c6(%rsp)
        movaps %xmm7,nb313_c12(%rsp)

        movq nb313_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb313_ixO(%rsp),%xmm0
        movss nb313_iyO(%rsp),%xmm1
        movss nb313_izO(%rsp),%xmm2
        movss nb313_ixH1(%rsp),%xmm3
        movss nb313_iyH1(%rsp),%xmm4
        movss nb313_izH1(%rsp),%xmm5
        unpcklps nb313_ixH2(%rsp),%xmm0         ## ixO ixH2 - -
        unpcklps nb313_iyH2(%rsp),%xmm1         ## iyO iyH2 - -
        unpcklps nb313_izH2(%rsp),%xmm2         ## izO izH2 - -
        unpcklps nb313_ixM(%rsp),%xmm3          ## ixH1 ixM - -
        unpcklps nb313_iyM(%rsp),%xmm4          ## iyH1 iyM - -
        unpcklps nb313_izM(%rsp),%xmm5          ## izH1 izM - -
        unpcklps %xmm3,%xmm0    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm4,%xmm1    ## same for y
        unpcklps %xmm5,%xmm2    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use M distances for storage
        movaps %xmm3,nb313_dxM(%rsp)
        movaps %xmm4,nb313_dyM(%rsp)
        movaps %xmm5,nb313_dzM(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb313_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb313_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb313_rinvM(%rsp)
        mulps  %xmm0,%xmm4      ## r

        mulps nb313_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm7
        movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movq nb313_VFtab(%rbp),%rsi
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        xorps  %xmm5,%xmm5
        movlps (%rsi,%rcx,4),%xmm3      ## data: Y3 F3  -  - 
        movhps (%rsi,%rbx,4),%xmm5      ## data:  0  0 Y2 F2
        movhps (%rsi,%rdx,4),%xmm3      ## data: Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## data:  0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## data:  0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## data:  0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3     ## data: G3 H3  -  - 
        movhps 8(%rsi,%rbx,4),%xmm7     ## data:  0  0 G2 H2
        movhps 8(%rsi,%rdx,4),%xmm3     ## data: G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## data:  0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## data:  0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## data:  0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4

        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        mulps  nb313_two(%rsp),%xmm7            ## two*Heps2 
        movaps nb313_qqM(%rsp),%xmm0
        addps  %xmm6,%xmm7
        addps  %xmm5,%xmm7 ## xmm7=FF 
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        mulps  %xmm7,%xmm0 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm0 fijC 
        ## increment vcoul - then we can get rid of mm5 
        addps  nb313_vctot(%rsp),%xmm5
        movaps %xmm5,nb313_vctot(%rsp)

        ## do nontable L-J  in first element only.
        movaps nb313_rinvM(%rsp),%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb313_c6(%rsp),%xmm1
        mulss  nb313_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb313_six(%rsp),%xmm1
        mulss  nb313_twelve(%rsp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb313_Vvdwtot(%rsp),%xmm3
        mulss  nb313_rinvM(%rsp),%xmm4
        ## add back coul stuff from memory, and work on all elements again
        mulps  nb313_tsc(%rsp),%xmm0
        subps  %xmm0,%xmm4
        movss %xmm3,nb313_Vvdwtot(%rsp)
        mulps  nb313_rinvM(%rsp),%xmm4

        movaps nb313_dxM(%rsp),%xmm0
        movaps nb313_dyM(%rsp),%xmm1
        movaps nb313_dzM(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb313_fixO(%rsp),%xmm3
        movss  nb313_fiyO(%rsp),%xmm4
        movss  nb313_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb313_fixO(%rsp)
        movss  %xmm4,nb313_fiyO(%rsp)
        movss  %xmm5,nb313_fizO(%rsp)   ## updated the O force now do the H's


        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixH1(%rsp),%xmm3
        addss  nb313_fiyH1(%rsp),%xmm4
        addss  nb313_fizH1(%rsp),%xmm5
        movss  %xmm3,nb313_fixH1(%rsp)
        movss  %xmm4,nb313_fiyH1(%rsp)
        movss  %xmm5,nb313_fizH1(%rsp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixH2(%rsp),%xmm3
        addss  nb313_fiyH2(%rsp),%xmm4
        addss  nb313_fizH2(%rsp),%xmm5
        movss  %xmm3,nb313_fixH2(%rsp)
        movss  %xmm4,nb313_fiyH2(%rsp)
        movss  %xmm5,nb313_fizH2(%rsp)          ## updated the H2 force 

        movq nb313_faction(%rbp),%rdi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb313_fixM(%rsp),%xmm3
        addss  nb313_fiyM(%rsp),%xmm4
        addss  nb313_fizM(%rsp),%xmm5
        movss  %xmm3,nb313_fixM(%rsp)
        movss  %xmm4,nb313_fiyM(%rsp)
        movss  %xmm5,nb313_fizM(%rsp)   ## updated the M force 

        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,%xmm0
        movaps  %xmm4,%xmm1
        movaps  %xmm5,%xmm2

        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2
        unpcklps %xmm1,%xmm0    ## x,y sum in xmm0, z sum in xmm2

        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb313_innerk(%rsp)
        jz    _nb_kernel313_x86_64_sse.nb313_updateouterdata
        jmp   _nb_kernel313_x86_64_sse.nb313_odd_loop
_nb_kernel313_x86_64_sse.nb313_updateouterdata: 
        movl  nb313_ii3(%rsp),%ecx
        movq  nb313_faction(%rbp),%rdi
        movq  nb313_fshift(%rbp),%rsi
        movl  nb313_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb313_fixO(%rsp),%xmm0
        movaps nb313_fiyO(%rsp),%xmm1
        movaps nb313_fizO(%rsp),%xmm2

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
        movaps nb313_fixH1(%rsp),%xmm0
        movaps nb313_fiyH1(%rsp),%xmm1
        movaps nb313_fizH1(%rsp),%xmm2

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
        movaps nb313_fixH2(%rsp),%xmm0
        movaps nb313_fiyH2(%rsp),%xmm1
        movaps nb313_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb313_fixM(%rsp),%xmm0
        movaps nb313_fiyM(%rsp),%xmm1
        movaps nb313_fizM(%rsp),%xmm2

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
        movss  36(%rdi,%rcx,4),%xmm3
        movss  40(%rdi,%rcx,4),%xmm4
        movss  44(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,36(%rdi,%rcx,4)
        movss  %xmm4,40(%rdi,%rcx,4)
        movss  %xmm5,44(%rdi,%rcx,4)

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
        movl nb313_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb313_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb313_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb313_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb313_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb313_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb313_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313_x86_64_sse.nb313_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313_n(%rsp)
        jmp _nb_kernel313_x86_64_sse.nb313_outer
_nb_kernel313_x86_64_sse.nb313_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313_x86_64_sse.nb313_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313_x86_64_sse.nb313_threadloop
_nb_kernel313_x86_64_sse.nb313_end: 
        movl nb313_nouter(%rsp),%eax
        movl nb313_ninner(%rsp),%ebx
        movq nb313_outeriter(%rbp),%rcx
        movq nb313_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1048,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel313nf_x86_64_sse
.globl _nb_kernel313nf_x86_64_sse
nb_kernel313nf_x86_64_sse:      
_nb_kernel313nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb313nf_fshift, 16
.set nb313nf_gid, 24
.set nb313nf_pos, 32
.set nb313nf_faction, 40
.set nb313nf_charge, 48
.set nb313nf_p_facel, 56
.set nb313nf_argkrf, 64
.set nb313nf_argcrf, 72
.set nb313nf_Vc, 80
.set nb313nf_type, 88
.set nb313nf_p_ntype, 96
.set nb313nf_vdwparam, 104
.set nb313nf_Vvdw, 112
.set nb313nf_p_tabscale, 120
.set nb313nf_VFtab, 128
.set nb313nf_invsqrta, 136
.set nb313nf_dvda, 144
.set nb313nf_p_gbtabscale, 152
.set nb313nf_GBtab, 160
.set nb313nf_p_nthreads, 168
.set nb313nf_count, 176
.set nb313nf_mtx, 184
.set nb313nf_outeriter, 192
.set nb313nf_inneriter, 200
.set nb313nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb313nf_ixO, 0
.set nb313nf_iyO, 16
.set nb313nf_izO, 32
.set nb313nf_ixH1, 48
.set nb313nf_iyH1, 64
.set nb313nf_izH1, 80
.set nb313nf_ixH2, 96
.set nb313nf_iyH2, 112
.set nb313nf_izH2, 128
.set nb313nf_ixM, 144
.set nb313nf_iyM, 160
.set nb313nf_izM, 176
.set nb313nf_iqM, 192
.set nb313nf_iqH, 208
.set nb313nf_qqM, 224
.set nb313nf_qqH, 240
.set nb313nf_rinvH1, 256
.set nb313nf_rinvH2, 272
.set nb313nf_rinvM, 288
.set nb313nf_rH1, 304
.set nb313nf_rH2, 320
.set nb313nf_rM, 336
.set nb313nf_tsc, 352
.set nb313nf_two, 368
.set nb313nf_c6, 384
.set nb313nf_c12, 400
.set nb313nf_vctot, 416
.set nb313nf_Vvdwtot, 432
.set nb313nf_half, 448
.set nb313nf_three, 464
.set nb313nf_is3, 480
.set nb313nf_ii3, 484
.set nb313nf_nri, 488
.set nb313nf_iinr, 496
.set nb313nf_jindex, 504
.set nb313nf_jjnr, 512
.set nb313nf_shift, 520
.set nb313nf_shiftvec, 528
.set nb313nf_facel, 536
.set nb313nf_innerjjnr, 544
.set nb313nf_ntia, 552
.set nb313nf_innerk, 556
.set nb313nf_n, 560
.set nb313nf_nn1, 564
.set nb313nf_nouter, 568
.set nb313nf_ninner, 572
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $584,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb313nf_nouter(%rsp)
        movl %eax,nb313nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb313nf_nri(%rsp)
        movq %rsi,nb313nf_iinr(%rsp)
        movq %rdx,nb313nf_jindex(%rsp)
        movq %rcx,nb313nf_jjnr(%rsp)
        movq %r8,nb313nf_shift(%rsp)
        movq %r9,nb313nf_shiftvec(%rsp)
        movq nb313nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb313nf_facel(%rsp)

        movq nb313nf_p_tabscale(%rbp),%rax
        movss (%rax),%xmm3
        shufps $0,%xmm3,%xmm3
        movaps %xmm3,nb313nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb313nf_half(%rsp)
        movss nb313nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb313nf_half(%rsp)
        movaps %xmm2,nb313nf_two(%rsp)
        movaps %xmm3,nb313nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb313nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb313nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb313nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb313nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb313nf_iqM(%rsp)
        movaps %xmm4,nb313nf_iqH(%rsp)

        movq  nb313nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb313nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb313nf_ntia(%rsp)
_nb_kernel313nf_x86_64_sse.nb313nf_threadloop: 
        movq  nb313nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel313nf_x86_64_sse.nb313nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel313nf_x86_64_sse.nb313nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb313nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb313nf_n(%rsp)
        movl %ebx,nb313nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel313nf_x86_64_sse.nb313nf_outerstart
        jmp _nb_kernel313nf_x86_64_sse.nb313nf_end

_nb_kernel313nf_x86_64_sse.nb313nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb313nf_nouter(%rsp),%ebx
        movl %ebx,nb313nf_nouter(%rsp)

_nb_kernel313nf_x86_64_sse.nb313nf_outer: 
        movq  nb313nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb313nf_is3(%rsp)            ## store is3 

        movq  nb313nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb313nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb313nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb313nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb313nf_ixO(%rsp)
        movaps %xmm4,nb313nf_iyO(%rsp)
        movaps %xmm5,nb313nf_izO(%rsp)
        movaps %xmm6,nb313nf_ixH1(%rsp)
        movaps %xmm7,nb313nf_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb313nf_izH1(%rsp)
        movaps %xmm0,nb313nf_ixH2(%rsp)
        movaps %xmm1,nb313nf_iyH2(%rsp)
        movaps %xmm2,nb313nf_izH2(%rsp)
        movaps %xmm3,nb313nf_ixM(%rsp)
        movaps %xmm4,nb313nf_iyM(%rsp)
        movaps %xmm5,nb313nf_izM(%rsp)

        ## clear vctot
        xorps %xmm4,%xmm4
        movaps %xmm4,nb313nf_vctot(%rsp)
        movaps %xmm4,nb313nf_Vvdwtot(%rsp)

        movq  nb313nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb313nf_pos(%rbp),%rsi
        movq  nb313nf_faction(%rbp),%rdi
        movq  nb313nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb313nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb313nf_ninner(%rsp),%ecx
        movl  %ecx,nb313nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb313nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel313nf_x86_64_sse.nb313nf_unroll_loop
        jmp   _nb_kernel313nf_x86_64_sse.nb313nf_odd_inner
_nb_kernel313nf_x86_64_sse.nb313nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb313nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb313nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb313nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb313nf_iqM(%rsp),%xmm3
        mulps  nb313nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb313nf_qqM(%rsp)
        movaps  %xmm4,nb313nf_qqH(%rsp)

        movq nb313nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb313nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb313nf_ntia(%rsp),%edi
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

        movaps %xmm4,nb313nf_c6(%rsp)
        movaps %xmm6,nb313nf_c12(%rsp)

        movq nb313nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
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
        movaps nb313nf_ixO(%rsp),%xmm4
        movaps nb313nf_iyO(%rsp),%xmm5
        movaps nb313nf_izO(%rsp),%xmm6

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
        movaps nb313nf_ixH1(%rsp),%xmm4
        movaps nb313nf_iyH1(%rsp),%xmm5
        movaps nb313nf_izH1(%rsp),%xmm6

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
        movaps nb313nf_ixH2(%rsp),%xmm3
        movaps nb313nf_iyH2(%rsp),%xmm4
        movaps nb313nf_izH2(%rsp),%xmm5

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

        ## move ixM-izM to xmm2-xmm4  
        movaps nb313nf_iyM(%rsp),%xmm3
        movaps nb313nf_izM(%rsp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb313nf_ixM(%rsp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313nf_three(%rsp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%rsp),%xmm0
        movaps  %xmm0,nb313nf_rinvH1(%rsp)      ## rinvH1 in xmm4 
        mulps   %xmm0,%xmm6
        movaps  %xmm6,nb313nf_rH1(%rsp)

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313nf_three(%rsp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%rsp),%xmm0
        movaps  %xmm0,nb313nf_rinvH2(%rsp)      ## rinvH2 in xmm4 
        mulps   %xmm0,%xmm5
        movaps  %xmm5,nb313nf_rH2(%rsp)

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb313nf_three(%rsp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb313nf_half(%rsp),%xmm0
        movaps  %xmm0,nb313nf_rinvM(%rsp)       ## rinvM in xmm5 
        mulps   %xmm0,%xmm4
        movaps  %xmm4,nb313nf_rM(%rsp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb313nf_two(%rsp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb313nf_c6(%rsp),%xmm0
        mulps   nb313nf_c12(%rsp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb313nf_Vvdwtot(%rsp),%xmm3
        movaps  %xmm3,nb313nf_Vvdwtot(%rsp)

        ## Do H1 interaction
        movq nb313nf_VFtab(%rbp),%rsi

        movaps nb313nf_rH1(%rsp),%xmm7
        mulps   nb313nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %eax,%mm0
        movd %ebx,%mm1
        movd %ecx,%mm2
        movd %edx,%mm3

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqH(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb313nf_vctot(%rsp)

        ## Done with H1, do H2 interactions 
        movaps nb313nf_rH2(%rsp),%xmm7
        mulps   nb313nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## shuffle 10001000
        shufps $221,%xmm7,%xmm5 ## shuffle 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## shuf 10001000
        shufps $221,%xmm3,%xmm7 ## shuf 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqH(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb313nf_vctot(%rsp)

        ## Done with H2, do M interactions 
        movaps nb313nf_rM(%rsp),%xmm7
        mulps   nb313nf_tsc(%rsp),%xmm7
        movhlps %xmm7,%xmm4
        cvttps2pi %xmm7,%mm6
        cvttps2pi %xmm4,%mm7    ## mm6/mm7 contain lu indices 

        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm4
        movlhps %xmm4,%xmm3

        subps %xmm3,%xmm7
        movaps %xmm7,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm6,%ebx
        movd %mm7,%edx

        movlps (%rsi,%rax,4),%xmm5
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm5
        movhps (%rsi,%rdx,4),%xmm7 ## got half coulomb table 

        movaps %xmm5,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm5 ## 11011101

        movlps 8(%rsi,%rax,4),%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3
        movhps 8(%rsi,%rbx,4),%xmm7
        movhps 8(%rsi,%rdx,4),%xmm3    ## other half of coulomb table  
        movaps %xmm7,%xmm6
        shufps $136,%xmm3,%xmm6 ## 10001000
        shufps $221,%xmm3,%xmm7 ## 11011101
        ## coulomb table ready, in xmm4-xmm7      

        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqM(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb313nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb313nf_innerk(%rsp)
        jl    _nb_kernel313nf_x86_64_sse.nb313nf_odd_inner
        jmp   _nb_kernel313nf_x86_64_sse.nb313nf_unroll_loop
_nb_kernel313nf_x86_64_sse.nb313nf_odd_inner: 
        addl $4,nb313nf_innerk(%rsp)
        jnz   _nb_kernel313nf_x86_64_sse.nb313nf_odd_loop
        jmp   _nb_kernel313nf_x86_64_sse.nb313nf_updateouterdata
_nb_kernel313nf_x86_64_sse.nb313nf_odd_loop: 
        movq  nb313nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb313nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb313nf_iqM(%rsp),%xmm4
        movq nb313nf_charge(%rbp),%rsi
        movhps nb313nf_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb313nf_qqM(%rsp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb313nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb313nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb313nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb313nf_c6(%rsp)
        movaps %xmm7,nb313nf_c12(%rsp)

        movq nb313nf_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb313nf_ixO(%rsp),%xmm3
        movss nb313nf_iyO(%rsp),%xmm4
        movss nb313nf_izO(%rsp),%xmm5
        movss nb313nf_ixH1(%rsp),%xmm0
        movss nb313nf_iyH1(%rsp),%xmm1
        movss nb313nf_izH1(%rsp),%xmm2
        unpcklps nb313nf_ixH2(%rsp),%xmm3       ## ixO ixH2 - -
        unpcklps nb313nf_iyH2(%rsp),%xmm4       ## iyO iyH2 - -
        unpcklps nb313nf_izH2(%rsp),%xmm5       ## izO izH2 - -
        unpcklps nb313nf_ixM(%rsp),%xmm0        ## ixH1 ixM - -
        unpcklps nb313nf_iyM(%rsp),%xmm1        ## iyH1 iyM - -
        unpcklps nb313nf_izM(%rsp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb313nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb313nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,nb313nf_rinvM(%rsp)
        mulps  %xmm0,%xmm4      ## r

        mulps nb313nf_tsc(%rsp),%xmm4
        movhlps %xmm4,%xmm7
        cvttps2pi %xmm4,%mm6
        cvttps2pi %xmm7,%mm7    ## mm6/mm7 contain lu indices 
        cvtpi2ps %mm6,%xmm3
        cvtpi2ps %mm7,%xmm7
        movlhps %xmm7,%xmm3

        subps   %xmm3,%xmm4
        movaps %xmm4,%xmm1      ## xmm1=eps 
        movaps %xmm1,%xmm2
        mulps  %xmm2,%xmm2      ## xmm2=eps2 
        pslld $2,%mm6
        pslld $2,%mm7

        movq nb313nf_VFtab(%rbp),%rsi
        psrlq $32,%mm6
        movd %mm6,%ebx
        movd %mm7,%ecx
        psrlq $32,%mm7
        movd %mm7,%edx

        xorps  %xmm5,%xmm5
        movlps (%rsi,%rcx,4),%xmm3      ## data: Y3 F3  -  - 
        movhps (%rsi,%rbx,4),%xmm5      ## data:  0  0 Y2 F2
        movhps (%rsi,%rdx,4),%xmm3      ## data: Y3 F3 Y4 F4 

        movaps %xmm5,%xmm4              ## data:  0  0 Y2 F2 
        shufps $0x88,%xmm3,%xmm4       ## data:  0 Y2 Y3 Y3
        shufps $0xDD,%xmm3,%xmm5       ## data:  0 F2 F3 F4 

        xorps  %xmm7,%xmm7
        movlps 8(%rsi,%rcx,4),%xmm3     ## data: G3 H3  -  - 
        movhps 8(%rsi,%rbx,4),%xmm7     ## data:  0  0 G2 H2
        movhps 8(%rsi,%rdx,4),%xmm3     ## data: G3 H3 G4 H4 

        movaps %xmm7,%xmm6              ## data:  0  0 G2 H2 
        shufps $0x88,%xmm3,%xmm6       ## data:  0 G2 G3 G3
        shufps $0xDD,%xmm3,%xmm7       ## data:  0 H2 H3 H4 

        ## xmm4 =  0  Y2 Y3 Y4
        ## xmm5 =  0  F2 F3 F4
        ## xmm6 =  0  G2 G3 G4
        ## xmm7 =  0  H2 H3 H4

        ## coulomb table ready, in xmm4-xmm7      
        mulps  %xmm1,%xmm6      ## xmm6=Geps 
        mulps  %xmm2,%xmm7      ## xmm7=Heps2 
        addps  %xmm6,%xmm5
        addps  %xmm7,%xmm5      ## xmm5=Fp        
        movaps nb313nf_qqM(%rsp),%xmm0
        mulps  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addps  %xmm4,%xmm5 ## xmm5=VV 
        mulps  %xmm0,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 
        ## increment vcoul 
        addps  nb313nf_vctot(%rsp),%xmm5
        movaps %xmm5,nb313nf_vctot(%rsp)

        ## do nontable L-J  in first element only.
        movaps nb313nf_rinvM(%rsp),%xmm2
        mulss  %xmm2,%xmm2
        movaps %xmm2,%xmm1
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb313nf_c6(%rsp),%xmm1
        mulss  nb313nf_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6
        addss  nb313nf_Vvdwtot(%rsp),%xmm3
        movss %xmm3,nb313nf_Vvdwtot(%rsp)

        decl nb313nf_innerk(%rsp)
        jz    _nb_kernel313nf_x86_64_sse.nb313nf_updateouterdata
        jmp   _nb_kernel313nf_x86_64_sse.nb313nf_odd_loop
_nb_kernel313nf_x86_64_sse.nb313nf_updateouterdata: 
        ## get n from stack
        movl nb313nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb313nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb313nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb313nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb313nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb313nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb313nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel313nf_x86_64_sse.nb313nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb313nf_n(%rsp)
        jmp _nb_kernel313nf_x86_64_sse.nb313nf_outer
_nb_kernel313nf_x86_64_sse.nb313nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb313nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel313nf_x86_64_sse.nb313nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel313nf_x86_64_sse.nb313nf_threadloop
_nb_kernel313nf_x86_64_sse.nb313nf_end: 
        movl nb313nf_nouter(%rsp),%eax
        movl nb313nf_ninner(%rsp),%ebx
        movq nb313nf_outeriter(%rbp),%rcx
        movq nb313nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $584,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





