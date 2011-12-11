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






.globl nb_kernel302_x86_64_sse2
.globl _nb_kernel302_x86_64_sse2
nb_kernel302_x86_64_sse2:       
_nb_kernel302_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb302_fshift, 16
.set nb302_gid, 24
.set nb302_pos, 32
.set nb302_faction, 40
.set nb302_charge, 48
.set nb302_p_facel, 56
.set nb302_argkrf, 64
.set nb302_argcrf, 72
.set nb302_Vc, 80
.set nb302_type, 88
.set nb302_p_ntype, 96
.set nb302_vdwparam, 104
.set nb302_Vvdw, 112
.set nb302_p_tabscale, 120
.set nb302_VFtab, 128
.set nb302_invsqrta, 136
.set nb302_dvda, 144
.set nb302_p_gbtabscale, 152
.set nb302_GBtab, 160
.set nb302_p_nthreads, 168
.set nb302_count, 176
.set nb302_mtx, 184
.set nb302_outeriter, 192
.set nb302_inneriter, 200
.set nb302_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb302_ixO, 0
.set nb302_iyO, 16
.set nb302_izO, 32
.set nb302_ixH1, 48
.set nb302_iyH1, 64
.set nb302_izH1, 80
.set nb302_ixH2, 96
.set nb302_iyH2, 112
.set nb302_izH2, 128
.set nb302_jxO, 144
.set nb302_jyO, 160
.set nb302_jzO, 176
.set nb302_jxH1, 192
.set nb302_jyH1, 208
.set nb302_jzH1, 224
.set nb302_jxH2, 240
.set nb302_jyH2, 256
.set nb302_jzH2, 272
.set nb302_dxOO, 288
.set nb302_dyOO, 304
.set nb302_dzOO, 320
.set nb302_dxOH1, 336
.set nb302_dyOH1, 352
.set nb302_dzOH1, 368
.set nb302_dxOH2, 384
.set nb302_dyOH2, 400
.set nb302_dzOH2, 416
.set nb302_dxH1O, 432
.set nb302_dyH1O, 448
.set nb302_dzH1O, 464
.set nb302_dxH1H1, 480
.set nb302_dyH1H1, 496
.set nb302_dzH1H1, 512
.set nb302_dxH1H2, 528
.set nb302_dyH1H2, 544
.set nb302_dzH1H2, 560
.set nb302_dxH2O, 576
.set nb302_dyH2O, 592
.set nb302_dzH2O, 608
.set nb302_dxH2H1, 624
.set nb302_dyH2H1, 640
.set nb302_dzH2H1, 656
.set nb302_dxH2H2, 672
.set nb302_dyH2H2, 688
.set nb302_dzH2H2, 704
.set nb302_qqOO, 720
.set nb302_qqOH, 736
.set nb302_qqHH, 752
.set nb302_two, 768
.set nb302_tsc, 784
.set nb302_vctot, 800
.set nb302_fixO, 816
.set nb302_fiyO, 832
.set nb302_fizO, 848
.set nb302_fixH1, 864
.set nb302_fiyH1, 880
.set nb302_fizH1, 896
.set nb302_fixH2, 912
.set nb302_fiyH2, 928
.set nb302_fizH2, 944
.set nb302_epsO, 960
.set nb302_epsH1, 976
.set nb302_epsH2, 992
.set nb302_fjxH1, 1008
.set nb302_fjyH1, 1024
.set nb302_fjzH1, 1040
.set nb302_fjxH2, 1056
.set nb302_fjyH2, 1072
.set nb302_fjzH2, 1088
.set nb302_half, 1104
.set nb302_three, 1120
.set nb302_rsqOO, 1136
.set nb302_rsqOH1, 1152
.set nb302_rsqOH2, 1168
.set nb302_rsqH1O, 1184
.set nb302_rsqH1H1, 1200
.set nb302_rsqH1H2, 1216
.set nb302_rsqH2O, 1232
.set nb302_rsqH2H1, 1248
.set nb302_rsqH2H2, 1264
.set nb302_rinvOO, 1280
.set nb302_rinvOH1, 1296
.set nb302_rinvOH2, 1312
.set nb302_rinvH1O, 1328
.set nb302_rinvH1H1, 1344
.set nb302_rinvH1H2, 1360
.set nb302_rinvH2O, 1376
.set nb302_rinvH2H1, 1392
.set nb302_rinvH2H2, 1408
.set nb302_is3, 1424
.set nb302_ii3, 1428
.set nb302_nri, 1432
.set nb302_iinr, 1440
.set nb302_jindex, 1448
.set nb302_jjnr, 1456
.set nb302_shift, 1464
.set nb302_shiftvec, 1472
.set nb302_facel, 1480
.set nb302_innerjjnr, 1488
.set nb302_innerk, 1496
.set nb302_n, 1500
.set nb302_nn1, 1504
.set nb302_nouter, 1508
.set nb302_ninner, 1512
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1528,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb302_nouter(%rsp)
        movl %eax,nb302_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb302_nri(%rsp)
        movq %rsi,nb302_iinr(%rsp)
        movq %rdx,nb302_jindex(%rsp)
        movq %rcx,nb302_jjnr(%rsp)
        movq %r8,nb302_shift(%rsp)
        movq %r9,nb302_shiftvec(%rsp)
        movq nb302_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb302_facel(%rsp)

        movq nb302_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb302_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb302_half(%rsp)
        movl %ebx,nb302_half+4(%rsp)
        movsd nb302_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb302_half(%rsp)
        movapd %xmm2,nb302_two(%rsp)
        movapd %xmm3,nb302_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb302_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb302_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb302_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb302_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302_qqOO(%rsp)
        movapd %xmm4,nb302_qqOH(%rsp)
        movapd %xmm5,nb302_qqHH(%rsp)

_nb_kernel302_x86_64_sse2.nb302_threadloop: 
        movq  nb302_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel302_x86_64_sse2.nb302_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302_x86_64_sse2.nb302_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302_n(%rsp)
        movl %ebx,nb302_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302_x86_64_sse2.nb302_outerstart
        jmp _nb_kernel302_x86_64_sse2.nb302_end

_nb_kernel302_x86_64_sse2.nb302_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302_nouter(%rsp),%ebx
        movl %ebx,nb302_nouter(%rsp)

_nb_kernel302_x86_64_sse2.nb302_outer: 
        movq  nb302_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb302_is3(%rsp)      ## store is3 

        movq  nb302_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb302_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb302_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb302_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302_ixO(%rsp)
        movapd %xmm4,nb302_iyO(%rsp)
        movapd %xmm5,nb302_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb302_ixH1(%rsp)
        movapd %xmm1,nb302_iyH1(%rsp)
        movapd %xmm2,nb302_izH1(%rsp)
        movapd %xmm3,nb302_ixH2(%rsp)
        movapd %xmm4,nb302_iyH2(%rsp)
        movapd %xmm5,nb302_izH2(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb302_vctot(%rsp)
        movapd %xmm4,nb302_fixO(%rsp)
        movapd %xmm4,nb302_fiyO(%rsp)
        movapd %xmm4,nb302_fizO(%rsp)
        movapd %xmm4,nb302_fixH1(%rsp)
        movapd %xmm4,nb302_fiyH1(%rsp)
        movapd %xmm4,nb302_fizH1(%rsp)
        movapd %xmm4,nb302_fixH2(%rsp)
        movapd %xmm4,nb302_fiyH2(%rsp)
        movapd %xmm4,nb302_fizH2(%rsp)

        movq  nb302_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb302_pos(%rbp),%rsi
        movq  nb302_faction(%rbp),%rdi
        movq  nb302_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb302_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb302_ninner(%rsp),%ecx
        movl  %ecx,nb302_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb302_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel302_x86_64_sse2.nb302_unroll_loop
        jmp   _nb_kernel302_x86_64_sse2.nb302_checksingle
_nb_kernel302_x86_64_sse2.nb302_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb302_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb302_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb302_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j O coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb302_ixO(%rsp),%xmm0
    subpd nb302_iyO(%rsp),%xmm1
    subpd nb302_izO(%rsp),%xmm2
    subpd nb302_ixH1(%rsp),%xmm3
    subpd nb302_iyH1(%rsp),%xmm4
    subpd nb302_izH1(%rsp),%xmm5
    subpd nb302_ixH2(%rsp),%xmm6
    subpd nb302_iyH2(%rsp),%xmm7
    subpd nb302_izH2(%rsp),%xmm8

        movapd %xmm0,nb302_dxOO(%rsp)
        movapd %xmm1,nb302_dyOO(%rsp)
        movapd %xmm2,nb302_dzOO(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1O(%rsp)
        movapd %xmm4,nb302_dyH1O(%rsp)
        movapd %xmm5,nb302_dzH1O(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb302_dxH2O(%rsp)
        movapd %xmm7,nb302_dyH2O(%rsp)
        movapd %xmm8,nb302_dzH2O(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb302_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb302_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOO 
        mulpd   %xmm15,%xmm10 ##   rinvH1O
    mulpd   %xmm15,%xmm11 ##   rinvH2O

        movapd  %xmm9,nb302_rinvOO(%rsp)
        movapd  %xmm10,nb302_rinvH1O(%rsp)
        movapd  %xmm11,nb302_rinvH2O(%rsp)

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb302_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb302_epsO(%rsp)
    movapd    %xmm3,nb302_epsH1(%rsp)
    movapd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb302_epsO(%rsp),%xmm12
    movapd nb302_epsH1(%rsp),%xmm13
    movapd nb302_epsH2(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb302_qqOO(%rsp),%xmm12
    movapd nb302_qqOH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm13,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm13,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb302_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb302_vctot(%rsp)

    movapd nb302_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb302_rinvOO(%rsp),%xmm4
    mulpd nb302_rinvH1O(%rsp),%xmm8
    mulpd nb302_rinvH2O(%rsp),%xmm11

    ## move j O forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb302_dxOO(%rsp),%xmm3
        mulpd nb302_dyOO(%rsp),%xmm4
        mulpd nb302_dzOO(%rsp),%xmm5
        mulpd nb302_dxH1O(%rsp),%xmm7
        mulpd nb302_dyH1O(%rsp),%xmm8
        mulpd nb302_dzH1O(%rsp),%xmm9
        mulpd nb302_dxH2O(%rsp),%xmm10
        mulpd nb302_dyH2O(%rsp),%xmm11
        mulpd nb302_dzH2O(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb302_fixO(%rsp),%xmm3
    addpd nb302_fiyO(%rsp),%xmm4
    addpd nb302_fizO(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb302_fixH1(%rsp),%xmm7
    addpd nb302_fiyH1(%rsp),%xmm8
    addpd nb302_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb302_fixH2(%rsp),%xmm10
    addpd nb302_fiyH2(%rsp),%xmm11
    addpd nb302_fizH2(%rsp),%xmm12

    movapd %xmm3,nb302_fixO(%rsp)
    movapd %xmm4,nb302_fiyO(%rsp)
    movapd %xmm5,nb302_fizO(%rsp)
    movapd %xmm7,nb302_fixH1(%rsp)
    movapd %xmm8,nb302_fiyH1(%rsp)
    movapd %xmm9,nb302_fizH1(%rsp)
    movapd %xmm10,nb302_fixH2(%rsp)
    movapd %xmm11,nb302_fiyH2(%rsp)
    movapd %xmm12,nb302_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## move j H1 coordinates to local temp variables 
    movq nb302_pos(%rbp),%rsi
    movlpd 24(%rsi,%rax,8),%xmm0
    movlpd 32(%rsi,%rax,8),%xmm1
    movlpd 40(%rsi,%rax,8),%xmm2
    movhpd 24(%rsi,%rbx,8),%xmm0
    movhpd 32(%rsi,%rbx,8),%xmm1
    movhpd 40(%rsi,%rbx,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb302_ixO(%rsp),%xmm0
    subpd nb302_iyO(%rsp),%xmm1
    subpd nb302_izO(%rsp),%xmm2
    subpd nb302_ixH1(%rsp),%xmm3
    subpd nb302_iyH1(%rsp),%xmm4
    subpd nb302_izH1(%rsp),%xmm5
    subpd nb302_ixH2(%rsp),%xmm6
    subpd nb302_iyH2(%rsp),%xmm7
    subpd nb302_izH2(%rsp),%xmm8

        movapd %xmm0,nb302_dxOH1(%rsp)
        movapd %xmm1,nb302_dyOH1(%rsp)
        movapd %xmm2,nb302_dzOH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1H1(%rsp)
        movapd %xmm4,nb302_dyH1H1(%rsp)
        movapd %xmm5,nb302_dzH1H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb302_dxH2H1(%rsp)
        movapd %xmm7,nb302_dyH2H1(%rsp)
        movapd %xmm8,nb302_dzH2H1(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb302_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb302_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH1
        mulpd   %xmm15,%xmm10 ##   rinvH1H1
    mulpd   %xmm15,%xmm11 ##   rinvH2H1

        movapd  %xmm9,nb302_rinvOH1(%rsp)
        movapd  %xmm10,nb302_rinvH1H1(%rsp)
        movapd  %xmm11,nb302_rinvH2H1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb302_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb302_epsO(%rsp)
    movapd    %xmm3,nb302_epsH1(%rsp)
    movapd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb302_epsO(%rsp),%xmm12
    movapd nb302_epsH1(%rsp),%xmm13
    movapd nb302_epsH2(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb302_qqOH(%rsp),%xmm12
    movapd nb302_qqHH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm13,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm13,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb302_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb302_vctot(%rsp)

    movapd nb302_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb302_rinvOH1(%rsp),%xmm4
    mulpd nb302_rinvH1H1(%rsp),%xmm8
    mulpd nb302_rinvH2H1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb302_dxOH1(%rsp),%xmm3
        mulpd nb302_dyOH1(%rsp),%xmm4
        mulpd nb302_dzOH1(%rsp),%xmm5
        mulpd nb302_dxH1H1(%rsp),%xmm7
        mulpd nb302_dyH1H1(%rsp),%xmm8
        mulpd nb302_dzH1H1(%rsp),%xmm9
        mulpd nb302_dxH2H1(%rsp),%xmm10
        mulpd nb302_dyH2H1(%rsp),%xmm11
        mulpd nb302_dzH2H1(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb302_fixO(%rsp),%xmm3
    addpd nb302_fiyO(%rsp),%xmm4
    addpd nb302_fizO(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb302_fixH1(%rsp),%xmm7
    addpd nb302_fiyH1(%rsp),%xmm8
    addpd nb302_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb302_fixH2(%rsp),%xmm10
    addpd nb302_fiyH2(%rsp),%xmm11
    addpd nb302_fizH2(%rsp),%xmm12

    movapd %xmm3,nb302_fixO(%rsp)
    movapd %xmm4,nb302_fiyO(%rsp)
    movapd %xmm5,nb302_fizO(%rsp)
    movapd %xmm7,nb302_fixH1(%rsp)
    movapd %xmm8,nb302_fiyH1(%rsp)
    movapd %xmm9,nb302_fizH1(%rsp)
    movapd %xmm10,nb302_fixH2(%rsp)
    movapd %xmm11,nb302_fiyH2(%rsp)
    movapd %xmm12,nb302_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb302_pos(%rbp),%rsi
    movlpd 48(%rsi,%rax,8),%xmm0
    movlpd 56(%rsi,%rax,8),%xmm1
    movlpd 64(%rsi,%rax,8),%xmm2
    movhpd 48(%rsi,%rbx,8),%xmm0
    movhpd 56(%rsi,%rbx,8),%xmm1
    movhpd 64(%rsi,%rbx,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb302_ixO(%rsp),%xmm0
    subpd nb302_iyO(%rsp),%xmm1
    subpd nb302_izO(%rsp),%xmm2
    subpd nb302_ixH1(%rsp),%xmm3
    subpd nb302_iyH1(%rsp),%xmm4
    subpd nb302_izH1(%rsp),%xmm5
    subpd nb302_ixH2(%rsp),%xmm6
    subpd nb302_iyH2(%rsp),%xmm7
    subpd nb302_izH2(%rsp),%xmm8

        movapd %xmm0,nb302_dxOH2(%rsp)
        movapd %xmm1,nb302_dyOH2(%rsp)
        movapd %xmm2,nb302_dzOH2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1H2(%rsp)
        movapd %xmm4,nb302_dyH1H2(%rsp)
        movapd %xmm5,nb302_dzH1H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb302_dxH2H2(%rsp)
        movapd %xmm7,nb302_dyH2H2(%rsp)
        movapd %xmm8,nb302_dzH2H2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb302_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb302_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb302_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvOH2
        mulpd   %xmm15,%xmm10 ##   rinvH1H2
    mulpd   %xmm15,%xmm11 ##   rinvH2H2


        movapd  %xmm9,nb302_rinvOH2(%rsp)
        movapd  %xmm10,nb302_rinvH1H2(%rsp)
        movapd  %xmm11,nb302_rinvH2H2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb302_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb302_epsO(%rsp)
    movapd    %xmm3,nb302_epsH1(%rsp)
    movapd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb302_epsO(%rsp),%xmm12
    movapd nb302_epsH1(%rsp),%xmm13
    movapd nb302_epsH2(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb302_qqOH(%rsp),%xmm12
    movapd nb302_qqHH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm13,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm13,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb302_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb302_vctot(%rsp)

    movapd nb302_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb302_rinvOH2(%rsp),%xmm4
    mulpd nb302_rinvH1H2(%rsp),%xmm8
    mulpd nb302_rinvH2H2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb302_dxOH2(%rsp),%xmm3
        mulpd nb302_dyOH2(%rsp),%xmm4
        mulpd nb302_dzOH2(%rsp),%xmm5
        mulpd nb302_dxH1H2(%rsp),%xmm7
        mulpd nb302_dyH1H2(%rsp),%xmm8
        mulpd nb302_dzH1H2(%rsp),%xmm9
        mulpd nb302_dxH2H2(%rsp),%xmm10
        mulpd nb302_dyH2H2(%rsp),%xmm11
        mulpd nb302_dzH2H2(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb302_fixO(%rsp),%xmm3
    addpd nb302_fiyO(%rsp),%xmm4
    addpd nb302_fizO(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb302_fixH1(%rsp),%xmm7
    addpd nb302_fiyH1(%rsp),%xmm8
    addpd nb302_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb302_fixH2(%rsp),%xmm10
    addpd nb302_fiyH2(%rsp),%xmm11
    addpd nb302_fizH2(%rsp),%xmm12

    movapd %xmm3,nb302_fixO(%rsp)
    movapd %xmm4,nb302_fiyO(%rsp)
    movapd %xmm5,nb302_fizO(%rsp)
    movapd %xmm7,nb302_fixH1(%rsp)
    movapd %xmm8,nb302_fiyH1(%rsp)
    movapd %xmm9,nb302_fizH1(%rsp)
    movapd %xmm10,nb302_fixH2(%rsp)
    movapd %xmm11,nb302_fiyH2(%rsp)
    movapd %xmm12,nb302_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb302_innerk(%rsp)
        jl    _nb_kernel302_x86_64_sse2.nb302_checksingle
        jmp   _nb_kernel302_x86_64_sse2.nb302_unroll_loop
_nb_kernel302_x86_64_sse2.nb302_checksingle: 
        movl  nb302_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel302_x86_64_sse2.nb302_dosingle
        jmp   _nb_kernel302_x86_64_sse2.nb302_updateouterdata
_nb_kernel302_x86_64_sse2.nb302_dosingle: 
        movq  nb302_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb302_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax


        ## move j O coordinates to local temp variables 
    movsd (%rsi,%rax,8),%xmm0
    movsd 8(%rsi,%rax,8),%xmm1
    movsd 16(%rsi,%rax,8),%xmm2

    ## xmm0 = Ox
    ## xmm1 = Oy
    ## xmm2 = Oz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb302_ixO(%rsp),%xmm0
    subsd nb302_iyO(%rsp),%xmm1
    subsd nb302_izO(%rsp),%xmm2
    subsd nb302_ixH1(%rsp),%xmm3
    subsd nb302_iyH1(%rsp),%xmm4
    subsd nb302_izH1(%rsp),%xmm5
    subsd nb302_ixH2(%rsp),%xmm6
    subsd nb302_iyH2(%rsp),%xmm7
    subsd nb302_izH2(%rsp),%xmm8

        movsd %xmm0,nb302_dxOO(%rsp)
        movsd %xmm1,nb302_dyOO(%rsp)
        movsd %xmm2,nb302_dzOO(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb302_dxH1O(%rsp)
        movsd %xmm4,nb302_dyH1O(%rsp)
        movsd %xmm5,nb302_dzH1O(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb302_dxH2O(%rsp)
        movsd %xmm7,nb302_dyH2O(%rsp)
        movsd %xmm8,nb302_dzH2O(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jO atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb302_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOO 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1O
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2O 

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb302_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOO 
        mulsd   %xmm15,%xmm10 ##   rinvH1O
    mulsd   %xmm15,%xmm11 ##   rinvH2O

        movsd  %xmm9,nb302_rinvOO(%rsp)
        movsd  %xmm10,nb302_rinvH1O(%rsp)
        movsd  %xmm11,nb302_rinvH2O(%rsp)

        ## O interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb302_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shll  $2,%r8d
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb302_epsO(%rsp)
    movsd    %xmm3,nb302_epsH1(%rsp)
    movsd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb302_epsO(%rsp),%xmm12
    movsd nb302_epsH1(%rsp),%xmm13
    movsd nb302_epsH2(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb302_qqOO(%rsp),%xmm12
    movsd nb302_qqOH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm13,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm13,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb302_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb302_vctot(%rsp)

    movsd nb302_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb302_rinvOO(%rsp),%xmm4
    mulsd nb302_rinvH1O(%rsp),%xmm8
    mulsd nb302_rinvH2O(%rsp),%xmm11

    ## move j O forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movsd (%rdi,%rax,8),%xmm0
        movsd 8(%rdi,%rax,8),%xmm1
        movsd 16(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb302_dxOO(%rsp),%xmm3
        mulsd nb302_dyOO(%rsp),%xmm4
        mulsd nb302_dzOO(%rsp),%xmm5
        mulsd nb302_dxH1O(%rsp),%xmm7
        mulsd nb302_dyH1O(%rsp),%xmm8
        mulsd nb302_dzH1O(%rsp),%xmm9
        mulsd nb302_dxH2O(%rsp),%xmm10
        mulsd nb302_dyH2O(%rsp),%xmm11
        mulsd nb302_dzH2O(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb302_fixO(%rsp),%xmm3
    addsd nb302_fiyO(%rsp),%xmm4
    addsd nb302_fizO(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb302_fixH1(%rsp),%xmm7
    addsd nb302_fiyH1(%rsp),%xmm8
    addsd nb302_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb302_fixH2(%rsp),%xmm10
    addsd nb302_fiyH2(%rsp),%xmm11
    addsd nb302_fizH2(%rsp),%xmm12

    movsd %xmm3,nb302_fixO(%rsp)
    movsd %xmm4,nb302_fiyO(%rsp)
    movsd %xmm5,nb302_fizO(%rsp)
    movsd %xmm7,nb302_fixH1(%rsp)
    movsd %xmm8,nb302_fiyH1(%rsp)
    movsd %xmm9,nb302_fizH1(%rsp)
    movsd %xmm10,nb302_fixH2(%rsp)
    movsd %xmm11,nb302_fiyH2(%rsp)
    movsd %xmm12,nb302_fizH2(%rsp)

    ## store back j O forces from xmm0-xmm2
        movsd %xmm0,(%rdi,%rax,8)
        movsd %xmm1,8(%rdi,%rax,8)
        movsd %xmm2,16(%rdi,%rax,8)

        ## move j H1 coordinates to local temp variables 
    movq nb302_pos(%rbp),%rsi
    movsd 24(%rsi,%rax,8),%xmm0
    movsd 32(%rsi,%rax,8),%xmm1
    movsd 40(%rsi,%rax,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb302_ixO(%rsp),%xmm0
    subsd nb302_iyO(%rsp),%xmm1
    subsd nb302_izO(%rsp),%xmm2
    subsd nb302_ixH1(%rsp),%xmm3
    subsd nb302_iyH1(%rsp),%xmm4
    subsd nb302_izH1(%rsp),%xmm5
    subsd nb302_ixH2(%rsp),%xmm6
    subsd nb302_iyH2(%rsp),%xmm7
    subsd nb302_izH2(%rsp),%xmm8

        movsd %xmm0,nb302_dxOH1(%rsp)
        movsd %xmm1,nb302_dyOH1(%rsp)
        movsd %xmm2,nb302_dzOH1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb302_dxH1H1(%rsp)
        movsd %xmm4,nb302_dyH1H1(%rsp)
        movsd %xmm5,nb302_dzH1H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb302_dxH2H1(%rsp)
        movsd %xmm7,nb302_dyH2H1(%rsp)
        movsd %xmm8,nb302_dzH2H1(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb302_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2OH1

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb302_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH1
        mulsd   %xmm15,%xmm10 ##   rinvH1H1
    mulsd   %xmm15,%xmm11 ##   rinvH2H1

        movsd  %xmm9,nb302_rinvOH1(%rsp)
        movsd  %xmm10,nb302_rinvH1H1(%rsp)
        movsd  %xmm11,nb302_rinvH2H1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb302_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shll  $2,%r8d
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb302_epsO(%rsp)
    movsd    %xmm3,nb302_epsH1(%rsp)
    movsd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb302_epsO(%rsp),%xmm12
    movsd nb302_epsH1(%rsp),%xmm13
    movsd nb302_epsH2(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb302_qqOH(%rsp),%xmm12
    movsd nb302_qqHH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm13,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm13,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb302_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb302_vctot(%rsp)

    movsd nb302_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb302_rinvOH1(%rsp),%xmm4
    mulsd nb302_rinvH1H1(%rsp),%xmm8
    mulsd nb302_rinvH2H1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb302_dxOH1(%rsp),%xmm3
        mulsd nb302_dyOH1(%rsp),%xmm4
        mulsd nb302_dzOH1(%rsp),%xmm5
        mulsd nb302_dxH1H1(%rsp),%xmm7
        mulsd nb302_dyH1H1(%rsp),%xmm8
        mulsd nb302_dzH1H1(%rsp),%xmm9
        mulsd nb302_dxH2H1(%rsp),%xmm10
        mulsd nb302_dyH2H1(%rsp),%xmm11
        mulsd nb302_dzH2H1(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb302_fixO(%rsp),%xmm3
    addsd nb302_fiyO(%rsp),%xmm4
    addsd nb302_fizO(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb302_fixH1(%rsp),%xmm7
    addsd nb302_fiyH1(%rsp),%xmm8
    addsd nb302_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb302_fixH2(%rsp),%xmm10
    addsd nb302_fiyH2(%rsp),%xmm11
    addsd nb302_fizH2(%rsp),%xmm12

    movsd %xmm3,nb302_fixO(%rsp)
    movsd %xmm4,nb302_fiyO(%rsp)
    movsd %xmm5,nb302_fizO(%rsp)
    movsd %xmm7,nb302_fixH1(%rsp)
    movsd %xmm8,nb302_fiyH1(%rsp)
    movsd %xmm9,nb302_fizH1(%rsp)
    movsd %xmm10,nb302_fixH2(%rsp)
    movsd %xmm11,nb302_fiyH2(%rsp)
    movsd %xmm12,nb302_fizH2(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movq nb302_pos(%rbp),%rsi
    movsd 48(%rsi,%rax,8),%xmm0
    movsd 56(%rsi,%rax,8),%xmm1
    movsd 64(%rsi,%rax,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb302_ixO(%rsp),%xmm0
    subsd nb302_iyO(%rsp),%xmm1
    subsd nb302_izO(%rsp),%xmm2
    subsd nb302_ixH1(%rsp),%xmm3
    subsd nb302_iyH1(%rsp),%xmm4
    subsd nb302_izH1(%rsp),%xmm5
    subsd nb302_ixH2(%rsp),%xmm6
    subsd nb302_iyH2(%rsp),%xmm7
    subsd nb302_izH2(%rsp),%xmm8

        movsd %xmm0,nb302_dxOH2(%rsp)
        movsd %xmm1,nb302_dyOH2(%rsp)
        movsd %xmm2,nb302_dzOH2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb302_dxH1H2(%rsp)
        movsd %xmm4,nb302_dyH1H2(%rsp)
        movsd %xmm5,nb302_dzH1H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb302_dxH2H2(%rsp)
        movsd %xmm7,nb302_dyH2H2(%rsp)
        movsd %xmm8,nb302_dzH2H2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb302_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvOH2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH1H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvH2H2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb302_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb302_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvOH2
        mulsd   %xmm15,%xmm10 ##   rinvH1H2
    mulsd   %xmm15,%xmm11 ##   rinvH2H2

        movsd  %xmm9,nb302_rinvOH2(%rsp)
        movsd  %xmm10,nb302_rinvH1H2(%rsp)
        movsd  %xmm11,nb302_rinvH2H2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb302_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shll  $2,%r8d
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb302_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb302_epsO(%rsp)
    movsd    %xmm3,nb302_epsH1(%rsp)
    movsd    %xmm6,nb302_epsH2(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb302_epsO(%rsp),%xmm12
    movsd nb302_epsH1(%rsp),%xmm13
    movsd nb302_epsH2(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb302_qqOH(%rsp),%xmm12
    movsd nb302_qqHH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm13,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm13,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb302_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb302_vctot(%rsp)

    movsd nb302_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb302_rinvOH2(%rsp),%xmm4
    mulsd nb302_rinvH1H2(%rsp),%xmm8
    mulsd nb302_rinvH2H2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb302_faction(%rbp),%rdi
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb302_dxOH2(%rsp),%xmm3
        mulsd nb302_dyOH2(%rsp),%xmm4
        mulsd nb302_dzOH2(%rsp),%xmm5
        mulsd nb302_dxH1H2(%rsp),%xmm7
        mulsd nb302_dyH1H2(%rsp),%xmm8
        mulsd nb302_dzH1H2(%rsp),%xmm9
        mulsd nb302_dxH2H2(%rsp),%xmm10
        mulsd nb302_dyH2H2(%rsp),%xmm11
        mulsd nb302_dzH2H2(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb302_fixO(%rsp),%xmm3
    addsd nb302_fiyO(%rsp),%xmm4
    addsd nb302_fizO(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb302_fixH1(%rsp),%xmm7
    addsd nb302_fiyH1(%rsp),%xmm8
    addsd nb302_fizH1(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb302_fixH2(%rsp),%xmm10
    addsd nb302_fiyH2(%rsp),%xmm11
    addsd nb302_fizH2(%rsp),%xmm12

    movsd %xmm3,nb302_fixO(%rsp)
    movsd %xmm4,nb302_fiyO(%rsp)
    movsd %xmm5,nb302_fizO(%rsp)
    movsd %xmm7,nb302_fixH1(%rsp)
    movsd %xmm8,nb302_fiyH1(%rsp)
    movsd %xmm9,nb302_fizH1(%rsp)
    movsd %xmm10,nb302_fixH2(%rsp)
    movsd %xmm11,nb302_fiyH2(%rsp)
    movsd %xmm12,nb302_fizH2(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

_nb_kernel302_x86_64_sse2.nb302_updateouterdata: 
        movl  nb302_ii3(%rsp),%ecx
        movq  nb302_faction(%rbp),%rdi
        movq  nb302_fshift(%rbp),%rsi
        movl  nb302_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb302_fixO(%rsp),%xmm0
        movapd nb302_fiyO(%rsp),%xmm1
        movapd nb302_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb302_fixH1(%rsp),%xmm0
        movapd nb302_fiyH1(%rsp),%xmm1
        movapd nb302_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb302_fixH2(%rsp),%xmm0
        movapd nb302_fiyH2(%rsp),%xmm1
        movapd nb302_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb302_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb302_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb302_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb302_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb302_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302_x86_64_sse2.nb302_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302_n(%rsp)
        jmp _nb_kernel302_x86_64_sse2.nb302_outer
_nb_kernel302_x86_64_sse2.nb302_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302_x86_64_sse2.nb302_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302_x86_64_sse2.nb302_threadloop
_nb_kernel302_x86_64_sse2.nb302_end: 
        movl nb302_nouter(%rsp),%eax
        movl nb302_ninner(%rsp),%ebx
        movq nb302_outeriter(%rbp),%rcx
        movq nb302_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1528,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel302nf_x86_64_sse2
.globl _nb_kernel302nf_x86_64_sse2
nb_kernel302nf_x86_64_sse2:     
_nb_kernel302nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb302nf_fshift, 16
.set nb302nf_gid, 24
.set nb302nf_pos, 32
.set nb302nf_faction, 40
.set nb302nf_charge, 48
.set nb302nf_p_facel, 56
.set nb302nf_argkrf, 64
.set nb302nf_argcrf, 72
.set nb302nf_Vc, 80
.set nb302nf_type, 88
.set nb302nf_p_ntype, 96
.set nb302nf_vdwparam, 104
.set nb302nf_Vvdw, 112
.set nb302nf_p_tabscale, 120
.set nb302nf_VFtab, 128
.set nb302nf_invsqrta, 136
.set nb302nf_dvda, 144
.set nb302nf_p_gbtabscale, 152
.set nb302nf_GBtab, 160
.set nb302nf_p_nthreads, 168
.set nb302nf_count, 176
.set nb302nf_mtx, 184
.set nb302nf_outeriter, 192
.set nb302nf_inneriter, 200
.set nb302nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb302nf_ixO, 0
.set nb302nf_iyO, 16
.set nb302nf_izO, 32
.set nb302nf_ixH1, 48
.set nb302nf_iyH1, 64
.set nb302nf_izH1, 80
.set nb302nf_ixH2, 96
.set nb302nf_iyH2, 112
.set nb302nf_izH2, 128
.set nb302nf_jxO, 144
.set nb302nf_jyO, 160
.set nb302nf_jzO, 176
.set nb302nf_jxH1, 192
.set nb302nf_jyH1, 208
.set nb302nf_jzH1, 224
.set nb302nf_jxH2, 240
.set nb302nf_jyH2, 256
.set nb302nf_jzH2, 272
.set nb302nf_qqOO, 288
.set nb302nf_qqOH, 304
.set nb302nf_qqHH, 320
.set nb302nf_tsc, 336
.set nb302nf_vctot, 352
.set nb302nf_half, 368
.set nb302nf_three, 384
.set nb302nf_rsqOO, 400
.set nb302nf_rsqOH1, 416
.set nb302nf_rsqOH2, 432
.set nb302nf_rsqH1O, 448
.set nb302nf_rsqH1H1, 464
.set nb302nf_rsqH1H2, 480
.set nb302nf_rsqH2O, 496
.set nb302nf_rsqH2H1, 512
.set nb302nf_rsqH2H2, 528
.set nb302nf_rinvOO, 544
.set nb302nf_rinvOH1, 560
.set nb302nf_rinvOH2, 576
.set nb302nf_rinvH1O, 592
.set nb302nf_rinvH1H1, 608
.set nb302nf_rinvH1H2, 624
.set nb302nf_rinvH2O, 640
.set nb302nf_rinvH2H1, 656
.set nb302nf_rinvH2H2, 672
.set nb302nf_is3, 688
.set nb302nf_ii3, 692
.set nb302nf_nri, 696
.set nb302nf_iinr, 704
.set nb302nf_jindex, 712
.set nb302nf_jjnr, 720
.set nb302nf_shift, 728
.set nb302nf_shiftvec, 736
.set nb302nf_facel, 744
.set nb302nf_innerjjnr, 752
.set nb302nf_innerk, 760
.set nb302nf_n, 764
.set nb302nf_nn1, 768
.set nb302nf_nouter, 772
.set nb302nf_ninner, 776
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $792,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb302nf_nouter(%rsp)
        movl %eax,nb302nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb302nf_nri(%rsp)
        movq %rsi,nb302nf_iinr(%rsp)
        movq %rdx,nb302nf_jindex(%rsp)
        movq %rcx,nb302nf_jjnr(%rsp)
        movq %r8,nb302nf_shift(%rsp)
        movq %r9,nb302nf_shiftvec(%rsp)
        movq nb302nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb302nf_facel(%rsp)

        movq nb302nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb302nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb302nf_half(%rsp)
        movl %ebx,nb302nf_half+4(%rsp)
        movsd nb302nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb302nf_half(%rsp)
        movapd %xmm3,nb302nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb302nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb302nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb302nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb302nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302nf_qqOO(%rsp)
        movapd %xmm4,nb302nf_qqOH(%rsp)
        movapd %xmm5,nb302nf_qqHH(%rsp)

_nb_kernel302nf_x86_64_sse2.nb302nf_threadloop: 
        movq  nb302nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel302nf_x86_64_sse2.nb302nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302nf_x86_64_sse2.nb302nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302nf_n(%rsp)
        movl %ebx,nb302nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302nf_x86_64_sse2.nb302nf_outerstart
        jmp _nb_kernel302nf_x86_64_sse2.nb302nf_end

_nb_kernel302nf_x86_64_sse2.nb302nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302nf_nouter(%rsp),%ebx
        movl %ebx,nb302nf_nouter(%rsp)

_nb_kernel302nf_x86_64_sse2.nb302nf_outer: 
        movq  nb302nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb302nf_is3(%rsp)            ## store is3 

        movq  nb302nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb302nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx            ## ebx =ii 

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb302nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb302nf_ii3(%rsp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%rax,%rbx,8),%xmm3
        addsd 8(%rax,%rbx,8),%xmm4
        addsd 16(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302nf_ixO(%rsp)
        movapd %xmm4,nb302nf_iyO(%rsp)
        movapd %xmm5,nb302nf_izO(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%rax,%rbx,8),%xmm0
        addsd 32(%rax,%rbx,8),%xmm1
        addsd 40(%rax,%rbx,8),%xmm2
        addsd 48(%rax,%rbx,8),%xmm3
        addsd 56(%rax,%rbx,8),%xmm4
        addsd 64(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb302nf_ixH1(%rsp)
        movapd %xmm1,nb302nf_iyH1(%rsp)
        movapd %xmm2,nb302nf_izH1(%rsp)
        movapd %xmm3,nb302nf_ixH2(%rsp)
        movapd %xmm4,nb302nf_iyH2(%rsp)
        movapd %xmm5,nb302nf_izH2(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb302nf_vctot(%rsp)

        movq  nb302nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb302nf_pos(%rbp),%rsi
        movq  nb302nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb302nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb302nf_ninner(%rsp),%ecx
        movl  %ecx,nb302nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb302nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel302nf_x86_64_sse2.nb302nf_unroll_loop
        jmp   _nb_kernel302nf_x86_64_sse2.nb302nf_checksingle
_nb_kernel302nf_x86_64_sse2.nb302nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb302nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb302nf_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb302nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movhpd (%rsi,%rbx,8),%xmm2
        movhpd 8(%rsi,%rbx,8),%xmm3
        movhpd 16(%rsi,%rbx,8),%xmm4
        movhpd 24(%rsi,%rbx,8),%xmm5
        movhpd 32(%rsi,%rbx,8),%xmm6
        movhpd 40(%rsi,%rbx,8),%xmm7
        movapd  %xmm2,nb302nf_jxO(%rsp)
        movapd  %xmm3,nb302nf_jyO(%rsp)
        movapd  %xmm4,nb302nf_jzO(%rsp)
        movapd  %xmm5,nb302nf_jxH1(%rsp)
        movapd  %xmm6,nb302nf_jyH1(%rsp)
        movapd  %xmm7,nb302nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 48(%rsi,%rbx,8),%xmm2
        movhpd 56(%rsi,%rbx,8),%xmm3
        movhpd 64(%rsi,%rbx,8),%xmm4
        movapd  %xmm2,nb302nf_jxH2(%rsp)
        movapd  %xmm3,nb302nf_jyH2(%rsp)
        movapd  %xmm4,nb302nf_jzH2(%rsp)

        movapd nb302nf_ixO(%rsp),%xmm0
        movapd nb302nf_iyO(%rsp),%xmm1
        movapd nb302nf_izO(%rsp),%xmm2
        movapd nb302nf_ixO(%rsp),%xmm3
        movapd nb302nf_iyO(%rsp),%xmm4
        movapd nb302nf_izO(%rsp),%xmm5
        subpd  nb302nf_jxO(%rsp),%xmm0
        subpd  nb302nf_jyO(%rsp),%xmm1
        subpd  nb302nf_jzO(%rsp),%xmm2
        subpd  nb302nf_jxH1(%rsp),%xmm3
        subpd  nb302nf_jyH1(%rsp),%xmm4
        subpd  nb302nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqOO(%rsp)
        movapd %xmm3,nb302nf_rsqOH1(%rsp)

        movapd nb302nf_ixO(%rsp),%xmm0
        movapd nb302nf_iyO(%rsp),%xmm1
        movapd nb302nf_izO(%rsp),%xmm2
        movapd nb302nf_ixH1(%rsp),%xmm3
        movapd nb302nf_iyH1(%rsp),%xmm4
        movapd nb302nf_izH1(%rsp),%xmm5
        subpd  nb302nf_jxH2(%rsp),%xmm0
        subpd  nb302nf_jyH2(%rsp),%xmm1
        subpd  nb302nf_jzH2(%rsp),%xmm2
        subpd  nb302nf_jxO(%rsp),%xmm3
        subpd  nb302nf_jyO(%rsp),%xmm4
        subpd  nb302nf_jzO(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqOH2(%rsp)
        movapd %xmm3,nb302nf_rsqH1O(%rsp)

        movapd nb302nf_ixH1(%rsp),%xmm0
        movapd nb302nf_iyH1(%rsp),%xmm1
        movapd nb302nf_izH1(%rsp),%xmm2
        movapd nb302nf_ixH1(%rsp),%xmm3
        movapd nb302nf_iyH1(%rsp),%xmm4
        movapd nb302nf_izH1(%rsp),%xmm5
        subpd  nb302nf_jxH1(%rsp),%xmm0
        subpd  nb302nf_jyH1(%rsp),%xmm1
        subpd  nb302nf_jzH1(%rsp),%xmm2
        subpd  nb302nf_jxH2(%rsp),%xmm3
        subpd  nb302nf_jyH2(%rsp),%xmm4
        subpd  nb302nf_jzH2(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqH1H1(%rsp)
        movapd %xmm3,nb302nf_rsqH1H2(%rsp)

        movapd nb302nf_ixH2(%rsp),%xmm0
        movapd nb302nf_iyH2(%rsp),%xmm1
        movapd nb302nf_izH2(%rsp),%xmm2
        movapd nb302nf_ixH2(%rsp),%xmm3
        movapd nb302nf_iyH2(%rsp),%xmm4
        movapd nb302nf_izH2(%rsp),%xmm5
        subpd  nb302nf_jxO(%rsp),%xmm0
        subpd  nb302nf_jyO(%rsp),%xmm1
        subpd  nb302nf_jzO(%rsp),%xmm2
        subpd  nb302nf_jxH1(%rsp),%xmm3
        subpd  nb302nf_jyH1(%rsp),%xmm4
        subpd  nb302nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb302nf_rsqH2O(%rsp)
        movapd %xmm4,nb302nf_rsqH2H1(%rsp)

        movapd nb302nf_ixH2(%rsp),%xmm0
        movapd nb302nf_iyH2(%rsp),%xmm1
        movapd nb302nf_izH2(%rsp),%xmm2
        subpd  nb302nf_jxH2(%rsp),%xmm0
        subpd  nb302nf_jyH2(%rsp),%xmm1
        subpd  nb302nf_jzH2(%rsp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb302nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb302nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH2H2(%rsp)
        movapd %xmm5,nb302nf_rinvH2H1(%rsp)

        movapd nb302nf_rsqOO(%rsp),%xmm0
        movapd nb302nf_rsqOH1(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb302nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb302nf_rinvOO(%rsp)
        movapd %xmm5,nb302nf_rinvOH1(%rsp)

        movapd nb302nf_rsqOH2(%rsp),%xmm0
        movapd nb302nf_rsqH1O(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb302nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvOH2(%rsp)
        movapd %xmm5,nb302nf_rinvH1O(%rsp)

        movapd nb302nf_rsqH1H1(%rsp),%xmm0
        movapd nb302nf_rsqH1H2(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb302nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH1H1(%rsp)
        movapd %xmm5,nb302nf_rinvH1H2(%rsp)

        movapd nb302nf_rsqH2O(%rsp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb302nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb302nf_rinvOO(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOO(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## O-H1 interaction 
        movapd nb302nf_rinvOH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOH1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## O-H2 interaction  
        movapd nb302nf_rinvOH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOH2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H1-O interaction 
        movapd nb302nf_rinvH1O(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1O(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H1-H1 interaction 
        movapd nb302nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H1-H2 interaction 
        movapd nb302nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H2-O interaction 
        movapd nb302nf_rinvH2O(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2O(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb302nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb302nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%rsp),%xmm5
    movapd %xmm5,nb302nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb302nf_innerk(%rsp)
        jl    _nb_kernel302nf_x86_64_sse2.nb302nf_checksingle
        jmp   _nb_kernel302nf_x86_64_sse2.nb302nf_unroll_loop
_nb_kernel302nf_x86_64_sse2.nb302nf_checksingle: 
        movl  nb302nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel302nf_x86_64_sse2.nb302nf_dosingle
        jmp   _nb_kernel302nf_x86_64_sse2.nb302nf_updateouterdata
_nb_kernel302nf_x86_64_sse2.nb302nf_dosingle: 
        movq  nb302nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb302nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## fetch j coordinates 
        movlpd (%rsi,%rax,8),%xmm2
        movlpd 8(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rax,8),%xmm4
        movlpd 24(%rsi,%rax,8),%xmm5
        movlpd 32(%rsi,%rax,8),%xmm6
        movlpd 40(%rsi,%rax,8),%xmm7
        movapd  %xmm2,nb302nf_jxO(%rsp)
        movapd  %xmm3,nb302nf_jyO(%rsp)
        movapd  %xmm4,nb302nf_jzO(%rsp)
        movapd  %xmm5,nb302nf_jxH1(%rsp)
        movapd  %xmm6,nb302nf_jyH1(%rsp)
        movapd  %xmm7,nb302nf_jzH1(%rsp)
        movlpd 48(%rsi,%rax,8),%xmm2
        movlpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movapd  %xmm2,nb302nf_jxH2(%rsp)
        movapd  %xmm3,nb302nf_jyH2(%rsp)
        movapd  %xmm4,nb302nf_jzH2(%rsp)

        movapd nb302nf_ixO(%rsp),%xmm0
        movapd nb302nf_iyO(%rsp),%xmm1
        movapd nb302nf_izO(%rsp),%xmm2
        movapd nb302nf_ixO(%rsp),%xmm3
        movapd nb302nf_iyO(%rsp),%xmm4
        movapd nb302nf_izO(%rsp),%xmm5
        subsd  nb302nf_jxO(%rsp),%xmm0
        subsd  nb302nf_jyO(%rsp),%xmm1
        subsd  nb302nf_jzO(%rsp),%xmm2
        subsd  nb302nf_jxH1(%rsp),%xmm3
        subsd  nb302nf_jyH1(%rsp),%xmm4
        subsd  nb302nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqOO(%rsp)
        movapd %xmm3,nb302nf_rsqOH1(%rsp)

        movapd nb302nf_ixO(%rsp),%xmm0
        movapd nb302nf_iyO(%rsp),%xmm1
        movapd nb302nf_izO(%rsp),%xmm2
        movapd nb302nf_ixH1(%rsp),%xmm3
        movapd nb302nf_iyH1(%rsp),%xmm4
        movapd nb302nf_izH1(%rsp),%xmm5
        subsd  nb302nf_jxH2(%rsp),%xmm0
        subsd  nb302nf_jyH2(%rsp),%xmm1
        subsd  nb302nf_jzH2(%rsp),%xmm2
        subsd  nb302nf_jxO(%rsp),%xmm3
        subsd  nb302nf_jyO(%rsp),%xmm4
        subsd  nb302nf_jzO(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqOH2(%rsp)
        movapd %xmm3,nb302nf_rsqH1O(%rsp)

        movapd nb302nf_ixH1(%rsp),%xmm0
        movapd nb302nf_iyH1(%rsp),%xmm1
        movapd nb302nf_izH1(%rsp),%xmm2
        movapd nb302nf_ixH1(%rsp),%xmm3
        movapd nb302nf_iyH1(%rsp),%xmm4
        movapd nb302nf_izH1(%rsp),%xmm5
        subsd  nb302nf_jxH1(%rsp),%xmm0
        subsd  nb302nf_jyH1(%rsp),%xmm1
        subsd  nb302nf_jzH1(%rsp),%xmm2
        subsd  nb302nf_jxH2(%rsp),%xmm3
        subsd  nb302nf_jyH2(%rsp),%xmm4
        subsd  nb302nf_jzH2(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302nf_rsqH1H1(%rsp)
        movapd %xmm3,nb302nf_rsqH1H2(%rsp)

        movapd nb302nf_ixH2(%rsp),%xmm0
        movapd nb302nf_iyH2(%rsp),%xmm1
        movapd nb302nf_izH2(%rsp),%xmm2
        movapd nb302nf_ixH2(%rsp),%xmm3
        movapd nb302nf_iyH2(%rsp),%xmm4
        movapd nb302nf_izH2(%rsp),%xmm5
        subsd  nb302nf_jxO(%rsp),%xmm0
        subsd  nb302nf_jyO(%rsp),%xmm1
        subsd  nb302nf_jzO(%rsp),%xmm2
        subsd  nb302nf_jxH1(%rsp),%xmm3
        subsd  nb302nf_jyH1(%rsp),%xmm4
        subsd  nb302nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb302nf_rsqH2O(%rsp)
        movapd %xmm4,nb302nf_rsqH2H1(%rsp)

        movapd nb302nf_ixH2(%rsp),%xmm0
        movapd nb302nf_iyH2(%rsp),%xmm1
        movapd nb302nf_izH2(%rsp),%xmm2
        subsd  nb302nf_jxH2(%rsp),%xmm0
        subsd  nb302nf_jyH2(%rsp),%xmm1
        subsd  nb302nf_jzH2(%rsp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb302nf_rsqH2H2(%rsp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb302nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH2H2(%rsp)
        movapd %xmm5,nb302nf_rinvH2H1(%rsp)

        movapd nb302nf_rsqOO(%rsp),%xmm0
        movapd nb302nf_rsqOH1(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb302nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%rsp),%xmm5   ## rinv
        movapd %xmm1,nb302nf_rinvOO(%rsp)
        movapd %xmm5,nb302nf_rinvOH1(%rsp)

        movapd nb302nf_rsqOH2(%rsp),%xmm0
        movapd nb302nf_rsqH1O(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb302nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvOH2(%rsp)
        movapd %xmm5,nb302nf_rinvH1O(%rsp)

        movapd nb302nf_rsqH1H1(%rsp),%xmm0
        movapd nb302nf_rsqH1H2(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb302nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH1H1(%rsp)
        movapd %xmm5,nb302nf_rinvH1H2(%rsp)

        movapd nb302nf_rsqH2O(%rsp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302nf_half(%rsp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb302nf_rinvH2O(%rsp)

        ## start with OO interaction 
        movapd nb302nf_rinvOO(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOO(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOO(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## O-H1 interaction 
        movapd nb302nf_rinvOH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOH1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## O-H2 interaction  
        movapd nb302nf_rinvOH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOH2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H1-O interaction 
        movapd nb302nf_rinvH1O(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1O(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H1-H1 interaction 
        movapd nb302nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H1-H2 interaction 
        movapd nb302nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H2-O interaction 
        movapd nb302nf_rinvH2O(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2O(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb302nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb302nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb302nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%rsp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%rsp)

_nb_kernel302nf_x86_64_sse2.nb302nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb302nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb302nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb302nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb302nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb302nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302nf_x86_64_sse2.nb302nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302nf_n(%rsp)
        jmp _nb_kernel302nf_x86_64_sse2.nb302nf_outer
_nb_kernel302nf_x86_64_sse2.nb302nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302nf_x86_64_sse2.nb302nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302nf_x86_64_sse2.nb302nf_threadloop
_nb_kernel302nf_x86_64_sse2.nb302nf_end: 
        movl nb302nf_nouter(%rsp),%eax
        movl nb302nf_ninner(%rsp),%ebx
        movq nb302nf_outeriter(%rbp),%rcx
        movq nb302nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $792,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



