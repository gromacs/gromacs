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







.globl nb_kernel103_x86_64_sse2
.globl _nb_kernel103_x86_64_sse2
nb_kernel103_x86_64_sse2:       
_nb_kernel103_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb103_fshift, 16
.set nb103_gid, 24
.set nb103_pos, 32
.set nb103_faction, 40
.set nb103_charge, 48
.set nb103_p_facel, 56
.set nb103_argkrf, 64
.set nb103_argcrf, 72
.set nb103_Vc, 80
.set nb103_type, 88
.set nb103_p_ntype, 96
.set nb103_vdwparam, 104
.set nb103_Vvdw, 112
.set nb103_p_tabscale, 120
.set nb103_VFtab, 128
.set nb103_invsqrta, 136
.set nb103_dvda, 144
.set nb103_p_gbtabscale, 152
.set nb103_GBtab, 160
.set nb103_p_nthreads, 168
.set nb103_count, 176
.set nb103_mtx, 184
.set nb103_outeriter, 192
.set nb103_inneriter, 200
.set nb103_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb103_ixH1, 0
.set nb103_iyH1, 16
.set nb103_izH1, 32
.set nb103_ixH2, 48
.set nb103_iyH2, 64
.set nb103_izH2, 80
.set nb103_ixM, 96
.set nb103_iyM, 112
.set nb103_izM, 128
.set nb103_iqM, 144
.set nb103_iqH, 160
.set nb103_dxH1, 176
.set nb103_dyH1, 192
.set nb103_dzH1, 208
.set nb103_dxH2, 224
.set nb103_dyH2, 240
.set nb103_dzH2, 256
.set nb103_dxM, 272
.set nb103_dyM, 288
.set nb103_dzM, 304
.set nb103_qqM, 320
.set nb103_qqH, 336
.set nb103_vctot, 352
.set nb103_fixM, 368
.set nb103_fiyM, 384
.set nb103_fizM, 400
.set nb103_fixH1, 416
.set nb103_fiyH1, 432
.set nb103_fizH1, 448
.set nb103_fixH2, 464
.set nb103_fiyH2, 480
.set nb103_fizH2, 496
.set nb103_fjx, 512
.set nb103_fjy, 528
.set nb103_fjz, 544
.set nb103_half, 560
.set nb103_three, 576
.set nb103_is3, 592
.set nb103_ii3, 596
.set nb103_nri, 600
.set nb103_innerjjnr, 608
.set nb103_iinr, 616
.set nb103_jindex, 624
.set nb103_jjnr, 632
.set nb103_shift, 640
.set nb103_shiftvec, 648
.set nb103_facel, 656
.set nb103_innerk, 664
.set nb103_n, 668
.set nb103_nn1, 672
.set nb103_nouter, 676
.set nb103_ninner, 680
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $696,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb103_nouter(%rsp)
        movl %eax,nb103_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb103_nri(%rsp)
        movq %rsi,nb103_iinr(%rsp)
        movq %rdx,nb103_jindex(%rsp)
        movq %rcx,nb103_jjnr(%rsp)
        movq %r8,nb103_shift(%rsp)
        movq %r9,nb103_shiftvec(%rsp)
        movq nb103_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb103_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb103_half(%rsp)
        movl %ebx,nb103_half+4(%rsp)
        movsd nb103_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb103_half(%rsp)
        movapd %xmm3,nb103_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb103_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb103_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb103_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb103_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb103_iqH(%rsp)
        movapd %xmm4,nb103_iqM(%rsp)

_nb_kernel103_x86_64_sse2.nb103_threadloop: 
        movq  nb103_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel103_x86_64_sse2.nb103_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103_x86_64_sse2.nb103_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103_n(%rsp)
        movl %ebx,nb103_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel103_x86_64_sse2.nb103_outerstart
        jmp _nb_kernel103_x86_64_sse2.nb103_end

_nb_kernel103_x86_64_sse2.nb103_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103_nouter(%rsp),%ebx
        movl %ebx,nb103_nouter(%rsp)

_nb_kernel103_x86_64_sse2.nb103_outer: 
        movq  nb103_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb103_is3(%rsp)      ## store is3 

        movq  nb103_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb103_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb103_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb103_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb103_ixH1(%rsp)
        movapd %xmm4,nb103_iyH1(%rsp)
        movapd %xmm5,nb103_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb103_ixH2(%rsp)
        movapd %xmm1,nb103_iyH2(%rsp)
        movapd %xmm2,nb103_izH2(%rsp)
        movapd %xmm3,nb103_ixM(%rsp)
        movapd %xmm4,nb103_iyM(%rsp)
        movapd %xmm5,nb103_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb103_vctot(%rsp)
        movapd %xmm4,nb103_fixM(%rsp)
        movapd %xmm4,nb103_fiyM(%rsp)
        movapd %xmm4,nb103_fizM(%rsp)
        movapd %xmm4,nb103_fixH1(%rsp)
        movapd %xmm4,nb103_fiyH1(%rsp)
        movapd %xmm4,nb103_fizH1(%rsp)
        movapd %xmm4,nb103_fixH2(%rsp)
        movapd %xmm4,nb103_fiyH2(%rsp)
        movapd %xmm4,nb103_fizH2(%rsp)

        movq  nb103_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx     ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb103_pos(%rbp),%rsi
        movq  nb103_faction(%rbp),%rdi
        movq  nb103_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb103_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb103_ninner(%rsp),%ecx
        movl  %ecx,nb103_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb103_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel103_x86_64_sse2.nb103_unroll_loop
        jmp   _nb_kernel103_x86_64_sse2.nb103_checksingle
_nb_kernel103_x86_64_sse2.nb103_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb103_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb103_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb103_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm6      ## jq A 
        movhpd (%rsi,%rbx,8),%xmm6      ## jq B 
        movapd nb103_iqM(%rsp),%xmm3
        movapd nb103_iqH(%rsp),%xmm4
        mulpd %xmm6,%xmm3               ## qqM 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103_qqM(%rsp)
        movapd  %xmm4,nb103_qqH(%rsp)

        movq nb103_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb103_ixH1(%rsp),%xmm0
    subpd nb103_iyH1(%rsp),%xmm1
    subpd nb103_izH1(%rsp),%xmm2
    subpd nb103_ixH2(%rsp),%xmm3
    subpd nb103_iyH2(%rsp),%xmm4
    subpd nb103_izH2(%rsp),%xmm5
    subpd nb103_ixM(%rsp),%xmm6
    subpd nb103_iyM(%rsp),%xmm7
    subpd nb103_izM(%rsp),%xmm8

        movapd %xmm0,nb103_dxH1(%rsp)
        movapd %xmm1,nb103_dyH1(%rsp)
        movapd %xmm2,nb103_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb103_dxH2(%rsp)
        movapd %xmm4,nb103_dyH2(%rsp)
        movapd %xmm5,nb103_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb103_dxM(%rsp)
        movapd %xmm7,nb103_dyM(%rsp)
        movapd %xmm8,nb103_dzM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
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

        movapd  nb103_three(%rsp),%xmm9
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

        movapd  nb103_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb103_three(%rsp),%xmm1
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

        movapd  nb103_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        ## interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb103_qqH(%rsp),%xmm0
    mulpd  nb103_qqH(%rsp),%xmm1
    mulpd  nb103_qqM(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb103_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb103_vctot(%rsp)

    ## move j forces to xmm0-xmm2
        movq  nb103_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

        mulpd nb103_dxH1(%rsp),%xmm7
        mulpd nb103_dyH1(%rsp),%xmm8
        mulpd nb103_dzH1(%rsp),%xmm9
        mulpd nb103_dxH2(%rsp),%xmm10
        mulpd nb103_dyH2(%rsp),%xmm11
        mulpd nb103_dzH2(%rsp),%xmm12
        mulpd nb103_dxM(%rsp),%xmm13
        mulpd nb103_dyM(%rsp),%xmm14
        mulpd nb103_dzM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb103_fixH1(%rsp),%xmm7
    addpd nb103_fiyH1(%rsp),%xmm8
    addpd nb103_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb103_fixH2(%rsp),%xmm10
    addpd nb103_fiyH2(%rsp),%xmm11
    addpd nb103_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb103_fixM(%rsp),%xmm13
    addpd nb103_fiyM(%rsp),%xmm14
    addpd nb103_fizM(%rsp),%xmm15

    movapd %xmm7,nb103_fixH1(%rsp)
    movapd %xmm8,nb103_fiyH1(%rsp)
    movapd %xmm9,nb103_fizH1(%rsp)
    movapd %xmm10,nb103_fixH2(%rsp)
    movapd %xmm11,nb103_fiyH2(%rsp)
    movapd %xmm12,nb103_fizH2(%rsp)
    movapd %xmm13,nb103_fixM(%rsp)
    movapd %xmm14,nb103_fiyM(%rsp)
    movapd %xmm15,nb103_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb103_innerk(%rsp)
        jl    _nb_kernel103_x86_64_sse2.nb103_checksingle
        jmp   _nb_kernel103_x86_64_sse2.nb103_unroll_loop
_nb_kernel103_x86_64_sse2.nb103_checksingle:    
        movl  nb103_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel103_x86_64_sse2.nb103_dosingle
        jmp    _nb_kernel103_x86_64_sse2.nb103_updateouterdata
_nb_kernel103_x86_64_sse2.nb103_dosingle: 
        movq  nb103_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb103_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%rsi,%rax,8),%xmm6      ## jq A 

        movapd nb103_iqM(%rsp),%xmm3
        movapd nb103_iqH(%rsp),%xmm4
        mulsd %xmm6,%xmm3               ## qqM
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103_qqM(%rsp)
        movapd  %xmm4,nb103_qqH(%rsp)

        movq nb103_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm4-xmm6 & xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb103_ixM(%rsp),%xmm4
        subsd nb103_iyM(%rsp),%xmm5
        subsd nb103_izM(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxM(%rsp)
        movapd %xmm5,nb103_dyM(%rsp)
        movapd %xmm6,nb103_dzM(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb103_ixH1(%rsp),%xmm4
        subsd nb103_iyH1(%rsp),%xmm5
        subsd nb103_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb103_dxH1(%rsp)
        movapd %xmm5,nb103_dyH1(%rsp)
        movapd %xmm6,nb103_dzH1(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move j coords to xmm3-xmm5 
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## calc dr 
        subsd nb103_ixH2(%rsp),%xmm3
        subsd nb103_iyH2(%rsp),%xmm4
        subsd nb103_izH2(%rsp),%xmm5

        ## store dr 
        movapd %xmm3,nb103_dxH2(%rsp)
        movapd %xmm4,nb103_dyH2(%rsp)
        movapd %xmm5,nb103_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        movapd  %xmm7,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm7=rinv, xmm4=rinvsq 
        mulsd  nb103_qqM(%rsp),%xmm7    ## xmm7=vcoul 

        mulsd  %xmm7,%xmm4      ## total fsM in xmm4 

        addsd  nb103_vctot(%rsp),%xmm7

        movlpd %xmm7,nb103_vctot(%rsp)

        movapd nb103_dxM(%rsp),%xmm0
        movapd nb103_dyM(%rsp),%xmm1
        movapd nb103_dzM(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb103_fixM(%rsp),%xmm3
        movapd nb103_fiyM(%rsp),%xmm4
        movapd nb103_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixM(%rsp)
        movlpd %xmm4,nb103_fiyM(%rsp)
        movlpd %xmm7,nb103_fizM(%rsp)
        ## update j forces with water M 
        movlpd %xmm0,nb103_fjx(%rsp)
        movlpd %xmm1,nb103_fjy(%rsp)
        movlpd %xmm2,nb103_fjz(%rsp)

        ## H1 interactions 
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb103_qqH(%rsp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb103_vctot(%rsp),%xmm6

        movapd nb103_dxH1(%rsp),%xmm0
        movapd nb103_dyH1(%rsp),%xmm1
        movapd nb103_dzH1(%rsp),%xmm2
        movlpd %xmm6,nb103_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb103_fixH1(%rsp),%xmm3
        movapd nb103_fiyH1(%rsp),%xmm4
        movapd nb103_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixH1(%rsp)
        movlpd %xmm4,nb103_fiyH1(%rsp)
        movlpd %xmm7,nb103_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb103_fjx(%rsp),%xmm0
        addsd  nb103_fjy(%rsp),%xmm1
        addsd  nb103_fjz(%rsp),%xmm2
        movsd %xmm0,nb103_fjx(%rsp)
        movsd %xmm1,nb103_fjy(%rsp)
        movsd %xmm2,nb103_fjz(%rsp)

        ## H2 interactions 
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb103_qqH(%rsp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb103_vctot(%rsp),%xmm5

        movapd nb103_dxH2(%rsp),%xmm0
        movapd nb103_dyH2(%rsp),%xmm1
        movapd nb103_dzH2(%rsp),%xmm2
        movlpd %xmm5,nb103_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb103_fixH2(%rsp),%xmm3
        movapd nb103_fiyH2(%rsp),%xmm4
        movapd nb103_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movlpd %xmm3,nb103_fixH2(%rsp)
        movlpd %xmm4,nb103_fiyH2(%rsp)
        movlpd %xmm7,nb103_fizH2(%rsp)

        movq nb103_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb103_fjx(%rsp),%xmm0
        addsd  nb103_fjy(%rsp),%xmm1
        addsd  nb103_fjz(%rsp),%xmm2

        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel103_x86_64_sse2.nb103_updateouterdata: 
        movl  nb103_ii3(%rsp),%ecx
        movq  nb103_faction(%rbp),%rdi
        movq  nb103_fshift(%rbp),%rsi
        movl  nb103_is3(%rsp),%edx

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixH1(%rsp),%xmm0
        movapd nb103_fiyH1(%rsp),%xmm1
        movapd nb103_fizH1(%rsp),%xmm2

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
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixH2(%rsp),%xmm0
        movapd nb103_fiyH2(%rsp),%xmm1
        movapd nb103_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

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

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb103_fixM(%rsp),%xmm0
        movapd nb103_fiyM(%rsp),%xmm1
        movapd nb103_fizM(%rsp),%xmm2

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
        movsd  72(%rdi,%rcx,8),%xmm3
        movsd  80(%rdi,%rcx,8),%xmm4
        movsd  88(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,72(%rdi,%rcx,8)
        movsd  %xmm4,80(%rdi,%rcx,8)
        movsd  %xmm5,88(%rdi,%rcx,8)

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
        movl nb103_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb103_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb103_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb103_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb103_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel103_x86_64_sse2.nb103_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103_n(%rsp)
        jmp _nb_kernel103_x86_64_sse2.nb103_outer
_nb_kernel103_x86_64_sse2.nb103_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103_x86_64_sse2.nb103_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103_x86_64_sse2.nb103_threadloop
_nb_kernel103_x86_64_sse2.nb103_end: 
        movl nb103_nouter(%rsp),%eax
        movl nb103_ninner(%rsp),%ebx
        movq nb103_outeriter(%rbp),%rcx
        movq nb103_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $696,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



.globl nb_kernel103nf_x86_64_sse2
.globl _nb_kernel103nf_x86_64_sse2
nb_kernel103nf_x86_64_sse2:     
_nb_kernel103nf_x86_64_sse2:    
.set nb103nf_fshift, 16
.set nb103nf_gid, 24
.set nb103nf_pos, 32
.set nb103nf_faction, 40
.set nb103nf_charge, 48
.set nb103nf_p_facel, 56
.set nb103nf_argkrf, 64
.set nb103nf_argcrf, 72
.set nb103nf_Vc, 80
.set nb103nf_type, 88
.set nb103nf_p_ntype, 96
.set nb103nf_vdwparam, 104
.set nb103nf_Vvdw, 112
.set nb103nf_p_tabscale, 120
.set nb103nf_VFtab, 128
.set nb103nf_invsqrta, 136
.set nb103nf_dvda, 144
.set nb103nf_p_gbtabscale, 152
.set nb103nf_GBtab, 160
.set nb103nf_p_nthreads, 168
.set nb103nf_count, 176
.set nb103nf_mtx, 184
.set nb103nf_outeriter, 192
.set nb103nf_inneriter, 200
.set nb103nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb103nf_ixM, 0
.set nb103nf_iyM, 16
.set nb103nf_izM, 32
.set nb103nf_ixH1, 48
.set nb103nf_iyH1, 64
.set nb103nf_izH1, 80
.set nb103nf_ixH2, 96
.set nb103nf_iyH2, 112
.set nb103nf_izH2, 128
.set nb103nf_iqM, 144
.set nb103nf_iqH, 160
.set nb103nf_qqM, 176
.set nb103nf_qqH, 192
.set nb103nf_vctot, 208
.set nb103nf_half, 224
.set nb103nf_three, 240
.set nb103nf_is3, 256
.set nb103nf_ii3, 260
.set nb103nf_nri, 264
.set nb103nf_iinr, 272
.set nb103nf_jindex, 280
.set nb103nf_jjnr, 288
.set nb103nf_shift, 296
.set nb103nf_shiftvec, 304
.set nb103nf_facel, 312
.set nb103nf_innerjjnr, 320
.set nb103nf_innerk, 328
.set nb103nf_n, 332
.set nb103nf_nn1, 336
.set nb103nf_nouter, 340
.set nb103nf_ninner, 344
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
        movl %eax,nb103nf_nouter(%rsp)
        movl %eax,nb103nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb103nf_nri(%rsp)
        movq %rsi,nb103nf_iinr(%rsp)
        movq %rdx,nb103nf_jindex(%rsp)
        movq %rcx,nb103nf_jjnr(%rsp)
        movq %r8,nb103nf_shift(%rsp)
        movq %r9,nb103nf_shiftvec(%rsp)
        movq nb103nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb103nf_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb103nf_half(%rsp)
        movl %ebx,nb103nf_half+4(%rsp)
        movsd nb103nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb103nf_half(%rsp)
        movapd %xmm3,nb103nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb103nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb103nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb103nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb103nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb103nf_iqH(%rsp)
        movapd %xmm4,nb103nf_iqM(%rsp)

_nb_kernel103nf_x86_64_sse2.nb103nf_threadloop: 
        movq  nb103nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel103nf_x86_64_sse2.nb103nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103nf_x86_64_sse2.nb103nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103nf_n(%rsp)
        movl %ebx,nb103nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel103nf_x86_64_sse2.nb103nf_outerstart
        jmp _nb_kernel103nf_x86_64_sse2.nb103nf_end

_nb_kernel103nf_x86_64_sse2.nb103nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103nf_nouter(%rsp),%ebx
        movl %ebx,nb103nf_nouter(%rsp)

_nb_kernel103nf_x86_64_sse2.nb103nf_outer: 
        movq  nb103nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb103nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb103nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb103nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb103nf_ii3(%rsp)

        addsd 24(%rax,%rbx,8),%xmm3
        addsd 32(%rax,%rbx,8),%xmm4
        addsd 40(%rax,%rbx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb103nf_ixH1(%rsp)
        movapd %xmm4,nb103nf_iyH1(%rsp)
        movapd %xmm5,nb103nf_izH1(%rsp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 48(%rax,%rbx,8),%xmm0
        addsd 56(%rax,%rbx,8),%xmm1
        addsd 64(%rax,%rbx,8),%xmm2
        addsd 72(%rax,%rbx,8),%xmm3
        addsd 80(%rax,%rbx,8),%xmm4
        addsd 88(%rax,%rbx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb103nf_ixH2(%rsp)
        movapd %xmm1,nb103nf_iyH2(%rsp)
        movapd %xmm2,nb103nf_izH2(%rsp)
        movapd %xmm3,nb103nf_ixM(%rsp)
        movapd %xmm4,nb103nf_iyM(%rsp)
        movapd %xmm5,nb103nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb103nf_vctot(%rsp)

        movq  nb103nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb103nf_pos(%rbp),%rsi
        movq  nb103nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb103nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb103nf_ninner(%rsp),%ecx
        movl  %ecx,nb103nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb103nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel103nf_x86_64_sse2.nb103nf_unroll_loop
        jmp   _nb_kernel103nf_x86_64_sse2.nb103nf_checksingle
_nb_kernel103nf_x86_64_sse2.nb103nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb103nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb103nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb103nf_charge(%rbp),%rsi     ## base of charge[] 


        movlpd (%rsi,%rax,8),%xmm6      ## jq A 
        movhpd (%rsi,%rbx,8),%xmm6      ## jq B 
        movapd nb103nf_iqM(%rsp),%xmm3
        movapd nb103nf_iqH(%rsp),%xmm4
        mulpd %xmm6,%xmm3               ## qqM 
        mulpd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103nf_qqM(%rsp)
        movapd  %xmm4,nb103nf_qqH(%rsp)

        movq nb103nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103nf_ixM(%rsp),%xmm4
        movapd nb103nf_iyM(%rsp),%xmm5
        movapd nb103nf_izM(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103nf_ixH1(%rsp),%xmm4
        movapd nb103nf_iyH1(%rsp),%xmm5
        movapd nb103nf_izH1(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb103nf_ixH2(%rsp),%xmm3
        movapd nb103nf_iyH2(%rsp),%xmm4
        movapd nb103nf_izH2(%rsp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtpd2ps %xmm7,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulpd   %xmm7,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subpd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subpd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulpd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulpd %xmm4,%xmm4       ## lu*lu 
        mulpd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subpd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulpd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        mulpd  nb103nf_qqM(%rsp),%xmm7          ## xmm7=vcoul 
        addpd  nb103nf_vctot(%rsp),%xmm7
        movapd %xmm7,nb103nf_vctot(%rsp)

        ## H1 interactions 
        mulpd  nb103nf_qqH(%rsp),%xmm6          ## xmm6=vcoul 
        addpd  nb103nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb103nf_vctot(%rsp)

        ## H2 interactions 
        mulpd  nb103nf_qqH(%rsp),%xmm5          ## xmm5=vcoul 
        addpd  nb103nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb103nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb103nf_innerk(%rsp)
        jl    _nb_kernel103nf_x86_64_sse2.nb103nf_checksingle
        jmp   _nb_kernel103nf_x86_64_sse2.nb103nf_unroll_loop
_nb_kernel103nf_x86_64_sse2.nb103nf_checksingle: 
        movl  nb103nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel103nf_x86_64_sse2.nb103nf_dosingle
        jmp   _nb_kernel103nf_x86_64_sse2.nb103nf_updateouterdata
_nb_kernel103nf_x86_64_sse2.nb103nf_dosingle: 
        movq  nb103nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb103nf_charge(%rbp),%rsi     ## base of charge[] 
        xorpd %xmm6,%xmm6
        movlpd (%rsi,%rax,8),%xmm6      ## jq A 

        movapd nb103nf_iqM(%rsp),%xmm3
        movapd nb103nf_iqH(%rsp),%xmm4
        mulsd %xmm6,%xmm3               ## qqM 
        mulsd %xmm6,%xmm4               ## qqH 

        movapd  %xmm3,nb103nf_qqM(%rsp)
        movapd  %xmm4,nb103nf_qqH(%rsp)

        movq nb103nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixM-izM to xmm4-xmm6 
        movapd nb103nf_ixM(%rsp),%xmm4
        movapd nb103nf_iyM(%rsp),%xmm5
        movapd nb103nf_izM(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqM in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb103nf_ixH1(%rsp),%xmm4
        movapd nb103nf_iyH1(%rsp),%xmm5
        movapd nb103nf_izH1(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb103nf_ixH2(%rsp),%xmm3
        movapd nb103nf_iyH2(%rsp),%xmm4
        movapd nb103nf_izH2(%rsp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

        ## start with rsqM - put seed in xmm2 
        cvtsd2ss %xmm7,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulsd   %xmm7,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm7       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subsd %xmm7,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm7     ## rinvM in xmm7 

        ## rsqH1 - seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm6       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subsd %xmm6,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm6     ## rinvH1 in xmm6 

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb103nf_three(%rsp),%xmm4
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulsd   nb103nf_half(%rsp),%xmm4   ## iter1 ( new lu) 

        movapd %xmm4,%xmm3
        mulsd %xmm4,%xmm4       ## lu*lu 
        mulsd %xmm4,%xmm5       ## rsq*lu*lu 
        movapd nb103nf_three(%rsp),%xmm4
        subsd %xmm5,%xmm4       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm4       ## lu*( 3-rsq*lu*lu) 
        mulsd nb103nf_half(%rsp),%xmm4   ## rinv 
        movapd  %xmm4,%xmm5     ## rinvH2 in xmm5 

        ## do M interactions 
        mulsd  nb103nf_qqM(%rsp),%xmm7          ## xmm7=vcoul 
        addsd  nb103nf_vctot(%rsp),%xmm7
        movlpd %xmm7,nb103nf_vctot(%rsp)

        ## H1 interactions 
        mulsd  nb103nf_qqH(%rsp),%xmm6          ## xmm6=vcoul 
        addsd  nb103nf_vctot(%rsp),%xmm6
        movlpd %xmm6,nb103nf_vctot(%rsp)

        ## H2 interactions 
        mulsd  nb103nf_qqH(%rsp),%xmm5          ## xmm5=vcoul 
        addsd  nb103nf_vctot(%rsp),%xmm5
        movlpd %xmm5,nb103nf_vctot(%rsp)

_nb_kernel103nf_x86_64_sse2.nb103nf_updateouterdata: 
        ## get n from stack
        movl nb103nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb103nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        movapd nb103nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 

        ## add earlier value from mem 
        movq  nb103nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb103nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel103nf_x86_64_sse2.nb103nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103nf_n(%rsp)
        jmp _nb_kernel103nf_x86_64_sse2.nb103nf_outer
_nb_kernel103nf_x86_64_sse2.nb103nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103nf_x86_64_sse2.nb103nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103nf_x86_64_sse2.nb103nf_threadloop
_nb_kernel103nf_x86_64_sse2.nb103nf_end: 
        movl nb103nf_nouter(%rsp),%eax
        movl nb103nf_ninner(%rsp),%ebx
        movq nb103nf_outeriter(%rbp),%rcx
        movq nb103nf_inneriter(%rbp),%rdx
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


