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







.globl nb_kernel310_x86_64_sse2
.globl _nb_kernel310_x86_64_sse2
nb_kernel310_x86_64_sse2:       
_nb_kernel310_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb310_fshift, 16
.set nb310_gid, 24
.set nb310_pos, 32
.set nb310_faction, 40
.set nb310_charge, 48
.set nb310_p_facel, 56
.set nb310_argkrf, 64
.set nb310_argcrf, 72
.set nb310_Vc, 80
.set nb310_type, 88
.set nb310_p_ntype, 96
.set nb310_vdwparam, 104
.set nb310_Vvdw, 112
.set nb310_p_tabscale, 120
.set nb310_VFtab, 128
.set nb310_invsqrta, 136
.set nb310_dvda, 144
.set nb310_p_gbtabscale, 152
.set nb310_GBtab, 160
.set nb310_p_nthreads, 168
.set nb310_count, 176
.set nb310_mtx, 184
.set nb310_outeriter, 192
.set nb310_inneriter, 200
.set nb310_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb310_ix, 0
.set nb310_iy, 16
.set nb310_iz, 32
.set nb310_iq, 48
.set nb310_dx, 64
.set nb310_dy, 80
.set nb310_dz, 96
.set nb310_two, 112
.set nb310_six, 128
.set nb310_twelve, 144
.set nb310_tsc, 160
.set nb310_qq, 176
.set nb310_c6, 192
.set nb310_c12, 208
.set nb310_eps, 224
.set nb310_vctot, 240
.set nb310_Vvdwtot, 256
.set nb310_fix, 272
.set nb310_fiy, 288
.set nb310_fiz, 304
.set nb310_half, 320
.set nb310_three, 336
.set nb310_nri, 352
.set nb310_iinr, 360
.set nb310_jindex, 368
.set nb310_jjnr, 376
.set nb310_shift, 384
.set nb310_shiftvec, 392
.set nb310_facel, 400
.set nb310_innerjjnr, 408
.set nb310_is3, 416
.set nb310_ii3, 420
.set nb310_ntia, 424
.set nb310_innerk, 428
.set nb310_n, 432
.set nb310_nn1, 436
.set nb310_ntype, 440
.set nb310_nouter, 444
.set nb310_ninner, 448
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $472,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb310_nouter(%rsp)
        movl %eax,nb310_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb310_nri(%rsp)
        movq %rsi,nb310_iinr(%rsp)
        movq %rdx,nb310_jindex(%rsp)
        movq %rcx,nb310_jjnr(%rsp)
        movq %r8,nb310_shift(%rsp)
        movq %r9,nb310_shiftvec(%rsp)
        movq nb310_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb310_ntype(%rsp)
        movq nb310_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb310_facel(%rsp)

        movq nb310_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb310_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb310_half(%rsp)
        movl %ebx,nb310_half+4(%rsp)
        movsd nb310_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb310_half(%rsp)
        movapd %xmm2,nb310_two(%rsp)
        movapd %xmm3,nb310_three(%rsp)
        movapd %xmm4,nb310_six(%rsp)
        movapd %xmm5,nb310_twelve(%rsp)

_nb_kernel310_x86_64_sse2.nb310_threadloop: 
        movq  nb310_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel310_x86_64_sse2.nb310_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310_x86_64_sse2.nb310_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310_n(%rsp)
        movl %ebx,nb310_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310_x86_64_sse2.nb310_outerstart
        jmp _nb_kernel310_x86_64_sse2.nb310_end

_nb_kernel310_x86_64_sse2.nb310_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310_nouter(%rsp),%ebx
        movl %ebx,nb310_nouter(%rsp)

_nb_kernel310_x86_64_sse2.nb310_outer: 
        movq  nb310_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb310_is3(%rsp)      ## store is3 

        movq  nb310_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb310_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb310_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb310_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb310_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb310_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb310_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb310_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb310_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb310_ix(%rsp)
        movapd %xmm1,nb310_iy(%rsp)
        movapd %xmm2,nb310_iz(%rsp)

        movl  %ebx,nb310_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm15,%xmm15
        movapd %xmm15,nb310_vctot(%rsp)
        movapd %xmm15,nb310_Vvdwtot(%rsp)
        movapd %xmm15,%xmm14
        movapd %xmm15,%xmm13

        movq  nb310_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb310_pos(%rbp),%rsi
        movq  nb310_faction(%rbp),%rdi
        movq  nb310_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb310_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb310_ninner(%rsp),%ecx
        movl  %ecx,nb310_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb310_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel310_x86_64_sse2.nb310_unroll_loop
        jmp   _nb_kernel310_x86_64_sse2.nb310_checksingle
_nb_kernel310_x86_64_sse2.nb310_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb310_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r10d
        movl  4(%rdx),%r11d
        addq $8,nb310_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb310_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r10,%r10,2),%rax     ## replace jnr with j3 
        lea  (%r11,%r11,2),%rbx

        ## move two coordinates to xmm4-xmm6    
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb310_ix(%rsp),%xmm4
        subpd nb310_iy(%rsp),%xmm5
        subpd nb310_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb310_charge(%rbp),%rsi     ## base of charge[] 
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4

        movlpd (%rsi,%r10,8),%xmm3
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        movhpd (%rsi,%r11,8),%xmm3

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb310_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb310_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        mulpd  nb310_iq(%rsp),%xmm3
        movapd %xmm3,nb310_qq(%rsp)

        movq nb310_type(%rbp),%rsi

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb310_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb310_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
    movapd %xmm0,%xmm1  ## xmm1=rinv

        movl (%rsi,%r10,4),%r12d
        movl (%rsi,%r11,4),%r13d

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb310_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,nb310_eps(%rsp)
        shll %r12d
        shll %r13d

        pslld $2,%mm6           ## idx *= 4 

        movq nb310_VFtab(%rbp),%rsi
        movd %mm6,%r8d
        psrlq $32,%mm6
        movd %mm6,%r9d
        movl nb310_ntia(%rsp),%edi
        addl %edi,%r12d
        addl %edi,%r13d

        movq nb310_vdwparam(%rbp),%rdi

        movapd (%rsi,%r8,8),%xmm4       ## Y1 F1        
        movapd (%rsi,%r9,8),%xmm8       ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm8,%xmm4    ## Y1 Y2 
        unpckhpd %xmm8,%xmm5    ## F1 F2 

        movlpd (%rdi,%r12,8),%xmm3
        movlpd 8(%rdi,%r12,8),%xmm7

    movapd %xmm1,%xmm0 ## rinv
    mulpd  %xmm0,%xmm0 ## rinvsq
    movapd %xmm0,%xmm2 ## rinvsq
    mulpd  %xmm2,%xmm2 ## rinv4
    mulpd  %xmm0,%xmm2 ## rinv6
    movapd %xmm2,%xmm12
    mulpd  %xmm12,%xmm12 ## rinv12

        movhpd (%rdi,%r13,8),%xmm3
        movhpd 8(%rdi,%r13,8),%xmm7

        movapd 16(%rsi,%r8,8),%xmm6     ## G1 H1        
        movapd 16(%rsi,%r9,8),%xmm8     ## G2 H2 

    mulpd  %xmm3,%xmm2   ## vvdw6=c6*rinv6
        mulpd  %xmm7,%xmm12  ## vvdw12=c12*rinv12     

        movapd %xmm12,%xmm0
        subpd  %xmm2,%xmm12     ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addpd  nb310_Vvdwtot(%rsp),%xmm12
    movapd %xmm12,nb310_Vvdwtot(%rsp)

        movapd %xmm6,%xmm7
        unpcklpd %xmm8,%xmm6    ## G1 G2 
        unpckhpd %xmm8,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            

        movapd nb310_eps(%rsp),%xmm3

    mulpd %xmm3,%xmm7  ## Heps
    mulpd  %xmm3,%xmm6 ## Geps
    mulpd %xmm3,%xmm7  ## Heps2

    addpd  %xmm6,%xmm5  ## F+Geps
    addpd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm7  ## 2*Heps2
    addpd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addpd  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulpd  %xmm3,%xmm5  ## eps*Fp
    addpd  %xmm4,%xmm5  ## VV
    mulpd  nb310_qq(%rsp),%xmm5     ## VV*qq=vcoul
    mulpd  nb310_qq(%rsp),%xmm7     ## FF*qq=fijC

        ## the fj's - start by accumulating forces from memory 
    movq nb310_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm6
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm6

    ## LJ forces
    mulpd  nb310_six(%rsp),%xmm2
    mulpd  nb310_twelve(%rsp),%xmm0
    subpd  %xmm2,%xmm0
    mulpd  %xmm1,%xmm0 ## (12*vnb12-6*vnb6)*rinv

    ## add potential to vctot 
        addpd  nb310_vctot(%rsp),%xmm5
    movapd %xmm5,nb310_vctot(%rsp)

    mulpd  nb310_tsc(%rsp),%xmm7
    subpd  %xmm7,%xmm0

    mulpd  %xmm1,%xmm0 ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulpd  %xmm0,%xmm9
        mulpd  %xmm0,%xmm10
        mulpd  %xmm0,%xmm11

        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm6

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm6,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm6,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb310_innerk(%rsp)
        jl    _nb_kernel310_x86_64_sse2.nb310_checksingle
        jmp   _nb_kernel310_x86_64_sse2.nb310_unroll_loop
_nb_kernel310_x86_64_sse2.nb310_checksingle: 
        movl  nb310_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel310_x86_64_sse2.nb310_dosingle
        jmp    _nb_kernel310_x86_64_sse2.nb310_updateouterdata
_nb_kernel310_x86_64_sse2.nb310_dosingle: 
        movq nb310_charge(%rbp),%rsi
        movq nb310_pos(%rbp),%rdi
        movq  nb310_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb310_charge(%rbp),%rsi     ## base of charge[] 
        movsd (%rsi,%rax,8),%xmm3
        mulsd  nb310_iq(%rsp),%xmm3
        movapd %xmm3,nb310_qq(%rsp)

        movq nb310_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb310_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb310_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb310_c6(%rsp)
        movapd %xmm6,nb310_c12(%rsp)

        movq nb310_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6    
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb310_ix(%rsp),%xmm4
        subsd nb310_iy(%rsp),%xmm5
        subsd nb310_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb310_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb310_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb310_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb310_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
    movapd %xmm0,%xmm1  ## xmm1=rinv

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb310_tsc(%rsp),%xmm4

        cvttsd2si %xmm4,%r8d    ## mm6 = lu idx 
        cvtsi2sd %r8d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm3      ## xmm3=eps 

        shll $2,%r8d            ## idx *= 4 

        movq nb310_VFtab(%rbp),%rsi

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm5
        movsd 16(%rsi,%r8,8),%xmm6
        movsd 24(%rsi,%r8,8),%xmm7

    movapd %xmm1,%xmm0 ## rinv
    mulsd  %xmm0,%xmm0 ## rinvsq
    movapd %xmm0,%xmm2 ## rinvsq
    mulsd  %xmm2,%xmm2 ## rinv4
    mulsd  %xmm0,%xmm2 ## rinv6
    movapd %xmm2,%xmm12
    mulsd  %xmm12,%xmm12 ## rinv12

    mulsd  nb310_c6(%rsp),%xmm2      ## vvdw6=c6*rinv6
        mulsd  nb310_c12(%rsp),%xmm12     ## vvdw12=c12*rinv12     

        movapd %xmm12,%xmm0
        subsd  %xmm2,%xmm12     ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addsd  nb310_Vvdwtot(%rsp),%xmm12
    movsd %xmm12,nb310_Vvdwtot(%rsp)

        ## coulomb table ready, in xmm4-xmm7            

    mulsd %xmm3,%xmm7  ## Heps
    mulsd  %xmm3,%xmm6 ## Geps
    mulsd %xmm3,%xmm7  ## Heps2

    addsd  %xmm6,%xmm5  ## F+Geps
    addsd  %xmm7,%xmm5  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm7  ## 2*Heps2
    addsd  %xmm6,%xmm7  ## 2*Heps2+Geps
    addsd  %xmm5,%xmm7  ## FF = Fp + 2*Heps2 + Geps
    mulsd  %xmm3,%xmm5  ## eps*Fp
    addsd  %xmm4,%xmm5  ## VV
    mulsd  nb310_qq(%rsp),%xmm5     ## VV*qq=vcoul
    mulsd  nb310_qq(%rsp),%xmm7     ## FF*qq=fijC

    ## LJ forces
    mulsd  nb310_six(%rsp),%xmm2
    mulsd  nb310_twelve(%rsp),%xmm0
    subsd  %xmm2,%xmm0
    mulsd  %xmm1,%xmm0 ## (12*vnb12-6*vnb6)*rinv

    ## add potential to vctot 
        addsd  nb310_vctot(%rsp),%xmm5
    movsd %xmm5,nb310_vctot(%rsp)

    mulsd  nb310_tsc(%rsp),%xmm7
    subsd  %xmm7,%xmm0

    mulsd  %xmm1,%xmm0 ## fscal

    ## calculate scalar force by multiplying dx/dy/dz with fscal
        mulsd  %xmm0,%xmm9
        mulsd  %xmm0,%xmm10
        mulsd  %xmm0,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15

        ## the fj's - start by accumulating forces from memory 
    movq nb310_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel310_x86_64_sse2.nb310_updateouterdata: 
        movl  nb310_ii3(%rsp),%ecx
        movq  nb310_faction(%rbp),%rdi
        movq  nb310_fshift(%rbp),%rsi
        movl  nb310_is3(%rsp),%edx

        ## accumulate i forces in xmm13, xmm14, xmm15
        movhlps %xmm13,%xmm3
        movhlps %xmm14,%xmm4
        movhlps %xmm15,%xmm5
        addsd  %xmm3,%xmm13
        addsd  %xmm4,%xmm14
        addsd  %xmm5,%xmm15 ## sum is in low xmm13-xmm15

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## increment fshift force  
        movsd  (%rsi,%rdx,8),%xmm3
        movsd  8(%rsi,%rdx,8),%xmm4
        movsd  16(%rsi,%rdx,8),%xmm5
        subsd  %xmm13,%xmm3
        subsd  %xmm14,%xmm4
        subsd  %xmm15,%xmm5
        movsd  %xmm3,(%rsi,%rdx,8)
        movsd  %xmm4,8(%rsi,%rdx,8)
        movsd  %xmm5,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb310_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb310_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb310_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb310_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb310_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb310_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb310_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310_x86_64_sse2.nb310_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310_n(%rsp)
        jmp _nb_kernel310_x86_64_sse2.nb310_outer
_nb_kernel310_x86_64_sse2.nb310_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310_x86_64_sse2.nb310_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310_x86_64_sse2.nb310_threadloop
_nb_kernel310_x86_64_sse2.nb310_end: 
        movl nb310_nouter(%rsp),%eax
        movl nb310_ninner(%rsp),%ebx
        movq nb310_outeriter(%rbp),%rcx
        movq nb310_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $472,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret





.globl nb_kernel310nf_x86_64_sse2
.globl _nb_kernel310nf_x86_64_sse2
nb_kernel310nf_x86_64_sse2:     
_nb_kernel310nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb310nf_fshift, 16
.set nb310nf_gid, 24
.set nb310nf_pos, 32
.set nb310nf_faction, 40
.set nb310nf_charge, 48
.set nb310nf_p_facel, 56
.set nb310nf_argkrf, 64
.set nb310nf_argcrf, 72
.set nb310nf_Vc, 80
.set nb310nf_type, 88
.set nb310nf_p_ntype, 96
.set nb310nf_vdwparam, 104
.set nb310nf_Vvdw, 112
.set nb310nf_p_tabscale, 120
.set nb310nf_VFtab, 128
.set nb310nf_invsqrta, 136
.set nb310nf_dvda, 144
.set nb310nf_p_gbtabscale, 152
.set nb310nf_GBtab, 160
.set nb310nf_p_nthreads, 168
.set nb310nf_count, 176
.set nb310nf_mtx, 184
.set nb310nf_outeriter, 192
.set nb310nf_inneriter, 200
.set nb310nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb310nf_ix, 0
.set nb310nf_iy, 16
.set nb310nf_iz, 32
.set nb310nf_iq, 48
.set nb310nf_tsc, 64
.set nb310nf_qq, 80
.set nb310nf_c6, 96
.set nb310nf_c12, 112
.set nb310nf_vctot, 128
.set nb310nf_Vvdwtot, 144
.set nb310nf_half, 160
.set nb310nf_three, 176
.set nb310nf_nri, 192
.set nb310nf_iinr, 200
.set nb310nf_jindex, 208
.set nb310nf_jjnr, 216
.set nb310nf_shift, 224
.set nb310nf_shiftvec, 232
.set nb310nf_facel, 240
.set nb310nf_innerjjnr, 248
.set nb310nf_is3, 256
.set nb310nf_ii3, 260
.set nb310nf_ntia, 264
.set nb310nf_innerk, 268
.set nb310nf_n, 272
.set nb310nf_nn1, 276
.set nb310nf_ntype, 280
.set nb310nf_nouter, 284
.set nb310nf_ninner, 288
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $312,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb310nf_nouter(%rsp)
        movl %eax,nb310nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb310nf_nri(%rsp)
        movq %rsi,nb310nf_iinr(%rsp)
        movq %rdx,nb310nf_jindex(%rsp)
        movq %rcx,nb310nf_jjnr(%rsp)
        movq %r8,nb310nf_shift(%rsp)
        movq %r9,nb310nf_shiftvec(%rsp)
        movq nb310nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb310nf_ntype(%rsp)
        movq nb310nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb310nf_facel(%rsp)

        movq nb310nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb310nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb310nf_half(%rsp)
        movl %ebx,nb310nf_half+4(%rsp)
        movsd nb310nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb310nf_half(%rsp)
        movapd %xmm3,nb310nf_three(%rsp)

_nb_kernel310nf_x86_64_sse2.nb310nf_threadloop: 
        movq  nb310nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel310nf_x86_64_sse2.nb310nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel310nf_x86_64_sse2.nb310nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb310nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb310nf_n(%rsp)
        movl %ebx,nb310nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel310nf_x86_64_sse2.nb310nf_outerstart
        jmp _nb_kernel310nf_x86_64_sse2.nb310nf_end

_nb_kernel310nf_x86_64_sse2.nb310nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb310nf_nouter(%rsp),%ebx
        movl %ebx,nb310nf_nouter(%rsp)

_nb_kernel310nf_x86_64_sse2.nb310nf_outer: 
        movq  nb310nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb310nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb310nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb310nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb310nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb310nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb310nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb310nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb310nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb310nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb310nf_ix(%rsp)
        movapd %xmm1,nb310nf_iy(%rsp)
        movapd %xmm2,nb310nf_iz(%rsp)

        movl  %ebx,nb310nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb310nf_vctot(%rsp)
        movapd %xmm4,nb310nf_Vvdwtot(%rsp)

        movq  nb310nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb310nf_pos(%rbp),%rsi
        movq  nb310nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb310nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb310nf_ninner(%rsp),%ecx
        movl  %ecx,nb310nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb310nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel310nf_x86_64_sse2.nb310nf_unroll_loop
        jmp   _nb_kernel310nf_x86_64_sse2.nb310nf_checksingle
_nb_kernel310nf_x86_64_sse2.nb310nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb310nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb310nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb310nf_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb310nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb310nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb310nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb310nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb310nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb310nf_c6(%rsp)
        movapd %xmm6,nb310nf_c12(%rsp)

        movq nb310nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb310nf_ix(%rsp),%xmm4
        movapd nb310nf_iy(%rsp),%xmm5
        movapd nb310nf_iz(%rsp),%xmm6

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
        ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb310nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb310nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb310nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb310nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        mulpd nb310nf_tsc(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movq nb310nf_VFtab(%rbp),%rsi
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
        movapd nb310nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        ## L-J 
        movapd %xmm0,%xmm4
        mulpd  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addpd  nb310nf_vctot(%rsp),%xmm5

        movapd %xmm4,%xmm6
        mulpd  %xmm4,%xmm6

        movapd %xmm5,nb310nf_vctot(%rsp)

        mulpd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb310nf_c6(%rsp),%xmm6
        mulpd  nb310nf_c12(%rsp),%xmm4
        movapd nb310nf_Vvdwtot(%rsp),%xmm7
        addpd  %xmm4,%xmm7
        subpd  %xmm6,%xmm7
        movapd %xmm7,nb310nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb310nf_innerk(%rsp)
        jl    _nb_kernel310nf_x86_64_sse2.nb310nf_checksingle
        jmp   _nb_kernel310nf_x86_64_sse2.nb310nf_unroll_loop
_nb_kernel310nf_x86_64_sse2.nb310nf_checksingle: 
        movl  nb310nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel310nf_x86_64_sse2.nb310nf_dosingle
        jmp    _nb_kernel310nf_x86_64_sse2.nb310nf_updateouterdata
_nb_kernel310nf_x86_64_sse2.nb310nf_dosingle: 
        movq nb310nf_charge(%rbp),%rsi
        movq nb310nf_pos(%rbp),%rdi
        movq  nb310nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        xorpd  %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd nb310nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb310nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movq nb310nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb310nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb310nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb310nf_c6(%rsp)
        movapd %xmm6,nb310nf_c12(%rsp)

        movq nb310nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb310nf_ix(%rsp),%xmm4
        movapd nb310nf_iy(%rsp),%xmm5
        movapd nb310nf_iz(%rsp),%xmm6

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
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb310nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb310nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb310nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb310nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        mulsd nb310nf_tsc(%rsp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movq nb310nf_VFtab(%rbp),%rsi

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
        movapd nb310nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        ## at this point mm5 contains vcoul 

        ## L-J 
        movapd %xmm0,%xmm4
        mulsd  %xmm0,%xmm4      ## xmm4=rinvsq 

        ## increment vcoul - then we can get rid of mm5 
        ## update vctot 
        addsd  nb310nf_vctot(%rsp),%xmm5

        movapd %xmm4,%xmm6
        mulsd  %xmm4,%xmm6

        movlpd %xmm5,nb310nf_vctot(%rsp)

        mulsd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb310nf_c6(%rsp),%xmm6
        mulsd  nb310nf_c12(%rsp),%xmm4
        movapd nb310nf_Vvdwtot(%rsp),%xmm7
        addsd  %xmm4,%xmm7
        subsd  %xmm6,%xmm7
        movlpd %xmm7,nb310nf_Vvdwtot(%rsp)

_nb_kernel310nf_x86_64_sse2.nb310nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb310nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb310nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb310nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb310nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb310nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb310nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb310nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel310nf_x86_64_sse2.nb310nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb310nf_n(%rsp)
        jmp _nb_kernel310nf_x86_64_sse2.nb310nf_outer
_nb_kernel310nf_x86_64_sse2.nb310nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb310nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel310nf_x86_64_sse2.nb310nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel310nf_x86_64_sse2.nb310nf_threadloop
_nb_kernel310nf_x86_64_sse2.nb310nf_end: 
        movl nb310nf_nouter(%rsp),%eax
        movl nb310nf_ninner(%rsp),%ebx
        movq nb310nf_outeriter(%rbp),%rcx
        movq nb310nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $312,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

