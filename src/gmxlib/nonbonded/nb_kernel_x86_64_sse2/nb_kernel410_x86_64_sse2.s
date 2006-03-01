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







.globl nb_kernel410_x86_64_sse2
.globl _nb_kernel410_x86_64_sse2
nb_kernel410_x86_64_sse2:       
_nb_kernel410_x86_64_sse2:      
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
        ## bottom of stack is cache-aligned for sse2 use 
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
.set nb410_ii, 496
.set nb410_is3, 500
.set nb410_ii3, 504
.set nb410_ntia, 508
.set nb410_innerk, 512
.set nb410_n, 516
.set nb410_nn1, 520
.set nb410_ntype, 524
.set nb410_nouter, 528
.set nb410_ninner, 532
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $552,%rsp          ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb410_facel(%rsp)

        movq nb410_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb410_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb410_half(%rsp)
        movl %ebx,nb410_half+4(%rsp)
        movsd nb410_half(%rsp),%xmm1
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
        movapd %xmm1,nb410_half(%rsp)
        movapd %xmm2,nb410_two(%rsp)
        movapd %xmm3,nb410_three(%rsp)
        movapd %xmm4,nb410_six(%rsp)
        movapd %xmm5,nb410_twelve(%rsp)

_nb_kernel410_x86_64_sse2.nb410_threadloop: 
        movq  nb410_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel410_x86_64_sse2.nb410_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410_x86_64_sse2.nb410_spinlock

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
        jg  _nb_kernel410_x86_64_sse2.nb410_outerstart
        jmp _nb_kernel410_x86_64_sse2.nb410_end

_nb_kernel410_x86_64_sse2.nb410_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410_nouter(%rsp),%ebx
        movl %ebx,nb410_nouter(%rsp)

_nb_kernel410_x86_64_sse2.nb410_outer: 
        movq  nb410_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb410_is3(%rsp)      ## store is3 

        movq  nb410_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb410_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb410_ii(%rsp)

        movq  nb410_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb410_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb410_invsqrta(%rbp),%rdx         ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movq  nb410_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb410_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb410_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb410_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb410_iq(%rsp)
        movapd %xmm4,nb410_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb410_ix(%rsp)
        movapd %xmm1,nb410_iy(%rsp)
        movapd %xmm2,nb410_iz(%rsp)

        movl  %ebx,nb410_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm13,%xmm13
        movapd %xmm13,%xmm12
        movapd %xmm13,nb410_Vvdwtot(%rsp)
        movapd %xmm13,nb410_dvdasum(%rsp)
        movapd %xmm13,%xmm14
        movapd %xmm13,%xmm15

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
        subl  $2,%edx
        addl  nb410_ninner(%rsp),%ecx
        movl  %ecx,nb410_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb410_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel410_x86_64_sse2.nb410_unroll_loop
        jmp   _nb_kernel410_x86_64_sse2.nb410_checksingle
_nb_kernel410_x86_64_sse2.nb410_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb410_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r14d
        movl  4(%rdx),%r15d
        addq $8,nb410_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        movq nb410_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%r14,%r14,2),%r10     ## replace jnr with j3 
        lea  (%r15,%r15,2),%r11

        ## move two coordinates to xmm4-xmm6    
        movlpd (%rsi,%r10,8),%xmm4
        movlpd 8(%rsi,%r10,8),%xmm5
        movlpd 16(%rsi,%r10,8),%xmm6
        movhpd (%rsi,%r11,8),%xmm4
        movhpd 8(%rsi,%r11,8),%xmm5
        movhpd 16(%rsi,%r11,8),%xmm6

        ## calc dr 
        subpd nb410_ix(%rsp),%xmm4
        subpd nb410_iy(%rsp),%xmm5
        subpd nb410_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb410_dx(%rsp)
        movapd %xmm5,nb410_dy(%rsp)
        movapd %xmm6,nb410_dz(%rsp)

        ## load isaj
        movq nb410_invsqrta(%rbp),%rsi

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 

        movlpd (%rsi,%r14,8),%xmm3
        movhpd (%rsi,%r15,8),%xmm3

        movq nb410_type(%rbp),%rdi
        movl (%rdi,%r14,4),%r8d
        movl (%rdi,%r15,4),%r9d

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        mulpd  nb410_isai(%rsp),%xmm3
        movapd %xmm3,nb410_isaprod(%rsp)

        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb410_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm3,%xmm6
        mulpd nb410_gbtsc(%rsp),%xmm6
        movapd %xmm6,nb410_gbscale(%rsp)

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb410_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

    mulpd  nb410_iq(%rsp),%xmm3
        movq nb410_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%r14,8),%xmm6
        movhpd (%rsi,%r15,8),%xmm6
        mulpd  %xmm3,%xmm6
        movapd %xmm6,nb410_qq(%rsp)

        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410_r(%rsp)
        mulpd nb410_gbscale(%rsp),%xmm4
        movl nb410_ntia(%rsp),%edi

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        shll %r8d
        shll %r9d
        addl %edi,%r8d
        addl %edi,%r9d

        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 
        movq nb410_vdwparam(%rbp),%rdi

        pslld $2,%mm6           ## idx *= 4 

        movq nb410_GBtab(%rbp),%rsi
        movd %mm6,%r12d
        psrlq $32,%mm6
        movd %mm6,%r13d         ## indices in r12/r13

        movlpd (%rdi,%r8,8),%xmm6
        movlpd 8(%rdi,%r8,8),%xmm7

    movapd %xmm0,%xmm9 ## rinv
    mulpd  %xmm9,%xmm9 ## rinvsq
    movapd %xmm9,%xmm10 ## rinvsq
    mulpd  %xmm10,%xmm10 ## rinv4
    mulpd  %xmm9,%xmm10 ## rinv6
    movapd %xmm10,%xmm11
    mulpd  %xmm11,%xmm11 ## rinv12


        movhpd (%rdi,%r9,8),%xmm6
        movhpd 8(%rdi,%r9,8),%xmm7

    ## load table data
        movapd (%rsi,%r12,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%r13,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

    mulpd  %xmm6,%xmm10  ## vvdw6=c6*rinv6
        mulpd  %xmm7,%xmm11  ## vvdw12=c12*rinv12     

        movapd %xmm11,%xmm9
        subpd  %xmm10,%xmm11    ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addpd  nb410_Vvdwtot(%rsp),%xmm11
    movapd %xmm11,nb410_Vvdwtot(%rsp)

        movapd 16(%rsi,%r12,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%r13,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            

        mulpd  %xmm1,%xmm7      ## xmm7=Heps
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm1,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb410_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb410_qq(%rsp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 

   ## LJ forces
    mulpd  nb410_six(%rsp),%xmm10
    mulpd  nb410_twelve(%rsp),%xmm9
    subpd  %xmm10,%xmm9
    mulpd  %xmm0,%xmm9 ## (12*vnb12-6*vnb6)*rinv

        movq nb410_dvda(%rbp),%rsi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulpd nb410_gbscale(%rsp),%xmm3
        movapd %xmm3,%xmm6
        mulpd  nb410_r(%rsp),%xmm6
        addpd  %xmm5,%xmm6

    ## update vctot 
        addpd  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subpd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        movq nb410_faction(%rbp),%rdi
        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%r10,8),%xmm2
        movlpd 8(%rdi,%r10,8),%xmm4
        movlpd 16(%rdi,%r10,8),%xmm5

        ## update dvdasum
        addpd  nb410_dvdasum(%rsp),%xmm7
        movapd %xmm7,nb410_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%rsi,%r14,8),%xmm6
        addsd  (%rsi,%r15,8),%xmm7
        movsd  %xmm6,(%rsi,%r14,8)
        movsd  %xmm7,(%rsi,%r15,8)

        movhpd (%rdi,%r11,8),%xmm2
        movhpd 8(%rdi,%r11,8),%xmm4
        movhpd 16(%rdi,%r11,8),%xmm5

    subpd  %xmm3,%xmm9
    mulpd  %xmm0,%xmm9 ## fscal

    movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

    mulpd   nb410_dx(%rsp),%xmm9
    mulpd   nb410_dy(%rsp),%xmm10
    mulpd   nb410_dz(%rsp),%xmm11

        addpd %xmm9,%xmm2
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5

        movlpd %xmm2,(%rdi,%r10,8)
        movlpd %xmm4,8(%rdi,%r10,8)
        movlpd %xmm5,16(%rdi,%r10,8)

        ## accumulate i forces
    addpd %xmm9,%xmm13
    addpd %xmm10,%xmm14
    addpd %xmm11,%xmm15

        movhpd %xmm2,(%rdi,%r11,8)
        movhpd %xmm4,8(%rdi,%r11,8)
        movhpd %xmm5,16(%rdi,%r11,8)

        ## should we do one more iteration? 
        subl $2,nb410_innerk(%rsp)
        jl    _nb_kernel410_x86_64_sse2.nb410_checksingle
        jmp   _nb_kernel410_x86_64_sse2.nb410_unroll_loop
_nb_kernel410_x86_64_sse2.nb410_checksingle: 
        movl  nb410_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410_x86_64_sse2.nb410_dosingle
        jmp    _nb_kernel410_x86_64_sse2.nb410_updateouterdata
_nb_kernel410_x86_64_sse2.nb410_dosingle: 
        movq nb410_charge(%rbp),%rsi
        movq nb410_invsqrta(%rbp),%rdx
        movq nb410_pos(%rbp),%rdi
        movq  nb410_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        ## load isaj
        movq nb410_invsqrta(%rbp),%rsi
        movsd (%rsi,%rax,8),%xmm2
        mulsd  nb410_isai(%rsp),%xmm2
        movapd %xmm2,nb410_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulsd nb410_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb410_gbscale(%rsp)

    mulsd nb410_iq(%rsp),%xmm2
        movq nb410_charge(%rbp),%rsi     ## base of charge[] 
        movsd (%rsi,%rax,8),%xmm3
        mulsd  %xmm2,%xmm3
        movapd %xmm3,nb410_qq(%rsp)

        movq nb410_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb410_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb410_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb410_c6(%rsp)
        movapd %xmm6,nb410_c12(%rsp)

        movq nb410_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%r10     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6    
        movsd (%rsi,%r10,8),%xmm4
        movsd 8(%rsi,%r10,8),%xmm5
        movsd 16(%rsi,%r10,8),%xmm6

        ## calc dr 
        subsd nb410_ix(%rsp),%xmm4
        subsd nb410_iy(%rsp),%xmm5
        subsd nb410_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb410_dx(%rsp)
        movapd %xmm5,nb410_dy(%rsp)
        movapd %xmm6,nb410_dz(%rsp)

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
        movapd nb410_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb410_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410_r(%rsp)
        mulsd nb410_gbscale(%rsp),%xmm4

        cvttsd2si %xmm4,%r12d   ## mm6 = lu idx 
        cvtsi2sd %r12d,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%r12d           ## idx *= 4 

        movq nb410_GBtab(%rbp),%rsi

    movapd %xmm0,%xmm9 ## rinv
    mulsd  %xmm9,%xmm9 ## rinvsq
    movapd %xmm9,%xmm10 ## rinvsq
    mulsd  %xmm10,%xmm10 ## rinv4
    mulsd  %xmm9,%xmm10 ## rinv6
    movapd %xmm10,%xmm11
    mulsd  %xmm11,%xmm11 ## rinv12

    ## load table data
        movapd (%rsi,%r12,8),%xmm4      ## Y1 F1        
    movhlps %xmm4,%xmm5

    mulsd  nb410_c6(%rsp),%xmm10      ## vvdw6=c6*rinv6
        mulsd  nb410_c12(%rsp),%xmm11     ## vvdw12=c12*rinv12     

        movapd %xmm11,%xmm9
        subsd  %xmm10,%xmm11    ## Vvdw=Vvdw12-Vvdw6

    ## add potential to vvdwtot 
        addsd  nb410_Vvdwtot(%rsp),%xmm11
    movsd %xmm11,nb410_Vvdwtot(%rsp)

        movapd 16(%rsi,%r12,8),%xmm6    ## G1 H1        
    movhlps %xmm6,%xmm7
        ## coulomb table ready, in xmm4-xmm7            

        mulsd  %xmm1,%xmm7      ## xmm7=Heps
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm1,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb410_two(%rsp),%xmm7    ## two*Heps2 
        movapd nb410_qq(%rsp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 

   ## LJ forces
    mulsd  nb410_six(%rsp),%xmm10
    mulsd  nb410_twelve(%rsp),%xmm9
    subsd  %xmm10,%xmm9
    mulsd  %xmm0,%xmm9 ## (12*vnb12-6*vnb6)*rinv

        movq nb410_dvda(%rbp),%rsi

        ## Calculate dVda
        xorpd %xmm7,%xmm7
        mulsd nb410_gbscale(%rsp),%xmm3
        movapd %xmm3,%xmm6
        mulsd  nb410_r(%rsp),%xmm6
        addsd  %xmm5,%xmm6

    ## update vctot 
        addsd  %xmm5,%xmm12

        ## xmm6=(vcoul+fijC*r)
        subsd  %xmm6,%xmm7
        movapd %xmm7,%xmm6

        ## update dvdasum
        addsd  nb410_dvdasum(%rsp),%xmm7
        movsd %xmm7,nb410_dvdasum(%rsp)

        ## update j atoms dvdaj
        movhlps %xmm6,%xmm7
        addsd  (%rsi,%rax,8),%xmm6
        addsd  (%rsi,%rbx,8),%xmm7
        movsd  %xmm6,(%rsi,%rax,8)
        movsd  %xmm7,(%rsi,%rbx,8)

    subsd  %xmm3,%xmm9
    mulsd  %xmm0,%xmm9 ## fscal

    movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

    mulsd   nb410_dx(%rsp),%xmm9
    mulsd   nb410_dy(%rsp),%xmm10
    mulsd   nb410_dz(%rsp),%xmm11

        ## accumulate i forces
    addsd %xmm9,%xmm13
    addsd %xmm10,%xmm14
    addsd %xmm11,%xmm15

        movq nb410_faction(%rbp),%rdi
        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%r10,8),%xmm9
        addsd 8(%rdi,%r10,8),%xmm10
        addsd 16(%rdi,%r10,8),%xmm11
        movsd %xmm9,(%rdi,%r10,8)
        movsd %xmm10,8(%rdi,%r10,8)
        movsd %xmm11,16(%rdi,%r10,8)

_nb_kernel410_x86_64_sse2.nb410_updateouterdata: 
        movl  nb410_ii3(%rsp),%ecx
        movq  nb410_faction(%rbp),%rdi
        movq  nb410_fshift(%rbp),%rsi
        movl  nb410_is3(%rsp),%edx

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
        movl nb410_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb410_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 has the sum now 

        ## add earlier value from mem 
        movq  nb410_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb410_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb410_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate dVda and update it 
        movapd nb410_dvdasum(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        movl nb410_ii(%rsp),%edx
        movq nb410_dvda(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb410_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410_x86_64_sse2.nb410_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410_n(%rsp)
        jmp _nb_kernel410_x86_64_sse2.nb410_outer
_nb_kernel410_x86_64_sse2.nb410_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410_x86_64_sse2.nb410_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410_x86_64_sse2.nb410_threadloop
_nb_kernel410_x86_64_sse2.nb410_end: 
        movl nb410_nouter(%rsp),%eax
        movl nb410_ninner(%rsp),%ebx
        movq nb410_outeriter(%rbp),%rcx
        movq nb410_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $552,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret









.globl nb_kernel410nf_x86_64_sse2
.globl _nb_kernel410nf_x86_64_sse2
nb_kernel410nf_x86_64_sse2:     
_nb_kernel410nf_x86_64_sse2:    
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
        ## bottom of stack is cache-aligned for sse2 use 
.set nb410nf_ix, 0
.set nb410nf_iy, 16
.set nb410nf_iz, 32
.set nb410nf_iq, 48
.set nb410nf_two, 64
.set nb410nf_gbtsc, 80
.set nb410nf_qq, 96
.set nb410nf_c6, 112
.set nb410nf_c12, 128
.set nb410nf_vctot, 144
.set nb410nf_Vvdwtot, 160
.set nb410nf_half, 176
.set nb410nf_three, 192
.set nb410nf_r, 208
.set nb410nf_isai, 224
.set nb410nf_isaprod, 240
.set nb410nf_gbscale, 256
.set nb410nf_nri, 272
.set nb410nf_iinr, 280
.set nb410nf_jindex, 288
.set nb410nf_jjnr, 296
.set nb410nf_shift, 304
.set nb410nf_shiftvec, 312
.set nb410nf_facel, 320
.set nb410nf_innerjjnr, 328
.set nb410nf_ii, 336
.set nb410nf_is3, 340
.set nb410nf_ii3, 344
.set nb410nf_ntia, 348
.set nb410nf_innerk, 352
.set nb410nf_n, 356
.set nb410nf_nn1, 360
.set nb410nf_ntype, 364
.set nb410nf_nouter, 368
.set nb410nf_ninner, 372
        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $392,%rsp          ## local variable stack space (n*16+8)

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
        movsd (%rsi),%xmm0
        movsd %xmm0,nb410nf_facel(%rsp)

        movq nb410nf_p_gbtabscale(%rbp),%rbx
        movsd (%rbx),%xmm4
        shufpd $0,%xmm4,%xmm4
        movapd %xmm4,nb410nf_gbtsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb410nf_half(%rsp)
        movl %ebx,nb410nf_half+4(%rsp)
        movsd nb410nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb410nf_half(%rsp)
        movapd %xmm2,nb410nf_two(%rsp)
        movapd %xmm3,nb410nf_three(%rsp)

_nb_kernel410nf_x86_64_sse2.nb410nf_threadloop: 
        movq  nb410nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel410nf_x86_64_sse2.nb410nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel410nf_x86_64_sse2.nb410nf_spinlock

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
        jg  _nb_kernel410nf_x86_64_sse2.nb410nf_outerstart
        jmp _nb_kernel410nf_x86_64_sse2.nb410nf_end

_nb_kernel410nf_x86_64_sse2.nb410nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb410nf_nouter(%rsp),%ebx
        movl %ebx,nb410nf_nouter(%rsp)

_nb_kernel410nf_x86_64_sse2.nb410nf_outer: 
        movq  nb410nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb410nf_is3(%rsp)            ## store is3 

        movq  nb410nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb410nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 
        movl  %ebx,nb410nf_ii(%rsp)

        movq  nb410nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb410nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb410nf_invsqrta(%rbp),%rdx       ## load invsqrta[ii]
        movsd (%rdx,%rbx,8),%xmm4
        shufpd $0,%xmm4,%xmm4

        movq  nb410nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb410nf_ntype(%rsp),%edx
        shll  %edx
    movl  %edx,nb410nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb410nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb410nf_iq(%rsp)
        movapd %xmm4,nb410nf_isai(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb410nf_ix(%rsp)
        movapd %xmm1,nb410nf_iy(%rsp)
        movapd %xmm2,nb410nf_iz(%rsp)

        movl  %ebx,nb410nf_ii3(%rsp)

        ## clear vctot and Vvdwtot
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb410nf_vctot(%rsp)
        movapd %xmm4,nb410nf_Vvdwtot(%rsp)

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
        subl  $2,%edx
        addl  nb410nf_ninner(%rsp),%ecx
        movl  %ecx,nb410nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb410nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel410nf_x86_64_sse2.nb410nf_unroll_loop
        jmp   _nb_kernel410nf_x86_64_sse2.nb410nf_checksingle
_nb_kernel410nf_x86_64_sse2.nb410nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb410nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb410nf_innerjjnr(%rsp)             ## advance pointer (unrolled 2) 

        ## load isaj
        movq nb410nf_invsqrta(%rbp),%rsi
        movlpd (%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm2
        mulpd  nb410nf_isai(%rsp),%xmm2
        movapd %xmm2,nb410nf_isaprod(%rsp)
        movapd %xmm2,%xmm1
        mulpd nb410nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb410nf_gbscale(%rsp)

        movq nb410nf_charge(%rbp),%rsi     ## base of charge[] 
        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        mulpd nb410nf_iq(%rsp),%xmm2
        mulpd  %xmm2,%xmm3
        movapd %xmm3,nb410nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb410nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb410nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb410nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb410nf_c6(%rsp)
        movapd %xmm6,nb410nf_c12(%rsp)

        movq nb410nf_pos(%rbp),%rsi        ## base of pos[] 

        movd  %eax,%mm2
        movd  %ebx,%mm3
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
        movapd nb410nf_ix(%rsp),%xmm4
        movapd nb410nf_iy(%rsp),%xmm5
        movapd nb410nf_iz(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square dr 
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
        movapd nb410nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb410nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        mulpd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410nf_r(%rsp)
        mulpd nb410nf_gbscale(%rsp),%xmm4

        cvttpd2pi %xmm4,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm5
        subpd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 

        movd %eax,%mm0
        movd %ebx,%mm1

        movq nb410nf_GBtab(%rbp),%rsi
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
        movapd nb410nf_qq(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb410nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb410nf_vctot(%rsp)

        ## L-J 
        movapd %xmm0,%xmm4
        mulpd  %xmm0,%xmm4      ## xmm4=rinvsq 

        movapd %xmm4,%xmm6
        mulpd  %xmm4,%xmm6

        mulpd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulpd  nb410nf_c6(%rsp),%xmm6
        mulpd  nb410nf_c12(%rsp),%xmm4
        movapd nb410nf_Vvdwtot(%rsp),%xmm7
        addpd  %xmm4,%xmm7
        subpd  %xmm6,%xmm7
        movapd %xmm7,nb410nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb410nf_innerk(%rsp)
        jl    _nb_kernel410nf_x86_64_sse2.nb410nf_checksingle
        jmp   _nb_kernel410nf_x86_64_sse2.nb410nf_unroll_loop
_nb_kernel410nf_x86_64_sse2.nb410nf_checksingle: 
        movl  nb410nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel410nf_x86_64_sse2.nb410nf_dosingle
        jmp    _nb_kernel410nf_x86_64_sse2.nb410nf_updateouterdata
_nb_kernel410nf_x86_64_sse2.nb410nf_dosingle: 
        movq nb410nf_charge(%rbp),%rsi
        movq nb410nf_invsqrta(%rbp),%rdx
        movq nb410nf_pos(%rbp),%rdi
        movq  nb410nf_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        xorpd  %xmm6,%xmm6
        movapd %xmm6,%xmm7
        movsd  (%rdx,%rax,8),%xmm7
        movlpd (%rsi,%rax,8),%xmm6      ## xmm6(0) has the charge
        mulsd  nb410nf_isai(%rsp),%xmm7
        movapd %xmm7,nb410nf_isaprod(%rsp)
        movapd %xmm7,%xmm1
        mulpd nb410nf_gbtsc(%rsp),%xmm1
        movapd %xmm1,nb410nf_gbscale(%rsp)

        mulsd  nb410nf_iq(%rsp),%xmm7
        mulsd  %xmm7,%xmm6
        movapd %xmm6,nb410nf_qq(%rsp)

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movq nb410nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb410nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb410nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb410nf_c6(%rsp)
        movapd %xmm6,nb410nf_c12(%rsp)

        movq nb410nf_pos(%rbp),%rsi        ## base of pos[]

        movd  %eax,%mm2
        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb410nf_ix(%rsp),%xmm4
        movapd nb410nf_iy(%rsp),%xmm5
        movapd nb410nf_iz(%rsp),%xmm6

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
        movapd nb410nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb410nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb410nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb410nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        mulsd %xmm0,%xmm4       ## xmm4=r 
        movapd %xmm4,nb410nf_r(%rsp)
        mulsd nb410nf_gbscale(%rsp),%xmm4

        movd %eax,%mm0
        cvttsd2si %xmm4,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm5
        subsd %xmm5,%xmm4
        movapd %xmm4,%xmm1      ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 

        movq nb410nf_GBtab(%rbp),%rsi

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
        movapd nb410nf_qq(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb410nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb410nf_vctot(%rsp)

        ## L-J 
        movapd %xmm0,%xmm4
        mulsd  %xmm0,%xmm4      ## xmm4=rinvsq 


        movapd %xmm4,%xmm6
        mulsd  %xmm4,%xmm6

        mulsd  %xmm4,%xmm6      ## xmm6=rinvsix 
        movapd %xmm6,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulsd  nb410nf_c6(%rsp),%xmm6
        mulsd  nb410nf_c12(%rsp),%xmm4
        movapd nb410nf_Vvdwtot(%rsp),%xmm7
        addsd  %xmm4,%xmm7
        subsd  %xmm6,%xmm7
        movlpd %xmm7,nb410nf_Vvdwtot(%rsp)

_nb_kernel410nf_x86_64_sse2.nb410nf_updateouterdata: 
        movl  nb410nf_ii3(%rsp),%ecx
        movl  nb410nf_is3(%rsp),%edx

        ## get n from stack
        movl nb410nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb410nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb410nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb410nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb410nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb410nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb410nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jecxz _nb_kernel410nf_x86_64_sse2.nb410nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb410nf_n(%rsp)
        jmp _nb_kernel410nf_x86_64_sse2.nb410nf_outer
_nb_kernel410nf_x86_64_sse2.nb410nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb410nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jecxz _nb_kernel410nf_x86_64_sse2.nb410nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel410nf_x86_64_sse2.nb410nf_threadloop
_nb_kernel410nf_x86_64_sse2.nb410nf_end: 
        movl nb410nf_nouter(%rsp),%eax
        movl nb410nf_ninner(%rsp),%ebx
        movq nb410nf_outeriter(%rbp),%rcx
        movq nb410nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $392,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret



