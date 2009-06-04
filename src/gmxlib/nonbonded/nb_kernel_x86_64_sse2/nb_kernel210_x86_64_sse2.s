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





.globl nb_kernel210_x86_64_sse2
.globl _nb_kernel210_x86_64_sse2
nb_kernel210_x86_64_sse2:       
_nb_kernel210_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb210_fshift, 16
.set nb210_gid, 24
.set nb210_pos, 32
.set nb210_faction, 40
.set nb210_charge, 48
.set nb210_p_facel, 56
.set nb210_argkrf, 64
.set nb210_argcrf, 72
.set nb210_Vc, 80
.set nb210_type, 88
.set nb210_p_ntype, 96
.set nb210_vdwparam, 104
.set nb210_Vvdw, 112
.set nb210_p_tabscale, 120
.set nb210_VFtab, 128
.set nb210_invsqrta, 136
.set nb210_dvda, 144
.set nb210_p_gbtabscale, 152
.set nb210_GBtab, 160
.set nb210_p_nthreads, 168
.set nb210_count, 176
.set nb210_mtx, 184
.set nb210_outeriter, 192
.set nb210_inneriter, 200
.set nb210_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb210_ix, 0
.set nb210_iy, 16
.set nb210_iz, 32
.set nb210_iq, 48
.set nb210_dx, 64
.set nb210_dy, 80
.set nb210_dz, 96
.set nb210_c6, 112
.set nb210_c12, 128
.set nb210_six, 144
.set nb210_twelve, 160
.set nb210_vctot, 176
.set nb210_Vvdwtot, 192
.set nb210_fix, 208
.set nb210_fiy, 224
.set nb210_fiz, 240
.set nb210_half, 256
.set nb210_three, 272
.set nb210_two, 288
.set nb210_krf, 304
.set nb210_crf, 320
.set nb210_nri, 336
.set nb210_iinr, 344
.set nb210_jindex, 352
.set nb210_jjnr, 360
.set nb210_shift, 368
.set nb210_shiftvec, 376
.set nb210_facel, 384
.set nb210_innerjjnr, 392
.set nb210_is3, 400
.set nb210_ii3, 404
.set nb210_ntia, 408
.set nb210_innerk, 412
.set nb210_n, 416
.set nb210_nn1, 420
.set nb210_ntype, 424
.set nb210_nouter, 428
.set nb210_ninner, 432
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $456,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb210_nouter(%rsp)
        movl %eax,nb210_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb210_nri(%rsp)
        movq %rsi,nb210_iinr(%rsp)
        movq %rdx,nb210_jindex(%rsp)
        movq %rcx,nb210_jjnr(%rsp)
        movq %r8,nb210_shift(%rsp)
        movq %r9,nb210_shiftvec(%rsp)
        movq nb210_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb210_ntype(%rsp)
        movq nb210_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb210_facel(%rsp)

        movq nb210_argkrf(%rbp),%rsi
        movq nb210_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb210_krf(%rsp)
        movapd %xmm2,nb210_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb210_half(%rsp)
        movl %ebx,nb210_half+4(%rsp)
        movsd nb210_half(%rsp),%xmm1
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
        movapd %xmm1,nb210_half(%rsp)
        movapd %xmm2,nb210_two(%rsp)
        movapd %xmm3,nb210_three(%rsp)
        movapd %xmm4,nb210_six(%rsp)
        movapd %xmm5,nb210_twelve(%rsp)

_nb_kernel210_x86_64_sse2.nb210_threadloop: 
        movq  nb210_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel210_x86_64_sse2.nb210_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210_x86_64_sse2.nb210_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210_n(%rsp)
        movl %ebx,nb210_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210_x86_64_sse2.nb210_outerstart
        jmp _nb_kernel210_x86_64_sse2.nb210_end

_nb_kernel210_x86_64_sse2.nb210_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210_nouter(%rsp),%ebx
        movl %ebx,nb210_nouter(%rsp)

_nb_kernel210_x86_64_sse2.nb210_outer: 
        movq  nb210_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb210_is3(%rsp)      ## store is3 

        movq  nb210_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb210_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb210_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb210_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb210_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb210_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb210_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb210_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb210_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb210_ix(%rsp)
        movapd %xmm1,nb210_iy(%rsp)
        movapd %xmm2,nb210_iz(%rsp)

        movl  %ebx,nb210_ii3(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm12,%xmm12
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15
        movapd %xmm12,nb210_vctot(%rsp)
        movapd %xmm12,nb210_Vvdwtot(%rsp)

        movq  nb210_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb210_pos(%rbp),%rsi
        movq  nb210_faction(%rbp),%rdi
        movq  nb210_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb210_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb210_ninner(%rsp),%ecx
        movl  %ecx,nb210_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb210_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel210_x86_64_sse2.nb210_unroll_loop
        jmp   _nb_kernel210_x86_64_sse2.nb210_checksingle
_nb_kernel210_x86_64_sse2.nb210_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb210_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r10d
        movl  4(%rdx),%r11d
        addq $8,nb210_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 


        movq nb210_pos(%rbp),%rsi        ## base of pos[] 

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
        subpd nb210_ix(%rsp),%xmm4
        subpd nb210_iy(%rsp),%xmm5
        subpd nb210_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb210_charge(%rbp),%rsi     ## base of charge[] 
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 
        movq nb210_type(%rbp),%rdi

        movlpd (%rsi,%r10,8),%xmm3
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        movhpd (%rsi,%r11,8),%xmm3
        movl (%rdi,%r10,4),%r8d
        movl (%rdi,%r11,4),%r9d

        movapd nb210_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb210_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        shll %r8d
        shll %r9d
        movl nb210_ntia(%rsp),%edi
        movapd nb210_half(%rsp),%xmm0
        movq nb210_vdwparam(%rbp),%rsi
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 
        mulpd nb210_iq(%rsp),%xmm3              ## qq 
        addl %edi,%r8d
        addl %edi,%r9d


        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb210_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210_half(%rsp),%xmm0

    movlpd  (%rsi,%r8,8),%xmm4
    movlpd  8(%rsi,%r8,8),%xmm6
    movhpd  (%rsi,%r9,8),%xmm4
    movhpd  8(%rsi,%r9,8),%xmm6
    movapd %xmm4,nb210_c6(%rsp)
    movapd %xmm6,nb210_c12(%rsp)

        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb210_crf(%rsp),%xmm6
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addpd  %xmm7,%xmm7
        mulpd  nb210_c6(%rsp),%xmm1
        mulpd  nb210_c12(%rsp),%xmm2

        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb210_Vvdwtot(%rsp),%xmm5
        mulpd  nb210_six(%rsp),%xmm1
        mulpd  nb210_twelve(%rsp),%xmm2
        subpd  %xmm1,%xmm2
        subpd  %xmm7,%xmm0
        mulpd  %xmm0,%xmm3
        addpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm4      ## xmm4=total fscal 
        addpd  nb210_vctot(%rsp),%xmm6

        ## the fj's - start by accumulating forces from memory 
        movq   nb210_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

        movapd %xmm6,nb210_vctot(%rsp)
        movapd %xmm5,nb210_Vvdwtot(%rsp)

        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        addpd %xmm9,%xmm0
        addpd %xmm10,%xmm1
        addpd %xmm11,%xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb210_innerk(%rsp)
        jl    _nb_kernel210_x86_64_sse2.nb210_checksingle
        jmp   _nb_kernel210_x86_64_sse2.nb210_unroll_loop

_nb_kernel210_x86_64_sse2.nb210_checksingle:    
        movl  nb210_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210_x86_64_sse2.nb210_dosingle
        jmp    _nb_kernel210_x86_64_sse2.nb210_updateouterdata
_nb_kernel210_x86_64_sse2.nb210_dosingle: 
        movq nb210_charge(%rbp),%rsi
        movq nb210_pos(%rbp),%rdi
        movq  nb210_innerjjnr(%rsp),%rcx

        movl  (%rcx),%eax

        movq nb210_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        mulsd nb210_iq(%rsp),%xmm3              ## qq 

        movq nb210_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb210_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb210_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm4
        movsd 8(%rsi,%r8,8),%xmm6
        movapd %xmm4,nb210_c6(%rsp)
        movapd %xmm6,nb210_c12(%rsp)

        movq nb210_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb210_ix(%rsp),%xmm4
        subsd nb210_iy(%rsp),%xmm5
        subsd nb210_iz(%rsp),%xmm6

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
        ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm5,%xmm2    ## lu in low xmm2 

        movapd nb210_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb210_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb210_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb210_crf(%rsp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addsd  %xmm7,%xmm7
        mulsd  nb210_c6(%rsp),%xmm1
        mulsd  nb210_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb210_Vvdwtot(%rsp),%xmm5
        mulsd  nb210_six(%rsp),%xmm1
        mulsd  nb210_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        subsd  %xmm7,%xmm0
        mulsd  %xmm0,%xmm3
        addsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm4      ## xmm4=total fscal 
        addsd  nb210_vctot(%rsp),%xmm6

        movsd %xmm6,nb210_vctot(%rsp)
        movsd %xmm5,nb210_Vvdwtot(%rsp)

        movq   nb210_faction(%rbp),%rdi
        mulsd  %xmm4,%xmm9
        mulsd  %xmm4,%xmm10
        mulsd  %xmm4,%xmm11

        ## now update f_i 
        addsd  %xmm9,%xmm13
        addsd  %xmm10,%xmm14
        addsd  %xmm11,%xmm15
        ## the fj's - start by accumulating forces from memory 
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)

_nb_kernel210_x86_64_sse2.nb210_updateouterdata: 
        movl  nb210_ii3(%rsp),%ecx
        movq  nb210_faction(%rbp),%rdi
        movq  nb210_fshift(%rbp),%rsi
        movl  nb210_is3(%rsp),%edx

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
        movl nb210_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb210_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb210_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb210_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb210_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb210_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb210_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210_x86_64_sse2.nb210_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210_n(%rsp)
        jmp _nb_kernel210_x86_64_sse2.nb210_outer
_nb_kernel210_x86_64_sse2.nb210_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210_x86_64_sse2.nb210_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210_x86_64_sse2.nb210_threadloop
_nb_kernel210_x86_64_sse2.nb210_end: 
        movl nb210_nouter(%rsp),%eax
        movl nb210_ninner(%rsp),%ebx
        movq nb210_outeriter(%rbp),%rcx
        movq nb210_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $456,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel210nf_x86_64_sse2
.globl _nb_kernel210nf_x86_64_sse2
nb_kernel210nf_x86_64_sse2:     
_nb_kernel210nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb210nf_fshift, 16
.set nb210nf_gid, 24
.set nb210nf_pos, 32
.set nb210nf_faction, 40
.set nb210nf_charge, 48
.set nb210nf_p_facel, 56
.set nb210nf_argkrf, 64
.set nb210nf_argcrf, 72
.set nb210nf_Vc, 80
.set nb210nf_type, 88
.set nb210nf_p_ntype, 96
.set nb210nf_vdwparam, 104
.set nb210nf_Vvdw, 112
.set nb210nf_p_tabscale, 120
.set nb210nf_VFtab, 128
.set nb210nf_invsqrta, 136
.set nb210nf_dvda, 144
.set nb210nf_p_gbtabscale, 152
.set nb210nf_GBtab, 160
.set nb210nf_p_nthreads, 168
.set nb210nf_count, 176
.set nb210nf_mtx, 184
.set nb210nf_outeriter, 192
.set nb210nf_inneriter, 200
.set nb210nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb210nf_ix, 0
.set nb210nf_iy, 16
.set nb210nf_iz, 32
.set nb210nf_iq, 48
.set nb210nf_c6, 64
.set nb210nf_c12, 80
.set nb210nf_vctot, 96
.set nb210nf_Vvdwtot, 112
.set nb210nf_half, 128
.set nb210nf_three, 144
.set nb210nf_krf, 160
.set nb210nf_crf, 176
.set nb210nf_nri, 192
.set nb210nf_iinr, 200
.set nb210nf_jindex, 208
.set nb210nf_jjnr, 216
.set nb210nf_shift, 224
.set nb210nf_shiftvec, 232
.set nb210nf_facel, 240
.set nb210nf_innerjjnr, 248
.set nb210nf_is3, 256
.set nb210nf_ii3, 260
.set nb210nf_ntia, 264
.set nb210nf_innerk, 268
.set nb210nf_n, 272
.set nb210nf_nn1, 276
.set nb210nf_ntype, 280
.set nb210nf_nouter, 284
.set nb210nf_ninner, 288
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
        movl %eax,nb210nf_nouter(%rsp)
        movl %eax,nb210nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb210nf_nri(%rsp)
        movq %rsi,nb210nf_iinr(%rsp)
        movq %rdx,nb210nf_jindex(%rsp)
        movq %rcx,nb210nf_jjnr(%rsp)
        movq %r8,nb210nf_shift(%rsp)
        movq %r9,nb210nf_shiftvec(%rsp)
        movq nb210nf_p_ntype(%rbp),%rdi
        movl (%rdi),%edi
        movl %edi,nb210nf_ntype(%rsp)
        movq nb210nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb210nf_facel(%rsp)

        movq nb210nf_argkrf(%rbp),%rsi
        movq nb210nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb210nf_krf(%rsp)
        movapd %xmm2,nb210nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb210nf_half(%rsp)
        movl %ebx,nb210nf_half+4(%rsp)
        movsd nb210nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb210nf_half(%rsp)
        movapd %xmm3,nb210nf_three(%rsp)

_nb_kernel210nf_x86_64_sse2.nb210nf_threadloop: 
        movq  nb210nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel210nf_x86_64_sse2.nb210nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel210nf_x86_64_sse2.nb210nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb210nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb210nf_n(%rsp)
        movl %ebx,nb210nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel210nf_x86_64_sse2.nb210nf_outerstart
        jmp _nb_kernel210nf_x86_64_sse2.nb210nf_end

_nb_kernel210nf_x86_64_sse2.nb210nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb210nf_nouter(%rsp),%ebx
        movl %ebx,nb210nf_nouter(%rsp)

_nb_kernel210nf_x86_64_sse2.nb210nf_outer: 
        movq  nb210nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb210nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb210nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb210nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb210nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movq  nb210nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%edx
        imull nb210nf_ntype(%rsp),%edx
        shll  %edx
        movl  %edx,nb210nf_ntia(%rsp)

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb210nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb210nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb210nf_ix(%rsp)
        movapd %xmm1,nb210nf_iy(%rsp)
        movapd %xmm2,nb210nf_iz(%rsp)

        movl  %ebx,nb210nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb210nf_vctot(%rsp)
        movapd %xmm4,nb210nf_Vvdwtot(%rsp)

        movq  nb210nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb210nf_pos(%rbp),%rsi
        movq  nb210nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb210nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb210nf_ninner(%rsp),%ecx
        movl  %ecx,nb210nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb210nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel210nf_x86_64_sse2.nb210nf_unroll_loop
        jmp   _nb_kernel210nf_x86_64_sse2.nb210nf_checksingle
_nb_kernel210nf_x86_64_sse2.nb210nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb210nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb210nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb210nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb210nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movq nb210nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb210nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb210nf_ntia(%rsp),%edi
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
        movapd %xmm4,nb210nf_c6(%rsp)
        movapd %xmm6,nb210nf_c12(%rsp)

        movq nb210nf_pos(%rbp),%rsi        ## base of pos[] 

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
        movapd nb210nf_ix(%rsp),%xmm4
        movapd nb210nf_iy(%rsp),%xmm5
        movapd nb210nf_iz(%rsp),%xmm6

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

        movapd nb210nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb210nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb210nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb210nf_crf(%rsp),%xmm6
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulpd  nb210nf_c6(%rsp),%xmm1
        mulpd  nb210nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb210nf_Vvdwtot(%rsp),%xmm5
        addpd  nb210nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb210nf_vctot(%rsp)
        movapd %xmm5,nb210nf_Vvdwtot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb210nf_innerk(%rsp)
        jl    _nb_kernel210nf_x86_64_sse2.nb210nf_checksingle
        jmp   _nb_kernel210nf_x86_64_sse2.nb210nf_unroll_loop

_nb_kernel210nf_x86_64_sse2.nb210nf_checksingle: 
        movl  nb210nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel210nf_x86_64_sse2.nb210nf_dosingle
        jmp    _nb_kernel210nf_x86_64_sse2.nb210nf_updateouterdata
_nb_kernel210nf_x86_64_sse2.nb210nf_dosingle: 
        movq nb210nf_charge(%rbp),%rsi
        movq nb210nf_pos(%rbp),%rdi
        movq  nb210nf_innerjjnr(%rsp),%rcx
        xorpd %xmm3,%xmm3
        movl  (%rcx),%eax

        movlpd (%rsi,%rax,8),%xmm3
        movapd nb210nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movq nb210nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb210nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb210nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movapd %xmm4,nb210nf_c6(%rsp)
        movapd %xmm6,nb210nf_c12(%rsp)

        movq nb210nf_pos(%rbp),%rsi        ## base of pos[] 

        lea (%rax,%rax,2),%rax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb210nf_ix(%rsp),%xmm4
        movapd nb210nf_iy(%rsp),%xmm5
        movapd nb210nf_iz(%rsp),%xmm6

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

        movapd nb210nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb210nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb210nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb210nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb210nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb210nf_crf(%rsp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        mulsd  nb210nf_c6(%rsp),%xmm1
        mulsd  nb210nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb210nf_Vvdwtot(%rsp),%xmm5
        addsd  nb210nf_vctot(%rsp),%xmm6
        movlpd %xmm6,nb210nf_vctot(%rsp)
        movlpd %xmm5,nb210nf_Vvdwtot(%rsp)

_nb_kernel210nf_x86_64_sse2.nb210nf_updateouterdata: 
        ## get n from stack
        movl nb210nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb210nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb210nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb210nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb210nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb210nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb210nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel210nf_x86_64_sse2.nb210nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb210nf_n(%rsp)
        jmp _nb_kernel210nf_x86_64_sse2.nb210nf_outer
_nb_kernel210nf_x86_64_sse2.nb210nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb210nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel210nf_x86_64_sse2.nb210nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel210nf_x86_64_sse2.nb210nf_threadloop
_nb_kernel210nf_x86_64_sse2.nb210nf_end: 
        movl nb210nf_nouter(%rsp),%eax
        movl nb210nf_ninner(%rsp),%ebx
        movq nb210nf_outeriter(%rbp),%rcx
        movq nb210nf_inneriter(%rbp),%rdx
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

