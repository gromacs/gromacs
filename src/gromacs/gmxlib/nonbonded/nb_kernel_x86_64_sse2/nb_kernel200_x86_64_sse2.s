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





.globl nb_kernel200_x86_64_sse2
.globl _nb_kernel200_x86_64_sse2
nb_kernel200_x86_64_sse2:       
_nb_kernel200_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb200_fshift, 16
.set nb200_gid, 24
.set nb200_pos, 32
.set nb200_faction, 40
.set nb200_charge, 48
.set nb200_p_facel, 56
.set nb200_argkrf, 64
.set nb200_argcrf, 72
.set nb200_Vc, 80
.set nb200_type, 88
.set nb200_p_ntype, 96
.set nb200_vdwparam, 104
.set nb200_Vvdw, 112
.set nb200_p_tabscale, 120
.set nb200_VFtab, 128
.set nb200_invsqrta, 136
.set nb200_dvda, 144
.set nb200_p_gbtabscale, 152
.set nb200_GBtab, 160
.set nb200_p_nthreads, 168
.set nb200_count, 176
.set nb200_mtx, 184
.set nb200_outeriter, 192
.set nb200_inneriter, 200
.set nb208_work, 200
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb200_ix, 0
.set nb200_iy, 16
.set nb200_iz, 32
.set nb200_iq, 48
.set nb200_dx, 64
.set nb200_dy, 80
.set nb200_dz, 96
.set nb200_vctot, 112
.set nb200_fix, 128
.set nb200_fiy, 144
.set nb200_fiz, 160
.set nb200_half, 176
.set nb200_three, 192
.set nb200_two, 208
.set nb200_krf, 224
.set nb200_crf, 240
.set nb200_is3, 256
.set nb200_ii3, 260
.set nb200_nri, 264
.set nb200_iinr, 272
.set nb200_jindex, 280
.set nb200_jjnr, 288
.set nb200_shift, 296
.set nb200_shiftvec, 304
.set nb200_facel, 312
.set nb200_innerjjnr, 320
.set nb200_innerk, 328
.set nb200_n, 332
.set nb200_nn1, 336
.set nb200_nouter, 340
.set nb200_ninner, 344
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
        movl %eax,nb200_nouter(%rsp)
        movl %eax,nb200_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb200_nri(%rsp)
        movq %rsi,nb200_iinr(%rsp)
        movq %rdx,nb200_jindex(%rsp)
        movq %rcx,nb200_jjnr(%rsp)
        movq %r8,nb200_shift(%rsp)
        movq %r9,nb200_shiftvec(%rsp)
        movq nb200_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb200_facel(%rsp)

        movq nb200_argkrf(%rbp),%rsi
        movq nb200_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb200_krf(%rsp)
        movapd %xmm2,nb200_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb200_half(%rsp)
        movl %ebx,nb200_half+4(%rsp)
        movsd nb200_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb200_half(%rsp)
        movapd %xmm2,nb200_two(%rsp)
        movapd %xmm3,nb200_three(%rsp)

_nb_kernel200_x86_64_sse2.nb200_threadloop: 
        movq  nb200_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel200_x86_64_sse2.nb200_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200_x86_64_sse2.nb200_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200_n(%rsp)
        movl %ebx,nb200_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200_x86_64_sse2.nb200_outerstart
        jmp _nb_kernel200_x86_64_sse2.nb200_end

_nb_kernel200_x86_64_sse2.nb200_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200_nouter(%rsp),%ebx
        movl %ebx,nb200_nouter(%rsp)

_nb_kernel200_x86_64_sse2.nb200_outer: 
        movq  nb200_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb200_is3(%rsp)      ## store is3 

        movq  nb200_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb200_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb200_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb200_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb200_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb200_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb200_ix(%rsp)
        movapd %xmm1,nb200_iy(%rsp)
        movapd %xmm2,nb200_iz(%rsp)

        movl  %ebx,nb200_ii3(%rsp)

        ## clear vctot (xmm12) and i forces (xmm13-xmm15)
        xorpd %xmm12,%xmm12
        movapd %xmm12,%xmm13
        movapd %xmm12,%xmm14
        movapd %xmm12,%xmm15

        movq  nb200_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb200_pos(%rbp),%rsi
        movq  nb200_faction(%rbp),%rdi
        movq  nb200_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb200_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb200_ninner(%rsp),%ecx
        movl  %ecx,nb200_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb200_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel200_x86_64_sse2.nb200_unroll_loop
        jmp   _nb_kernel200_x86_64_sse2.nb200_checksingle
_nb_kernel200_x86_64_sse2.nb200_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb200_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        addq $8,nb200_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        lea  (%r8,%r8,2),%rax     ## replace jnr with j3 
        lea  (%r9,%r9,2),%rbx

        movq nb200_pos(%rbp),%rsi        ## base of pos[] 

        ## move two coordinates to xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
        movhpd (%rsi,%rbx,8),%xmm4
        movhpd 8(%rsi,%rbx,8),%xmm5
        movhpd 16(%rsi,%rbx,8),%xmm6

        ## calc dr 
        subpd nb200_ix(%rsp),%xmm4
        subpd nb200_iy(%rsp),%xmm5
        subpd nb200_iz(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,%xmm9
        movapd %xmm5,%xmm10
        movapd %xmm6,%xmm11

        movq nb200_charge(%rbp),%rsi     ## base of charge[] 
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        ## rsq in xmm4 


        movlpd (%rsi,%r8,8),%xmm3

        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm5,%xmm2    ## lu in low xmm2 

        movhpd (%rsi,%r9,8),%xmm3
        movapd nb200_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb200_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        mulpd nb200_iq(%rsp),%xmm3              ## qq 
        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb200_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb200_crf(%rsp),%xmm6
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 

        movq   nb200_faction(%rbp),%rdi

        addpd  %xmm7,%xmm7

        subpd  %xmm7,%xmm0
        mulpd  %xmm0,%xmm3
        mulpd  %xmm3,%xmm4      ## xmm4=total fscal 

    ## increment vctot
        addpd  %xmm6,%xmm12

        mulpd  %xmm4,%xmm9
        mulpd  %xmm4,%xmm10
        mulpd  %xmm4,%xmm11

        ## the fj's - start by accumulating forces from memory 
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5

        ## now update f_i 
        addpd  %xmm9,%xmm13
        addpd  %xmm10,%xmm14
        addpd  %xmm11,%xmm15

        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5

        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb200_innerk(%rsp)
        jl    _nb_kernel200_x86_64_sse2.nb200_checksingle
        jmp   _nb_kernel200_x86_64_sse2.nb200_unroll_loop

_nb_kernel200_x86_64_sse2.nb200_checksingle:    
        movl  nb200_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200_x86_64_sse2.nb200_dosingle
        jmp    _nb_kernel200_x86_64_sse2.nb200_updateouterdata
_nb_kernel200_x86_64_sse2.nb200_dosingle: 
        movq nb200_charge(%rbp),%rsi
        movq nb200_pos(%rbp),%rdi
        movq  nb200_innerjjnr(%rsp),%rcx
        movl  (%rcx),%eax

        movq nb200_charge(%rbp),%rsi     ## base of charge[] 

        movsd (%rsi,%rax,8),%xmm3
        mulsd nb200_iq(%rsp),%xmm3              ## qq 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        movq nb200_pos(%rbp),%rsi        ## base of pos[] 

        ## move two coordinates to xmm4-xmm6
        movsd (%rsi,%rax,8),%xmm4
        movsd 8(%rsi,%rax,8),%xmm5
        movsd 16(%rsi,%rax,8),%xmm6

        ## calc dr 
        subsd nb200_ix(%rsp),%xmm4
        subsd nb200_iy(%rsp),%xmm5
        subsd nb200_iz(%rsp),%xmm6

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

        movapd nb200_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb200_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb200_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb200_crf(%rsp),%xmm6
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 

        addsd  %xmm7,%xmm7

        subsd  %xmm7,%xmm0
        mulsd  %xmm0,%xmm3
        mulsd  %xmm3,%xmm4      ## xmm4=total fscal 

    ## increment vctot
        addsd  %xmm6,%xmm12

        movq   nb200_faction(%rbp),%rdi
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

_nb_kernel200_x86_64_sse2.nb200_updateouterdata: 
        movl  nb200_ii3(%rsp),%ecx
        movq  nb200_faction(%rbp),%rdi
        movq  nb200_fshift(%rbp),%rsi
        movl  nb200_is3(%rsp),%edx

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
        movl nb200_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb200_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total coulomb energy and update it 
        movhlps %xmm12,%xmm6
        addsd  %xmm6,%xmm12     ## low xmm12 have the sum now 

        ## add earlier value from mem 
        movq  nb200_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm12
        ## move back to mem 
        movsd %xmm12,(%rax,%rdx,8)

        ## finish if last 
        movl nb200_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200_x86_64_sse2.nb200_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200_n(%rsp)
        jmp _nb_kernel200_x86_64_sse2.nb200_outer
_nb_kernel200_x86_64_sse2.nb200_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200_x86_64_sse2.nb200_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200_x86_64_sse2.nb200_threadloop
_nb_kernel200_x86_64_sse2.nb200_end: 
        movl nb200_nouter(%rsp),%eax
        movl nb200_ninner(%rsp),%ebx
        movq nb200_outeriter(%rbp),%rcx
        movq nb200_inneriter(%rbp),%rdx
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




.globl nb_kernel200nf_x86_64_sse2
.globl _nb_kernel200nf_x86_64_sse2
nb_kernel200nf_x86_64_sse2:     
_nb_kernel200nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb200nf_fshift, 16
.set nb200nf_gid, 24
.set nb200nf_pos, 32
.set nb200nf_faction, 40
.set nb200nf_charge, 48
.set nb200nf_p_facel, 56
.set nb200nf_argkrf, 64
.set nb200nf_argcrf, 72
.set nb200nf_Vc, 80
.set nb200nf_type, 88
.set nb200nf_p_ntype, 96
.set nb200nf_vdwparam, 104
.set nb200nf_Vvdw, 112
.set nb200nf_p_tabscale, 120
.set nb200nf_VFtab, 128
.set nb200nf_invsqrta, 136
.set nb200nf_dvda, 144
.set nb200nf_p_gbtabscale, 152
.set nb200nf_GBtab, 160
.set nb200nf_p_nthreads, 168
.set nb200nf_count, 176
.set nb200nf_mtx, 184
.set nb200nf_outeriter, 192
.set nb200nf_inneriter, 200
.set nb208nf_work, 200
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb200nf_ix, 0
.set nb200nf_iy, 16
.set nb200nf_iz, 32
.set nb200nf_iq, 48
.set nb200nf_vctot, 64
.set nb200nf_half, 80
.set nb200nf_three, 96
.set nb200nf_krf, 112
.set nb200nf_crf, 128
.set nb200nf_is3, 144
.set nb200nf_ii3, 148
.set nb200nf_nri, 152
.set nb200nf_iinr, 160
.set nb200nf_jindex, 168
.set nb200nf_jjnr, 176
.set nb200nf_shift, 184
.set nb200nf_shiftvec, 192
.set nb200nf_facel, 200
.set nb200nf_innerjjnr, 208
.set nb200nf_innerk, 216
.set nb200nf_n, 220
.set nb200nf_nn1, 224
.set nb200nf_nouter, 228
.set nb200nf_ninner, 232
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $248,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb200nf_nouter(%rsp)
        movl %eax,nb200nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb200nf_nri(%rsp)
        movq %rsi,nb200nf_iinr(%rsp)
        movq %rdx,nb200nf_jindex(%rsp)
        movq %rcx,nb200nf_jjnr(%rsp)
        movq %r8,nb200nf_shift(%rsp)
        movq %r9,nb200nf_shiftvec(%rsp)
        movq nb200nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb200nf_facel(%rsp)

        movq nb200nf_argkrf(%rbp),%rsi
        movq nb200nf_argcrf(%rbp),%rdi
        movsd (%rsi),%xmm1
        movsd (%rdi),%xmm2
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        movapd %xmm1,nb200nf_krf(%rsp)
        movapd %xmm2,nb200nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb200nf_half(%rsp)
        movl %ebx,nb200nf_half+4(%rsp)
        movsd nb200nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb200nf_half(%rsp)
        movapd %xmm3,nb200nf_three(%rsp)

_nb_kernel200nf_x86_64_sse2.nb200nf_threadloop: 
        movq  nb200nf_count(%rbp),%rsi          ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel200nf_x86_64_sse2.nb200nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel200nf_x86_64_sse2.nb200nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb200nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb200nf_n(%rsp)
        movl %ebx,nb200nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel200nf_x86_64_sse2.nb200nf_outerstart
        jmp _nb_kernel200nf_x86_64_sse2.nb200nf_end

_nb_kernel200nf_x86_64_sse2.nb200nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb200nf_nouter(%rsp),%ebx
        movl %ebx,nb200nf_nouter(%rsp)

_nb_kernel200nf_x86_64_sse2.nb200nf_outer: 
        movq  nb200nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 

        movq  nb200nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb200nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movq  nb200nf_charge(%rbp),%rdx
        movsd (%rdx,%rbx,8),%xmm3
        mulsd nb200nf_facel(%rsp),%xmm3
        shufpd $0,%xmm3,%xmm3

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb200nf_pos(%rbp),%rax      ## rax = base of pos[]  

        addsd (%rax,%rbx,8),%xmm0
        addsd 8(%rax,%rbx,8),%xmm1
        addsd 16(%rax,%rbx,8),%xmm2

        movapd %xmm3,nb200nf_iq(%rsp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb200nf_ix(%rsp)
        movapd %xmm1,nb200nf_iy(%rsp)
        movapd %xmm2,nb200nf_iz(%rsp)

        movl  %ebx,nb200nf_ii3(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb200nf_vctot(%rsp)

        movq  nb200nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb200nf_pos(%rbp),%rsi
        movq  nb200nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb200nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb200nf_ninner(%rsp),%ecx
        movl  %ecx,nb200nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb200nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel200nf_x86_64_sse2.nb200nf_unroll_loop
        jmp   _nb_kernel200nf_x86_64_sse2.nb200nf_checksingle
_nb_kernel200nf_x86_64_sse2.nb200nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb200nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        addq $8,nb200nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb200nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3

        movapd nb200nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        movq nb200nf_pos(%rbp),%rsi        ## base of pos[] 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200nf_ix(%rsp),%xmm4
        movapd nb200nf_iy(%rsp),%xmm5
        movapd nb200nf_iz(%rsp),%xmm6

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

        movapd nb200nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulpd %xmm2,%xmm2       ## lu*lu 
        movapd nb200nf_three(%rsp),%xmm1
        mulpd %xmm4,%xmm7       ## krsq 
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200nf_half(%rsp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb200nf_three(%rsp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200nf_half(%rsp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subpd  nb200nf_crf(%rsp),%xmm6
        mulpd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addpd  nb200nf_vctot(%rsp),%xmm6
        movapd %xmm6,nb200nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb200nf_innerk(%rsp)
        jl    _nb_kernel200nf_x86_64_sse2.nb200nf_checksingle
        jmp   _nb_kernel200nf_x86_64_sse2.nb200nf_unroll_loop

_nb_kernel200nf_x86_64_sse2.nb200nf_checksingle: 
        movl  nb200nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz    _nb_kernel200nf_x86_64_sse2.nb200nf_dosingle
        jmp    _nb_kernel200nf_x86_64_sse2.nb200nf_updateouterdata
_nb_kernel200nf_x86_64_sse2.nb200nf_dosingle: 
        movq nb200nf_charge(%rbp),%rsi
        movq nb200nf_pos(%rbp),%rdi
        movq  nb200nf_innerjjnr(%rsp),%rcx

        xorpd %xmm3,%xmm3
        movl  (%rcx),%eax

        movlpd (%rsi,%rax,8),%xmm3
        movapd nb200nf_iq(%rsp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movq nb200nf_pos(%rbp),%rsi        ## base of pos[] 

        lea (%rax,%rax,2),%rax    ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb200nf_ix(%rsp),%xmm4
        movapd nb200nf_iy(%rsp),%xmm5
        movapd nb200nf_iz(%rsp),%xmm6

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

        movapd nb200nf_krf(%rsp),%xmm7
        ## lookup seed in xmm2 
        movapd %xmm2,%xmm5      ## copy of lu 
        mulsd %xmm2,%xmm2       ## lu*lu 
        movapd nb200nf_three(%rsp),%xmm1
        mulsd %xmm4,%xmm7       ## krsq 
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb200nf_half(%rsp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb200nf_three(%rsp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb200nf_half(%rsp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 
        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm0,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        movapd %xmm4,%xmm1
        subsd  nb200nf_crf(%rsp),%xmm6
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm3,%xmm6      ## xmm6=vcoul=qq*(rinv+ krsq) 
        addsd  nb200nf_vctot(%rsp),%xmm6
        movlpd %xmm6,nb200nf_vctot(%rsp)

_nb_kernel200nf_x86_64_sse2.nb200nf_updateouterdata: 
        ## get n from stack
        movl nb200nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb200nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb200nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb200nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb200nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel200nf_x86_64_sse2.nb200nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb200nf_n(%rsp)
        jmp _nb_kernel200nf_x86_64_sse2.nb200nf_outer
_nb_kernel200nf_x86_64_sse2.nb200nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb200nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel200nf_x86_64_sse2.nb200nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel200nf_x86_64_sse2.nb200nf_threadloop
_nb_kernel200nf_x86_64_sse2.nb200nf_end: 
        movl nb200nf_nouter(%rsp),%eax
        movl nb200nf_ninner(%rsp),%ebx
        movq nb200nf_outeriter(%rbp),%rcx
        movq nb200nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $248,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




