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


.globl nb_kernel010_ia32_sse2
.globl _nb_kernel010_ia32_sse2
nb_kernel010_ia32_sse2: 
_nb_kernel010_ia32_sse2:        
.set nb010_p_nri, 8
.set nb010_iinr, 12
.set nb010_jindex, 16
.set nb010_jjnr, 20
.set nb010_shift, 24
.set nb010_shiftvec, 28
.set nb010_fshift, 32
.set nb010_gid, 36
.set nb010_pos, 40
.set nb010_faction, 44
.set nb010_charge, 48
.set nb010_p_facel, 52
.set nb010_argkrf, 56
.set nb010_argcrf, 60
.set nb010_Vc, 64
.set nb010_type, 68
.set nb010_p_ntype, 72
.set nb010_vdwparam, 76
.set nb010_Vvdw, 80
.set nb010_p_tabscale, 84
.set nb010_VFtab, 88
.set nb010_invsqrta, 92
.set nb010_dvda, 96
.set nb010_p_gbtabscale, 100
.set nb010_GBtab, 104
.set nb010_p_nthreads, 108
.set nb010_count, 112
.set nb010_mtx, 116
.set nb010_outeriter, 120
.set nb010_inneriter, 124
.set nb010_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010_ix, 0
.set nb010_iy, 16
.set nb010_iz, 32
.set nb010_dx, 48
.set nb010_dy, 64
.set nb010_dz, 80
.set nb010_two, 96
.set nb010_c6, 112
.set nb010_c12, 128
.set nb010_six, 144
.set nb010_twelve, 160
.set nb010_Vvdwtot, 176
.set nb010_fix, 192
.set nb010_fiy, 208
.set nb010_fiz, 224
.set nb010_half, 240
.set nb010_three, 256
.set nb010_is3, 272
.set nb010_ii3, 276
.set nb010_ntia, 280
.set nb010_innerjjnr, 284
.set nb010_innerk, 288
.set nb010_n, 292
.set nb010_nn1, 296
.set nb010_nri, 300
.set nb010_ntype, 304
.set nb010_nouter, 308
.set nb010_ninner, 312
.set nb010_salign, 316
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $320,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb010_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb010_p_nri(%ebp),%ecx
        movl nb010_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb010_nri(%esp)
        movl %edi,nb010_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010_nouter(%esp)
        movl %eax,nb010_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 2.0 IEEE (hex)
        movl $0x40000000,%ebx
        movl %eax,nb010_two(%esp)
        movl %ebx,nb010_two+4(%esp)
        movsd nb010_two(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm2
        addpd  %xmm1,%xmm2      ## 4.0
        addpd  %xmm1,%xmm2      ## 6.0
        movapd %xmm2,%xmm3
        addpd  %xmm3,%xmm3      ## 12.0
        movapd %xmm1,nb010_two(%esp)
        movapd %xmm2,nb010_six(%esp)
        movapd %xmm3,nb010_twelve(%esp)

_nb_kernel010_ia32_sse2.nb010_threadloop: 
        movl  nb010_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010_ia32_sse2.nb010_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010_ia32_sse2.nb010_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010_n(%esp)
        movl %ebx,nb010_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010_ia32_sse2.nb010_outerstart
        jmp _nb_kernel010_ia32_sse2.nb010_end

_nb_kernel010_ia32_sse2.nb010_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010_nouter(%esp),%ebx
        movl %ebx,nb010_nouter(%esp)

_nb_kernel010_ia32_sse2.nb010_outer: 
        movl  nb010_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb010_is3(%esp)      ## store is3 

        movl  nb010_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb010_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb010_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb010_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb010_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb010_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb010_ix(%esp)
        movapd %xmm1,nb010_iy(%esp)
        movapd %xmm2,nb010_iz(%esp)

        movl  %ebx,nb010_ii3(%esp)

        ## clear Vvdwtot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb010_Vvdwtot(%esp)
        movapd %xmm4,nb010_fix(%esp)
        movapd %xmm4,nb010_fiy(%esp)
        movapd %xmm4,nb010_fiz(%esp)

        movl  nb010_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb010_pos(%ebp),%esi
        movl  nb010_faction(%ebp),%edi
        movl  nb010_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010_ninner(%esp),%ecx
        movl  %ecx,nb010_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb010_innerk(%esp)      ## number of innerloop atoms 

        jge   _nb_kernel010_ia32_sse2.nb010_unroll_loop
        jmp   _nb_kernel010_ia32_sse2.nb010_checksingle
_nb_kernel010_ia32_sse2.nb010_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb010_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl  $8,nb010_innerjjnr(%esp)              ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb010_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb010_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb010_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx

        movapd %xmm4,nb010_c6(%esp)
        movapd %xmm6,nb010_c12(%esp)

        movl nb010_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    

        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010_ix(%esp),%xmm4
        movapd nb010_iy(%esp),%xmm5
        movapd nb010_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb010_dx(%esp)
        movapd %xmm5,nb010_dy(%esp)
        movapd %xmm6,nb010_dz(%esp)
        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4       ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm6
        rcpps %xmm6,%xmm6
        cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010_two(%esp),%xmm0
        movapd %xmm4,%xmm5
        mulpd %xmm6,%xmm4       ## lu*rsq 
        subpd %xmm4,%xmm0       ## 2-lu*rsq 
        mulpd %xmm0,%xmm6       ## (new lu) 

        movapd nb010_two(%esp),%xmm0
        mulpd %xmm6,%xmm5       ## lu*rsq 
        subpd %xmm5,%xmm0       ## 2-lu*rsq 
        mulpd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulpd  nb010_c6(%esp),%xmm1
        mulpd  nb010_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb010_Vvdwtot(%esp),%xmm5
        mulpd  nb010_six(%esp),%xmm1
        mulpd  nb010_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm2,%xmm0      ## xmm4=total fscal 
        movapd %xmm0,%xmm4

        movapd nb010_dx(%esp),%xmm0
        movapd nb010_dy(%esp),%xmm1
        movapd nb010_dz(%esp),%xmm2

        movapd %xmm5,nb010_Vvdwtot(%esp)

        movl   nb010_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb010_fix(%esp),%xmm3
        movapd nb010_fiy(%esp),%xmm4
        movapd nb010_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb010_fix(%esp)
        movapd %xmm4,nb010_fiy(%esp)
        movapd %xmm5,nb010_fiz(%esp)

        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        movhpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%ebx,8),%xmm4
        movhpd 16(%edi,%ebx,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)
        movhpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm4,8(%edi,%ebx,8)
        movhpd %xmm5,16(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl  $2,nb010_innerk(%esp)
        jl    _nb_kernel010_ia32_sse2.nb010_checksingle
        jmp   _nb_kernel010_ia32_sse2.nb010_unroll_loop
_nb_kernel010_ia32_sse2.nb010_checksingle:      
        movl  nb010_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010_ia32_sse2.nb010_dosingle
        jmp    _nb_kernel010_ia32_sse2.nb010_updateouterdata
_nb_kernel010_ia32_sse2.nb010_dosingle: 
        movl nb010_pos(%ebp),%edi
        movl  nb010_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage    
        movl nb010_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb010_vdwparam(%ebp),%esi
        shll %eax
        movl nb010_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb010_c6(%esp)
        movapd %xmm6,nb010_c12(%esp)

        movl nb010_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        

        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010_ix(%esp),%xmm4
        movapd nb010_iy(%esp),%xmm5
        movapd nb010_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb010_dx(%esp)
        movapd %xmm5,nb010_dy(%esp)
        movapd %xmm6,nb010_dz(%esp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4       ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm6
        rcpss %xmm6,%xmm6
        cvtss2sd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010_two(%esp),%xmm0
        movapd %xmm4,%xmm5
        mulsd %xmm6,%xmm4       ## lu*rsq 
        subsd %xmm4,%xmm0       ## 2-lu*rsq 
        mulsd %xmm0,%xmm6       ## (new lu) 

        movapd nb010_two(%esp),%xmm0
        mulsd %xmm6,%xmm5       ## lu*rsq 
        subsd %xmm5,%xmm0       ## 2-lu*rsq 
        mulsd %xmm6,%xmm0       ## xmm0=rinvsq 
        movapd %xmm0,%xmm4

        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulsd  nb010_c6(%esp),%xmm1
        mulsd  nb010_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb010_Vvdwtot(%esp),%xmm5
        mulsd  nb010_six(%esp),%xmm1
        mulsd  nb010_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm2,%xmm4      ## xmm4=total fscal 

        movapd nb010_dx(%esp),%xmm0
        movapd nb010_dy(%esp),%xmm1
        movapd nb010_dz(%esp),%xmm2

        movlpd %xmm5,nb010_Vvdwtot(%esp)

        movl   nb010_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb010_fix(%esp),%xmm3
        movlpd nb010_fiy(%esp),%xmm4
        movlpd nb010_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb010_fix(%esp)
        movlpd %xmm4,nb010_fiy(%esp)
        movlpd %xmm5,nb010_fiz(%esp)

        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel010_ia32_sse2.nb010_updateouterdata: 
        movl  nb010_ii3(%esp),%ecx
        movl  nb010_faction(%ebp),%edi
        movl  nb010_fshift(%ebp),%esi
        movl  nb010_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb010_fix(%esp),%xmm0
        movapd nb010_fiy(%esp),%xmm1
        movapd nb010_fiz(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## increment fshift force  
        movsd  (%esi,%edx,8),%xmm3
        movsd  8(%esi,%edx,8),%xmm4
        movsd  16(%esi,%edx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%esi,%edx,8)
        movsd  %xmm4,8(%esi,%edx,8)
        movsd  %xmm5,16(%esi,%edx,8)

        ## get n from stack
        movl nb010_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb010_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb010_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb010_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010_ia32_sse2.nb010_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010_n(%esp)
        jmp _nb_kernel010_ia32_sse2.nb010_outer
_nb_kernel010_ia32_sse2.nb010_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010_ia32_sse2.nb010_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010_ia32_sse2.nb010_threadloop
_nb_kernel010_ia32_sse2.nb010_end: 
        emms

        movl nb010_nouter(%esp),%eax
        movl nb010_ninner(%esp),%ebx
        movl nb010_outeriter(%ebp),%ecx
        movl nb010_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb010_salign(%esp),%eax
        addl %eax,%esp
        addl $320,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel010nf_ia32_sse2
.globl _nb_kernel010nf_ia32_sse2
nb_kernel010nf_ia32_sse2:       
_nb_kernel010nf_ia32_sse2:      
.set nb010nf_p_nri, 8
.set nb010nf_iinr, 12
.set nb010nf_jindex, 16
.set nb010nf_jjnr, 20
.set nb010nf_shift, 24
.set nb010nf_shiftvec, 28
.set nb010nf_fshift, 32
.set nb010nf_gid, 36
.set nb010nf_pos, 40
.set nb010nf_faction, 44
.set nb010nf_charge, 48
.set nb010nf_p_facel, 52
.set nb010nf_argkrf, 56
.set nb010nf_argcrf, 60
.set nb010nf_Vc, 64
.set nb010nf_type, 68
.set nb010nf_p_ntype, 72
.set nb010nf_vdwparam, 76
.set nb010nf_Vvdw, 80
.set nb010nf_p_tabscale, 84
.set nb010nf_VFtab, 88
.set nb010nf_invsqrta, 92
.set nb010nf_dvda, 96
.set nb010nf_p_gbtabscale, 100
.set nb010nf_GBtab, 104
.set nb010nf_p_nthreads, 108
.set nb010nf_count, 112
.set nb010nf_mtx, 116
.set nb010nf_outeriter, 120
.set nb010nf_inneriter, 124
.set nb010nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb010nf_ix, 0
.set nb010nf_iy, 16
.set nb010nf_iz, 32
.set nb010nf_two, 48
.set nb010nf_c6, 64
.set nb010nf_c12, 80
.set nb010nf_Vvdwtot, 96
.set nb010nf_half, 112
.set nb010nf_three, 128
.set nb010nf_is3, 144
.set nb010nf_ii3, 148
.set nb010nf_ntia, 152
.set nb010nf_innerjjnr, 156
.set nb010nf_innerk, 160
.set nb010nf_n, 164
.set nb010nf_nn1, 168
.set nb010nf_nri, 172
.set nb010nf_ntype, 176
.set nb010nf_nouter, 180
.set nb010nf_ninner, 184
.set nb010nf_salign, 188
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $192,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb010nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb010nf_p_nri(%ebp),%ecx
        movl nb010nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movl (%edi),%edi
        movl %ecx,nb010nf_nri(%esp)
        movl %edi,nb010nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb010nf_nouter(%esp)
        movl %eax,nb010nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 2.0 IEEE (hex)
        movl $0x40000000,%ebx
        movl %eax,nb010nf_two(%esp)
        movl %ebx,nb010nf_two+4(%esp)
        movsd nb010nf_two(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,nb010nf_two(%esp)

_nb_kernel010nf_ia32_sse2.nb010nf_threadloop: 
        movl  nb010nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel010nf_ia32_sse2.nb010nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel010nf_ia32_sse2.nb010nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb010nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb010nf_n(%esp)
        movl %ebx,nb010nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel010nf_ia32_sse2.nb010nf_outerstart
        jmp _nb_kernel010nf_ia32_sse2.nb010nf_end

_nb_kernel010nf_ia32_sse2.nb010nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb010nf_nouter(%esp),%ebx
        movl %ebx,nb010nf_nouter(%esp)

_nb_kernel010nf_ia32_sse2.nb010nf_outer: 
        movl  nb010nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb010nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb010nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb010nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb010nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb010nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb010nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb010nf_ix(%esp)
        movapd %xmm1,nb010nf_iy(%esp)
        movapd %xmm2,nb010nf_iz(%esp)

        movl  %ebx,nb010nf_ii3(%esp)

        ## clear Vvdwtot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb010nf_Vvdwtot(%esp)

        movl  nb010nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb010nf_pos(%ebp),%esi
        movl  nb010nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb010nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb010nf_ninner(%esp),%ecx
        movl  %ecx,nb010nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb010nf_innerk(%esp)      ## number of innerloop atoms 

        jge   _nb_kernel010nf_ia32_sse2.nb010nf_unroll_loop
        jmp   _nb_kernel010nf_ia32_sse2.nb010nf_checksingle
_nb_kernel010nf_ia32_sse2.nb010nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb010nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl  $8,nb010nf_innerjjnr(%esp)              ## advance pointer (unrolled 2) 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb010nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb010nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb010nf_ntia(%esp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movlpd (%esi,%ebx,8),%xmm7      ## c6b
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        movhpd 8(%esi,%ebx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx

        movapd %xmm4,nb010nf_c6(%esp)
        movapd %xmm6,nb010nf_c12(%esp)

        movl nb010nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move two coordinates to xmm0-xmm2    

        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2
        movhpd (%esi,%ebx,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm1
        movhpd 16(%esi,%ebx,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010nf_ix(%esp),%xmm4
        movapd nb010nf_iy(%esp),%xmm5
        movapd nb010nf_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4       ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm6
        rcpps %xmm6,%xmm6
        cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010nf_two(%esp),%xmm0
        movapd %xmm4,%xmm5
        mulpd %xmm6,%xmm4       ## lu*rsq 
        subpd %xmm4,%xmm0       ## 2-lu*rsq 
        mulpd %xmm0,%xmm6       ## (new lu) 

        movapd nb010nf_two(%esp),%xmm0
        mulpd %xmm6,%xmm5       ## lu*rsq 
        subpd %xmm5,%xmm0       ## 2-lu*rsq 
        mulpd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulpd  nb010nf_c6(%esp),%xmm1
        mulpd  nb010nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb010nf_Vvdwtot(%esp),%xmm5
        movapd %xmm5,nb010nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl  $2,nb010nf_innerk(%esp)
        jl    _nb_kernel010nf_ia32_sse2.nb010nf_checksingle
        jmp   _nb_kernel010nf_ia32_sse2.nb010nf_unroll_loop
_nb_kernel010nf_ia32_sse2.nb010nf_checksingle:  
        movl  nb010nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel010nf_ia32_sse2.nb010nf_dosingle
        jmp    _nb_kernel010nf_ia32_sse2.nb010nf_updateouterdata
_nb_kernel010nf_ia32_sse2.nb010nf_dosingle: 
        movl nb010nf_pos(%ebp),%edi
        movl  nb010nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        movd  %eax,%mm0         ## use mmx registers as temp storage    
        movl nb010nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb010nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb010nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb010nf_c6(%esp)
        movapd %xmm6,nb010nf_c12(%esp)

        movl nb010nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2        

        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb010nf_ix(%esp),%xmm4
        movapd nb010nf_iy(%esp),%xmm5
        movapd nb010nf_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4       ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm6
        rcpss %xmm6,%xmm6
        cvtss2sd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb010nf_two(%esp),%xmm0
        movapd %xmm4,%xmm5
        mulsd %xmm6,%xmm4       ## lu*rsq 
        subsd %xmm4,%xmm0       ## 2-lu*rsq 
        mulsd %xmm0,%xmm6       ## (new lu) 

        movapd nb010nf_two(%esp),%xmm0
        mulsd %xmm6,%xmm5       ## lu*rsq 
        subsd %xmm5,%xmm0       ## 2-lu*rsq 
        mulsd %xmm6,%xmm0       ## xmm0=rinvsq 
        movapd %xmm0,%xmm4

        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulsd  nb010nf_c6(%esp),%xmm1
        mulsd  nb010nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb010nf_Vvdwtot(%esp),%xmm5
        movlpd %xmm5,nb010nf_Vvdwtot(%esp)

_nb_kernel010nf_ia32_sse2.nb010nf_updateouterdata: 
        ## get n from stack
        movl nb010nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb010nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total lj energy and update it 
        movapd nb010nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 have the sum now 

        ## add earlier value from mem 
        movl  nb010nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb010nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel010nf_ia32_sse2.nb010nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb010nf_n(%esp)
        jmp _nb_kernel010nf_ia32_sse2.nb010nf_outer
_nb_kernel010nf_ia32_sse2.nb010nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb010nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel010nf_ia32_sse2.nb010nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel010nf_ia32_sse2.nb010nf_threadloop
_nb_kernel010nf_ia32_sse2.nb010nf_end: 
        emms

        movl nb010nf_nouter(%esp),%eax
        movl nb010nf_ninner(%esp),%ebx
        movl nb010nf_outeriter(%ebp),%ecx
        movl nb010nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb010nf_salign(%esp),%eax
        addl %eax,%esp
        addl $192,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret


