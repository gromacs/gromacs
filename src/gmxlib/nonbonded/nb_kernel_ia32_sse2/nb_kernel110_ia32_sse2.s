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



.globl nb_kernel110_ia32_sse2
.globl _nb_kernel110_ia32_sse2
nb_kernel110_ia32_sse2: 
_nb_kernel110_ia32_sse2:        
.set nb110_p_nri, 8
.set nb110_iinr, 12
.set nb110_jindex, 16
.set nb110_jjnr, 20
.set nb110_shift, 24
.set nb110_shiftvec, 28
.set nb110_fshift, 32
.set nb110_gid, 36
.set nb110_pos, 40
.set nb110_faction, 44
.set nb110_charge, 48
.set nb110_p_facel, 52
.set nb110_argkrf, 56
.set nb110_argcrf, 60
.set nb110_Vc, 64
.set nb110_type, 68
.set nb110_p_ntype, 72
.set nb110_vdwparam, 76
.set nb110_Vvdw, 80
.set nb110_p_tabscale, 84
.set nb110_VFtab, 88
.set nb110_invsqrta, 92
.set nb110_dvda, 96
.set nb110_p_gbtabscale, 100
.set nb110_GBtab, 104
.set nb110_p_nthreads, 108
.set nb110_count, 112
.set nb110_mtx, 116
.set nb110_outeriter, 120
.set nb110_inneriter, 124
.set nb110_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb110_ix, 0
.set nb110_iy, 16
.set nb110_iz, 32
.set nb110_iq, 48
.set nb110_dx, 64
.set nb110_dy, 80
.set nb110_dz, 96
.set nb110_c6, 112
.set nb110_c12, 128
.set nb110_six, 144
.set nb110_twelve, 160
.set nb110_vctot, 176
.set nb110_Vvdwtot, 192
.set nb110_fix, 208
.set nb110_fiy, 224
.set nb110_fiz, 240
.set nb110_half, 256
.set nb110_three, 272
.set nb110_facel, 288
.set nb110_is3, 304
.set nb110_ii3, 308
.set nb110_ntia, 312
.set nb110_innerjjnr, 316
.set nb110_innerk, 320
.set nb110_n, 324
.set nb110_nn1, 328
.set nb110_nri, 332
.set nb110_ntype, 336
.set nb110_nouter, 340
.set nb110_ninner, 344
.set nb110_salign, 348

        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $352,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb110_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb110_p_nri(%ebp),%ecx
        movl nb110_p_facel(%ebp),%esi
        movl nb110_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb110_nri(%esp)
        movsd %xmm7,nb110_facel(%esp)
        movl %edi,nb110_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110_nouter(%esp)
        movl %eax,nb110_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb110_half(%esp)
        movl %ebx,nb110_half+4(%esp)
        movsd nb110_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb110_half(%esp)
        movapd %xmm3,nb110_three(%esp)
        movapd %xmm4,nb110_six(%esp)
        movapd %xmm5,nb110_twelve(%esp)

_nb_kernel110_ia32_sse2.nb110_threadloop: 
        movl  nb110_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110_ia32_sse2.nb110_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110_ia32_sse2.nb110_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110_n(%esp)
        movl %ebx,nb110_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110_ia32_sse2.nb110_outerstart
        jmp _nb_kernel110_ia32_sse2.nb110_end

_nb_kernel110_ia32_sse2.nb110_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110_nouter(%esp),%ebx
        movl %ebx,nb110_nouter(%esp)

_nb_kernel110_ia32_sse2.nb110_outer: 
        movl  nb110_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb110_is3(%esp)      ## store is3 

        movl  nb110_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb110_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb110_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb110_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb110_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb110_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb110_ix(%esp)
        movapd %xmm1,nb110_iy(%esp)
        movapd %xmm2,nb110_iz(%esp)

        movl  %ebx,nb110_ii3(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb110_vctot(%esp)
        movapd %xmm4,nb110_Vvdwtot(%esp)
        movapd %xmm4,nb110_fix(%esp)
        movapd %xmm4,nb110_fiy(%esp)
        movapd %xmm4,nb110_fiz(%esp)

        movl  nb110_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110_pos(%ebp),%esi
        movl  nb110_faction(%ebp),%edi
        movl  nb110_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb110_ninner(%esp),%ecx
        movl  %ecx,nb110_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb110_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel110_ia32_sse2.nb110_unroll_loop
        jmp   _nb_kernel110_ia32_sse2.nb110_checksingle
_nb_kernel110_ia32_sse2.nb110_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb110_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb110_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb110_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb110_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb110_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb110_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb110_ntia(%esp),%edi
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
        movapd %xmm4,nb110_c6(%esp)
        movapd %xmm6,nb110_c12(%esp)

        movl nb110_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb110_ix(%esp),%xmm4
        movapd nb110_iy(%esp),%xmm5
        movapd nb110_iz(%esp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb110_dx(%esp)
        movapd %xmm5,nb110_dy(%esp)
        movapd %xmm6,nb110_dz(%esp)
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
        movapd nb110_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb110_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  nb110_c6(%esp),%xmm1
        mulpd  nb110_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb110_Vvdwtot(%esp),%xmm5
        mulpd  nb110_six(%esp),%xmm1
        mulpd  nb110_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        addpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm4      ## xmm4=total fscal 
        addpd  nb110_vctot(%esp),%xmm3

        movapd nb110_dx(%esp),%xmm0
        movapd nb110_dy(%esp),%xmm1
        movapd nb110_dz(%esp),%xmm2

        movapd %xmm3,nb110_vctot(%esp)
        movapd %xmm5,nb110_Vvdwtot(%esp)

        movl   nb110_faction(%ebp),%edi
        mulpd  %xmm4,%xmm0
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movapd nb110_fix(%esp),%xmm3
        movapd nb110_fiy(%esp),%xmm4
        movapd nb110_fiz(%esp),%xmm5
        addpd  %xmm0,%xmm3
        addpd  %xmm1,%xmm4
        addpd  %xmm2,%xmm5
        movapd %xmm3,nb110_fix(%esp)
        movapd %xmm4,nb110_fiy(%esp)
        movapd %xmm5,nb110_fiz(%esp)
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
        subl $2,nb110_innerk(%esp)
        jl    _nb_kernel110_ia32_sse2.nb110_checksingle
        jmp   _nb_kernel110_ia32_sse2.nb110_unroll_loop
_nb_kernel110_ia32_sse2.nb110_checksingle: 
        movl  nb110_innerk(%esp),%edx
        andl  $1,%edx
        jnz    _nb_kernel110_ia32_sse2.nb110_dosingle
        jmp    _nb_kernel110_ia32_sse2.nb110_updateouterdata
_nb_kernel110_ia32_sse2.nb110_dosingle: 
        movl nb110_charge(%ebp),%esi
        movl nb110_pos(%ebp),%edi
        movl nb110_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3

        movapd nb110_iq(%esp),%xmm5
        mulsd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb110_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb110_vdwparam(%ebp),%esi
        shll %eax
        movl nb110_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 
        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb110_c6(%esp)
        movapd %xmm6,nb110_c12(%esp)

        movl nb110_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb110_ix(%esp),%xmm4
        movapd nb110_iy(%esp),%xmm5
        movapd nb110_iz(%esp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## store dr 
        movapd %xmm4,nb110_dx(%esp)
        movapd %xmm5,nb110_dy(%esp)
        movapd %xmm6,nb110_dz(%esp)
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
        movapd nb110_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb110_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  nb110_c6(%esp),%xmm1
        mulsd  nb110_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb110_Vvdwtot(%esp),%xmm5
        mulsd  nb110_six(%esp),%xmm1
        mulsd  nb110_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        addsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm4      ## xmm4=total fscal 
        addsd  nb110_vctot(%esp),%xmm3

        movapd nb110_dx(%esp),%xmm0
        movapd nb110_dy(%esp),%xmm1
        movapd nb110_dz(%esp),%xmm2

        movlpd %xmm3,nb110_vctot(%esp)
        movlpd %xmm5,nb110_Vvdwtot(%esp)

        movl   nb110_faction(%ebp),%edi
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        ## now update f_i 
        movlpd nb110_fix(%esp),%xmm3
        movlpd nb110_fiy(%esp),%xmm4
        movlpd nb110_fiz(%esp),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movlpd %xmm3,nb110_fix(%esp)
        movlpd %xmm4,nb110_fiy(%esp)
        movlpd %xmm5,nb110_fiz(%esp)
        ## the fj's - start by accumulating forces from memory 
        movlpd (%edi,%eax,8),%xmm3
        movlpd 8(%edi,%eax,8),%xmm4
        movlpd 16(%edi,%eax,8),%xmm5
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        movlpd %xmm3,(%edi,%eax,8)
        movlpd %xmm4,8(%edi,%eax,8)
        movlpd %xmm5,16(%edi,%eax,8)

_nb_kernel110_ia32_sse2.nb110_updateouterdata: 
        movl  nb110_ii3(%esp),%ecx
        movl  nb110_faction(%ebp),%edi
        movl  nb110_fshift(%ebp),%esi
        movl  nb110_is3(%esp),%edx

        ## accumulate i forces in xmm0, xmm1, xmm2 
        movapd nb110_fix(%esp),%xmm0
        movapd nb110_fiy(%esp),%xmm1
        movapd nb110_fiz(%esp),%xmm2

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
        movl nb110_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb110_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb110_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb110_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb110_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

       ## finish if last 
        movl nb110_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110_ia32_sse2.nb110_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110_n(%esp)
        jmp _nb_kernel110_ia32_sse2.nb110_outer
_nb_kernel110_ia32_sse2.nb110_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110_ia32_sse2.nb110_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110_ia32_sse2.nb110_threadloop
_nb_kernel110_ia32_sse2.nb110_end: 
        emms

        movl nb110_nouter(%esp),%eax
        movl nb110_ninner(%esp),%ebx
        movl nb110_outeriter(%ebp),%ecx
        movl nb110_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb110_salign(%esp),%eax
        addl %eax,%esp
        addl $352,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret






.globl nb_kernel110nf_ia32_sse2
.globl _nb_kernel110nf_ia32_sse2
nb_kernel110nf_ia32_sse2:       
_nb_kernel110nf_ia32_sse2:      
.set nb110nf_p_nri, 8
.set nb110nf_iinr, 12
.set nb110nf_jindex, 16
.set nb110nf_jjnr, 20
.set nb110nf_shift, 24
.set nb110nf_shiftvec, 28
.set nb110nf_fshift, 32
.set nb110nf_gid, 36
.set nb110nf_pos, 40
.set nb110nf_faction, 44
.set nb110nf_charge, 48
.set nb110nf_p_facel, 52
.set nb110nf_argkrf, 56
.set nb110nf_argcrf, 60
.set nb110nf_Vc, 64
.set nb110nf_type, 68
.set nb110nf_p_ntype, 72
.set nb110nf_vdwparam, 76
.set nb110nf_Vvdw, 80
.set nb110nf_p_tabscale, 84
.set nb110nf_VFtab, 88
.set nb110nf_invsqrta, 92
.set nb110nf_dvda, 96
.set nb110nf_p_gbtabscale, 100
.set nb110nf_GBtab, 104
.set nb110nf_p_nthreads, 108
.set nb110nf_count, 112
.set nb110nf_mtx, 116
.set nb110nf_outeriter, 120
.set nb110nf_inneriter, 124
.set nb110nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb110nf_ix, 0
.set nb110nf_iy, 16
.set nb110nf_iz, 32
.set nb110nf_iq, 48
.set nb110nf_c6, 64
.set nb110nf_c12, 80
.set nb110nf_vctot, 96
.set nb110nf_Vvdwtot, 112
.set nb110nf_half, 128
.set nb110nf_three, 144
.set nb110nf_is3, 160
.set nb110nf_ii3, 164
.set nb110nf_ntia, 168
.set nb110nf_innerjjnr, 172
.set nb110nf_innerk, 176
.set nb110nf_n, 180
.set nb110nf_nn1, 184
.set nb110nf_nri, 188
.set nb110nf_facel, 192                       ## uses 8 bytes
.set nb110nf_ntype, 200
.set nb110nf_nouter, 204
.set nb110nf_ninner, 208
.set nb110nf_salign, 212
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
        movl %eax,nb110nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb110nf_p_nri(%ebp),%ecx
        movl nb110nf_p_facel(%ebp),%esi
        movl nb110nf_p_ntype(%ebp),%edi
        movl (%ecx),%ecx
        movsd (%esi),%xmm7
        movl (%edi),%edi
        movl %ecx,nb110nf_nri(%esp)
        movsd %xmm7,nb110nf_facel(%esp)
        movl %edi,nb110nf_ntype(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb110nf_nouter(%esp)
        movl %eax,nb110nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb110nf_half(%esp)
        movl %ebx,nb110nf_half+4(%esp)
        movsd nb110nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb110nf_half(%esp)
        movapd %xmm3,nb110nf_three(%esp)

_nb_kernel110nf_ia32_sse2.nb110nf_threadloop: 
        movl  nb110nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel110nf_ia32_sse2.nb110nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel110nf_ia32_sse2.nb110nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb110nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb110nf_n(%esp)
        movl %ebx,nb110nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel110nf_ia32_sse2.nb110nf_outerstart
        jmp _nb_kernel110nf_ia32_sse2.nb110nf_end

_nb_kernel110nf_ia32_sse2.nb110nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb110nf_nouter(%esp),%ebx
        movl %ebx,nb110nf_nouter(%esp)

_nb_kernel110nf_ia32_sse2.nb110nf_outer: 
        movl  nb110nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb110nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb110nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movl  nb110nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        mulsd nb110nf_facel(%esp),%xmm3
        shufpd $0,%xmm3,%xmm3

        movl  nb110nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%edx
        imull nb110nf_ntype(%esp),%edx
        shll  %edx
        movl  %edx,nb110nf_ntia(%esp)

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb110nf_pos(%ebp),%eax      ## eax = base of pos[]  

        addsd (%eax,%ebx,8),%xmm0
        addsd 8(%eax,%ebx,8),%xmm1
        addsd 16(%eax,%ebx,8),%xmm2

        movapd %xmm3,nb110nf_iq(%esp)

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2

        movapd %xmm0,nb110nf_ix(%esp)
        movapd %xmm1,nb110nf_iy(%esp)
        movapd %xmm2,nb110nf_iz(%esp)

        movl  %ebx,nb110nf_ii3(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb110nf_vctot(%esp)
        movapd %xmm4,nb110nf_Vvdwtot(%esp)

        movl  nb110nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb110nf_pos(%ebp),%esi
        movl  nb110nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb110nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb110nf_ninner(%esp),%ecx
        movl  %ecx,nb110nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb110nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel110nf_ia32_sse2.nb110nf_unroll_loop
        jmp   _nb_kernel110nf_ia32_sse2.nb110nf_checksingle
_nb_kernel110nf_ia32_sse2.nb110nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb110nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx
        addl $8,nb110nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb110nf_charge(%ebp),%esi     ## base of charge[] 

        movlpd (%esi,%eax,8),%xmm3
        movhpd (%esi,%ebx,8),%xmm3

        movapd nb110nf_iq(%esp),%xmm5
        mulpd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movl nb110nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl (%esi,%ebx,4),%ebx
        movl nb110nf_vdwparam(%ebp),%esi
        shll %eax
        shll %ebx
        movl nb110nf_ntia(%esp),%edi
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
        movapd %xmm4,nb110nf_c6(%esp)
        movapd %xmm6,nb110nf_c12(%esp)

        movl nb110nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd nb110nf_ix(%esp),%xmm4
        movapd nb110nf_iy(%esp),%xmm5
        movapd nb110nf_iz(%esp),%xmm6

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
        movapd nb110nf_three(%esp),%xmm1
        mulpd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110nf_half(%esp),%xmm0
        subpd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm1
        mulpd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulpd %xmm1,%xmm1       ## lu*lu 
        movapd nb110nf_three(%esp),%xmm2
        mulpd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110nf_half(%esp),%xmm0
        subpd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulpd %xmm5,%xmm2
        mulpd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulpd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulpd  %xmm4,%xmm1
        mulpd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulpd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulpd  nb110nf_c6(%esp),%xmm1
        mulpd  nb110nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addpd  nb110nf_Vvdwtot(%esp),%xmm5
        addpd  nb110nf_vctot(%esp),%xmm3
        movapd %xmm3,nb110nf_vctot(%esp)
        movapd %xmm5,nb110nf_Vvdwtot(%esp)

        ## should we do one more iteration? 
        subl $2,nb110nf_innerk(%esp)
        jl    _nb_kernel110nf_ia32_sse2.nb110nf_checksingle
        jmp   _nb_kernel110nf_ia32_sse2.nb110nf_unroll_loop
_nb_kernel110nf_ia32_sse2.nb110nf_checksingle: 
        movl  nb110nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel110nf_ia32_sse2.nb110nf_dosingle
        jmp   _nb_kernel110nf_ia32_sse2.nb110nf_updateouterdata
_nb_kernel110nf_ia32_sse2.nb110nf_dosingle: 
        movl nb110nf_charge(%ebp),%esi
        movl nb110nf_pos(%ebp),%edi
        movl nb110nf_innerjjnr(%esp),%ecx
        movl  (%ecx),%eax

        xorpd %xmm3,%xmm3
        movlpd (%esi,%eax,8),%xmm3

        movapd nb110nf_iq(%esp),%xmm5
        mulsd %xmm5,%xmm3               ## qq 

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movl nb110nf_type(%ebp),%esi
        movl (%esi,%eax,4),%eax
        movl nb110nf_vdwparam(%ebp),%esi
        shll %eax
        movl nb110nf_ntia(%esp),%edi
        addl %edi,%eax

        movlpd (%esi,%eax,8),%xmm6      ## c6a
        movhpd 8(%esi,%eax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax

        movapd %xmm4,nb110nf_c6(%esp)
        movapd %xmm6,nb110nf_c12(%esp)

        movl nb110nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move two coordinates to xmm0-xmm2    
        movlpd (%esi,%eax,8),%xmm0
        movlpd 8(%esi,%eax,8),%xmm1
        movlpd 16(%esi,%eax,8),%xmm2

        ## move ix-iz to xmm4-xmm6 
        movapd nb110nf_ix(%esp),%xmm4
        movapd nb110nf_iy(%esp),%xmm5
        movapd nb110nf_iz(%esp),%xmm6

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
        movapd nb110nf_three(%esp),%xmm1
        mulsd %xmm4,%xmm2       ## rsq*lu*lu                    
        movapd nb110nf_half(%esp),%xmm0
        subsd %xmm2,%xmm1       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm1
        mulsd %xmm0,%xmm1       ## xmm0=iter1 of rinv (new lu) 

        movapd %xmm1,%xmm5      ## copy of lu 
        mulsd %xmm1,%xmm1       ## lu*lu 
        movapd nb110nf_three(%esp),%xmm2
        mulsd %xmm4,%xmm1       ## rsq*lu*lu                    
        movapd nb110nf_half(%esp),%xmm0
        subsd %xmm1,%xmm2       ## 30-rsq*lu*lu 
        mulsd %xmm5,%xmm2
        mulsd %xmm2,%xmm0       ## xmm0=rinv 

        movapd %xmm0,%xmm4
        mulsd  %xmm4,%xmm4      ## xmm4=rinvsq 
        movapd %xmm4,%xmm1
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  %xmm0,%xmm3      ## xmm3=vcoul 
        mulsd  nb110nf_c6(%esp),%xmm1
        mulsd  nb110nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        addsd  nb110nf_Vvdwtot(%esp),%xmm5
        addsd  nb110nf_vctot(%esp),%xmm3
        movlpd %xmm3,nb110nf_vctot(%esp)
        movlpd %xmm5,nb110nf_Vvdwtot(%esp)

_nb_kernel110nf_ia32_sse2.nb110nf_updateouterdata: 
        ## get n from stack
        movl nb110nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb110nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb110nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb110nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb110nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb110nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb110nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel110nf_ia32_sse2.nb110nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb110nf_n(%esp)
        jmp _nb_kernel110nf_ia32_sse2.nb110nf_outer
_nb_kernel110nf_ia32_sse2.nb110nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb110nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel110nf_ia32_sse2.nb110nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel110nf_ia32_sse2.nb110nf_threadloop
_nb_kernel110nf_ia32_sse2.nb110nf_end: 
        emms

        movl nb110nf_nouter(%esp),%eax
        movl nb110nf_ninner(%esp),%ebx
        movl nb110nf_outeriter(%ebp),%ecx
        movl nb110nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb110nf_salign(%esp),%eax
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



